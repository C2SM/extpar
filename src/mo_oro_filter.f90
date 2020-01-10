!+  Fortran module with routines and settings for orography smoothing
!
! Description:
! Fortran module with routines and settings for orography smoothing
!
! Current Code Owner: MeteoSwiss, Anne Roches
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_4         2011/04/21 Anne Roches
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with routines and settings for orography smoothing
!> \author Anne Roches
!>
MODULE mo_oro_filter

  !###
  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp, i8, i4

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar,      ONLY: abort_extpar,           &
       &                              extend_field2D,         &
       &                              horizontal_filtering
  USE mo_logging

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_namelists_extpar_orosmooth, do_orosmooth

CONTAINS
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !> subroutine to read namelist for orography smoothing in EXTPAR
  SUBROUTINE read_namelists_extpar_orosmooth(namelist_file,        &
       lfilter_oro,          &
       ilow_pass_oro,        &
       numfilt_oro,          &
       eps_filter,           &
       ifill_valley,         &
       rfill_valley,         &
       ilow_pass_xso,        &
       numfilt_xso,          &
       lxso_first,           &
       rxso_mask)

    USE mo_utilities_extpar, ONLY: free_un, & ! function to get free unit number
         abort_extpar

    CHARACTER (len=*), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings
    ! orography smoothing

    LOGICAL,          INTENT(OUT) :: lfilter_oro              !< oro smoothing to be performed? (TRUE/FALSE)
    INTEGER(KIND=i4), INTENT(OUT) :: ilow_pass_oro            !< type of oro smoothing and 
    !  stencil width (1,4,5,6,8)
    INTEGER(KIND=i4), INTENT(OUT) :: numfilt_oro              !< number of applications of the filter
    REAL(KIND=wp),    INTENT(OUT) :: eps_filter               !< smoothing param ("strength" of the filtering)
    INTEGER(KIND=i4), INTENT(OUT) :: ifill_valley             !< fill valleys before or after oro smoothing 
    !  (1: before, 2: after)
    REAL(KIND=wp),    INTENT(OUT) :: rfill_valley             !< mask for valley filling (threshold value)
    INTEGER(KIND=i4), INTENT(OUT) :: ilow_pass_xso            !< type of oro eXtra SmOothing for steep
    !  orography and stencil width (1,4,5,6,8)
    INTEGER(KIND=i4), INTENT(OUT) :: numfilt_xso              !< number of applications of the eXtra filter
    LOGICAL,          INTENT(OUT) :: lxso_first               !< eXtra SmOothing before or after oro
    !  smoothing? (TRUE/FALSE)
    REAL(KIND=wp),    INTENT(OUT) :: rxso_mask                !< mask for eXtra SmOothing (threshold value)

    !> local variables
    INTEGER           :: nuin     !< unit number
    INTEGER (KIND=i4) :: ierr     !< error flag

    !> variables for default values
    LOGICAL           :: lfilter_oro_d,     &
         &               lxso_first_d

    INTEGER(KIND=i4)  :: ilow_pass_oro_d,   & 
         &               numfilt_oro_d,     &
         &               ifill_valley_d,    &
         &               ilow_pass_xso_d,   &
         &               numfilt_xso_d

    REAL(KIND=wp)     :: eps_filter_d,      &
         &               rfill_valley_d,    &   
         &               rxso_mask_d     

    !> define the namelist group
    NAMELIST /orography_smoothing/ &
         lfilter_oro, ilow_pass_oro, numfilt_oro, eps_filter, ifill_valley, &
         rfill_valley, ilow_pass_xso, numfilt_xso, lxso_first, rxso_mask

    !> initialization
    ierr     = 0

    !> default values definition
    lfilter_oro_d   = .FALSE.
    lxso_first_d    = .FALSE.
    ilow_pass_oro_d =  0
    numfilt_oro_d   =  1
    ifill_valley_d  =  1
    ilow_pass_xso_d =  1  
    numfilt_xso_d   =  1
    eps_filter_d    = 10
    rfill_valley_d  =  0  
    rxso_mask_d     =  0

    !> default values attribution
    lfilter_oro   = lfilter_oro_d 
    lxso_first    = lxso_first_d
    ilow_pass_oro = ilow_pass_oro_d
    numfilt_oro   = numfilt_oro_d
    ifill_valley  = ifill_valley_d
    ilow_pass_xso = ilow_pass_xso_d 
    numfilt_xso   = numfilt_xso_d
    eps_filter    = eps_filter_d
    rfill_valley  = rfill_valley_d  
    rxso_mask     = rxso_mask_d

    !> read namelist  
    nuin = free_un()  ! function free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL abort_extpar('read_namelists_extpar_orosmooth: cannot open file')
    ENDIF
    READ(nuin, NML=orography_smoothing, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL abort_extpar('read_namelists_extpar_orosmooth: cannot read file')
    ENDIF
    CLOSE(nuin)

    !> check values for consistency
    IF ( lfilter_oro ) THEN
      SELECT CASE( ilow_pass_oro )
      CASE( 1 )
        ! Default Raymond filter
      CASE( 3:6, 8 )
        ! Low-pass filter with predefined filter weights
        IF ( (numfilt_oro < 1) .OR. (numfilt_oro > 10) ) THEN
          WRITE(logging%fileunit,*)' WARNING:  *** numfilt_oro has to be between  1 and 10 *** '
          WRITE(logging%fileunit,*)'          *** set numfilt_oro = 1 (default value)!    *** '
          numfilt_oro = 1
        ENDIF
      CASE default
        WRITE(logging%fileunit,*)' WARNING:  *** ilow_pass_oro has to be 1, 3, 4, 5, 6 or 8 *** '
        WRITE(logging%fileunit,*)'          *** set ilow_pass_oro = 1 (default value)!   *** '
        ilow_pass_oro = 1
      END SELECT
      IF ( ilow_pass_xso >= ilow_pass_oro ) THEN
        SELECT CASE( ilow_pass_xso )
        CASE( 3:6, 8 )
          ! Low-pass filter with predefined filter weights
          IF ( (numfilt_xso < 1) .OR. (numfilt_xso > 10) ) THEN
            WRITE(logging%fileunit,*)' WARNING:  *** numfilt_xso has to be between  1 and 10 *** '
            WRITE(logging%fileunit,*)'          *** set numfilt_xso = 1 (default value)!    *** '
            numfilt_xso = 1
          ENDIF
        CASE default
          WRITE(logging%fileunit,*)' WARNING:  *** ilow_pass_xso has to be 3, 4, 5, 6 or 8 *** '
          WRITE(logging%fileunit,*)'          *** and to be  greater/equal ilow_pass_oro *** '
          WRITE(logging%fileunit,*)'          *** set ilow_pass_xso = 0 (default value)! *** '
          ilow_pass_xso = 0
        END SELECT
      ELSE
        IF ( ilow_pass_xso /= 0 ) THEN
          WRITE(logging%fileunit,*)' WARNING:  *** ilow_pass_xso has to be 3, 4, 5, 6 or 8 *** '
          WRITE(logging%fileunit,*)'          *** and to be  greater/equal ilow_pass_oro *** '
          WRITE(logging%fileunit,*)'          *** set ilow_pass_xso = 0 (default value)! *** '
          ilow_pass_xso = 0
        ENDIF
      ENDIF
      IF ((ifill_valley < 1) .OR. (ifill_valley > 2)) THEN
        WRITE(logging%fileunit,*)' WARNING:  *** ifill valley has to be 1 or 2 *** '
        WRITE(logging%fileunit,*)'          *** set ifill valley = 1 (default value)! *** '
        ifill_valley = 1
      ENDIF
    ENDIF

  END SUBROUTINE read_namelists_extpar_orosmooth
  !---------------------------------------------------------------------------
  !> subroutine to handle orography smoothing
  SUBROUTINE do_orosmooth   (tg,                                 &
       &                                      hh_target,        &
       &                                      fr_land_topo,    &
       &                                      lfilter_oro,      &
       &                                      ilow_pass_oro,    &
       &                                      numfilt_oro,      &
       &                                      eps_filter,       &
       &                                      ifill_valley,     &
       &                                      rfill_valley,     &
       &                                      ilow_pass_xso,    &
       &                                      numfilt_xso,      &
       &                                      lxso_first,       &
       &                                      rxso_mask ,       &
       &                                      hsmooth            ) 

    USE mo_grid_structures, ONLY: target_grid_def  !< Definition of data type with target grid definition

    TYPE(target_grid_def), INTENT(IN)      :: tg              !< !< structure with target grid description
    REAL(KIND=wp), INTENT(IN)   :: hh_target(1:tg%ie,1:tg%je,1:tg%ke)  !< mean height of target grid element
    REAL(KIND=wp), INTENT(IN)   :: fr_land_topo(1:tg%ie,1:tg%je,1:tg%ke)  !< mean height of target grid element
    LOGICAL, INTENT(IN) :: lfilter_oro  !< oro smoothing to be performed? (TRUE/FALSE) 
    INTEGER(KIND=i4), INTENT(IN) :: ilow_pass_oro            !< type of oro smoothing and 
    !  stencil width (1,4,5,6,8)
    INTEGER(KIND=i4), INTENT(IN) :: numfilt_oro              !< number of applications of the filter
    REAL(KIND=wp),    INTENT(IN) :: eps_filter               !< smoothing param ("strength" of the filtering)
    INTEGER(KIND=i4), INTENT(IN) :: ifill_valley             !< fill valleys before or after oro smoothing 
    !  (1: before, 2: after)
    REAL(KIND=wp),    INTENT(IN) :: rfill_valley             !< mask for valley filling (threshold value)
    INTEGER(KIND=i4), INTENT(IN) :: ilow_pass_xso            !< type of oro eXtra SmOothing for steep
    !  orography and stencil width (1,4,5,6,8)
    INTEGER(KIND=i4), INTENT(IN) :: numfilt_xso              !< number of applications of the eXtra filter
    LOGICAL,          INTENT(IN) :: lxso_first               !< eXtra SmOothing before or after oro
    !  smoothing? (TRUE/FALSE)
    REAL(KIND=wp),    INTENT(IN) :: rxso_mask                !< mask for eXtra SmOothing (threshold value)
    REAL(KIND=wp), INTENT(INOUT) :: hsmooth(1:tg%ie,1:tg%je,1:tg%ke)   !< mean smoothed height of target grid element

    !> local variables
    INTEGER(KIND=i8)             :: ndim, ie, je, ie_ext_hf, je_ext_hf, &
         ile, iri, jlo, jup
    INTEGER(KIND=i4)             :: hfwidth, hfw_m_nb, errorcode,       &
         n, nit, nfrl
    REAL(KIND=wp)                :: zdh_x1, zdh_x2, zdh_y1, zdh_y2,&
         zdh_max, zdh_xy1, zdh_xy2,     &
         zdh_xy3, zdh_xy4 
    REAL(KIND=wp), ALLOCATABLE   :: xy_vec(:), ci(:), cj(:), ck(:),&
         cl(:), cm(:), a(:), b(:),      &
         c(:), d(:), e(:), f(:),        &
         ff_filt(:,:,:)
    LOGICAL, ALLOCATABLE         :: hfx_mask(:,:),  &
         hfy_mask(:,:)
    LOGICAL                      :: ldhsurf_xy

    REAL(KIND=wp),PARAMETER      :: zhmax_sea = 1.0


    ALLOCATE( ff_filt(tg%ie,tg%je,tg%ke), STAT = errorcode )
    IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate ff_filt')

    ff_filt(:,:,:) = hh_target (:,:,:)


    !> filling of valleys before orography smoothing
    IF ( ifill_valley == 1 ) CALL fill_valleys (rfill_valley, tg,         &
         ifill_valley, ff_filt) 

    !> orography smoothing
    SELECT CASE( ilow_pass_oro )

      !>> global Raymond filtering
    CASE( 1 )

      DO nit = 1, numfilt_oro
        !> standard orography smoothing performed
        !  before eXtra SmOothing of steep oro.
        IF ( .NOT.lxso_first ) THEN

          !> filtering in x-direction

          ! Set the dimension
          ndim = tg%ie   
          ! allocate the necessary fields for gaussian elimination
          ALLOCATE (xy_vec(ndim), ci(ndim), cj(ndim), ck(ndim), cl(ndim), cm(ndim),&
               a(ndim),  b(ndim),  c(ndim),  d(ndim),  e(ndim),  f(ndim),STAT=errorcode)
          IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate xy_vec, ci, and so on')
          ! Set the above variables for the gaussian elimination
          CALL set_gauss (a, b, c, d, e, f, ci, cj, ck, cl, cm, ndim, eps_filter)
          ! Apply the filter for every line
          DO je = 1, tg%je
            xy_vec(:) = ff_filt(:,je,1)
            CALL low_pass_filter (xy_vec, a, b, c, d, e, f, ci, cj, ck, cl, cm,   &
                 ndim, eps_filter)
            ff_filt(:,je,1) = xy_vec(:)
          ENDDO
          ! Release memory
          DEALLOCATE (xy_vec, ci, cj, ck, cl, cm, a, b, c, d, e, f,STAT=errorcode)
          IF(errorcode/=0) CALL abort_extpar('Cant deallocate xy_vec, ci, and so on')

          !> filtering in y-direction

          ! Set the dimension
          ndim = tg%je
          ! allocate the necessary fields for gaussian elimination
          ALLOCATE (xy_vec(ndim), ci(ndim), cj(ndim), ck(ndim), cl(ndim), cm(ndim),&
               a(ndim),  b(ndim),  c(ndim),  d(ndim),  e(ndim),  f(ndim),STAT=errorcode) 
          IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate xy_vec, ci, and so on')
          ! Set the above variables for the gaussian elimination
          CALL set_gauss (a, b, c, d, e, f, ci, cj, ck, cl, cm, ndim, eps_filter)
          ! Apply the filter for every column
          DO ie = 1, tg%ie
            xy_vec(:) = ff_filt(ie,:,1)
            CALL low_pass_filter (xy_vec, a, b, c, d, e, f, ci, cj, ck, cl, cm,   &
                 ndim, eps_filter)
            ff_filt(ie,:,1) = xy_vec(:)
          ENDDO
          ! Release memory
          DEALLOCATE (xy_vec, ci, cj, ck, cl, cm, a, b, c, d, e, f,STAT=errorcode)
          IF(errorcode/=0) CALL abort_extpar('Cant deallocate xy_vec, ci, and so on')

        ENDIF ! .NOT.lxso_first

        !> eXtra SmOothing of steep oro. 
        ldhsurf_xy = .TRUE. !steepness of oro. in diago. dir. taken into account
        IF ( ilow_pass_xso >= ilow_pass_oro ) THEN
          IF ( rxso_mask > 0.0_wp ) THEN           
            ! set width of the stencil for the horizontal filter
            SELECT CASE( ilow_pass_xso )
            CASE( 3, 4, 6 )
              hfwidth = 4
            CASE( 5, 8 )
              hfwidth = 6
            END SELECT
            hfw_m_nb = hfwidth 
            ie_ext_hf = tg%ie + 2*hfw_m_nb
            je_ext_hf = tg%je + 2*hfw_m_nb  
            ALLOCATE( hfx_mask(ie_ext_hf,je_ext_hf),  &
                 hfy_mask(ie_ext_hf,je_ext_hf),  &
                 STAT = errorcode )
            IF(errorcode/=0) CALL abort_extpar('Cant allocate hfx_mask, hfy_mask')

            DO n = 1, numfilt_xso        
              hfx_mask(:,:) = .FALSE.
              hfy_mask(:,:) = .FALSE.
              ! set mask for extra smoothing
              DO je = 2, tg%je-1
                DO ie = 2, tg%ie-1         
                  zdh_x1 = ABS( ff_filt(ie-1,je,1) - ff_filt(ie,je,1) )
                  zdh_x2 = ABS( ff_filt(ie+1,je,1) - ff_filt(ie,je,1) )
                  zdh_max = MAX( zdh_x1, zdh_x2 )
                  IF ( zdh_max > rxso_mask ) THEN
                    hfx_mask(ie+hfw_m_nb,je+hfw_m_nb) = .TRUE.
                  ENDIF
                  zdh_y1 = ABS( ff_filt(ie,je-1,1) - ff_filt(ie,je,1) )
                  zdh_y2 = ABS( ff_filt(ie,je+1,1) - ff_filt(ie,je,1) )
                  zdh_max = MAX( zdh_y1, zdh_y2 )
                  IF ( zdh_max > rxso_mask ) THEN
                    hfy_mask(ie+hfw_m_nb,je+hfw_m_nb) = .TRUE.
                  ENDIF
                  IF ( ldhsurf_xy ) THEN
                    zdh_xy1 = ABS( ff_filt(ie-1,je-1,1) - ff_filt(ie,je,1) )
                    zdh_xy2 = ABS( ff_filt(ie+1,je+1,1) - ff_filt(ie,je,1) )
                    zdh_xy3 = ABS( ff_filt(ie-1,je+1,1) - ff_filt(ie,je,1) )
                    zdh_xy4 = ABS( ff_filt(ie+1,je-1,1) - ff_filt(ie,je,1) )
                    zdh_max = MAX( zdh_xy1, zdh_xy2, zdh_xy3, zdh_xy4 )
                    IF ( zdh_max > SQRT(2.0_wp)*rxso_mask ) THEN
                      hfx_mask(ie+hfw_m_nb,je+hfw_m_nb) = .TRUE.
                      hfy_mask(ie+hfw_m_nb,je+hfw_m_nb) = .TRUE.
                    ENDIF
                  ENDIF

                ENDDO
              ENDDO
              CALL hfilter_orography( ncutoff=ilow_pass_xso, lhf_mask=.TRUE.,&
                   hfx_mask=hfx_mask, hfy_mask=hfy_mask,  &
                   tg=tg, ie_ext_hf=ie_ext_hf,            &
                   je_ext_hf=je_ext_hf, hfw_m_nb=hfw_m_nb,&
                   hfwidth=hfwidth, field=ff_filt)  
            ENDDO

            DEALLOCATE( hfx_mask, hfy_mask, STAT = errorcode )
            IF(errorcode/=0) CALL abort_extpar('Cant deallocate hfx_mask, hfy_mask')

          ENDIF ! ilow_pass_xso >= ilow_pass_oro
        ENDIF ! rxso_mask > 0.0_wp

        IF ( lxso_first ) THEN

          !> filtering in x-direction

          ! Set the dimension
          ndim = tg%ie   
          ! allocate the necessary fields for gaussian elimination
          ALLOCATE (xy_vec(ndim), ci(ndim), cj(ndim), ck(ndim), cl(ndim), cm(ndim),&
               a(ndim),  b(ndim),  c(ndim),  d(ndim),  e(ndim),  f(ndim),STAT=errorcode)
          IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate xy_vec, ci, and so on')
          ! Set the above variables for the gaussian elimination
          CALL set_gauss (a, b, c, d, e, f, ci, cj, ck, cl, cm, ndim, eps_filter)
          ! Apply the filter for every line
          DO je = 1, tg%je
            xy_vec(:) = ff_filt(:,je,1)
            CALL low_pass_filter (xy_vec, a, b, c, d, e, f, ci, cj, ck, cl, cm,   &
                 ndim, eps_filter)
            ff_filt(:,je,1) = xy_vec(:)
          ENDDO
          ! Release memory
          DEALLOCATE (xy_vec, ci, cj, ck, cl, cm, a, b, c, d, e, f,STAT=errorcode)
          IF(errorcode/=0) CALL abort_extpar('Cant deallocate xy_vec, ci, and so on')

          !> filtering in y-direction

          ! Set the dimension
          ndim = tg%je
          ! allocate the necessary fields for gaussian elimination
          ALLOCATE (xy_vec(ndim), ci(ndim), cj(ndim), ck(ndim), cl(ndim), cm(ndim),&
               a(ndim),  b(ndim),  c(ndim),  d(ndim),  e(ndim),  f(ndim),STAT=errorcode) 
          IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate xy_vec, ci, and so on')
          ! Set the above variables for the gaussian elimination
          CALL set_gauss (a, b, c, d, e, f, ci, cj, ck, cl, cm, ndim, eps_filter)
          ! Apply the filter for every column
          DO ie = 1, tg%ie
            xy_vec(:) = ff_filt(ie,:,1)
            CALL low_pass_filter (xy_vec, a, b, c, d, e, f, ci, cj, ck, cl, cm,   &
                 ndim, eps_filter)
            ff_filt(ie,:,1) = xy_vec(:)
          ENDDO
          ! Release memory
          DEALLOCATE (xy_vec, ci, cj, ck, cl, cm, a, b, c, d, e, f,STAT=errorcode)
          IF(errorcode/=0) CALL abort_extpar('Cant deallocate xy_vec, ci, and so on')

        ENDIF ! lxso_first
      ENDDO ! numfilt_oro

      !>> local filtering (weights by J. Foerstner, DWD)
    CASE( 3:6, 8 )
      !> standard orography smoothing performed
      !  before eXtra SmOothing of steep oro.
      IF ( .NOT.lxso_first ) THEN
        ! set width of the stencil for the horizontal filter
        SELECT CASE( ilow_pass_oro )
        CASE( 3, 4, 6 )
          hfwidth = 4
        CASE( 5, 8 )
          hfwidth = 6
        END SELECT
        hfw_m_nb = hfwidth 
        ie_ext_hf = INT(tg%ie + 2*hfw_m_nb  ,i8)
        je_ext_hf = INT(tg%je + 2*hfw_m_nb ,i8)
        DO n = 1, numfilt_oro
          CALL hfilter_orography( ncutoff=ilow_pass_oro, lhf_mask=.FALSE., &
               tg=tg, ie_ext_hf=ie_ext_hf,              &
               je_ext_hf=je_ext_hf,hfw_m_nb=hfw_m_nb,   &
               hfwidth=hfwidth, field=ff_filt)
        ENDDO
      ENDIF

      !> eXtra SmOothing of steep oro. 
      ldhsurf_xy = .TRUE. !steepness of oro. in diago. dir. taken into account
      IF ( ilow_pass_xso >= ilow_pass_oro ) THEN
        IF ( rxso_mask > 0.0_wp ) THEN           
          ! set width of the stencil for the horizontal filter
          SELECT CASE( ilow_pass_xso )
          CASE( 3, 4, 6 )
            hfwidth = 4
          CASE( 5, 8 )
            hfwidth = 6
          END SELECT
          hfw_m_nb = hfwidth 
          ie_ext_hf = tg%ie + 2*hfw_m_nb
          je_ext_hf = tg%je + 2*hfw_m_nb  
          ALLOCATE( hfx_mask(ie_ext_hf,je_ext_hf),  &
               hfy_mask(ie_ext_hf,je_ext_hf),  &
               STAT = errorcode )
          IF(errorcode/=0) CALL abort_extpar('Cant allocate hfx_mask, hfy_mask')

          DO n = 1, numfilt_xso        
            hfx_mask(:,:) = .FALSE.
            hfy_mask(:,:) = .FALSE.
            ! set mask for extra smoothing
            DO je = 2, tg%je-1
              DO ie = 2, tg%ie-1         
                zdh_x1 = ABS( ff_filt(ie-1,je,1) - ff_filt(ie,je,1) )
                zdh_x2 = ABS( ff_filt(ie+1,je,1) - ff_filt(ie,je,1) )
                zdh_max = MAX( zdh_x1, zdh_x2 )
                IF ( zdh_max > rxso_mask ) THEN
                  hfx_mask(ie+hfw_m_nb,je+hfw_m_nb) = .TRUE.
                ENDIF
                zdh_y1 = ABS( ff_filt(ie,je-1,1) - ff_filt(ie,je,1) )
                zdh_y2 = ABS( ff_filt(ie,je+1,1) - ff_filt(ie,je,1) )
                zdh_max = MAX( zdh_y1, zdh_y2 )
                IF ( zdh_max > rxso_mask ) THEN
                  hfy_mask(ie+hfw_m_nb,je+hfw_m_nb) = .TRUE.
                ENDIF
                IF ( ldhsurf_xy ) THEN
                  zdh_xy1 = ABS( ff_filt(ie-1,je-1,1) - ff_filt(ie,je,1) )
                  zdh_xy2 = ABS( ff_filt(ie+1,je+1,1) - ff_filt(ie,je,1) )
                  zdh_xy3 = ABS( ff_filt(ie-1,je+1,1) - ff_filt(ie,je,1) )
                  zdh_xy4 = ABS( ff_filt(ie+1,je-1,1) - ff_filt(ie,je,1) )
                  zdh_max = MAX( zdh_xy1, zdh_xy2, zdh_xy3, zdh_xy4 )
                  IF ( zdh_max > SQRT(2.0_wp)*rxso_mask ) THEN
                    hfx_mask(ie+hfw_m_nb,je+hfw_m_nb) = .TRUE.
                    hfy_mask(ie+hfw_m_nb,je+hfw_m_nb) = .TRUE.
                  ENDIF
                ENDIF

              ENDDO
            ENDDO
            CALL hfilter_orography( ncutoff=ilow_pass_xso, lhf_mask=.TRUE.,&
                 hfx_mask=hfx_mask, hfy_mask=hfy_mask,  &
                 tg=tg, ie_ext_hf=ie_ext_hf,            &
                 je_ext_hf=je_ext_hf, hfw_m_nb=hfw_m_nb,&
                 hfwidth=hfwidth, field=ff_filt)  
          ENDDO

          DEALLOCATE( hfx_mask, hfy_mask, STAT = errorcode )
          IF(errorcode/=0) CALL abort_extpar('Cant deallocate hfx_mask, hfy_mask')


        ENDIF
      ENDIF
      !> standard orography smoothing performed
      !  after eXtra SmOothing of steep oro.
      IF ( lxso_first ) THEN
        ! set width of the stencil for the horizontal filter
        SELECT CASE( ilow_pass_oro )
        CASE( 3, 4, 6 )
          hfwidth = 4
        CASE( 5, 8 )
          hfwidth = 6
        END SELECT
        hfw_m_nb = hfwidth 
        ie_ext_hf = tg%ie + 2*hfw_m_nb
        je_ext_hf = tg%je + 2*hfw_m_nb
        DO n = 1, numfilt_oro
          CALL hfilter_orography( ncutoff=ilow_pass_oro, lhf_mask=.FALSE., &
               tg=tg, ie_ext_hf=ie_ext_hf,              &
               je_ext_hf=je_ext_hf,hfw_m_nb=hfw_m_nb,   &
               hfwidth=hfwidth, field=ff_filt)
        ENDDO
      ENDIF

    END SELECT

    !> correction of the filtered oro. by setting the height to the original
    !  value if:
    !     - it is a sea point surrounded by at least 4 other sea-points
    !     - the sign of the filtered oro. does not coincide with the
    !       original sign
    DO je = 1, tg%je
      DO ie = 1, tg%ie
        IF ((fr_land_topo(ie,je,1) < 0.5) .AND.  &
             (hh_target(ie,je,1) <= zhmax_sea)) THEN
          ile = MAX (1_i8,ie-1_i8)
          iri = MIN (ie+1_i8,tg%ie)
          jlo = MAX (1_i8,je-1_i8)
          jup = MIN (je+1_i8,tg%je)
          nfrl = 0
          IF (fr_land_topo(ile,jlo,1 ) < 0.5) nfrl = nfrl + 1
          IF (fr_land_topo(ie ,jlo,1 ) < 0.5) nfrl = nfrl + 1
          IF (fr_land_topo(iri,jlo,1 ) < 0.5) nfrl = nfrl + 1
          IF (fr_land_topo(ile,je ,1 ) < 0.5) nfrl = nfrl + 1
          IF (fr_land_topo(iri,je ,1 ) < 0.5) nfrl = nfrl + 1
          IF (fr_land_topo(ile,jup,1 ) < 0.5) nfrl = nfrl + 1
          IF (fr_land_topo(ie ,jup,1 ) < 0.5) nfrl = nfrl + 1
          IF (fr_land_topo(iri,jup,1 ) < 0.5) nfrl = nfrl + 1
          IF (nfrl >= 4) THEN
            ff_filt(ie,je,1) = hh_target(ie,je,1)
          ENDIF
        ENDIF

        IF (ff_filt(ie,je,1)*hh_target(ie,je,1) < 0.0) THEN
          ff_filt(ie,je,1) = hh_target(ie,je,1)
        ENDIF
      ENDDO
    ENDDO

    !> filling of valleys after orography smoothing
    IF ( ifill_valley == 2 ) CALL fill_valleys (rfill_valley, tg,         &
         ifill_valley, ff_filt)

    !> save back the filtered oro. in the variable hsmooth 
    !  and deallocate the tmp variable ff_filt

    hsmooth(:,:,1) = ff_filt (:,:,1)

    DEALLOCATE( ff_filt, STAT = errorcode )
    IF (errorcode /= 0 ) CALL abort_extpar('Cant deallocate ff_filt')


  END SUBROUTINE do_orosmooth
  !---------------------------------------------------------------------------

  SUBROUTINE fill_valleys (rfill_valley,                                    &
       &                  tg,                                               &
       &                  ifill_valley,                                   &
       &                  hh_filt)

    USE mo_grid_structures, ONLY: target_grid_def  !< Definition of data type with target grid definition


    TYPE(target_grid_def), INTENT(IN)      :: tg              !< !< structure with target grid description

    INTEGER(KIND=i4) , INTENT(IN):: ifill_valley             !< fill valleys before or after oro smoothing 
    !  (1: before, 2: after)
    REAL(KIND=wp) , INTENT(IN) :: rfill_valley             !< mask for valley filling (threshold value)


    REAL(KIND=wp), INTENT(INOUT) :: hh_filt(1:tg%ie,1:tg%je,1:tg%ke)   !< mean smoothed height of target grid element



    ! local variables
    INTEGER :: errorcode !< error status variable
    INTEGER (KIND=i8) :: ie, je, ke  ! indices for grid elements

    REAL(KIND=wp) :: zdh_x1, zdh_x2, zdh_y1, zdh_y2, zdh_max               

    REAL (KIND=wp), ALLOCATABLE ::  &
         ff_tmp(:,:,:)



    ! Check if dh <= rfill_valley in a V-valley
    ! otherwise fill valley with maximum value of neighborhood
    IF ( rfill_valley > 0.0_wp ) THEN

      ! Allocate and set temporary field
      ALLOCATE( ff_tmp(tg%ie,tg%je,tg%ke), STAT = errorcode )
      IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate ff_tmp')

      ff_tmp(:,:,:) = hh_filt(:,:,:)


      DO ke =1, tg%ke
        DO je = 2, tg%je-1
          DO ie = 2, tg%ie-1
            zdh_x1 = ff_tmp(ie-1,je,ke) - ff_tmp(ie,je,ke)
            zdh_x2 = ff_tmp(ie+1,je,ke) - ff_tmp(ie,je,ke)
            zdh_y1 = ff_tmp(ie,je-1,ke) - ff_tmp(ie,je,ke)
            zdh_y2 = ff_tmp(ie,je+1,ke) - ff_tmp(ie,je,ke)

            IF ( zdh_x1 > 0.0_wp .AND. zdh_x2 > 0.0_wp .AND.      &
                 zdh_y1 > 0.0_wp .AND. zdh_y2 > 0.0_wp ) THEN
              zdh_max = MAX( zdh_x1, zdh_x2, zdh_y1, zdh_y2 )

              IF ( zdh_max > rfill_valley ) THEN
                ! 1-3: BEFORE  / 4-6: AFTER filtering of orography
                SELECT CASE( ifill_valley )
                CASE( 1, 2 )
                  ! MIN value      of 4-point neighborhood
                  hh_filt(ie,je,ke) = MIN( ff_tmp(ie-1,je,ke), ff_tmp(ie+1,je,ke),    &
                       ff_tmp(ie,je-1,ke), ff_tmp(ie,je+1,ke) )
                END SELECT
              ENDIF

            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE( ff_tmp, STAT = errorcode )
      IF(errorcode/=0) CALL abort_extpar('Cant deallocate ff_tmp')

    ENDIF

  END SUBROUTINE fill_valleys
  !------------------------------------------------------------------------------

  SUBROUTINE set_gauss (a, b, c, d, e, f, ci, cj, ck, cl, cm, ndim, eps)

    !------------------------------------------------------------------------------
    !
    ! Description:
    !
    ! Method:
    !
    !------------------------------------------------------------------------------

    ! Parameterlist

    INTEGER (KIND=i8), INTENT(IN)    ::   &
         ndim         ! dimension of xy_vec

    REAL    (KIND=wp)   , INTENT(IN)    ::   &
         eps                   ! filter parameter

    REAL    (KIND=wp)   , INTENT(OUT)   ::   &
         a(ndim),  b(ndim),  c(ndim),  d(ndim),  e(ndim),  f(ndim),    &
         ci(ndim), cj(ndim), ck(ndim), cl(ndim), cm(ndim)

    ! Local variables
    REAL    (KIND=wp)                   ::   &
         o(ndim),  p(ndim),  q(ndim),  r(ndim),  s(ndim),  t(ndim),    &
         u(ndim),  v(ndim),  w(ndim),  x(ndim)

    INTEGER (KIND=i8)                :: n

    !------------------------------------------------------------------------------
    !- End of header -
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE set_gauss
    !------------------------------------------------------------------------------

    ! Setup up the gaussian elimination step
    ! coefficients of matrix

    a(:) = 0.0;  b(:) = 0.0;  c(:) = 0.0;  d(:) = 0.0;  e(:) = 0.0; f(:) = 0.0;
    ci(:) = 0.0; cj(:) = 0.0; ck(:) = 0.0; cl(:) = 0.0; cm(:) = 0.0;
    o(:) = 0.0;  p(:) = 0.0;  q(:) = 0.0;  r(:) = 0.0;  s(:) = 0.0;  
    t(:) = 0.0;  u(:) = 0.0;  v(:) = 0.0;  w(:) = 0.0;  x(:) = 0.0;

    s(     1) = 1.0
    s(ndim  ) = s(1)
    s(     2) = 2.0*(1.0+eps)
    s(ndim-1) = s(2)
    s(     3) = 6.0*(1.0+eps)
    s(ndim-2) = s(3)
    s(     4) =20.0*(1.0+eps)
    s(ndim-3) = s(4)
    s(     5) =70.0*(1.0+eps)
    s(ndim-4) = s(5)

    t(     2) = 1.0-eps
    t(ndim-1) = t(2)
    t(     3) = 4.0*(1.0-eps)
    t(ndim-2) = t(3)
    t(     4) =15.0*(1.0-eps)
    t(ndim-3) = t(4)
    t(     5) =56.0*(1.0-eps)
    t(ndim-4) = t(5)

    u(     3) = 1.0+eps
    u(ndim-2) = u(3)
    u(     4) = 6.0*(1.0+eps)
    u(ndim-3) = u(4)
    u(     5) =28.0*(1.0+eps)
    u(ndim-4) = u(5)

    v(     4) = 1.0-eps
    v(ndim-3) = v(4)
    v(     5) = 8.0*(1.0-eps)
    v(ndim-4) = v(5)

    w(     5) = 1.0+eps
    w(ndim-4) = w(5)

    DO n = 6, ndim-5
      x(n) =  1.0-eps
      w(n) = 10.0*(1.0+eps)
      v(n) = 45.0*(1.0-eps)
      u(n) =120.0*(1.0+eps)
      t(n) =210.0*(1.0-eps)
      s(n) =252.0*(1.0+eps)
    ENDDO

    r(:) = t(:)
    q(:) = u(:)
    p(:) = v(:)
    o(:) = w(:)
    a(:) = x(:)

    !-------------Matrix inversion--------------

    ! Step I

    f (1) = s(1)
    ci(1) = t(1)/f(1)
    e (2) = r(2)
    f (2) = s(2)-e(2)*ci(1)
    cj(1) = u(1)/f(1)
    ci(2) = (t(2)-e(2)*cj(1))/f(2)
    d (3) = q(3)
    e (3) = r(3)-d(3)*ci(1)
    f (3) = s(3)-e(3)*ci(2)-d(3)*cj(1)
    ck(1) = v(1)/f(1)
    cj(2) = (u(2)-e(2)*ck(1))/f(2)
    ci(3) = (t(3)-d(3)*ck(1)-e(3)*cj(2))/f(3)
    c (4) = p(4)
    d (4) = q(4)-c(4)*ci(1)
    e (4) = r(4)-d(4)*ci(2)-c(4)*cj(1)
    f (4) = s(4)-e(4)*ci(3)-d(4)*cj(2)-c(4)*ck(1)
    cl(1) = w(1)/f(1)
    ck(2) = (v(2)-e(2)*cl(1))/f(2)
    cj(3) = (u(3)-d(3)*cl(1)-e(3)*ck(2))/f(3)
    ci(4) = (t(4)-c(4)*cl(1)-d(4)*ck(2)-e(4)*cj(3))/f(4)
    b (5) = o(5)
    c (5) = p(5)-b(5)*ci(1)
    d (5) = q(5)-c(5)*ci(2)-b(5)*cj(1)
    e (5) = r(5)-d(5)*ci(3)-c(5)*cj(2)-b(5)*ck(1)
    f (5) = s(5)-e(5)*ci(4)-d(5)*cj(3)-c(5)*ck(2)-b(5)*cl(1) 

    ! Step II
    DO n = 6, ndim
      cm(n-5) = x(n-5)/f(n-5)
      cl(n-4) = (w(n-4)-e(n-4)*cm(n-5))/f(n-4)
      ck(n-3) = (v(n-3)-d(n-3)*cm(n-5)-e(n-3)*cl(n-4))/f(n-3)
      cj(n-2) = (u(n-2)-c(n-2)*cm(n-5)-d(n-2)*cl(n-4)-e(n-2)*ck(n-3))/f(n-2)
      ci(n-1) = (t(n-1)-b(n-1)*cm(n-5)-c(n-1)*cl(n-4)-d(n-1)*ck(n-3)-  &
           &     e(n-1)*cj(n-2))/f(n-1)
      b (n  ) = o(n)-a(n)*ci(n-5)
      c (n  ) = p(n)-b(n)*ci(n-4)-a(n)*cj(n-5)
      d (n  ) = q(n)-c(n)*ci(n-3)-b(n)*cj(n-4)-a(n)*ck(n-5)
      e (n  ) = r(n)-d(n)*ci(n-2)-c(n)*cj(n-3)-b(n)*ck(n-4)-a(n)*cl(n-5)
      f (n  ) = s(n)-e(n)*ci(n-1)-d(n)*cj(n-2)-c(n)*ck(n-3)-b(n)*cl(n-4)- &
           &    a(n)*cm(n-5)
    ENDDO

  END SUBROUTINE set_gauss

  !------------------------------------------------------------------------------

  SUBROUTINE low_pass_filter (xy_vec, a, b, c, d, e, f, ci, cj, ck, cl, cm,   &
       &                      ndim, eps)

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   The filter of Raymond (Raymond, W.H., 1988: High-Order Low-Pass Implicit
    !   Tangent Filters for Use in Finite Area Calculations, MWR 116, 2132-2141)
    !   is used for filtering orography. We use fifth order filter. Computations
    !   are as proposed in the appendix of the paper. Some corrections to the
    !   filtered field are done afterwards.
    !
    ! Method:
    !
    !------------------------------------------------------------------------------

    ! Parameterlist

    INTEGER (KIND=i8), INTENT(IN)    :: ndim         ! dimension of xy_vec

    REAL    (KIND=wp), INTENT(INOUT) :: xy_vec(ndim) ! one dimensional vector to be filtered

    REAL    (KIND=wp), INTENT(IN)    :: eps,  &      ! filter parameter
                                                     ! variables for gaussian elimination:
         &                              a(ndim),  b(ndim),  c(ndim),  d(ndim),  e(ndim),  f(ndim), &
         &                              ci(ndim), cj(ndim), ck(ndim), cl(ndim), cm(ndim)

    ! Local variables
    REAL    (KIND=wp)                :: rhs(ndim), & ! RHS of equation
         &                              phi(ndim), & ! filter increment
         &                              h(ndim)      ! for intermediate storage

    INTEGER (KIND=i8)                :: n

    !> Specify RHS of equation

    rhs(     1) = 0
    rhs(ndim  ) = 0

    rhs(     2) = eps*(xy_vec(     1)-2.0*xy_vec(     2)+xy_vec(     3))
    rhs(ndim-1) = eps*(xy_vec(ndim-2)-2.0*xy_vec(ndim-1)+xy_vec(ndim  )) 

    rhs(     3) = eps*(-1.0*(xy_vec(     1)+xy_vec(     5))     & 
         &             +4.0*(xy_vec(     2)+xy_vec(     4))     &
         &             -6.0* xy_vec(     3)               )

    rhs(ndim-2) = eps*(-1.0*(xy_vec(ndim  )+xy_vec(ndim-4))     &
         &             +4.0*(xy_vec(ndim-1)+xy_vec(ndim-3))     &
         &             -6.0* xy_vec(ndim-2)         )

    rhs(     4) = eps*(      xy_vec(   1)+xy_vec(   7)          &
         &             -6.0*(xy_vec(   2)+xy_vec(   6))         &
         &            +15.0*(xy_vec(   3)+xy_vec(   5))         &
         &            -20.0* xy_vec(   4)             )

    rhs(ndim-3) = eps*(      xy_vec(ndim-6)+xy_vec(  ndim)      &
         &             -6.0*(xy_vec(ndim-5)+xy_vec(ndim-1))     &
         &            +15.0*(xy_vec(ndim-4)+xy_vec(ndim-2))     &
         &            -20.0* xy_vec(ndim-3)               )

    rhs(     5) = eps*(-1.0*(xy_vec(   1)+xy_vec(   9))         &
         &             +8.0*(xy_vec(   2)+xy_vec(   8))         &
         &            -28.0*(xy_vec(   3)+xy_vec(   7))         &
         &            +56.0*(xy_vec(   4)+xy_vec(   6))         &
         &            -70.0* xy_vec(   5)             )

    rhs(ndim-4) = eps*(-1.0*(xy_vec(ndim-8)+xy_vec(  ndim))     &
         &             +8.0*(xy_vec(ndim-7)+xy_vec(ndim-1))     &
         &            -28.0*(xy_vec(ndim-6)+xy_vec(ndim-2))     &
         &            +56.0*(xy_vec(ndim-5)+xy_vec(ndim-3))     &
         &            -70.0* xy_vec(ndim-4)               )

    DO n = 6, ndim-5
      rhs(n) = eps*(         xy_vec( n-5)+xy_vec( n+5)          &
         &            -10.0*(xy_vec( n-4)+xy_vec( n+4))         &
         &            +45.0*(xy_vec( n-3)+xy_vec( n+3))         &
         &           -120.0*(xy_vec( n-2)+xy_vec( n+2))         &
         &           +210.0*(xy_vec( n-1)+xy_vec( n+1))         &
         &           -252.0* xy_vec(   n)             )
    ENDDO

    !> Compute filter increment and filtered variables

    ! Step III

    h(1) =  rhs(1)/f(1)
    h(2) = (rhs(2)-h(1)*e(2))/f(2)
    h(3) = (rhs(3)-h(1)*d(3)-h(2)*e(3))/f(3)
    h(4) = (rhs(4)-h(1)*c(4)-h(2)*d(4)-h(3)*e(4))/f(4)
    h(5) = (rhs(5)-h(1)*b(5)-h(2)*c(5)-h(3)*d(5)-h(4)*e(5))/f(5)

    DO n = 6, ndim
      h(n) = (rhs(n)-h(n-5)*a(n)-h(n-4)*b(n)-h(n-3)*c(n)-h(n-2)*d(n)- &
           h(n-1)*e(n))/f(n)
    ENDDO

    ! Step IV
    ! Determine filter increment phi 
    phi(ndim  ) = h(ndim)
    phi(ndim-1) = h(ndim-1)-ci(ndim-1)*phi(ndim)
    phi(ndim-2) = h(ndim-2)-ci(ndim-2)*phi(ndim-1)-cj(ndim-2)*phi(ndim)
    phi(ndim-3) = h(ndim-3)-ci(ndim-3)*phi(ndim-2)-cj(ndim-3)*phi(ndim-1)-   &
         &                                         ck(ndim-3)*phi(ndim)
    phi(ndim-4) = h(ndim-4)-ci(ndim-4)*phi(ndim-3)-cj(ndim-4)*phi(ndim-2)-   &
         &                  ck(ndim-4)*phi(ndim-1)-cl(ndim-4)*phi(ndim)

    DO n = ndim-5,1,-1
      phi(n) = h(n)-ci(n)*phi(n+1)-cj(n)*phi(n+2)-ck(n)*phi(n+3)-      &
         &          cl(n)*phi(n+4)-cm(n)*phi(n+5)
    ENDDO

    DO n = 1, ndim
      xy_vec(n) = xy_vec(n) + phi(n)
    ENDDO

  END SUBROUTINE low_pass_filter

  !------------------------------------------------------------------------------

  SUBROUTINE hfilter_orography( ncutoff, lhf_mask,hfx_mask, hfy_mask, tg, ie_ext_hf,    &
       je_ext_hf, hfw_m_nb, hfwidth, field)

    USE mo_grid_structures, ONLY: target_grid_def  !< Definition of data type with target grid definition

    INTEGER (KIND=i4), INTENT(IN)     :: ncutoff
    LOGICAL, INTENT(IN)               :: lhf_mask
    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i8), INTENT(IN)     :: ie_ext_hf
    INTEGER (KIND=i8), INTENT(IN)     :: je_ext_hf
    INTEGER (KIND=i4), INTENT(IN)     :: hfw_m_nb
    INTEGER (KIND=i4), INTENT(IN)     :: hfwidth

    LOGICAL, INTENT(IN), OPTIONAL     :: hfx_mask(ie_ext_hf,je_ext_hf) 
    LOGICAL, INTENT(IN), OPTIONAL     ::  hfy_mask(ie_ext_hf,je_ext_hf)


    REAL(KIND=wp), INTENT(INOUT)      :: field(1:tg%ie,1:tg%je,1:tg%ke)   

    ! Local variables
    REAL (KIND=wp), ALLOCATABLE       :: ff_tmp(:,:)
    INTEGER (KIND=i4)                 :: errorcode
    INTEGER (KIND=i8)                 :: i, j

    ! Allocate and set temporary field
    ALLOCATE( ff_tmp(ie_ext_hf,je_ext_hf), STAT = errorcode )
    IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate ff_tmp')

    CALL extend_field2D ( field(:,:,1), tg%ie, tg%je,                &
         &                ff_tmp(:,:), ie_ext_hf, je_ext_hf,         &
         &                hfw_m_nb)

    IF ( lhf_mask ) THEN
      CALL horizontal_filtering( ff_tmp(:,:), ie_ext_hf, je_ext_hf,  &
           &                     hfwidth, ncutoff,         &
           &                     hfx_mask, hfy_mask )
    ELSE
      CALL horizontal_filtering( ff_tmp(:,:), ie_ext_hf, je_ext_hf,  &
           &                     hfwidth, ncutoff)
    ENDIF

    DO j = 1, tg%je
      DO i = 1, tg%ie
        field(i,j,1) = ff_tmp(i+hfw_m_nb,j+hfw_m_nb)
      ENDDO
    ENDDO

    DEALLOCATE( ff_tmp, STAT = errorcode )
    IF(errorcode/=0) CALL abort_extpar('Cant deallocate ff_tmp')

  END SUBROUTINE hfilter_orography

  !---------------------------------------------------------------------------

END MODULE mo_oro_filter

