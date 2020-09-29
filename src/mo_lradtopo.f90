!+  Fortran module with routines and settings for lradtopo parameters
!
! Description:
! Fortran module with routines and settings for the computation of external
! parameters related to orography correction of the radiation
!
! Current Code Owner: MeteoSwiss, Anne Roches
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! @VERSION@    @DATE@     Anne Roches
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with routines and settings for lradtopo
!> \author Anne Roches
!>
MODULE mo_lradtopo
  
  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4
  USE mo_utilities_extpar,      ONLY: &
       &                              phi2phirot, &
       &                              rla2rlarot, &
       &                              free_un

  USE mo_cosmo_grid,            ONLY: nborder, &
       &                              res_in, &
       &                              cosmo_grid, &
       &                              lon_rot, &
       &                              lat_rot

  USE mo_icon_grid_data,        ONLY: icon_grid, & !< structure which contains the definition of the ICON grid
                                      icon_grid_region

  USE mo_grid_structures,       ONLY: target_grid_def

  USE mo_search_icongrid,       ONLY: find_nc

  USE mo_icon_grid_data,        ONLY: icon_grid_region, &
       &                              nvertex_per_cell, &
       &                              icon_dom_def

  USE mo_base_geometry,         ONLY: geographical_coordinates, &
       &                              cartesian_coordinates, &
       &                              vector_product
  USE mo_topo_tg_fields,        ONLY: add_parameters_domain
  !USE mo_math_constans          ONLY: pi

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_namelists_extpar_lradtopo, &
       &    compute_lradtopo, &
       &    lradtopo_icon, &
       &    deghor, pi, rad2deg, deg2rad

  REAL(KIND=wp)            :: deghor !< number of degrees per sector [deg]

  REAL(KIND=wp), PARAMETER :: &  
       &                      pi             = 3.14159265359_wp, & !< pi
       &                      rad2deg        = 180._wp/pi,       & !< radians to degrees
       &                      deg2rad        = pi/180._wp          !< degrees to radians

  CONTAINS

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  !> subroutine to read namelist for lradtopo in EXTPAR
  SUBROUTINE read_namelists_extpar_lradtopo(namelist_file,        &
       &                                    lradtopo,             &
       &                                    nhori,                &
       &                                    radius,               &
       &                                    min_circ_cov,         &
       &                                    max_missing)



    CHARACTER (len=*), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings
    !  lradtopo

    LOGICAL,          INTENT(OUT) :: lradtopo        !< parameters for lradtopo to be computed? (TRUE/FALSE)

    INTEGER(KIND=i4), INTENT(OUT) :: nhori,  &       !< number of sectors for the horizon computation 
         &                           radius, &       !< radius [m] considered for horizon computation
         &                           min_circ_cov    !< only consider every min_circ_cov point at circumference 

    REAL(KIND=wp), INTENT(OUT)    :: max_missing     !< max missing values per nhori for each cell        

    !> local variables
    INTEGER(KIND=i4)  :: nuin, &     !< unit number
         &               ierr, &     !< error flag
         &               nhori_d, &  
         &               radius_d,&
         &               min_circ_cov_d

    REAL(KIND=wp)     :: max_missing_d

    LOGICAL           :: lradtopo_d

    !> define the namelist group
    NAMELIST /radtopo/ lradtopo, nhori, radius, min_circ_cov, max_missing

    !> initialization
    ierr            = 0

    !> default values definition
    lradtopo_d      = .FALSE.
    nhori_d         = 24
    radius_d        = 40000
    min_circ_cov_d  = 1
    max_missing_d   = 0.9_wp

    !> default values attribution
    lradtopo        = lradtopo_d 
    nhori           = nhori_d
    radius          = radius_d
    min_circ_cov    = min_circ_cov_d
    max_missing     = max_missing_d


    !> read namelist  
    nuin = free_un()  ! function free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=radtopo, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist radtopo',__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

    !> check values for consistency
    IF ( nhori > 24_i4 ) THEN
      CALL logging%warning('nhori larger than 24 is not recommended')
    ENDIF

  END SUBROUTINE read_namelists_extpar_lradtopo

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  !> subroutine to compute the lradtopo parameters in EXTPAR
  SUBROUTINE compute_lradtopo(nhori,tg,hh_topo,slope_asp,slope_ang,horizon,skyview)

    !> \author Buzzi, Luethi

    ! History:
    ! Version    Date       Name
    ! ---------- ---------- ----
    !    -        2006        Matteo Buzzi:  initial version based
    !                                        on an IDL code
    !             2012        Anne Roches:   cleanup, new interpolation
    !                                        routine written by O. Fuhrer
    !             2013        Daniel Luethi: almost complete rewrite changing
    !                                        the grid representation paradigme       
    !             2013        Anne Roches:   cleanup and documentation
    !---------------------------------------------------------------------------

    INTEGER(KIND=i4),      INTENT(IN) :: nhori                !< number of sectors for the horizon computation 
    TYPE(target_grid_def), INTENT(IN) :: tg                   !< structure with target grid description
    REAL(KIND=wp),         INTENT(IN) :: hh_topo(:,:,:)       !< mean height 
    REAL(KIND=wp),         INTENT(OUT):: slope_asp(:,:,:), &
         &                               slope_ang(:,:,:), &
         &                               horizon  (:,:,:,:), &
         &                               skyview  (:,:,:)

    !> local variables
    INTEGER (KIND=i4)                 :: i, j, errorcode, &
         &                               nsec                !< number of gridpoints per sector (in both x and y dir.)
    REAL(KIND=wp)                     :: rlon_np, rlat_np, & !< location of true North Pole in rot. coord.
                                         rdx, rdy,         & !< distance from the sector center in x and y dir.
                                         dhdx, dhdy,       & !< slope in x and y directions resp.
                                         coslat              !< cosine of rotated latitude

    REAL(KIND=wp), ALLOCATABLE        :: zhh(:,:),                   & !< local mean height
         &                               zslope_asp(:,:),            & !< local slope aspect
         &                               zslope_ang(:,:),            & !< local slope angle
         &                               zhorizon  (:,:,:),          & !< local horizon
         &                               zskyview  (:,:),            & !< local skyview
         &                               slope_x(:,:), slope_y(:,:), & !< slope in x and y directions resp.
         &                               h_hres(:,:),                & !< corrected height above see level of target grid
         &                               h_corr(:,:),                & !< height correction (Earth's curvature)
         &                               dist(:,:),                  & !< distance from sector center (local search grid)
         &                               aberr(:,:)                    !< angle between geographic meridian and grid meridian

    !> parameters
    REAL(KIND=wp), PARAMETER          :: semimaj = 6378137.0         !< semimajor radius WGS 84

    !---------------------------------------------------------------------------
    CALL logging%info('Enter routine: compute_lradtopo')
    errorcode = 0
    nsec   = 1 + 2 * nborder
    deghor = 360.0_wp / nhori

    !> allocations and initializations
    ALLOCATE( h_hres(tg%ie,tg%je), STAT=errorcode )
    IF ( errorcode /= 0 ) CALL logging%error( 'Cant allocate the array h_hres',__FILE__,__LINE__)

    ALLOCATE( zhh(tg%ie,tg%je), STAT=errorcode )
    IF ( errorcode /= 0 ) CALL logging%error( 'Cant allocate the array zhh',__FILE__,__LINE__ )

    ALLOCATE( h_corr(nsec,nsec), dist(nsec,nsec), STAT=errorcode )
    IF ( errorcode /= 0 ) CALL logging%error( 'Cant allocate the arrays h_corr and dist',__FILE__,__LINE__ )

    ALLOCATE( zslope_asp(tg%ie-2*nborder,tg%je-2*nborder),        &
         &    zslope_ang(tg%ie-2*nborder,tg%je-2*nborder),        &
         &    zhorizon  (tg%ie-2*nborder,tg%je-2*nborder,nhori),  &
         &    zskyview  (tg%ie-2*nborder,tg%je-2*nborder),        &
         &    slope_x   (tg%ie-2*nborder,tg%je-2*nborder),        &
         &    slope_y   (tg%ie-2*nborder,tg%je-2*nborder),        &
         &    aberr     (tg%ie-2*nborder,tg%je-2*nborder),        &
         &    STAT = errorcode )
    IF ( errorcode /= 0 ) CALL logging%error( 'Cant allocate the lradtopo arrays',__FILE__,__LINE__ )

    h_hres      (:,:)   = 0.0_wp
    h_corr      (:,:)   = 0.0_wp
    zslope_asp  (:,:)   = 0.0_wp
    zslope_ang  (:,:)   = 0.0_wp
    zhorizon    (:,:,:) = 0.0_wp
    zskyview    (:,:)   = 0.0_wp
    slope_x     (:,:)   = 0.0_wp
    slope_y     (:,:)   = 0.0_wp
    aberr       (:,:)   = 0.0_wp

    ! copy the orography in a 2D local variable
    zhh(:,:) = hh_topo(:,:,1)

    !> position of true North Pole in rotated coordinates
    rlat_np = phi2phirot( 90.0_wp, 0.0_wp, cosmo_grid%pollat, cosmo_grid%pollon )
    rlon_np = rla2rlarot( 90.0_wp, 0.0_wp, cosmo_grid%pollat, cosmo_grid%pollon, cosmo_grid%polgam  )

    !> distance matrix for all points within one sector 
    !        (distance from the sector center)
    DO j = 1, nsec
      rdy = REAL(j-nborder-1) !_br 04.04.14
      DO i = 1, nsec
        rdx = REAL(i-nborder-1) !_br 04.04.14
        dist(i,j) = SQRT(rdx*rdx + rdy*rdy)
      ENDDO
    ENDDO
    h_corr(:,:) = semimaj * (1.0_wp / COS(dist(:,:)*res_in/semimaj) - 1.0_wp)

    !> loop over all target gridpoints: computation of the output quantities
    DO j = 1, tg%je - 2 * nborder
      coslat = COS( deg2rad * lat_rot(j) )
      ! -> this is not working correct -> !$OMP PARALLEL DO PRIVATE(i,dhdx,dhdy)
      DO i = 1, tg%ie - 2 * nborder

        ! compute corrected height within the search radius of horizon
        h_hres(i:i+2*nborder,j:j+2*nborder) = MAX(0.0_wp, zhh(i:i+2*nborder,j:j+2*nborder) - h_corr(:,:))

        ! compute angle between y axis and true north in rotated grid
        aberr(i,j) = rad2deg * atan2(rlon_np - lon_rot(i+nborder), rlat_np - lat_rot(j+nborder))

        ! compute slope in both x and y direction
        dhdx = 0.5_wp * ( zhh(i+nborder+1,j+nborder  ) - zhh(i+nborder-1,j+nborder) ) / (res_in * coslat)
        dhdy = 0.5_wp * ( zhh(i+nborder  ,j+nborder+1) - zhh(i+nborder,j+nborder-1) ) / res_in

        slope_x(i,j) = dhdx * COS( deg2rad * aberr(i,j) ) - dhdy * SIN( deg2rad * aberr(i,j) )
        slope_y(i,j) = dhdy * COS( deg2rad * aberr(i,j) ) + dhdx * SIN( deg2rad * aberr(i,j) )

        ! compute slope angle 
        zslope_ang(i,j) = rad2deg * ATAN( SQRT( dhdx**2 + dhdy**2 ) )

        ! compute slope aspect
        IF ( ( dhdx == 0._wp ) .AND. ( dhdy == 0._wp ) ) THEN
          zslope_asp(i,j) = 0.0_wp 
        ELSE
          zslope_asp(i,j) = 90._wp - rad2deg * ATAN2( -dhdy, -dhdx )
        ENDIF
        IF ( zslope_ang(i,j) >  0.001_wp) THEN
          zslope_asp(i,j) = zslope_asp(i,j) - aberr(i,j)
        ENDIF
        IF ( zslope_asp(i,j) >= 360.0_wp) THEN
          zslope_asp(i,j) = zslope_asp(i,j) - 360.0_wp
        ENDIF
        IF ( zslope_asp(i,j) <    0.0_wp) THEN
          zslope_asp(i,j) = zslope_asp(i,j) + 360.0_wp
        ENDIF

        ! compute horizon
        CALL comp_horiz( h_hres(i:i+2*nborder,j:j+2*nborder), res_in,             &
             slope_x(i,j), slope_y(i,j), nhori, nsec, nsec,         &
             zhorizon(i,j,:), aberr(i,j) )

        ! compute skyview factor
        CALL compute_skyview( zslope_ang(i,j), zslope_asp(i,j), zhorizon(i,j,:),   &
             nhori, zskyview(i,j) )

      ENDDO
      ! -> this is not working correct -> !$OMP END PARALLEL DO
    ENDDO

    !> transformation in correct units
    zslope_asp(:,:) = deg2rad * zslope_asp(:,:)
    zslope_ang(:,:) = deg2rad * zslope_ang(:,:)

    !> copy 2D local output variables into the 3D output variables
    slope_asp(nborder+1:nborder+tg%ie-2*nborder,nborder+1:nborder+tg%je-2*nborder,1)   = zslope_asp(:,:)
    slope_ang(nborder+1:nborder+tg%ie-2*nborder,nborder+1:nborder+tg%je-2*nborder,1)   = zslope_ang(:,:)
    horizon  (nborder+1:nborder+tg%ie-2*nborder,nborder+1:nborder+tg%je-2*nborder,1,:) = zhorizon(:,:,:)
    skyview  (nborder+1:nborder+tg%ie-2*nborder,nborder+1:nborder+tg%je-2*nborder,1)   = zskyview(:,:)

    !> deallocations 
    DEALLOCATE( h_hres, STAT=errorcode )
    IF ( errorcode /= 0 ) CALL logging%error( 'Cant deallocate the array h_hres',__FILE__,__LINE__ )

    DEALLOCATE( zhh, STAT=errorcode )
    IF ( errorcode /= 0 ) CALL logging%error( 'Cant deallocate the array zhh',__FILE__,__LINE__ )

    DEALLOCATE( h_corr, dist, STAT=errorcode )
    IF ( errorcode /= 0 ) CALL logging%error( 'Cant deallocate the arrays h_corr and dist',__FILE__,__LINE__ )

    DEALLOCATE(zslope_asp,                         &
         &     zslope_ang,                         &
         &     zhorizon  ,                         &
         &     zskyview  ,                         &
         &     slope_x   ,                         &
         &     slope_y   ,                         &
         &     aberr     ,                         &
         STAT = errorcode )
    IF ( errorcode /= 0 ) CALL logging%error( 'Cant deallocate the lradtopo arrays',__FILE__,__LINE__ )

    CALL logging%info('Exit routine: compute_lradtopo')

  END SUBROUTINE compute_lradtopo

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  !> subroutine to compute the lradtopo parameters in EXTPAR for Icon
  SUBROUTINE lradtopo_icon(nhori,radius, min_circ_cov, tg,hh_topo,horizon, skyview, max_missing)

    !  This routine computes the horizon (assuming plain surfaces) for 
    !  unrotated Icon grids. Additionally a simple geometric skyview-factor is
    !  computed too. For more robust results each nhori-sector is further
    !  subdivided. The azimuth-angle rotates clockwise, starting at N.

  
    INTEGER(KIND=i4),      INTENT(IN) :: nhori, &             !< number of sectors for the horizon computation 
         &                               radius, &            !< search radius [m]
         &                               min_circ_cov         !< factor to reduce coverage at outermost points

    TYPE(target_grid_def), INTENT(IN) :: tg                   !< structure with target grid description
    REAL(KIND=wp),         INTENT(IN) :: hh_topo(:,:,:),   &  !< mean height 
         &                               max_missing          !< max missing values per nhori for each cell

    REAL(KIND=wp),         INTENT(OUT):: horizon(:,:,:,:), &  !< horizon angle [deg]
         &                               skyview(:,:,:)       !< skyview-factor [-]

    !> local variables
    INTEGER (KIND=i4)                 :: i, j, k,nh, errorcode,&
         &                               size_radius,       & !< number of gridcells along radius
         &                               start_cell_id,     & !< id of cell to start find_nc 
         &                               nearest_cell_id,   & !< return id of find_nc
         &                               refine_factor,     & !< refinement factor for nhorinhori_iter, 
         &                               nhori_iter           !< number of "true" nhori for computation

    REAL(KIND=wp)                     :: deghor,            & !< azimut increment
         &                               icon_resolution,   & !< aprox. icon grid resolution
         &                               domain_center_lat, & !< center of Icon grid
         &                               domain_center_lon, & !< center of Icon grid
         &                               rlat_cell,rlon_cell,&!< radians of target cell 
         &                               angle_deg,         & !< azimuth angle for each nhori
         &                               skyview_sum,       & !< helper
         &                               percentage_of_fails  !< % of failed find_nc

    REAL(KIND=wp), ALLOCATABLE        :: zhh(:),                   & !< local mean height
         &                               zhorizon  (:,:),          & !< local horizon
         &                               h_corr(:),                & !< height correction (Earth's curvature)
         &                               dlat(:,:),                & !< dlat across search_radius
         &                               dlon(:,:),                & !< dlat across search_radius
         &                               dh(:),                    & !< horizontal difference [m]
         &                               dz(:),                    & !< vertical difference [m]
         &                               angle(:),                 & !< azimuth angle for each nhori [deg]
         &                               rates(:),                 & !< dz/dh across search_radius
         &                               zskyview(:),              & !< local skyview
         &                               z_search_radius(:),       & !< debug field to visualize search_radius
         &                               z_missing_data(:,:)         !< field to store fraction of missing data
    
    TYPE(geographical_coordinates)    :: target_co_rad,            & !< rad. coord. of target cell
         &                               target_co_deg               ! geog. coord. of target cell

    !> parameters
    REAL(KIND=wp), PARAMETER          :: semimaj = 6378137.0         !< semimajor radius WGS 84

    !---------------------------------------------------------------------------
    CALL logging%info('Enter routine: lradtopo_icon')
    errorcode = 0

    ! aprox. horizontal resolution
    icon_resolution = 5050.e3_wp/(icon_grid%grid_root*2**icon_grid%grid_level)

    size_radius = NINT(radius/icon_resolution)

    CALL determine_domain_center(tg, domain_center_lat, domain_center_lon, size_radius, icon_resolution, semimaj)

    ! subdivide nhori for more robust results:
    refine_factor = INT( (2 * pi * size_radius / (min_circ_cov * nhori) ) )

    nhori_iter = refine_factor * nhori

    deghor = 360._wp/nhori_iter

    ! info prints
    WRITE(message_text, '(A,F6.1,A)'), ' Grid resolution for lradtopo-computation is ', &
         & icon_resolution, ' m'
    CALL logging%info(message_text)

    IF ( refine_factor > 1 ) THEN
      WRITE(message_text, '(A,I3)'), ' Subdivide each nhori-sector further by ', refine_factor
      CALL logging%info(message_text)
    ENDIF

    !> allocations and initializations
    ALLOCATE( zhh            (tg%ie),  &
         &    angle          (tg%ie),  &
         &    z_search_radius(tg%ie),  &
         &    zskyview       (tg%ie),  &
         &    STAT = errorcode )
    IF ( errorcode /= 0 ) THEN
      CALL logging%error( 'Cant allocate the lradtopo arrays with dim(tg%ie)',__FILE__,__LINE__ )
    ENDIF

    ALLOCATE( zhorizon       (tg%ie,nhori),  &
         &    z_missing_data (tg%ie, nhori), &
         &    STAT = errorcode )
    IF ( errorcode /= 0 ) THEN
      CALL logging%error( 'Cant allocate the lradtopo arrays with dim(tg%ie,nhori)',__FILE__,__LINE__ )
    ENDIF

    ALLOCATE( dh     (size_radius),          &
         &    dz     (size_radius),          &
         &    rates  (size_radius),          &
         &    h_corr (size_radius),          &
         &    STAT = errorcode )
    IF ( errorcode /= 0 ) THEN
      CALL logging%error( 'Cant allocate the lradtopo arrays with dim(size_radius)',__FILE__,__LINE__ )
    ENDIF

    ALLOCATE( dlat(size_radius, nhori_iter),          &
         &    dlon(size_radius, nhori_iter),          &
         &    STAT = errorcode )
    IF ( errorcode /= 0 ) THEN
      CALL logging%error( 'Cant allocate the lradtopo arrays with dim(size_radius,nhori_iter)',__FILE__,__LINE__ )
    ENDIF

    ! copy the orography in a 1D local variable
    zhh(:) = hh_topo(:,1,1)

    h_corr(:)     = 0.0_wp
    zhorizon(:,:) = 0.0_wp
    z_missing_data(:,:) = 0.0_wp

    ! debug fields
    z_search_radius(:) = -5.0_wp

    ! calculate dlat/dlon on vector from center of circle for each nhori_iter
    DO j=1, nhori_iter
      angle_deg = 0 + deghor * (j-1)

      CALL distance_relative_to_cell(domain_center_lat, &
           &                         domain_center_lon, &
           &                         angle_deg,         &
           &                         size_radius,       &
           &                         dlat(:,j),         &
           &                         dlon(:,j),         &
           &                         icon_resolution,   &
           &                         semimaj)
    ENDDO
    
    ! height correction because of earth's curvature 
    DO i= 1, size_radius
      dh(i) = i * icon_resolution
      h_corr(i) = SQRT(semimaj**2 + dh(i)**2) -semimaj
    ENDDO

    ! loop over all cells
    DO i=1,tg%ie 

      !get_coordinates of cell in radians
      rlat_cell = icon_grid_region%cells%center(i)%lat
      rlon_cell = icon_grid_region%cells%center(i)%lon

      nh = 0
      ! iterate over refined horizons clockwise
      DO j= 1, nhori_iter

        ! increase nhori after "refine_factor"-times
        IF (MODULO(j,refine_factor ) == 1 .OR. nhori_iter == nhori) THEN
          nh = nh + 1
        ENDIF

        ! reset start_cell index
        start_cell_id = i

        DO k = 1, size_radius

          ! lat/lon of point on radius in radians
          target_co_rad%lat = rlat_cell - dlat(k,j) * deg2rad
          target_co_rad%lon = rlon_cell - dlon(k,j) * deg2rad

          ! convert to degress
          target_co_deg%lat = rad2deg * target_co_rad%lat
          target_co_deg%lon = rad2deg * target_co_rad%lon

          
          ! check if point on radius is outside of domain
          IF (   target_co_deg%lat > tg%maxlat .OR. &
               & target_co_deg%lat < tg%minlat .OR. &
               & target_co_deg%lon > tg%maxlon .OR. &
               & target_co_deg%lon < tg%minlon) THEN

             ! mark points as outside domain for later handling in the code
             z_missing_data(i, nh) = z_missing_data(i, nh) + 1

             ! set height difference to center cell to 0
             dz(k) = 0.0_wp
             
          ELSE ! point inside domain

            ! walk to nearest cell for given "first try" start_cell_id
            CALL find_nc(target_co_rad,    &
              &          nvertex_per_cell, &
              &          icon_dom_def,     &
              &          icon_grid_region, &
              &          start_cell_id,    &
              &          nearest_cell_id)

            
            ! if no nearest cell could be determined index returned is 0
            IF ( nearest_cell_id == 0) THEN
              z_missing_data(i, nh) = z_missing_data(i, nh) + 1

              ! set height difference to center cell to 0
              dz(k) = 0.0_wp

              ! reset start_cell_id
              start_cell_id = i

            ELSE ! nearest_cell_id found
              
              ! height difference to center cell with correction of earth's
              ! curvature
              dz(k)=  MAX( (zhh(nearest_cell_id)-h_corr(k)) - zhh(i), 0.0_wp)

              ! visualize search-radius (for debugging)
              IF (i == (tg%ie/2) ) THEN
                IF (z_search_radius(nearest_cell_id) <= -4._wp)THEN
                  z_search_radius(nearest_cell_id) = 1
                ELSE
                  z_search_radius(nearest_cell_id) = z_search_radius(nearest_cell_id) + 1
                ENDIF
              ENDIF

            ENDIF

          ENDIF ! point inside domain

        ENDDO ! size_radius

        ! rates along radius
        rates(:) = dz(:) / dh(:)

        ! zhorizon unfiltered in degrees
        zhorizon(i,nh) = zhorizon(i,nh) + rad2deg * ATAN(MAXVAL(rates))
      ENDDO ! nhori_iter

      ! average horizon
      zhorizon(i,:) = zhorizon(i,:) / refine_factor

      ! norm missingness by number of cell considered for one horizon value
      z_missing_data(i,:) = z_missing_data(i,:) / REAL(size_radius * refine_factor)

    ENDDO ! tg%ie

    ! set all horizon angles to 0 for missingness-fraction above max_missing
    DO nh=1, nhori
      DO i = 1, tg%ie 
        IF ( z_missing_data(i,nh) > max_missing) THEN
          zhorizon(i,nh)= 0.0_wp 
        ENDIF
      ENDDO
    ENDDO

    WRITE(message_text,*), ' Total missing points for nhori:'
    CALL logging%info(message_text)

    DO nh= 1,nhori
      percentage_of_fails = 100 * SUM( z_missing_data(:,nh) )/ REAL(tg%ie)
      WRITE(message_text,'(A,I3,A,F6.1,A)'), '  ', nh, ': ', percentage_of_fails, '%'
      CALL logging%info(message_text)
    ENDDO
      
    !---------------------------------------------------------------------
    ! compute skyview-factor (experimental)
    DO i=1, tg%ie
      skyview_sum = 0.0_wp
      DO nh=1, nhori

        ! pure geometric skyview-factor
        skyview_sum = skyview_sum + (1 - SIN(zhorizon(i,nh) * deg2rad))

        ! geometric scaled with sin(horizon)
        !skyview_sum = skyview_sum + (1 - (SIN(zhorizon(i,nh)*deg2rad) * SIN(zhorizon(i,nh)*deg2rad)))

        ! geometric scaled with sin(horizon)**2
        !skyview_sum = skyview_sum + (1 - (SIN(zhorizon(i,nh)*deg2rad) * SIN(zhorizon(i,nh)*deg2rad)**2))

        ! pure cos(horizon)
        !skyview_sum = skyview_sum + COS(zhorizon(i,nh)*deg2rad) 

        ! pure cos(horizon)**2
        !skyview_sum = skyview_sum + COS(zhorizon(i,nh)*deg2rad)**2 
      ENDDO
      zskyview(i) = skyview_sum / REAL(nhori)
    ENDDO

    ! assign local fields to output fields
    horizon(:,1,1,:) = zhorizon(:,:)
    skyview(:,1,1) = zskyview(:)

    !> deallocations 
    DEALLOCATE( zhh          ,  &
         &    angle          ,  &
         &    z_search_radius,  &
         &    zskyview       ,  &
         &    STAT = errorcode )
    IF ( errorcode /= 0 ) THEN
      CALL logging%error( 'Cant deallocate the lradtopo arrays with dim(tg%ie)',__FILE__,__LINE__ )
    ENDIF

    DEALLOCATE( zhorizon     , &
         &    z_missing_data , &
         &    STAT = errorcode )
    IF ( errorcode /= 0 ) THEN
      CALL logging%error( 'Cant deallocate the lradtopo arrays with dim(tg%ie,nhori)',__FILE__,__LINE__ )
    ENDIF

    DEALLOCATE( dh           , &
         &    dz             , &
         &    rates          , &
         &    h_corr         , &
         &    STAT = errorcode )
    IF ( errorcode /= 0 )THEN
      CALL logging%error( 'Cant deallocate the lradtopo arrays with dim(size_radius)',__FILE__,__LINE__ )
    ENDIF

    DEALLOCATE( dlat,          &
         &      dlon,          &
         &    STAT = errorcode )
    IF ( errorcode /= 0 ) THEN
      CALL logging%error( 'Cant deallocate the lradtopo arrays with dim(size_radius,nhori_iter)',__FILE__,__LINE__ )
    ENDIF

    CALL logging%info('Exit routine: lradtopo_icon')

  END SUBROUTINE lradtopo_icon

  !> subroutine to compute the horizon
  SUBROUTINE comp_horiz( h_hres, dhres, dhdx, dhdy, nhori, nx, ny, hor, rot_ang )

    ! Compute horizon angles (considering self-shading) 
    ! inside the domain given by h_hres
    ! Matteo Buzzi, 2009, Meteoswiss
    ! Daniel Luethi, 2012, IAC
    ! AnneRoches, 2013, C2SM

    !> arguments
    INTEGER(KIND=i4), INTENT(IN) :: nx, ny, &     ! gridpoints number in the x and y dir
         &                          nhori       !  number of sectors

    REAL   (KIND=wp), INTENT(IN) :: h_hres(:,:), & ! topography [m]
         &                          dhres,                 & ! resolution [m] 
         &                          dhdx , dhdy,           & ! slope in x and y dir resp. [-]
         &                          rot_ang                  ! rotation angle (=angle between grid meridian and geo. meridian) [deg]
    REAL   (KIND=wp), INTENT(OUT):: hor(:)               ! horizon angle [deg]

    !> local variables
    INTEGER(KIND=i4)             :: x0, y0,    & ! center of the sector in x and y dir resp.
                                    i,  j , k, & ! loop indices
                                    ngp,       & ! number of gridpoints at each side of the sector center
                                    nzstat       ! allocation status
                                 
    REAL(KIND=wp)                :: tmp,       & ! temporary horizon [deg]
                                    zshift,    & ! angle zshift [deg]
                                    azi,       & ! angle between the geo North and the sector start  [deg]
                                    ssa,       & ! self-shading [deg]
                                    ho           ! horizon (maximum elevation angle in a portion of the sector)
                                     ! without self-shading [deg]

    INTEGER(KIND=i4),ALLOCATABLE :: xcart(:), ycart(:)

    REAL(KIND=wp)   ,ALLOCATABLE :: pcr(:,:),    & ! polar coordinates (index1: angle, index2: radius)
         &                          rcr(:,:),    & ! rectangular coord.(index1: x    , index2: y     )
         &                          dh(:),dn(:), & ! distance between the sector center and the other points of 
                                                    ! interest in the sector in [m] in the vertical and in the horizontal resp.
         &                          rates(:)       ! slopes between the sector center and the other points of
                                                    ! interest in the sector

    INTEGER(KIND=i4), PARAMETER  :: niter = 10 ! number of iterations

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    !> preparations
    x0     = nx/2 +1
    y0     = ny/2 +1
    ngp    = nx/2

    !> allocations
    ALLOCATE( pcr(2,ngp), rcr(2,ngp), xcart(ngp), ycart(ngp), dh(ngp), dn(ngp), rates(ngp), STAT=nzstat )
    IF ( nzstat /= 0 ) CALL logging%error( 'Cant allocate fields',__FILE__,__LINE__)

    !> radius for polar coordinates
    DO k = 1, ngp
      pcr(2,k) = k
    ENDDO

    !> loop over all sectors to compute the horizon
    DO i = 1, nhori 
      tmp = 0.0_wp
      ! iterations by zshifting slightly the sector in order
      ! to get more robust results 
      DO j = 1, niter+1
        zshift = -deghor/2.0_wp + (j-1) * deghor/niter
        azi   = deghor * (i-1) + rot_ang + zshift 
        ! convert relative polar coordinates to relative
        ! rectangular (cartesian) coordinates
        ! (these are the coordinates relative to the sector center)
        pcr(1,1:ngp) = deg2rad      * (90.0_wp - azi)
        rcr(1,1:ngp) = pcr(2,1:ngp) * COS(pcr(1,1:ngp)) 
        rcr(2,1:ngp) = pcr(2,1:ngp) * SIN(pcr(1,1:ngp))
        ! absolute cartesian coordinates
        ! (these are the coordinates of the points relative
        ! to the lower left point of the sector)
        xcart(:) = NINT(rcr(1,1:ngp)) + x0 
        ycart(:) = NINT(rcr(2,1:ngp)) + y0 

        ! compute vertical and horizontal differences
        ! between sector center and all the points located
        ! in the sector part of interest 
        DO k = 1, ngp
          dh(k) = h_hres(xcart(k),ycart(k)) - h_hres(x0,y0)
          dn(k) = dhres * SQRT(REAL( (xcart(k)-x0)**2 + (ycart(k)-y0)**2 ))
        ENDDO
        rates(:) = dh(:) / dn(:)

        ! compute maximal horizon angle
        ho = rad2deg * ATAN(MAXVAL(rates))

        ! compute self-shading
        ssa = rad2deg * ATAN( SIN(deg2rad*azi)*dhdx + COS(deg2rad*azi)*dhdy )
        ssa = MAX(ssa, 0.0_wp)

        ! take the maximum between the maximal terrain angle and self-shading
        tmp = tmp + MAX( ssa, ho )
      ENDDO
      hor(i) = tmp / (niter+1)
    ENDDO

    !> cleanup
    DEALLOCATE( pcr, rcr, xcart, ycart, dh, dn, rates, STAT=nzstat )
    IF ( nzstat /= 0_i4 ) STOP

    IF ( nzstat /= 0 ) CALL logging%error( 'Cant deallocate fields',__FILE__,__LINE__)

  END SUBROUTINE comp_horiz

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  !> subroutine to compute the skyview
  SUBROUTINE compute_skyview( slope_ang, slope_asp, horizon, nhori, skyview )

    ! This procedure calculates the skyview factor for a tilted
    ! plane surrounded by a non negligible horizon.
    ! The horizon  elevations are assumed to be constant within each of the 
    ! nhori equally sized sectors but they can vary from sector to sector.
    ! The sectors start in N and increase clockwise over E to S and W.
    !
    ! Antoine Zelenka, 2005, MeteoSwiss
    ! Matteo Buzzi, 2009, MeteoSwiss
    ! Anne Roches, 2013, C2SM

    !> arguments
    REAL(KIND=wp), INTENT(IN)                   :: slope_ang, & ! slope angle [deg]
         &                                         slope_asp, &    ! slope aspect [deg]
         &                                         horizon(:)   ! horizon [deg]   

    INTEGER(KIND=i4), INTENT(IN)                :: nhori        ! number of sectors 

    REAL(KIND=wp), INTENT(OUT)                  :: skyview      ! skyview factor [-]

    !> local variables
    INTEGER(KIND=i4)                            :: i, i0, k, ip, count_clip_low, count_clip_high
                                               
    REAL(KIND=wp)                               :: rslope_ang,             & ! slope angle [rad]
         &                                         rdeghor,                & ! number of radians per sector [rad]
         &                                         cosa, cosah, sinah,     & ! cosine, half cosine, half sine of the slope angle resp. [-] 
         &                                         sf,                     & !
         &                                         azi,                    & ! azimut of the trailing trace of the slope on the horizontal
                                                                             ! plane [deg]
         &                                         slg                       !

    REAL(KIND=wp), DIMENSION(nhori+1)          :: ang_start,            & ! angle of the sector start  [rad]
         &                                        hop,                  & !
         &                                        sl,                   & !
         &                                        hpl

    REAL(KIND=wp), DIMENSION(2*nhori+1)        :: drhorizon             ! duplicated horizon [rad]

    REAL(KIND=wp), PARAMETER                   :: zepsilon = 2.0_wp * pi * 0.01_wp ! security parameter 

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    !> preparations
    count_clip_low = 0
    count_clip_high = 0

    ! convert all variables from degrees to radians
    rslope_ang                 = ( slope_ang + 0.00001 ) * deg2rad
    drhorizon(1:nhori)         = horizon(:)              * deg2rad
    drhorizon(nhori+1:2*nhori) = horizon(:)              * deg2rad
    drhorizon(2*nhori+1)       = horizon(1)              * deg2rad
    rdeghor                    = deghor                  * deg2rad

    ! angle of the sector start
    DO i = 1, nhori+1 
      ang_start(i) = rdeghor * (i-1)
    ENDDO

    ! angle between the grid North and the main slope
    azi = slope_asp - 90.0_wp 
    IF ( azi < 0.0_wp ) azi = azi + 360.0_wp
    azi = azi * deg2rad

    ! selection of the sector containing the main slope
    i  = FLOOR(azi/rdeghor)
    ! end of next sector
    i0 = i+2

    ! angle between the end of the next sector and the main slope
    sl(1) = ang_start(i0) - azi
    DO k = 2, nhori+1
      sl(k) = sl(k-1) + rdeghor
    ENDDO

    hop(1:nhori+1) = drhorizon(i0:i0+nhori)

    cosa  = COS(rslope_ang) 
    cosah = 0.5_wp * cosa
    sinah = 0.5_wp * SIN(rslope_ang)

    ! Treat the vertical plane separately because of cos(rslope_ang)=0 in the
    ! formula for HPL
    IF ( slope_ang /= 90.0_wp ) THEN
      hpl     = ATAN( -SIN(sl) * sinah/cosah )
      skyview = 0.0
      DO i =1,nhori
        ip = i + 1
        IF ( hop(i) > hpl(i) ) THEN 
          IF ( hop(i) >= hpl(ip) ) THEN 
            skyview = skyview + sffront(cosah, sinah, hop(i), sl(ip), sl(i)) 
          ELSE 
            slg = ASIN(TAN(hop(i))/TAN(rslope_ang)) + pi
            ! SLG too high
            IF ( (slg > sl(ip)) .AND. (slg > sl(ip)+zepsilon) ) THEN 
              WRITE(message_text,*) 'slope_ang: ', slope_ang, ' slope_asp: ', slope_asp, ' horizon: ', horizon
              CALL logging%warning(message_text)
              WRITE(message_text,*) 'slg:  ', slg, ' sl(ip): ', sl(ip), ' sl(i): ', sl(i), ' hop(i): ', hop(i), ' hpl(i): ', hpl(i)
              CALL logging%warning(message_text)
              WRITE(message_text,*)  ' hpl(ip): ', hpl(ip), ' i:  ',i, ' ip: ', ip
              CALL logging%warning(message_text)
              CALL logging%error('SLG too high',__FILE__,__LINE__)
            ENDIF
            ! SLG too low
            IF ((slg < sl(i)) .AND. (slg < sl(i)-zepsilon)) THEN 
              WRITE(message_text,*) 'slope_ang: ', slope_ang, ' slope_asp: ', slope_asp, ' horizon: ', horizon
              CALL logging%warning(message_text)
              WRITE(message_text,*) 'slg:  ', slg, ' sl(ip): ', sl(ip), ' sl(i): ', sl(i), ' hop(i): ', hop(i), ' hpl(i): ', hpl(i)
              CALL logging%warning(message_text)
              WRITE(message_text,*)  ' hpl(ip): ', hpl(ip), ' i:  ',i, ' ip: ', ip
              CALL logging%warning(message_text)
              CALL logging%error('SLG too low',__FILE__,__LINE__)
            ENDIF
            ! SLG in tolerance zone -> clipped to lower limit
            IF ( (slg > sl(ip)) .AND. (slg < sl(ip)+zepsilon) ) THEN 
              count_clip_low = count_clip_low + 1
              WRITE(message_text,*) 'SLG clipped to lower limit','slg before: ', slg, ' slg clipped: ', sl(ip)
              slg = sl(ip)
            ENDIF
            ! SLG in tolerance zone -> clipped to upper limit
            IF ( (slg < sl(i)) .AND. (slg > sl(i)-zepsilon) ) THEN
              count_clip_high = count_clip_high + 1
              WRITE(message_text,*) '!! SLG clipped to upper limit !!','slg before: ', slg, ' slg clipped: ', sl(i)
              slg = sl(i)
            ENDIF

            sf      = sffront( cosah, sinah, hop(i), slg, sl(i) )
            sf      = sf + sfback( cosa, sinah, hpl(ip), hop(i), sl(ip), slg )
            skyview = skyview + sf
          ENDIF
        ELSE
          IF ( hop(i) <= hpl(ip) ) THEN 
            skyview = skyview + sfback( cosa, sinah, hpl(ip), hpl(i), sl(ip), sl(i) ) 
          ELSE 
            slg = 2.0_wp * pi - ASIN(TAN(hop(i))/TAN(rslope_ang))
            ! slG too high
            IF ( (slg > sl(ip)) .AND. (slg > sl(ip)+zepsilon) ) THEN 
              WRITE(message_text,*) 'slope_ang: ', slope_ang, ' slope_asp: ', slope_asp, ' horizon: ', horizon
              CALL logging%warning(message_text)
              WRITE(message_text,*) 'slg:  ', slg, ' sl(ip): ', sl(ip), ' sl(i): ', sl(i), ' hop(i): ', hop(i), ' hpl(i): ', hpl(i)
              CALL logging%warning(message_text)
              WRITE(message_text,*)  ' hpl(ip): ', hpl(ip), ' i:  ',i, ' ip: ', ip
              CALL logging%warning(message_text)
              CALL logging%error('SLG too high',__FILE__,__LINE__)
            ENDIF
            ! slg too low
            IF ( (slg < sl(i)) .AND. (slg < sl(i)-zepsilon) ) THEN 
              WRITE(message_text,*) 'slope_ang: ', slope_ang, ' slope_asp: ', slope_asp, ' horizon: ', horizon
              CALL logging%warning(message_text)
              WRITE(message_text,*) 'slg:  ', slg, ' sl(ip): ', sl(ip), ' sl(i): ', sl(i), ' hop(i): ', hop(i), ' hpl(i): ', hpl(i)
              CALL logging%warning(message_text)
              WRITE(message_text,*)  ' hpl(ip): ', hpl(ip), ' i:  ',i, ' ip: ', ip
              CALL logging%warning(message_text)
              CALL logging%error('SLG too low',__FILE__,__LINE__)
            ENDIF
            ! slg in tolerance zone -> clipped to lower limit
            IF ( (slg > sl(ip)) .AND. (slg < sl(ip)+zepsilon) ) THEN 
              count_clip_low = count_clip_low + 1
              WRITE(message_text,*) 'slg clipped to lower limit','slg before: ', slg, ' slg clipped: ', sl(ip)
              slg = sl(ip)
            ENDIF
            ! slg in tolerance zone -> clipped to upper limit
            IF ( (slg < sl(i)) .AND. (slg > sl(i)-zepsilon) ) THEN 
              count_clip_high = count_clip_high + 1
              WRITE(message_text,*) '!! slg clipped to upper limit !!','slg before: ', slg, ' slg clipped: ', sl(i)
              slg = sl(i)
            ENDIF
            sf      = sfback( cosa, sinah, hop(i), hpl(i), slg, sl(i) )
            sf      = sf + sffront( cosah, sinah, hop(i), sl(ip), slg )
            skyview = skyview + sf
          ENDIF
        ENDIF
      ENDDO
      skyview = skyview / pi
    ELSE
      IF ( slope_ang == 90.0_wp ) THEN 
        skyview = 0.5_wp * pi
        !      skyview = skyview + sffront( 0.0_wp, sinah, hop(nhori-1), sl(0), 0.0_wp )!_br 21.02.14 sl field starts at index 1
        skyview = skyview + sffront( 0.0_wp, sinah, hop(nhori-1), sl(1), 0.0_wp ) !_br 21.02.14
        i = 1
        DO WHILE ( sl(i+1) < pi )
          ip = i + 1
          skyview = skyview + sffront( 0.0_wp, sinah, hop(i), sl(ip), sl(i) )
          i = ip
        ENDDO
        skyview = skyview + sffront( 0.0_wp, sinah, hop(i), pi, sl(i) )
        skyview = skyview / pi
      ENDIF
    ENDIF

    IF (count_clip_high > 0) THEN
      WRITE(message_text,*)'Number of SLG clipped to upper limit: ' , count_clip_high 
      CALL logging%info(message_text)
    ENDIF

    IF (count_clip_low > 0) THEN
      WRITE(message_text,*)'Number of SLG clipped to lower limit: ' , count_clip_low 
      CALL logging%info(message_text)
    ENDIF

  END SUBROUTINE compute_skyview

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  !> function to ???
  FUNCTION atan22( psi, cosalfa )

    REAL (KIND=wp),INTENT(IN) :: psi, cosalfa

    REAL (KIND=wp)            :: atan22
    REAL (KIND=wp)            :: pih, pi3h, x

    IF( cosalfa == 0.0_wp ) THEN 
      WRITE(logging%fileunit,*) 'WRONG argument in atan22 -> STOP'
    ENDIF
    pih  = 0.5*pi
    pi3h = 1.5*pi
    IF ( psi == pih ) THEN 
      atan22 = pih
    ENDIF
    IF ( psi == pi3h ) THEN 
      atan22 = pi3h
    ENDIF
    !Regular cases
    x = TAN(psi)/cosalfa
    IF ( psi < pih ) THEN 
      atan22 = ATAN(x) 
    ELSE
      IF ( psi < pi3h ) THEN 
        atan22 = pi + ATAN(x) 
      ELSE 
        atan22 = 2.0_wp * pi + ATAN(x)
      ENDIF
    ENDIF

  END FUNCTION atan22

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  !> function to ???
  FUNCTION sffront( c, s, h, so, su )

    REAL(KIND=wp),INTENT(IN) :: c, s, h, so, su

    REAL(KIND=wp)            :: sffront
    REAL(KIND=wp)            :: coshi, sinhi 

    coshi = cos(h)
    sinhi = sin(h)
    sffront = c * coshi**2 * (so-su) 
    sffront = sffront - s * ( h + sinhi*coshi ) * ( -COS(so) + COS(su) )

  END FUNCTION sffront

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  !> function to ???
  FUNCTION sfback( c, s, hpip, hpi, so, su )

    REAL(KIND=wp),INTENT(IN):: c, s, hpip, hpi, so, su

    REAL(KIND=wp)           :: sfback

    sfback = s * ( hpip * COS(so) - hpi * COS(su) ) 
    sfback = sfback + 0.5_wp * ( atan22(so,c) - atan22(su,c) ) 

  END FUNCTION sfback

  ! calculate radial dlat/dlon across vector for a specific bearing angle
  SUBROUTINE distance_relative_to_cell(lat_celld,lon_celld,angled,length,dlat,dlon,icon_resolution, semimaj)

    INTEGER(KIND=i4), INTENT(IN):: length            !< length of vector  
    REAL(KIND=wp),INTENT(IN)    :: lat_celld, &      !< lat of cell [deg] 
         &                         lon_celld, &      !< lon of cell [deg]
         &                         angled,    &      !< bearing angle [deg]
         &                         icon_resolution, &!< grid resolution of Icon
         &                         semimaj           !< earth-radius [m]

    REAL(KIND=wp), INTENT(OUT)  :: dlat(length),    &!< dlat [deg] along vector
         &                         dlon(length)      !< dlon [deg] along vector
    
    ! local variables
    INTEGER(KIND=i4)            :: i
    REAL(KIND=wp)               :: lat_tmp, lon_tmp,&!< helpers
         &                         rlat_cell,       &!< lat of cell [rad]
         &                         rlon_cell,       &!< lon of cell [rad]
         &                         angle             !< bearing angle [rad]

    !conversion to radians
    rlat_cell = deg2rad * lat_celld
    rlon_cell = deg2rad * lon_celld
    angle    = deg2rad * angled
    
    DO i=1,length
      
      lat_tmp = ASIN(SIN(rlat_cell) * COS(icon_resolution*i/semimaj) &
           &    + COS(rlat_cell) * SIN(icon_resolution*i/semimaj) * COS(angle))

      lon_tmp = rlon_cell + ATAN2(SIN(angle)*SIN(icon_resolution*i/semimaj) * COS(rlat_cell), &
           &   COS(icon_resolution*i/semimaj) - SIN(rlat_cell) * SIN(lat_tmp))

      dlat(i) = lat_celld - lat_tmp * rad2deg
      dlon(i) = lon_celld- lon_tmp * rad2deg
    ENDDO

  END SUBROUTINE distance_relative_to_cell

  ! determine domain center and additionaly check if latitudional extent of
  ! domain still allows a correct function of the search-algorithm for the
  ! compuation of the horizon field
  SUBROUTINE determine_domain_center(tg, center_lat, center_lon, size_radius, icon_resolution, semimaj)

    !  CAUTION:
    !  The search-algorithm assumes constant values of dlon across
    !  the entire domain, its behaviour on larger domains is not tested.
    !  If the domain extent has that size Extpar prints out a warning.

    TYPE(target_grid_def), INTENT(IN) :: tg                   !< structure with target grid description
    REAL(KIND=wp), INTENT(IN)         :: icon_resolution, &   !< hor. resolution
         &                               semimaj              !< radius earth

    INTEGER(KIND=i4), INTENT(IN)      :: size_radius          !< number of gridcells along radius
 
    REAL(KIND=wp), INTENT(OUT)        :: center_lat,&         !< lat at center of domain [deg]
      &                                  center_lon           !< lon at center of domain [deg]

    ! local variables
    REAL(KIND=wp)                     :: dlat(size_radius,2),&!< dlat [deg] 
         &                               dlon(size_radius,2),&!< dlon [deg]
         &                               dist_a, dist_b, &    !< helpers
         &                               diff_at_latmax, &    !< diff [m]
         &                               diff_at_latmin       !< diff [m]

      ! get dlat/dlon for corner with maxlat/maxlon
      CALL distance_relative_to_cell(tg%maxlat, &
           &                         tg%maxlon, &
           &                         90.0_wp,         &
           &                         size_radius,       &
           &                         dlat(:,1),         &
           &                         dlon(:,1),         &
           &                         icon_resolution,   &
           &                         semimaj)

      ! get dlat/dlon for corner with minlat/minlon
      CALL distance_relative_to_cell(tg%minlat, &
           &                         tg%minlon, &
           &                         90.0_wp,         &
           &                         size_radius,       &
           &                         dlat(:,2),         &
           &                         dlon(:,2),         &
           &                         icon_resolution,   &
           &                         semimaj)

      ! calculate distance for a "full" radius at that latitudes
      CALL haversine(semimaj, tg%maxlat, tg%maxlon, &
           &         tg%maxlat+dlat(size_radius,1), &
           &         tg%maxlon+dlon(size_radius,1), &
           &         dist_a)

      CALL haversine(semimaj, tg%minlat, tg%minlon, &
           &         tg%minlat+dlat(size_radius,2), &
           &         tg%minlon+dlon(size_radius,2), &
           &         dist_b)

      diff_at_latmax = ABS(dist_a - (size_radius * icon_resolution) )
      diff_at_latmin = ABS(dist_b - (size_radius * icon_resolution) )

      IF (diff_at_latmin > icon_resolution/2.0_wp .OR. &
          diff_at_latmax > icon_resolution/2.0_wp) THEN

        WRITE(message_text,*), 'The extent of the domain in latitudional direction is very big! ' , &
             &                 'Some gridcells may be skipped by the search-alogrithm!'
        CALL logging%warning(message_text)
      ENDIF

    ! determine domain center later used for search-algorithm
    center_lat = ( tg%maxlat + tg%minlat ) / 2.0_wp
    center_lon = ( tg%maxlon + tg%minlon ) / 2.0_wp

  END SUBROUTINE determine_domain_center

  ! calculate great-circle distance on the globe using the haversine formula
  SUBROUTINE haversine(semimaj,start_lat, start_lon, end_lat, end_lon, distance)

    REAL(KIND=wp), INTENT(IN) :: start_lat, & !< [deg]
         &                       start_lon, & !< [deg]
         &                       end_lat,   & !< [deg]
         &                       end_lon,   & !< [deg]
         &                       semimaj      !< [m]

    REAL(KIND=wp), INTENT(OUT):: distance !< [m]

    ! local variables
    REAL(KIND=wp)             :: a,c, &       !< helpers
         &                       rstart_lat, &!< [rad]
         &                       rstart_lon, &!< [rad]
         &                       rend_lat,   &!< [rad]
         &                       rend_lon     !< [rad]
    
      ! convert degree to radians
      rstart_lat=start_lat * deg2rad
      rstart_lon=start_lon * deg2rad
      rend_lat=end_lat * deg2rad
      rend_lon=end_lon * deg2rad

      ! distance between two point with haversine formula
      a = SIN(rend_lat - rstart_lat)* SIN(rend_lon - rstart_lon) + &
          &   COS(rend_lat) * COS(rstart_lat) * SIN((rend_lon-rstart_lon) / 2.0_wp) * &
          &   SIN((rend_lon-rstart_lon) / 2.0_wp)

      c = 2* ATAN2(SQRT(a),SQRT(1-a))

      distance = semimaj * c

  END SUBROUTINE haversine
  !---------------------------------------------------------------------------

END MODULE mo_lradtopo

