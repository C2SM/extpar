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
  
  USE mo_kind,             ONLY: wp, i8, i4
  USE mo_utilities_extpar, ONLY: abort_extpar, phi2phirot, rla2rlarot, free_un
  USE mo_cosmo_grid,       ONLY: nborder, res_in, cosmo_grid, lon_rot, lat_rot
  USE mo_grid_structures,  ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_namelists_extpar_lradtopo
  PUBLIC :: compute_lradtopo
  PUBLIC :: deghor, pi, rad2deg, deg2rad

  REAL(wp)     :: deghor               !< number of degrees per sector [deg]

  REAL(wp), PARAMETER :: &  
       & pi             = 3.14159265359_wp,          & !< pi
       & rad2deg        = 180._wp/pi,                & !< radians to degrees
       & deg2rad        = pi/180._wp                   !< degrees to radians

CONTAINS

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  !> subroutine to read namelist for lradtopo in EXTPAR
  SUBROUTINE read_namelists_extpar_lradtopo(namelist_file,        &
       &                                    lradtopo,             &
       &                                    nhori)



    CHARACTER (len=*), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings
    !  lradtopo

    LOGICAL,          INTENT(OUT) :: lradtopo                 !< parameters for lradtopo to be computed? (TRUE/FALSE)
    INTEGER(KIND=i4), INTENT(OUT) :: nhori                    !< number of sectors for the horizon computation 

    !> local variables
    INTEGER           :: nuin     !< unit number
    INTEGER (KIND=i4) :: ierr     !< error flag

    !> variables for default values
    LOGICAL           :: lradtopo_d

    INTEGER(KIND=i4)  :: nhori_d

    !> define the namelist group
    NAMELIST /radtopo/ lradtopo, nhori

    !> initialization
    ierr            = 0

    !> default values definition
    lradtopo_d      = .FALSE.
    nhori_d         = 24

    !> default values attribution
    lradtopo        = lradtopo_d 
    nhori           = nhori_d

    !> read namelist  
    nuin = free_un()  ! function free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL abort_extpar('read_namelists_extpar_lradtopo: cannot open file')
    ENDIF
    READ(nuin, NML=radtopo, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL abort_extpar('read_namelists_extpar_lradtopo: cannot read file')
    ENDIF
    CLOSE(nuin)

    !> check values for consistency
    IF ( nhori > 24_i8 ) THEN
      PRINT *,' Warning *** nhori larger than 24 is not recommended  *** '
    ENDIF

    !> print values
    PRINT *, 'lradtopo:', lradtopo
    PRINT *, 'nhori: '  , nhori

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
    INTEGER           :: errorcode
    INTEGER (KIND=i8) :: i, j
    INTEGER (KIND=i8) :: nsec                !< number of gridpoints per sector (in both x and y dir.)
    REAL(KIND=wp)     :: rlon_np, rlat_np, & !< location of true North Pole in rot. coord.
                         rdx, rdy,         & !< distance from the sector center in x and y dir.
                         dhdx, dhdy,       & !< slope in x and y directions resp.
                         coslat              !< cosine of rotated latitude

    REAL(KIND=wp), ALLOCATABLE :: zhh(:,:),                   & !< local mean height
         &                        zslope_asp(:,:),            & !< local slope aspect
         &                        zslope_ang(:,:),            & !< local slope angle
         &                        zhorizon  (:,:,:),          & !< local horizon
         &                        zskyview  (:,:),            & !< local skyview
         &                        slope_x(:,:), slope_y(:,:), & !< slope in x and y directions resp.
         &                        h_hres(:,:),                & !< corrected height above see level of target grid
         &                        h_corr(:,:),                & !< height correction (Earth's curvature)
         &                        dist(:,:),                  & !< distance from sector center (local search grid)
         &                        aberr(:,:)                    !< angle between geographic meridian and grid meridian

    !> parameters
    REAL(KIND=wp), PARAMETER :: semimaj = 6378137.0         !< semimajor radius WGS 84

    !---------------------------------------------------------------------------

    errorcode = 0

    nsec   = 1 + 2 * nborder
    deghor = 360.0_wp / nhori

    !> allocations and initializations
    ALLOCATE( h_hres(tg%ie,tg%je), STAT=errorcode )
    IF ( errorcode /= 0 ) CALL abort_extpar( 'Cant allocate the array h_hres' )

    ALLOCATE( zhh(tg%ie,tg%je), STAT=errorcode )
    IF ( errorcode /= 0 ) CALL abort_extpar( 'Cant allocate the array zhh' )

    ALLOCATE( h_corr(nsec,nsec), dist(nsec,nsec), STAT=errorcode )
    IF ( errorcode /= 0 ) CALL abort_extpar( 'Cant allocate the arrays h_corr and dist' )

    ALLOCATE( zslope_asp(tg%ie-2*nborder,tg%je-2*nborder),        &
         &    zslope_ang(tg%ie-2*nborder,tg%je-2*nborder),        &
         &    zhorizon  (tg%ie-2*nborder,tg%je-2*nborder,nhori),  &
         &    zskyview  (tg%ie-2*nborder,tg%je-2*nborder),        &
         &    slope_x   (tg%ie-2*nborder,tg%je-2*nborder),        &
         &    slope_y   (tg%ie-2*nborder,tg%je-2*nborder),        &
         &    aberr     (tg%ie-2*nborder,tg%je-2*nborder),        &
         &    STAT = errorcode )
    IF ( errorcode /= 0 ) CALL abort_extpar( 'Cant allocate the lradtopo arrays' )

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
        !        IF (ldebug) PRINT*, ' treating point (',i,',',j,')'

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
    IF ( errorcode /= 0 ) CALL abort_extpar( 'Cant deallocate the array h_hres' )

    DEALLOCATE( zhh, STAT=errorcode )
    IF ( errorcode /= 0 ) CALL abort_extpar( 'Cant deallocate the array zhh' )

    DEALLOCATE( h_corr, dist, STAT=errorcode )
    IF ( errorcode /= 0 ) CALL abort_extpar( 'Cant deallocate the arrays h_corr and dist' )

    DEALLOCATE(zslope_asp,                         &
         &     zslope_ang,                         &
         &     zhorizon  ,                         &
         &     zskyview  ,                         &
         &     slope_x   ,                         &
         &     slope_y   ,                         &
         &     aberr     ,                         &
         STAT = errorcode )
    IF ( errorcode /= 0 ) CALL abort_extpar( 'Cant deallocate the lradtopo arrays' )

  END SUBROUTINE compute_lradtopo

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  !> subroutine to compute the horizon
  SUBROUTINE comp_horiz( h_hres, dhres, dhdx, dhdy, nhori, nx, ny, hor, rot_ang )

    ! Compute horizon angles (considering self-shading) 
    ! inside the domain given by h_hres
    ! Matteo Buzzi, 2009, Meteoswiss
    ! Daniel Luethi, 2012, IAC
    ! AnneRoches, 2013, C2SM

    !> arguments
    INTEGER(KIND=i8), INTENT(IN) :: nx, ny     ! gridpoints number in the x and y dir
    INTEGER(KIND=i4), INTENT(IN) :: nhori       !  number of sectors
    REAL   (KIND=wp), INTENT(IN) :: h_hres(:,:) ! topography [m]
    REAL   (KIND=wp), INTENT(IN) :: dhres,                 & ! resolution [m] 
                                    dhdx , dhdy,           & ! slope in x and y dir resp. [-]
                                    rot_ang                  ! rotation angle (=angle between grid meridian and geo. meridian) [deg]
    REAL   (KIND=wp), INTENT(OUT):: hor(:)               ! horizon angle [deg]

    !> local variables
    INTEGER(KIND=i8) :: x0, y0,    & ! center of the sector in x and y dir resp.
                        i,  j , k, & ! loop indices
                        ngp,       & ! number of gridpoints at each side of the sector center
                        nzstat       ! allocation status

    REAL(KIND=wp)    :: tmp,       & ! temporary horizon [deg]
                        zshift,    & ! angle zshift [deg]
                        azi,       & ! angle between the geo North and the sector start  [deg]
                        ssa,       & ! self-shading [deg]
                        ho           ! horizon (maximum elevation angle in a portion of the sector)
                                     ! without self-shading [deg]

    INTEGER(KIND=i8), ALLOCATABLE :: xcart(:), ycart(:)

    REAL(KIND=wp)   , ALLOCATABLE :: pcr(:,:),    & ! polar coordinates (index1: angle, index2: radius)
         &                           rcr(:,:),    & ! rectangular coord.(index1: x    , index2: y     )
         &                           dh(:),dn(:), & ! distance between the sector center and the other points of 
                                                    ! interest in the sector in [m] in the vertical and in the horizontal resp.
         &                           rates(:)       ! slopes between the sector center and the other points of
                                                    ! interest in the sector

    INTEGER(KIND=i8), PARAMETER :: niter = 10 ! number of iterations

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    !> preparations
    x0     = nx/2 +1
    y0     = ny/2 +1
    ngp    = nx/2

    !> allocations
    ALLOCATE( pcr(2,ngp), rcr(2,ngp), xcart(ngp), ycart(ngp), dh(ngp), dn(ngp), rates(ngp), STAT=nzstat )
    IF ( nzstat /= 0_i8 ) STOP

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
    IF ( nzstat /= 0_i8 ) STOP


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
                                                   slope_asp    ! slope aspect [deg]

    INTEGER(KIND=i4), INTENT(IN)                :: nhori        ! number of sectors 
    REAL(KIND=wp), INTENT(IN)                   :: horizon(:)   ! horizon [deg]   
    REAL(KIND=wp), INTENT(OUT)                  :: skyview      ! skyview factor [-]

    !> local variables
    INTEGER(KIND=i8)               :: i, i0, k, ip

    REAL(KIND=wp)                  :: rslope_ang,             & ! slope angle [rad]
         &                            rdeghor,                & ! number of radians per sector [rad]
         &                            cosa, cosah, sinah,     & ! cosine, half cosine, half sine of the slope angle resp. [-] 
         &                            sf,                     & !
         &                            azi,                    & ! azimut of the trailing trace of the slope on the horizontal
                                                                ! plane [deg]
         &                            slg                       !

    REAL(KIND=wp), DIMENSION(nhori+1):: ang_start,            & ! angle of the sector start  [rad]
         &                              hop,                  & !
         &                              sl,                   & !
         &                              hpl

    REAL(KIND=wp), DIMENSION(2*nhori+1):: drhorizon             ! duplicated horizon [rad]


    REAL(KIND=wp), PARAMETER :: zepsilon = 2.0_wp * pi * 0.01_wp ! security parameter 

    LOGICAL:: ldebug=.FALSE.

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    !> preparations

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
              PRINT*, '!! SLG TOO HIGH    !!','slope_ang: ', slope_ang, ' slope_asp: ', slope_asp, ' horizon: ', horizon
              PRINT*, 'slg:  ', slg, ' sl(ip): ', sl(ip), ' sl(i): ', sl(i), ' hop(i): ', hop(i), ' hpl(i): ', hpl(i)
              PRINT*,  ' hpl(ip): ', hpl(ip), ' i:  ',i, ' ip: ', ip
              STOP
            ENDIF
            ! SLG too low
            IF ((slg < sl(i)) .AND. (slg < sl(i)-zepsilon)) THEN 
              PRINT*, '!! SLG TOO LOW !!','slope_ang: ', slope_ang, ' slope_asp: ', slope_asp, ' horizon: ', horizon
              PRINT*, 'slg:  ', slg, ' sl(ip): ', sl(ip), ' sl(i): ', sl(i), ' hop(i): ', hop(i), ' hpl(i): ', hpl(i)
              PRINT*,  ' hpl(ip): ', hpl(ip), ' i:  ',i, ' ip: ', ip
              STOP
            ENDIF
            ! SLG in tolerance zone -> clipped to lower limit
            IF ( (slg > sl(ip)) .AND. (slg < sl(ip)+zepsilon) ) THEN 
              IF ( ldebug ) THEN
                PRINT*, 'SLG clipped to lower limit','slg before: ', slg, ' slg clipped: ', sl(ip)
              ENDIF
              slg = sl(ip)
            ENDIF
            ! SLG in tolerance zone -> clipped to upper limit
            IF ( (slg < sl(i)) .AND. (slg > sl(i)-zepsilon) ) THEN
              IF ( ldebug ) THEN 
                PRINT*, '!! SLG clipped to upper limit !!','slg before: ', slg, ' slg clipped: ', sl(i)
              ENDIF
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
              PRINT*, '!! slg TOO HIGH    !!','slope_ang: ', slope_ang, ' slope_asp: ', slope_asp, ' horizon: ', horizon
              PRINT*, 'slg:  ', slg, ' sl(ip): ', sl(ip), ' sl(i): ', sl(i), ' hop(i): ', hop(i), ' hpl(i): ', hpl(i)
              PRINT*,  ' hpl(ip): ', hpl(ip), ' i:  ',i, ' ip: ', ip
              STOP
            ENDIF
            ! slg too low
            IF ( (slg < sl(i)) .AND. (slg < sl(i)-zepsilon) ) THEN 
              PRINT*, '!! slg TOO LOW !!','slope_ang: ', slope_ang, ' slope_asp: ', slope_asp, ' horizon: ', horizon
              PRINT*, 'slg:  ', slg, ' sl(ip): ', sl(ip), ' sl(i): ', sl(i), ' hop(i): ', hop(i), ' hpl(i): ', hpl(i)
              PRINT*,  ' hpl(ip): ', hpl(ip), ' i:  ',i, ' ip: ', ip
              STOP
            ENDIF
            ! slg in tolerance zone -> clipped to lower limit
            IF ( (slg > sl(ip)) .AND. (slg < sl(ip)+zepsilon) ) THEN 
              IF ( ldebug ) THEN
                PRINT*, 'slg clipped to lower limit','slg before: ', slg, ' slg clipped: ', sl(ip)
              ENDIF
              slg = sl(ip)
            ENDIF
            ! slg in tolerance zone -> clipped to upper limit
            IF ( (slg < sl(i)) .AND. (slg > sl(i)-zepsilon) ) THEN 
              IF ( ldebug ) THEN
                PRINT*, '!! slg clipped to upper limit !!','slg before: ', slg, ' slg clipped: ', sl(i)
              ENDIF
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

  END SUBROUTINE compute_skyview

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  !> function to ???
  FUNCTION atan22( psi, cosalfa )

    REAL (KIND=wp),INTENT(IN) :: psi, cosalfa

    REAL (KIND=wp)            :: atan22
    REAL (KIND=wp)            :: pih, pi3h, x

    IF( cosalfa == 0.0_wp ) THEN 
      PRINT*, 'WRONG argument in atan22 -> STOP'
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

  !---------------------------------------------------------------------------

END MODULE mo_lradtopo

