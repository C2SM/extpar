!+ Fortran module with routines for the COSMO target grid 
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio
!  small bug fixes accroding to Fortran compiler warnings
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with routines for the COSMO target grid 
!> definition of coordinates, input/ouptut routines
!> \author Hermann Asensio
MODULE mo_cosmo_grid

  !> kind parameters are defined in MODULE data_parameters
  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4

  USE mo_grid_structures,       ONLY: rotated_lonlat_grid, &
    &                                 target_grid_def, &
    &                                 igrid_cosmo

  USE mo_utilities_extpar,      ONLY: free_un, & 
    &                                 rlarot2rla, &
    &                                 phirot2phi

  USE mo_io_units,              ONLY: filename_max

  USE mo_target_grid_data,      ONLY: lon_geo, &
    &                                lat_geo

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_cosmo_domain_namelist, &
    &       calculate_cosmo_domain_coordinates, &
    &       get_cosmo_grid_info, &
    &       cosmo_grid, &
    &       lon_rot, &
    &       lat_rot, &
    &       res_in, &
    &       nborder, &
    &       allocate_cosmo_rc, &
            calculate_cosmo_target_coordinates


  REAL (KIND=wp), ALLOCATABLE  :: lon_rot(:), &          !< longitide coordinates of the COSMO grid in the rotated system 
    &                             lat_rot(:)          !< latitude coordinates of the COSMO grid in the rotated system
  REAL (KIND=wp)               :: res_in              !< horizontal resolution at the equator [m]
  INTEGER (KIND=i4)            :: nborder             !< number of grid points required at the border for horizon computation
  TYPE (rotated_lonlat_grid)   :: cosmo_grid  !< structure which contains the definition of the COSMO grid

  CONTAINS

  !> allocate the variables for the rotated coordinates of the COSMO grid
  SUBROUTINE allocate_cosmo_rc(ie_tot,je_tot)

    INTEGER (KIND=i4), INTENT(IN) :: ie_tot, & !< number of COSMO grid elements in zonal direction (in the rotated system) 
      &                              je_tot !< number of COSMO grid elements in meridional direction (in the rotated system) 

    INTEGER (KIND=i4)             :: errorcode
   
    ALLOCATE (lon_rot(1:ie_tot), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lon_rot', __FILE__,__LINE__)
    lon_rot = 0.0

    ALLOCATE (lat_rot(1:je_tot), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lat_rot',__FILE__,__LINE__)
    lat_rot = 0.0

  END SUBROUTINE allocate_cosmo_rc

  !> read namelist with settings for COSMO target grid
  SUBROUTINE read_cosmo_domain_namelist(input_namelist_file, &
    &                                   lrad, COSMO_grid)

    CHARACTER (len=filename_max), INTENT(IN) :: input_namelist_file !< file with input namelist with COSMO grid definition
    LOGICAL,                      INTENT(IN) :: lrad

    ! Variables for domain
    TYPE(rotated_lonlat_grid), INTENT(OUT)   :: COSMO_grid !< structure which contains the definition of the COSMO grid

     ! Variables for default values
    REAL  (KIND=wp)                          :: pollon_d, & !< longitude of the rotated north pole (in degrees, E>0)
      &                                         pollat_d, & !< latitude of the rotated north pole (in degrees, N>0)
      &                                         polgam_d, & !< angle between the north poles of the systems
      &                                         dlon_d, & !< grid point distance in zonal direction (in degrees)
      &                                         dlat_d, & !< grid point distance in meridional direction (in degrees)
      &                                         startlon_tot_d, & !< transformed longitude of lower left grid point for total domain
      &                                         startlat_tot_d, &  !< transformed latitude of lower left grid point for total domaim
      ! Variables for domain
      &                                         pollon, & !< longitude of the rotated north pole (in degrees, E>0)
      &                                         pollat, & !< latitude of the rotated north pole (in degrees, N>0)
      &                                         polgam, & !< angle between the north poles of the systems
      &                                         dlon, & !< grid point distance in zonal direction (in degrees)
      &                                         dlat, & !< grid point distance in meridional direction (in degrees)
      &                                         startlon_tot, & !< transformed longitude of lower left grid point for total domain
      &                                         startlat_tot !< transformed latitude of lower left grid point for total domain

    INTEGER  (KIND=i4)                       :: ie_tot_d, je_tot_d,ke_tot_d, &
      &                                         ie_tot, je_tot, ke_tot, &
      &                                         ierr, nuin


    REAL(KIND=wp), PARAMETER                 :: circum_earth = 40075160.0_wp, & ! Earth circumference at equator [m]
      &                                         horizon_radius = 40000.0_wp ! Radius used for horizon calculation [m]
    INTEGER(KIND=i4), PARAMETER              :: securi = 4 ! Minimum number of points required for horizon computation

    !> Define the namelist group 
    NAMELIST /lmgrid/ pollon, pollat, polgam, dlon, dlat,          &   
                 startlon_tot,startlat_tot, ie_tot, je_tot, ke_tot

    ! Comment HA: namelist input and initialization like in COSMO model, src_setup.f90
    !------------------------------------------------------------------------------
    !- Section 1: Initialize the default variables
    !------------------------------------------------------------------------------

    pollon_d       = -170.0_wp
    pollat_d       =   32.5_wp
    polgam_d       =    0.0_wp
    dlon_d         =    0.008_wp
    dlat_d         =    0.008_wp
    startlon_tot_d =   -1.252_wp
    startlat_tot_d =   -7.972_wp
    ie_tot_d    = 51
    je_tot_d    = 51
    ke_tot_d    = 0
        
    !------------------------------------------------------------------------------
    !- Section 2: Initialize variables with default
    !------------------------------------------------------------------------------

    pollon       = pollon_d
    pollat       = pollat_d
    polgam       = polgam_d
    dlon         = dlon_d
    dlat         = dlat_d
    startlon_tot = startlon_tot_d
    startlat_tot = startlat_tot_d
    ie_tot    = ie_tot_d 
    je_tot    = je_tot_d 
    ke_tot    = ke_tot_d 


    !------------------------------------------------------------------------------
    !- Section 3: Input of the namelist values
    !------------------------------------------------------------------------------

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    open(nuin,FILE=input_namelist_file, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(input_namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    read(nuin, NML=lmgrid, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist lmgrid',__FILE__, __LINE__) 
    ENDIF
    close(nuin)

    ! checks
    IF (lrad .AND. (dlon /= dlat)) THEN
      CALL logging%error('dlon must equal dlat in case of lradtopo',__FILE__,__LINE__)
    ENDIF

    ! compute number of additional grid points required in case of lradtopo
    IF (lrad) THEN
      res_in  = dlon * (circum_earth/360.0_wp) ![m]
      nborder = MAX(INT(horizon_radius/res_in,i4), securi) 
    ENDIF

    ! put values to data structure COSMO_grid
    IF (lrad) THEN
      COSMO_grid%pollon = pollon
      COSMO_grid%pollat = pollat
      COSMO_grid%polgam = polgam
      COSMO_grid%startlon_rot = startlon_tot - nborder * dlon
      COSMO_grid%startlat_rot = startlat_tot - nborder * dlat
      COSMO_grid%dlon_rot = dlon
      COSMO_grid%dlat_rot = dlat
      COSMO_grid%nlon_rot = ie_tot + 2 * nborder
      COSMO_grid%nlat_rot = je_tot + 2 * nborder
      COSMO_grid%ke_tot = ke_tot
    ELSE
      COSMO_grid%pollon = pollon
      COSMO_grid%pollat = pollat
      COSMO_grid%polgam = polgam
      COSMO_grid%startlon_rot = startlon_tot
      COSMO_grid%startlat_rot = startlat_tot
      COSMO_grid%dlon_rot = dlon
      COSMO_grid%dlat_rot = dlat
      COSMO_grid%nlon_rot = ie_tot
      COSMO_grid%nlat_rot = je_tot
      COSMO_grid%ke_tot = ke_tot
    ENDIF

  END SUBROUTINE read_cosmo_domain_namelist

  !> calculate coordinates of COSMO domain grid
  SUBROUTINE calculate_cosmo_domain_coordinates(tg,COSMO_grid)

    TYPE(target_grid_def), INTENT(IN)      :: tg !< structure with target grid description
    TYPE(rotated_lonlat_grid), INTENT(IN)  :: cosmo_grid !< structure which contains the definition of the COSMO grid

    ! local variables
    REAL (KIND=wp) :: pollat, pollon, polgam
    INTEGER :: i,j,k  ! counters

    k = 1 ! for a COSMO grid the third dimension of target arrays is 1
    pollon = COSMO_grid%pollon
    pollat = COSMO_grid%pollat
    polgam = COSMO_grid%polgam
    
    ! coordinates in the rotated system
    DO i=1, tg%ie
      lon_rot(i) = COSMO_grid%startlon_rot + (i-1) * COSMO_grid%dlon_rot
    ENDDO
    DO j=1, tg%je
      lat_rot(j) =  COSMO_grid%startlat_rot + (j-1) * COSMO_grid%dlat_rot
    ENDDO
       
    ! convert coordinates from the rotated system to the geographical system
    DO j=1, tg%je
      DO i=1, tg%ie
        lon_geo(i,j,k) = rlarot2rla(lat_rot(j),lon_rot(i),pollat,pollon,polgam)
        lat_geo(i,j,k) = phirot2phi(lat_rot(j),lon_rot(i),pollat,polgam)
      ENDDO
    ENDDO
    
  END SUBROUTINE calculate_cosmo_domain_coordinates

!> calculate coordinates of COSMO target grid
!! determin the coordintas for the rotated system and 
!! transform for each target grid element the coordinates into the geographical system  
  SUBROUTINE calculate_cosmo_target_coordinates(tg,cosmo_grid,lon_geo,lat_geo,lon_rot,lat_rot)

    TYPE(target_grid_def), INTENT(in)      :: tg              !< !< structure with target grid description
    TYPE(rotated_lonlat_grid), INTENT(in)  :: cosmo_grid !< structure which contains the definition of the COSMO grid

    REAL (KIND=wp), INTENT(inout)          :: lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
      &                                       lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
      &                                       lon_rot(:), & !< longitide coordinates of the COSMO grid in the rotated system 
                                              lat_rot(:) !< latitude coordinates of the COSMO grid in the rotated system
    ! local variables
    INTEGER :: i,j,k  ! counters

    CALL logging%info('Enter routine: calculate_cosmo_target_grid_coordinates')

    k = 1 ! for a COSMO grid the third dimension of target arrays is 1
    ! coordinates in the rotated system
    DO i=1, tg%ie
      lon_rot(i) = cosmo_grid%startlon_rot + (i-1) * cosmo_grid%dlon_rot
    ENDDO
    DO j=1, tg%je
      lat_rot(j) =  cosmo_grid%startlat_rot + (j-1) * cosmo_grid%dlat_rot
    ENDDO

    ! convert coordinates from the rotated system to the geographical system
    DO j=1, tg%je
      DO i=1, tg%ie
        lon_geo(i,j,k) = rlarot2rla(lat_rot(j),lon_rot(i), &
          &                         cosmo_grid%pollat,cosmo_grid%pollon,cosmo_grid%polgam)
        lat_geo(i,j,k) = phirot2phi(lat_rot(j),lon_rot(i), &
         &                         cosmo_grid%pollat,cosmo_grid%polgam)
      ENDDO
    ENDDO

  END SUBROUTINE calculate_cosmo_target_coordinates

!> get Information for COSMO_grid from namelist INPUT_COSMO_GRID
  SUBROUTINE get_cosmo_grid_info(input_namelist_file,tg,cosmo_grid,lrad)

   CHARACTER (len=filename_max), INTENT(IN) :: input_namelist_file !< file with input namelist with COSMO grid definition
   LOGICAL                     , INTENT(IN) :: lrad
   TYPE(target_grid_def), INTENT(OUT)       :: tg              !< !< structure with target grid description
   TYPE(rotated_lonlat_grid), INTENT(OUT)   :: cosmo_grid !< structure which contains the definition of the COSMO grid

   CALL  read_cosmo_domain_namelist(input_namelist_file,     &
     &                              lrad, cosmo_grid)

   !describe the target grid
   tg%igrid_type = igrid_cosmo ! igrid_cosmo=2 for the COSMO grid
   tg%ie = cosmo_grid%nlon_rot
   tg%je = cosmo_grid%nlat_rot
   tg%ke = 1              ! third dimension with length 1     
      
  END SUBROUTINE get_cosmo_grid_info
 
END Module mo_cosmo_grid
