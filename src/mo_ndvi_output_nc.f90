!+ Fortran module for netcdf output of NDVI data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_2         2011/03/25 Hermann Asensio
!  small bug fix: remove some USE statements for unused global variables
! V1_3         2011/04/19 Hermann Asensio
!  change netcdf output:  time variable
! V1_4         2011/04/21 Anne Roche
!  clean up
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for netcdf output of NDVI data
!> \author Hermann Asensio
MODULE mo_ndvi_output_nc

  
  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  !> data type structures form module GRID_structures
  USE mo_grid_structures, ONLY: rotated_lonlat_grid
  USE mo_grid_structures, ONLY: icosahedral_triangular_grid
  USE mo_grid_structures, ONLY: target_grid_def

  USE mo_io_utilities, ONLY: netcdf_attributes
  USE mo_io_utilities, ONLY: dim_meta_info
  USE mo_io_utilities, ONLY: netcdf_put_var
  USE mo_io_utilities, ONLY: open_new_netcdf_file
  USE mo_io_utilities, ONLY: close_netcdf_file
  USE mo_io_utilities, ONLY: netcdf_def_grid_mapping
  USE mo_io_utilities, ONLY: set_date_mm_extpar_field

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_netcdf_buffer_ndvi
  PUBLIC :: write_netcdf_cosmo_grid_ndvi
  PUBLIC :: write_netcdf_icon_grid_ndvi

  PUBLIC :: read_netcdf_buffer_ndvi

  CONTAINS

  SUBROUTINE write_netcdf_buffer_ndvi(netcdf_filename,  &
   &                                     tg,         &
   &                                     ntime, &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     ndvi_max,  &
   &                                     ndvi_field_mom,&
   &                                     ndvi_ratio_mom)


    USE mo_var_meta_data, ONLY: dim_3d_tg, &
      &                         def_dimension_info_buffer

    USE mo_var_meta_data, ONLY: lon_geo_meta, &
      &                         lat_geo_meta, &
      &                         def_com_target_fields_meta  

    USE mo_var_meta_data, ONLY: ndvi_max_meta, &
      &                         ndvi_field_mom_meta, &
      &                         ndvi_ratio_mom_meta,&
      &                         def_ndvi_meta



    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN) :: ntime !< number of times of ndvi data (12 monthly mean values)
    REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
    INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
    REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
    REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
    REAL (KIND=wp), INTENT(IN) :: ndvi_max(:,:,:) !< field for ndvi maximum
    REAL (KIND=wp), INTENT(IN) :: ndvi_field_mom(:,:,:,:) !< field for monthly mean ndvi data (12 months)
    REAL (KIND=wp), INTENT(IN) :: ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)

    ! local variables
    ! local variables
    REAL (KIND=wp),ALLOCATABLE :: time(:) !< time variable
    INTEGER (KIND=i8) :: dataDate  !< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
   INTEGER (KIND=i8) :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm

    INTEGER :: ndims  
    INTEGER :: ncid

    TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
    
    INTEGER, PARAMETER :: nglob_atts=6
    TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

    INTEGER :: errorcode !< error status variable

    INTEGER :: n !< counter

    PRINT *,'ENTER write_netcdf_buffer_ndvi'

    PRINT *,'set_global_att_ndvi'

    !-------------------------------------------------------------
    ! define global attributes
    CALL set_global_att_ndvi(global_attributes)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    PRINT *,'def_com_target_fields_meta'
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various NDVI data related variables for netcdf output
    CALL def_ndvi_meta(tg,ntime,dim_3d_tg)
    ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta
    
    ALLOCATE(time(1:ntime),STAT=errorcode)
    IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array time')
    DO n=1,ntime
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO

    ! set up dimensions for netcdf output 
    ndims = 4
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')
      
    dim_list(1) = dim_3d_tg(1) ! ie
    dim_list(2) = dim_3d_tg(2) ! je
    dim_list(3) = dim_3d_tg(3) ! ke
    dim_list(4)%dimname = 'time' ! ntime
    dim_list(4)%dimsize = ntime 

    !-----------------------------------------------------------------

    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       time=time,          &
      &                       ncid=ncid)

    ! lon
    CALL netcdf_put_var(ncid,lon_geo,lon_geo_meta,undefined)

    ! lat
    CALL netcdf_put_var(ncid,lat_geo,lat_geo_meta,undefined)

    ! ndvi_max
    CALL netcdf_put_var(ncid,ndvi_max,ndvi_max_meta,undefined)

    ! ndvi_field_mom
    CALL netcdf_put_var(ncid,ndvi_field_mom,ndvi_field_mom_meta,undefined)

    ! ndvi_ratio_mom
    CALL netcdf_put_var(ncid,ndvi_ratio_mom,ndvi_ratio_mom_meta,undefined)

    CALL close_netcdf_file(ncid)


   END SUBROUTINE write_netcdf_buffer_ndvi
   !-----------------------------------------------------------------
   !-----------------------------------------------------------------
   !-----------------------------------------------------------------



   SUBROUTINE write_netcdf_cosmo_grid_ndvi(netcdf_filename,  &
   &                                     cosmo_grid,         &
   &                                     tg,         &
   &                                     ntime, &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     ndvi_max,  &
   &                                     ndvi_field_mom,&
   &                                     ndvi_ratio_mom)

    
    USE mo_var_meta_data, ONLY: nc_grid_def_cosmo, &
    &                           set_nc_grid_def_cosmo

    USE mo_var_meta_data, ONLY: dim_rlon_cosmo, &
    &                         dim_rlat_cosmo, &
    &                         dim_2d_cosmo,   &
    &                         rlon_meta,      &
    &                         rlat_meta,      &
    &                         def_dimension_info_cosmo

    USE mo_cosmo_grid, ONLY: lon_rot, lat_rot

    USE mo_var_meta_data, ONLY: def_dimension_info_buffer
    USE mo_var_meta_data, ONLY: def_com_target_fields_meta  
                               
    USE mo_var_meta_data, ONLY: ndvi_max_meta, &
      &                         ndvi_field_mom_meta, &
      &                         ndvi_ratio_mom_meta,&
      &                         def_ndvi_meta



    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(rotated_lonlat_grid), INTENT(IN)  :: COSMO_grid      !< structure which contains the definition of the COSMO grid
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN) :: ntime !< number of times of ndvi data (12 monthly mean values)
    REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
    INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
    REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
    REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
    REAL (KIND=wp), INTENT(IN) :: ndvi_max(:,:,:) !< field for ndvi maximum
    REAL (KIND=wp), INTENT(IN) :: ndvi_field_mom(:,:,:,:) !< field for monthly mean ndvi data (12 months)
    REAL (KIND=wp), INTENT(IN) :: ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)


    ! local variables
     REAL (KIND=wp),ALLOCATABLE :: time(:) !< time variable
    INTEGER (KIND=i8) :: dataDate  !< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
   INTEGER (KIND=i8) :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm



    INTEGER :: ndims  
    INTEGER :: ncid
    INTEGER :: varid

    TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
    
    INTEGER, PARAMETER :: nglob_atts=6
    TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

    CHARACTER (len=80):: grid_mapping !< netcdf attribute grid mapping
    CHARACTER (len=80):: coordinates  !< netcdf attribute coordinates


    INTEGER :: errorcode !< error status variable

    INTEGER :: n !< counter

    PRINT *,'ENTER write_netcdf_buffer_ndvi'

    PRINT *,'set_global_att_ndvi'

    !-------------------------------------------------------------
    ! define global attributes
    CALL set_global_att_ndvi(global_attributes)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    
    !set up dimensions for COSMO grid
    CALL def_dimension_info_cosmo(cosmo_grid)
    ! dim_rlon_cosmo, dim_rlat_cosmo, dim_2d_cosmo, rlon_meta, rlat_meta

    ! set mapping parameters for netcdf
    grid_mapping="rotated_pole"
    coordinates="lon lat"
    CALL set_nc_grid_def_cosmo(cosmo_grid,grid_mapping)
    ! nc_grid_def_cosmo
    PRINT *,'def_com_target_fields_meta'
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_2d_cosmo,coordinates,grid_mapping)
    ! lon_geo_meta and lat_geo_meta


    !define meta information for various NDVI data related variables for netcdf output
    CALL def_ndvi_meta(tg,ntime,dim_2d_cosmo,coordinates,grid_mapping)
    ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta

    ALLOCATE(time(1:ntime),STAT=errorcode)
    IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array time')
    DO n=1,ntime
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO

    ! set up dimensions for netcdf output 
    ndims = 3
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')
    dim_list(1) = dim_rlon_cosmo(1) ! rlon
    dim_list(2) = dim_rlat_cosmo(1) ! rlat
    dim_list(3)%dimname = 'time'
    dim_list(3)%dimsize = ntime
    
   !-----------------------------------------------------------------
    PRINT *,' CALL open_new_netcdf_file'
    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       time=time,          &
      &                       ncid=ncid)
    !-----------------------------------------------------------------

    ! rlon
    !HA debug
    PRINT *,'HA debug: put rlon to netcdf'
    CALL netcdf_put_var(ncid,lon_rot(1:cosmo_grid%nlon_rot),rlon_meta,undefined)

    PRINT *,'HA debug: put rlat to netcdf'
    ! rlat
    CALL netcdf_put_var(ncid,lat_rot(1:cosmo_grid%nlat_rot),rlat_meta,undefined)

    ! ndvi_max
    CALL netcdf_put_var(ncid,ndvi_max(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
      &                 ndvi_max_meta,undefined)  

    ! ndvi_field_mom
    CALL netcdf_put_var(ncid,&
                       & ndvi_field_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntime), &
                       & ndvi_field_mom_meta, &
                       & undefined)

    ! ndvi_ratio_mom
    CALL netcdf_put_var(ncid,&
                       & ndvi_ratio_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntime), &
                       & ndvi_ratio_mom_meta, &
                       & undefined)

    !-----------------------------------------------------------------
    CALL netcdf_def_grid_mapping(ncid, nc_grid_def_cosmo, varid)

    CALL close_netcdf_file(ncid)


   END SUBROUTINE write_netcdf_cosmo_grid_ndvi
   !-----------------------------------------------------------------
   !-----------------------------------------------------------------
   !-----------------------------------------------------------------


   SUBROUTINE write_netcdf_icon_grid_ndvi(netcdf_filename,  &
   &                                     icon_grid,         &
   &                                     tg,         &
   &                                     ntime, &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     ndvi_max,  &
   &                                     ndvi_field_mom,&
   &                                     ndvi_ratio_mom)


    USE mo_var_meta_data, ONLY:  dim_icon, &
     &                          def_dimension_info_icon

    USE mo_var_meta_data, ONLY: set_nc_grid_def_icon

    USE mo_var_meta_data, ONLY: def_dimension_info_buffer
                               
    USE mo_var_meta_data, ONLY: lon_geo_meta, &
      &                         lat_geo_meta, &
      &                         def_com_target_fields_meta  

    USE mo_var_meta_data, ONLY: ndvi_max_meta, &
      &                         ndvi_field_mom_meta, &
      &                         ndvi_ratio_mom_meta,&
      &                         def_ndvi_meta



    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(icosahedral_triangular_grid), INTENT(IN)  :: icon_grid      !< structure which contains the definition of the ICON grid
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN) :: ntime !< number of times of ndvi data (12 monthly mean values)
    REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
    INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
    REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
    REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
    REAL (KIND=wp), INTENT(IN) :: ndvi_max(:,:,:) !< field for ndvi maximum
    REAL (KIND=wp), INTENT(IN) :: ndvi_field_mom(:,:,:,:) !< field for monthly mean ndvi data (12 months)
    REAL (KIND=wp), INTENT(IN) :: ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)


    ! local variables
     REAL (KIND=wp),ALLOCATABLE :: time(:) !< time variable
    INTEGER (KIND=i8) :: dataDate  !< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
   INTEGER (KIND=i8) :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm



    INTEGER :: ndims 
    INTEGER :: ncid

    TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
    TYPE(dim_meta_info), TARGET :: dim_1d_icon(1:1)
    
    INTEGER, PARAMETER :: nglob_atts=6
    TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

    CHARACTER (len=80):: grid_mapping !< netcdf attribute grid mapping

    INTEGER :: errorcode !< error status variable

    INTEGER :: n !< counter

    PRINT *,'ENTER write_netcdf_icon_grid_ndvi'

    PRINT *,'set_global_att_ndvi'

    !-------------------------------------------------------------
    ! define global attributes
    CALL set_global_att_ndvi(global_attributes)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    

    !set up dimensions for ICON grid
    CALL def_dimension_info_icon(icon_grid)
    ! dim_icon
    dim_1d_icon = dim_icon(1) ! cell

    
    ! set mapping parameters for netcdf
    grid_mapping="lon_lat_on_sphere"
    CALL set_nc_grid_def_icon(grid_mapping)
    ! nc_grid_def_icon
    PRINT *,'def_soil_meta'

    
    PRINT *,'def_com_target_fields_meta'
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_1d_icon)
    ! lon_geo_meta and lat_geo_meta



    !define meta information for various NDVI data related variables for netcdf output
    CALL def_ndvi_meta(tg,ntime,dim_1d_icon)
    ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta

    ALLOCATE(time(1:ntime),STAT=errorcode)
    IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array time')
    DO n=1,ntime
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO


    ! set up dimensions for netcdf output 
    ndims = 2
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')
    dim_list(1) =  dim_icon(1) ! cell
    dim_list(2)%dimname = 'time'
    dim_list(2)%dimsize = ntime

     !-----------------------------------------------------------------
    PRINT *,' CALL open_new_netcdf_file'
    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
        &                       dim_list=dim_list,                  &
        &                       global_attributes=global_attributes, &
        &                       time=time,          &
        &                       ncid=ncid)
    !-----------------------------------------------------------------

    ! lon
    CALL netcdf_put_var(ncid,lon_geo(1:icon_grid%ncell,1,1),lon_geo_meta,undefined)

    ! lat
    CALL netcdf_put_var(ncid,lat_geo(1:icon_grid%ncell,1,1),lat_geo_meta,undefined)

    ! ndvi_max
    CALL netcdf_put_var(ncid,ndvi_max(1:icon_grid%ncell,1,1),ndvi_max_meta,undefined)

    ! ndvi_field_mom
    CALL netcdf_put_var(ncid,ndvi_field_mom(1:icon_grid%ncell,1,1,1:ntime), &
       &                 ndvi_field_mom_meta, undefined)

    ! ndvi_ratio_mom
    CALL netcdf_put_var(ncid,ndvi_ratio_mom(1:icon_grid%ncell,1,1,1:ntime), &
       &                 ndvi_ratio_mom_meta, undefined)

    CALL close_netcdf_file(ncid)

   END SUBROUTINE write_netcdf_icon_grid_ndvi

   !----------------------------------------------------------------------- 
   !-----------------------------------------------------------------
   !-----------------------------------------------------------------------
  !> set global attributes for netcdf with NDVI data
  SUBROUTINE set_global_att_ndvi(global_attributes)
    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:6)

    !local variables
    CHARACTER(len=10) :: ydate
    CHARACTER(len=10) :: ytime
    CHARACTER(len=2)  :: cc
    CHARACTER(len=2)  :: yy
    CHARACTER(len=2)  :: mm
    CHARACTER(len=2)  :: dd
    CHARACTER(len=2)  :: hh
    CHARACTER(len=2)  :: minute

    ! define global attributes
    
    global_attributes(1)%attname = 'title'
    global_attributes(1)%attributetext='NDVI data '
    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='Deutscher Wetterdienst'

    global_attributes(3)%attname = 'source'
    global_attributes(3)%attributetext='NASA/GSFC SeaWiFS'

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' ndvi_to_buffer'

    global_attributes(5)%attname = 'references'
    global_attributes(5)%attributetext='http://oceancolor.gsfc.nasa.gov/PRODUCTS/'

    global_attributes(6)%attname = 'comment'
    global_attributes(6)%attributetext=''

  END SUBROUTINE set_global_att_ndvi
  !-----------------------------------------------------------------------

  SUBROUTINE read_netcdf_buffer_ndvi(netcdf_filename,  &
   &                                     tg,         &
   &                                     ntime, &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     ndvi_max,  &
   &                                     ndvi_field_mom,&
   &                                     ndvi_ratio_mom)

    USE mo_var_meta_data, ONLY: dim_3d_tg, &
      &                         def_dimension_info_buffer

    USE mo_var_meta_data, ONLY: def_com_target_fields_meta
    USE mo_var_meta_data, ONLY: ndvi_max_meta, &
      &                         ndvi_field_mom_meta, &
      &                         ndvi_ratio_mom_meta,&
      &                         def_ndvi_meta

    USE mo_io_utilities, ONLY: netcdf_get_var

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(INOUT) :: ntime !< number of times of ndvi data (12 monthly mean values)
    REAL(KIND=wp), INTENT(OUT)          :: undefined       !< value to indicate undefined grid elements 
    INTEGER, INTENT(OUT)                :: undef_int       !< value to indicate undefined grid elements
    REAL (KIND=wp), INTENT(OUT) :: ndvi_max(:,:,:) !< field for ndvi maximum
    REAL (KIND=wp), INTENT(OUT) :: ndvi_field_mom(:,:,:,:) !< field for monthly mean ndvi data (12 months)
    REAL (KIND=wp), INTENT(OUT) :: ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)

    ! local variables
    INTEGER, PARAMETER :: nglob_atts=6

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    PRINT *,'def_com_target_fields_meta'
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various NDVI data related variables for netcdf output
    CALL def_ndvi_meta(tg,ntime,dim_3d_tg)
    ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta

    PRINT *,'CALL read netcdf data NDVI'

    CALL netcdf_get_var(TRIM(netcdf_filename),ndvi_max_meta,ndvi_max)
    PRINT *,'ndvi_max read'

    CALL netcdf_get_var(TRIM(netcdf_filename),ndvi_field_mom_meta,ndvi_field_mom)
    PRINT *,'ndvi_field_mom read'

    CALL netcdf_get_var(TRIM(netcdf_filename),ndvi_ratio_mom_meta,ndvi_ratio_mom)
    PRINT *,'ndvi_ratio_mom read'



   END SUBROUTINE read_netcdf_buffer_ndvi
   !-----------------------------------------------------------------






END Module mo_ndvi_output_nc

