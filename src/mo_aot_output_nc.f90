!+  Fortran module for aerosol input output routines
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_3         2011/04/19 Hermann Asensio
!  change netcdf output:  time variable
! V1_4         2011/04/21 Hermann Asensio
!  clean up
! V2_0         2013/08/18 Daniel Luethi
!  added support for alternative aerosol climatologies AEROCOM and MNACC-II
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for aerosol input output routines
!> \author Hermann Asensio
MODULE mo_aot_output_nc

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  !> data type structures form module GRID_structures
  USE mo_grid_structures, ONLY: reg_lonlat_grid
  USE mo_grid_structures, ONLY: rotated_lonlat_grid
  USE mo_grid_structures, ONLY: icosahedral_triangular_grid
  USE mo_grid_structures, ONLY: target_grid_def

  USE mo_io_utilities, ONLY: var_meta_info

  USE mo_io_utilities, ONLY: check_netcdf
  USE mo_io_utilities, ONLY: netcdf_attributes

  USE mo_io_utilities, ONLY: netcdf_put_var
  USE mo_io_utilities, ONLY: open_new_netcdf_file
  USE mo_io_utilities, ONLY: close_netcdf_file
  USE mo_io_utilities, ONLY: netcdf_def_grid_mapping

  USE mo_io_utilities, ONLY: get_date_const_field
  USE mo_io_utilities, ONLY: set_date_mm_extpar_field

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_netcdf_buffer_aot
  PUBLIC :: write_netcdf_cosmo_grid_aot
  PUBLIC :: write_netcdf_icon_grid_aot

  PUBLIC :: read_netcdf_buffer_aot


  CONTAINS

  !> create a netcdf file for the AOT data in the buffer
  SUBROUTINE write_netcdf_buffer_aot(netcdf_filename,    &
   &                                     tg,             &
   &                                     undefined,      &
   &                                     undef_int,      &
   &                                     lon_geo,        &
   &                                     lat_geo,        &
   &                                     ntype,          &
   &                                     ntime,          &
   &                                     aot_tg,         &
   &                                     iaot_type)


  USE mo_io_utilities, ONLY: netcdf_grid_mapping, &
    &                        netcdf_char_attributes, &
    &                        netcdf_real_attributes

  USE mo_io_utilities, ONLY: dim_meta_info

  USE mo_io_utilities, ONLY: vartype_int 
  USE mo_io_utilities, ONLY: vartype_real
  USE mo_io_utilities, ONLY: vartype_char

  USE mo_grid_structures, ONLY: target_grid_def

  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer

  USE mo_var_meta_data, ONLY: dim_aot_tg, &
    &                         aot_tg_meta, &
    &                         def_aot_tg_meta

  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  INTEGER (KIND=i8), INTENT(IN) :: ntype !< number of types of aerosols
  INTEGER (KIND=i8), INTENT(IN) :: ntime !< number of times
  INTEGER (KIND=i4), INTENT(IN) :: iaot_type !< ID of aeorosol raw data

  REAL (KIND=wp), INTENT(IN)  :: aot_tg(:,:,:,:,:) !< aerosol optical thickness, aot_tg(ie,je,ke,ntype,ntime)

  ! local variables
  REAL (KIND=wp),ALLOCATABLE :: time(:) !< time variable
  INTEGER (KIND=i8) :: dataDate  !< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
  INTEGER (KIND=i8) :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm

  INTEGER, PARAMETER :: nglob_atts=5
  TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

  INTEGER :: ndims  
  INTEGER :: ncid
  TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file

  TYPE(dim_meta_info), TARGET :: dim_3d_buffer(1:3)
  TYPE(dim_meta_info), TARGET :: dim_5d_buffer(1:5)


  INTEGER :: errorcode !< error status variable
  CHARACTER (len=80) :: attname
  CHARACTER (len=255) :: attributetext

  INTEGER :: n !< counter

  !-------------------------------------------------------------
  ! define global attributes
  IF (iaot_type == 1 ) THEN
     CALL set_global_att_aot(global_attributes)
  ELSEIF(iaot_type == 2 ) THEN
     CALL set_global_att_aot_aero(global_attributes)
  ELSEIF(iaot_type == 3 ) THEN
     CALL set_global_att_aot_macc(global_attributes)
  ELSE
     PRINT *, 'UNKNOWN AOT DATA OPTION: '
     STOP 1
  ENDIF

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg
  
  ! define dimensions and meta information for variable aot_tg for netcdf output
  CALL def_aot_tg_meta(tg,ntime,ntype,dim_3d_tg)
  ! dim_aot_tg and aot_tg_meta
  
  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta

  ALLOCATE(time(1:ntime),STAT=errorcode)
    IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array time')
    DO n=1,ntime
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO

  ! set up dimensions for netcdf output 
  ndims = 5
  ALLOCATE(dim_list(1:ndims),STAT=errorcode)
  IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')

  dim_list(1) = dim_3d_tg(1) ! ie
  dim_list(2) = dim_3d_tg(2) ! je
  dim_list(3) = dim_3d_tg(3) ! ke
  dim_list(4)%dimname = 'ntype'
  dim_list(4)%dimsize = ntype
  dim_list(5)%dimname = 'time'
  dim_list(5)%dimsize = ntime


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

  ! aot_tg
  CALL netcdf_put_var(ncid,aot_tg,aot_tg_meta,undefined)

  CALL close_netcdf_file(ncid)



  END SUBROUTINE write_netcdf_buffer_aot
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------


  !> create a netcdf file for the AOT data in the COSMO grid
  SUBROUTINE write_netcdf_cosmo_grid_aot(netcdf_filename,  &
   &                                     cosmo_grid,       &
   &                                     tg,               &
   &                                     undefined,        &
   &                                     undef_int,        &
   &                                     lon_geo,          &
   &                                     lat_geo,          &
   &                                     ntype,            &
   &                                     ntime,            &
   &                                     aot_tg,           &
   &                                     iaot_type)


  USE mo_io_utilities, ONLY: netcdf_grid_mapping, &
    &                        netcdf_char_attributes, &
    &                        netcdf_real_attributes, &
    &                        netcdf_def_grid_mapping

  USE mo_io_utilities, ONLY: dim_meta_info

  USE mo_io_utilities, ONLY: vartype_int 
  USE mo_io_utilities, ONLY: vartype_real
  USE mo_io_utilities, ONLY: vartype_char

  USE mo_grid_structures, ONLY: rotated_lonlat_grid
  USE mo_grid_structures, ONLY: igrid_cosmo

  USE mo_cosmo_grid, ONLY: lon_rot, lat_rot

  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer

  USE mo_var_meta_data, ONLY: dim_aot_tg, &
    &                         aot_tg_meta, &
    &                         def_aot_tg_meta

  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  

  USE mo_var_meta_data, ONLY: dim_rlon_cosmo, &
    &                         dim_rlat_cosmo, &
    &                         dim_2d_cosmo,   &
    &                         rlon_meta,      &
    &                         rlat_meta,      &
    &                         def_dimension_info_cosmo


  USE mo_var_meta_data, ONLY: nc_grid_def_cosmo, &
    &                         set_nc_grid_def_cosmo


  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(rotated_lonlat_grid), INTENT(IN) :: cosmo_grid !< structure which contains the definition of the COSMO grid
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  INTEGER (KIND=i8), INTENT(IN) :: ntype !< number of types of aerosols
  INTEGER (KIND=i8), INTENT(IN) :: ntime !< number of times
  INTEGER (KIND=i4), INTENT(IN) :: iaot_type !< ID of aeorosol raw data

  REAL (KIND=wp), INTENT(IN)  :: aot_tg(:,:,:,:,:) !< aerosol optical thickness, aot_tg(ie,je,ke,ntype,ntime)


  ! local variables
  REAL (KIND=wp),ALLOCATABLE :: time(:) !< time variable
  INTEGER (KIND=i8) :: dataDate  !< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
  INTEGER (KIND=i8) :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm


  INTEGER, PARAMETER :: nglob_atts=5
  TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

  TYPE(netcdf_grid_mapping) :: nc_grid_def !< mapping parameters for netcdf


  INTEGER :: ndims  
  INTEGER :: ncid
  INTEGER :: varid
  TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
  TYPE(dim_meta_info), TARGET :: dim_4d_buffer(1:4)
  TYPE(dim_meta_info), TARGET :: dim_2d_buffer(1:2)


  INTEGER :: n_1d_real = 0 !< number of 1D real variables
  INTEGER :: n_2d_real = 0 !< number of 2D real variables
  INTEGER :: n_4d_real = 0 !< number of 4D real variables

  INTEGER :: errorcode !< error status variable

  CHARACTER (len=80):: grid_mapping !< netcdf attribute grid mapping
  CHARACTER (len=80):: coordinates  !< netcdf attribute coordinates

  CHARACTER (len=1), PARAMETER :: c_undef = "-" !< default charcter for undefined string
  CHARACTER (len=80) :: attname
  CHARACTER (len=255) :: attributetext

  INTEGER :: n !< counter


  !-------------------------------------------------------------
  ! define global attributes
  IF (iaot_type == 1 ) THEN
     CALL set_global_att_aot(global_attributes)
  ELSEIF(iaot_type == 2 ) THEN
     CALL set_global_att_aot_aero(global_attributes)
  ELSEIF(iaot_type == 3 ) THEN
     CALL set_global_att_aot_macc(global_attributes)
  ELSE
     PRINT *, 'UNKNOWN AOT DATA OPTION: '
     STOP 1
  ENDIF

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

  !set up dimensions for aot_tg
  CALL def_aot_tg_meta(tg,ntime,ntype,dim_2d_cosmo,coordinates,grid_mapping)
  ! dim_aot_tg and aot_tg_meta

  CALL def_com_target_fields_meta(dim_2d_cosmo,coordinates,grid_mapping)
  ! lon_geo_meta and lat_geo_meta

  
  ALLOCATE(time(1:ntime),STAT=errorcode)
    IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array time')
    DO n=1,ntime
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO



  !set up dimensions for buffer
  ndims = 4
  ALLOCATE(dim_list(1:ndims),STAT=errorcode)
  IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')

  dim_list(1) = dim_rlon_cosmo(1) ! rlon
  dim_list(2) = dim_rlat_cosmo(1) ! rlat
  dim_list(3)%dimname = 'ntype'
  dim_list(3)%dimsize = ntype
  dim_list(4)%dimname = 'time' 
  dim_list(4)%dimsize = ntime 


  dim_2d_buffer(1:2) = dim_list(1:2)
  dim_4d_buffer(1:4) = dim_list(1:4)

  CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       time=time,                           &
      &                       ncid=ncid)
    !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  ! start with real 1d variables
       
  ! rlon
  CALL netcdf_put_var(ncid,lon_rot(1:cosmo_grid%nlon_rot),rlon_meta,undefined)

  ! rlat
  CALL netcdf_put_var(ncid,lat_rot(1:cosmo_grid%nlat_rot),rlat_meta,undefined)

  !-----------------------------------------------------------------
  ! lon
  CALL netcdf_put_var(ncid,lon_geo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 lon_geo_meta,undefined)
  ! lat
  CALL netcdf_put_var(ncid,lat_geo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 lat_geo_meta,undefined)

  ! aot_tg
   CALL netcdf_put_var(ncid,&
                       &  aot_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntype,1:ntime), &
                       & aot_tg_meta, &
                       & undefined)
  !-----------------------------------------------------------------
  CALL netcdf_def_grid_mapping(ncid, nc_grid_def_cosmo, varid)

  CALL close_netcdf_file(ncid)

  END SUBROUTINE write_netcdf_cosmo_grid_aot
  !----------------------------------------------------------------------------

  !> create a netcdf file for the AOT data in the ICON grid
  SUBROUTINE write_netcdf_icon_grid_aot(netcdf_filename,  &
   &                                     icon_grid,       &
   &                                     tg,              &
   &                                     undefined,       &
   &                                     undef_int,       &
   &                                     lon_geo,         &
   &                                     lat_geo,         &
   &                                     ntype,           &
   &                                     ntime,           &
   &                                     aot_tg,          &
   &                                     iaot_type)

  USE mo_io_utilities, ONLY: netcdf_grid_mapping, &
    &                        netcdf_char_attributes, &
    &                        netcdf_real_attributes, &
    &                        netcdf_def_grid_mapping

  USE mo_io_utilities, ONLY: dim_meta_info

  USE mo_io_utilities, ONLY: vartype_int 
  USE mo_io_utilities, ONLY: vartype_real
  USE mo_io_utilities, ONLY: vartype_char

  USE mo_physical_constants, ONLY: re !< av. radius of the earth [m]

  USE mo_grid_structures, ONLY: icosahedral_triangular_grid
 
  USE mo_grid_structures, ONLY: igrid_icon

  
  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer

  USE mo_var_meta_data, ONLY: dim_aot_tg, &
    &                         aot_tg_meta, &
    &                         def_aot_tg_meta

  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  


  USE mo_var_meta_data, ONLY:  dim_icon, &
    &                          def_dimension_info_icon

  USE mo_var_meta_data, ONLY: nc_grid_def_icon, &
    &                         set_nc_grid_def_icon


  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(icosahedral_triangular_grid), INTENT(IN) :: icon_grid !< structure which contains the definition of the ICON grid
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  INTEGER (KIND=i8), INTENT(IN) :: ntype !< number of types of aerosols
  INTEGER (KIND=i8), INTENT(IN) :: ntime !< number of times
  INTEGER (KIND=i4), INTENT(IN) :: iaot_type !< ID of aeorosol raw data

  REAL (KIND=wp), INTENT(IN)  :: aot_tg(:,:,:,:,:) !< aerosol optical thickness, aot_tg(ie,je,ke,ntype,ntime)


  ! local variables
  REAL (KIND=wp),ALLOCATABLE :: time(:) !< time variable
  INTEGER (KIND=i8) :: dataDate  !< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
  INTEGER (KIND=i8) :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm



  INTEGER, PARAMETER :: nglob_atts=5
  TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

  TYPE(netcdf_grid_mapping) :: nc_grid_def !< mapping parameters for netcdf
  INTEGER :: ndims 
  INTEGER :: ncid
  INTEGER :: varid

  TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
  TYPE(dim_meta_info), TARGET :: dim_1d_icon(1:1)
  TYPE(dim_meta_info), TARGET :: dim_3d_buffer(1:3)


  INTEGER :: n_1d_real = 0 !< number of 1D real variables
  INTEGER :: n_3d_real = 0 !< number of 3D real variables

  INTEGER :: errorcode !< error status variable

  CHARACTER (len=80):: grid_mapping !< netcdf attribute grid mapping
  CHARACTER (len=80):: coordinates  !< netcdf attribute coordinates

  CHARACTER (len=80) :: attname
  CHARACTER (len=255) :: attributetext

  INTEGER :: n !< counter

  !-------------------------------------------------------------
  ! define global attributes
  IF (iaot_type == 1 ) THEN
     CALL set_global_att_aot(global_attributes)
  ELSEIF(iaot_type == 2 ) THEN
     CALL set_global_att_aot_aero(global_attributes)
  ELSEIF(iaot_type == 3 ) THEN
     CALL set_global_att_aot_macc(global_attributes)
  ELSE
     PRINT *, 'UNKNOWN AOT DATA OPTION: '
     STOP 1
  ENDIF

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

  !set up dimensions for ICON grid
  CALL def_dimension_info_icon(icon_grid)
  ! dim_icon
   dim_1d_icon = dim_icon(1) ! cell

  
  !set up dimensions for aot_tg
  CALL def_aot_tg_meta(tg,ntime,ntype,dim_1d_icon)
  ! dim_aot_tg and aot_tg_meta

  CALL def_com_target_fields_meta(dim_1d_icon)
  ! lon_geo_meta and lat_geo_meta

  ! set mapping parameters for netcdf
  grid_mapping="lon_lat_on_sphere"
  coordinates="lon lat"

  CALL set_nc_grid_def_icon(grid_mapping)
  ! nc_grid_def_icon

  
  ALLOCATE(time(1:ntime),STAT=errorcode)
    IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array time')
    DO n=1,ntime
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO


  !set up dimensions for buffer
  ndims = 4
  ALLOCATE(dim_list(1:ndims),STAT=errorcode)
  IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')

  dim_list(1) = dim_icon(1) ! cell
  dim_list(2) = dim_icon(2) ! nv,  number of vertices per cell
  dim_list(3)%dimname = 'time'
  dim_list(3)%dimsize = ntime
  dim_list(4)%dimname = 'ntype'
  dim_list(4)%dimsize = ntype
  
  dim_1d_icon(1) = dim_icon(1)

  dim_3d_buffer(1) = dim_list(1)
  dim_3d_buffer(2) = dim_list(3)
  dim_3d_buffer(3) = dim_list(4)

  ! set variable lists for output
  n_1d_real = 2 ! for ICON write out 3 1D real variables ( lon, lat)
  n_3d_real = 1  ! for ICON write out 1 3D real variable (aot_tg)

  !-----------------------------------------------------------------
    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
        &                       dim_list=dim_list,                  &
        &                       global_attributes=global_attributes, &
        &                       time=time,                           &
        &                       ncid=ncid)

  !-----------------------------------------------------------------
  ! start with real 1d variables
     

  ! lon
  CALL netcdf_put_var(ncid,lon_geo(1:icon_grid%ncell,1,1),lon_geo_meta,undefined)

  ! lat
  CALL netcdf_put_var(ncid,lat_geo(1:icon_grid%ncell,1,1),lat_geo_meta,undefined)
  !-----------------------------------------------------------------

  ! aot_tg
  CALL netcdf_put_var(ncid,aot_tg(1:icon_grid%ncell,1,1,1:ntype,1:ntime), &
       &                 aot_tg_meta, undefined)
     
  CALL close_netcdf_file(ncid)


  END SUBROUTINE write_netcdf_icon_grid_aot

  !----------------------------------------------------------------------------


  !> set global attributes for netcdf with aerosol optical thickness data
  !global climatology from Ina Tegen (Tegen et al. 1997)
  SUBROUTINE set_global_att_aot(global_attributes)
    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:5)

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
    global_attributes(1)%attributetext='Aerosol Optical Thickness'

    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='DWD'


    global_attributes(3)%attname = 'source'
    global_attributes(3)%attributetext='Global Aerosol Climatology Project'

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' aot_to_buffer'

    global_attributes(5)%attname = 'references'
    global_attributes(5)%attributetext='http://gacp.giss.nasa.gov/data_sets/transport/'

  END SUBROUTINE set_global_att_aot
!----------------------------------------------------------------------------

!> set global attributes for netcdf with aerosol optical thickness data    AeroCom1
!gs_21.03.12
    SUBROUTINE set_global_att_aot_aero(global_attributes)
    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:5)

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
    global_attributes(1)%attributetext='14 model median (AeroCom1): LO,LS,UL,SP,CT,MI,NF,OT,OG,IM,GM,GO,GI,GR'

    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='MPI_MET'


    global_attributes(3)%attname = 'source'
    global_attributes(3)%attributetext='Global Aerosol Climatology Project'

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute)

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' aot_to_buffer'

    global_attributes(5)%attname = 'references'
    global_attributes(5)%attributetext='created by stefan_kinne in 2010_11'

  END SUBROUTINE set_global_att_aot_aero

!>
!----------------------------------------------------------------------------

!> set global attributes for netcdf with aerosol optical thickness data 
!  from ECMWF-MACC II dataset

    SUBROUTINE set_global_att_aot_macc(global_attributes)
    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:5)

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
    global_attributes(1)%attributetext='monthly mean climatology of AOD compiled from MACC dataset 2003-2012'

    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='ECMWF'


    global_attributes(3)%attname = 'source'
    global_attributes(3)%attributetext='Monitoring atmospheric composition and climate(MACC) project'

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute)

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' aot_to_buffer'

    global_attributes(5)%attname = 'references'
    global_attributes(5)%attributetext='http://www.gmes-atmosphere.eu/'

  END SUBROUTINE set_global_att_aot_macc

!>
!----------------------------------------------------------------------------
!> read netcdf file for the AOT data in the buffer
  SUBROUTINE read_netcdf_buffer_aot(netcdf_filename,     &
   &                                     tg,             &
   &                                     ntype,          &
   &                                     ntime,          &
   &                                     aot_tg)

  USE mo_grid_structures, ONLY: target_grid_def

  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer

  USE mo_var_meta_data, ONLY: dim_aot_tg, &
    &                         aot_tg_meta, &
    &                         def_aot_tg_meta

  USE mo_io_utilities, ONLY: netcdf_get_var

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  INTEGER (KIND=i8), INTENT(IN) :: ntype !< number of types of aerosols
  INTEGER (KIND=i8), INTENT(IN) :: ntime !< number of times

  REAL (KIND=wp), INTENT(OUT)  :: aot_tg(:,:,:,:,:) !< aerosol optical thickness, aot_tg(ie,je,ke,ntype,ntime)

  ! local variables

  INTEGER :: n !< counter

  !-------------------------------------------------------------
  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg
  
  ! define dimensions and meta information for variable aot_tg for netcdf output
  CALL def_aot_tg_meta(tg,ntime,ntype,dim_3d_tg)
  ! dim_aot_tg and aot_tg_meta
  

  PRINT *,' read netcdf data aot'

  CALL netcdf_get_var(TRIM(netcdf_filename),aot_tg_meta,aot_tg)
  PRINT *,'aot_tg read'

  


  END SUBROUTINE read_netcdf_buffer_aot
  !----------------------------------------------------------------------------



  END MODULE mo_aot_output_nc

