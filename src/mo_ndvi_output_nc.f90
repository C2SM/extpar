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

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: rotated_lonlat_grid, &
       &                              icosahedral_triangular_grid, &
       &                              target_grid_def

  USE mo_io_utilities,          ONLY: netcdf_attributes, &
       &                              dim_meta_info, &
       &                              netcdf_put_var, &
       &                              open_new_netcdf_file, &
       &                              close_netcdf_file, &
       &                              netcdf_get_var, &
       &                              netcdf_def_grid_mapping, &
       &                              set_date_mm_extpar_field

  USE mo_var_meta_data,        ONLY:  dim_3d_tg, &
       &                              def_dimension_info_buffer, &
       &                              lon_geo_meta, &
       &                              lat_geo_meta, &
       &                              def_com_target_fields_meta, &  
       &                              ndvi_max_meta, &
       &                              ndvi_field_mom_meta, &
       &                              ndvi_ratio_mom_meta,&
       &                              def_ndvi_meta, &
       &                              nc_grid_def_cosmo, &
       &                              set_nc_grid_def_cosmo, &
       &                              dim_rlon_cosmo, &
       &                              dim_rlat_cosmo, &
       &                              dim_2d_cosmo,   &
       &                              rlon_meta,      &
       &                              rlat_meta,      &
       &                              def_dimension_info_cosmo, &
       &                              dim_icon, &
       &                              def_dimension_info_icon, &
       &                              set_nc_grid_def_icon
                               
  USE mo_cosmo_grid,            ONLY: lon_rot, lat_rot

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_netcdf_buffer_ndvi, &
       &    write_netcdf_cosmo_grid_ndvi, &
       &    write_netcdf_icon_grid_ndvi, &
       &    read_netcdf_buffer_ndvi

  CONTAINS

  SUBROUTINE write_netcdf_buffer_ndvi(netcdf_filename,  &
       &                              tg,         &
       &                              ntime, &
       &                              undefined, &
       &                              lon_geo,     &
       &                              lat_geo, &
       &                              ndvi_max,  &
       &                              ndvi_field_mom,&
       &                              ndvi_ratio_mom)


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file

    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description

    INTEGER (KIND=i4), INTENT(IN)      :: ntime !< number of times of ndvi data (12 monthly mean values)

    REAL(KIND=wp), INTENT(IN)          :: undefined, &       !< value to indicate undefined grid elements 
         &                                lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
         &                                lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
         &                                ndvi_max(:,:,:), & !< field for ndvi maximum
         &                                ndvi_field_mom(:,:,:,:), & !< field for monthly mean ndvi data (12 months)
         &                                ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)

    ! local variables
    REAL (KIND=wp),ALLOCATABLE        :: time(:) !< time variable
    INTEGER (KIND=i4)                 :: dataDate, &
         &                               dataTime, &  !< time, for edition independent use GRIB_API dataTime in the format hhmm
         &                               ndims, &  
         &                               ncid, &
         &                               errorcode, & !< error status variable
         &                               n !< counter

    INTEGER(KINd=i4), PARAMETER       :: nglob_atts=6

    TYPE(netcdf_attributes)           :: global_attributes(nglob_atts)
    TYPE(dim_meta_info), ALLOCATABLE  :: dim_list(:) !< dimensions for netcdf file

    CALL logging%info('Enter routine: write_netcdf_buffer_ndvi')

    ! define global attributes
    CALL set_global_att_ndvi(global_attributes)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various NDVI data related variables for netcdf output
    CALL def_ndvi_meta(ntime,dim_3d_tg)
    ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta
    
    ALLOCATE(time(1:ntime),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array time',__FILE__,__LINE__)
    DO n=1,ntime
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO

    ! set up dimensions for netcdf output 
    ndims = 4
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)
      
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

    CALL logging%info('Exit routine: write_netcdf_buffer_ndvi')

  END SUBROUTINE write_netcdf_buffer_ndvi

  SUBROUTINE write_netcdf_cosmo_grid_ndvi(netcdf_filename,  &
       &                                  cosmo_grid,         &
       &                                  tg,         &
       &                                  ntime, &
       &                                  undefined, &
       &                                  ndvi_max,  &
       &                                  ndvi_field_mom,&
       &                                  ndvi_ratio_mom)


    CHARACTER (len=*), INTENT(IN)          :: netcdf_filename !< filename for the netcdf file

    TYPE(rotated_lonlat_grid), INTENT(IN)  :: COSMO_grid      !< structure which contains the definition of the COSMO grid
    TYPE(target_grid_def), INTENT(IN)      :: tg !< structure with target grid description

    INTEGER (KIND=i4), INTENT(IN)          :: ntime !< number of times of ndvi data (12 monthly mean values)

    REAL(KIND=wp), INTENT(IN)              :: undefined, &       !< value to indicate undefined grid elements 
         &                                    ndvi_max(:,:,:), & !< field for ndvi maximum
         &                                    ndvi_field_mom(:,:,:,:), & !< field for monthly mean ndvi data (12 months)
         &                                    ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)

    ! local variables
    REAL (KIND=wp),ALLOCATABLE             :: time(:) !< time variable

    INTEGER (KIND=i4)                      :: dataDate, &
         &                                    dataTime, &
         &                                    ndims, &  
         &                                    ncid, &
         &                                    varid, n, &
         &                                    errorcode !< error status variable

    INTEGER(KIND=i4), PARAMETER           :: nglob_atts=6

    TYPE(dim_meta_info), ALLOCATABLE      :: dim_list(:) !< dimensions for netcdf file
    
    TYPE(netcdf_attributes)               :: global_attributes(nglob_atts)

    CHARACTER (len=80)                    :: grid_mapping, & !< netcdf attribute grid mapping
         &                                   coordinates  !< netcdf attribute coordinates

    CALL logging%info('Enter routine: write_netcdf_buffer_ndvi')

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
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_2d_cosmo,coordinates,grid_mapping)
    ! lon_geo_meta and lat_geo_meta

    !define meta information for various NDVI data related variables for netcdf output
    CALL def_ndvi_meta(ntime,dim_2d_cosmo,coordinates,grid_mapping)
    ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta

    ALLOCATE(time(1:ntime),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array time',__FILE__,__LINE__)
    DO n=1,ntime
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO

    ! set up dimensions for netcdf output 
    ndims = 3
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)
    dim_list(1) = dim_rlon_cosmo(1) ! rlon
    dim_list(2) = dim_rlat_cosmo(1) ! rlat
    dim_list(3)%dimname = 'time'
    dim_list(3)%dimsize = ntime
    
    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       time=time,          &
      &                       ncid=ncid)
    !-----------------------------------------------------------------

    ! rlon
    CALL netcdf_put_var(ncid,lon_rot(1:cosmo_grid%nlon_rot),rlon_meta,undefined)

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

    CALL logging%info('Exit routine: write_netcdf_buffer_ndvi')

   END SUBROUTINE write_netcdf_cosmo_grid_ndvi

   SUBROUTINE write_netcdf_icon_grid_ndvi(netcdf_filename,  &
        &                                 icon_grid,         &
        &                                 tg,         &
        &                                 ntime, &
        &                                 undefined, &
        &                                 lon_geo,     &
        &                                 lat_geo, &
        &                                 ndvi_max,  &
        &                                 ndvi_field_mom,&
        &                                 ndvi_ratio_mom)

    CHARACTER (len=*), INTENT(IN)                  :: netcdf_filename !< filename for the netcdf file

    TYPE(icosahedral_triangular_grid), INTENT(IN)  :: icon_grid      !< structure which contains the definition of the ICON grid
    TYPE(target_grid_def), INTENT(IN)              :: tg !< structure with target grid description

    INTEGER (KIND=i4), INTENT(IN)                  :: ntime !< number of times of ndvi data (12 monthly mean values)

    REAL(KIND=wp), INTENT(IN)                      :: undefined, &       !< value to indicate undefined grid elements 
         &                                            lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
         &                                            lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
         &                                            ndvi_max(:,:,:), & !< field for ndvi maximum
         &                                            ndvi_field_mom(:,:,:,:), & !< field for monthly mean ndvi data (12 months)
         &                                            ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)

    ! local variables
    REAL (KIND=wp),ALLOCATABLE                    :: time(:) !< time variable
    INTEGER (KIND=i4)                             :: dataDate, &
         &                                           dataTime, &
         &                                           ndims, errorcode, n, &
         &                                           ncid

    INTEGER(KINd=i4), PARAMETER                  :: nglob_atts=6

    TYPE(dim_meta_info), ALLOCATABLE             :: dim_list(:) !< dimensions for netcdf file
    TYPE(dim_meta_info), TARGET                  :: dim_1d_icon(1:1)
    TYPE(netcdf_attributes)                      :: global_attributes(nglob_atts)

    CHARACTER (len=80)                           :: grid_mapping !< netcdf attribute grid mapping

    CALL logging%info('Enter routine: write_netcdf_icon_grid_ndvi')

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

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_1d_icon)
    ! lon_geo_meta and lat_geo_meta

    !define meta information for various NDVI data related variables for netcdf output
    CALL def_ndvi_meta(ntime,dim_1d_icon)
    ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta

    ALLOCATE(time(1:ntime),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array time',__FILE__,__LINE__)
    DO n=1,ntime
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO

    ! set up dimensions for netcdf output 
    ndims = 2
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)
    dim_list(1) =  dim_icon(1) ! cell
    dim_list(2)%dimname = 'time'
    dim_list(2)%dimsize = ntime

    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
        &                       dim_list=dim_list,                  &
        &                       global_attributes=global_attributes, &
        &                       time=time,          &
        &                       ncid=ncid)

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

    CALL logging%info('Exit routine: write_netcdf_icon_grid_ndvi')

  END SUBROUTINE write_netcdf_icon_grid_ndvi

  !> set global attributes for netcdf with NDVI data
  SUBROUTINE set_global_att_ndvi(global_attributes)

    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:6)

    !local variables
    CHARACTER(len=10)                      :: ydate, &
         &                                    ytime

    CHARACTER(len=2)                       :: cc, &
         &                                    yy, &
         &                                    mm, &
         &                                    dd, &
         &                                    hh, &
         &                                    minute

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
       &                             tg,         &
       &                             ntime, &
       &                             ndvi_max,  &
       &                             ndvi_field_mom,&
       &                             ndvi_ratio_mom)


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(INOUT)   :: ntime !< number of times of ndvi data (12 monthly mean values)

    REAL (KIND=wp), INTENT(OUT)        :: ndvi_max(:,:,:), & !< field for ndvi maximum
         &                                ndvi_field_mom(:,:,:,:), & !< field for monthly mean ndvi data (12 months)
         &                                ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)

    ! local variables
    INTEGER(KIND=i4), PARAMETER       :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_ndvi')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various NDVI data related variables for netcdf output
    CALL def_ndvi_meta(ntime,dim_3d_tg)
    ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),ndvi_max_meta,ndvi_max)

    CALL netcdf_get_var(TRIM(netcdf_filename),ndvi_field_mom_meta,ndvi_field_mom)

    CALL netcdf_get_var(TRIM(netcdf_filename),ndvi_ratio_mom_meta,ndvi_ratio_mom)
    
    CALL logging%info('Exit routine: read_netcdf_buffer_ndvi')

   END SUBROUTINE read_netcdf_buffer_ndvi

END MODULE mo_ndvi_output_nc
