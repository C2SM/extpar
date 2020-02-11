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
! V4_0         2016/08/17 authors from RHM and Daniel Lthi
!  added support for MACv2 climatological aerosol fields
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for aerosol input output routines
!> \author Hermann Asensio
MODULE mo_aot_output_nc

  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4

  USE mo_grid_structures,       ONLY: target_grid_def, &
    &                                 rotated_lonlat_grid, &
    &                                 icosahedral_triangular_grid

  USE mo_io_utilities,          ONLY: netcdf_attributes, &
    &                                 netcdf_put_var, &
    &                                 open_new_netcdf_file, &
    &                                 close_netcdf_file, &
    &                                 set_date_mm_extpar_field, &
    &                                 netcdf_def_grid_mapping, &
    &                                 netcdf_get_var, &
    &                                 dim_meta_info

  USE mo_var_meta_data,         ONLY: dim_2d_tg, &
    &                                 dim_3d_tg, &
    &                                 def_dimension_info_buffer, &
    &                                 aot_tg_MAC_meta,&
    &                                 ssa_tg_MAC_meta,&
    &                                 asy_tg_MAC_meta,&
    &                                 def_aot_tg_meta, &
    &                                 lat_geo_meta, &
    &                                 def_com_target_fields_meta, &
    &                                 lon_geo_meta, &
    &                                 dim_rlon_cosmo, &
    &                                 dim_rlat_cosmo, &
    &                                 dim_2d_cosmo,   &
    &                                 rlon_meta,      &
    &                                 rlat_meta,      &
    &                                 def_dimension_info_cosmo, &
    &                                 nc_grid_def_cosmo, &
    &                                 set_nc_grid_def_cosmo, &
    &                                 dim_icon, &
    &                                 def_dimension_info_icon, &
    &                                 set_nc_grid_def_icon, &
    &                                 aot_tg_meta

  USE mo_cosmo_grid,            ONLY: lon_rot, lat_rot

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_netcdf_buffer_aot, &
    &       write_netcdf_cosmo_grid_aot, &
    &       write_netcdf_icon_grid_aot, &
    &       read_netcdf_buffer_aot, &
    &       read_netcdf_buffer_aot_MAC


  CONTAINS

  !> create a netcdf file for the AOT data in the buffer
  SUBROUTINE write_netcdf_buffer_aot(netcdf_filename,  &
   &                                     tg,         &
   &                                     undefined, &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     ntype,           &
   &                                     ntime,        &
   &                                     n_spectr,       & !new
   &                                     aot_tg,         &
   &                                     MAC_aot_tg,     &
   &                                     MAC_ssa_tg,     &
   &                                     MAC_asy_tg,     &
   &                                     iaot_type)


    CHARACTER (len=*),     INTENT(IN):: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN):: tg !< structure with target grid description
    REAL(KIND=wp),         INTENT(IN):: undefined, &
      &                                 lon_geo(:,:,:), &
      &                                 lat_geo(:,:,:), &
      &                                 aot_tg(:,:,:,:,:), & 
      &                                 MAC_aot_tg(:,:,:,:), &
      &                                 MAC_ssa_tg(:,:,:,:), &
      &                                 MAC_asy_tg(:,:,:,:)

    INTEGER (KIND=i4), INTENT(IN)    :: ntype, & !< number of types of aerosols
      &                                 ntime, & !< number of times
      &                                 n_spectr, & !< number of spectral bands
      &                                 iaot_type !< ID of aeorosol raw data

    
    ! local variables
    REAL (KIND=wp),ALLOCATABLE :: time(:) !< time variable
    INTEGER (KIND=i4) :: dataDate, &  !< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
      &                  dataTime, & !< time, for edition independent use GRIB_API dataTime in the format hhmm
      &                  ndims,ncid, errorcode, &
      &                  n

    INTEGER, PARAMETER :: nglob_atts=5
    TYPE(netcdf_attributes) :: global_attributes(nglob_atts)
    TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file

    CALL logging%info('Enter routine: write_netcdf_buffer_aot')
    !-------------------------------------------------------------
    ! define global attributes
    IF (iaot_type == 1 ) THEN
      CALL set_global_att_aot(global_attributes)
    ELSEIF(iaot_type == 2 ) THEN
       CALL set_global_att_aot_aero(global_attributes)
    ELSEIF(iaot_type == 3 ) THEN
       CALL set_global_att_aot_macc(global_attributes)
    ELSEIF(iaot_type == 4 ) THEN
       CALL set_global_att_aot_MACv2(global_attributes)
    ELSE
       CALL logging%error('Unknown aot data option',__FILE__,__LINE__)
    ENDIF
    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    
    ! define dimensions and meta information for variable aot_tg for netcdf output
    IF (iaot_type == 4) THEN
      CALL def_aot_tg_meta(ntime,ntype,dim_3d_tg,n_spectr=n_spectr) !new

      ! dim_aot_tg and aot_tg_metar
      ! define meta information for target field variables lon_geo, lat_geo 
      CALL def_com_target_fields_meta(dim_3d_tg)
      ! lon_geo_meta and lat_geo_meta

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

      dim_list(1) = dim_2d_tg(1) ! ie
      dim_list(2) = dim_2d_tg(2) ! je
      dim_list(3)%dimname = 'spectr'
      dim_list(3)%dimsize = n_spectr
      dim_list(4)%dimname = 'time'
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

      CALL netcdf_put_var(ncid,MAC_aot_tg,aot_tg_MAC_meta,undefined)
      CALL netcdf_put_var(ncid,MAC_ssa_tg,ssa_tg_MAC_meta,undefined)
      CALL netcdf_put_var(ncid,MAC_asy_tg,asy_tg_MAC_meta,undefined)

    ELSE
      CALL def_aot_tg_meta(ntime,ntype,dim_3d_tg)
      ! dim_aot_tg and aot_tg_meta
      
      ! define meta information for target field variables lon_geo, lat_geo 
      CALL def_com_target_fields_meta(dim_3d_tg)
      ! lon_geo_meta and lat_geo_meta

      ALLOCATE(time(1:ntime),STAT=errorcode)
      IF (errorcode /= 0 ) CALL logging%error('Cant allocate array time',__FILE__,__LINE__)
      DO n=1,ntime
        CALL set_date_mm_extpar_field(n,dataDate,dataTime)
        time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
      ENDDO

      ! set up dimensions for netcdf output 
      ndims = 5
      ALLOCATE(dim_list(1:ndims),STAT=errorcode)
      IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)

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
    ENDIF

    CALL close_netcdf_file(ncid)

    CALL logging%info('Exit routine: write_netcdf_buffer_aot')

  END SUBROUTINE write_netcdf_buffer_aot
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------


  !> create a netcdf file for the AOT data in the COSMO grid
  SUBROUTINE write_netcdf_cosmo_grid_aot(netcdf_filename,  &
   &                                     cosmo_grid,       &
   &                                     tg,               &
   &                                     undefined, &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     ntype,           &
   &                                     ntime,        &
   &                                     n_spectr,         & !new
   &                                     aot_tg,           &
   &                                     MAC_aot_tg,       &
   &                                     MAC_ssa_tg,       &
   &                                     MAC_asy_tg,       &
   &                                     iaot_type)




    CHARACTER (len=*), INTENT(IN)         :: netcdf_filename !< filename for the netcdf file
    TYPE(rotated_lonlat_grid), INTENT(IN) :: cosmo_grid !< structure which contains the definition of the COSMO grid
    TYPE(target_grid_def), INTENT(IN)     :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN)         :: ntype, & !< number of types of aerosols
      &                                      ntime, & !< number of times
      &                                      n_spectr, & !< number of spectral bands
      &                                      iaot_type !< ID of aeorosol raw data

    REAL(KIND=wp), INTENT(IN)             :: undefined, &       !< value to indicate undefined grid elements 
      &                                      aot_tg(:,:,:,:,:), & !< aerosol optical thickness, aot_tg(ie,je,ke,ntype,ntime)
      &                                      lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
      &                                      lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
      &                                      MAC_aot_tg(:,:,:,:), &
      &                                      MAC_ssa_tg(:,:,:,:), &
      &                                      MAC_asy_tg(:,:,:,:)

    ! local variables
    REAL (KIND=wp),ALLOCATABLE            :: time(:) !< time variable
    INTEGER, PARAMETER                    :: nglob_atts=5

    TYPE(netcdf_attributes)               :: global_attributes(nglob_atts)
    TYPE(dim_meta_info), ALLOCATABLE      :: dim_list(:) !< dimensions for netcdf file

    INTEGER(KIND=i4)                      :: ndims, ncid, varid, n, &
      &                                      errorcode,&
      &                                      dataDate, &
      &                                      dataTime


    CHARACTER (len=80)                    :: grid_mapping !< netcdf attribute grid mapping
    CHARACTER (len=80)                    :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=1), PARAMETER          :: c_undef = "-" !< default charcter for undefined string


    CALL logging%info('Enter routine: write_netcdf_cosmo_grid_aot')
    !-------------------------------------------------------------
    ! define global attributes
    IF (iaot_type == 1) THEN
      CALL set_global_att_aot(global_attributes)
    ELSEIF (iaot_type == 2) THEN
      CALL set_global_att_aot_aero(global_attributes)
    ELSEIF (iaot_type == 3) THEN
      CALL set_global_att_aot_macc(global_attributes)
    ELSEIF (iaot_type == 4) THEN
      CALL set_global_att_aot_MACv2(global_attributes)
    ELSE
      CALL logging%error('Unknown aot data option',__FILE__,__LINE__)
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
    CALL def_aot_tg_meta(ntime,ntype,dim_2d_cosmo,coordinates,grid_mapping)
    ! dim_aot_tg and aot_tg_meta

    CALL def_com_target_fields_meta(dim_2d_cosmo,coordinates,grid_mapping)
    ! lon_geo_meta and lat_geo_meta

    
    ALLOCATE(time(1:ntime),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array time',__FILE__,__LINE__)
    DO n=1,ntime
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO

    !set up dimensions for buffer
    ndims = 4
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)

    dim_list(1) = dim_rlon_cosmo(1) ! rlon
    dim_list(2) = dim_rlat_cosmo(1) ! rlat
    IF (iaot_type == 4) THEN
      dim_list(3)%dimname = 'spectr' 
      dim_list(3)%dimsize = n_spectr 
    ELSE
      dim_list(3)%dimname = 'ntype'
      dim_list(3)%dimsize = ntype
    ENDIF
    dim_list(4)%dimname = 'time' 
    dim_list(4)%dimsize = ntime 

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

    IF (iaot_type == 4) THEN
      CALL netcdf_put_var(ncid, MAC_aot_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:n_spectr,1:ntime), & 
                         & aot_tg_MAC_meta, &
                         & undefined)
      CALL netcdf_put_var(ncid, MAC_ssa_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:n_spectr,1:ntime), & 
                         & ssa_tg_MAC_meta, &
                         & undefined)
      CALL netcdf_put_var(ncid, MAC_asy_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:n_spectr,1:ntime), & 
                         & asy_tg_MAC_meta, &
                         & undefined)
    ELSE
    ! aot_tg
     CALL netcdf_put_var(ncid,&
                         &  aot_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntype,1:ntime), &
                         & aot_tg_meta, &
                         & undefined)
    ENDIF
    !-----------------------------------------------------------------
    CALL netcdf_def_grid_mapping(ncid, nc_grid_def_cosmo, varid)

    CALL close_netcdf_file(ncid)

    CALL logging%info('Exit routine: write_netcdf_cosmo_grid_aot')

  END SUBROUTINE write_netcdf_cosmo_grid_aot
  !----------------------------------------------------------------------------

  !> create a netcdf file for the AOT data in the ICON grid
  SUBROUTINE write_netcdf_icon_grid_aot(netcdf_filename,  &
   &                                     icon_grid,       &
   &                                     tg,              &
   &                                     undefined, &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     ntype,           &
   &                                     ntime,        &
   &                                     n_spectr,        & !new
   &                                     aot_tg,          &
   &                                     iaot_type)


    CHARACTER (len=*), INTENT(IN)                :: netcdf_filename !< filename for the netcdf file
    TYPE(icosahedral_triangular_grid), INTENT(IN):: icon_grid !< structure which contains the definition of the ICON grid
    TYPE(target_grid_def), INTENT(IN)            :: tg !< structure with target grid description

    REAL(KIND=wp), INTENT(IN)                    :: undefined, &       !< value to indicate undefined grid elements 
      &                                             lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
      &                                             lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
      &                                             aot_tg(:,:,:,:,:) !< aerosol optical thickness, aot_tg(ie,je,ke,ntype,ntime)

    INTEGER (KIND=i4), INTENT(IN)                :: ntype, & !< number of types of aerosols
      &                                             ntime, & !< number of times
      &                                             n_spectr, & !< number of times new  
      &                                             iaot_type !< ID of aeorosol raw data


    ! local variables
    REAL (KIND=wp),ALLOCATABLE                   :: time(:) !< time variable
    INTEGER (KIND=i4)                            :: dataDate, dataTime, ndims, ncid, n, &
      &                                             errorcode

    INTEGER, PARAMETER                           :: nglob_atts=5
    TYPE(netcdf_attributes)                      :: global_attributes(nglob_atts)
    TYPE(dim_meta_info), ALLOCATABLE             :: dim_list(:) !< dimensions for netcdf file
    TYPE(dim_meta_info), TARGET                  :: dim_1d_icon(1:1)

    CHARACTER (len=80)                           :: grid_mapping !< netcdf attribute grid mapping

    CALL logging%info('Enter routine: write_netcdf_icon_grid_aot')
  !-------------------------------------------------------------
    ! define global attributes
    IF (iaot_type == 1) THEN
    CALL set_global_att_aot(global_attributes)
    ELSEIF (iaot_type == 2) THEN
       CALL set_global_att_aot_aero(global_attributes)
    ELSEIF (iaot_type == 3) THEN
       CALL set_global_att_aot_macc(global_attributes)
    ELSEIF (iaot_type == 4) THEN
       CALL set_global_att_aot_MACv2(global_attributes)
    ELSE
      CALL logging%error('Unknown aot data option',__FILE__,__LINE__)
    ENDIF

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg

    !set up dimensions for ICON grid
    CALL def_dimension_info_icon(icon_grid)
    ! dim_icon
     dim_1d_icon = dim_icon(1) ! cell

    !set up dimensions for aot_tg
    CALL def_aot_tg_meta(ntime,ntype,dim_1d_icon,n_spectr=n_spectr)
    ! dim_aot_tg and aot_tg_meta

    CALL def_com_target_fields_meta(dim_1d_icon)
    ! lon_geo_meta and lat_geo_meta

    ! set mapping parameters for netcdf
    grid_mapping="lon_lat_on_sphere"

    CALL set_nc_grid_def_icon(grid_mapping)
    ! nc_grid_def_icon

    ALLOCATE(time(1:ntime),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array time',__FILE__,__LINE__)
    DO n=1,ntime
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO

    !set up dimensions for buffer
    ndims = 4
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)

    dim_list(1) = dim_icon(1) ! cell
    dim_list(2) = dim_icon(2) ! nv,  number of vertices per cell
    dim_list(3)%dimname = 'time'
    dim_list(3)%dimsize = ntime
    dim_list(4)%dimname = 'ntype'
    dim_list(4)%dimsize = ntype
    
    dim_1d_icon(1) = dim_icon(1)


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

    CALL logging%info('Exit routine: write_netcdf_icon_grid_aot')

  END SUBROUTINE write_netcdf_icon_grid_aot

   !----------------------------------------------------------------------------
  !> set global attributes for netcdf with aerosol optical thickness data
  !global climatology from MACv2(Kinne et al. 1997)
  SUBROUTINE set_global_att_aot_MACv2(global_attributes)
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
    global_attributes(1)%attributetext='Aerosol Optical Thickness from MACv2'

    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='MPI'


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
    global_attributes(5)%attributetext='created by stefan_kinne in 2015..'

  END SUBROUTINE set_global_att_aot_MACv2
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

  !> read netcdf file for the AOT data in the buffer
  SUBROUTINE read_netcdf_buffer_aot(netcdf_filename,  &
   &                                     tg,         &
   &                                     ntype,           &
   &                                     ntime,        &
   &                                     aot_tg)


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN)      :: ntype, & !< number of types of aerosols
         &                                ntime !< number of times

    REAL (KIND=wp), INTENT(OUT)        :: aot_tg(:,:,:,:,:) !< aerosol optical thickness, aot_tg(ie,je,ke,ntype,ntime)

    CALL logging%info('Enter routine: read_netcdf_buffer_aot')
    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    
    ! define dimensions and meta information for variable aot_tg for netcdf output
    CALL def_aot_tg_meta(ntime,ntype,dim_3d_tg)
    ! dim_aot_tg and aot_tg_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),aot_tg_meta,aot_tg)

    CALL logging%info('Exit routine: read_netcdf_buffer_aot')

  END SUBROUTINE read_netcdf_buffer_aot
  
  !> read netcdf file for the AOT data in the buffer
  SUBROUTINE read_netcdf_buffer_aot_MAC (netcdf_filename,     &
   &                                     tg,             &
   &                                     ntype,          &
   &                                     ntime,          &
   &                                     n_spectr,       & 
   &                                     MAC_aot_tg,     &
   &                                     MAC_ssa_tg,     &
   &                                     MAC_asy_tg)


  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description

  INTEGER (KIND=i4), INTENT(IN)      :: ntype, & !< number of types of aerosols
       &                                ntime, & !< number of times
       &                                n_spectr !< number of times new  
     
  REAL (KIND=wp), INTENT(OUT)        :: MAC_aot_tg(:,:,:,:), & !< aerosol optical thickness
       &                                MAC_ssa_tg(:,:,:,:), & !< single scattering albedo
       &                                MAC_asy_tg(:,:,:,:) !< factor asymmetry

    CALL logging%info('Enter routine: read_netcdf_buffer_aot_MAC')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg

    ! define dimensions and meta information for variable aot_tg for netcdf output
    CALL def_aot_tg_meta(ntime,ntype,dim_2d_tg,n_spectr=n_spectr)
    ! dim_aot_tg and aot_tg_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),aot_tg_MAC_meta,MAC_aot_tg)
    CALL netcdf_get_var(TRIM(netcdf_filename),ssa_tg_MAC_meta,MAC_ssa_tg)
    CALL netcdf_get_var(TRIM(netcdf_filename),asy_tg_MAC_meta,MAC_asy_tg)

    CALL logging%info('Exit routine: read_netcdf_buffer_aot_MAC')

  END SUBROUTINE read_netcdf_buffer_aot_MAC

END MODULE mo_aot_output_nc
