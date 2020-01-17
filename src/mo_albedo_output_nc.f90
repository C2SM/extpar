!+ Fortran module for netcdf output of Albedo data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_8         2013/03/12 Frank Brenner
!  introduced MODIS albedo dataset(s) as new external parameter(s) 
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for netcdf output of Albedo data
!> \author Hermann Asensio, Frank Brenner
MODULE mo_albedo_output_nc

  
  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  !> data type structures form module GRID_structures
  USE mo_grid_structures, ONLY: target_grid_def

  USE mo_io_utilities, ONLY: netcdf_attributes

  USE mo_io_utilities, ONLY: dim_meta_info

  USE mo_io_utilities, ONLY: netcdf_put_var
  USE mo_io_utilities, ONLY: open_new_netcdf_file
  USE mo_io_utilities, ONLY: close_netcdf_file

  USE mo_io_utilities, ONLY: set_date_mm_extpar_field

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_netcdf_buffer_alb
  PUBLIC :: read_netcdf_buffer_alb

  CONTAINS

  SUBROUTINE write_netcdf_buffer_alb(netcdf_filename,  &
   &                                     tg,         &
   &                                     ntime, &
   &                                     undefined, &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     alb_field_mom, &
   &                                     alnid_field_mom, &
   &                                     aluvd_field_mom, &
   &                                     alb_dry, &
   &                                     alb_sat)


    USE mo_var_meta_data, ONLY: dim_3d_tg, &
      &                         def_dimension_info_buffer

    USE mo_var_meta_data, ONLY: lon_geo_meta, &
      &                         lat_geo_meta, &
      &                         def_com_target_fields_meta  

    USE mo_var_meta_data, ONLY: alb_field_mom_meta, &
      &                         alnid_field_mom_meta, &
      &                         aluvd_field_mom_meta, &
      &                         alb_dry_meta, &
      &                         alb_sat_meta, &
      &                         def_alb_meta


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN) :: ntime !< number of times of input data (12 monthly mean values)
    REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
    REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
    REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
    REAL (KIND=wp), INTENT(IN), OPTIONAL :: alb_field_mom(:,:,:,:) !< field for monthly mean albedo data (12 months)
    REAL (KIND=wp), INTENT(IN), OPTIONAL :: alnid_field_mom(:,:,:,:)
    REAL (KIND=wp), INTENT(IN), OPTIONAL :: aluvd_field_mom(:,:,:,:)
    REAL (KIND=wp), INTENT(IN), OPTIONAL :: alb_dry(:,:,:)
    REAL (KIND=wp), INTENT(IN), OPTIONAL :: alb_sat(:,:,:)

    ! local variables
    ! local variables
    REAL (KIND=wp),ALLOCATABLE :: time(:) !< time variable
    INTEGER (KIND=i8) :: dataDate  
!< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
   INTEGER (KIND=i8) :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm

    INTEGER :: ndims  
    INTEGER :: ncid

    TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
    
    INTEGER, PARAMETER :: nglob_atts=6
    TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

    INTEGER :: errorcode !< error status variable

    INTEGER :: n !< counter

    PRINT *,'ENTER write_netcdf_buffer_alb'

    PRINT *,'set_global_att_alb'

    !-------------------------------------------------------------
    ! define global attributes
    CALL set_global_att_alb(global_attributes)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    PRINT *,'def_com_target_fields_meta'
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define albedo meta information, related variables for netcdf output
    CALL def_alb_meta(ntime,dim_3d_tg)

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

    PRINT *, ' ============= In Routine: write_netcdf_buffer_alb =============='
    PRINT *, TRIM(netcdf_filename)

    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       time=time,          &
      &                       ncid=ncid)

    ! lon
    CALL netcdf_put_var(ncid,lon_geo,lon_geo_meta,undefined)

    ! lat
    CALL netcdf_put_var(ncid,lat_geo,lat_geo_meta,undefined)

 
    ! alb_field_mom
    IF (PRESENT(alb_field_mom)) THEN
      CALL netcdf_put_var(ncid,alb_field_mom,alb_field_mom_meta,undefined)
    ENDIF
    IF (PRESENT(alnid_field_mom)) THEN
      CALL netcdf_put_var(ncid,alnid_field_mom,alnid_field_mom_meta,undefined)
    ENDIF
    IF (PRESENT(aluvd_field_mom)) THEN
      CALL netcdf_put_var(ncid,aluvd_field_mom,aluvd_field_mom_meta,undefined)
    ENDIF

    IF (PRESENT(alb_dry)) THEN
      CALL netcdf_put_var(ncid,alb_dry,alb_dry_meta,undefined)
    ENDIF
    IF (PRESENT(alb_sat)) THEN
      CALL netcdf_put_var(ncid,alb_sat,alb_sat_meta,undefined)
    ENDIF


    CALL close_netcdf_file(ncid)

   END SUBROUTINE write_netcdf_buffer_alb
   !-----------------------------------------------------------------
   !-----------------------------------------------------------------
   !-----------------------------------------------------------------

  !> set global attributes for netcdf with albedo data
  SUBROUTINE set_global_att_alb(global_attributes)

    USE mo_albedo_data, ONLY : ialb_type

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

    IF (ialb_type == 2) THEN
      global_attributes(1)%attname = 'title'
      global_attributes(1)%attributetext='soil albedo data '
      global_attributes(2)%attname = 'institution'
      global_attributes(2)%attributetext='ETH Zurich'
      global_attributes(3)%attname = 'source'
      global_attributes(3)%attributetext='CESM-CLM 4.0'
      global_attributes(5)%attname = 'references'
      global_attributes(5)%attributetext='https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/'// & !_br 21.02.14 spitted too long line
           'inputdata/lnd/clm2/rawdata/mksrf_soilcol.081008.nc & JGR:Lawrence et al. (2007)'  !_br 21.02.14
    ELSE    
      global_attributes(1)%attname = 'title'
      global_attributes(1)%attributetext='albedo data '
      global_attributes(2)%attname = 'institution'
      global_attributes(2)%attributetext='Deutscher Wetterdienst'
      global_attributes(3)%attname = 'source'
      global_attributes(3)%attributetext='NASA MODIS'
      global_attributes(5)%attname = 'references'
      global_attributes(5)%attributetext='http://www-modis.bu.edu/brdf/'
    ENDIF


    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' alb_to_buffer'

    global_attributes(6)%attname = 'comment'
    global_attributes(6)%attributetext=''

  END SUBROUTINE set_global_att_alb

  SUBROUTINE read_netcdf_buffer_alb(netcdf_filename,  &
   &                                     tg,         &
   &                                     ntime, &
   &                                     alb_field_mom, &
   &                                     alnid_field_mom, &
   &                                     aluvd_field_mom, &
   &                                     alb_dry, &
   &                                     alb_sat)

    USE mo_var_meta_data, ONLY: dim_3d_tg, &
      &                         def_dimension_info_buffer

    USE mo_var_meta_data, ONLY: def_com_target_fields_meta  

    USE mo_var_meta_data, ONLY: alb_field_mom_meta, &
      &                         alnid_field_mom_meta, &
      &                         aluvd_field_mom_meta, &
      &                         alb_dry_meta, &
      &                         alb_sat_meta, &
      &                         def_alb_meta

    USE mo_io_utilities, ONLY: netcdf_get_var


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(OUT) :: ntime !< number of times of input data (12 monthly mean values)
    REAL (KIND=wp), INTENT(OUT), OPTIONAL :: alb_field_mom(:,:,:,:) !< field for monthly mean albedo data (12 months)
    REAL (KIND=wp), INTENT(OUT), OPTIONAL :: alnid_field_mom(:,:,:,:)
    REAL (KIND=wp), INTENT(OUT), OPTIONAL :: aluvd_field_mom(:,:,:,:)
    REAL (KIND=wp), INTENT(OUT), OPTIONAL :: alb_dry(:,:,:)
    REAL (KIND=wp), INTENT(OUT), OPTIONAL :: alb_sat(:,:,:)

    ! local variables
    
    INTEGER, PARAMETER :: nglob_atts=6

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! dim_3d_tg
    PRINT *,'def_com_target_fields_meta'
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta

    !define albedo meta information, related variables for netcdf output
    CALL def_alb_meta(ntime,dim_3d_tg)

    PRINT *,'CALL read netcdf data ALB'
    PRINT *, TRIM(netcdf_filename)

    IF (PRESENT(alb_field_mom)) THEN
      CALL netcdf_get_var(TRIM(netcdf_filename),alb_field_mom_meta,alb_field_mom)
      PRINT *,'alb_field_mom read'
    ENDIF
    IF (PRESENT(alnid_field_mom)) THEN
      CALL netcdf_get_var(TRIM(netcdf_filename),alnid_field_mom_meta,alnid_field_mom)
      PRINT *,'alnid_field_mom read'
    ENDIF
    IF (PRESENT(aluvd_field_mom)) THEN
      CALL netcdf_get_var(TRIM(netcdf_filename),aluvd_field_mom_meta,aluvd_field_mom)
      PRINT *,'aluvd_field_mom read'
    ENDIF

    IF (PRESENT(alb_dry)) THEN
      CALL netcdf_get_var(TRIM(netcdf_filename),alb_dry_meta,alb_dry)
      PRINT *,'alb_dry read'
    ENDIF
    IF (PRESENT(alb_sat)) THEN
      CALL netcdf_get_var(TRIM(netcdf_filename),alb_sat_meta,alb_sat)
      PRINT *,'alb_sat read'
    ENDIF


   END SUBROUTINE read_netcdf_buffer_alb

END Module mo_albedo_output_nc
