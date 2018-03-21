!+ <A one line description of this module>
!
!
! Description:
! <Say what this module is for>
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_10        2018-02-19 Juergen Helmert 
!  lsubtract_mean_slope, ERA-I surface temp for land points
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!
MODULE mo_era_output_nc
!
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
  USE mo_io_utilities, ONLY: netcdf_attributes

  USE mo_io_utilities, ONLY: dim_meta_info

  USE mo_io_utilities, ONLY: netcdf_put_var
  USE mo_io_utilities, ONLY: open_new_netcdf_file
  USE mo_io_utilities, ONLY: close_netcdf_file
  USE mo_io_utilities, ONLY: netcdf_def_grid_mapping

  USE mo_io_utilities, ONLY: get_date_const_field
  USE mo_io_utilities, ONLY: set_date_mm_extpar_field



  USE mo_io_utilities, ONLY: vartype_int 
  USE mo_io_utilities, ONLY: vartype_real
  USE mo_io_utilities, ONLY: vartype_char

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  IMPLICIT NONE 

  PRIVATE



  PUBLIC :: read_netcdf_buffer_sst
  PUBLIC :: read_netcdf_buffer_t2m

  CONTAINS


  !> set global attributes for netcdf with SST data
  SUBROUTINE set_global_att_sst(global_attributes)
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
    global_attributes(1)%attributetext='SST data '
    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='Deutscher Wetterdienst'

    global_attributes(3)%attname = 'source'
    global_attributes(3)%attributetext='ECMWF'

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' sst_to_buffer'

    global_attributes(5)%attname = 'references'
    global_attributes(5)%attributetext='ERA-I'

    global_attributes(6)%attname = 'comment'
    global_attributes(6)%attributetext=''

  END SUBROUTINE set_global_att_sst
  !-----------------------------------------------------------------------

  SUBROUTINE read_netcdf_buffer_sst(netcdf_filename,  &
   &                                     tg,         &
   &                                     ntime, &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     sst_field,&
   &                                     wsnow_field)

    USE mo_var_meta_data, ONLY: dim_3d_tg, dim_icon, &
      &                         def_dimension_info_buffer

    USE mo_var_meta_data, ONLY: lon_geo_meta, &
      &                         lat_geo_meta, &
      &                         no_raw_data_pixel_meta, &
      &                         def_com_target_fields_meta  

    USE mo_var_meta_data, ONLY: dim_era_tg
    USE mo_var_meta_data, ONLY: sst_field_meta, &
      &                         wsnow_field_meta,&
      &                         def_era_meta


    USE mo_io_utilities, ONLY: netcdf_get_var


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(INOUT) :: ntime !< number of times of sst data (12 monthly mean values)
    REAL(KIND=wp), INTENT(OUT)          :: undefined       !< value to indicate undefined grid elements 
    INTEGER, INTENT(OUT)                :: undef_int       !< value to indicate undefined grid elements
    REAL (KIND=wp), INTENT(OUT) :: sst_field(:,:,:,:) !< field for monthly mean sst data (12 months)
    REAL (KIND=wp), INTENT(OUT) :: wsnow_field(:,:,:,:) !< field for monthly sst ratio (12 months)


    ! local variables
    INTEGER :: ndims  
    TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
    
    INTEGER, PARAMETER :: nglob_atts=6
    TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

    INTEGER :: errorcode !< error status variable

    INTEGER :: n !< counter





    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    PRINT *,'def_com_target_fields_meta'
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various SST data related variables for netcdf output
    CALL def_era_meta(tg,ntime,dim_3d_tg)
    ! dim_sst_tg, sst_field_meta, wsnow_field_meta

    PRINT *,'CALL read netcdf data SST'


    CALL netcdf_get_var(TRIM(netcdf_filename),sst_field_meta,sst_field)
    PRINT *,'sst_field read'

    CALL netcdf_get_var(TRIM(netcdf_filename),wsnow_field_meta,wsnow_field)
    PRINT *,'wsnow_field read'



   END SUBROUTINE read_netcdf_buffer_sst
   !-----------------------------------------------------------------

  SUBROUTINE read_netcdf_buffer_t2m(netcdf_filename,  &
   &                                     tg,         &
   &                                     ntime, &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     t2m_field,&
   &                                     hsurf_field)

    USE mo_var_meta_data, ONLY: dim_3d_tg, dim_icon, &
      &                         def_dimension_info_buffer

    USE mo_var_meta_data, ONLY: lon_geo_meta, &
      &                         lat_geo_meta, &
      &                         no_raw_data_pixel_meta, &
      &                         def_com_target_fields_meta  

    USE mo_var_meta_data, ONLY: dim_era_tg
    USE mo_var_meta_data, ONLY: t2m_field_meta, &
      &                         hsurf_field_meta,&
      &                         def_era_meta


    USE mo_io_utilities, ONLY: netcdf_get_var


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(INOUT) :: ntime !< number of times of t2m data (12 monthly mean values)
    REAL(KIND=wp), INTENT(OUT)          :: undefined       !< value to indicate undefined grid elements 
    INTEGER, INTENT(OUT)                :: undef_int       !< value to indicate undefined grid elements
    REAL (KIND=wp), INTENT(OUT) :: t2m_field(:,:,:,:) !< field for monthly mean t2m data (12 months)
    REAL (KIND=wp), INTENT(OUT) :: hsurf_field(:,:,:) !< field for hsurf from ERA-I


    ! local variables
    INTEGER :: ndims  
    TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
    
    INTEGER, PARAMETER :: nglob_atts=6
    TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

    INTEGER :: errorcode !< error status variable

    INTEGER :: n !< counter





    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    PRINT *,'def_com_target_fields_meta'
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various T2M data related variables for netcdf output
    CALL def_era_meta(tg,ntime,dim_3d_tg)
    ! dim_sst_tg, sst_field_meta, hsurf_field_meta

    PRINT *,'CALL read netcdf data T2M'


    CALL netcdf_get_var(TRIM(netcdf_filename),t2m_field_meta,t2m_field)
    PRINT *,'ERA-I t2m_field read'

    CALL netcdf_get_var(TRIM(netcdf_filename),hsurf_field_meta,hsurf_field)
    PRINT *,'ERA-I hsurf_field read'



   END SUBROUTINE read_netcdf_buffer_t2m
   !-----------------------------------------------------------------





END Module mo_era_output_nc

