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
MODULE mo_era_sst_output_nc
!
  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: reg_lonlat_grid, &
       &                              rotated_lonlat_grid, &
       &                              icosahedral_triangular_grid, &
       &                              target_grid_def

  USE mo_var_meta_data,         ONLY: dim_3d_tg, dim_icon, &
       &                              def_dimension_info_buffer, &
       &                              lon_geo_meta, &
       &                              lat_geo_meta, &
       &                              no_raw_data_pixel_meta, &
       &                              def_com_target_fields_meta, &  
       &                              dim_era_tg, &
       &                              sst_field_meta, &
       &                              wsnow_field_meta,&
       &                              def_era_meta, &
       &                              t2m_field_meta, &
       &                              hsurf_field_meta,&
       &                              def_era_meta

  USE mo_io_utilities,          ONLY: var_meta_info, &
       &                              netcdf_attributes, &
       &                              dim_meta_info, &
       &                              netcdf_put_var, &
       &                              open_new_netcdf_file, &
       &                              close_netcdf_file, &
       &                              netcdf_def_grid_mapping, &
       &                              get_date_const_field, &
       &                              set_date_mm_extpar_field, &
       &                              vartype_int, & 
       &                              vartype_real, &
       &                              netcdf_get_var, &
       &                              vartype_char

  IMPLICIT NONE 

  PRIVATE

  PUBLIC :: read_netcdf_buffer_sst, &
       &    read_netcdf_buffer_t2m

  CONTAINS

  !> set global attributes for netcdf with SST data
  SUBROUTINE set_global_att_sst(global_attributes)

    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:6)

    !local variables
    CHARACTER(len=10) :: ydate, ytime
    CHARACTER(len=2)  :: cc, &
         &               yy, &
         &               mm, &
         &               dd, &
         &               hh, &
         &               minute

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

  SUBROUTINE read_netcdf_buffer_sst(netcdf_filename,  &
   &                                     tg,         &
   &                                     ntime, &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     sst_field,&
   &                                     wsnow_field)



    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(INOUT)   :: ntime !< number of times of sst data (12 monthly mean values)
    REAL(KIND=wp), INTENT(OUT)         :: undefined       !< value to indicate undefined grid elements 
    INTEGER, INTENT(OUT)               :: undef_int       !< value to indicate undefined grid elements

    REAL (KIND=wp), INTENT(OUT)        :: sst_field(:,:,:,:), & !< field for monthly mean sst data (12 months)
         &                                wsnow_field(:,:,:,:) !< field for monthly sst ratio (12 months)

    ! local variables
    INTEGER(KIND=i4)                   :: ndims, errorcode, n
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    TYPE(dim_meta_info), ALLOCATABLE   :: dim_list(:) !< dimensions for netcdf file
    TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various SST data related variables for netcdf output
    CALL def_era_meta(tg,ntime,dim_3d_tg)
    ! dim_sst_tg, sst_field_meta, wsnow_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),sst_field_meta,sst_field)

    CALL netcdf_get_var(TRIM(netcdf_filename),wsnow_field_meta,wsnow_field)

  END SUBROUTINE read_netcdf_buffer_sst

  SUBROUTINE read_netcdf_buffer_t2m(netcdf_filename,  &
   &                                     tg,         &
   &                                     ntime, &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     t2m_field,&
   &                                     hsurf_field)



    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(INOUT)   :: ntime !< number of times of t2m data (12 monthly mean values)
    REAL(KIND=wp), INTENT(OUT)         :: undefined       !< value to indicate undefined grid elements 
    INTEGER, INTENT(OUT)               :: undef_int       !< value to indicate undefined grid elements

    REAL (KIND=wp), INTENT(OUT)        :: t2m_field(:,:,:,:), & !< field for monthly mean t2m data (12 months)
         &                                hsurf_field(:,:,:) !< field for hsurf from ERA-I


    ! local variables
    INTEGER(KIND=i4)                   :: ndims, errorcode, n  
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    TYPE(dim_meta_info), ALLOCATABLE   :: dim_list(:) !< dimensions for netcdf file
    TYPE(netcdf_attributes)            :: global_attributes(nglob_atts)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various T2M data related variables for netcdf output
    CALL def_era_meta(tg,ntime,dim_3d_tg)
    ! dim_sst_tg, sst_field_meta, hsurf_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),t2m_field_meta,t2m_field)

    CALL netcdf_get_var(TRIM(netcdf_filename),hsurf_field_meta,hsurf_field)

   END SUBROUTINE read_netcdf_buffer_t2m

END MODULE mo_era_sst_output_nc
