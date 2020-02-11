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
  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  USE mo_var_meta_data,         ONLY: dim_3d_tg, &
       &                              def_dimension_info_buffer, &
       &                              def_com_target_fields_meta, &  
       &                              sst_field_meta, &
       &                              wsnow_field_meta,&
       &                              t2m_field_meta, & 
       &                              hsurf_field_meta,&
       &                              def_era_meta


  USE mo_io_utilities,          ONLY: netcdf_get_var

  IMPLICIT NONE 

  PRIVATE

  PUBLIC :: read_netcdf_buffer_sst, &
       &    read_netcdf_buffer_t2m

  CONTAINS

  SUBROUTINE read_netcdf_buffer_sst(netcdf_filename,  &
       &                            tg,         &
       &                            ntime, &
       &                            sst_field,&
       &                            wsnow_field)



    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(INOUT)   :: ntime !< number of times of sst data (12 monthly mean values)

    REAL (KIND=wp), INTENT(OUT)        :: sst_field(:,:,:,:), & !< field for monthly mean sst data (12 months)
         &                                wsnow_field(:,:,:,:) !< field for monthly sst ratio (12 months)

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6
    
    CALL logging%info('Enter routine: read_netcdf_buffer_sst')
    
    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various SST data related variables for netcdf output
    CALL def_era_meta(ntime,dim_3d_tg)
    ! dim_sst_tg, sst_field_meta, wsnow_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),sst_field_meta,sst_field)

    CALL netcdf_get_var(TRIM(netcdf_filename),wsnow_field_meta,wsnow_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_sst')

  END SUBROUTINE read_netcdf_buffer_sst

  SUBROUTINE read_netcdf_buffer_t2m(netcdf_filename,  &
       &                            tg,         &
       &                            ntime, &
       &                            t2m_field,&
       &                            hsurf_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(INOUT)   :: ntime !< number of times of t2m data (12 monthly mean values)

    REAL (KIND=wp), INTENT(OUT)        :: t2m_field(:,:,:,:), & !< field for monthly mean t2m data (12 months)
         &                                hsurf_field(:,:,:) !< field for hsurf from ERA-I

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6
    
    CALL logging%info('Enter routine: read_netcdf_buffer_t2m')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various T2M data related variables for netcdf output
    CALL def_era_meta(ntime,dim_3d_tg)
    ! dim_sst_tg, sst_field_meta, hsurf_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),t2m_field_meta,t2m_field)

    CALL netcdf_get_var(TRIM(netcdf_filename),hsurf_field_meta,hsurf_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_t2m')
    
  END SUBROUTINE read_netcdf_buffer_t2m

END MODULE mo_era_output_nc
