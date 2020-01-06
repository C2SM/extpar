!+ Fortran module with netcdf routines for land_sea_mask on target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_11        2013/04/16 Juergen Helmert
!  Adaptions for using external land-sea mask  
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with netcdf output routines for FLAKE data on the target grid
!> ouptut routines
!> \author Hermann Asensio
MODULE mo_lsm_output_nc

!wgrib -s lm_d5_07000_965x773.sso.mol.albproz_api_1_9_18.stf | 
!egrep "(:FR_LAND)" | wgrib -i -grib  lm_d5_07000_965x773.sso.mol.albproz_api_1_9_18.stf 
!-o lm_d5_07000_965x773_fr_land.g1

!cdo -f nc copy lm_d5_07000_965x773_fr_land.g1 lm_d5_07000_965x773_fr_land.nc
!replace lon and lat from buffer nc lona and lat fields using ncdump and ncgen

 

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

  USE mo_io_utilities, ONLY: vartype_int 
  USE mo_io_utilities, ONLY: vartype_real
  USE mo_io_utilities, ONLY: vartype_char

  USE mo_io_utilities, ONLY: netcdf_put_var
  USE mo_io_utilities, ONLY: open_new_netcdf_file
  USE mo_io_utilities, ONLY: close_netcdf_file
  USE mo_io_utilities, ONLY: netcdf_def_grid_mapping

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_netcdf_buffer_lsm


  CONTAINS



  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  !> read Land_Sea_Mask derived from ROUTI field
    SUBROUTINE read_netcdf_buffer_lsm(netcdf_filename,  &
    &                                     tg,         &
    &                                     fr_land_mask )


  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: def_com_target_fields_meta  
  
  USE mo_var_meta_data, ONLY: def_lsm_fields_meta

  USE mo_var_meta_data, ONLY: fr_land_mask_meta



  !USE mo_io_utilities, ONLY: netcdf_get_var_real_3d, netcdf_get_var_real_4d
  !USE mo_io_utilities, ONLY: netcdf_get_var_int_3d, netcdf_get_var_int_4d

  USE mo_io_utilities, ONLY: netcdf_get_var

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  REAL (KIND=wp), INTENT(OUT)  :: fr_land_mask(:,:,:)     !< fraction of fresh water (lakes)

  PRINT *,'ENTER read_netcdf_buffer_land_sea_mask'


  !set up dimensions for buffer
  !CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

  
  ! define meta information for various land use related variables (FLAKE) for netcdf output
  CALL def_lsm_fields_meta(dim_3d_tg)

 ! fr_land_mask_meta


  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta

  CALL netcdf_get_var(TRIM(netcdf_filename),fr_land_mask_meta,fr_land_mask)
  PRINT *,'Land-Sea-Mask from external file read'


  END SUBROUTINE read_netcdf_buffer_lsm
  !-----------------------------------------------------------------------

    
 
END Module mo_lsm_output_nc

