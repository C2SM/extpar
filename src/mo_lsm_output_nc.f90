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

  USE mo_logging
  USE mo_kind,                  ONLY: wp

  USE mo_var_meta_data,         ONLY: dim_3d_tg, &
       &                              def_com_target_fields_meta, &  
       &                              def_lsm_fields_meta, &
       &                              fr_land_mask_meta

  USE mo_io_utilities,          ONLY: netcdf_get_var

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_netcdf_buffer_lsm


  CONTAINS

  !> read Land_Sea_Mask derived from ROUTI field
    SUBROUTINE read_netcdf_buffer_lsm(netcdf_filename,  &
    &                                     fr_land_mask )

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  REAL (KIND=wp), INTENT(OUT)        :: fr_land_mask(:,:,:)     !< fraction of fresh water (lakes)

  CALL logging%info('Enter routine: read_netcdf_buffer_land_sea_mask')

  ! define meta information for various land use related variables (FLAKE) for netcdf output
  CALL def_lsm_fields_meta(dim_3d_tg)

  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta

  CALL netcdf_get_var(TRIM(netcdf_filename),fr_land_mask_meta,fr_land_mask)

  CALL logging%info('Exit routine: read_netcdf_buffer_land_sea_mask')

  END SUBROUTINE read_netcdf_buffer_lsm
                                                                          
END MODULE mo_lsm_output_nc
