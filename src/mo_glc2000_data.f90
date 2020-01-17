!+ Fortran Module with data fields for the GLC2000 data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran Module with data fields for the GLC2000 data
!> \author Hermann Asensio
!!
MODULE mo_glc2000_data

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, &
                   i4

!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar

USE mo_grid_structures, ONLY: reg_lonlat_grid
                           
IMPLICIT NONE

PRIVATE

PUBLIC :: glc2000_grid, &
          lon_glc2000,  &
          lat_glc2000,  &
          allocate_raw_glc2000_fields, &
          deallocate_glc2000_fields
          



TYPE(reg_lonlat_grid) :: glc2000_grid !< structure with defenition of the raw data grid for the whole GLC2000 dataset

REAL (KIND=wp), ALLOCATABLE    :: lon_glc2000(:) !< longitude of glc2000 raw data
REAL (KIND=wp), ALLOCATABLE    :: lat_glc2000(:) !< latitude of glc2000 raw data


CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_glc2000_fields(nrows,ncolumns)
  IMPLICIT NONE
  INTEGER (KIND=i4), INTENT(IN) :: nrows !< number of rows
  INTEGER (KIND=i4), INTENT(IN) :: ncolumns !< number of columns

  INTEGER :: errorcode !< error status variable


    ALLOCATE (lon_glc2000(1:ncolumns), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lon_glc2000')
    lon_glc2000 = 0.0

     ALLOCATE (lat_glc2000(1:nrows), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lat_glc2000')
    lat_glc2000 = 0.0



  END  SUBROUTINE allocate_raw_glc2000_fields

  SUBROUTINE deallocate_glc2000_fields()

 
  USE mo_lu_tg_fields, ONLY: fr_land_lu,       &
  &                          ice_lu,           &
  &                          z0_lu,            &
  &                          root_lu,          &
  &                          plcov_mn_lu,      &
  &                          plcov_mx_lu,      &
  &                          lai_mn_lu,        &
  &                          lai_mx_lu,        &
  &                          rs_min_lu,        &
  &                          urban_lu,         &
  &                          for_d_lu,         &
  &                          for_e_lu,         &
  &                          emissivity_lu,    &
  &                          lu_class_fraction,&
  &                          lu_class_npixel,  &
  &                          lu_tot_npixel                            

     IMPLICIT NONE


     INTEGER :: errorcode


       DEALLOCATE (lat_glc2000, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lat_glc2000')
       DEALLOCATE (lon_glc2000, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lon_glc2000')
       DEALLOCATE (ice_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector ice_lu')
       DEALLOCATE (z0_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector z0_lu')
       DEALLOCATE (root_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector root_lu')
       DEALLOCATE (plcov_mn_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector plcov_mn_lu')
       DEALLOCATE (plcov_mx_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector plcov_mx_lu')
       DEALLOCATE (lai_mn_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lai_mn_lu')
       DEALLOCATE (lai_mx_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lai_mx_lu')
       DEALLOCATE (rs_min_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector rs_min_lu')
       DEALLOCATE (urban_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector urban_lu')
       DEALLOCATE (for_d_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector for_d_lu')
       DEALLOCATE (for_e_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector for_e_lu')
       DEALLOCATE (emissivity_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector emissivity_lu')
       DEALLOCATE (lu_class_fraction, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector glc2000_class_fraction')
       DEALLOCATE (lu_class_npixel, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector glc2000_class_npixel')
       DEALLOCATE (lu_tot_npixel, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector glc2000_tot_npixel')
       DEALLOCATE (fr_land_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector fr_land_glc2000')


     END SUBROUTINE deallocate_glc2000_fields

END MODULE mo_glc2000_data
