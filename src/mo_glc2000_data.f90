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

  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4

  USE mo_grid_structures,       ONLY: reg_lonlat_grid
  
  USE mo_lu_tg_fields,          ONLY: fr_land_lu,       &
       &                              ice_lu,           &
       &                              z0_lu,            &
       &                              root_lu,          &
       &                              plcov_mn_lu,      &
       &                              plcov_mx_lu,      &
       &                              lai_mn_lu,        &
       &                              lai_mx_lu,        &
       &                              rs_min_lu,        &
       &                              urban_lu,         &
       &                              for_d_lu,         &
       &                              for_e_lu,         &
       &                              emissivity_lu,    &
       &                              lu_class_fraction,&
       &                              lu_class_npixel,  &
       &                              lu_tot_npixel                            
                           
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: glc2000_grid, &
       &    lon_glc2000,  &
       &    lat_glc2000,  &
       &    allocate_raw_glc2000_fields, &
       &    deallocate_glc2000_fields

  TYPE(reg_lonlat_grid)          :: glc2000_grid !< structure with defenition of the raw data grid for the whole GLC2000 dataset

  REAL (KIND=wp), ALLOCATABLE    :: lon_glc2000(:), &  !< longitude of glc2000 raw data
       &                            lat_glc2000(:) !< latitude of glc2000 raw data

  CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_glc2000_fields(nrows,ncolumns)

    IMPLICIT NONE

    INTEGER (KIND=i4), INTENT(IN) :: nrows, &  !< number of rows
         &                           ncolumns !< number of columns

    INTEGER(KIND=i4)              :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_raw_glc2000_fields')

    ALLOCATE (lon_glc2000(1:ncolumns), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lon_glc2000',__FILE__,__LINE__)
    lon_glc2000 = 0.0

    ALLOCATE (lat_glc2000(1:nrows), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lat_glc2000',__FILE__,__LINE__)
    lat_glc2000 = 0.0

  END  SUBROUTINE allocate_raw_glc2000_fields

  SUBROUTINE deallocate_glc2000_fields()

    IMPLICIT NONE

    INTEGER(KIND=i4) :: errorcode

    CALL logging%info('Enter routine: deallocate_glc2000_fields')

    DEALLOCATE (lat_glc2000, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lat_glc2000',__FILE__,__LINE__)
    DEALLOCATE (lon_glc2000, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lon_glc2000',__FILE__,__LINE__)
    DEALLOCATE (ice_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector ice_lu',__FILE__,__LINE__)
    DEALLOCATE (z0_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector z0_lu',__FILE__,__LINE__)
    DEALLOCATE (root_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector root_lu',__FILE__,__LINE__)
    DEALLOCATE (plcov_mn_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector plcov_mn_lu',__FILE__,__LINE__)
    DEALLOCATE (plcov_mx_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector plcov_mx_lu',__FILE__,__LINE__)
    DEALLOCATE (lai_mn_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lai_mn_lu',__FILE__,__LINE__)
    DEALLOCATE (lai_mx_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lai_mx_lu',__FILE__,__LINE__)
    DEALLOCATE (rs_min_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector rs_min_lu',__FILE__,__LINE__)
    DEALLOCATE (urban_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector urban_lu',__FILE__,__LINE__)
    DEALLOCATE (for_d_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector for_d_lu',__FILE__,__LINE__)
    DEALLOCATE (for_e_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector for_e_lu',__FILE__,__LINE__)
    DEALLOCATE (emissivity_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector emissivity_lu',__FILE__,__LINE__)
    DEALLOCATE (lu_class_fraction, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector glc2000_class_fraction',__FILE__,__LINE__)
    DEALLOCATE (lu_class_npixel, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector glc2000_class_npixel',__FILE__,__LINE__)
    DEALLOCATE (lu_tot_npixel, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector glc2000_tot_npixel',__FILE__,__LINE__)
    DEALLOCATE (fr_land_lu, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector fr_land_glc2000',__FILE__,__LINE__)

  END SUBROUTINE deallocate_glc2000_fields

END MODULE mo_glc2000_data
