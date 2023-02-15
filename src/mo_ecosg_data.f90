!+ Fortran Module with data fields for the ECOSG data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2023/02/05 Andrzej Wyszogrodzki
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran Module with data fields for the ECOSG data
!> \author Hermann Asensio
!!
MODULE mo_ecosg_data

  USE mo_logging
  USE mo_kind,                    ONLY: wp, i4

  USE mo_grid_structures,         ONLY: reg_lonlat_grid

  USE mo_ecosg_tg_fields,          ONLY: fr_land_ecosg, &
       &                                ecosg_class_fraction,    &
       &                                ecosg_class_npixel, &
       &                                ecosg_tot_npixel, &
       &                                ice_ecosg, &
       &                                z0_ecosg, &
       &                                root_ecosg, &
       &                                plcov_mn_ecosg, &
       &                                plcov_mx_ecosg, &
       &                                lai_mn_ecosg, &
       &                                lai_mx_ecosg, &
       &                                rs_min_ecosg, &
       &                                urban_ecosg,  &
       &                                for_d_ecosg,  &
       &                                for_e_ecosg, &
       &                                skinc_ecosg, &
       &                                emissivity_ecosg
                           
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ecosg_grid, &
   &         lon_ecosg,  &
   &         lat_ecosg,  &
   &         allocate_raw_ecosg_fields,&
   &         deallocate_ecosg_fields
            
  TYPE(reg_lonlat_grid)          :: ecosg_grid !< structure with defenition of the raw data grid for the whole ECOSG dataset

  REAL (KIND=wp), ALLOCATABLE    :: lon_ecosg(:), & !< longitude of ecosg raw data
       &                            lat_ecosg(:) !< latitude of ecosg raw data

  CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_ecosg_fields(nrows,ncolumns)

  IMPLICIT NONE

  INTEGER (KIND=i4), INTENT(IN) :: nrows, & !< number of rows
       &                           ncolumns !< number of columns

  INTEGER(KIND=i4)              :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_raw_ecosg_fields')

    ALLOCATE (lon_ecosg(1:ncolumns), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lon_ecosg',__FILE__,__LINE__)
    lon_ecosg = 0.0

    ALLOCATE (lat_ecosg(1:nrows), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lat_ecosg',__FILE__,__LINE__)
    lat_ecosg = 0.0

    CALL logging%info('Exit routine: allocate_raw_ecosg_fields')

  END  SUBROUTINE allocate_raw_ecosg_fields

  SUBROUTINE deallocate_ecosg_fields()

    IMPLICIT NONE

    INTEGER(KIND=i4) :: errorcode
    
    CALL logging%info('Enter routine: deallocate_ecosg_fields')

    DEALLOCATE (lat_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lat_ecosg',__FILE__,__LINE__)
    DEALLOCATE (lon_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lon_ecosg',__FILE__,__LINE__)

    DEALLOCATE (fr_land_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector fr_land_ecosg',__FILE__,__LINE__)
    DEALLOCATE (ice_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector ice_ecosg',__FILE__,__LINE__)
    DEALLOCATE (z0_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector z0_ecosg',__FILE__,__LINE__)
    DEALLOCATE (root_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector root_ecosg',__FILE__,__LINE__)
    DEALLOCATE (plcov_mn_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector plcov_mn_ecosg',__FILE__,__LINE__)
    DEALLOCATE (plcov_mx_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector plcov_mx_ecosg',__FILE__,__LINE__)
    DEALLOCATE (lai_mn_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lai_mn_ecosg',__FILE__,__LINE__)
    DEALLOCATE (lai_mx_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lai_mx_ecosg',__FILE__,__LINE__)
    DEALLOCATE (rs_min_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector rs_min_ecosg',__FILE__,__LINE__)
    DEALLOCATE (urban_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector urban_ecosg',__FILE__,__LINE__)
    DEALLOCATE (for_d_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector for_d_ecosg',__FILE__,__LINE__)
    DEALLOCATE (for_e_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector for_e_ecosg',__FILE__,__LINE__)
    DEALLOCATE (skinc_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector skinc_ecosg',__FILE__,__LINE__)
    DEALLOCATE (emissivity_ecosg, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector emissivity_ecosg',__FILE__,__LINE__)
    DEALLOCATE (ecosg_class_fraction, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector ecosg_class_fraction',__FILE__,__LINE__)
    DEALLOCATE (ecosg_class_npixel, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector ecosg_class_npixel',__FILE__,__LINE__)
    DEALLOCATE (ecosg_tot_npixel, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector ecosg_tot_npixel',__FILE__,__LINE__)

    CALL logging%info('Exit routine: deallocate_ecosg_fields')
    
  END SUBROUTINE deallocate_ecosg_fields

END MODULE mo_ecosg_data
