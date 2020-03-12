!+ Fortran Module with data fields for the GLCC data
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
!> Fortran Module with data fields for the GLCC data
!> \author Hermann Asensio
!!
MODULE mo_glcc_data

  USE mo_logging
  USE mo_kind,                    ONLY: wp, i4

  USE mo_grid_structures,         ONLY: reg_lonlat_grid

  USE mo_glcc_tg_fields,          ONLY: fr_land_glcc, &
       &                                glcc_class_fraction,    &
       &                                glcc_class_npixel, &
       &                                glcc_tot_npixel, &
       &                                ice_glcc, &
       &                                z0_glcc, &
       &                                root_glcc, &
       &                                plcov_mn_glcc, &
       &                                plcov_mx_glcc, &
       &                                lai_mn_glcc, &
       &                                lai_mx_glcc, &
       &                                rs_min_glcc, &
       &                                urban_glcc,  &
       &                                for_d_glcc,  &
       &                                for_e_glcc, &
       &                                emissivity_glcc
                           
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: glcc_grid, &
   &         lon_glcc,  &
   &         lat_glcc,  &
   &         allocate_raw_glcc_fields,&
   &         deallocate_glcc_fields
            
  TYPE(reg_lonlat_grid)          :: glcc_grid !< structure with defenition of the raw data grid for the whole GLCC dataset

  REAL (KIND=wp), ALLOCATABLE    :: lon_glcc(:), & !< longitude of glcc raw data
       &                            lat_glcc(:) !< latitude of glcc raw data

  CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_glcc_fields(nrows,ncolumns)

  IMPLICIT NONE

  INTEGER (KIND=i4), INTENT(IN) :: nrows, & !< number of rows
       &                           ncolumns !< number of columns

  INTEGER(KIND=i4)              :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_raw_glcc_fields')

    ALLOCATE (lon_glcc(1:ncolumns), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lon_glcc',__FILE__,__LINE__)
    lon_glcc = 0.0

    ALLOCATE (lat_glcc(1:nrows), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lat_glcc',__FILE__,__LINE__)
    lat_glcc = 0.0

    CALL logging%info('Exit routine: allocate_raw_glcc_fields')

  END  SUBROUTINE allocate_raw_glcc_fields

  SUBROUTINE deallocate_glcc_fields()

    IMPLICIT NONE

    INTEGER(KIND=i4) :: errorcode
    
    CALL logging%info('Enter routine: deallocate_glcc_fields')

    DEALLOCATE (lat_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lat_glcc',__FILE__,__LINE__)
    DEALLOCATE (lon_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lon_glcc',__FILE__,__LINE__)

    DEALLOCATE (fr_land_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector fr_land_glcc',__FILE__,__LINE__)
    DEALLOCATE (ice_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector ice_glcc',__FILE__,__LINE__)
    DEALLOCATE (z0_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector z0_glcc',__FILE__,__LINE__)
    DEALLOCATE (root_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector root_glcc',__FILE__,__LINE__)
    DEALLOCATE (plcov_mn_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector plcov_mn_glcc',__FILE__,__LINE__)
    DEALLOCATE (plcov_mx_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector plcov_mx_glcc',__FILE__,__LINE__)
    DEALLOCATE (lai_mn_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lai_mn_glcc',__FILE__,__LINE__)
    DEALLOCATE (lai_mx_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lai_mx_glcc',__FILE__,__LINE__)
    DEALLOCATE (rs_min_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector rs_min_glcc',__FILE__,__LINE__)
    DEALLOCATE (urban_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector urban_glcc',__FILE__,__LINE__)
    DEALLOCATE (for_d_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector for_d_glcc',__FILE__,__LINE__)
    DEALLOCATE (for_e_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector for_e_glcc',__FILE__,__LINE__)
    DEALLOCATE (emissivity_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector emissivity_glcc',__FILE__,__LINE__)
    DEALLOCATE (glcc_class_fraction, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector glcc_class_fraction',__FILE__,__LINE__)
    DEALLOCATE (glcc_class_npixel, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector glcc_class_npixel',__FILE__,__LINE__)
    DEALLOCATE (glcc_tot_npixel, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector glcc_tot_npixel',__FILE__,__LINE__)

    CALL logging%info('Exit routine: deallocate_glcc_fields')
    
  END SUBROUTINE deallocate_glcc_fields

END MODULE mo_glcc_data
