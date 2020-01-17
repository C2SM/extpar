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

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, &
                   i4

!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar

USE mo_grid_structures, ONLY: reg_lonlat_grid
                           
IMPLICIT NONE

PRIVATE

PUBLIC :: glcc_grid, &
 &         lon_glcc,  &
 &         lat_glcc,  &
 &         allocate_raw_glcc_fields,&
 &         deallocate_glcc_fields
          
TYPE(reg_lonlat_grid) :: glcc_grid !< structure with defenition of the raw data grid for the whole GLCC dataset

REAL (KIND=wp), ALLOCATABLE    :: lon_glcc(:) !< longitude of glcc raw data
REAL (KIND=wp), ALLOCATABLE    :: lat_glcc(:) !< latitude of glcc raw data


CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_glcc_fields(nrows,ncolumns)
  IMPLICIT NONE
  INTEGER (KIND=i4), INTENT(IN) :: nrows !< number of rows
  INTEGER (KIND=i4), INTENT(IN) :: ncolumns !< number of columns

  INTEGER :: errorcode !< error status variable


    ALLOCATE (lon_glcc(1:ncolumns), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lon_glcc')
    lon_glcc = 0.0

     ALLOCATE (lat_glcc(1:nrows), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lat_glcc')
    lat_glcc = 0.0

  END  SUBROUTINE allocate_raw_glcc_fields

  SUBROUTINE deallocate_glcc_fields()

    USE mo_glcc_tg_fields, ONLY : fr_land_glcc, &
                         &        glcc_class_fraction,    &
                         &        glcc_class_npixel, &
                         &        glcc_tot_npixel, &
                         &        ice_glcc, &
                         &        z0_glcc, &
                         &        root_glcc, &
                         &        plcov_mn_glcc, &
                         &        plcov_mx_glcc, &
                         &        lai_mn_glcc, &
                         &        lai_mx_glcc, &
                         &        rs_min_glcc, &
                         &        urban_glcc,  &
                         &        for_d_glcc,  &
                         &        for_e_glcc, &
                         &        emissivity_glcc

    IMPLICIT NONE


    INTEGER :: errorcode
    
    
    DEALLOCATE (lat_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lat_glcc')
    DEALLOCATE (lon_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lon_glcc')

    DEALLOCATE (fr_land_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector fr_land_glcc')
    DEALLOCATE (ice_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector ice_glcc')
    DEALLOCATE (z0_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector z0_glcc')
    DEALLOCATE (root_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector root_glcc')
    DEALLOCATE (plcov_mn_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector plcov_mn_glcc')
    DEALLOCATE (plcov_mx_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector plcov_mx_glcc')
    DEALLOCATE (lai_mn_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lai_mn_glcc')
    DEALLOCATE (lai_mx_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lai_mx_glcc')
    DEALLOCATE (rs_min_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector rs_min_glcc')
    DEALLOCATE (urban_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector urban_glcc')
    DEALLOCATE (for_d_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector for_d_glcc')
    DEALLOCATE (for_e_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector for_e_glcc')
    DEALLOCATE (emissivity_glcc, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector emissivity_glcc')
    DEALLOCATE (glcc_class_fraction, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector glcc_class_fraction')
    DEALLOCATE (glcc_class_npixel, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector glcc_class_npixel')
    DEALLOCATE (glcc_tot_npixel, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector glcc_tot_npixel')
    
  END SUBROUTINE deallocate_glcc_fields


END MODULE mo_glcc_data
