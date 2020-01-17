!+ Fortran Module with data fields for the globcover data
!
! Current Code Owner: DWD, Hermann Asensio
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_3         2011/04/19 Hermann Asensio
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!! ECOCLIMAP option gs_08.03.12  programm not operational yet
MODULE mo_ecoclimap_data

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, &
                   i4, &
                   i4

!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar

USE mo_grid_structures, ONLY: reg_lonlat_grid

IMPLICIT NONE

PRIVATE

PUBLIC :: ecoclimap_grid, &
          lon_ecoclimap,  &
          lat_ecoclimap,  &
          ntime_ecoclimap, &
          allocate_raw_ecoclimap_fields, &
          deallocate_ecoclimap_fields

TYPE(reg_lonlat_grid) :: ecoclimap_grid !< structure with defenition of the raw data grid for the whole globcover dataset

REAL (KIND=wp), ALLOCATABLE    :: lon_ecoclimap(:) !< longitude of ecoclimap raw data
REAL (KIND=wp), ALLOCATABLE    :: lat_ecoclimap(:) !< latitude of ecoclimap raw data
INTEGER (KIND=i4)              :: ntime_ecoclimap = 12 !< number of timesteps (12 for monthly mean values)

CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_ecoclimap_fields(nrows,ncolumns)
  IMPLICIT NONE
  INTEGER (KIND=i4), INTENT(IN) :: nrows !< number of rows
  INTEGER (KIND=i4), INTENT(IN) :: ncolumns !< number of columns

  INTEGER :: errorcode !< error status variable

    ALLOCATE (lon_ecoclimap(1:ncolumns), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lon_ecoclimap')
    lon_ecoclimap = 0.0

     ALLOCATE (lat_ecoclimap(1:nrows), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lat_ecoclimap')
    lat_ecoclimap = 0.0

  END  SUBROUTINE allocate_raw_ecoclimap_fields

  SUBROUTINE deallocate_ecoclimap_fields()

    IMPLICIT NONE

    INTEGER :: errorcode

    IF (ALLOCATED(lat_ecoclimap)) DEALLOCATE (lat_ecoclimap, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lat_ecoclimap')
    IF (ALLOCATED(lon_ecoclimap)) DEALLOCATE (lon_ecoclimap, STAT = errorcode)
    IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lon_ecoclimap')

  END SUBROUTINE deallocate_ecoclimap_fields
  


END MODULE mo_ecoclimap_data

