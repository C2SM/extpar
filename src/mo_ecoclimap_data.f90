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

  USE mo_logging
  USE mo_kind, ONLY: wp, i4

  USE mo_grid_structures, ONLY: reg_lonlat_grid

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ecoclimap_grid, &
       &    lon_ecoclimap,  &
       &    lat_ecoclimap,  &
       &    ntime_ecoclimap, &
       &    allocate_raw_ecoclimap_fields, &
       &    deallocate_ecoclimap_fields

  TYPE(reg_lonlat_grid)          :: ecoclimap_grid !< structure with defenition of the raw data grid for the whole globcover dataset

  REAL (KIND=wp), ALLOCATABLE    :: lon_ecoclimap(:), &  !< longitude of ecoclimap raw data
       &                            lat_ecoclimap(:) !< latitude of ecoclimap raw data

  INTEGER (KIND=i4)              :: ntime_ecoclimap = 12 !< number of timesteps (12 for monthly mean values)

  CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_ecoclimap_fields(nrows,ncolumns)

    IMPLICIT NONE

    INTEGER (KIND=i4), INTENT(IN) :: nrows, &  !< number of rows
         &                           ncolumns !< number of columns

    INTEGER(KIND=i4)              :: errorcode !< error status variable

    CALL logging%info('allocate_raw_ecoclimap_fields')

    ALLOCATE (lon_ecoclimap(1:ncolumns), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lon_ecoclimap',__FILE__,__LINE__)
    lon_ecoclimap = 0.0

    ALLOCATE (lat_ecoclimap(1:nrows), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lat_ecoclimap',__FILE__,__LINE__)
    lat_ecoclimap = 0.0

  END  SUBROUTINE allocate_raw_ecoclimap_fields

  SUBROUTINE deallocate_ecoclimap_fields()

    IMPLICIT NONE

    INTEGER(KIND=i4) :: errorcode

    CALL logging%info('deallocate_ecoclimap_fields')

    IF (ALLOCATED(lat_ecoclimap)) DEALLOCATE (lat_ecoclimap, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lat_ecoclimap',__FILE__,__LINE__)
    IF (ALLOCATED(lon_ecoclimap)) DEALLOCATE (lon_ecoclimap, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lon_ecoclimap',__FILE__,__LINE__)

  END SUBROUTINE deallocate_ecoclimap_fields

END MODULE mo_ecoclimap_data
