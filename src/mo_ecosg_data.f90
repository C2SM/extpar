MODULE mo_ecosg_data

  USE mo_logging
  USE mo_kind,                    ONLY: wp, i4

  USE mo_grid_structures,         ONLY: reg_lonlat_grid

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ecosg_grid, &
   &         lon_ecosg,  &
   &         lat_ecosg,  &
   &         allocate_raw_ecosg_fields

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

END MODULE mo_ecosg_data
