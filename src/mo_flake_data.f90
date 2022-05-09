!+ Fortran Module with data fields for the flake data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_2         2011/03/25 Hermann Asensio
!  Add parameter DWD_min_lake_depth for minimal lake depth         
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran Module with data fields for the flake data
!> \author Hermann Asensio
!!
MODULE mo_flake_data

  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4

  USE mo_grid_structures,       ONLY: reg_lonlat_grid
  USE mo_flake_tg_fields,       ONLY: lake_depth, &
       &                              fr_lake,    &
       &                              flake_tot_npixel
                           
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: flake_grid, &
     &      lon_flake,  &
     &      lat_flake,  &
     &      allocate_raw_flake_fields, &
     &      deallocate_raw_flake_fields

  PUBLIC :: flake_depth_undef, flake_depth_default, DWD_max_lake_depth, DWD_min_lake_depth
          
  TYPE(reg_lonlat_grid) :: flake_grid !< structure with defenition of the raw data grid for the whole GLCC dataset

  REAL (KIND=wp), ALLOCATABLE    :: lon_flake(:), & !< longitude of flake raw data
       &                            lat_flake(:) !< latitude of flake raw data

  REAL (KIND=wp), PARAMETER      :: flake_depth_undef = -1., & !< default value for undefined lake depth
       &                            flake_depth_default = 10.0, & !< default value for default lake depth, 10 [m]
       &                            DWD_max_lake_depth = 50.0, & !< Maximum lake depth in [m] for FLAKE
       &                            DWD_min_lake_depth = 1.0 !< Minimal lake depth in [m] for FLAKE


  CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_flake_fields(nrows,ncolumns)
    
    IMPLICIT NONE

    INTEGER (KIND=i4), INTENT(IN) :: nrows, & !< number of rows
         &                           ncolumns !< number of columns

    INTEGER                       :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_raw_flake_fields')
    ALLOCATE (lon_flake(1:ncolumns+1), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lon_flake',__FILE__,__LINE__)
    lon_flake = 0.0

    ALLOCATE (lat_flake(1:nrows), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lat_flake',__FILE__,__LINE__)
    lat_flake = 0.0

  END  SUBROUTINE allocate_raw_flake_fields

  SUBROUTINE deallocate_raw_flake_fields()


    IMPLICIT NONE

    INTEGER :: errorcode
    
    CALL logging%info('Enter routine: deallocate_raw_flake_fields')

    DEALLOCATE (lon_flake, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array lon_flake',__FILE__,__LINE__)
    DEALLOCATE (lat_flake, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array lat_flake',__FILE__,__LINE__)
    DEALLOCATE (fr_lake, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array fr_lake',__FILE__,__LINE__)
    DEALLOCATE (lake_depth, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array lake_depth',__FILE__,__LINE__)
    DEALLOCATE (flake_tot_npixel, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array flake_tot_npixel',__FILE__,__LINE__)
    
  END SUBROUTINE deallocate_raw_flake_fields

END MODULE mo_flake_data
