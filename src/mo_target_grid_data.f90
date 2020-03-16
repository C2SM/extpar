!+ Fortran module for data specification of target grid for external parameters
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_7         2013/01/25 Guenther Zaengl 
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON  
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for data specification of target grid for external parameters
!> \author Hermann Asensio
MODULE mo_target_grid_data

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lon_geo, &
       &    lat_geo, &
       &    no_raw_data_pixel, &
       &    tg, &
       &    allocate_com_target_fields, &
       &    search_res

  REAL (KIND=wp), ALLOCATABLE  :: lon_geo(:,:,:), &    !< longitude coordinates of the target grid in the geographical system 
       &                          lat_geo(:,:,:)          !< latitude coordinates of the target grid in the geographical system

  TYPE(target_grid_def)        :: tg !< structure with target grid description

  INTEGER(KIND=i4), PARAMETER  :: search_res = 4 ! resolution of ICON search index list is 1/search_res in units of degrees

  INTEGER(KIND=i4),ALLOCATABLE :: no_raw_data_pixel(:,:,:) !< number of raw data pixel inside the target grid element

  CONTAINS

  !> allocate common fields for target grid
  SUBROUTINE allocate_com_target_fields(tg)

    TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description 

    !local variables
    INTEGER(KIND=i4)                  :: errorcode
     
    CALL logging%info('Enter routine: allocate_com_target_fields')

    ALLOCATE (lon_geo(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lon_geo',__FILE__,__LINE__)
    lon_geo = 0.0

    ALLOCATE (lat_geo(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lat_geo',__FILE__,__LINE__)
    lat_geo = 0.0


    ALLOCATE (no_raw_data_pixel(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array no_raw_data_pixel',__FILE__,__LINE__)
    no_raw_data_pixel = 0

    CALL logging%info('Exit routine: allocate_com_target_fields')

  END SUBROUTINE  allocate_com_target_fields

END MODULE mo_target_grid_data

