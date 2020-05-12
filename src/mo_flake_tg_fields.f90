!+ Fortran module for FLake data specification on target grid for external Parameters
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
!> Fortran module for FLake data specification on target grid for external Parameters
!> \author Hermann Asensio
!
MODULE mo_flake_tg_fields

  !> kind parameters are defined in MODULE data_parameters
  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4
  USE mo_array_cache,           ONLY: allocate_cached
  USE mo_grid_structures,       ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lake_depth, &
      &     fr_lake,    &
      &     flake_tot_npixel, &
      &     allocate_flake_target_fields

  INTEGER (KIND=i4), POINTER :: flake_tot_npixel(:,:,:) !< total number of flake raw data pixels on target grid
  REAL (KIND=wp), POINTER  :: lake_depth(:,:,:), &!< lake depth [m]
       &                      fr_lake(:,:,:)      !< lake fraction

  CONTAINS

!! the target grid for the GME has 3 dimension (ie,je,jd),
!! the target grid for the COSMO model has 2 dimension (ie,je)
!! the target grid for the ICON model has 1 dimension (ne)
!! depending of the target model the second and third dimension of the target fields should be
!! allocated with the length 1
  SUBROUTINE allocate_flake_target_fields(tg, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    LOGICAL, INTENT(in)               :: l_use_array_cache      

    INTEGER                           :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_flake_target_fields')

if (l_use_array_cache) then
   call allocate_cached('lake_depth', lake_depth, [tg%ie,tg%je,tg%ke])
else
   allocate(lake_depth(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lake_depth',__FILE__,__LINE__)
    lake_depth = 0.0

if (l_use_array_cache) then
   call allocate_cached('fr_lake', fr_lake, [tg%ie,tg%je,tg%ke])
else
   allocate(fr_lake(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_lake',__FILE__,__LINE__)
    fr_lake = 0.0

if (l_use_array_cache) then
   call allocate_cached('flake_tot_npixel', flake_tot_npixel, [tg%ie,tg%je,tg%ke])
else
   allocate(flake_tot_npixel(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array flake_tot_npixel',__FILE__,__LINE__)
    flake_tot_npixel = 0

  END SUBROUTINE allocate_flake_target_fields

END MODULE mo_flake_tg_fields

