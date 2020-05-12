!+ Fortran module for AHF data on target grid for external parameters
!
!
! Description:
! Fortran module for AHF data on target grid for external parameters
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_4         2015-05-21 Hendrik Wouters
!  Initial release
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!> Fortran module for AHF data on target grid for external parameters
!> \author Hermann Asensio
MODULE mo_ahf_tg_fields

  USE mo_logging
  USE mo_kind,                  ONLY: wp
  USE mo_array_cache,           ONLY: allocate_cached
  USE mo_grid_structures,       ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ahf_field, &
       &    allocate_ahf_target_fields

  REAL(KIND=wp), POINTER :: ahf_field(:,:,:) !< field for ahf data

  CONTAINS

  !> allocate fields for GLOBE target data
  SUBROUTINE allocate_ahf_target_fields(tg, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    LOGICAL, INTENT(in)               :: l_use_array_cache
    INTEGER                           :: errorcode !< error status variable

if (l_use_array_cache) then
   call allocate_cached('ahf_field', ahf_field, [tg%ie,tg%je,tg%ke])
else
   allocate(ahf_field(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ahf_field',__FILE__,__LINE__)
    ahf_field = 0.0

  END SUBROUTINE allocate_ahf_target_fields

END Module mo_ahf_tg_fields
