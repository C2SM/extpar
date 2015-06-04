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

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4
  USE mo_kind, ONLY: i8

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_grid_structures, ONLY: target_grid_def


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ahf_field, &
    &        allocate_ahf_target_fields

         REAL(KIND=wp), ALLOCATABLE  :: ahf_field(:,:,:) !< field for ahf data


  CONTAINS

  !> allocate fields for GLOBE target data 
    SUBROUTINE allocate_ahf_target_fields(tg)
      IMPLICIT NONE

      TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

      INTEGER :: errorcode !< error status variable
        
      ALLOCATE (ahf_field(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
          IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array ahf_field')
      ahf_field = 0.0

    END SUBROUTINE allocate_ahf_target_fields


END Module mo_ahf_tg_fields

