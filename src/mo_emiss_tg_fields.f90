!+ Fortran module for EMISS data on target grid for external parameters
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
!> Fortran module for EMISS data on target grid for external parameters 
!> \author Hermann Asensio
MODULE mo_emiss_tg_fields

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def


  IMPLICIT NONE

  PRIVATE

  PUBLIC ::  emiss_field, &
    &        emiss_max, &
    &        emiss_field_mom, &
    &        emiss_ratio_mom, &
    &        allocate_emiss_target_fields


  REAL(KIND=wp), ALLOCATABLE  :: emiss_field(:,:,:), & !< field for emiss data
       &                         emiss_max(:,:,:), & !< field for emiss maximum
       &                         emiss_field_mom(:,:,:,:), & !< field for monthly mean emiss data (12 months)
       &                         emiss_ratio_mom(:,:,:,:) !< field for monthly emiss ratio (12 months)

  CONTAINS

  !> allocate fields for GLOBE target data 
  SUBROUTINE allocate_emiss_target_fields(tg,nt)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(in)     :: nt !< number of timesteps (12 for monthly mean values)

    INTEGER(KIND=i4)                  :: errorcode !< error status variable
      
    ALLOCATE (emiss_field(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emiss_field',__FILE__,__LINE__)
    emiss_field = 0.0

    ALLOCATE (emiss_max(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emiss_max',__FILE__,__LINE__)
    emiss_max = 0.0

    ALLOCATE (emiss_field_mom(1:tg%ie,1:tg%je,1:tg%ke,1:nt), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emiss_field_mom',__FILE__,__LINE__)
    emiss_field_mom = 0.0

    ALLOCATE (emiss_ratio_mom(1:tg%ie,1:tg%je,1:tg%ke,1:nt), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emiss_ratio_mom',__FILE__,__LINE__)
    emiss_ratio_mom = 0.0

  END SUBROUTINE allocate_emiss_target_fields

END MODULE mo_emiss_tg_fields
