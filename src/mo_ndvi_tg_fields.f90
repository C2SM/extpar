!+ Fortran module for NDVI data on target grid for external parameters
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
!> Fortran module for NDVI data on target grid for external parameters 
!> \author Hermann Asensio
MODULE mo_ndvi_tg_fields

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_grid_structures, ONLY: target_grid_def


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ndvi_field, &
    &        ndvi_max, &
    &        ndvi_field_mom, &
    &        ndvi_ratio_mom, &
    &        allocate_ndvi_target_fields


         REAL(KIND=wp), ALLOCATABLE  :: ndvi_field(:,:,:) !< field for ndvi data

         REAL(KIND=wp), ALLOCATABLE :: ndvi_max(:,:,:) !< field for ndvi maximum

         REAL(KIND=wp), ALLOCATABLE  :: ndvi_field_mom(:,:,:,:) !< field for monthly mean ndvi data (12 months)

         REAL(KIND=wp), ALLOCATABLE  :: ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)


  CONTAINS

  !> allocate fields for GLOBE target data 
    SUBROUTINE allocate_ndvi_target_fields(tg,nt)
      IMPLICIT NONE

      TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
      INTEGER (KIND=i4), INTENT(in) :: nt !< number of timesteps (12 for monthly mean values)

      INTEGER :: errorcode !< error status variable
        
      ALLOCATE (ndvi_field(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
          IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array ndvi_field')
      ndvi_field = 0.0

       ALLOCATE (ndvi_max(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
          IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array ndvi_max')
      ndvi_max = 0.0

       ALLOCATE (ndvi_field_mom(1:tg%ie,1:tg%je,1:tg%ke,1:nt), STAT=errorcode)
          IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array ndvi_field_mom')
      ndvi_field_mom = 0.0

       ALLOCATE (ndvi_ratio_mom(1:tg%ie,1:tg%je,1:tg%ke,1:nt), STAT=errorcode)
          IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array ndvi_ratio_mom')
      ndvi_ratio_mom = 0.0

    END SUBROUTINE allocate_ndvi_target_fields

END Module mo_ndvi_tg_fields

