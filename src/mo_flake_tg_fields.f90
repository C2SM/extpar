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
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar
  USE mo_grid_structures, ONLY: target_grid_def


IMPLICIT NONE

PRIVATE

PUBLIC :: lake_depth, &
  &       fr_lake,    &
  &       flake_tot_npixel, &
  &       allocate_flake_target_fields

       INTEGER (KIND=i8), ALLOCATABLE :: flake_tot_npixel(:,:,:)  
                                         !< total number of flake raw data pixels on target grid (dimension (ie,je,ke))
       REAL (KIND=wp), ALLOCATABLE  :: lake_depth(:,:,:)     !< lake depth [m]
       REAL (KIND=wp), ALLOCATABLE  :: fr_lake(:,:,:)      !< lake fraction

CONTAINS


!> allocate fields for TARGET grid
!!
!! the target grid for the GME has 3 dimension (ie,je,jd),
!! the target grid for the COSMO model has 2 dimension (ie,je)
!! the target grid for the ICON model has 1 dimension (ne)
!! depending of the target model the second and third dimension of the target fields should be 
!! allocated with the length 1
  SUBROUTINE allocate_flake_target_fields(tg)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER :: errorcode !< error status variable
   
    ALLOCATE (lake_depth(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lake_depth')
    lake_depth = 0.0

     ALLOCATE (fr_lake(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_lake')
    fr_lake = 0.0

    ALLOCATE (flake_tot_npixel(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array flake_tot_npixel')
    flake_tot_npixel = 0


  END SUBROUTINE allocate_flake_target_fields



END MODULE mo_flake_tg_fields

