!+ <A one line description of this module>
!
!
! Description:
! <Say what this module is for>
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_10        2018-02-19 Juergen Helmert 
!  lsubtract_mean_slope, ERA-I surface temp for land points
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!
MODULE mo_era_tg_fields


  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4
  USE mo_kind, ONLY: i8

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_grid_structures, ONLY: target_grid_def


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sst_field, &
    &        wsnow_field, &
    &        t2m_field, &
    &        hsurf_field, &
    &        allocate_era_target_fields


         REAL(KIND=wp), ALLOCATABLE  :: sst_field(:,:,:,:) !< field for sst data (12 months)
         REAL(KIND=wp), ALLOCATABLE  :: wsnow_field(:,:,:,:) !< field for wsnow data (12 months)
         REAL(KIND=wp), ALLOCATABLE  :: t2m_field(:,:,:,:) !< field for wsnow data (12 months)
         REAL(KIND=wp), ALLOCATABLE  :: hsurf_field(:,:,:) !< field for wsnow data (12 months)


  CONTAINS

  !> allocate fields for GLOBE target data 
    SUBROUTINE allocate_era_target_fields(tg,nt)
      IMPLICIT NONE

      TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
      INTEGER (KIND=i4), INTENT(in) :: nt !< number of timesteps (12 for monthly mean values)

      INTEGER :: errorcode !< error status variable
        
 
       ALLOCATE (sst_field(1:tg%ie,1:tg%je,1:tg%ke,1:nt), STAT=errorcode)
          IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array sst_field')
      sst_field = 0.0

       ALLOCATE (wsnow_field(1:tg%ie,1:tg%je,1:tg%ke,1:nt), STAT=errorcode)
          IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array wsnow_field')
      wsnow_field = 0.0

       ALLOCATE (t2m_field(1:tg%ie,1:tg%je,1:tg%ke,1:nt), STAT=errorcode)
          IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array t2m_field')
      t2m_field = 0.0

       ALLOCATE (hsurf_field(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
          IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array hsurf_field')
      hsurf_field = 0.0


    END SUBROUTINE allocate_era_target_fields


  

  


END Module mo_era_tg_fields

