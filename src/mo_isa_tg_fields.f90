!+ Fortran module for land use data on target grid for external parameters 
!
!
! Description:
! Fortran module for land use data on target grid for external parameters
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
!> Fortran module for land use data on target grid for external parameters 
!> \author Hermann Asensio
!
MODULE mo_isa_tg_fields

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: isa_field, &
       &    isa_tot_npixel


  PUBLIC :: allocate_isa_target_fields, allocate_add_isa_fields

  REAL (KIND=wp), ALLOCATABLE    :: isa_field(:,:,:) !< fraction land due to land use raw data

  INTEGER (KIND=i4), ALLOCATABLE :: isa_tot_npixel(:,:,:)  

  CONTAINS

  !> allocate fields for TARGET grid
  !!
  !! the target grid for the GME has 3 dimension (ie,je,jd),
  !! the target grid for the COSMO model has 2 dimension (ie,je)
  !! the target grid for the ICON model has 1 dimension (ne)
  !! depending of the target model the second and third dimension of the target fields should be 
  !! allocated with the length 1
  SUBROUTINE allocate_isa_target_fields(tg)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER(KIND=i4)                   :: errorcode !< error status variable
   
    CALL logging%info('Enter routine: allocate_isa_target_fields')

    ALLOCATE (isa_field(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array isa_field',__FILE__,__LINE__)
    isa_field = 0.0

  END SUBROUTINE allocate_isa_target_fields

  !> allocate additional land use target fields
  SUBROUTINE allocate_add_isa_fields(tg)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

    INTEGER(KIND=i4)                   :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_add_isa_fields')

    ALLOCATE (isa_tot_npixel(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array isa_tot_npixel',__FILE__,__LINE__)
    isa_tot_npixel = 0

  END SUBROUTINE allocate_add_isa_fields

END MODULE mo_isa_tg_fields
