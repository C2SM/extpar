!+ Fortran module for Aerosol optical thickness data, specification of the target grid fields
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_3         2011/04/19 Hermann Asensio
!  cleanup of code
! V4_0         2016/08/17 authors from RHM and Daniel Lthi
!  added support for MACv2 spectrally stratified monthly aerosol fields
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for Aerosol optical thickness data, specification of the target grid fields 
!> \author Hermann Asensio
!

MODULE mo_aot_target_fields

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_grid_structures, ONLY: target_grid_def


IMPLICIT NONE

PRIVATE

PUBLIC :: allocate_aot_target_fields
PUBLIC :: aot_tg
PUBLIC :: MAC_aot_tg, MAC_ssa_tg, MAC_asy_tg !new

REAL (KIND=wp), ALLOCATABLE :: aot_tg(:,:,:,:,:) !< aerosol optical thickness, aot_tg(ie,je,ke,,ntype,ntime) 

REAL (KIND=wp), ALLOCATABLE :: MAC_aot_tg(:,:,:,:), MAC_ssa_tg(:,:,:,:), MAC_asy_tg(:,:,:,:) 



CONTAINS

!> allocate fields for TARGET grid
!!
!! the target grid for the GME has 3 dimension (ie,je,jd),
!! the target grid for the COSMO model has 2 dimension (ie,je)
!! the target grid for the ICON model has 1 dimension (ne)
!! depending of the target model the second and third dimension of the target fields should be 
!! allocated with the length 1
  SUBROUTINE allocate_aot_target_fields(tg, iaot_type, ntime, ntype, n_spectr)
  
    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN) :: iaot_type !< type of data source
    INTEGER (KIND=i8), INTENT(IN) :: ntime !< number of times
    INTEGER (KIND=i8), INTENT(IN) :: ntype !< number of types of aerosol
    INTEGER (KIND=i8), INTENT(IN) :: n_spectr !< number of spectral intervals

    INTEGER :: errorcode !< error status variable
    IF (iaot_type == 4) THEN
      ALLOCATE (MAC_aot_tg(1:tg%ie,1:tg%je,1:n_spectr,1:ntime), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array MAC_aot_tg')
      MAC_aot_tg = 0.0

      ALLOCATE (MAC_ssa_tg(1:tg%ie,1:tg%je,1:n_spectr,1:ntime), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array MAC_ssa_tg')
      MAC_ssa_tg = 0.0
   
      ALLOCATE (MAC_asy_tg(1:tg%ie,1:tg%je,1:n_spectr,1:ntime), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array MAC_asy_tg')
      MAC_asy_tg = 0.0
    ELSE
    ALLOCATE (aot_tg(1:tg%ie,1:tg%je,1:tg%ke,1:ntype,1:ntime), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array aot_tg')
    aot_tg = 0.0
    ENDIF

  END SUBROUTINE allocate_aot_target_fields


END MODULE mo_aot_target_fields

