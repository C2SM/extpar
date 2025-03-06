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

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4
  USE mo_array_cache,           ONLY: allocate_cached
  USE mo_grid_structures,       ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: allocate_aot_target_fields
  PUBLIC :: aot_tg

  !< aerosol optical thickness, aot_tg(ie,je,ke,,ntype,ntime)
  REAL(KIND=wp), POINTER :: aot_tg(:,:,:,:,:)

  CONTAINS

  !> allocate fields for TARGET grid
  !!
  !! the target grid for the GME has 3 dimension (ie,je,jd),
  !! the target grid for the COSMO model has 2 dimension (ie,je)
  !! the target grid for the ICON model has 1 dimension (ne)
  !! depending of the target model the second and third dimension of the target fields should be
  !! allocated with the length 1
  SUBROUTINE allocate_aot_target_fields(tg, iaot_type, ntime, ntype,l_use_array_cache)

   TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
   INTEGER (KIND=i4), INTENT(IN)     :: iaot_type, & !< type of data source
         &                               ntime, & !< number of times
         &                               ntype !< number of types of aerosol

   LOGICAL, INTENT(in)               :: l_use_array_cache

   INTEGER(KIND=i4)                  :: errorcode !< error status variable

   errorcode = 0
    
   CALL logging%info('Enter routine: allocate_aot_target_fields')

   IF (l_use_array_cache) then
      CALL allocate_cached('aot_tg', aot_tg, [tg%ie,tg%je,tg%ke,ntype,ntime])
   ELSE
      ALLOCATE(aot_tg(tg%ie,tg%je,tg%ke,ntype,ntime), stat=errorcode)
   ENDIF
   IF(errorcode /= 0) CALL logging%error('Cant allocate the array aot_tg',__FILE__,__LINE__)
   aot_tg = 0.0_wp

   CALL logging%info('Exit routine: allocate_aot_target_fields')

END SUBROUTINE allocate_aot_target_fields

END MODULE mo_aot_target_fields
