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
  PUBLIC :: MAC_aot_tg, MAC_ssa_tg, MAC_asy_tg !new
  PUBLIC :: CAMS_tg

  !< aerosol optical thickness, aot_tg(ie,je,ke,,ntype,ntime)
  REAL(KIND=wp), POINTER :: aot_tg(:,:,:,:,:), &
       &                    MAC_aot_tg(:,:,:,:),&
       &                    MAC_ssa_tg(:,:,:,:),&
       &                    MAC_asy_tg(:,:,:,:),&
       &                    CAMS_tg(:,:,:,:,:)  
	   
  CONTAINS

  !> allocate fields for TARGET grid
  !!
  !! the target grid for the GME has 3 dimension (ie,je,jd),
  !! the target grid for the COSMO model has 2 dimension (ie,je)
  !! the target grid for the ICON model has 1 dimension (ne)
  !! depending of the target model the second and third dimension of the target fields should be
  !! allocated with the length 1
  SUBROUTINE allocate_aot_target_fields(tg, iaot_type, ntime, ntype, n_spectr, nlevel_cams, ntype_cams, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN)     :: iaot_type, & !< type of data source
         &                               ntime, & !< number of times
         &                               ntype, & !< number of types of aerosol
         &                               n_spectr, & !< number of spectral intervals
         &                               ntype_cams, & !< number of types of aerosols in CAMS
         &                               nlevel_cams   !< number of level in CAMS
		 
    LOGICAL, INTENT(in)               :: l_use_array_cache

    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_aot_target_fields')

    IF (iaot_type == 4) THEN
if (l_use_array_cache) then
   call allocate_cached('MAC_aot_tg', MAC_aot_tg, [tg%ie,tg%je,n_spectr,ntime])
else
   allocate(MAC_aot_tg(tg%ie,tg%je,n_spectr,ntime), stat=errorcode)
endif
      IF(errorcode /= 0) CALL logging%error('Cant allocate the array MAC_aot_tg',__FILE__,__LINE__)
      MAC_aot_tg = 0.0_wp

if (l_use_array_cache) then
   call allocate_cached('MAC_ssa_tg', MAC_ssa_tg, [tg%ie,tg%je,n_spectr,ntime])
else
   allocate(MAC_ssa_tg(tg%ie,tg%je,n_spectr,ntime), stat=errorcode)
endif
      IF(errorcode /= 0) CALL logging%error('Cant allocate the array MAC_ssa_tg',__FILE__,__LINE__)
      MAC_ssa_tg = 0.0_wp

if (l_use_array_cache) then
   call allocate_cached('MAC_asy_tg', MAC_asy_tg, [tg%ie,tg%je,n_spectr,ntime])
else
   allocate(MAC_asy_tg(tg%ie,tg%je,n_spectr,ntime), stat=errorcode)
endif
      IF(errorcode /= 0) CALL logging%error('Cant allocate the array MAC_asy_tg',__FILE__,__LINE__)
      MAC_asy_tg = 0.0_wp

if (l_use_array_cache) then
   call allocate_cached('aot_tg', aot_tg, [0,0,0,0,0])
else
   allocate(aot_tg(0,0,0,0,0), stat=errorcode)
endif
    
    ELSEIF (iaot_type == 5) THEN
if (l_use_array_cache) then
   call allocate_cached('CAMS_tg', CAMS_tg, [tg%ie,tg%je,nlevel_cams,ntime,ntype_cams])
else
   allocate(CAMS_tg(tg%ie,tg%je,nlevel_cams,ntime,ntype_cams), stat=errorcode)
endif
      IF(errorcode /= 0) CALL logging%error('Cant allocate the array CAMS_tg',__FILE__,__LINE__)
      CAMS_tg = 0.0_wp
if (l_use_array_cache) then
   call allocate_cached('aot_tg', aot_tg, [0,0,0,0,0])
else
   allocate(aot_tg(0,0,0,0,0), stat=errorcode)
endif

    ELSE
if (l_use_array_cache) then
   call allocate_cached('aot_tg', aot_tg, [tg%ie,tg%je,tg%ke,ntype,ntime])
else
   allocate(aot_tg(tg%ie,tg%je,tg%ke,ntype,ntime), stat=errorcode)
endif
      IF(errorcode /= 0) CALL logging%error('Cant allocate the array aot_tg',__FILE__,__LINE__)
      aot_tg = 0.0_wp

if (l_use_array_cache) then
   call allocate_cached('MAC_aot_tg', MAC_aot_tg, [0,0,0,0])
else
   allocate(MAC_aot_tg(0,0,0,0), stat=errorcode)
endif
if (l_use_array_cache) then
   call allocate_cached('MAC_ssa_tg', MAC_ssa_tg, [0,0,0,0])
else
   allocate(MAC_ssa_tg(0,0,0,0), stat=errorcode)
endif
if (l_use_array_cache) then
   call allocate_cached('MAC_asy_tg', MAC_asy_tg, [0,0,0,0])
else
   allocate(MAC_asy_tg(0,0,0,0), stat=errorcode)
endif
if (l_use_array_cache) then
   call allocate_cached('CAMS_tg', CAMS_tg, [0,0,0,0,0])
else
   allocate(CAMS_tg(0,0,0,0,0), stat=errorcode)
endif
    ENDIF

    CALL logging%info('Exit routine: allocate_aot_target_fields')

  END SUBROUTINE allocate_aot_target_fields

END MODULE mo_aot_target_fields
