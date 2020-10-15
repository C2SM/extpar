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

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4
  USE mo_array_cache,           ONLY: allocate_cached
  USE mo_grid_structures,       ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sst_field, &
    &        wsnow_field, &
    &        t2m_field, &
    &        hsurf_field, &
    &        allocate_era_target_fields


  REAL(KIND=wp), POINTER :: sst_field(:,:,:,:), & !< field for sst data (12 months)
       &                    wsnow_field(:,:,:,:), & !< field for wsnow data (12 months)
       &                    t2m_field(:,:,:,:), & !< field for wsnow data (12 months)
       &                    hsurf_field(:,:,:) !< field for wsnow data (12 months)

  CONTAINS

  !> allocate fields for GLOBE target data
  SUBROUTINE allocate_era_target_fields(tg,nt, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(in)     :: nt !< number of timesteps (12 for monthly mean values)
    LOGICAL, INTENT(in)               :: l_use_array_cache

    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0
    
    CALL logging%info('Enter routine: allocate_era_target_fields')

if (l_use_array_cache) then
   call allocate_cached('sst_field', sst_field, [tg%ie,tg%je,tg%ke,nt])
else
   allocate(sst_field(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array sst_field',__FILE__,__LINE__)
    sst_field = 0.0

if (l_use_array_cache) then
   call allocate_cached('wsnow_field', wsnow_field, [tg%ie,tg%je,tg%ke,nt])
else
   allocate(wsnow_field(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array wsnow_field',__FILE__,__LINE__)
    wsnow_field = 0.0

if (l_use_array_cache) then
   call allocate_cached('t2m_field', t2m_field, [tg%ie,tg%je,tg%ke,nt])
else
   allocate(t2m_field(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array t2m_field',__FILE__,__LINE__)
    t2m_field = 0.0

if (l_use_array_cache) then
   call allocate_cached('hsurf_field', hsurf_field, [tg%ie,tg%je,tg%ke])
else
   allocate(hsurf_field(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hsurf_field',__FILE__,__LINE__)
    hsurf_field = 0.0

    CALL logging%info('Exit routine: allocate_era_target_fields')

  END SUBROUTINE allocate_era_target_fields

END MODULE mo_era_tg_fields

