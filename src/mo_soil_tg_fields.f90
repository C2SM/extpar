!+ Fortran module for soil data on target grid for external parameters
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V2_0         1013/06/04 Martina Messmer
!  adaptations in a way that HWSD data set can be used (top- and subsoil)
!  Code received from Juergen Helmert
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for soil data on target grid for external parameters
!> \author Hermann Asensio
MODULE mo_soil_tg_fields

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4
  USE mo_array_cache,           ONLY: allocate_cached
  USE mo_grid_structures,       ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fr_land_soil
  PUBLIC :: soiltype_fao
  PUBLIC :: soiltype_hwsd
  PUBLIC :: fr_sand,fr_silt,fr_clay,fr_oc,fr_bd,fr_dm

  PUBLIC :: allocate_soil_target_fields


  INTEGER(KIND=i4), POINTER ::     soiltype_fao(:,:,:), & !< soiltype due to FAO Digital Soil map of the World
       &                           soiltype_hwsd(:,:,:) !< soiltype due to HWSD

  REAL(KIND=wp), POINTER        :: fr_land_soil(:,:,:), & !< fraction land due to FAO Digital Soil map of the World
       &                           fr_sand(:,:,:), & !< fraction sand due to HWSD
       &                           fr_silt(:,:,:), & !< fraction silt due to HWSD
       &                           fr_clay(:,:,:), & !< fraction clay due to HWSD
       &                           fr_oc(:,:,:), & !< fraction oc due to HWSD
       &                           fr_bd(:,:,:), & !< fraction bd due to HWSD
       &                           fr_dm(:,:,:) !< dummy of HWSD

  CONTAINS

  !> allocate fields for GLOBE target data
  SUBROUTINE allocate_soil_target_fields(tg, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    LOGICAL, INTENT(in)               :: l_use_array_cache 

    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0
    
    CALL logging%info('Enter routine: allocate_soil_target_fields')

if (l_use_array_cache) then
   call allocate_cached('fr_land_soil', fr_land_soil, [tg%ie,tg%je,tg%ke])
else
   allocate(fr_land_soil(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_soil',__FILE__,__LINE__)
    fr_land_soil = 0.0

if (l_use_array_cache) then
   call allocate_cached('soiltype_fao', soiltype_fao, [tg%ie,tg%je,tg%ke])
else
   allocate(soiltype_fao(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array soiltype_fao',__FILE__,__LINE__)
    soiltype_fao = 3  ! default value for soiltype is 'sand' (3)

if (l_use_array_cache) then
   call allocate_cached('soiltype_hwsd', soiltype_hwsd, [tg%ie,tg%je,tg%ke])
else
   allocate(soiltype_hwsd(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array soiltype_hwsd',__FILE__,__LINE__)

    soiltype_hwsd = 0

if (l_use_array_cache) then
   call allocate_cached('fr_sand', fr_sand, [tg%ie,tg%je,tg%ke])
else
   allocate(fr_sand(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_sand',__FILE__,__LINE__)
    fr_sand = -1.0

if (l_use_array_cache) then
   call allocate_cached('fr_silt', fr_silt, [tg%ie,tg%je,tg%ke])
else
   allocate(fr_silt(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_silt',__FILE__,__LINE__)
    fr_silt = -1.0

if (l_use_array_cache) then
   call allocate_cached('fr_clay', fr_clay, [tg%ie,tg%je,tg%ke])
else
   allocate(fr_clay(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_clay',__FILE__,__LINE__)
    fr_clay = -1.0

if (l_use_array_cache) then
   call allocate_cached('fr_oc', fr_oc, [tg%ie,tg%je,tg%ke])
else
   allocate(fr_oc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_oc',__FILE__,__LINE__)
    fr_oc = -1.0

if (l_use_array_cache) then
   call allocate_cached('fr_bd', fr_bd, [tg%ie,tg%je,tg%ke])
else
   allocate(fr_bd(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_bd',__FILE__,__LINE__)
    fr_bd = -1.0

if (l_use_array_cache) then
   call allocate_cached('fr_dm', fr_dm, [tg%ie,tg%je,tg%ke])
else
   allocate(fr_dm(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_dm',__FILE__,__LINE__)
    fr_dm = -1.0

  END SUBROUTINE allocate_soil_target_fields

END MODULE mo_soil_tg_fields

