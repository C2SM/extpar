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

  USE mo_grid_structures,       ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fr_land_soil
  PUBLIC :: soiltype_fao
  PUBLIC :: soiltype_hwsd
  PUBLIC :: soiltype_hwsd_s
  PUBLIC :: soiltype_deep
  PUBLIC :: fr_sand,fr_silt,fr_clay,fr_oc,fr_bd,fr_dm
  PUBLIC :: fr_sand_deep,fr_silt_deep,fr_clay_deep, &
            fr_oc_deep,fr_bd_deep,fr_dm_deep

  PUBLIC :: allocate_soil_target_fields


  INTEGER(KIND=i4), ALLOCATABLE :: soiltype_fao(:,:,:), & !< soiltype due to FAO Digital Soil map of the World
       &                           soiltype_hwsd(:,:,:), & !< soiltype due to HWSD
       &                           soiltype_hwsd_s(:,:,:), & !< soiltype due to HWSD SUBSOIL 
       &                           soiltype_deep(:,:,:) !< deep soiltype due to HWSD data

  REAL(KIND=wp), ALLOCATABLE    :: fr_land_soil(:,:,:), & !< fraction land due to FAO Digital Soil map of the World
       &                           fr_sand(:,:,:), & !< fraction sand due to HWSD
       &                           fr_silt(:,:,:), & !< fraction silt due to HWSD
       &                           fr_clay(:,:,:), & !< fraction clay due to HWSD
       &                           fr_oc(:,:,:), & !< fraction oc due to HWSD
       &                           fr_bd(:,:,:), & !< fraction bd due to HWSD
       &                           fr_dm(:,:,:), & !< dummy of HWSD
       &                           fr_sand_deep(:,:,:), & !< fraction sand due to HWSD
       &                           fr_silt_deep(:,:,:), & !< fraction silt due to HWSD
       &                           fr_clay_deep(:,:,:), & !< fraction clay due to HWSD
       &                           fr_oc_deep(:,:,:), & !< fraction oc due to HWSD
       &                           fr_bd_deep(:,:,:), & !< fraction bd due to HWSD
       &                           fr_dm_deep(:,:,:) !< dummy of HWSD

  INTEGER(KIND=i4)              :: size_ie, size_je, size_ke

  CONTAINS

  !> allocate fields for GLOBE target data 
  SUBROUTINE allocate_soil_target_fields(tg, ldeep_soil)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    LOGICAL,               INTENT(IN) :: ldeep_soil !< logical switch for deep soil data
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_soil_target_fields')

    ALLOCATE (fr_land_soil(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_soil',__FILE__,__LINE__)
    fr_land_soil = 0.0
        
    ALLOCATE (soiltype_fao(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array soiltype_fao',__FILE__,__LINE__)
    soiltype_fao = 3  ! default value for soiltype is 'sand' (3)

    ALLOCATE (soiltype_hwsd(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array soiltype_hwsd',__FILE__,__LINE__)
    
    soiltype_hwsd = 0

    ALLOCATE (soiltype_hwsd_s(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array soiltype_hwsd',__FILE__,__LINE__)

    soiltype_hwsd_s = 0

    ALLOCATE (fr_sand(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_sand',__FILE__,__LINE__)
    fr_sand = -1.0

    ALLOCATE (fr_silt(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_silt',__FILE__,__LINE__)
    fr_silt = -1.0

    ALLOCATE (fr_clay(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_clay',__FILE__,__LINE__)
    fr_clay = -1.0

    ALLOCATE (fr_oc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_oc',__FILE__,__LINE__)
    fr_oc = -1.0

    ALLOCATE (fr_bd(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_bd',__FILE__,__LINE__)
    fr_bd = -1.0

    ALLOCATE (fr_dm(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_dm',__FILE__,__LINE__)
    fr_dm = -1.0

    ! The following arrays are only conditionally used
    IF (ldeep_soil) THEN
      size_ie = tg%ie
      size_je = tg%je
      size_ke = tg%ke
    ELSE
      size_ie = 0
      size_je = 0
      size_ke = 0
    ENDIF

    ALLOCATE (soiltype_deep(1:size_ie,1:size_je,1:size_ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array soiltype_deep',__FILE__,__LINE__)
    soiltype_deep = 3  ! default value for soiltype is 'sand' (3)

    ALLOCATE (fr_sand_deep(1:size_ie,1:size_je,1:size_ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_sand_deep',__FILE__,__LINE__)
    fr_sand_deep = -1.0

    ALLOCATE (fr_silt_deep(1:size_ie,1:size_je,1:size_ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_silt_deep',__FILE__,__LINE__)
    fr_silt_deep = -1.0

    ALLOCATE (fr_clay_deep(1:size_ie,1:size_je,1:size_ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_clay_deep',__FILE__,__LINE__)
    fr_clay_deep = -1.0
      
    ALLOCATE (fr_oc_deep(1:size_ie,1:size_je,1:size_ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_oc_deep',__FILE__,__LINE__)
    fr_oc_deep = -1.0

    ALLOCATE (fr_bd_deep(1:size_ie,1:size_je,1:size_ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_bd_deep',__FILE__,__LINE__)
    fr_bd_deep = -1.0

    ALLOCATE (fr_dm_deep(1:size_ie,1:size_je,1:size_ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_dm_deep',__FILE__,__LINE__)
    fr_dm_deep = -1.0

  END SUBROUTINE allocate_soil_target_fields

END MODULE mo_soil_tg_fields

