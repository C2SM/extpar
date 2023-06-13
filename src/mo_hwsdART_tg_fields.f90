MODULE mo_hwsdART_tg_fields
  USE mo_logging
  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp, &
  & i4


  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_grid_structures, ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fr_heavy_clay, &
   & fr_silty_clay, &
   & fr_light_clay, &
   & fr_silty_clay_loam, &
   & fr_clay_loam, &
   & fr_silt, &
   & fr_silt_loam, &
   & fr_sandy_clay, &
   & fr_loam, &
   & fr_sandy_clay_loam, &
   & fr_sandy_loam, &
   & fr_loamy_sand, &
   & fr_sand, &
   & fr_undef, &
   & allocate_hwsdART_target_fields


  
REAL(KIND=wp), ALLOCATABLE :: fr_heavy_clay(:,:,:)  , &     !< fraction of heavy clay
   &  fr_silty_clay(:,:,:)       , &   !< fraction of silty clay
   &  fr_light_clay(:,:,:)       , &   !< fraction of light clay
   &  fr_silty_clay_loam(:,:,:)  , &   !< fraction of silty clay loam
   &  fr_clay_loam(:,:,:)        , &   !< fraction of clay loam
   &  fr_silt(:,:,:)             , &   !< fraction of silt
   &  fr_silt_loam(:,:,:)        , &   !< fraction of silt loam
   &  fr_sandy_clay(:,:,:)       , &   !< fraction of sandy clay
   &  fr_loam(:,:,:)             , &   !< fraction of loam
   &  fr_sandy_clay_loam(:,:,:)  , &   !< fraction of sandy clay loam
   &  fr_sandy_loam(:,:,:)       , &   !< fraction of sandy loam
   &  fr_loamy_sand(:,:,:)       , &   !< fraction of loamy sand
   &  fr_sand(:,:,:)             , &   !< fraction of sand
   &  fr_undef(:,:,:)                  !< fraction of undef
  

  CONTAINS


  !> allocate fields for GLOBE target data 
  SUBROUTINE allocate_hwsdART_target_fields(tg)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER :: errorcode !< error status variable

    ALLOCATE (fr_heavy_clay(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_heavy_clay')
    fr_heavy_clay = 0.0_wp
    
    ALLOCATE (fr_silty_clay(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_silty_clay')
    fr_silty_clay = 0.0_wp
    
    ALLOCATE (fr_light_clay(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_light_clay')
    fr_light_clay = 0.0_wp
    
    ALLOCATE (fr_silty_clay_loam(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_silty_clay_loam')
    fr_silty_clay_loam = 0.0_wp
    
    ALLOCATE (fr_clay_loam(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_clay_loam')
    fr_clay_loam = 0.0_wp
    
    ALLOCATE (fr_silt(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_silt')
    fr_silt = 0.0_wp
    
    ALLOCATE (fr_silt_loam(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_silt_loam')
    fr_silt_loam = 0.0_wp
    
    ALLOCATE (fr_sandy_clay(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_sandy_clay')
    fr_sandy_clay = 0.0_wp
    
    ALLOCATE (fr_loam(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_loam')
    fr_loam = 0.0_wp
    
    ALLOCATE (fr_sandy_clay_loam(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_sandy_clay_loam')
    fr_sandy_clay_loam = 0.0_wp
    
    ALLOCATE (fr_sandy_loam(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_sandy_loam')
    fr_sandy_loam = 0.0_wp
    
    ALLOCATE (fr_loamy_sand(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_loamy_sand')
    fr_loamy_sand = 0.0_wp
    
    ALLOCATE (fr_sand(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_sand')
    fr_sand = 0.0_wp
        
    ALLOCATE (fr_undef(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_undef')
    fr_undef = 0.0_wp

  END SUBROUTINE allocate_hwsdART_target_fields

END MODULE mo_hwsdART_tg_fields

