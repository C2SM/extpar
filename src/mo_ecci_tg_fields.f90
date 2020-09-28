!+ Fortran module for ecci data specification on target grid for external Parameters

! Description:
! The ecci dataset contains the following land use classification scheme
!

! class no. value         description
!
!  01   11  'Post-flooding or irrigated croplands  ' 
!  02   14  'Rainfed croplands' 
!  03   20  'Mosaic Cropland (50-70%) / Vegetation (grassland, shrubland, forest) (20-50%)' 
!  04   30  'Mosaic Vegetation (grassland, shrubland, forest) (50-70%) / Cropland (20-50%) '
!  05   40  'Closed to open (>15%) broadleaved evergreen and/or semi-deciduous forest (>5m)   ' 
!  06   50  'Closed (>40%) broadleaved deciduous forest (>5m)           ' 
!  07   60  'Open (15-40%) broadleaved deciduous forest (>5m)  ' 
!  08   70  'Closed (>40%) needleleaved evergreen forest (>5m) ' 
!  09   90  'Open (15-40%) needleleaved deciduous or evergreen forest (>5m) 
!  10   100 'Closed (>40%) needleleaved evergreen forest (>5m) '       '
!  11   110 'Mosaic Forest/Shrubland (50-70%) / Grassland (20-50%)' 
!  12   120 'Mosaic Grassland (50-70%) / Forest/Shrubland (20-50%) ' 
!  13   130 'Closed to open (>15%) shrubland (<5m)' 
!  14   140 'Closed to open (>15%) grassland' 
!  15   150 'Sparse (>15%) vegetation (woody vegetation, shrubs, grassland)  ' 
!  16   160 'Closed (>40%) broadleaved forest regularly flooded - Fresh water  ' 
!  17   170 'Closed (>40%) broadleaved semi-deciduous and/or evergreen forest regularly flooded - Saline water' 
!  18   180 'Closed to open (>15%) vegetation (grassland, shrubland, woody vegetation) on regularly flooded or &
! &          waterlogged soil - Fresh, brackish or saline water '
!  19   190 ' regularly flooded or waterlogged soil - Fresh, brackish or saline water    '
!  20   200 'bare areas                      '
!  21   210 'water bodies              ' 
!  22   220 'artificial surfaces          ' 
!  23   230 'undefined                    '
!
!

MODULE mo_ecci_tg_fields

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  USE mo_ecci_lookup_tables,    ONLY: nclass_ecci

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fr_land_ecci, &
       &    ecci_class_fraction,    &
       &    ecci_class_npixel, &
       &    ecci_tot_npixel, &
       &    ice_ecci, &
       &    z0_ecci, &
       &    root_ecci, &
       &    plcov_mn_ecci, &
       &    plcov_mx_ecci, &
       &    lai_mn_ecci, &
       &    lai_mx_ecci, &
       &    rs_min_ecci, &
       &    urban_ecci,  &
       &    for_d_ecci,  &
       &    for_e_ecci, &
       &    skinc_ecci, &
       &    emissivity_ecci, &
       &    allocate_ecci_target_fields

  INTEGER (KIND=i4), ALLOCATABLE :: ecci_class_npixel(:,:,:,:), & !< number of raw data pixels for each &
       &                            ecci_tot_npixel(:,:,:)  !< total number of ecci raw data pixels &

  REAL (KIND=wp), ALLOCATABLE  :: fr_land_ecci(:,:,:), & !< fraction land due to ecci raw data
      &                          ice_ecci(:,:,:), &     !< fraction of ice due to ecci raw data
      &                          z0_ecci(:,:,:), &      !< roughness length due to ecci land use data
      &                          root_ecci(:,:,:), &    !< root depth due to ecci land use data
      &                          plcov_mx_ecci(:,:,:), &!< plant cover maximum due to ecci land use data
      &                          plcov_mn_ecci(:,:,:), &!< plant cover minimum due to ecci land use data
      &                          lai_mx_ecci(:,:,:), &  !< Leaf Area Index maximum due to ecci land use data
      &                          lai_mn_ecci(:,:,:), &  !< Leaf Area Index minimum due to ecci land use data
      &                          rs_min_ecci(:,:,:), &  !< minimal stomata resistance due to ecci land use data
      &                          urban_ecci(:,:,:), &   !< urban fraction due to ecci land use data
      &                          for_d_ecci(:,:,:), &   !< deciduous forest (fraction) due to ecci land use data
      &                          for_e_ecci(:,:,:), &   !< evergreen forest (fraction) due to ecci land use data
      &                          skinc_ecci(:,:,:), &   !< skin conductivity due to ecci land use data
      &                          emissivity_ecci(:,:,:), & !< longwave emissivity due to ecci land use data
      &                          ecci_class_fraction(:,:,:,:)  !< fraction for each ecci class &

  CONTAINS

  SUBROUTINE allocate_ecci_target_fields(tg)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    ALLOCATE (fr_land_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_ecci',__FILE__,__LINE__)
    fr_land_ecci = 0.0

    ALLOCATE (ecci_tot_npixel(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ecci_tot_npixel',__FILE__,__LINE__)
    ecci_tot_npixel = 0

    ALLOCATE (ecci_class_fraction(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_ecci), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ecci_class_fraction',__FILE__,__LINE__)
    ecci_class_fraction = 0.0

    ALLOCATE (ecci_class_npixel(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_ecci), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ecci_class_npixel',__FILE__,__LINE__)
    ecci_class_npixel = 0

    ALLOCATE (ice_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ice_ecci',__FILE__,__LINE__)
    ice_ecci = 0.0

    ALLOCATE (z0_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0_ecci',__FILE__,__LINE__)
    z0_ecci = 0.0

    ALLOCATE (root_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array root_ecci',__FILE__,__LINE__)
    root_ecci = 0.0

    ALLOCATE (plcov_mx_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mx_ecci',__FILE__,__LINE__)
    plcov_mx_ecci = 0.0

    ALLOCATE (plcov_mn_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mn_ecci',__FILE__,__LINE__)
    plcov_mn_ecci = 0.0

    ALLOCATE (lai_mx_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mx_ecci',__FILE__,__LINE__)
    lai_mx_ecci = 0.0

    ALLOCATE (lai_mn_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mn_ecci',__FILE__,__LINE__)
    lai_mn_ecci = 0.0

    ALLOCATE (rs_min_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array rs_min_ecci',__FILE__,__LINE__)
    rs_min_ecci = 0.0

    ALLOCATE (urban_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array urban_ecci',__FILE__,__LINE__)
    urban_ecci = 0.0

    ALLOCATE (for_d_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_d_ecci',__FILE__,__LINE__)
    for_d_ecci = 0.0

    ALLOCATE (for_e_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_e_ecci',__FILE__,__LINE__)
    for_e_ecci = 0.0

    ALLOCATE (skinc_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array skinc_ecci',__FILE__,__LINE__)
    skinc_ecci = 0.0

    ALLOCATE (emissivity_ecci(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emissivity_ecci',__FILE__,__LINE__)
    emissivity_ecci = 0.0

  END SUBROUTINE allocate_ecci_target_fields

END MODULE mo_ecci_tg_fields
