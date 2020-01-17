!+ Fortran module for globcover data specification on target grid for external Parameters
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_3         2011/04/19 Hermann Asensio
!  Initial release
!
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for globcover data specification on target grid for external Parameters
!> \author Hermann Asensio
!
! Description:
! The globcover dataset contains the following land use classification scheme
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

MODULE mo_globcover_tg_fields

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4


!> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_grid_structures, ONLY: target_grid_def

  USE mo_globcover_lookup_tables, ONLY: nclass_globcover

IMPLICIT NONE

PRIVATE

PUBLIC :: fr_land_globcover, &
          globcover_class_fraction,    &
          globcover_class_npixel, &
          globcover_tot_npixel, &
          ice_globcover, &
          z0_globcover, &
          root_globcover, &
          plcov_mn_globcover, &
          plcov_mx_globcover, &
          lai_mn_globcover, &
          lai_mx_globcover, &
          rs_min_globcover, &
          urban_globcover,  &
          for_d_globcover,  &
          for_e_globcover, &
          skinc_globcover, &
          emissivity_globcover, &
          allocate_globcover_target_fields


       REAL (KIND=wp), ALLOCATABLE  :: globcover_class_fraction(:,:,:,:)  !< fraction for each globcover class &
! &                                    on target grid (dimension (ie,je,ke,nclass_globcover))

       INTEGER (KIND=i4), ALLOCATABLE :: globcover_class_npixel(:,:,:,:) !< number of raw data pixels for each &
! &                                      globcover class on target grid (dimension (ie,je,ke,nclass_globcover))


       INTEGER (KIND=i4), ALLOCATABLE :: globcover_tot_npixel(:,:,:)  !< total number of globcover raw data pixels &
! &                                      on target grid (dimension (ie,je,ke))


       REAL (KIND=wp), ALLOCATABLE  :: fr_land_globcover(:,:,:) !< fraction land due to globcover raw data
       REAL (KIND=wp), ALLOCATABLE  :: ice_globcover(:,:,:)     !< fraction of ice due to globcover raw data
       REAL (KIND=wp), ALLOCATABLE  :: z0_globcover(:,:,:)      !< roughness length due to globcover land use data
       REAL (KIND=wp), ALLOCATABLE  :: root_globcover(:,:,:)    !< root depth due to globcover land use data
       REAL (KIND=wp), ALLOCATABLE  :: plcov_mx_globcover(:,:,:)!< plant cover maximum due to globcover land use data
       REAL (KIND=wp), ALLOCATABLE  :: plcov_mn_globcover(:,:,:)!< plant cover minimum due to globcover land use data
       REAL (KIND=wp), ALLOCATABLE  :: lai_mx_globcover(:,:,:)  !< Leaf Area Index maximum due to globcover land use data
       REAL (KIND=wp), ALLOCATABLE  :: lai_mn_globcover(:,:,:)  !< Leaf Area Index minimum due to globcover land use data
       REAL (KIND=wp), ALLOCATABLE  :: rs_min_globcover(:,:,:)  !< minimal stomata resistance due to globcover land use data
       REAL (KIND=wp), ALLOCATABLE  :: urban_globcover(:,:,:)   !< urban fraction due to globcover land use data
       REAL (KIND=wp), ALLOCATABLE  :: for_d_globcover(:,:,:)   !< deciduous forest (fraction) due to globcover land use data
       REAL (KIND=wp), ALLOCATABLE  :: for_e_globcover(:,:,:)   !< evergreen forest (fraction) due to globcover land use data
       REAL (KIND=wp), ALLOCATABLE  :: skinc_globcover(:,:,:)   !< skin conductivity due to globcover land use data
       REAL (KIND=wp), ALLOCATABLE  :: emissivity_globcover(:,:,:) !< longwave emissivity due to globcover land use data


CONTAINS





!> allocate fields for TARGET grid
!!
!! the target grid for the GME has 3 dimension (ie,je,jd),
!! the target grid for the COSMO model has 2 dimension (ie,je)
!! the target grid for the ICON model has 1 dimension (ne)
!! depending of the target model the second and third dimension of the target fields should be
!! allocated with the length 1
  SUBROUTINE allocate_globcover_target_fields(tg)


    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description


    INTEGER :: errorcode !< error status variable


    ALLOCATE (fr_land_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_land_globcover')
    fr_land_globcover = 0.0

    ALLOCATE (globcover_tot_npixel(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array globcover_tot_npixel')
    globcover_tot_npixel = 0

     ALLOCATE (globcover_class_fraction(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_globcover), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array globcover_class_fraction')
    globcover_class_fraction = 0.0


     ALLOCATE (globcover_class_npixel(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_globcover), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array globcover_class_npixel')
    globcover_class_npixel = 0

     ALLOCATE (ice_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array ice_globcover')
    ice_globcover = 0.0

     ALLOCATE (z0_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array z0_globcover')
    z0_globcover = 0.0

     ALLOCATE (root_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array root_globcover')
    root_globcover = 0.0

     ALLOCATE (plcov_mx_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array plcov_mx_globcover')
    plcov_mx_globcover = 0.0


     ALLOCATE (plcov_mn_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array plcov_mn_globcover')
    plcov_mn_globcover = 0.0

     ALLOCATE (lai_mx_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lai_mx_globcover')
    lai_mx_globcover = 0.0

     ALLOCATE (lai_mn_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lai_mn_globcover')
    lai_mn_globcover = 0.0

     ALLOCATE (rs_min_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array rs_min_globcover')
    rs_min_globcover = 0.0

     ALLOCATE (urban_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array urban_globcover')
    urban_globcover = 0.0


    ALLOCATE (for_d_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array for_d_globcover')
    for_d_globcover = 0.0

    ALLOCATE (for_e_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array for_e_globcover')
    for_e_globcover = 0.0

    ALLOCATE (skinc_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array skinc_globcover')
    skinc_globcover = 0.0

    ALLOCATE (emissivity_globcover(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array emissivity_globcover')
    emissivity_globcover = 0.0

  END SUBROUTINE allocate_globcover_target_fields



END Module mo_globcover_tg_fields


