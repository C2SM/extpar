!+ Fortran module for GLC2000 data specification on target grid for external Parameters
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for GLC2000 data specification on target grid for external Parameters 
!> \author Hermann Asensio
!
! Description:
! The GLC2000 dataset contains the following land use classification scheme
!

! class no.          description               !        IGBP-correspondence
!
!  01   'evergreen broadleaf tree     ' , &    !        evergreen broadleaf forest
!  02   'deciduous broadleaf tree clos' , &    !        deciduous broadleaf forest
!  03   'deciduous broadleaf tree open' , &    !        woody savannas & savanna
!  04   'evergreen needleleaf tree    ' , &    !        evergreen needleleaf forest
!  05   'deciduous needleleaf tree    ' , &    !        deciduous needleleaf forest
!  06   'mixed leaf tree              ' , &    !        mixed forest
!  07   'fresh water flooded tree     ' , &    !        evergreen broadleaf
!  08   'saline water flooded tree    ' , &    !        evergreen broadleaf
!  09   'mosaic tree / other nat.veg. ' , &    !        shrubland ?
!  10   'burnt tree cover             ' , &    ! 
!  11   'evergreen shrubs closed-open ' , &    !        shrubland open-closed
!  12   'deciduous shrubs closed-open ' , &    !        savannas
!  13   'herbaceous cover closed-open ' , &    !        grassland
!  14   'sparse herbaceous or grass   ' , &    !
!  15   'flooded shrub or herbaceous  ' , &    !        wetlands
!  16   'cultivated & managed areas   ' , &    !        croplands
!  17   'mosaic crop/tree/natural veg.' , &    !        cropland/other vegetation mosaic
!  18   'mosaic crop/shrub or grass   ' , &    !        cropland/other vegetation mosaic
!  19   'bare areas                   ' , &    !        barren or sparsely vegetated
!  20   'water bodies                 ' , &    !        water
!  21   'snow & ice                   ' , &    !        snow & ice
!  22   'artificial surfaces          ' , &    !        urban and built-up areas
!  23   'undefined                    ' /      ! 
!
!

MODULE mo_glc2000_tg_fields

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8


!> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_grid_structures, ONLY: target_grid_def

  USE mo_glc2000_lookup_tables, ONLY: nclass_glc2000

IMPLICIT NONE

PRIVATE

PUBLIC :: fr_land_glc2000, &
          glc2000_class_fraction,    &
          glc2000_class_npixel, &
          glc2000_tot_npixel, &
          ice_glc2000, &
          z0_glc2000, &
          root_glc2000, &
          plcov_mn_glc2000, &
          plcov_mx_glc2000, &
          lai_mn_glc2000, &
          lai_mx_glc2000, &
          rs_min_glc2000, &
          urban_glc2000,  &
          for_d_glc2000,  &
          for_e_glc2000, &
          emissivity_glc2000, &
          allocate_glc2000_target_fields


       REAL (KIND=wp), ALLOCATABLE  :: glc2000_class_fraction(:,:,:,:)  !< fraction for each glc2000 class &
! & on target grid (dimension (ie,je,ke,nclass_glc2000))

       INTEGER (KIND=i8), ALLOCATABLE :: glc2000_class_npixel(:,:,:,:) !< number of raw data pixels for each &
! & glc2000 class on target grid (dimension (ie,je,ke,nclass_glc2000))


       INTEGER (KIND=i8), ALLOCATABLE :: glc2000_tot_npixel(:,:,:)  !< total number of glc2000 raw data pixels &
! & on target grid (dimension (ie,je,ke))


       REAL (KIND=wp), ALLOCATABLE  :: fr_land_glc2000(:,:,:) !< fraction land due to glc2000 raw data
       REAL (KIND=wp), ALLOCATABLE  :: ice_glc2000(:,:,:)     !< fraction of ice due to glc2000 raw data
       REAL (KIND=wp), ALLOCATABLE  :: z0_glc2000(:,:,:)      !< roughness length due to glc2000 land use data
       REAL (KIND=wp), ALLOCATABLE  :: root_glc2000(:,:,:)    !< root depth due to glc2000 land use data
       REAL (KIND=wp), ALLOCATABLE  :: plcov_mx_glc2000(:,:,:)!< plant cover maximum due to glc2000 land use data
       REAL (KIND=wp), ALLOCATABLE  :: plcov_mn_glc2000(:,:,:)!< plant cover minimum due to glc2000 land use data
       REAL (KIND=wp), ALLOCATABLE  :: lai_mx_glc2000(:,:,:)  !< Leaf Area Index maximum due to glc2000 land use data
       REAL (KIND=wp), ALLOCATABLE  :: lai_mn_glc2000(:,:,:)  !< Leaf Area Index minimum due to glc2000 land use data
       REAL (KIND=wp), ALLOCATABLE  :: rs_min_glc2000(:,:,:)  !< minimal stomata resistance due to glc2000 land use data
       REAL (KIND=wp), ALLOCATABLE  :: urban_glc2000(:,:,:)   !< urban fraction due to glc2000 land use data
       REAL (KIND=wp), ALLOCATABLE  :: for_d_glc2000(:,:,:)   !< deciduous forest (fraction) due to glc2000 land use data
       REAL (KIND=wp), ALLOCATABLE  :: for_e_glc2000(:,:,:)   !< evergreen forest (fraction) due to glc2000 land use data
       REAL (KIND=wp), ALLOCATABLE  :: emissivity_glc2000(:,:,:) !< longwave emissivity due to glc2000 land use data


CONTAINS





!> allocate fields for TARGET grid
!!
!! the target grid for the GME has 3 dimension (ie,je,jd),
!! the target grid for the COSMO model has 2 dimension (ie,je)
!! the target grid for the ICON model has 1 dimension (ne)
!! depending of the target model the second and third dimension of the target fields should be 
!! allocated with the length 1
  subroutine allocate_glc2000_target_fields(tg)
  

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description


    INTEGER :: errorcode !< error status variable

   
    ALLOCATE (fr_land_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_land_glc2000')
    fr_land_glc2000 = 0.0

    allocate (glc2000_tot_npixel(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array glc2000_tot_npixel')
    glc2000_tot_npixel = 0

     allocate (glc2000_class_fraction(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_glc2000), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array glc2000_class_fraction')
    glc2000_class_fraction = 0.0


     allocate (glc2000_class_npixel(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_glc2000), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array glc2000_class_npixel')
    glc2000_class_npixel = 0

     allocate (ice_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array ice_glc2000')
    ice_glc2000 = 0.0

     allocate (z0_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array z0_glc2000')
    z0_glc2000 = 0.0

     allocate (root_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array root_glc2000')
    root_glc2000 = 0.0

     allocate (plcov_mx_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array plcov_mx_glc2000')
    plcov_mx_glc2000 = 0.0


     allocate (plcov_mn_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array plcov_mn_glc2000')
    plcov_mn_glc2000 = 0.0

     allocate (lai_mx_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lai_mx_glc2000')
    lai_mx_glc2000 = 0.0

     allocate (lai_mn_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lai_mn_glc2000')
    lai_mn_glc2000 = 0.0

     allocate (rs_min_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array rs_min_glc2000')
    rs_min_glc2000 = 0.0

     allocate (urban_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array urban_glc2000')
    urban_glc2000 = 0.0


    allocate (for_d_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array for_d_glc2000')
    for_d_glc2000 = 0.0

    allocate (for_e_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array for_e_glc2000')
    for_e_glc2000 = 0.0

    allocate (emissivity_glc2000(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array emissivity_glc2000')
    emissivity_glc2000 = 0.0

  end subroutine allocate_glc2000_target_fields



END Module mo_glc2000_tg_fields

