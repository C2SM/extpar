!+ Fortran module for land use data on target grid for external parameters 
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_3         2011/04/19 Hermann Asensio
! introduce Globcover 2009 land use data set for external parameters
! V1_11        2013/04/16 Juergen Helmert
!  Adaptions for using external land-sea mask 
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for land use data on target grid for external parameters 
!> \author Hermann Asensio
!
MODULE mo_lu_tg_fields

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fr_land_lu, &
       &    fr_land_mask, &
       &    ice_lu, &
       &    z0_lu, &
       &    z0_tot, &
       &    root_lu, &
       &    plcov_mn_lu, &
       &    plcov_mx_lu, &
       &    lai_mn_lu, &
       &    lai_mx_lu, &
       &    rs_min_lu, &
       &    urban_lu,  &
       &    for_d_lu,  &
       &    for_e_lu, &
       &    skinc_lu, &
       &    emissivity_lu, &
       &    fr_ocean_lu, &
       &    lu_class_fraction,    &
       &    lu_class_npixel, &
       &    lu_tot_npixel, &
       &    z012_tot,      &
       &    z012_lu,       &
       &    lai12_lu,         &
       &    plcov12_lu


  PUBLIC :: allocate_lu_target_fields, allocate_add_lu_fields
  PUBLIC :: i_lu_globcover, i_lu_glc2000, i_lu_glcc, i_lu_ecoclimap

  PUBLIC :: fr_land, &
       &    ice, &
       &    z0, &
       &    root, &
       &    plcov_mn, &
       &    plcov_mx, &
       &    lai_mn, &
       &    lai_mx, &
       &    rs_min, &
       &    urban,  &
       &    for_d,  &
       &    for_e, &
       &    skinc, &
       &    emissivity, &
       &    fr_ocean, &
       &    z012, &
       &    z012tot, & 
       &    lai12, &
       &    plcov12

  PUBLIC :: allocate_lu_ds_target_fields

  INTEGER (KIND=i4), ALLOCATABLE :: lu_class_npixel(:,:,:,:), &  
       &                            lu_tot_npixel(:,:,:)  

  INTEGER(KIND=i4), PARAMETER    :: i_lu_globcover = 1, &  !< id for landuse data set Globcover 2009
       &                            i_lu_glc2000   = 2, &  !< id for landuse data set GLC2000
       &                            i_lu_glcc      = 3, &  !< id for landuse data set GLCC
       &                            i_lu_ecoclimap = 4  !< id for landuse data set ecoclimap

  REAL (KIND=wp), ALLOCATABLE    :: fr_land_lu(:,:,:), &  !< fraction land due to land use raw data
       &                            fr_land_mask(:,:,:), &  !< fraction land due to external target data
       &                            ice_lu(:,:,:), &      !< fraction of ice due to land use raw data
       &                            z0_lu(:,:,:), &       !< roughness length due to land use land use data
       &                            z0_tot(:,:,:), &       !< total roughness length 
       &                            root_lu(:,:,:), &     !< root depth due to land use land use data
       &                            plcov_mx_lu(:,:,:), & !< plant cover maximum due to land use land use data
       &                            plcov_mn_lu(:,:,:), & !< plant cover minimum due to land use land use data
       &                            lai_mx_lu(:,:,:), &   !< Leaf Area Index maximum due to land use land use data
       &                            lai_mn_lu(:,:,:), &   !< Leaf Area Index minimum due to land use land use data
       &                            rs_min_lu(:,:,:), &   !< minimal stomata resistance due to land use land use data
       &                            urban_lu(:,:,:), &    !< urban fraction due to land use land use data
       &                            for_d_lu(:,:,:), &    !< deciduous forest (fraction) due to land use land use data
       &                            for_e_lu(:,:,:), &    !< evergreen forest (fraction) due to land use land use data
       &                            skinc_lu(:,:,:), &    !< skin conductivity due to land use data
       &                            emissivity_lu(:,:,:), &  !< longwave emissivity due to land use land use data
       &                            z012_lu(:,:,:,:), &  !< z0 veget. ecoclomap
       &                            z012_tot(:,:,:,:), &  !< z0 ecoclomap 
       &                            lai12_lu(:,:,:,:), &  ! <  lai12 ecoclimap
       &                            plcov12_lu(:,:,:,:), &  !<  plcov ecoclimap
       &                            fr_ocean_lu(:,:,:), &  !< fraction ocean due to land use raw data
       &                            lu_class_fraction(:,:,:,:), &   
       &                            fr_land(:,:,:,:), &  !< fraction land due to land use raw data
       &                            ice(:,:,:,:), &      !< fraction of ice due to land use raw data
       &                            z0(:,:,:,:), &       !< roughness length due to land use land use data
       &                            root(:,:,:,:), &     !< root depth due to land use data
       &                            plcov_mx(:,:,:,:), & !< plant cover maximum due to land use data
       &                            plcov_mn(:,:,:,:), & !< plant cover minimum due to land use data
       &                            lai_mx(:,:,:,:), &   !< Leaf Area Index maximum due to land use data
       &                            lai_mn(:,:,:,:), &   !< Leaf Area Index minimum due to land use data
       &                            rs_min(:,:,:,:), &   !< minimal stomata resistance due to land use data
       &                            urban(:,:,:,:), &    !< urban fraction due to land use data
       &                            for_d(:,:,:,:), &    !< deciduous forest (fraction) due to land use data
       &                            for_e(:,:,:,:), &    !< evergreen forest (fraction) due to land use data
       &                            skinc(:,:,:,:), &    !< skin conductivity due to land use data
       &                            emissivity(:,:,:,:), &  !< longwave emissivity due to land use data
       &                            fr_ocean(:,:,:,:), &  !< fraction ocean due to land use raw data
       &                            z012(:,:,:,:), &  !< z0 ecoclomap
       &                            z012tot(:,:,:,:), &  !< z0 ecoclomap 
       &                            lai12(:,:,:,:), &  ! <  lai12 ecoclimap
       &                            plcov12(:,:,:,:) !<  plcov ecoclimap

  CONTAINS

  !> allocate fields for TARGET grid
  !! the target grid for the GME has 3 dimension (ie,je,jd),
  !! the target grid for the COSMO model has 2 dimension (ie,je)
  !! the target grid for the ICON model has 1 dimension (ne)
  !! depending of the target model the second and third dimension of the target fields should be 
  !! allocated with the length 1
  SUBROUTINE allocate_lu_target_fields(tg)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER(KIND=i4)                  :: errorcode !< error status variable
   
    CALL logging%info('Enter routine: allocate_lu_target_fields')

    ALLOCATE (fr_land_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_lu',__FILE__,__LINE__)
    fr_land_lu = 0.0

    ALLOCATE (fr_land_mask(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_mask',__FILE__,__LINE__)
    fr_land_mask = 0.0

    ALLOCATE (ice_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ice_lu',__FILE__,__LINE__)
    ice_lu = 0.0

    ALLOCATE (z0_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0_lu',__FILE__,__LINE__)
    z0_lu = 0.0

    ALLOCATE (z0_tot(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0_tot',__FILE__,__LINE__)
    z0_tot = 0.0

    ALLOCATE (root_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array root_lu',__FILE__,__LINE__)
    root_lu = 0.0

    ALLOCATE (plcov_mx_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mx_lu',__FILE__,__LINE__)
    plcov_mx_lu = 0.0

    ALLOCATE (plcov_mn_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mn_lu',__FILE__,__LINE__)
    plcov_mn_lu = 0.0
        
    ALLOCATE (lai_mx_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mx_lu',__FILE__,__LINE__)
    lai_mx_lu = 0.0

    ALLOCATE (lai_mn_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mn_lu',__FILE__,__LINE__)
    lai_mn_lu = 0.0

    ALLOCATE (rs_min_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array rs_min_lu',__FILE__,__LINE__)
    rs_min_lu = 0.0

    ALLOCATE (urban_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array urban_lu',__FILE__,__LINE__)
    urban_lu = 0.0

    ALLOCATE (for_d_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_d_lu',__FILE__,__LINE__)
    for_d_lu = 0.0

    ALLOCATE (for_e_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_e_lu',__FILE__,__LINE__)
    for_e_lu = 0.0

    ALLOCATE (skinc_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array skinc_lu',__FILE__,__LINE__)
    skinc_lu = 0.0

    ALLOCATE (emissivity_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emissivity_lu',__FILE__,__LINE__)
    emissivity_lu = 0.0

    ALLOCATE (fr_ocean_lu(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_ocean_lu',__FILE__,__LINE__)
    fr_ocean_lu = 0.0

    ALLOCATE (lai12_lu(1:tg%ie,1:tg%je,1:tg%ke,1:12), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai12_lu',__FILE__,__LINE__)
    lai12_lu = 0.0

    ALLOCATE (plcov12_lu(1:tg%ie,1:tg%je,1:tg%ke,1:12), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov12_lu',__FILE__,__LINE__)
    plcov12_lu = 0.0

    ALLOCATE (z012_lu(1:tg%ie,1:tg%je,1:tg%ke,1:12), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z012',__FILE__,__LINE__)
    z012_lu = 0.0

    ALLOCATE (z012_tot(1:tg%ie,1:tg%je,1:tg%ke,1:12), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z012_tot',__FILE__,__LINE__)
    z012_tot = 0.0

  END SUBROUTINE allocate_lu_target_fields


  !> allocate additional land use target fields
  SUBROUTINE allocate_add_lu_fields(tg,nclass_lu)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER, INTENT(IN)               :: nclass_lu !< number of land use classes
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_add_lu_fields')

    ALLOCATE (lu_tot_npixel(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lu_tot_npixel',__FILE__,__LINE__)
    lu_tot_npixel = 0

    ALLOCATE (lu_class_fraction(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_lu), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lu_class_fraction',__FILE__,__LINE__)
    lu_class_fraction = 0.0

    ALLOCATE (lu_class_npixel(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_lu), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lu_class_npixel',__FILE__,__LINE__)
    lu_class_npixel = 0

  END SUBROUTINE allocate_add_lu_fields

  !> allocate land use fields for TARGET grid
  !!
  !! the target grid for the GME has 3 dimension (ie,je,jd),
  !! the target grid for the COSMO model has 2 dimension (ie,je)
  !! the target grid for the ICON model has 1 dimension (ne)
  !! depending of the target model the second and third dimension of the target fields should be 
  !! allocated with the length 1
  SUBROUTINE allocate_lu_ds_target_fields(tg,n_data)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER(KIND=i4), INTENT(IN)      :: n_data !< number of datasets for land use
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_lu_ds_target_fields')

    ALLOCATE (fr_land(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land',__FILE__,__LINE__)
    fr_land = 0.0

    ALLOCATE (ice(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ice',__FILE__,__LINE__)
    ice = 0.0

    ALLOCATE (z0(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0',__FILE__,__LINE__)
    z0 = 0.0

    ALLOCATE (root(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array root',__FILE__,__LINE__)
    root = 0.0

    ALLOCATE (plcov_mx(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mx',__FILE__,__LINE__)
    plcov_mx = 0.0

    ALLOCATE (plcov_mn(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mn',__FILE__,__LINE__)
    plcov_mn = 0.0

    ALLOCATE (lai_mx(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mx',__FILE__,__LINE__)
    lai_mx = 0.0

    ALLOCATE (lai_mn(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mn',__FILE__,__LINE__)
    lai_mn = 0.0

    ALLOCATE (rs_min(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array rs_min',__FILE__,__LINE__)
    rs_min = 0.0

    ALLOCATE (urban(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array urban',__FILE__,__LINE__)
    urban = 0.0

    ALLOCATE (for_d(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_d',__FILE__,__LINE__)
    for_d = 0.0

    ALLOCATE (for_e(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_e',__FILE__,__LINE__)
    for_e = 0.0

    ALLOCATE (skinc(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array skinc',__FILE__,__LINE__)
    skinc = 0.0

    ALLOCATE (emissivity(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emissivity',__FILE__,__LINE__)
    emissivity = 0.0

    ALLOCATE (fr_ocean(1:tg%ie,1:tg%je,1:tg%ke,1:n_data), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_ocean',__FILE__,__LINE__)
    fr_ocean = 0.0

  END SUBROUTINE allocate_lu_ds_target_fields

END MODULE mo_lu_tg_fields
