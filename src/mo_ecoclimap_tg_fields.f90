!+ Fortran module for ecoclimap data specification on target grid for external Parameters
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
!> Fortran module for ecoclimap data specification on target grid for external Parameters
!> \author Hermann Asensio
!! ECOCLIMAP option gs_08.03.12  programm not operational yet

! Description:
! The ecoclimap dataset contains the following land use classification scheme
!

! class no. value         description
!

MODULE mo_ecoclimap_tg_fields

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  USE mo_ecoclimap_lookup_tables, ONLY: nclass_ecoclimap

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fr_land_ecoclimap, &
       &    ecoclimap_class_fraction,    &
       &    ecoclimap_class_npixel, &
       &    ecoclimap_tot_npixel, &
       &    ice_ecoclimap, &
       &    z012_ecoclimap, &
       &    z012tot_ecoclimap, &
       &    root_ecoclimap, &
       &    plcov12_ecoclimap, &
       &    lai12_ecoclimap, &
       &    rs_min_ecoclimap, &
       &    urban_ecoclimap,  &
       &    for_d_ecoclimap,  &
       &    for_e_ecoclimap, &
       &    emissivity_ecoclimap, &
       &    allocate_ecoclimap_target_fields

  !< fraction for each ecoclimap class on target grid (dimension (ie,je,ke,nclass_ecoclimap))
  REAL (KIND=wp), ALLOCATABLE  :: ecoclimap_class_fraction(:,:,:,:)  

  INTEGER (KIND=i4), ALLOCATABLE :: ecoclimap_class_npixel(:,:,:,:), & 
       &                            ecoclimap_tot_npixel(:,:,:)  


  REAL (KIND=wp), ALLOCATABLE  :: fr_land_ecoclimap(:,:,:), & !< fraction land due to ecoclimap raw data
       &                          ice_ecoclimap(:,:,:), &     !< fraction of ice due to ecoclimap raw data
       &                          z012_ecoclimap(:,:,:,:), &      !< roughness length due to ecoclimap land use data
       &                          z012tot_ecoclimap(:,:,:,:), &      !< roughness length 
       &                          root_ecoclimap(:,:,:), &    !< root depth due to ecoclimap land use data
       &                          plcov12_ecoclimap(:,:,:,:), & !< plant cover  to ecoclimap land use data
       &                          lai12_ecoclimap(:,:,:,:), &   !< Leaf Area Index  due to ecoclimap land use data
       &                          rs_min_ecoclimap(:,:,:), &  !< minimal stomata resistance due to ecoclimap land use data
       &                          urban_ecoclimap(:,:,:), &   !< urban fraction due to ecoclimap land use data
       &                          for_d_ecoclimap(:,:,:), &   !< deciduous forest (fraction) due to ecoclimap land use data
       &                          for_e_ecoclimap(:,:,:), &   !< evergreen forest (fraction) due to ecoclimap land use data
       &                          emissivity_ecoclimap(:,:,:) !< longwave emissivity due to ecoclimap land use data


  CONTAINS

  !> allocate fields for TARGET grid
  !!
  !! the target grid for the GME has 3 dimension (ie,je,jd),
  !! the target grid for the COSMO model has 2 dimension (ie,je)
  !! the target grid for the ICON model has 1 dimension (ne)
  !! depending of the target model the second and third dimension of the target fields should be
  !! allocated with the length 1
  SUBROUTINE allocate_ecoclimap_target_fields(tg)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description


    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_ecoclimap_target_fields')

    ALLOCATE (fr_land_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_ecoclimap',__FILE__,__LINE__)
    fr_land_ecoclimap = 0.0

    ALLOCATE (ecoclimap_tot_npixel(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ecoclimap_tot_npixel',__FILE__,__LINE__)
    ecoclimap_tot_npixel = 0

    ALLOCATE (ecoclimap_class_fraction(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_ecoclimap), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ecoclimap_class_fraction',__FILE__,__LINE__)
    ecoclimap_class_fraction = 0.0

    ALLOCATE (ecoclimap_class_npixel(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_ecoclimap), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ecoclimap_class_npixel',__FILE__,__LINE__)
    ecoclimap_class_npixel = 0

    ALLOCATE (ice_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ice_ecoclimap',__FILE__,__LINE__)
    ice_ecoclimap = 0.0

    ALLOCATE (z012_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke,1:12), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z012_ecoclimap',__FILE__,__LINE__)
    z012_ecoclimap = 0.0

    ALLOCATE (z012tot_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke,1:12), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z012tot_ecoclimap',__FILE__,__LINE__)
    z012tot_ecoclimap = 0.0


    ALLOCATE (root_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array root_ecoclimap',__FILE__,__LINE__)
    root_ecoclimap = 0.0

    ALLOCATE (plcov12_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke,1:12), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov12_ecoclimap',__FILE__,__LINE__)
    plcov12_ecoclimap = 0.0

    ALLOCATE (lai12_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke,1:12), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai12_ecoclimap',__FILE__,__LINE__)
    lai12_ecoclimap = 0.0


    ALLOCATE (rs_min_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array rs_min_ecoclimap',__FILE__,__LINE__)
    rs_min_ecoclimap = 0.0

    ALLOCATE (urban_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array urban_ecoclimap',__FILE__,__LINE__)
    urban_ecoclimap = 0.0


    ALLOCATE (for_d_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_d_ecoclimap',__FILE__,__LINE__)
    for_d_ecoclimap = 0.0

    ALLOCATE (for_e_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_e_ecoclimap',__FILE__,__LINE__)
    for_e_ecoclimap = 0.0

    ALLOCATE (emissivity_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emissivity_ecoclimap',__FILE__,__LINE__)
    emissivity_ecoclimap = 0.0

    CALL logging%info('Exit routine: allocate_ecoclimap_target_fields')

  END SUBROUTINE allocate_ecoclimap_target_fields

END MODULE mo_ecoclimap_tg_fields
