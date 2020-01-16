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

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8


  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_grid_structures, ONLY: target_grid_def

  USE mo_ecoclimap_lookup_tables, ONLY: nclass_ecoclimap

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fr_land_ecoclimap, &
       ecoclimap_class_fraction,    &
       ecoclimap_class_npixel, &
       ecoclimap_tot_npixel, &
       ice_ecoclimap, &
       z012_ecoclimap, &
       z012tot_ecoclimap, &
       root_ecoclimap, &
       plcov12_ecoclimap, &
       lai12_ecoclimap, &
       rs_min_ecoclimap, &
       urban_ecoclimap,  &
       for_d_ecoclimap,  &
       for_e_ecoclimap, &
       emissivity_ecoclimap, &
       allocate_ecoclimap_target_fields

  !< fraction for each ecoclimap class on target grid (dimension (ie,je,ke,nclass_ecoclimap))
  REAL (KIND=wp), ALLOCATABLE  :: ecoclimap_class_fraction(:,:,:,:)  
  !< number of raw data pixels for each ecoclimap class on target grid (dimension (ie,je,ke,nclass_ecoclimap))
  INTEGER (KIND=i8), ALLOCATABLE :: ecoclimap_class_npixel(:,:,:,:) 
  !< total number of ecoclimap raw data pixels on target grid (dimension (ie,je,ke))
  INTEGER (KIND=i8), ALLOCATABLE :: ecoclimap_tot_npixel(:,:,:)  


  REAL (KIND=wp), ALLOCATABLE  :: fr_land_ecoclimap(:,:,:) !< fraction land due to ecoclimap raw data
  REAL (KIND=wp), ALLOCATABLE  :: ice_ecoclimap(:,:,:)     !< fraction of ice due to ecoclimap raw data
  REAL (KIND=wp), ALLOCATABLE  :: z012_ecoclimap(:,:,:,:)      !< roughness length due to ecoclimap land use data
  REAL (KIND=wp), ALLOCATABLE  :: z012tot_ecoclimap(:,:,:,:)      !< roughness length 
  REAL (KIND=wp), ALLOCATABLE  :: root_ecoclimap(:,:,:)    !< root depth due to ecoclimap land use data
  REAL (KIND=wp), ALLOCATABLE  :: plcov12_ecoclimap(:,:,:,:) !< plant cover  to ecoclimap land use data
  REAL (KIND=wp), ALLOCATABLE  :: lai12_ecoclimap(:,:,:,:)   !< Leaf Area Index  due to ecoclimap land use data
  REAL (KIND=wp), ALLOCATABLE  :: rs_min_ecoclimap(:,:,:)  !< minimal stomata resistance due to ecoclimap land use data
  REAL (KIND=wp), ALLOCATABLE  :: urban_ecoclimap(:,:,:)   !< urban fraction due to ecoclimap land use data
  REAL (KIND=wp), ALLOCATABLE  :: for_d_ecoclimap(:,:,:)   !< deciduous forest (fraction) due to ecoclimap land use data
  REAL (KIND=wp), ALLOCATABLE  :: for_e_ecoclimap(:,:,:)   !< evergreen forest (fraction) due to ecoclimap land use data
  REAL (KIND=wp), ALLOCATABLE  :: emissivity_ecoclimap(:,:,:) !< longwave emissivity due to ecoclimap land use data


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


    INTEGER :: errorcode !< error status variable


    ALLOCATE (fr_land_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_land_ecoclimap')
    fr_land_ecoclimap = 0.0

    ALLOCATE (ecoclimap_tot_npixel(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array ecoclimap_tot_npixel')
    ecoclimap_tot_npixel = 0

    ALLOCATE (ecoclimap_class_fraction(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_ecoclimap), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array ecoclimap_class_fraction')
    ecoclimap_class_fraction = 0.0


    ALLOCATE (ecoclimap_class_npixel(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_ecoclimap), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array ecoclimap_class_npixel')
    ecoclimap_class_npixel = 0

    ALLOCATE (ice_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array ice_ecoclimap')
    ice_ecoclimap = 0.0

    ALLOCATE (z012_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke,1:12), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array z012_ecoclimap')
    z012_ecoclimap = 0.0

    ALLOCATE (z012tot_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke,1:12), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array z012tot_ecoclimap')
    z012tot_ecoclimap = 0.0


    ALLOCATE (root_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array root_ecoclimap')
    root_ecoclimap = 0.0

    ALLOCATE (plcov12_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke,1:12), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array plcov12_ecoclimap')
    plcov12_ecoclimap = 0.0

    ALLOCATE (lai12_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke,1:12), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lai12_ecoclimap')
    lai12_ecoclimap = 0.0


    ALLOCATE (rs_min_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array rs_min_ecoclimap')
    rs_min_ecoclimap = 0.0

    ALLOCATE (urban_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array urban_ecoclimap')
    urban_ecoclimap = 0.0


    ALLOCATE (for_d_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array for_d_ecoclimap')
    for_d_ecoclimap = 0.0

    ALLOCATE (for_e_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array for_e_ecoclimap')
    for_e_ecoclimap = 0.0

    ALLOCATE (emissivity_ecoclimap(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array emissivity_ecoclimap')
    emissivity_ecoclimap = 0.0

  END SUBROUTINE allocate_ecoclimap_target_fields

END Module mo_ecoclimap_tg_fields


