!+ Fortran module for GLCC data specification on target grid for external Parameters
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
!> Fortran module for GLCC data specification on target grid for external Parameters 
!> \author Hermann Asensio
!

MODULE mo_glcc_tg_fields

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  USE mo_glcc_lookup_tables,    ONLY: nclass_glcc

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fr_land_glcc, &
       &    glcc_class_fraction,    &
       &    glcc_class_npixel, &
       &    glcc_tot_npixel, &
       &    ice_glcc, &
       &    z0_glcc, &
       &    root_glcc, &
       &    plcov_mn_glcc, &
       &    plcov_mx_glcc, &
       &    lai_mn_glcc, &
       &    lai_mx_glcc, &
       &    rs_min_glcc, &
       &    urban_glcc,  &
       &    for_d_glcc,  &
       &    for_e_glcc, &
       &    emissivity_glcc, &
       &    allocate_glcc_target_fields

  INTEGER (KIND=i4), ALLOCATABLE :: glcc_class_npixel(:,:,:,:), &
       &                            glcc_tot_npixel(:,:,:)  

  REAL (KIND=wp), ALLOCATABLE  :: glcc_class_fraction(:,:,:,:), &   
       &                          fr_land_glcc(:,:,:), &  !< fraction land due to glcc raw data
       &                          ice_glcc(:,:,:), &      !< fraction of ice due to glcc raw data
       &                          z0_glcc(:,:,:), &       !< roughness length due to glcc land use data
       &                          root_glcc(:,:,:), &     !< root depth due to glcc land use data
       &                          plcov_mx_glcc(:,:,:), & !< plant cover maximum due to glcc land use data
       &                          plcov_mn_glcc(:,:,:), & !< plant cover minimum due to glcc land use data
       &                          lai_mx_glcc(:,:,:), &   !< Leaf Area Index maximum due to glcc land use data
       &                          lai_mn_glcc(:,:,:), &   !< Leaf Area Index minimum due to glcc land use data
       &                          rs_min_glcc(:,:,:), &   !< minimal stomata resistance due to glcc land use data
       &                          urban_glcc(:,:,:), &    !< urban fraction due to glcc land use data
       &                          for_d_glcc(:,:,:), &    !< deciduous forest (fraction) due to glcc land use data
       &                          for_e_glcc(:,:,:), &    !< evergreen forest (fraction) due to glcc land use data
       &                          emissivity_glcc(:,:,:) !< longwave emissivity due to glcc land use data

  CONTAINS

  !> allocate fields for TARGET grid
  !!
  !! the target grid for the GME has 3 dimension (ie,je,jd),
  !! the target grid for the COSMO model has 2 dimension (ie,je)
  !! the target grid for the ICON model has 1 dimension (ne)
  !! depending of the target model the second and third dimension of the target fields should be 
  !! allocated with the length 1
  SUBROUTINE allocate_glcc_target_fields(tg)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_glcc_target_fields')
   
    ALLOCATE (fr_land_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_glcc',__FILE__,__LINE__)
    fr_land_glcc = 0.0

    ALLOCATE (glcc_tot_npixel(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array glcc_tot_npixel',__FILE__,__LINE__)
    glcc_tot_npixel = 0

    ALLOCATE (glcc_class_fraction(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_glcc), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array glcc_class_fraction',__FILE__,__LINE__)
    glcc_class_fraction = 0.0


    ALLOCATE (glcc_class_npixel(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_glcc), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array glcc_class_npixel',__FILE__,__LINE__)
    glcc_class_npixel = 0

    ALLOCATE (ice_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ice_glcc',__FILE__,__LINE__)
    ice_glcc = 0.0

    ALLOCATE (z0_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0_glcc',__FILE__,__LINE__)
    z0_glcc = 0.0

    ALLOCATE (root_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array root_glcc',__FILE__,__LINE__)
    root_glcc = 0.0

    ALLOCATE (plcov_mx_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mx_glcc',__FILE__,__LINE__)
    plcov_mx_glcc = 0.0


    ALLOCATE (plcov_mn_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mn_glcc',__FILE__,__LINE__)
    plcov_mn_glcc = 0.0

    ALLOCATE (lai_mx_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mx_glcc',__FILE__,__LINE__)
    lai_mx_glcc = 0.0

    ALLOCATE (lai_mn_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mn_glcc',__FILE__,__LINE__)
    lai_mn_glcc = 0.0

    ALLOCATE (rs_min_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array rs_min_glcc',__FILE__,__LINE__)
    rs_min_glcc = 0.0

    ALLOCATE (urban_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array urban_glcc',__FILE__,__LINE__)
    urban_glcc = 0.0

    ALLOCATE (for_d_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_d_glcc',__FILE__,__LINE__)
    for_d_glcc = 0.0

    ALLOCATE (for_e_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_e_glcc',__FILE__,__LINE__)
    for_e_glcc = 0.0

    ALLOCATE (emissivity_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emissivity_glcc',__FILE__,__LINE__)
    emissivity_glcc = 0.0

  END SUBROUTINE allocate_glcc_target_fields

END Module mo_glcc_tg_fields
