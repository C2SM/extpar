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

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8


!> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_grid_structures, ONLY: target_grid_def

  USE mo_glcc_lookup_tables, ONLY: nclass_glcc

IMPLICIT NONE

PRIVATE

PUBLIC :: fr_land_glcc, &
  &       glcc_class_fraction,    &
  &        glcc_class_npixel, &
  &        glcc_tot_npixel, &
  &        ice_glcc, &
  &        z0_glcc, &
  &        root_glcc, &
  &        plcov_mn_glcc, &
  &        plcov_mx_glcc, &
  &        lai_mn_glcc, &
  &        lai_mx_glcc, &
  &        rs_min_glcc, &
  &        urban_glcc,  &
  &        for_d_glcc,  &
  &        for_e_glcc, &
  &        emissivity_glcc, &
  &        allocate_glcc_target_fields


       REAL (KIND=wp), ALLOCATABLE  :: glcc_class_fraction(:,:,:,:)  
!< fraction for each glcc class on target grid (dimension (ie,je,ke,nclass_glcc))

       INTEGER (KIND=i8), ALLOCATABLE :: glcc_class_npixel(:,:,:,:) 
!< number of raw data pixels for each glcc class on target grid (dimension (ie,je,ke,nclass_glcc))


       INTEGER (KIND=i8), ALLOCATABLE :: glcc_tot_npixel(:,:,:)  
!< total number of glcc raw data pixels on target grid (dimension (ie,je,ke))


       REAL (KIND=wp), ALLOCATABLE  :: fr_land_glcc(:,:,:) !< fraction land due to glcc raw data
       REAL (KIND=wp), ALLOCATABLE  :: ice_glcc(:,:,:)     !< fraction of ice due to glcc raw data
       REAL (KIND=wp), ALLOCATABLE  :: z0_glcc(:,:,:)      !< roughness length due to glcc land use data
       REAL (KIND=wp), ALLOCATABLE  :: root_glcc(:,:,:)    !< root depth due to glcc land use data
       REAL (KIND=wp), ALLOCATABLE  :: plcov_mx_glcc(:,:,:)!< plant cover maximum due to glcc land use data
       REAL (KIND=wp), ALLOCATABLE  :: plcov_mn_glcc(:,:,:)!< plant cover minimum due to glcc land use data
       REAL (KIND=wp), ALLOCATABLE  :: lai_mx_glcc(:,:,:)  !< Leaf Area Index maximum due to glcc land use data
       REAL (KIND=wp), ALLOCATABLE  :: lai_mn_glcc(:,:,:)  !< Leaf Area Index minimum due to glcc land use data
       REAL (KIND=wp), ALLOCATABLE  :: rs_min_glcc(:,:,:)  !< minimal stomata resistance due to glcc land use data
       REAL (KIND=wp), ALLOCATABLE  :: urban_glcc(:,:,:)   !< urban fraction due to glcc land use data
       REAL (KIND=wp), ALLOCATABLE  :: for_d_glcc(:,:,:)   !< deciduous forest (fraction) due to glcc land use data
       REAL (KIND=wp), ALLOCATABLE  :: for_e_glcc(:,:,:)   !< evergreen forest (fraction) due to glcc land use data
       REAL (KIND=wp), ALLOCATABLE  :: emissivity_glcc(:,:,:) !< longwave emissivity due to glcc land use data


CONTAINS





!> allocate fields for TARGET grid
!!
!! the target grid for the GME has 3 dimension (ie,je,jd),
!! the target grid for the COSMO model has 2 dimension (ie,je)
!! the target grid for the ICON model has 1 dimension (ne)
!! depending of the target model the second and third dimension of the target fields should be 
!! allocated with the length 1
  subroutine allocate_glcc_target_fields(tg)
  

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description


    INTEGER :: errorcode !< error status variable

   
    ALLOCATE (fr_land_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array fr_land_glcc')
    fr_land_glcc = 0.0

    allocate (glcc_tot_npixel(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array glcc_tot_npixel')
    glcc_tot_npixel = 0

     allocate (glcc_class_fraction(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_glcc), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array glcc_class_fraction')
    glcc_class_fraction = 0.0


     allocate (glcc_class_npixel(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_glcc), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array glcc_class_npixel')
    glcc_class_npixel = 0

     allocate (ice_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array ice_glcc')
    ice_glcc = 0.0

     allocate (z0_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array z0_glcc')
    z0_glcc = 0.0

     allocate (root_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array root_glcc')
    root_glcc = 0.0

     allocate (plcov_mx_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array plcov_mx_glcc')
    plcov_mx_glcc = 0.0


     allocate (plcov_mn_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array plcov_mn_glcc')
    plcov_mn_glcc = 0.0

     allocate (lai_mx_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lai_mx_glcc')
    lai_mx_glcc = 0.0

     allocate (lai_mn_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lai_mn_glcc')
    lai_mn_glcc = 0.0

     allocate (rs_min_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array rs_min_glcc')
    rs_min_glcc = 0.0

     allocate (urban_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array urban_glcc')
    urban_glcc = 0.0


    allocate (for_d_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array for_d_glcc')
    for_d_glcc = 0.0

    allocate (for_e_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array for_e_glcc')
    for_e_glcc = 0.0

    allocate (emissivity_glcc(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array emissivity_glcc')
    emissivity_glcc = 0.0

  end subroutine allocate_glcc_target_fields



END Module mo_glcc_tg_fields

