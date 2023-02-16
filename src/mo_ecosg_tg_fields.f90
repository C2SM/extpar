!+ Fortran module for ecosg data specification on target grid for external Parameters
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2022/09/01 Andrzej Wyszogrodzki
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for ecosg data specification on target grid for external Parameters
!> \author AndrzeJ Wyszogrodzki
!
MODULE mo_ecosg_tg_fields

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  USE mo_ecosg_lookup_tables, ONLY: nclass_ecosg

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fr_land_ecosg, &
       &    ecosg_class_fraction,    &
       &    ecosg_class_npixel, &
       &    ecosg_tot_npixel, &
       &    ice_ecosg, &
       &    z0_ecosg, &
       &    root_ecosg, &
       &    plcov_mn_ecosg, &
       &    plcov_mx_ecosg, &
       &    lai_mn_ecosg, &
       &    lai_mx_ecosg, &
       &    rs_min_ecosg, &
       &    urban_ecosg,  &
       &    for_d_ecosg,  &
       &    for_e_ecosg, &
       &    skinc_ecosg, &
       &    emissivity_ecosg, &
       &    allocate_ecosg_target_fields



       INTEGER (KIND=i4), ALLOCATABLE :: ecosg_class_npixel(:,:,:,:), & !< number of raw data pixels for each &
  ! &                                    ecosg class on target grid (dimension (ie,je,ke,nclass_ecosg))
            &                            ecosg_tot_npixel(:,:,:)  !< total number of ecosg raw data pixels &
  ! &                                    on target grid (dimension (ie,je,ke))


       REAL (KIND=wp), ALLOCATABLE  :: ecosg_class_fraction(:,:,:,:), &  !< fraction for each ecosg class &
  ! &                                  on target grid (dimension (ie,je,ke,nclass_ecosg))
            &                          fr_land_ecosg(:,:,:), & !< fraction land due to ecosg raw data
            &                          ice_ecosg(:,:,:), &     !< fraction of ice due to ecosg raw data
            &                          z0_ecosg(:,:,:), &      !< roughness length due to ecosg land use data
            &                          root_ecosg(:,:,:), &    !< root depth due to ecosg land use data
            &                          plcov_mx_ecosg(:,:,:), &!< plant cover maximum due to ecosg land use data
            &                          plcov_mn_ecosg(:,:,:), &!< plant cover minimum due to ecosg land use data
            &                          lai_mx_ecosg(:,:,:), &  !< Leaf Area Index maximum due to ecosg land use data
            &                          lai_mn_ecosg(:,:,:), &  !< Leaf Area Index minimum due to ecosg land use data
            &                          rs_min_ecosg(:,:,:), &  !< minimal stomata resistance due to ecosg land use data
            &                          urban_ecosg(:,:,:), &   !< urban fraction due to ecosg land use data
            &                          for_d_ecosg(:,:,:), &   !< deciduous forest (fraction) due to ecosg land use data
            &                          for_e_ecosg(:,:,:), &   !< evergreen forest (fraction) due to ecosg land use data
            &                          skinc_ecosg(:,:,:), &   !< skin conductivity due to ecosg land use data
            &                          emissivity_ecosg(:,:,:) !< longwave emissivity due to ecosg land use data


  CONTAINS

  !> allocate fields for TARGET grid
  !!
  !! the target grid for the GME has 3 dimension (ie,je,jd),
  !! the target grid for the COSMO model has 2 dimension (ie,je)
  !! the target grid for the ICON model has 1 dimension (ne)
  !! depending of the target model the second and third dimension of the target fields should be
  !! allocated with the length 1
  SUBROUTINE allocate_ecosg_target_fields(tg)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description


    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0

    CALL logging%info('Enter routine: allocate_ecosg_target_fields')

    ALLOCATE (fr_land_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_ecosg',__FILE__,__LINE__)
    fr_land_ecosg = 0.0

    ALLOCATE (ecosg_tot_npixel(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ecosg_tot_npixel',__FILE__,__LINE__)
    ecosg_tot_npixel = 0

    ALLOCATE (ecosg_class_fraction(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_ecosg), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ecosg_class_fraction',__FILE__,__LINE__)
    ecosg_class_fraction = 0.0


    ALLOCATE (ecosg_class_npixel(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_ecosg), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ecosg_class_npixel',__FILE__,__LINE__)
    ecosg_class_npixel = 0

    ALLOCATE (ice_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ice_ecosg',__FILE__,__LINE__)
    ice_ecosg = 0.0

    ALLOCATE (z0_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0_ecosg',__FILE__,__LINE__)
    z0_ecosg = 0.0

    ALLOCATE (root_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array root_ecosg',__FILE__,__LINE__)
    root_ecosg = 0.0

    ALLOCATE (plcov_mx_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mx_ecosg',__FILE__,__LINE__)
    plcov_mx_ecosg = 0.0


    ALLOCATE (plcov_mn_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mn_ecosg',__FILE__,__LINE__)
    plcov_mn_ecosg = 0.0

    ALLOCATE (lai_mx_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mx_ecosg',__FILE__,__LINE__)
    lai_mx_ecosg = 0.0

    ALLOCATE (lai_mn_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mn_ecosg',__FILE__,__LINE__)
    lai_mn_ecosg = 0.0

    ALLOCATE (rs_min_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array rs_min_ecosg',__FILE__,__LINE__)
    rs_min_ecosg = 0.0

    ALLOCATE (urban_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array urban_ecosg',__FILE__,__LINE__)
    urban_ecosg = 0.0

    ALLOCATE (for_d_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_d_ecosg',__FILE__,__LINE__)
    for_d_ecosg = 0.0

    ALLOCATE (for_e_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_e_ecosg',__FILE__,__LINE__)
    for_e_ecosg = 0.0

    ALLOCATE (skinc_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array skinc_ecosg',__FILE__,__LINE__)
    skinc_ecosg = 0.0

    ALLOCATE (emissivity_ecosg(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emissivity_ecosg',__FILE__,__LINE__)
    emissivity_ecosg = 0.0

    CALL logging%info('Exit routine: allocate_ecosg_target_fields')

  END SUBROUTINE allocate_ecosg_target_fields

END MODULE mo_ecosg_tg_fields
