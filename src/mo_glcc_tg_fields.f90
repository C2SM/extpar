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
  USE mo_array_cache,           ONLY: allocate_cached
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

  INTEGER (KIND=i4), POINTER :: glcc_class_npixel(:,:,:,:), &
       &                        glcc_tot_npixel(:,:,:)

  REAL (KIND=wp), POINTER    :: glcc_class_fraction(:,:,:,:), &
       &                        fr_land_glcc(:,:,:), &  !< fraction land due to glcc raw data
       &                        ice_glcc(:,:,:), &      !< fraction of ice due to glcc raw data
       &                        z0_glcc(:,:,:), &       !< roughness length due to glcc land use data
       &                        root_glcc(:,:,:), &     !< root depth due to glcc land use data
       &                        plcov_mx_glcc(:,:,:), & !< plant cover maximum due to glcc land use data
       &                        plcov_mn_glcc(:,:,:), & !< plant cover minimum due to glcc land use data
       &                        lai_mx_glcc(:,:,:), &   !< Leaf Area Index maximum due to glcc land use data
       &                        lai_mn_glcc(:,:,:), &   !< Leaf Area Index minimum due to glcc land use data
       &                        rs_min_glcc(:,:,:), &   !< minimal stomata resistance due to glcc land use data
       &                        urban_glcc(:,:,:), &    !< urban fraction due to glcc land use data
       &                        for_d_glcc(:,:,:), &    !< deciduous forest (fraction) due to glcc land use data
       &                        for_e_glcc(:,:,:), &    !< evergreen forest (fraction) due to glcc land use data
       &                        emissivity_glcc(:,:,:) !< longwave emissivity due to glcc land use data

  CONTAINS

  !> allocate fields for TARGET grid
  !!
  !! the target grid for the GME has 3 dimension (ie,je,jd),
  !! the target grid for the COSMO model has 2 dimension (ie,je)
  !! the target grid for the ICON model has 1 dimension (ne)
  !! depending of the target model the second and third dimension of the target fields should be
  !! allocated with the length 1
  SUBROUTINE allocate_glcc_target_fields(tg, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    LOGICAL, INTENT(in)               :: l_use_array_cache 

    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_glcc_target_fields')

if (l_use_array_cache) then
   call allocate_cached('fr_land_glcc', fr_land_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(fr_land_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_glcc',__FILE__,__LINE__)
    fr_land_glcc = 0.0

if (l_use_array_cache) then
   call allocate_cached('glcc_tot_npixel', glcc_tot_npixel, [tg%ie,tg%je,tg%ke])
else
   allocate(glcc_tot_npixel(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array glcc_tot_npixel',__FILE__,__LINE__)
    glcc_tot_npixel = 0

if (l_use_array_cache) then
   call allocate_cached('glcc_class_fraction', glcc_class_fraction, [tg%ie,tg%je,tg%ke,nclass_glcc])
else
   allocate(glcc_class_fraction(tg%ie,tg%je,tg%ke,nclass_glcc), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array glcc_class_fraction',__FILE__,__LINE__)
    glcc_class_fraction = 0.0


if (l_use_array_cache) then
   call allocate_cached('glcc_class_npixel', glcc_class_npixel, [tg%ie,tg%je,tg%ke,nclass_glcc])
else
   allocate(glcc_class_npixel(tg%ie,tg%je,tg%ke,nclass_glcc), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array glcc_class_npixel',__FILE__,__LINE__)
    glcc_class_npixel = 0

if (l_use_array_cache) then
   call allocate_cached('ice_glcc', ice_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(ice_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ice_glcc',__FILE__,__LINE__)
    ice_glcc = 0.0

if (l_use_array_cache) then
   call allocate_cached('z0_glcc', z0_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(z0_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0_glcc',__FILE__,__LINE__)
    z0_glcc = 0.0

if (l_use_array_cache) then
   call allocate_cached('root_glcc', root_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(root_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array root_glcc',__FILE__,__LINE__)
    root_glcc = 0.0

if (l_use_array_cache) then
   call allocate_cached('plcov_mx_glcc', plcov_mx_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(plcov_mx_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mx_glcc',__FILE__,__LINE__)
    plcov_mx_glcc = 0.0


if (l_use_array_cache) then
   call allocate_cached('plcov_mn_glcc', plcov_mn_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(plcov_mn_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mn_glcc',__FILE__,__LINE__)
    plcov_mn_glcc = 0.0

if (l_use_array_cache) then
   call allocate_cached('lai_mx_glcc', lai_mx_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(lai_mx_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mx_glcc',__FILE__,__LINE__)
    lai_mx_glcc = 0.0

if (l_use_array_cache) then
   call allocate_cached('lai_mn_glcc', lai_mn_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(lai_mn_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mn_glcc',__FILE__,__LINE__)
    lai_mn_glcc = 0.0

if (l_use_array_cache) then
   call allocate_cached('rs_min_glcc', rs_min_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(rs_min_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array rs_min_glcc',__FILE__,__LINE__)
    rs_min_glcc = 0.0

if (l_use_array_cache) then
   call allocate_cached('urban_glcc', urban_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(urban_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array urban_glcc',__FILE__,__LINE__)
    urban_glcc = 0.0

if (l_use_array_cache) then
   call allocate_cached('for_d_glcc', for_d_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(for_d_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_d_glcc',__FILE__,__LINE__)
    for_d_glcc = 0.0

if (l_use_array_cache) then
   call allocate_cached('for_e_glcc', for_e_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(for_e_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_e_glcc',__FILE__,__LINE__)
    for_e_glcc = 0.0

if (l_use_array_cache) then
   call allocate_cached('emissivity_glcc', emissivity_glcc, [tg%ie,tg%je,tg%ke])
else
   allocate(emissivity_glcc(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emissivity_glcc',__FILE__,__LINE__)
    emissivity_glcc = 0.0

  END SUBROUTINE allocate_glcc_target_fields

END Module mo_glcc_tg_fields
