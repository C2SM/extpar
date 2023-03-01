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
  USE mo_array_cache,           ONLY: allocate_cached
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
       &    lu_tot_npixel


  PUBLIC :: allocate_lu_target_fields, allocate_add_lu_fields
  PUBLIC :: i_lu_globcover, i_lu_glc2000, i_lu_glcc, i_lu_ecci

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

  INTEGER (KIND=i4), POINTER     :: lu_class_npixel(:,:,:,:), &  
       &                            lu_tot_npixel(:,:,:)  

  INTEGER(KIND=i4), PARAMETER    :: i_lu_globcover = 1, &  !< id for landuse data set Globcover 2009
       &                            i_lu_glc2000   = 2, &  !< id for landuse data set GLC2000
       &                            i_lu_glcc      = 3, &  !< id for landuse data set GLCC
       &                            i_lu_ecci      = 5     !< id for landuse data set ESA CCI

 REAL (KIND=wp), POINTER         :: fr_land_lu(:,:,:), &  !< fraction land due to land use raw data
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
  SUBROUTINE allocate_lu_target_fields(tg, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    LOGICAL, INTENT(in)               :: l_use_array_cache
    
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0
    
    CALL logging%info('Enter routine: allocate_lu_target_fields')

if (l_use_array_cache) then
   call allocate_cached('fr_land_lu', fr_land_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache fr_land_lu')      
else
   allocate(fr_land_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_lu',__FILE__,__LINE__)
    fr_land_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('fr_land_mask', fr_land_mask, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache fr_land_mask')      
else
   allocate(fr_land_mask(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_mask',__FILE__,__LINE__)
    fr_land_mask = 0.0

if (l_use_array_cache) then
   call allocate_cached('ice_lu', ice_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache ice_lu')        
else
   allocate(ice_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ice_lu',__FILE__,__LINE__)
    ice_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('z0_lu', z0_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache z0_lu')      
else
   allocate(z0_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0_lu',__FILE__,__LINE__)
    z0_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('z0_tot', z0_tot, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache z0_tot')      
else
   allocate(z0_tot(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0_tot',__FILE__,__LINE__)
    z0_tot = 0.0

if (l_use_array_cache) then
   call allocate_cached('root_lu', root_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache root_lu')      
else
   allocate(root_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array root_lu',__FILE__,__LINE__)
    root_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('plcov_mx_lu', plcov_mx_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache plcov_mx_lu')      
else
   allocate(plcov_mx_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mx_lu',__FILE__,__LINE__)
    plcov_mx_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('plcov_mn_lu', plcov_mn_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache plcov_mn_lu')      
else
   allocate(plcov_mn_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mn_lu',__FILE__,__LINE__)
    plcov_mn_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('lai_mx_lu', lai_mx_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache lai_mx_lu')      
else
   allocate(lai_mx_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mx_lu',__FILE__,__LINE__)
    lai_mx_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('lai_mn_lu', lai_mn_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache lai_mn_lu')      
else
   allocate(lai_mn_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mn_lu',__FILE__,__LINE__)
    lai_mn_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('rs_min_lu', rs_min_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache rs_min_lu')      
else
   allocate(rs_min_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array rs_min_lu',__FILE__,__LINE__)
    rs_min_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('urban_lu', urban_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache urban_lu')      
else
   allocate(urban_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array urban_lu',__FILE__,__LINE__)
    urban_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('for_d_lu', for_d_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache for_d_lu')      
else
   allocate(for_d_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_d_lu',__FILE__,__LINE__)
    for_d_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('for_e_lu', for_e_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache for_e_lu')      
else
   allocate(for_e_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_e_lu',__FILE__,__LINE__)
    for_e_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('skinc_lu', skinc_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache skinc_lu')      
else
   allocate(skinc_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array skinc_lu',__FILE__,__LINE__)
    skinc_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('emissivity_lu', emissivity_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache emissivity_lu')      
else
   allocate(emissivity_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emissivity_lu',__FILE__,__LINE__)
    emissivity_lu = 0.0

if (l_use_array_cache) then
   call allocate_cached('fr_ocean_lu', fr_ocean_lu, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache fr_ocean_lu')      
else
   allocate(fr_ocean_lu(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_ocean_lu',__FILE__,__LINE__)
    fr_ocean_lu = 0.0

  END SUBROUTINE allocate_lu_target_fields


  !> allocate additional land use target fields
  SUBROUTINE allocate_add_lu_fields(tg, nclass_lu, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER, INTENT(in)               :: nclass_lu !< number of land use classes
    LOGICAL, INTENT(in)               :: l_use_array_cache 
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0
    
    CALL logging%info('Enter routine: allocate_add_lu_fields')
    write(message_text,'(a,l3)') '    Allocation scheme (false:allocate|true:cache): ', l_use_array_cache
    call logging%info(message_text)
    write(message_text,'(a,4(2x,i0))') '     Dimensions: ', tg%ie, tg%je, tg%ke, nclass_lu
    call logging%info(message_text)
    
if (l_use_array_cache) then
   call allocate_cached('lu_tot_npixel', lu_tot_npixel, [tg%ie,tg%je,tg%ke])
   CALL logging%info('cache lu_tot_npixel')
else
   allocate(lu_tot_npixel(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lu_tot_npixel',__FILE__,__LINE__)
    lu_tot_npixel = 0

if (l_use_array_cache) then
   call allocate_cached('lu_class_fraction', lu_class_fraction, [tg%ie,tg%je,tg%ke,nclass_lu])
   CALL logging%info('cache lu_class_fraction')
else
   allocate(lu_class_fraction(tg%ie,tg%je,tg%ke,nclass_lu), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lu_class_fraction',__FILE__,__LINE__)
    lu_class_fraction = 0.0

if (l_use_array_cache) then
   call allocate_cached('lu_class_npixel', lu_class_npixel, [tg%ie,tg%je,tg%ke,nclass_lu])
   CALL logging%info('cache lu_class_npixel')
else
   allocate(lu_class_npixel(tg%ie,tg%je,tg%ke,nclass_lu), stat=errorcode)
endif
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
  SUBROUTINE allocate_lu_ds_target_fields(tg,n_data, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER(KIND=i4), INTENT(IN)      :: n_data !< number of datasets for land use
    LOGICAL, INTENT(in)               :: l_use_array_cache 
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0
    
    CALL logging%info('Enter routine: allocate_lu_ds_target_fields')

if (l_use_array_cache) then
   call allocate_cached('fr_land', fr_land, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(fr_land(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land',__FILE__,__LINE__)
    fr_land = 0.0

if (l_use_array_cache) then
   call allocate_cached('ice', ice, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(ice(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ice',__FILE__,__LINE__)
    ice = 0.0

if (l_use_array_cache) then
   call allocate_cached('z0', z0, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(z0(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0',__FILE__,__LINE__)
    z0 = 0.0

if (l_use_array_cache) then
   call allocate_cached('root', root, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(root(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array root',__FILE__,__LINE__)
    root = 0.0

if (l_use_array_cache) then
   call allocate_cached('plcov_mx', plcov_mx, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(plcov_mx(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mx',__FILE__,__LINE__)
    plcov_mx = 0.0

if (l_use_array_cache) then
   call allocate_cached('plcov_mn', plcov_mn, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(plcov_mn(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array plcov_mn',__FILE__,__LINE__)
    plcov_mn = 0.0

if (l_use_array_cache) then
   call allocate_cached('lai_mx', lai_mx, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(lai_mx(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mx',__FILE__,__LINE__)
    lai_mx = 0.0

if (l_use_array_cache) then
   call allocate_cached('lai_mn', lai_mn, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(lai_mn(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lai_mn',__FILE__,__LINE__)
    lai_mn = 0.0

if (l_use_array_cache) then
   call allocate_cached('rs_min', rs_min, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(rs_min(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array rs_min',__FILE__,__LINE__)
    rs_min = 0.0

if (l_use_array_cache) then
   call allocate_cached('urban', urban, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(urban(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array urban',__FILE__,__LINE__)
    urban = 0.0

if (l_use_array_cache) then
   call allocate_cached('for_d', for_d, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(for_d(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_d',__FILE__,__LINE__)
    for_d = 0.0

if (l_use_array_cache) then
   call allocate_cached('for_e', for_e, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(for_e(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array for_e',__FILE__,__LINE__)
    for_e = 0.0

if (l_use_array_cache) then
   call allocate_cached('skinc', skinc, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(skinc(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array skinc',__FILE__,__LINE__)
    skinc = 0.0

if (l_use_array_cache) then
   call allocate_cached('emissivity', emissivity, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(emissivity(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emissivity',__FILE__,__LINE__)
    emissivity = 0.0

if (l_use_array_cache) then
   call allocate_cached('fr_ocean', fr_ocean, [tg%ie,tg%je,tg%ke,n_data])
else
   allocate(fr_ocean(tg%ie,tg%je,tg%ke,n_data), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_ocean',__FILE__,__LINE__)
    fr_ocean = 0.0

  END SUBROUTINE allocate_lu_ds_target_fields

END MODULE mo_lu_tg_fields
