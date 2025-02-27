!+ Fortran module for GLOBE data on target grid for external parameters
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_2         2011/03/25 Hermann Asensio
!  Update doxygen documetation (comments)
! V2_0         2013/06/04 Anne Roches
!  introduction of the topographical corrected radiation parameters
! V2_0         2013/06/04 Martina Messmer
!  renaming of all the variables that contained a 'globe' into 'topo'
! V1_14        2014-07-18 Juergen Helmert
!  Combined COSMO Release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for GLOBE data on target grid for external parameters
!> \author Hermann Asensio
MODULE mo_topo_tg_fields

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4
  USE mo_array_cache,           ONLY: allocate_cached
  USE mo_grid_structures,       ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fr_land_topo,       &
       &    hh_topo,            &
       &    hh_topo_max,        &
       &    hh_topo_min,        &
       &    stdh_topo,          &
       &    theta_topo,         &
       &    aniso_topo,         &
       &    slope_topo,         &
       &    z0_topo,            &
       &    slope_asp_topo,     &
       &    slope_ang_topo,     &
       &    horizon_topo,       &
       &    skyview_topo,       &
       &    sgsl,               &
       &    allocate_topo_target_fields, &
       &    fr_land_topo_globe,       &
       &    hh_topo_globe,            &
       &    hh_topo_max_globe,        &
       &    hh_topo_min_globe,        &
       &    stdh_topo_globe,          &
       &    theta_topo_globe,         &
       &    aniso_topo_globe,         &
       &    slope_topo_globe,         &
       &    z0_topo_globe,            &
       &    slope_asp_topo_globe,     &
       &    slope_ang_topo_globe,     &
       &    horizon_topo_globe,       &
       &    skyview_topo_globe,       &
       &    sgsl_globe,               &
       &    allocate_topo_target_fields_globe

 

  REAL(KIND=wp), POINTER  :: hh_topo(:,:,:), &      !< mean height
       &                     hh_topo_max(:,:,:), &  !< maximum height
       &                     hh_topo_min(:,:,:), &  !< minimum height
       &                     stdh_topo(:,:,:), &    !< standard deviation of subgrid scale orographic height
       &                     theta_topo(:,:,:), & !< sso parameter, angle of principal axis
       &                     aniso_topo(:,:,:), & !< sso parameter, anisotropie factor
       &                     slope_topo(:,:,:), & !< sso parameter, mean slope
       &                     fr_land_topo(:,:,:), & !< fraction land due to GLOBE raw data
       &                     z0_topo(:,:,:), & !< roughness length due to orography
       &                     slope_asp_topo(:,:,:), &   !< lradtopo parameter, slope aspect
       &                     slope_ang_topo(:,:,:), &   !< lradtopo parameter, slope angle
       &                     horizon_topo  (:,:,:,:), & !< lradtopo parameter, horizon
       &                     skyview_topo  (:,:,:), &   !< lradtopo parameter, skyview
       &                     sgsl(:,:,:), &   !< subgrid-scale slopes
       &                     hh_topo_globe(:,:,:), &      !< mean height
       &                     hh_topo_max_globe(:,:,:), &  !< maximum height
       &                     hh_topo_min_globe(:,:,:), &  !< minimum height
       &                     stdh_topo_globe(:,:,:), &    !< standard deviation of subgrid scale orographic height
       &                     theta_topo_globe(:,:,:), & !< sso parameter, angle of principal axis
       &                     aniso_topo_globe(:,:,:), & !< sso parameter, anisotropie factor
       &                     slope_topo_globe(:,:,:), & !< sso parameter, mean slope
       &                     fr_land_topo_globe(:,:,:), & !< fraction land due to GLOBE raw data
       &                     z0_topo_globe(:,:,:), & !< roughness length due to orography
       &                     slope_asp_topo_globe(:,:,:), &   !< lradtopo parameter, slope aspect
       &                     slope_ang_topo_globe(:,:,:), &   !< lradtopo parameter, slope angle
       &                     horizon_topo_globe  (:,:,:,:), & !< lradtopo parameter, horizon
       &                     skyview_topo_globe  (:,:,:), &   !< lradtopo parameter, skyview
       &                     sgsl_globe(:,:,:) !< subgrid-scale slopes

  CONTAINS

  !> allocate fields for GLOBE target data
  SUBROUTINE allocate_topo_target_fields(tg,nhori, lcompute_sgsl, l_use_array_cache)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN)     :: nhori
    LOGICAL, INTENT(IN)               :: lcompute_sgsl
    LOGICAL, INTENT(in)               :: l_use_array_cache
    
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0 
    
    CALL logging%info('Enter routine: allocate_topo_target_fields')

if (l_use_array_cache) then
   call allocate_cached('fr_land_topo', fr_land_topo, [tg%ie,tg%je,tg%ke])
else
   allocate(fr_land_topo(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_topo',__FILE__,__LINE__)
    fr_land_topo = 0.0

if (l_use_array_cache) then
   call allocate_cached('hh_topo', hh_topo, [tg%ie,tg%je,tg%ke])
else
   allocate(hh_topo(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hh_topo',__FILE__,__LINE__)
    hh_topo = 0.0

if (l_use_array_cache) then
   call allocate_cached('hh_topo_max', hh_topo_max, [tg%ie,tg%je,tg%ke])
else
   allocate(hh_topo_max(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hh_topo_max',__FILE__,__LINE__)
    hh_topo_max = 0.0

if (l_use_array_cache) then
   call allocate_cached('hh_topo_min', hh_topo_min, [tg%ie,tg%je,tg%ke])
else
   allocate(hh_topo_min(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hh_topo_min',__FILE__,__LINE__)
    hh_topo_min = 0.0

if (l_use_array_cache) then
   call allocate_cached('stdh_topo', stdh_topo, [tg%ie,tg%je,tg%ke])
else
   allocate(stdh_topo(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array stdh_topo',__FILE__,__LINE__)
    stdh_topo = 0.0

if (l_use_array_cache) then
   call allocate_cached('theta_topo', theta_topo, [tg%ie,tg%je,tg%ke])
else
   allocate(theta_topo(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array theta_topo',__FILE__,__LINE__)
    theta_topo = 0.0

if (l_use_array_cache) then
   call allocate_cached('aniso_topo', aniso_topo, [tg%ie,tg%je,tg%ke])
else
   allocate(aniso_topo(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array aniso_topo',__FILE__,__LINE__)
    aniso_topo = 0.0

if (l_use_array_cache) then
   call allocate_cached('slope_topo', slope_topo, [tg%ie,tg%je,tg%ke])
else
   allocate(slope_topo(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array slope_topo',__FILE__,__LINE__)
    slope_topo = 0.0


if (l_use_array_cache) then
   call allocate_cached('z0_topo', z0_topo, [tg%ie,tg%je,tg%ke])
else
   allocate(z0_topo(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0_topo',__FILE__,__LINE__)
    z0_topo = 0.0

if (l_use_array_cache) then
   call allocate_cached('slope_asp_topo', slope_asp_topo, [tg%ie,tg%je,tg%ke])
else
   allocate(slope_asp_topo(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array slope_asp_topo',__FILE__,__LINE__)
    slope_asp_topo = 0.0

if (l_use_array_cache) then
   call allocate_cached('slope_ang_topo', slope_ang_topo, [tg%ie,tg%je,tg%ke])
else
   allocate(slope_ang_topo(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array slope_ang_topo',__FILE__,__LINE__)
    slope_ang_topo = 0.0

if (l_use_array_cache) then
   call allocate_cached('horizon_topo', horizon_topo, [tg%ie,tg%je,tg%ke,nhori])
else
   allocate(horizon_topo(tg%ie,tg%je,tg%ke,nhori), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array horizon_topo',__FILE__,__LINE__)
    horizon_topo = 0.0

if (l_use_array_cache) then
   call allocate_cached('skyview_topo', skyview_topo, [tg%ie,tg%je,tg%ke])
else
   allocate(skyview_topo(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array skyview_topo',__FILE__,__LINE__)
    skyview_topo = 0.0

    IF (lcompute_sgsl) THEN
if (l_use_array_cache) then
   call allocate_cached('sgsl', sgsl, [tg%ie,tg%je,tg%ke])
else
   allocate(sgsl(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array sgsl',__FILE__,__LINE__)
      sgsl = 0.0
    ENDIF

    CALL logging%info('Exit routine: allocate_topo_target_fields')

  END SUBROUTINE allocate_topo_target_fields

  !> allocate fields for GLOBE target data
  SUBROUTINE allocate_topo_target_fields_globe(tg,nhori, lcompute_sgsl, l_use_array_cache)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN)     :: nhori
    LOGICAL, INTENT(IN)               :: lcompute_sgsl
    LOGICAL, INTENT(in)               :: l_use_array_cache
    
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0 
    
    CALL logging%info('Enter routine: allocate_topo_target_fields_globe')

if (l_use_array_cache) then
   call allocate_cached('fr_land_topo_globe', fr_land_topo_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(fr_land_topo_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_topo_globe',__FILE__,__LINE__)
    fr_land_topo_globe = 0.0

if (l_use_array_cache) then
   call allocate_cached('hh_topo_globe', hh_topo_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(hh_topo_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hh_topo_globe',__FILE__,__LINE__)
    hh_topo_globe = 0.0

if (l_use_array_cache) then
   call allocate_cached('hh_topo_max_globe', hh_topo_max_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(hh_topo_max_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hh_topo_max_globe',__FILE__,__LINE__)
    hh_topo_max_globe = 0.0

if (l_use_array_cache) then
   call allocate_cached('hh_topo_min_globe', hh_topo_min_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(hh_topo_min_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hh_topo_min_globe',__FILE__,__LINE__)
    hh_topo_min_globe = 0.0

if (l_use_array_cache) then
   call allocate_cached('stdh_topo_globe', stdh_topo_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(stdh_topo_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array stdh_topo_globe',__FILE__,__LINE__)
    stdh_topo_globe = 0.0

if (l_use_array_cache) then
   call allocate_cached('theta_topo_globe', theta_topo_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(theta_topo_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array theta_topo_globe',__FILE__,__LINE__)
    theta_topo_globe = 0.0

if (l_use_array_cache) then
   call allocate_cached('aniso_topo_globe', aniso_topo_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(aniso_topo_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array aniso_topo_globe',__FILE__,__LINE__)
    aniso_topo_globe= 0.0

if (l_use_array_cache) then
   call allocate_cached('slope_topo_globe', slope_topo_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(slope_topo_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array slope_topo_globe',__FILE__,__LINE__)
    slope_topo_globe = 0.0


if (l_use_array_cache) then
   call allocate_cached('z0_topo_globe', z0_topo_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(z0_topo_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0_topo_globe',__FILE__,__LINE__)
    z0_topo_globe = 0.0

if (l_use_array_cache) then
   call allocate_cached('slope_asp_topo_globe', slope_asp_topo_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(slope_asp_topo_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array slope_asp_topo_globe',__FILE__,__LINE__)
    slope_asp_topo_globe = 0.0

if (l_use_array_cache) then
   call allocate_cached('slope_ang_topo_globe', slope_ang_topo_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(slope_ang_topo_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array slope_ang_topo_globe',__FILE__,__LINE__)
    slope_ang_topo_globe = 0.0

if (l_use_array_cache) then
   call allocate_cached('horizon_topo_globe', horizon_topo_globe, [tg%ie,tg%je,tg%ke,nhori])
else
   allocate(horizon_topo_globe(tg%ie,tg%je,tg%ke,nhori), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array horizon_topo_globe',__FILE__,__LINE__)
    horizon_topo_globe = 0.0

if (l_use_array_cache) then
   call allocate_cached('skyview_topo_globe', skyview_topo_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(skyview_topo_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array skyview_topo_globe',__FILE__,__LINE__)
    skyview_topo_globe = 0.0

    IF (lcompute_sgsl) THEN
if (l_use_array_cache) then
   call allocate_cached('sgsl_globe', sgsl_globe, [tg%ie,tg%je,tg%ke])
else
   allocate(sgsl_globe(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array sgsl_globe',__FILE__,__LINE__)
      sgsl_globe = 0.0
    ENDIF

    CALL logging%info('Exit routine: allocate_topo_target_fields_globe')

  END SUBROUTINE allocate_topo_target_fields_globe
  

END MODULE mo_topo_tg_fields
