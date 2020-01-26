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
       &    allocate_topo_target_fields
            

  PUBLIC ::   add_parameters_domain, &
       &       vertex_param, &
       &       allocate_additional_hh_param


       
  REAL(KIND=wp), ALLOCATABLE  :: hh_topo(:,:,:), &      !< mean height
       &                         hh_topo_max(:,:,:), &  !< maximum height
       &                         hh_topo_min(:,:,:), &  !< minimum height 
       &                         stdh_topo(:,:,:), &    !< standard deviation of subgrid scale orographic height
       &                         theta_topo(:,:,:), & !< sso parameter, angle of principal axis
       &                         aniso_topo(:,:,:), & !< sso parameter, anisotropie factor
       &                         slope_topo(:,:,:), & !< sso parameter, mean slope
       &                         fr_land_topo(:,:,:), & !< fraction land due to GLOBE raw data
       &                         z0_topo(:,:,:), & !< roughness length due to orography
       &                         slope_asp_topo(:,:,:), &   !< lradtopo parameter, slope aspect
       &                         slope_ang_topo(:,:,:), &   !< lradtopo parameter, slope angle
       &                         horizon_topo  (:,:,:,:), & !< lradtopo parameter, horizon
       &                         skyview_topo  (:,:,:)   !< lradtopo parameter, skyview

  !> data structure for parameters on vertices of Icon grid
  TYPE add_parameters_domain
     REAL(KIND=wp), ALLOCATABLE     :: hh_vert(:,:,:)   !< height on vertex
     INTEGER (KIND=i4), ALLOCATABLE :: npixel_vert(:,:,:) !< number of raw data pixel corresponding to vertex
  END TYPE add_parameters_domain

  TYPE(add_parameters_domain) :: vertex_param  !< additional external parameters for ICON domain

  CONTAINS

  !> allocate fields for GLOBE target data 
  SUBROUTINE allocate_topo_target_fields(tg,nhori)
  
    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN)     :: nhori

    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_topo_target_fields')
      
    ALLOCATE (fr_land_topo(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array fr_land_topo',__FILE__,__LINE__)
    fr_land_topo = 0.0

    ALLOCATE (hh_topo(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hh_topo',__FILE__,__LINE__)
    hh_topo = 0.0

    ALLOCATE (hh_topo_max(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hh_topo_max',__FILE__,__LINE__)
    hh_topo_max = 0.0

    ALLOCATE (hh_topo_min(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hh_topo_min',__FILE__,__LINE__)
    hh_topo_min = 0.0
    
    ALLOCATE (stdh_topo(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array stdh_topo',__FILE__,__LINE__)
    stdh_topo = 0.0

    ALLOCATE (theta_topo(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array theta_topo',__FILE__,__LINE__)
    theta_topo = 0.0

    ALLOCATE (aniso_topo(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array aniso_topo',__FILE__,__LINE__)
    aniso_topo = 0.0

    ALLOCATE (slope_topo(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array slope_topo',__FILE__,__LINE__)
    slope_topo = 0.0

    
    ALLOCATE (z0_topo(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array z0_topo',__FILE__,__LINE__)
    z0_topo = 0.0

    ALLOCATE (slope_asp_topo(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array slope_asp_topo',__FILE__,__LINE__)
    slope_asp_topo = 0.0
    ALLOCATE (slope_ang_topo(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array slope_ang_topo',__FILE__,__LINE__)
    slope_ang_topo = 0.0
    ALLOCATE (horizon_topo(1:tg%ie,1:tg%je,1:tg%ke,nhori), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array horizon_topo',__FILE__,__LINE__)
    horizon_topo = 0.0
    ALLOCATE (skyview_topo(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array skyview_topo',__FILE__,__LINE__)
    skyview_topo = 0.0

    CALL logging%info('Exit routine: allocate_topo_target_fields')

  END SUBROUTINE allocate_topo_target_fields

  !> allocate additional parameters which correspond to the vertex
  !!
  !! the target grid has the dimension nvertex
  !! for future developments (optimizations, other code structure) the target grid is 
  !! defined as a 3-dimensional matrix, but the dimension are set to (nvertex,1,1) in this case
  SUBROUTINE allocate_additional_hh_param(nvertex)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvertex  !< number of vertices in target domains

    INTEGER, PARAMETER  :: je = 1, ke = 1

    INTEGER(KIND=i4)    :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_additional_hh_param')

    ALLOCATE(vertex_param%hh_vert(1:nvertex,1:je,1:ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the vertex_param%hh_vert(nvertex,je,ke',__FILE__,__LINE__)
    vertex_param%hh_vert = 0.0

    ALLOCATE(vertex_param%npixel_vert(1:nvertex,1:je,1:ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the vertex_param%npixel_vert(nvertex,je,ke)',__FILE__,__LINE__)
    vertex_param%npixel_vert = 0

    CALL logging%info('Exit routine: allocate_additional_hh_param')

  END SUBROUTINE allocate_additional_hh_param



END Module mo_topo_tg_fields

