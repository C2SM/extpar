!+ Fortran Module to search index in target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_2         2011/03/25 Hermann Asensio
!  update to support ICON refinement grids
! V1_7         2013/01/25 Guenther Zaengl 
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON  
! V1_14        2014-07-18 Juergen Helmert
!  Combined COSMO Release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran Module to search index in target grid
!! Fortran module to find nearest grid element index in target grid for external parameters
!! 
!! Depending of the type of the grid, the nearest grid point search for the ICON grid or the 
!! Cosmo-grid are called
!!
MODULE mo_search_target_grid


  USE mo_kind,            ONLY: wp, i4, i8
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: find_nearest_target_grid_element


  CONTAINS

  SUBROUTINE find_nearest_target_grid_element(  point_lon_geo, &
    &                                           point_lat_geo, &
    &                                           tg,            &
    &                                           start_cell_id, &
    &                                           tg_el_ie,      &
    &                                           tg_el_je,      &
    &                                           tg_el_ke)

  USE mo_grid_structures, ONLY: target_grid_def
  USE mo_grid_structures, ONLY: reg_lonlat_grid
  USE mo_grid_structures, ONLY: target_grid_def

  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo

  USE mo_search_ll_grid, ONLY: find_reg_lonlat_grid_element_index
  USE mo_search_ll_grid, ONLY: find_rotated_lonlat_grid_element_index

  ! USE structure which contains the definition of the COSMO grid
  USE  mo_cosmo_grid, ONLY: COSMO_grid !< structure which contains the definition of the COSMO grid

  ! USE structure which contains the definition of the ICON grid
  USE  mo_icon_grid_data, ONLY: ICON_grid !< structure which contains the definition of the ICON grid

  ! USE icon domain structure wich contains the ICON coordinates (and parent-child indices etc)
   USE mo_icon_grid_data, ONLY:  icon_grid_region
   USE mo_icon_grid_data, ONLY: nvertex_per_cell
   USE mo_icon_grid_data, ONLY: icon_dom_def

  ! USE modules to search in ICON grid
  USE mo_search_icongrid, ONLY: find_nc

  USE mo_icon_domain,     ONLY: icon_domain

  USE mo_base_geometry,   ONLY: geographical_coordinates
  USE mo_base_geometry,   ONLY: cartesian_coordinates

  USE mo_additional_geometry,   ONLY: gc2cc

  USE mo_math_constants, ONLY: deg2rad, pi

  
  REAL (KIND=wp), INTENT(in) :: point_lon_geo       !< longitude coordinate in geographical system of input point 
  REAL (KIND=wp), INTENT(in) :: point_lat_geo       !< latitude coordinate in geographical system of input point
  TYPE(target_grid_def), INTENT(IN) :: tg           !< structure with target grid description

  INTEGER (KIND=i8), INTENT(INOUT) :: start_cell_id !< ID of starting cell
  INTEGER (KIND=i8), INTENT(OUT) :: tg_el_ie        !< Index tg_el_ie of target grid element nearest to the input point
  INTEGER (KIND=i8), INTENT(OUT) :: tg_el_je        !< Index tg_el_je of target grid element nearest to the input point
  INTEGER (KIND=i8), INTENT(OUT) :: tg_el_ke        !< Index tg_el_ke of target grid element nearest to the input point

  ! local variables
   TYPE(geographical_coordinates) :: target_geo_co  !< structure for geographical coordinates of raw data pixel
   TYPE(cartesian_coordinates)  :: target_cc_co     
!< coordinates in cartesian system of point for which the nearest ICON grid cell is to be determined

   ! variables for GME search
   INTEGER :: nip1 ! grid mesh dimension 
   REAL (KIND=wp)  :: zx,zy,zz ! cartesian coordinates of point
   REAL (KIND=wp), SAVE  :: spd_t = 1. ! threshold value for scalar product 
   INTEGER, SAVE :: kd = 1  ! diamond containing point
   INTEGER, SAVE :: kj1 = 0 ! nodal indices of nearest grid point
   INTEGER, SAVE :: kj2 = 1 ! on entry, kj1 and kj2 are first guess values
   REAL (KIND=wp), SAVE  :: sp =1.! scalar product between point and nearest GME nodal point
   LOGICAL :: ldebug=.FALSE.

       tg_el_ie = 0_i8 ! default settings
       tg_el_je = 0_i8
       tg_el_ke = 0_i8

       SELECT CASE(tg%igrid_type)
       CASE(igrid_icon)  ! ICON GRID

         target_geo_co%lon = point_lon_geo * deg2rad ! note that the ICON coordinates do not have the unit degree but radians
         target_geo_co%lat = point_lat_geo * deg2rad

         !target_cc_co = gc2cc(target_geo_co) ! transform the geographical coordinates of the point to cartesian coordinates
         CALL find_nc(target_geo_co,    &
           &          nvertex_per_cell, &
           &          icon_dom_def,     &
           &          icon_grid_region, &
           &          start_cell_id,    &
           &          tg_el_ie)

         tg_el_je = 1_i8
         tg_el_ke = 1_i8

       CASE(igrid_cosmo)  ! COSMO GRID

         CALL find_rotated_lonlat_grid_element_index(point_lon_geo, &
           &                                    point_lat_geo,     &
           &                                    COSMO_grid,        &
           &                                    tg_el_ie,          &
           &                                    tg_el_je)
         tg_el_ke = 1_i8

       END SELECT

  END  SUBROUTINE find_nearest_target_grid_element


END MODULE mo_search_target_grid



