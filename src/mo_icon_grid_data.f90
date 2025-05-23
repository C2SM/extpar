!+ Fortran module for icon grid data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_2         2011/03/25 Hermann Asensio 
!  update to support ICON refinement grids
! V1_3         2011/04/19 Hermann Asensio
!  introduce variables clon, clat, clon_vertices, clat_vertices
! V1_7         2013/01/25 Guenther Zaengl 
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON  
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for icon grid data 
!> \author Hermann Asensio
MODULE mo_icon_grid_data

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_icon_domain,           ONLY: icon_domain

  USE mo_grid_structures,       ONLY: icosahedral_triangular_grid, &
       &                              icon_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: icon_grid, &
       &    icon_grid_region, &
       &    nvertex_per_cell, &
       &    icon_dom_def, &
       &    clon, clat, &
       &    clon_vertices, clat_vertices, &
       &    allocate_icon_coor

  TYPE(icosahedral_triangular_grid) :: icon_grid          !< structure which contains the definition of the ICON grid
  TYPE(icon_grid_def)               :: icon_dom_def       !< structure with the definition of the various ICON domains
  TYPE(icon_domain)                 :: icon_grid_region

  INTEGER (KIND=i4)                 :: nvertex_per_cell   !< number of verices per Icon grid cell

  REAL (KIND=wp), ALLOCATABLE       :: clon(:), &            !< longitude of Icon grid cells, (1:ncell) in radians
       &                               clat(:), &            !< latitude of Icon grid cells, (1:ncell) in radians
       &                               clon_vertices(:,:), & !< longitude of Icon grid vertices, (1:ncell,1:nvertex_per_cell) in radians
       &                               clat_vertices(:,:)    !< latitude of Icon grid vertices, (1:ncell,1:nvertex_per_cell) in radians

CONTAINS

  !> allocate icon coordinate arrays
  SUBROUTINE allocate_icon_coor(ncell, nvertex_per_cell)

    INTEGER (KIND=i4), INTENT(IN)  :: ncell, &            !< number of cells
                                      nvertex_per_cell !< number of vertices per cell
    ! local variables
    INTEGER(KIND=i4)               :: istat, ist

    ist = 0
    istat = 0

    ALLOCATE(clon(1:ncell),STAT=istat)
    ist = ist+istat

    ALLOCATE(clat(1:ncell),STAT=istat)
    ist = ist+istat

    ALLOCATE(clon_vertices(1:ncell,1:nvertex_per_cell),STAT=istat)
    ist = ist+istat

    ALLOCATE(clat_vertices(1:ncell,1:nvertex_per_cell),STAT=istat)
    ist = ist+istat

    IF (istat /= 0) THEN
      CALL logging%error('allocate_icon_coor',__FILE__,__LINE__)
    ENDIF

  END SUBROUTINE allocate_icon_coor

END MODULE mo_icon_grid_data

