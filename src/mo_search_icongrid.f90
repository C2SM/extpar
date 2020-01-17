!+ Fortran module to find grid element index in ICON grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_2         2011/03/25 Hermann Asensio
!  add point in polygon test for ICON grid search
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
!>
!! Fortran module to find grid element index in ICON grid
!!
!! @author Hermann Asensio, DWD
!!
!! @par Revision History
!! Initial realease by Hermann Asensio (2009-07-31)
!!
MODULE mo_search_icongrid

  USE mo_kind,                ONLY: wp, i4

  USE mo_additional_geometry, ONLY: gc2cc,                  &
       &                            scal_pro,               &
       &                            point_in_grid_element

  USE mo_base_geometry,       ONLY: geographical_coordinates
  USE mo_base_geometry,       ONLY: cartesian_coordinates

  USE mo_grid_structures,     ONLY: icon_grid_def
  USE mo_icon_domain,         ONLY: icon_domain

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id: mo_search_icongrid.f90,v 1.11 2013-04-16 11:09:20 for0adm Exp $'

  PUBLIC :: walk_to_nc
  PUBLIC :: find_nc
  PUBLIC :: find_nearest_vert

CONTAINS

  !> Find the nearest grid cell index in the ICON grid for given (geographical) target coordinates. 
  !!
  !! Search for nearest grid point for given (geographical) target coordinates.
  !! Start at the highest level (level = 0) and find the nearest grid element at this level.
  !! Proceed to the next level and check there for the four child cell which is the nearest grid cell.
  !! Repeat the last step until the first domain level (start_lev of the Namelist input for the grid generator) is reached.
  !! If further refinement domains are given, 
  !! search also in the refinement domains and put result in array nearest_cell_ids(ndom).
  SUBROUTINE find_nc( target_geo_co,               &
       &              nvertex_per_cell,            &
       &              icon_dom_def,                &
       &              icon_grid_region,            &
       &              start_cell_id,               &
       &              nearest_cell_id )

    TYPE(geographical_coordinates), INTENT(IN)  :: target_geo_co    
    !< target coordinates in geographical system of point for which the nearest ICON grid cell is to be determined
    INTEGER, INTENT(IN)                         :: nvertex_per_cell   !< number of vertices per cell
    TYPE(icon_grid_def), INTENT(IN)             :: icon_dom_def       !< structure which contains the definition of the ICON grid
    TYPE(icon_domain), INTENT(IN), TARGET       :: icon_grid_region   !< Data structure with ICON domain with refinement domain,
                                                                      !< dimension (1:ndom)
    INTEGER (KIND=i4), INTENT(INOUT)            :: start_cell_id      !< id of starting cell

    INTEGER (KIND=i4), INTENT(OUT)              :: nearest_cell_id    !< id of nearest cell

    ! local variables

    TYPE(cartesian_coordinates)  :: target_cc_co     
    !< coordinates in cartesian system of point for which the nearest ICON grid cell is to be determined

    INTEGER :: undefined

    !!HA debug
    ! print *,'Entering subroutine find_nc'
    undefined = 0
    nearest_cell_id = undefined ! set the nearest_cell_id to "undefined" as default
    ! transform geographical coordinates to cartesian coordinates of target point
    target_cc_co = gc2cc(target_geo_co)

    CALL walk_to_nc( icon_grid_region,               &
         &           target_cc_co,                   &
         &           start_cell_id,                  &
         &           nvertex_per_cell,               &
         &           icon_dom_def%nedges_per_vertex, &
         &           nearest_cell_id )

  END SUBROUTINE find_nc

  !-----------------------------------------------------------------------------

  !> Go to the nearest grid cell in the ICON grid 
  !!
  !! Search for nearest grid point for given (geographical) target coordinates.
  !! Walk from a (given) starting cell in the direction of the target coordintes
  !! until the neighbour cells have a larger distance to the current cell.
  !! This algorithm works on global domains, for regional domains the search
  !! might get stuck at the boundaries of the domains.
  SUBROUTINE walk_to_nc ( grid,              &
       &                  target_cc_co,      &
       &                  start_cell_id,     &
       &                  nvertex_per_cell,  &
       &                  ncells_per_vertex, &
       &                  nearest_cell_id )

    TYPE(icon_domain), INTENT(IN)               :: grid              !> Data structure with ICON grid
    TYPE(cartesian_coordinates), INTENT(IN)     :: target_cc_co     
    !>  target coordinates in cartesian system of point for which the nearest ICON grid cell is to be determined
    INTEGER (KIND=i4), INTENT(INOUT)            :: start_cell_id     !> id of starting cell
    INTEGER, INTENT(IN)                         :: nvertex_per_cell  !< number of vertices per cell
    INTEGER, INTENT(IN)                         :: ncells_per_vertex !< number of cells per vertex
    INTEGER (KIND=i4), INTENT(OUT)              :: nearest_cell_id   !> id of nearest cell

    ! local variables
    TYPE(cartesian_coordinates)  :: cell_cc          !> coordinates of cell centre in cartesian system 
    TYPE(cartesian_coordinates)  :: neighbour_cc     !> coordinates of a neighbour cell centre in cartesian system
    INTEGER                      :: nb_cell_id       !> neighbour cell id

    INTEGER   (KIND=i4)          :: current_cell_id
    INTEGER   (KIND=i4)          :: next_cell_id

    REAL(KIND=wp)                :: sp               
    !> cos arc length of  of geodesic arc with endpoints x0,x1 (normalized scalar product of the two points)
    REAL(KIND=wp)                :: sp_max           !> maximum of the scalar product of two points (minimal distance)

    LOGICAL                      :: searching

    INTEGER                      :: nj               !< counter

    !PRINT *,'entering walk_to_nc'
    searching = .TRUE.   ! set searching to "true"

    nearest_cell_id = start_cell_id ! initial setting
    current_cell_id = start_cell_id ! initial setting
    next_cell_id    = start_cell_id ! initial setting 

    ! cartesian coordinates of start cell centre
    cell_cc =  grid%cells%cc_center(start_cell_id)

    ! calculate a measure for the distance to target point
    sp = scal_pro(target_cc_co, cell_cc)
    sp_max = sp

    DO WHILE(searching)
      searching = .false. ! abort condition
      DO nj=1, nvertex_per_cell
        nb_cell_id = grid%cells%neighbor_index(current_cell_id,nj) ! get cell id of neighbour cells
        IF (nb_cell_id > 0 ) THEN                                  ! 0 is the "undefined" value for the cell id
          neighbour_cc = grid%cells%cc_center(nb_cell_id)          ! get cartesian coordinates of neighbour cell
          sp = scal_pro(target_cc_co,neighbour_cc)                 ! calculate measure for distance to target point
          IF (sp > sp_max) THEN                                    ! if neighbour cell is nearer to target point than the "old" cell
            sp_max = sp                                            ! save new distance measure
            next_cell_id = nb_cell_id                              ! save cell id
            searching = .true.                                     ! continue with search loop
          ENDIF
        ENDIF
      ENDDO
      current_cell_id = next_cell_id       ! move one cell toward target point
    ENDDO
    nearest_cell_id = current_cell_id      ! set nearest_cell_id to the cell id
                                           ! which has smallest distance (i.e. largest sp) to target point

    ! check with a point in polygon test
    CALL control_nc( grid,              &
         &           target_cc_co,      &
         &           ncells_per_vertex, &
         &           nvertex_per_cell,  &
         &           nearest_cell_id )

    ! save the start cell ID for next search point
    ! if it is zero, it will be reinitialized from the search index list in the calling program
    start_cell_id = nearest_cell_id

  END SUBROUTINE walk_to_nc

  !-----------------------------------------------------------------------------

  !> Find the nearest vertex for a given grid cell in the ICON grid 
  !!
  !! Check the distance of the target point to the vertices
  !! give out the id of the nearest vertex
  SUBROUTINE find_nearest_vert( grid,             &
       &                        target_cc_co,     &
       &                        cell_id,          &
       &                        nvertex_per_cell, &
       &                        nearest_vert_id )

    TYPE(icon_domain), INTENT(IN)               :: grid             !> Data structure with ICON grid
    TYPE(cartesian_coordinates), INTENT(IN)     :: target_cc_co     
    !>  target coordinates in cartesian system of point for which the nearest ICON grid cell is to be determined
    INTEGER (KIND=i4), INTENT(IN)               :: cell_id          !> id of cell
    INTEGER, INTENT(IN)                         :: nvertex_per_cell !< number of vertices per cell
    INTEGER (KIND=i4), INTENT(OUT)              :: nearest_vert_id  !> id of nearest cell

    ! local variables
    INTEGER (KIND=i4)           :: vert_id_vec(1:nvertex_per_cell) !< indices of cells vertices
    REAL(KIND=wp)               :: sp_vec(1:nvertex_per_cell)      !< cos arc length of  of geodesic arc 
    !! is used as a distance measure  
    INTEGER :: max_index

    vert_id_vec(:) = grid%cells%vertex_index(cell_id,:) ! get the indices of the cells vertices

    ! calculate a measure for the distance to target point,
    ! the sp value is between [-1,1], with 1 for identical points and -1 for antipodes (on the sphere!)
    sp_vec(:) = scal_pro(target_cc_co, grid%verts%cc_vertex(vert_id_vec)) 

    max_index = MAXLOC(sp_vec,DIM=1)  ! the 
    nearest_vert_id = vert_id_vec(max_index)


  END SUBROUTINE find_nearest_vert

  !-----------------------------------------------------------------------------

  !> control with a point in polygon test, if nearest cell contains search point
  !! if not, find the correct the cell id
  SUBROUTINE control_nc( grid,              &
       &                 target_cc_co,      &
       &                 ncells_per_vertex, &
       &                 nvertex_per_cell,  &
       &                 nearest_cell_id )

    TYPE(icon_domain), INTENT(IN)               :: grid              !< Data structure with ICON grid
    TYPE(cartesian_coordinates), INTENT(IN)     :: target_cc_co      
    !<  target coordinates in cartesian system of point for which the nearest ICON grid cell is to be determined
    INTEGER, INTENT(IN)                         :: ncells_per_vertex !< number of cells per vertex
    INTEGER, INTENT(IN)                         :: nvertex_per_cell  !< number of vertices per cell
    INTEGER (KIND=i4), INTENT(INOUT)            :: nearest_cell_id   !> id of nearest cell

    ! local variables
    TYPE(cartesian_coordinates)  :: cc_vertices(1:nvertex_per_cell) 
    ! cartesian coordinates of vertices of grid element for point in polygon test
    INTEGER  (KIND=i4)           :: ivert            !< counter
    INTEGER  (KIND=i4)           :: vert_index       !< index
    INTEGER                      :: inflag
    INTEGER (KIND=i4)            :: n_vert_id        !< vertex id
    INTEGER (KIND=i4)            :: nb_cell_id       !< cell id
    INTEGER (KIND=i4)            :: nev              !< counter

    DO ivert=1,nvertex_per_cell
      vert_index = grid%cells%vertex_index(nearest_cell_id,ivert)
      cc_vertices(ivert) = grid%verts%cc_vertex(vert_index)
    ENDDO

    CALL point_in_grid_element(target_cc_co,nvertex_per_cell,cc_vertices,inflag)
    !PRINT *,'--HA debug SUBROUTINE control_nc --'
    !PRINT *,'inflag: ',inflag
    !PRINT *,'--HA debug SUBROUTINE control_nc --'

    IF (inflag > 0) THEN ! point in grid element
      RETURN
    ELSE
      CALL find_nearest_vert( grid,             &
           &                  target_cc_co,     &
           &                  nearest_cell_id,  &
           &                  nvertex_per_cell, &
           &                  n_vert_id)

      DO nev=1,ncells_per_vertex
        nb_cell_id = grid%verts%cell_index(n_vert_id,nev) ! neighbour cell index
        IF (nb_cell_id > 0) THEN
          DO ivert=1,nvertex_per_cell
            vert_index = grid%cells%vertex_index(nb_cell_id,ivert)
            cc_vertices(ivert) = grid%verts%cc_vertex(vert_index)
          ENDDO
          CALL point_in_grid_element(target_cc_co,nvertex_per_cell,cc_vertices,inflag)
          IF (inflag > 0) THEN ! point in grid element
            nearest_cell_id = nb_cell_id
            RETURN
          ENDIF
        ENDIF
      ENDDO

      ! If we arrive here, none of the points is within the cell
      nearest_cell_id = 0

    ENDIF
  END SUBROUTINE control_nc

END MODULE mo_search_icongrid
