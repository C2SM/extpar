!+ Fortran module with structures to describe various grid types
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_2         2011/03/25 Hermann Asensio
!  change definition of data type icon_grid_def
!  update doxygen documentation
! V1_7         2013/01/25 Guenther Zaengl 
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON 
! V2_0_3       2015-01-12 Juergen Helmert
!  Combined COSMO Release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with structures to describe various grid types
!!
!! data type for the description of
!!  - rotated longitude latitude grid 
!!  - regular longitude latitude grid
!!  - icosahedral triangular grid 
!> \author Hermann Asensio
MODULE mo_grid_structures

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp
USE mo_kind, ONLY: i4
USE mo_kind, ONLY: i4
USE mo_io_units, ONLY: filename_max


IMPLICIT NONE
PRIVATE

PUBLIC :: reg_lonlat_grid
PUBLIC :: rotated_lonlat_grid
PUBLIC :: icosahedral_triangular_grid
PUBLIC :: icon_grid_def
PUBLIC :: target_grid_def
PUBLIC :: igrid_icon
PUBLIC :: igrid_cosmo

!> Definition of data type with target grid definition 
TYPE :: target_grid_def
  INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO
  INTEGER (KIND=i4) :: ie !< number of grid elements in first dimension of target grid
  INTEGER (KIND=i4) :: je !< number of grid elements in second dimension of target grid
  INTEGER (KIND=i4) :: ke !< number of grid elements in third dimension of target grid
  REAL    (KIND=wp) :: minlon, maxlon, minlat, maxlat !< minimum and maximum longitudes and latitudes
  REAL    (KIND=wp) :: minlon_s, maxlon_s !< shifted minimum and maximum longitudes (range 0 - 360 deg)
  INTEGER (KIND=i4), ALLOCATABLE :: search_index(:,:) !< list for ICON search start index
END TYPE target_grid_def

INTEGER (KIND=i4), PARAMETER :: igrid_icon = 1 !< parameter to identify ICON grid
INTEGER (KIND=i4), PARAMETER :: igrid_cosmo = 2 !< parameter to identify COSMO grid

!> Definition of data type with ICON grid definition 
TYPE :: icosahedral_triangular_grid
  INTEGER (KIND=i4) :: nvertex_per_cell  !< number of vertices per cell
  INTEGER (KIND=i4) :: nchilds_per_cell  !< number of child cells per cell
  INTEGER (KIND=i4) :: ncells_per_edge   !< number of cells per edge
  INTEGER (KIND=i4) :: nedges_per_vertex !< number of edges per vertex
  INTEGER (KIND=i4) :: ncell             !< number of cells
  INTEGER (KIND=i4) :: nvertex           !< number of vertices
  INTEGER (KIND=i4) :: nedge             !< number of edges
  CHARACTER (LEN=filename_max) :: nc_grid_file !< filename of grid file which contains the Icon grid coordinates 
  INTEGER           :: number_of_grid_used     !< number of grid used, for grib2
  CHARACTER(LEN=36) :: uuidOfHGrid          ! ICON grid ID
  INTEGER           :: grid_root             ! root division ratio
  INTEGER           :: grid_level            ! grid level
END TYPE icosahedral_triangular_grid

!> Definition of data type with ICON grid definition for multiple model domains
TYPE :: icon_grid_def
  INTEGER (KIND=i4) :: nvertex_per_cell  !< number of vertices per cell
  INTEGER (KIND=i4) :: nchilds_per_cell  !< number of child cells per cell
  INTEGER (KIND=i4) :: ncells_per_edge   !< number of cells per edge
  INTEGER (KIND=i4) :: nedges_per_vertex !< number of edges per vertex
  INTEGER (KIND=i4) :: ncell             !< number of cells
  INTEGER (KIND=i4) :: nvertex           !< number of vertices
  INTEGER (KIND=i4) :: nedge             !< number of edges
  CHARACTER (LEN=filename_max) :: grid_file !< filename of Icon grid coordinates
END TYPE icon_grid_def

!> Definition of data type to describe a rotated lonlat grid
TYPE :: rotated_lonlat_grid
  REAL  (KIND=wp) :: pollon !< longitude of the rotated north pole (in degrees, E>0)
  REAL  (KIND=wp) :: pollat !< latitude of the rotated north pole (in degrees, N>0)
  REAL  (KIND=wp) :: polgam !< longitude (in the rotated system) of the geographical north pole
  REAL  (KIND=wp)  :: startlon_rot !< transformed longitude of the lower left grid point of the total domain (in degrees, E>0)
  REAL  (KIND=wp)  :: startlat_rot   !< transformed latitude of the lower left grid point of the total domain (in degrees, N>0)
  REAL  (KIND=wp)  :: dlon_rot !< grid point distance in zonal direction (in degrees, rotated system)
  REAL  (KIND=wp)  :: dlat_rot !< grid point distance in meridional direction (in degrees, rotated system)
  INTEGER  (KIND=i4) :: nlon_rot !< number of grid points in zonal direction (rotated system)
  INTEGER  (KIND=i4) :: nlat_rot !< number of grid points in meridional direction (rotated system)
  INTEGER  (KIND=i4) :: ke_tot !< number of grid points in vertical direction
END TYPE rotated_lonlat_grid

!> Definition of data type to describe a regular lonlat grid 
TYPE :: reg_lonlat_grid
       REAL (KIND=wp) :: start_lon_reg !< start longitude for regular lon-lat grid (usually western boundary)
       REAL (KIND=wp) :: end_lon_reg   !< end longitude for regular lon-lat grid (usually eastern boundary)

       REAL (KIND=wp) :: start_lat_reg !< start latitude for regular lon-lat grid 
!(usually southern boundary, but northern boundary possible of dlat_reg has negative increment (grid order from north to south)
       REAL (KIND=wp) :: end_lat_reg !< end latitude for regular lon-lat grid  
!(usually northern boundary, but southern boundary possible of dlat_reg has negative increment (grid order from north to south)


       REAL (KIND=wp) :: dlon_reg !< grid point distance in zonal direction (in degrees) for regular lon-lat grid
       REAL (KIND=wp) :: dlat_reg !< grid point distance in meridional direction (in degrees) 
!for regular lon-lat grid (negative increment for grid order from north to south)
       INTEGER (KIND=i4) :: nlon_reg !< number of grid elements in zonal direction for regular lon-lat grid
       INTEGER (KIND=i4) :: nlat_reg !< number of grid elements in meridional direction for regular lon-lat grid
END TYPE reg_lonlat_grid

END MODULE mo_grid_structures


