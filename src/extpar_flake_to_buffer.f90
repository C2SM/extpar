!+ Fortran main program to aggregate lake depth data to a target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio 
!  small bug fixes accroding to Fortran compiler warnings         
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
!> Fortran main program to aggregate lake depth data to a target grid
!!
!! @par extpar_flake_to_buffer 
!!
!! 
!! @author
!!     Hermann Asensio
!!     (DWD)
!!
!!
PROGRAM extpar_flake_to_buffer

  !Load the library information data:
  USE info_extpar, ONLY: info_define, info_print

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  USE mo_grid_structures, ONLY: target_grid_def,   &
    &                            reg_lonlat_grid,   &
    &                            rotated_lonlat_grid
  
  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo
  USE mo_grid_structures, ONLY: igrid_gme

  USE mo_target_grid_data, ONLY: lon_geo, &
    &                            lat_geo, &
    &                            no_raw_data_pixel, &
    &                            allocate_com_target_fields
  
  USE mo_target_grid_data, ONLY: tg
  
  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_icon_grid_data, ONLY: ICON_grid  !< structure which contains the definition of the ICON grid
 
  USE  mo_cosmo_grid, ONLY: COSMO_grid, &
    &                       lon_rot, &
    &                       lat_rot, &
    &                       allocate_cosmo_rc, &
    &                       get_cosmo_grid_info, &
    &                       calculate_cosmo_domain_coordinates

  USE mo_base_geometry,    ONLY:  geographical_coordinates, &
    &                             cartesian_coordinates

  USE mo_icon_domain,          ONLY: icon_domain, &
    &                             grid_cells,               &
    &                             grid_vertices,            &
    &                             construct_icon_domain,    &
    &                             destruct_icon_domain

  USE mo_io_units,          ONLY: filename_max

  USE mo_exception,         ONLY: message_text, message, finish

  USE mo_utilities_extpar, ONLY: abort_extpar
  
  USE mo_additional_geometry,   ONLY: cc2gc,                  &
    &                            gc2cc,                  &
    &                            arc_length,             &
    &                            cos_arc_length,         &
    &                            inter_section,          &
    &                            vector_product,         &
    &                            point_in_polygon_sp


  USE mo_math_constants,  ONLY: pi, pi_2, dbl_eps,rad2deg

  USE mo_flake_routines, ONLY: read_namelists_extpar_flake

  USE mo_flake_routines, ONLY:  get_dimension_flake_data, &
    &                             get_lonlat_flake_data
  

  USE mo_flake_data, ONLY: flake_grid, &
 &         lon_flake,  &
 &         lat_flake,  &
 &         allocate_raw_flake_fields,  &
 &         deallocate_raw_flake_fields

  USE mo_flake_tg_fields, ONLY: fr_lake, &
  &       lake_depth,    &
  &       flake_tot_npixel, &
  &       allocate_flake_target_fields

  USE mo_agg_flake, ONLY : agg_flake_data_to_target_grid


  USE mo_flake_output_nc, ONLY: write_netcdf_buffer_flake, &
    &                             write_netcdf_cosmo_grid_flake, &
    &                             write_netcdf_icon_grid_flake

  
  IMPLICIT NONE
  
  CHARACTER(len=filename_max) :: filename
  CHARACTER(len=filename_max) :: netcdf_filename

  CHARACTER(len=filename_max) :: input_namelist_file
  CHARACTER(len=filename_max) :: input_namelist_cosmo_grid !< file with input namelist with COSMO grid definition

  CHARACTER(len=filename_max) :: namelist_grid_def

  CHARACTER (len=filename_max) :: namelist_topo_data_input !< file with input namelist with GLOBE data information
  CHARACTER(len=filename_max) :: input_flake_namelist_file 
  CHARACTER(len=filename_max) :: flake_file

  CHARACTER (len=filename_max) :: raw_data_flake_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_flake_filename !< filename flake raw data

  CHARACTER (len=filename_max) :: flake_buffer_file !< name for flake buffer file
  CHARACTER (len=filename_max) :: flake_output_file !< name for flake output file


  INTEGER :: i, ip, ic, in


  INTEGER                      :: i_nc       !< number of cells
  INTEGER                      :: i_ne       !< number of edges
  INTEGER                      :: i_nv       !< number of vertices
  INTEGER                      :: nc_p_e     !< number of cells per edge
  INTEGER                      :: nv_p_c     !< number of vertices per cell
  INTEGER                      :: ne_p_v     !< number of edges per vertex

  TYPE(icon_domain) , ALLOCATABLE, TARGET :: icon_grid_all(:)

  TYPE(geographical_coordinates) :: tpoint

  INTEGER :: start_id 
  INTEGER :: nearest_cell_id

  INTEGER :: nj
  INTEGER :: nb_cell_id
  TYPE(cartesian_coordinates)  :: neighbour_cc     !> coordinates of a neighbour cell centre in cartesian system
  REAL(KIND=wp)                :: sp               !> cos arc length of  of geodesic arc with endpoints x0,x1 
                                                   !> (normalized scalar product of the two points)
  REAL(KIND=wp)                :: sp_max
  TYPE(geographical_coordinates) :: target_geo_co  !> target coordinates in geographical system of point for which 
                                                   !> the nearest ICON grid cell is to be determined
  TYPE(cartesian_coordinates)  :: target_cc_co     !> target coordinates in cartesian system of point for which 
                                                   !> the nearest ICON grid cell is to be determined

  INTEGER, ALLOCATABLE :: nearest_cell_ids(:)    !< array with ids of nearest cell for the domains
  TYPE(cartesian_coordinates), ALLOCATABLE :: polygon(:)
  TYPE(cartesian_coordinates)              :: point
  TYPE(cartesian_coordinates)              :: out_point
  TYPE(geographical_coordinates)           :: out_point_geo
  TYPE(geographical_coordinates), ALLOCATABLE :: poly_geo(:)

  INTEGER                                  :: inflag

  INTEGER                                  :: vert_index
  INTEGER                                  :: ivert

  TYPE(cartesian_coordinates), ALLOCATABLE :: test_poly(:)
  TYPE(cartesian_coordinates)              :: test_point
  TYPE(geographical_coordinates)           :: test_point_geo
  TYPE(cartesian_coordinates)              :: test_out_point
  TYPE(geographical_coordinates)           :: test_out_point_geo
  TYPE(geographical_coordinates), ALLOCATABLE :: test_poly_geo(:)

  INTEGER :: j,k !< counter
  INTEGER (KIND=i8) :: icell

  REAL (KIND=wp) :: undefined
  INTEGER :: undef_int


  INTEGER (KIND=i8) :: nlon_flake !< number of grid elements in zonal direction for flake data
  INTEGER (KIND=i8) :: nlat_flake !< number of grid elements in meridional direction for flake data

  !--------------------------------------------------------------------------------------

  INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid

 ! Print the default information to stdout:
  CALL info_define ('flake_to_buffer')      ! Pre-define the program name as binary name
  CALL info_print ()                     ! Print the information to stdout



  namelist_grid_def = 'INPUT_grid_org'
  CALL init_target_grid(namelist_grid_def)

  PRINT *,'target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat

  igrid_type = tg%igrid_type

  CALL allocate_flake_target_fields(tg)

  print *,'Grid defined, target fields allocated'

  !------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------


  ! get information about FLAE data

  ! get info on raw data file
  input_flake_namelist_file = 'INPUT_FLAKE'

  !---------------------------------------------------------------------------
  CALL read_namelists_extpar_flake(input_flake_namelist_file, &
    &                                      raw_data_flake_path, &
                                           raw_data_flake_filename, &
                                           flake_buffer_file, &
                                           flake_output_file)


  PRINT *,'raw_data_flake_filename: ',TRIM(raw_data_flake_filename)
  flake_file = TRIM(raw_data_flake_path) // TRIM(raw_data_flake_filename)

  PRINT *,'flake file: ', TRIM(flake_file)

  CALL get_dimension_flake_data(flake_file, &
    &                                  nlon_flake, &
    &                                  nlat_flake)

  print *,'nlon_flake: ',nlon_flake
  print *,'nlat_flake: ',nlat_flake


   ! &                             get_lonlat_flake_data
  CALL allocate_raw_flake_fields(nlat_flake,nlon_flake)

  CALL  get_lonlat_flake_data(flake_file, &
                                      nlon_flake, &
                                      nlat_flake, &
                                      lon_flake,  &
                                      lat_flake,  &
                                      flake_grid)

  PRINT *,'MINVAL(lat_flake) :', MINVAL(lat_flake)
  PRINT *,'MAXVAL(lat_flake) :', MAXVAL(lat_flake)

  PRINT *,'flake_grid: ', flake_grid
  PRINT *,'MINVAL(lon_flake) :', MINVAL(lon_flake)
  PRINT *,'MAXVAL(lon_flake) :', MAXVAL(lon_flake)

  PRINT *,'MINVAL(lat_flake) :', MINVAL(lat_flake)
  PRINT *,'MAXVAL(lat_flake) :', MAXVAL(lat_flake)

  undefined = 0.0_wp


  PRINT *,'CALL agg_flake_data_to_target_grid'

  CALL agg_flake_data_to_target_grid(flake_file, &
    &                                      undefined,  &
    &                                      tg,         &
    &                                      lake_depth, &
    &                                      fr_lake,    &
    &                                      flake_tot_npixel)

  PRINT *,'agg_flake_data_to_target_grid done'

  !--------------------------------------------------------------------------------
  ! output

   undefined = -999.0_wp
   undef_int = -999


   netcdf_filename = TRIM(flake_buffer_file)
   PRINT *, 'FLake data buffer filename: ',TRIM(netcdf_filename)

   CALL write_netcdf_buffer_flake(TRIM(netcdf_filename),  &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     lake_depth, &
    &                                     fr_lake,    &
    &                                     flake_tot_npixel)



    SELECT CASE(igrid_type)

      CASE(igrid_icon) ! ICON GRID
        
        netcdf_filename = TRIM(flake_output_file)
        PRINT *,'write out ', TRIM(netcdf_filename)

        CALL write_netcdf_icon_grid_flake(TRIM(netcdf_filename),  &
    &                                     icon_grid,       &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     lake_depth, &
    &                                     fr_lake,    &
    &                                     flake_tot_npixel)




         
      CASE(igrid_cosmo) ! COSMO grid

         netcdf_filename = TRIM(flake_output_file)
         PRINT *,'write out ', TRIM(netcdf_filename)

         CALL write_netcdf_cosmo_grid_flake(TRIM(netcdf_filename), &
    &                                     cosmo_grid,       &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     lake_depth, &
    &                                     fr_lake,    &
    &                                     flake_tot_npixel)



      CASE(igrid_gme) ! GME grid   

    END SELECT

    CALL deallocate_raw_flake_fields()


  PRINT *,'============= flake_to_buffer done ==============='


END PROGRAM extpar_flake_to_buffer
