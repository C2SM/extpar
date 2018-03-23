!+  Fortran main program to aggregate the CRU near surface climatology to target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio
!  small bug fixes accroding to Fortran compiler warnings
! V1_7         2013/01/25 Guenther Zaengl
!   Parallel threads for ICON and COSMO using Open-MP,
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON         
! V2_0         2013/06/04 Martina Messmer
!  introduction of a finer temperature climatology (CRU) and 
!  CRU temperature elevation (CLM Community)
!  there is a switch to choose between the fine and coarse data set
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!>  Fortran main program to aggregate the CRU near surface climatology to target grid
!>
!! @par extpar_cru_to_buffer 
!!
!!  This program interpolates the CRU near surface climatology to a given target grid (COSMO/ICON).
!!
!! @author
!!     Hermann Asensio
!!     (DWD)
!!
PROGRAM extpar_cru_to_buffer

! Load the library information data:
  USE info_extpar, ONLY: info_define, info_print


  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  USE mo_grid_structures, ONLY: target_grid_def,   &
                                reg_lonlat_grid,   &
                                rotated_lonlat_grid
  
  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo
  USE mo_grid_structures, ONLY: igrid_gme


  USE mo_target_grid_data, ONLY: lon_geo, &
    &                            lat_geo, &
    &                            no_raw_data_pixel, &
    &                            allocate_com_target_fields
  
  USE mo_target_grid_data, ONLY: tg

  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_icon_grid_data, ONLY: ICON_grid !< structure which contains the definition of the ICON grid
 
  USE  mo_cosmo_grid, ONLY: COSMO_grid, &
    &                       lon_rot, &
    &                       lat_rot, &
    &                       allocate_cosmo_rc, &
    &                       get_cosmo_grid_info, &
    &                       calculate_cosmo_domain_coordinates

  USE mo_base_geometry,    ONLY:  geographical_coordinates, &
    &                                cartesian_coordinates

                                  
  
  USE mo_icon_domain,          ONLY: icon_domain, &
    &                             grid_cells,               &
    &                             grid_vertices,            &
    &                             construct_icon_domain,    &
    &                             destruct_icon_domain

  USE mo_io_units,          ONLY: filename_max

  USE mo_exception,         ONLY: message_text, message, finish

  USE mo_utilities_extpar, ONLY: abort_extpar


  USE mo_additional_geometry,   ONLY: cc2gc,             &
    &                            gc2cc,                  &
    &                            arc_length,             &
    &                            cos_arc_length,         &
    &                            inter_section,          &
    &                            vector_product,         &
    &                            point_in_polygon_sp


  USE mo_math_constants,  ONLY: pi, pi_2, dbl_eps,rad2deg

  USE mo_cru_data, ONLY: read_namelists_extpar_t_clim

  USE mo_cru_data, ONLY : allocate_cru_data, &
    &                     deallocate_cru_data, &
    &                     read_cru_data_input_namelist, &
    &                     get_dimension_cru_data, &
    &                     get_cru_grid_and_data, &
    &                     lon_cru, &
    &                     lat_cru, &
    &                     cru_raw_data, &
    &                     cru_grid
 
  USE mo_cru_target_fields, ONLY: allocate_cru_target_fields,&
    &                              crutemp, &
    &                              meta_crutemp, &
    &                              cruelev,      &
    &                              meta_cruelev, &
    &                              i_t_cru_fine, &
    &                              i_t_cru_coarse

  USE mo_agg_cru, ONLY: agg_cru_data_to_target_grid

  USE mo_cru_output_nc, ONLY: write_netcdf_buffer_cru, &
    &                         write_netcdf_cosmo_grid_cru, &
    &                         write_netcdf_icon_grid_cru  

  IMPLICIT NONE
  
  CHARACTER(len=filename_max) :: filename
  CHARACTER(len=filename_max) :: netcdf_filename
  CHARACTER(len=filename_max) :: input_namelist_file
  CHARACTER(len=filename_max) :: input_namelist_cosmo_grid !< file with input namelist with COSMO grid definition

  CHARACTER(len=filename_max) :: namelist_file
  CHARACTER(len=filename_max) :: namelist_grid_def



  CHARACTER(len=filename_max) :: input_glc2000_namelist_file 
  CHARACTER(len=filename_max) :: glc2000_file

  INTEGER (KIND=i8):: it_cl_type !< integer switch to choose a land use raw data set
                                        !! 1 CRU (fine), 2 CRU (coarse)

  CHARACTER (len=filename_max) :: raw_data_t_clim_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_t_clim_filename !< filename temperature climatology raw data


  CHARACTER (len=filename_max) :: t_clim_buffer_file !< name for temperature climatology buffer
  CHARACTER (len=filename_max) :: t_clim_output_file !< name for temperature climatology output file



  INTEGER :: i, ip, ic, in


  INTEGER                      :: i_nc       !< number of cells
  INTEGER                      :: i_ne       !< number of edges
  INTEGER                      :: i_nv       !< number of vertices
  INTEGER                      :: nc_p_e     !< number of cells per edge
  INTEGER                      :: nv_p_c     !< number of vertices per cell
  INTEGER                      :: ne_p_v     !< number of edges per vertex



  TYPE(icon_domain) , ALLOCATABLE, TARGET :: icon_grid_all(:)

  TYPE(geographical_coordinates) :: tpoint

  INTEGER :: nearest_cell_id

  INTEGER :: nj
  INTEGER :: nb_cell_id
  TYPE(cartesian_coordinates)  :: neighbour_cc     !> coordinates of a neighbour cell centre in cartesian system
  REAL(KIND=wp)                :: sp               !> cos arc length of  of geodesic arc with endpoints x0,x1
  REAL(KIND=wp)                :: sp_max
  TYPE(geographical_coordinates) :: target_geo_co    !> target coordinates in geographical system of point 
  TYPE(cartesian_coordinates)  :: target_cc_co     !>  target coordinates in cartesian system of point 

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
  INTEGER :: l,m               !< counter
  INTEGER (KIND=i8) :: icell

  REAL (KIND=wp) :: undefined
  INTEGER :: undef_int


  !--------------------------------------------------------------------------------------
  INTEGER (KIND=i8) :: nlon_glc2000 !< number of grid elements in zonal direction for glc2000 data
  INTEGER (KIND=i8) :: nlat_glc2000 !< number of grid elements in meridional direction for glc2000 data
  !--------------------------------------------------------------------------------------

  
  INTEGER (KIND=i8) :: ntype !< number of types of aerosols
  INTEGER (KIND=i8) :: nrows !< number of rows
  INTEGER (KIND=i8) :: ncolumns !< number of columns
  INTEGER (KIND=i8) :: ntime !< number of times


  INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid

  ! Print the default information to stdout:
  CALL info_define ('cru_to_buffer')      ! Pre-define the program name as binary name
  CALL info_print ()                     ! Print the information to stdout

  !--------------------------------------------------------------------------------------------------------
  ! get information on target grid, allocate target fields with coordinates and determin the coordinates 
  ! for th target grid

  namelist_grid_def = 'INPUT_grid_org'
  CALL  init_target_grid(namelist_grid_def)

  PRINT *,'target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat

  igrid_type = tg%igrid_type


  !------------------------------------------------------------------------------------


  ! get information about temperature climatology data

  namelist_file = 'INPUT_TCLIM'
  CALL  read_namelists_extpar_t_clim(namelist_file, &
                                         it_cl_type, &
                                         raw_data_t_clim_path, &
                                         raw_data_t_clim_filename, &
                                         t_clim_buffer_file, &
                                         t_clim_output_file)

  filename = TRIM(raw_data_t_clim_path) // TRIM(raw_data_t_clim_filename)

  PRINT *,'filename ',TRIM(filename)

  ! inquire dimensions

  CALL  get_dimension_cru_data(filename,     &
    &                          nrows,        &
    &                          ncolumns,     &
    &                          ntime)

    PRINT *, 'nrows: ',nrows
    PRINT *, 'ncolumns: ',ncolumns
    PRINT *, 'ntime: ',ntime

    CALL allocate_cru_data(nrows,ncolumns,ntime)


    ! read in aot raw data
    CALL get_cru_grid_and_data(filename,     &
         &                     it_cl_type,   &
         &                     nrows,        &
         &                     ncolumns,     &
         &                     ntime)
                                     
    PRINT *, 'cru_grid: ', cru_grid

    PRINT *,'target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat


    ! allocate target grid fields for aerosol optical thickness

    CALL allocate_cru_target_fields(tg)


    undefined = -999.0_wp
    undef_int = -999

    crutemp  =  undefined  ! set target grid values to undefined
    cruelev  = undefined

    PRINT *,'agg_cru_data_to_target_grid'
    CALL  agg_cru_data_to_target_grid(nrows,ncolumns,ntime,it_cl_type)

    PRINT *,'aggregation done'


    !write out data
    filename = TRIM(t_clim_output_file)
    

    SELECT CASE(igrid_type)

      CASE(igrid_icon) ! ICON GRID
        
        netcdf_filename = TRIM(t_clim_output_file)
        PRINT *,'write out ', TRIM(netcdf_filename)

        SELECT CASE (it_cl_type)
        CASE(i_t_cru_coarse)
         
        CALL write_netcdf_icon_grid_cru(netcdf_filename,  &
    &                                     icon_grid,       &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     crutemp)
  
        CASE(i_t_cru_fine)

        CALL write_netcdf_icon_grid_cru(netcdf_filename,  &
    &                                     icon_grid,       &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     crutemp, &
    &                                     cruelev=cruelev)

      END SELECT

      CASE(igrid_cosmo) ! COSMO grid

         netcdf_filename = TRIM(t_clim_output_file)
         PRINT *,'write out ', TRIM(netcdf_filename)

        SELECT CASE (it_cl_type)
        CASE(i_t_cru_coarse)
        CALL write_netcdf_cosmo_grid_cru(netcdf_filename,  &
    &                                     cosmo_grid,       &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     crutemp)
     

       CASE(i_t_cru_fine)
         CALL write_netcdf_cosmo_grid_cru(netcdf_filename,  &
    &                                     cosmo_grid,       &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     crutemp, &
    &                                     cruelev=cruelev)

      END SELECT


      CASE(igrid_gme) ! GME grid   

    END SELECT

    !write out data
    netcdf_filename = TRIM(t_clim_buffer_file)
    PRINT *,'write out ', TRIM(netcdf_filename)

    SELECT CASE (it_cl_type)
    CASE(i_t_cru_coarse)
    CALL write_netcdf_buffer_cru(netcdf_filename,  &
   &                                     tg,         &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     crutemp)

    CASE(i_t_cru_fine)
      CALL write_netcdf_buffer_cru(netcdf_filename,  &
   &                                     tg,         &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     crutemp, &
   &                                     cruelev=cruelev)

    END SELECT



    CALL deallocate_cru_data()



  PRINT *,'============= cru_to_buffer done ==============='
 
  

END PROGRAM extpar_cru_to_buffer


