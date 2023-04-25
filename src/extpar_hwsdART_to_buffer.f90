!+ Fortran main program to read in hwsdART data and aggregate to target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2014/04/24 Daniel Rieger
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran main program to read in hwsdART data and aggregate to target grid
!>  
!> \author Daniel Rieger
PROGRAM extpar_hwsdART_to_buffer

! Load the library information data:
!USE info_extpar, ONLY: info_define, info_readnl, info_print


USE mo_kind, ONLY: wp
USE mo_kind, ONLY: i4
USE mo_kind, ONLY: i8

USE mo_target_grid_data, ONLY: no_raw_data_pixel
USE mo_target_grid_data, ONLY: lon_geo
USE mo_target_grid_data, ONLY: lat_geo

USE mo_target_grid_data, ONLY: allocate_com_target_fields
USE mo_target_grid_data, ONLY: tg
 
USE mo_grid_structures, ONLY: rotated_lonlat_grid
USE mo_grid_structures, ONLY: icosahedral_triangular_grid
USE mo_grid_structures, ONLY: target_grid_def

USE mo_grid_structures, ONLY: igrid_icon
USE mo_grid_structures, ONLY: igrid_cosmo
!USE mo_grid_structures, ONLY: igrid_gme

USE  mo_cosmo_grid, ONLY: COSMO_grid, &
  &                         lon_rot, &
  &                         lat_rot, &
  &                         allocate_cosmo_rc, &
  &                         get_cosmo_grid_info, &
  &                         calculate_cosmo_domain_coordinates

  
USE  mo_cosmo_grid, ONLY: calculate_cosmo_target_coordinates


USE  mo_icon_grid_data, ONLY: ICON_grid !< structure which contains the definition of the ICON grid


USE mo_base_geometry,    ONLY:  geographical_coordinates, &
                                 cartesian_coordinates
  
USE mo_additional_geometry,   ONLY: cc2gc,                  &
                              gc2cc,                  &
                              arc_length,             &
                              cos_arc_length,         &
                              inter_section,          &
                              vector_product,         &
                              point_in_polygon_sp

                              

USE mo_icon_domain,          ONLY: icon_domain, &
                              grid_cells,               &
                              grid_vertices,            &
                              construct_icon_domain,    &
                              destruct_icon_domain

USE mo_io_units,          ONLY: filename_max

USE mo_exception,         ONLY: message_text, message, finish

USE mo_utilities_extpar, ONLY: abort_extpar

USE mo_math_constants,  ONLY: pi, pi_2, dbl_eps,rad2deg

USE mo_agg_hwsdART, ONLY: agg_hwsdART_data_to_target_grid
!USE mo_agg_soil, ONLY: agg_soil_data_to_target_grid, &
!                       nearest_soil_data_to_target_grid

USE mo_read_extpar_namelists, ONLY: read_namelists_extpar_grid_def

USE mo_hwsdART_routines, ONLY: read_namelists_extpar_hwsdART, &
                               get_hwsdART_data, &
                               get_dimension_hwsdART_data
                               
USE mo_hwsdART_data, ONLY: define_hwsdARTtype, &
                           hwsdART_data,       &
                           hwsdART_grid,       &
                           undef_hwsdARTtype,  &
                           default_hwsdARTtype,&
                           no_data,            &
                           type_clay_heavy,    &
                           type_silty_clay,    &
                           type_clay_light,    &
                           type_silty_clay_loam, &
                           type_clay_loam,     &
                           type_silt,          &
                           type_silt_loam,     &
                           type_sandy_clay,    &
                           type_loam,          &
                           type_sandy_clay_loam, &
                           type_sandy_loam,    &
                           type_loamy_sand,    &
                           type_sand,          &
                           lon_hwsdART,        &
                           lat_hwsdART,        &
                           hwsdART_soil_unit,  &
                           allocate_raw_hwsdART_fields
                           

USE mo_hwsdART_tg_fields, ONLY: allocate_hwsdART_target_fields

USE mo_hwsdART_output_nc, ONLY: write_netcdf_hwsdART_icon_grid
USE mo_hwsdART_output_nc, ONLY: write_netcdf_hwsdART_cosmo_grid

USE mo_target_grid_routines, ONLY: init_target_grid

  IMPLICIT NONE


      INTEGER  (KIND=i4) :: nlon_hwsdART  !< number of grid elements in zonal direction for hwsdART raw dataset
      INTEGER  (KIND=i4) :: nlat_hwsdART  !< number of grid elements in meridional direction for hwsdART raw dataset

      CHARACTER(len=filename_max) :: netcdf_filename
 
      CHARACTER (len=filename_max) :: input_namelist_cosmo_grid !< file with input namelist with COSMO grid definition
      CHARACTER (len=filename_max) :: namelist_hwsdART_data_input !< file with input namelist with hwsdART data information
    
      CHARACTER (len=filename_max) :: raw_data_path        !< path to raw data
      CHARACTER (len=filename_max) :: path_hwsdART_file      !< filename with path for hwsdART raw data     
      CHARACTER (len=filename_max) :: path_deep_hwsdART_file      !< filename with path for hwsdART raw data
      CHARACTER (len=filename_max) :: netcdf_out_filename      !< filename for netcdf file with hwsdART data on COSMO grid
      CHARACTER (len=filename_max) :: hwsdART_buffer_file  !< name for hwsdART buffer file
      CHARACTER (len=filename_max) :: hwsdART_output_file  !< name for hwsdART output file

      CHARACTER (len=filename_max) :: raw_data_hwsdART_path        !< path to raw data
      CHARACTER (len=filename_max) :: raw_data_hwsdART_filename !< filename hwsdART raw data
      CHARACTER (len=filename_max) :: raw_data_deep_hwsdART_filename !< filename deep hwsdART raw data

      CHARACTER (len=filename_max) :: namelist_grid_def !< filename with namelists for grid settings for EXTPAR
      CHARACTER (len=filename_max) :: domain_def_namelist !< namelist file with domain definition

      
      REAL (KIND=wp) :: point_lon_geo !< longitude of a point in geographical system
      REAL (KIND=wp) :: point_lat_geo !< latitude of a point in geographical system


      REAL(KIND=wp) :: undefined !< value to indicate undefined grid elements in cosmo_ndvi_field
      INTEGER (KIND=i4) :: undefined_integer   !< value for undefined integer

      INTEGER (KIND=i4) :: default_value !< default value


      INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid

      INTEGER :: idom  !< ICON Domain Number

      INTEGER :: i,j,k !< counters
      INTEGER :: errorcode




      ! Print the default information to stdout:
!      CALL info_define ('hwsdART_to_buffer')      ! Pre-define the program name as binary name
!      CALL info_print ()                     ! Print the information to stdout
      !--------------------------------------------------------------------------------------------------------

      undefined_integer = 0 ! set undefined to zero
      undefined = -99.0 ! undef vlaue
      default_value =  3 ! default value
      path_deep_hwsdART_file = "" !default name

      !--------------------------------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------------------------------
      ! get information on target grid, allocate target fields with coordinates and determin the coordinates 
      ! for th target grid
      
      namelist_grid_def = 'INPUT_grid_org'
      CALL  init_target_grid(namelist_grid_def)
      
      !--------------------------------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------------------------------
      
      
      ! get information on hwsdART raw data
      !--------------------------------------------------------------------------------------------------------
        
      ! read namelist with hwsdART data information (path, filename)
      
      namelist_hwsdART_data_input = 'INPUT_hwsdART'
      CALL read_namelists_extpar_hwsdART(namelist_hwsdART_data_input,        &
                                           raw_data_hwsdART_path,         &
                                           raw_data_hwsdART_filename,     &
                                           hwsdART_buffer_file,           &
                                           hwsdART_output_file)

                                           

           

      print *,'raw_data_hwsdART_path: ', TRIM(raw_data_hwsdART_path)
      print *,'raw_data_hwsdART_filename: ', TRIM(raw_data_hwsdART_filename)
      print *,'hwsdART_output_file: ', TRIM(hwsdART_output_file)

      path_hwsdART_file = TRIM(raw_data_hwsdART_path) // TRIM(raw_data_hwsdART_filename)

      print *, 'path_hwsdART_file: ', TRIM(path_hwsdART_file)

      ! inquire dimensions from raw data file

      CALL  get_dimension_hwsdART_data(path_hwsdART_file,  &
                                          nlon_hwsdART, &
                                          nlat_hwsdART)


      print *, 'nlon_hwsdART', nlon_hwsdART
      print *, 'nlat_hwsdART', nlat_hwsdART
      

      ! define value of global variables hwsdART types
      !--------------------------------------------------------------------------------------------------------
      CALL define_hwsdARTtype()
      print*, 'define_hwsdARTtype done'
      CALL allocate_raw_hwsdART_fields(nlon_hwsdART, nlat_hwsdART)
      print*, 'allocate_raw_hwsdART_fields done'
      CALL get_hwsdART_data(path_hwsdART_file)
      print*, 'get_hwsdART_data'
      CALL allocate_hwsdART_target_fields(tg)
      print*, 'allocate_hwsdART_target_fields done'
      
      print *,'hwsdART read from file ', TRIM(path_hwsdART_file)
      
      ! aggregate hwsdART data to target grid
      print *,'aggregate hwsdART data to target grid'
      undefined = 0.0_wp
      
      
      CALL agg_hwsdART_data_to_target_grid(tg,              &
                  &                   hwsdART_soil_unit, &
                  &                   hwsdART_grid,      &
                  &                   lon_hwsdART,       &
                  &                   lat_hwsdART)


      
      print *,'MAXVAL(no_raw_data_pixel): ', MAXVAL(no_raw_data_pixel)
      print *,'MINVAL(no_raw_data_pixel): ', MINVAL(no_raw_data_pixel)

      DEALLOCATE (hwsdART_soil_unit, STAT = errorcode)
      IF (errorcode /= 0) print*, 'Cant deallocate hwsdART_soil_unit'

!      PRINT *,'Start buffer output'


!      netcdf_filename=  TRIM(hwsdART_buffer_file)

!      undefined = -999.0_wp
!      undefined_integer= 999

!        CALL write_netcdf_soil_buffer(netcdf_filename,   &
!   &                                   tg,               &
!   &                                   isoil_data,       &
!   &                                   ldeep_soil,       &
!   &                                   undefined,        &
!   &                                   undefined_integer,&
!   &                                   lon_geo,          &
!   &                                   lat_geo,          &
!   &                                   fr_land_soil,     &
!   &                                   soiltype_fao,     &
!   &                                   soiltype_deep = soiltype_deep)

!      PRINT *,'buffer output done'
      PRINT *,'Start target grid output'

      undefined = -999._wp

      SELECT CASE(tg%igrid_type)
       !-----------------------------------------------------------------
       CASE(igrid_icon) ! ICON GRID

         netcdf_filename= TRIM(hwsdART_output_file)

         CALL write_netcdf_hwsdART_icon_grid(netcdf_filename,  &
                                          & icon_grid,         &
                                          & tg,                &
                                          & undefined,         &
                                          & undefined_integer, &
                                          & lon_geo,           &
                                          & lat_geo)
       CASE(igrid_cosmo) ! COSMO grid

         netcdf_filename= TRIM(hwsdART_output_file)

          CALL write_netcdf_hwsdART_cosmo_grid(netcdf_filename, &
                                          & cosmo_grid,         &
                                          & tg,                 &
                                          & undefined,          &
                                          & undefined_integer,  &
                                          & lon_geo,            &
                                          & lat_geo)
    END SELECT
                           
  PRINT *,'============= hwsdART_to_buffer done ==============='

        

END PROGRAM extpar_hwsdART_to_buffer
