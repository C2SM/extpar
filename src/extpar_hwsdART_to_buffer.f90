PROGRAM extpar_hwsdART_to_buffer

USE mo_logging
USE info_extpar,         ONLY: info_print
USE mo_io_units,         ONLY: filename_max

USE mo_kind, ONLY: wp, &
             &     i4
 
USE mo_target_grid_data, ONLY: no_raw_data_pixel, &
 & allocate_com_target_fields, &
 & tg
 
USE mo_grid_structures, ONLY: rotated_lonlat_grid, &
  &   icosahedral_triangular_grid, &
  &   target_grid_def, &
  &   igrid_icon, &
  &   igrid_cosmo 


USE  mo_cosmo_grid, ONLY: COSMO_grid, &
  &                         allocate_cosmo_rc, &
  &                         get_cosmo_grid_info, &
  &                         calculate_cosmo_domain_coordinates

  
USE mo_cosmo_grid,     ONLY: calculate_cosmo_target_coordinates


USE mo_icon_grid_data, ONLY: ICON_grid !< structure which contains the definition of the ICON grid


USE mo_base_geometry,  ONLY: geographical_coordinates, &
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

USE mo_math_constants,    ONLY: pi, pi_2, dbl_eps,rad2deg

USE mo_agg_hwsdART,        ONLY: agg_hwsdART_data_to_target_grid


USE mo_read_extpar_namelists, ONLY: read_namelists_extpar_grid_def

USE mo_hwsdART_routines, ONLY: read_namelists_extpar_hwsdART, &
                           &    get_hwsdART_data, &
                           &    get_dimension_hwsdART_data
                               
USE mo_hwsdART_data, ONLY: define_hwsdARTtype, &
                         & hwsdART_grid,       &
                         &  lon_hwsdART,        &
                         &  lat_hwsdART,        &
                         &  hwsdART_soil_unit,  &
                         &  allocate_raw_hwsdART_fields
                           

USE mo_hwsdART_tg_fields, ONLY: allocate_hwsdART_target_fields


USE mo_hwsdART_output_nc, ONLY: write_netcdf_hwsdART_icon_grid, & 
 & write_netcdf_hwsdART_cosmo_grid


USE mo_target_grid_routines, ONLY: init_target_grid

  IMPLICIT NONE


      INTEGER  (KIND=i4) :: nlon_hwsdART  !< number of grid elements in zonal direction for hwsdART raw dataset
      INTEGER  (KIND=i4) :: nlat_hwsdART  !< number of grid elements in meridional direction for hwsdART raw dataset

      CHARACTER(len=filename_max) :: netcdf_filename , &
       &  namelist_hwsdART_data_input , &!< file with input namelist with hwsdART data information
       &  path_hwsdART_file           , &  !< filename with path for hwsdART raw data     
       &  path_deep_hwsdART_file      , &   !< filename with path for hwsdART raw data
       &  hwsdART_buffer_file         , &!< name for hwsdART buffer file
       &  hwsdART_output_file         , & !< name for hwsdART output file
       &  raw_data_hwsdART_path       , &    !< path to raw data
       &  raw_data_hwsdART_filename   , &!< filename hwsdART raw data
       &  namelist_grid_def                 !< filename with namelists for grid settings for EXTPAR



      REAL(KIND=wp) :: undefined !< value to indicate undefined grid elements in cosmo_ndvi_field
      INTEGER (KIND=i4) :: undefined_integer   !< value for undefined integer
      INTEGER (KIND=i4) :: default_value !< default value

      INTEGER :: errorcode

      undefined_integer = 0 ! set undefined to zero
      undefined = -99.0 ! undef vlaue
      default_value =  3 ! default value

 
  CALL initialize_logging("extpar_hwsdART_to_buffer.log")
  CALL info_print ()
  !--------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= start hwsdART_to_buffer =============')
  CALL logging%info( '')

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= read namelist and init grid ======')
  CALL logging%info( '')

  !--------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------

      
      namelist_grid_def = 'INPUT_grid_org'
      CALL  init_target_grid(namelist_grid_def)
      
      namelist_hwsdART_data_input = 'INPUT_hwsdART'
      CALL read_namelists_extpar_hwsdART(namelist_hwsdART_data_input,        &
                                           raw_data_hwsdART_path,         &
                                           raw_data_hwsdART_filename,     &
                                           hwsdART_output_file)

                                           

      WRITE(message_text,*) 'raw_data_hwsdART_path: ', TRIM(raw_data_hwsdART_path)
      CALL logging%info(message_text) 

      WRITE(message_text,*) 'raw_data_hwsdART_filename: ', TRIM(raw_data_hwsdART_filename)
      CALL logging%info(message_text) 

      WRITE(message_text,*) 'hwsdART_output_file: ', TRIM(hwsdART_output_file)
      CALL logging%info(message_text) 

      path_hwsdART_file = TRIM(raw_data_hwsdART_path) // TRIM(raw_data_hwsdART_filename)

      WRITE(message_text,*)  'path_hwsdART_file: ', TRIM(path_hwsdART_file)
      CALL logging%info(message_text) 
      ! inquire dimensions from raw data file

      CALL  get_dimension_hwsdART_data(path_hwsdART_file,  &
                                          nlon_hwsdART, &
                                          nlat_hwsdART)


      WRITE(message_text,*)  'nlon_hwsdART', nlon_hwsdART
      CALL logging%info(message_text) 
      WRITE(message_text,*)  'nlat_hwsdART', nlat_hwsdART
      CALL logging%info(message_text) 

      ! define value of global variables hwsdART types
      !--------------------------------------------------------------------------------------------------------
      CALL define_hwsdARTtype()
      WRITE(message_text,*)   'define_hwsdARTtype done'
      CALL logging%info(message_text)  

      CALL allocate_raw_hwsdART_fields(nlon_hwsdART, nlat_hwsdART)
      WRITE(message_text,*)   'allocate_raw_hwsdART_fields done'
      CALL logging%info(message_text)  

      CALL get_hwsdART_data(path_hwsdART_file)
      WRITE(message_text,*)   'get_hwsdART_data'
      CALL logging%info(message_text)  

      CALL allocate_hwsdART_target_fields(tg)
      WRITE(message_text,*)   'allocate_hwsdART_target_fields done'
      CALL logging%info(message_text)  
      
      WRITE(message_text,*) 'hwsdART read from file ', TRIM(path_hwsdART_file)
      CALL logging%info(message_text)
  
      
      ! aggregate hwsdART data to target grid
      WRITE(message_text,*) 'aggregate hwsdART data to target grid'
      CALL logging%info(message_text)  
      undefined = 0.0_wp
      
      
      CALL agg_hwsdART_data_to_target_grid(tg,              &
                  &                   hwsdART_soil_unit, &
                  &                   hwsdART_grid,      &
                  &                   lon_hwsdART,       &
                  &                   lat_hwsdART)


      
      WRITE(message_text,*) 'MAXVAL(no_raw_data_pixel): ', MAXVAL(no_raw_data_pixel)
      CALL logging%info(message_text)

      WRITE(message_text,*) 'MINVAL(no_raw_data_pixel): ', MINVAL(no_raw_data_pixel)
      CALL logging%info(message_text)

      DEALLOCATE (hwsdART_soil_unit, STAT = errorcode)
      IF (errorcode /= 0) WRITE(message_text,*)  'Cant deallocate hwsdART_soil_unit'


      WRITE(MESSAGE_TEXT,*) 'Start target grid output'
      CALL logging%info(message_text)  

      undefined = -999._wp

      SELECT CASE(tg%igrid_type)
       !-----------------------------------------------------------------
       CASE(igrid_icon) ! ICON GRID

         netcdf_filename= TRIM(hwsdART_output_file)

         CALL write_netcdf_hwsdART_icon_grid(netcdf_filename,  &
                                          & icon_grid,         &
                                          & tg,                &
                                          & undefined          )


       CASE(igrid_cosmo) ! COSMO grid

         netcdf_filename= TRIM(hwsdART_output_file)

          CALL write_netcdf_hwsdART_cosmo_grid(netcdf_filename, &
                                          & cosmo_grid,         &
                                          & tg,                 &
                                          & undefined           )
    END SELECT
                           
    CALL logging%info('============= hwsdART_to_buffer done ===============')
        

END PROGRAM extpar_hwsdART_to_buffer
