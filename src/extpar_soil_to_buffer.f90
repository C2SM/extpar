!+ Fortran main program to read in soil data and aggregate to target grid
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
! V2_0         2013/06/04 Martina Messmer
!   introduction of the HWSD data sets (topsoil and subsoil) for the 
!   external parameters (code for the topsoil obtained from Juergen Helmert)
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran main program to read in soil data and aggregate to target grid
!>  
!> \author Hermann Asensio
PROGRAM extpar_soil_to_buffer

! Load the library information data:
USE info_extpar, ONLY: info_define, info_print


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
USE mo_grid_structures, ONLY: igrid_gme

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

USE mo_agg_soil, ONLY: agg_soil_data_to_target_grid, &
                       nearest_soil_data_to_target_grid

USE mo_read_extpar_namelists, ONLY: read_namelists_extpar_grid_def
USE mo_soil_routines, ONLY: read_namelists_extpar_soil

                      
USE mo_soil_routines, ONLY: get_soil_data, &
                            get_deep_soil_data, &
                            get_dimension_soil_data

USE mo_soil_routines, ONLY: &
        nlon_soil,          &
        nlat_soil

USE mo_soil_data, ONLY: allocate_raw_soil_fields, &
        allocate_raw_deep_soil_fields,            &
        define_soiltype,    &
        dsmw_legend,        &
        soil_texslo,        &
        soil_texslo_deep,   & 
        dsmw_soil_unit,     &
        dsmw_deep_soil_unit,&
        n_unit,             &
        dsmw_grid,          &
        lon_soil,           &
        lat_soil

USE mo_soil_data, ONLY: FAO_data, HWSD_data, HWSD_map, soil_data

USE mo_soil_data, ONLY:     &
       lon_full,            &
       lat_full

USE   mo_soil_tg_fields, ONLY:  fr_land_soil
USE   mo_soil_tg_fields, ONLY:  soiltype_fao,soiltype_hwsd, soiltype_deep,soiltype_hwsd_s
USE   mo_soil_tg_fields, ONLY:  allocate_soil_target_fields

USE mo_soil_output_nc, ONLY: write_netcdf_soil_cosmo_grid
USE mo_soil_output_nc, ONLY: write_netcdf_soil_icon_grid
USE mo_soil_output_nc, ONLY: write_netcdf_soil_buffer



USE mo_target_grid_routines, ONLY: init_target_grid

  IMPLICIT NONE

      CHARACTER(len=filename_max) :: netcdf_filename
 
      CHARACTER (len=filename_max) :: input_namelist_cosmo_grid !< file with input namelist with COSMO grid definition
      CHARACTER (len=filename_max) :: namelist_soil_data_input !< file with input namelist with soil data information
    
      CHARACTER (len=filename_max) :: raw_data_path        !< path to raw data
      CHARACTER (len=filename_max) :: path_soil_file      !< filename with path for soil raw data     
      CHARACTER (len=filename_max) :: path_deep_soil_file      !< filename with path for soil raw data
      CHARACTER (len=filename_max) :: netcdf_out_filename      !< filename for netcdf file with soil data on COSMO grid
      CHARACTER (len=filename_max) :: soil_buffer_file  !< name for soil buffer file
      CHARACTER (len=filename_max) :: soil_output_file  !< name for soil output file
      CHARACTER (len=filename_max) :: soil_buffer_file_consistent !< name for soil buffer file after consistency check
      CHARACTER (len=filename_max) :: soil_output_file_consistent !< name for soil output file after consistency check

      CHARACTER (len=filename_max) :: raw_data_soil_path        !< path to raw data
      CHARACTER (len=filename_max) :: raw_data_soil_filename !< filename soil raw data
      CHARACTER (len=filename_max) :: raw_data_deep_soil_filename !< filename deep soil raw data

      CHARACTER (len=filename_max) :: namelist_grid_def !< filename with namelists for grid settings for EXTPAR
      CHARACTER (len=filename_max) :: domain_def_namelist !< namelist file with domain definition

      
      REAL (KIND=wp) :: point_lon_geo !< longitude of a point in geographical system
      REAL (KIND=wp) :: point_lat_geo !< latitude of a point in geographical system


      REAL(KIND=wp) :: undefined !< value to indicate undefined grid elements in cosmo_ndvi_field
      INTEGER (KIND=i4) :: undefined_integer   !< value for undefined integer

      INTEGER (KIND=i4) :: default_value !< default value


      INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid
      INTEGER (KIND=i4) :: isoil_data  !< soil data, 1 for FAO raw data, 
                                       !             2 for HWSD raw data,
                                       !             3 for HWSD terra mapping

      INTEGER (KIND=i4) :: undef_soiltype
      INTEGER (KIND=i4) :: default_soiltype
      INTEGER (KIND=i4) :: soiltype_ice
      INTEGER (KIND=i4) :: soiltype_water

      INTEGER :: idom  !< ICON Domain Number

      INTEGER :: i,j,k !< counters
      INTEGER :: errorcode

      LOGICAL :: ldeep_soil            !< switch to decide weather the deep soil layer is desired or not

      INTEGER  (KIND=i4) :: nlon_full
      INTEGER  (KIND=i4) :: nlat_full
      INTEGER  (KIND=i4) :: lon_low, lat_low, lon_hig, lat_hig
      INTEGER  (KIND=i4) :: start(2)


      ! Print the default information to stdout:
      CALL info_define ('soil_to_buffer')      ! Pre-define the program name as binary name
      CALL info_print ()                     ! Print the information to stdout
      !--------------------------------------------------------------------------------------------------------

      undefined_integer = 0 ! set undefined to zero
      undefined = -99.0 ! undef vlaue
      default_value =  3 ! default value
      path_deep_soil_file = "" !default name

      !--------------------------------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------------------------------
      ! get information on target grid, allocate target fields with coordinates and determin the coordinates 
      ! for th target grid
      
      namelist_grid_def = 'INPUT_grid_org'
      CALL  init_target_grid(namelist_grid_def)
      
      !--------------------------------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------------------------------
      
      
      ! get information on soil raw data
      !--------------------------------------------------------------------------------------------------------
        
      ! read namelist with soil data information (path, filename)
      
      namelist_soil_data_input = 'INPUT_SOIL'
      CALL read_namelists_extpar_soil(namelist_soil_data_input,        &
                                           isoil_data,                 &
                                           ldeep_soil,                 &
                                           raw_data_soil_path,         &
                                           raw_data_soil_filename,     &
                                           raw_data_deep_soil_filename,&
                                           soil_buffer_file,           &
                                           soil_output_file,           &
                                           soil_buffer_file_consistent,&
                                           soil_output_file_consistent)
      
       
      
      
      IF (ldeep_soil .AND. isoil_data /= HWSD_data) THEN !_br 21.02.14 replace eq by eqv
        ldeep_soil = .FALSE.
        print*, '********* you can only use the deep soil if HWSD data is used *********'
        print*, '********* ldeep_soil is set to FALSE *********'
      ENDIF


      !HA debug
      print *,'raw_data_soil_path: ', TRIM(raw_data_soil_path)
      print *,'raw_data_soil_filename: ', TRIM(raw_data_soil_filename)
      IF(ldeep_soil) print *,'raw_data_deep_soil_filename: ', TRIM(raw_data_deep_soil_filename) 
      print *,'soil_output_file: ', TRIM(soil_output_file)

      path_soil_file = TRIM(raw_data_soil_path) // TRIM(raw_data_soil_filename)
      IF (ldeep_soil) THEN
        path_deep_soil_file = TRIM(raw_data_soil_path) // TRIM(raw_data_deep_soil_filename)
      ENDIF

      !HA debug
      print *, 'path_soil_file: ', TRIM(path_soil_file)
      IF (ldeep_soil) print *, 'path_deep_soil_file: ', TRIM(path_deep_soil_file)

      ! inquire dimensions from raw data file

      CALL  get_dimension_soil_data(path_soil_file,  &
                                          nlon_full, &
                                          nlat_full, &
                                          n_unit)

      
      ! determine section of full data covered by target domain

      IF (MAXVAL(lon_geo) > 180.0) THEN
        nlon_soil = nlon_full
        lon_low = 1
        lon_hig = nlon_full
      ELSE
        lon_low = MAX(1,INT((MINVAL(lon_geo)-dsmw_grid%dlon_reg-0.5_wp-lon_full(1))/dsmw_grid%dlon_reg))
        lon_hig = MIN(nlon_full,INT((MAXVAL(lon_geo)+dsmw_grid%dlon_reg+0.5_wp-lon_full(1))/dsmw_grid%dlon_reg))
        nlon_soil = lon_hig + 1 - lon_low
      ENDIF

!     latitude runs from north to south in the raw data
      lat_low = MAX(1,INT((lat_full(1)-MAXVAL(lat_geo)-0.5_wp+dsmw_grid%dlat_reg)/ABS(dsmw_grid%dlat_reg)))
      lat_hig = MIN(nlat_full,INT((lat_full(1)-MINVAL(lat_geo)+0.5_wp-dsmw_grid%dlat_reg)/ABS(dsmw_grid%dlat_reg)))
      nlat_soil = lat_hig + 1 - lat_low
      start(1) = lon_low
      start(2) = lat_low

      !HA debug
      print *, 'nlon_soil', nlon_soil
      print *, 'nlat_soil', nlat_soil
      print *, 'n_unit', n_unit



      ! get coordinates and legend and data from raw data file
      ! define value of global variables soil_raw_data_grid, lon_reg, lat_reg, soil_texslo, dsmw_soil_unit
      !--------------------------------------------------------------------------------------------------------
      CALL define_soiltype(isoil_data, ldeep_soil, &
                           undef_soiltype,         &
                           default_soiltype,       &
                           soiltype_ice,           &
                           soiltype_water,         &
                           soil_data)
      print*, 'define_soiltype done'
      CALL allocate_raw_soil_fields(nlon_soil, nlat_soil, n_unit)
      print*, 'allocate_raw_soil_fields done'
      CALL get_soil_data(path_soil_file, start)
      lon_soil = lon_full(lon_low:lon_hig)
      lat_soil = lat_full(lat_low:lat_hig)

      CALL get_soil_data(path_soil_file,start)
      print*, 'get_soil_data'
      CALL allocate_soil_target_fields(tg, ldeep_soil)
      print*, 'allocate_soil_target_fields done'

      !--------------------------------------------------------------------------------------------------------

      !HA debug
      SELECT CASE(isoil_data)
      CASE(FAO_data)
        print *,'FAO DSMW read from file ', TRIM(path_soil_file)
      CASE(HWSD_data, HWSD_map)
        print *,'HWSD read from file ', TRIM(path_soil_file)
      END SELECT


      ! aggregate soil data to target grid

      print *,'aggregate soil data to target grid'

      undefined = 0.0

      CALL agg_soil_data_to_target_grid(tg,              &
                  &                   undefined,         &
                  &                   soil_texslo,       &
                  &                   dsmw_soil_unit,    &
                  &                   n_unit,            &
                  &                   dsmw_grid,         &
                  &                   lon_soil,          &
                  &                   lat_soil,          &
                  &                   soiltype_fao,      &
                  &                   soiltype_hwsd,     &
                  &                   fr_land_soil)

 
      
      print *,'MAXVAL(no_raw_data_pixel): ', MAXVAL(no_raw_data_pixel)
      print *,'MINVAL(no_raw_data_pixel): ', MINVAL(no_raw_data_pixel)

      print *,'MAXVAL(cosmo_soiltyp): ', MAXVAL(soiltype_fao)
      print *,'MINVAL(cosmo_soiltyp): ', MINVAL(soiltype_fao)

      print *,'MAXVAL(fr_land_soil): ', MAXVAL(fr_land_soil)
      print *,'MINVAL(fr_land_soil): ', MINVAL(fr_land_soil)



      print *,'Start filling of undefined target grid elements with nearest grid point raw data'

      CALL nearest_soil_data_to_target_grid(tg,         &
                  &                   undefined,        &
                  &                   soil_texslo,      &
                  &                   dsmw_soil_unit,   &
                  &                   n_unit,           &
                  &                   dsmw_grid,        &
                  &                   lon_soil,         &
                  &                   lat_soil,         & 
                  &                   soiltype_fao,     &
                  &                   soiltype_hwsd,    &
                  &                   fr_land_soil)


      print *,'Filling of undefined target grid elements with nearest grid point raw data done.'


      print *,'MAXVAL(no_raw_data_pixel): ', MAXVAL(no_raw_data_pixel)
      print *,'MINVAL(no_raw_data_pixel): ', MINVAL(no_raw_data_pixel)

      print *,'MAXVAL(cosmo_soiltyp): ', MAXVAL(soiltype_fao)
      print *,'MINVAL(cosmo_soiltyp): ', MINVAL(soiltype_fao)

      print *,'MAXVAL(fr_land_soil): ', MAXVAL(fr_land_soil)
      print *,'MINVAL(fr_land_soil): ', MINVAL(fr_land_soil)


      DEALLOCATE (dsmw_soil_unit, STAT = errorcode)
      IF (errorcode /= 0) print*, 'Cant deallocate dsmw_soil_unit'
      DEALLOCATE (soil_texslo, STAT = errorcode)
      IF (errorcode /= 0) print*, 'Cant deallocate soil_texslo'

      IF (ldeep_soil) THEN
        CALL allocate_raw_deep_soil_fields(nlon_soil, nlat_soil, n_unit)
        CALL get_deep_soil_data(path_deep_soil_file,start)

        print *,'HWSD deep soil read from file ', TRIM(path_deep_soil_file)

        print *,'aggregate deep_soil data to target grid'

        CALL agg_soil_data_to_target_grid(tg,             &
                  &                   undefined,          &
                  &                   soil_texslo_deep,   &
                  &                   dsmw_deep_soil_unit,&
                  &                   n_unit,             &
                  &                   dsmw_grid,          &
                  &                   lon_soil,           &
                  &                   lat_soil,           &
                  &                   soiltype_deep,      &
                  &                   soiltype_hwsd_s,       &
                  &                   fr_land_soil)

        print *,'MAXVAL(no_raw_data_pixel): ', MAXVAL(no_raw_data_pixel)
        print *,'MINVAL(no_raw_data_pixel): ', MINVAL(no_raw_data_pixel)

        print *,'MAXVAL(cosmo_deep_soiltyp): ', MAXVAL(soiltype_deep)
        print *,'MINVAL(cosmo_deep_soiltyp): ', MINVAL(soiltype_deep)

        print *,'Start filling of undefined deep soil target grid elements with nearest grid point raw data'

        CALL nearest_soil_data_to_target_grid(tg,         &
                  &                   undefined,          &
                  &                   soil_texslo_deep,   &
                  &                   dsmw_deep_soil_unit,&
                  &                   n_unit,             &
                  &                   dsmw_grid,          &
                  &                   lon_soil,           &
                  &                   lat_soil,           &   
                  &                   soiltype_deep,      &
                  &                   soiltype_hwsd_s,      & 
                  &                   fr_land_soil)

        print *,'Filling of undefined deep soil target grid elements with nearest grid point raw data done.'


        print *,'MAXVAL(no_raw_data_pixel): ', MAXVAL(no_raw_data_pixel)
        print *,'MINVAL(no_raw_data_pixel): ', MINVAL(no_raw_data_pixel)

        print *,'MAXVAL(cosmo_deep_soiltyp): ', MAXVAL(soiltype_deep)
        print *,'MINVAL(cosmo_deep_soiltyp): ', MINVAL(soiltype_deep)

        print *,'MAXVAL(cosmo_deep_soiltyp HWSD): ', MAXVAL(soiltype_hwsd_s)
        print *,'MINVAL(cosmo_deep_soiltyp HWSD): ', MINVAL(soiltype_hwsd_s)

        DEALLOCATE (dsmw_deep_soil_unit, STAT = errorcode)
        IF (errorcode /= 0) print*, 'Cant deallocate dsmw_deep_soil_unit'
        DEALLOCATE (soil_texslo_deep, STAT = errorcode)
        IF (errorcode /= 0) print*, 'Cant deallocate soil_texslo_deep'
      ENDIF

      PRINT *,'Start buffer output'

      netcdf_filename=  TRIM(soil_buffer_file)

      undefined = -999.0
      undefined_integer= 999

!roa bug fix soiltype_deep>
      IF (ldeep_soil) THEN
        CALL write_netcdf_soil_buffer(netcdf_filename,   &
   &                                   tg,               &
   &                                   isoil_data,       &
   &                                   ldeep_soil,       &
   &                                   undefined,        &
   &                                   undefined_integer,&
   &                                   lon_geo,          &
   &                                   lat_geo,          &
   &                                   fr_land_soil,     &
   &                                   soiltype_fao,     &
   &                                   soiltype_hwsd,     &
   &                                   soiltype_fao_deep = soiltype_deep,&
   &                                   soiltype_hwsd_deep= soiltype_hwsd_s   )
      ELSE
        CALL write_netcdf_soil_buffer(netcdf_filename,   &
   &                                   tg,               &
   &                                   isoil_data,       &
   &                                   ldeep_soil,       &
   &                                   undefined,        &
   &                                   undefined_integer,&
   &                                   lon_geo,          &
   &                                   lat_geo,          &
   &                                   fr_land_soil,     &
   &                                   soiltype_fao,     &
   &                                   soiltype_hwsd     )
      ENDIF
!roa bug fix soiltype_deep<

      PRINT *,'buffer output done'
      PRINT *,'Start target grid output'


      SELECT CASE(tg%igrid_type)
       !-----------------------------------------------------------------
       CASE(igrid_icon) ! ICON GRID

         netcdf_filename= TRIM(soil_output_file)

         undefined = -999.0
         undefined_integer= -999

!roa bug fix soiltype_deep>
      IF (ldeep_soil) THEN
         CALL write_netcdf_soil_icon_grid(netcdf_filename,  &
   &                                     icon_grid,         &
   &                                     isoil_data,        &
   &                                     tg,                &
   &                                     ldeep_soil,        &
   &                                     undefined,         &
   &                                     undefined_integer, &
   &                                     lon_geo,           &
   &                                     lat_geo,           &
   &                                     fr_land_soil,      &
   &                                     soiltype_fao,      &
   &                                     soiltype_deep = soiltype_deep)
      ELSE
         CALL write_netcdf_soil_icon_grid(netcdf_filename,  &
   &                                     icon_grid,         &
   &                                     isoil_data,        &
   &                                     tg,                &
   &                                     ldeep_soil,        &
   &                                     undefined,         &
   &                                     undefined_integer, &
   &                                     lon_geo,           &
   &                                     lat_geo,           &
   &                                     fr_land_soil,      &
   &                                     soiltype_fao)
      ENDIF
!roa bug fix soiltype_deep<

       CASE(igrid_cosmo) ! COSMO grid

         netcdf_filename= TRIM(soil_output_file)

         undefined = -999.0
         undefined_integer= -999

!roa bug fix soiltype_deep>
      IF (ldeep_soil) THEN
          CALL write_netcdf_soil_cosmo_grid(netcdf_filename, &
   &                                     cosmo_grid,         &
   &                                     tg,                 &
   &                                     isoil_data,         &
   &                                     ldeep_soil,         &
   &                                     undefined,          &
   &                                     undefined_integer,  &
   &                                     lon_geo,            &
   &                                     lat_geo,            &
   &                                     fr_land_soil,       &
   &                                     soiltype_fao,       &
   &                                     soiltype_deep = soiltype_deep)
      ELSE
          CALL write_netcdf_soil_cosmo_grid(netcdf_filename, &
   &                                     cosmo_grid,         &
   &                                     tg,                 &
   &                                     isoil_data,         &
   &                                     ldeep_soil,         &
   &                                     undefined,          &
   &                                     undefined_integer,  &
   &                                     lon_geo,            &
   &                                     lat_geo,            &
   &                                     fr_land_soil,       &
   &                                     soiltype_fao)
      ENDIF
!roa bug fix soiltype_deep<

       CASE(igrid_gme) ! GME grid   
    END SELECT


  PRINT *,'============= soil_to_buffer done ==============='

        

END PROGRAM extpar_soil_to_buffer
