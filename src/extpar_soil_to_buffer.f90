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

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4
  USE info_extpar,              ONLY: info_print
  USE mo_io_units,              ONLY: filename_max

  USE mo_target_grid_data,      ONLY: &
       &                              lon_geo, &
       &                              tg, &
       &                              lat_geo

  USE mo_target_grid_routines,  ONLY: init_target_grid

  USE mo_agg_soil,              ONLY: agg_soil_data_to_target_grid, &
       &                              nearest_soil_data_to_target_grid

  USE mo_soil_routines,         ONLY: read_namelists_extpar_soil, &
       &                              get_soil_data, &
       &                              get_deep_soil_data, &
       &                              nlon_soil, nlat_soil, &
       &                              get_dimension_soil_data

  USE mo_soil_data,             ONLY: allocate_raw_soil_fields, &
       &                              allocate_raw_deep_soil_fields, &
       &                              define_soiltype,    &
       &                              soil_texslo,        &
       &                              soil_texslo_deep,   & 
       &                              dsmw_soil_unit,     &
       &                              dsmw_deep_soil_unit,&
       &                              n_unit,             &
       &                              dsmw_grid,          &
       &                              FAO_data, HWSD_data, HWSD_map, &
       &                              soil_data, &
       &                              lon_soil,  &
       &                              lon_full, &
       &                              lat_full, &
       &                              lat_soil

  USE mo_soil_tg_fields,       ONLY:  fr_land_soil, &
       &                              soiltype_fao,soiltype_hwsd, soiltype_deep,soiltype_hwsd_s, &
       &                              allocate_soil_target_fields

  USE mo_soil_output_nc,       ONLY:  write_netcdf_soil_buffer

  IMPLICIT NONE

  CHARACTER(len=filename_max) :: netcdf_filename, &
       &                         namelist_soil_data_input, & !< file with input namelist with soil data information
       &                         path_soil_file, &      !< filename with path for soil raw data     
       &                         path_deep_soil_file, &      !< filename with path for soil raw data
       &                         soil_buffer_file, &  !< name for soil buffer file
       &                         soil_buffer_file_consistent, & !< name for soil buffer file after consistency check
       &                         soil_output_file_consistent, & !< name for soil output file after consistency check
       &                         raw_data_soil_path, &        !< path to raw data
       &                         raw_data_soil_filename, & !< filename soil raw data
       &                         raw_data_deep_soil_filename, & !< filename deep soil raw data
       &                         namelist_grid_def !< filename with namelists for grid settings for EXTPAR

  REAL(KIND=wp)               :: undefined !< value to indicate undefined grid elements in cosmo_ndvi_field

  INTEGER (KIND=i4)           :: undefined_integer, &   !< value for undefined integer
       &                         isoil_data, &  !< soil data, 1 for FAO raw data, 
       &                         undef_soiltype, &
       &                         default_soiltype, &
       &                         soiltype_ice, &
       &                         soiltype_water, &
       &                         nlon_full, &
       &                         nlat_full, &
       &                         lon_low, lat_low, lon_hig, lat_hig, &
       &                         errorcode, &
       &                         start(2)

  LOGICAL                     :: ldeep_soil            !< switch to decide weather the deep soil layer is desired or not

  undefined_integer = 0 ! set undefined to zero
  undefined = -99.0 ! undef vlaue

  path_deep_soil_file          = "" !default name
  namelist_grid_def            = 'INPUT_grid_org'
  namelist_soil_data_input     = 'INPUT_SOIL'

  CALL initialize_logging("extpar_soil_to_buffer.log")
  CALL info_print ()
  !--------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= start soil_to_buffer =============')
  CALL logging%info( '')

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= read namelist and init grid ======')
  CALL logging%info( '')


  !--------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
  ! get information on target grid, allocate target fields with coordinates and determin the coordinates 
  ! for th target grid
  
  CALL  init_target_grid(namelist_grid_def)
  
  !--------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
  
  
  ! get information on soil raw data
  !--------------------------------------------------------------------------------------------------------
    
  ! read namelist with soil data information (path, filename)
  CALL read_namelists_extpar_soil(namelist_soil_data_input,        &
                                       isoil_data,                 &
                                       ldeep_soil,                 &
                                       raw_data_soil_path,         &
                                       raw_data_soil_filename,     &
                                       raw_data_deep_soil_filename,&
                                       soil_buffer_file,           &
                                       soil_buffer_file_consistent,&
                                       soil_output_file_consistent)
  
   
  
  
  IF (ldeep_soil .AND. isoil_data /= HWSD_data) THEN
    ldeep_soil = .FALSE.
    CALL logging%warning( 'One can only use the deep soil if HWSD data is used => ldeep_soil is set to FALSE!')
  ENDIF

  IF(ldeep_soil) THEN 
    WRITE(message_text,*)'raw_data_deep_soil_filename: ', TRIM(raw_data_deep_soil_filename) 
    CALL logging%info(message_text)
  ENDIF

  path_soil_file = TRIM(raw_data_soil_path) // TRIM(raw_data_soil_filename)
  IF (ldeep_soil) THEN
    path_deep_soil_file = TRIM(raw_data_soil_path) // TRIM(raw_data_deep_soil_filename)
  ENDIF

  IF (ldeep_soil)THEN
    WRITE(message_text,*)'path_deep_soil_file: ', TRIM(path_deep_soil_file)
    CALL logging%info(message_text)
  ENDIF

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

! latitude runs from north to south in the raw data
  lat_low = MAX(1,INT((lat_full(1)-MAXVAL(lat_geo)-0.5_wp+dsmw_grid%dlat_reg)/ABS(dsmw_grid%dlat_reg)))
  lat_hig = MIN(nlat_full,INT((lat_full(1)-MINVAL(lat_geo)+0.5_wp-dsmw_grid%dlat_reg)/ABS(dsmw_grid%dlat_reg)))
  nlat_soil = lat_hig + 1 - lat_low
  start(1) = lon_low
  start(2) = lat_low

  ! get coordinates and legend and data from raw data file
  ! define value of global variables soil_raw_data_grid, lon_reg, lat_reg, soil_texslo, dsmw_soil_unit
  !--------------------------------------------------------------------------------------------------------
  CALL define_soiltype(isoil_data, ldeep_soil, &
                       undef_soiltype,         &
                       default_soiltype,       &
                       soiltype_ice,           &
                       soiltype_water,         &
                       soil_data)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info( '============= allocate fields ==================')
  CALL logging%info( '')

  CALL allocate_raw_soil_fields(nlon_soil, nlat_soil, n_unit)

  CALL get_soil_data(path_soil_file, start)
  lon_soil = lon_full(lon_low:lon_hig)
  lat_soil = lat_full(lat_low:lat_hig)

  CALL get_soil_data(path_soil_file,start)
  CALL allocate_soil_target_fields(tg, ldeep_soil, l_use_array_cache=.FALSE.)

  !--------------------------------------------------------------------------------------------------------

  SELECT CASE(isoil_data)
    CASE(FAO_data)
      WRITE(message_text,*)'FAO DSMW read from file ', TRIM(path_soil_file)
      CALL logging%info(message_text)
    CASE(HWSD_data, HWSD_map)
      WRITE(message_text,*)'HWSD read from file ', TRIM(path_soil_file)
      CALL logging%info(message_text)
    END SELECT


  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= start aggregation ================'
  WRITE(logging%fileunit,*) ''

  undefined = 0.0

  CALL agg_soil_data_to_target_grid(tg,              &
              &                   undefined,         &
              &                   soil_texslo,       &
              &                   dsmw_soil_unit,    &
              &                   dsmw_grid,         &
              &                   lon_soil,          &
              &                   lat_soil,          &
              &                   soiltype_fao,      &
              &                   soiltype_hwsd,     &
              &                   fr_land_soil)

 
  
  CALL logging%info('Fill undefined target grid elements with nearest grid point raw data')

  CALL nearest_soil_data_to_target_grid(tg,         &
              &                   undefined,        &
              &                   soil_texslo,      &
              &                   dsmw_soil_unit,   &
              &                   dsmw_grid,        &
              &                   soiltype_fao,     &
              &                   soiltype_hwsd,    &
              &                   fr_land_soil)


  DEALLOCATE (dsmw_soil_unit, STAT = errorcode)
  IF (errorcode /= 0) CALL logging%error('Cant deallocate dsmw_soil_unit',__FILE__,__LINE__)
  DEALLOCATE (soil_texslo, STAT = errorcode)
  IF (errorcode /= 0) CALL logging%error('Cant deallocate soil_texslo',__FILE__,__LINE__)

  IF (ldeep_soil) THEN
    CALL allocate_raw_deep_soil_fields(nlon_soil, nlat_soil, n_unit)
    CALL get_deep_soil_data(path_deep_soil_file,start)

    CALL agg_soil_data_to_target_grid(tg,             &
              &                   undefined,          &
              &                   soil_texslo_deep,   &
              &                   dsmw_deep_soil_unit,&
              &                   dsmw_grid,          &
              &                   lon_soil,           &
              &                   lat_soil,           &
              &                   soiltype_deep,      &
              &                   soiltype_hwsd_s,       &
              &                   fr_land_soil)

    CALL logging%info('Fill undefined target grid elements with nearest grid point raw data')

    CALL nearest_soil_data_to_target_grid(tg,         &
              &                   undefined,          &
              &                   soil_texslo_deep,   &
              &                   dsmw_deep_soil_unit,&
              &                   dsmw_grid,          &
              &                   soiltype_deep,      &
              &                   soiltype_hwsd_s,      & 
              &                   fr_land_soil)

    DEALLOCATE (dsmw_deep_soil_unit, STAT = errorcode)
    IF (errorcode /= 0) CALL logging%error('Cant deallocate dsmw_deep_soil_unit',__FILE__,__LINE__)
    DEALLOCATE (soil_texslo_deep, STAT = errorcode)
    IF (errorcode /= 0) CALL logging%error('Cant deallocate soil_texslo_deep',__FILE__,__LINE__)
  ENDIF

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= write data to netcdf==============')
  CALL logging%info( '')

  netcdf_filename=  TRIM(soil_buffer_file)

  undefined = -999.0
  undefined_integer= 999

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

  CALL logging%info( '')
  CALL logging%info('============= soil_to_buffer done ===============')

END PROGRAM extpar_soil_to_buffer
