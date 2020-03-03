!+ Fortran main program to read in DEM slope data and aggregate to target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_0         2016/07/28 Daniel Luethi
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran main program to read in DEM slope data and aggregate to target grid
!>  
!! @par extpar_sgsl_to_buffer 
!!
!! This program reads in the GLOBE/ASTER slope data set and aggregates it to the target grid 
!!
!> Purpose: read in GLOBE/ASTER slope data and aggregate to COSMO/ICON grid
!> \author Daniel Luethi
PROGRAM extpar_sgsl_to_buffer

  USE mo_logging
  USE info_extpar,              ONLY: info_print

  USE mo_kind,                  ONLY: wp, i4

  USE mo_target_grid_data,      ONLY: lon_geo, &
       &                              lat_geo, &
       &                              no_raw_data_pixel

  USE mo_target_grid_data,      ONLY: tg  !< structure with target grid description

  USE mo_target_grid_routines,  ONLY: init_target_grid

  USE mo_grid_structures,       ONLY: igrid_icon

  USE mo_icon_grid_data,        ONLY: icon_grid_region

  USE mo_io_units,              ONLY: filename_max

  USE mo_sgsl_routines,         ONLY: read_namelists_extpar_sg_slope, &
       &                              det_sgsl_grid, &
       &                              det_sgsl_tiles_grid

  USE mo_sgsl_tg_fields,        ONLY: sgsl, &
       &                              allocate_sgsl_target_fields, &
       &                              allocate_additional_sgsl_param
                                
  USE mo_sgsl_data,             ONLY: dem_aster,        &
       &                              dem_gl,           &
       &                              idem_type,        &    
       &                              sgsl_tiles_grid,   &
       &                              sgsl_grid,         &
       &                              ntiles,            &
       &                              max_tiles,         &
       &                              nc_tot,            &
       &                              nr_tot,            &
       &                              nc_tile,           &
       &                              tiles_lon_min,     &
       &                              tiles_lon_max,     &
       &                              tiles_lat_min,     &
       &                              tiles_lat_max,     &
       &                              demraw_lat_min,    &
       &                              demraw_lat_max,    &
       &                              demraw_lon_min,    &
       &                              demraw_lon_max,    &
       &                              num_tiles,         &
       &                              allocate_sgsl_data,&
       &                              fill_sgsl_data,    &
       &                              deallocate_sgsl_fields


  USE mo_agg_sgsl,              ONLY: agg_sgsl_data_to_target_grid

  USE mo_sgsl_output_nc,        ONLY: write_netcdf_buffer_sgsl
  USE mo_preproc_for_sgsl,       ONLY: prepare_preproc

  IMPLICIT NONE

  CHARACTER(len=filename_max)    :: netcdf_filename, &
       &                            namelist_preproc, &
       &                            namelist_grid_def, &
       &                            namelist_sgsl_data_input, & !< file with input namelist with GLOBE data information
       &                            sgsl_files(1:max_tiles), &  !< filenames globe raw data
       &                            sgsl_buffer_file, & !< name for orography buffer file
       &                            raw_data_sgsl_path        !< path to raw data

  REAL(KIND=wp)                  :: undefined !< value to indicate undefined grid elements 
                               
  INTEGER (KIND=i4)              :: k, nvertex, &
       &                            ntiles_column, &        !< number of tile columns in total domain
       &                            ntiles_row, &           !< number of tile rows in total domain
       &                            igrid_type  !< target grid type, 1 for ICON, 2 for COSMO

  INTEGER (KIND=i4), ALLOCATABLE :: sgsl_startrow(:), &    !< startrow indices for each DEM tile
       &                            sgsl_endrow(:), &      !< endrow indices for each DEM tile
       &                            sgsl_startcolumn(:), &  !< starcolumn indices for each DEM tile
       &                            sgsl_endcolumn(:)   !< endcolumn indices for each DEM tile



  namelist_preproc         = 'INPUT_ORO'
  namelist_grid_def        = 'INPUT_grid_org'
  namelist_sgsl_data_input = 'INPUT_SGSL'

  ALLOCATE (sgsl_startrow(1:ntiles), sgsl_endrow(1:ntiles),sgsl_startcolumn(1:ntiles),sgsl_endcolumn(1:ntiles))

  CALL initialize_logging("extpar_sgsl_to_buffer.log")
  CALL info_print ()

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= start sgsl_to_buffer =============')
  CALL logging%info( '')

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= read namelist and init grid ======')
  CALL logging%info( '')
 

  CALL init_target_grid(namelist_grid_def)

  igrid_type = tg%igrid_type
  ! get information on target grid
  ! get GLOBE raw data information

  ! read namelist with globe data information
  CALL read_namelists_extpar_sg_slope(namelist_sgsl_data_input, &
    &                                  raw_data_sgsl_path,   &
    &                                  sgsl_files,                &  !mes>
    &                                  ntiles_column,             &
    &                                  ntiles_row,                &
    &                                  idem_type,                &
    &                                  sgsl_buffer_file)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= preprocess raw oro data === ======')
  CALL logging%info( '')

  CALL prepare_preproc(namelist_preproc, &
       &               sgsl_files)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= read namelist and init grid ======')
  CALL logging%info( '')

  CALL num_tiles(ntiles_column, ntiles_row,ntiles)        
   ! gives back the number of tiles that are available 16 for GLOBE or 36 for ASTER
  
  CALL allocate_sgsl_data(ntiles)                  ! allocates the data using ntiles

  CALL fill_sgsl_data(raw_data_sgsl_path,sgsl_files, &! the allocated vectors need to be filled with the respective value.
                                           tiles_lon_min, &
                                           tiles_lon_max, &    
                                           tiles_lat_min, &
                                           tiles_lat_max, &
                                           nc_tot,        &
                                           nr_tot,        &
                                           nc_tile)

  SELECT CASE(idem_type)
    CASE(dem_aster, dem_gl)
      WRITE(message_text,*)'Edges of raw data domain: ', demraw_lon_min,' ', demraw_lon_max,' ', demraw_lat_min,' ',demraw_lat_max
      CALL logging%info(message_text)
  END SELECT

  SELECT CASE (idem_type)
    CASE (dem_aster, dem_gl)
      IF (lon_geo (tg%ie,tg%je,tg%ke) > demraw_lon_max .OR. lon_geo(1,1,1) < demraw_lon_min) THEN
        WRITE(message_text,*) 'raw data min lon is: ', demraw_lon_min, ' and raw data max lon is: ', demraw_lon_max
        CALL logging%warning(message_text)
        CALL logging%error('The chosen longitude edges are not within the ASTER domain.',__FILE__,__LINE__)
      END IF
      IF (lat_geo(tg%ie,tg%je,tg%ke) > demraw_lat_max .OR. lat_geo(1,1,1) < demraw_lat_min) THEN
        WRITE(message_text,*) 'raw data min lat is: ', demraw_lat_min, ' and raw data max lat is: ', demraw_lat_max
        CALL logging%warning(message_text)
        CALL logging%error('The chosen latitude edges are not within the raw data domain.',__FILE__,__LINE__)
      END IF
  END SELECT

  CALL det_sgsl_tiles_grid(sgsl_tiles_grid)
  DO k=1,ntiles
    WRITE(message_text,*)'sgsl files: ', TRIM(sgsl_files(k))
    CALL logging%info(message_text)
  ENDDO

  CALL det_sgsl_grid(sgsl_grid)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info( '============= allocate fields ==================')
  CALL logging%info( '')

  CALL allocate_sgsl_target_fields(tg)

  ! allocate additional fields for icon grid
  !--------------------------------------------------------------------------------------------------------

  SELECT CASE(igrid_type)
    CASE(igrid_icon) ! ICON GRID
      ! allocate addtional target fields
      nvertex = icon_grid_region%nverts
      CALL  allocate_additional_sgsl_param(nvertex)
  END SELECT

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= start aggregation ================')
  CALL logging%info( '')

  CALL agg_sgsl_data_to_target_grid(sgsl_tiles_grid, &
  &                                sgsl_grid,        &
  &                                tg,               &
  &                                sgsl_files,       &
  &                                sgsl,         &
  &                                no_raw_data_pixel,    &
  &                                raw_data_sgsl_path=raw_data_sgsl_path)

  ! if the target domain has a higher resolution of than the GLOBE data set (30'') some grid elements might not
  ! be set by the routine agg_sgsl_data_to_target_grid, (no_raw_data_pixel(ie,je,ke) == 0 in this case
  ! loop over all grid elements to check and perform a bilinear interplation if necessary
  k = 0
  undefined = -999.9

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= write data to netcdf==============')
  CALL logging%info( '')

  undefined = -999.9

  netcdf_filename = TRIM(sgsl_buffer_file)

  SELECT CASE(igrid_type)
    CASE(igrid_icon) ! ICON GRID
      CALL write_netcdf_buffer_sgsl(netcdf_filename,  &
       &                                tg,            &
       &                                undefined,     &
       &                                igrid_type,    &
       &                                lon_geo,       &
       &                                lat_geo,       &
       &                                sgsl)

    CASE DEFAULT

      CALL write_netcdf_buffer_sgsl(netcdf_filename,     &
       &                                tg,              &
       &                                undefined,       &
       &                                igrid_type,      &
       &                                lon_geo,         &
       &                                lat_geo,         &
       &                                sgsl)

  END SELECT

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= deallocate fields =================')
  CALL logging%info( '')

   CALL deallocate_sgsl_fields()

  CALL logging%info( '')
  CALL logging%info('============= sgsl_to_buffer done ===============')

END PROGRAM extpar_sgsl_to_buffer
