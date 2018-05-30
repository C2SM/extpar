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

  ! Load the library information data:
  USE info_extpar, ONLY: info_define, info_print


  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4


  USE mo_target_grid_data, ONLY: lon_geo, &
    &                            lat_geo, &
    &                            no_raw_data_pixel

  USE mo_target_grid_data, ONLY: tg  !< structure with target grid description

  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_grid_structures, ONLY: igrid_icon

  USE mo_icon_grid_data, ONLY: icon_grid_region

  USE mo_io_units,          ONLY: filename_max

  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_sgsl_routines, ONLY: read_namelists_extpar_sg_slope

  USE mo_sgsl_tg_fields, ONLY:  sgsl

  USE mo_sgsl_tg_fields, ONLY:  allocate_sgsl_target_fields

  USE mo_sgsl_tg_fields, ONLY:  allocate_additional_sgsl_param

! mes > -------------------------------------------------------------
  USE mo_sgsl_data,      ONLY:  dem_aster,        &
    &                           dem_gl,           &
    &                           idem_type,        &    
    &                           sgsl_tiles_grid,   &
    &                           sgsl_grid,         &
    &                           ntiles,            &
    &                           max_tiles,         &
    &                           nc_tot,            &
    &                           nr_tot,            &
    &                           nc_tile,           &
    &                           tiles_lon_min,     &
    &                           tiles_lon_max,     &
    &                           tiles_lat_min,     &
    &                           tiles_lat_max,     &
    &                           demraw_lat_min,    &
    &                           demraw_lat_max,    &
    &                           demraw_lon_min,    &
    &                           demraw_lon_max,    &
    &                           num_tiles,         &
    &                           allocate_sgsl_data,&
    &                           fill_sgsl_data,    &
    &                           deallocate_sgsl_fields

! mes < -------------------------------------------------------------

  USE mo_sgsl_routines, ONLY:   det_sgsl_tiles_grid,           &
    &                           det_sgsl_grid

  USE mo_agg_sgsl, ONLY: agg_sgsl_data_to_target_grid

  USE mo_sgsl_output_nc, ONLY: write_netcdf_buffer_sgsl

  IMPLICIT NONE

  CHARACTER(len=filename_max) :: netcdf_filename
  CHARACTER(len=filename_max) :: namelist_grid_def

  CHARACTER (len=filename_max) :: namelist_sgsl_data_input !< file with input namelist with GLOBE data information
    
  CHARACTER (LEN=filename_max) :: sgsl_files(1:max_tiles)  !< filenames globe raw data

  CHARACTER (len=filename_max) :: sgsl_buffer_file !< name for orography buffer file
  CHARACTER (len=filename_max) :: raw_data_sgsl_path        !< path to raw data

  REAL(KIND=wp) :: undefined !< value to indicate undefined grid elements 
  INTEGER (KIND=i4) :: undef_int   !< value for undefined integer

  INTEGER (KIND=i4), ALLOCATABLE :: sgsl_startrow(:)    !< startrow indices for each DEM tile
  INTEGER (KIND=i4), ALLOCATABLE :: sgsl_endrow(:)      !< endrow indices for each DEM tile
  INTEGER (KIND=i4), ALLOCATABLE :: sgsl_startcolumn(:)  !< starcolumn indices for each DEM tile
  INTEGER (KIND=i4), ALLOCATABLE :: sgsl_endcolumn(:)   !< endcolumn indices for each DEM tile

  INTEGER :: k !< counter

  INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO

  ! variables for the ICON grid 
  INTEGER :: nvertex  !< total number of vertices

 !mes > -------------------------
 INTEGER (KIND=i4) :: ntiles_column        !< number of tile columns in total domain
 INTEGER (KIND=i4) :: ntiles_row           !< number of tile rows in total domain


  ALLOCATE (sgsl_startrow(1:ntiles), sgsl_endrow(1:ntiles),sgsl_startcolumn(1:ntiles),sgsl_endcolumn(1:ntiles))  !_br 21.02.14
!_br 21.02.14 for clean programming this should be deallocated somewhere

 ! Print the default information to stdout:
  CALL info_define ()
  CALL info_print ()
  !--------------------------------------------------------------------------------------------------------
 
  namelist_grid_def = 'INPUT_grid_org'

  CALL init_target_grid(namelist_grid_def)

  PRINT *,' target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
  igrid_type = tg%igrid_type
  ! get information on target grid
  ! get GLOBE raw data information

  !--------------------------------------------------------------------------------------------------------

  ! read namelist with globe data information

  namelist_sgsl_data_input = 'INPUT_SGSL'
  CALL read_namelists_extpar_sg_slope(namelist_sgsl_data_input, &
    &                                  raw_data_sgsl_path,   &
    &                                  sgsl_files,                &  !mes>
    &                                  ntiles_column,             &
    &                                  ntiles_row,                &
    &                                  idem_type,                &
    &                                  sgsl_buffer_file)

  CALL num_tiles(idem_type,ntiles_column, ntiles_row,ntiles,idem_type)        
 ! gives back the number of tiles that are available 16 for GLOBE or 36 for ASTER
  
!mes <

! mes >
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
      PRINT*,'edges of raw data domain: ', demraw_lon_min,' ', demraw_lon_max,' ', demraw_lat_min,' ',demraw_lat_max
    END SELECT
! mes <

!mes >
  PRINT*,' idem_type :', idem_type

  SELECT CASE (idem_type)
  CASE (dem_aster, dem_gl)
   IF (lon_geo (tg%ie,tg%je,tg%ke) > demraw_lon_max .OR. lon_geo(1,1,1) < demraw_lon_min) THEN
   PRINT*, 'raw data min lon is: ', demraw_lon_min, ' and raw data max lon is: ', demraw_lon_max
   CALL abort_extpar('The chosen longitude edges are not within the ASTER domain.')
   END IF
   IF (lat_geo(tg%ie,tg%je,tg%ke) > demraw_lat_max .OR. lat_geo(1,1,1) < demraw_lat_min) THEN
   PRINT*, 'raw data min lat is: ', demraw_lat_min, ' and raw data max lat is: ', demraw_lat_max
   CALL abort_extpar('The chosen latitude edges are not within the raw data domain.')
   END IF
  END SELECT


  CALL det_sgsl_tiles_grid(sgsl_tiles_grid)
  !HA debug
  DO k=1,ntiles
    print *,'sgsl files: ', TRIM(sgsl_files(k))
  ENDDO
    print *,'sgsl_tiles_grid(1): ', sgsl_tiles_grid(1)

  CALL det_sgsl_grid(sgsl_grid)
  !HA debug
  print *,'sgsl_grid: ', sgsl_grid


  ! allocate globe fields for target grid
  !--------------------------------------------------------------------------------------------------------
  CALL allocate_sgsl_target_fields(tg)

  ! allocate additional fields for icon grid
  !--------------------------------------------------------------------------------------------------------

  SELECT CASE(igrid_type)
    CASE(igrid_icon) ! ICON GRID
    ! allocate addtional target fields
    nvertex = icon_grid_region%nverts
    CALL  allocate_additional_sgsl_param(nvertex)
  END SELECT

  ! call the aggregation routine
  !--------------------------------------------------------------------------------------------------------
  PRINT *,'CALL agg_sgsl_data_to_target_grid'
  CALL agg_sgsl_data_to_target_grid(sgsl_tiles_grid, &
   &                                sgsl_grid,        &
   &                                tg,               &
   &                                sgsl_files,       &
   &                                sgsl,         &
   &                                no_raw_data_pixel,    &
   &                                raw_data_sgsl_path=raw_data_sgsl_path) !_br 17.09.14)

   ! if the target domain has a higher resolution of than the GLOBE data set (30'') some grid elements might not
   ! be set by the routine agg_sgsl_data_to_target_grid, (no_raw_data_pixel(ie,je,ke) == 0 in this case
   ! loop over all grid elements to check and perform a bilinear interplation if necessary
   k = 0
   undefined = -999.9

   PRINT *,'Max number of sgsl raw data pixel in a target grid element: '
   PRINT *,'MAXVAL(no_raw_data_pixel): ', MAXVAL(no_raw_data_pixel)

   PRINT *,'Min number of sgsl raw data pixel in a target grid element: '
   PRINT *,'MINVAL(no_raw_data_pixel): ', MINVAL(no_raw_data_pixel)

   PRINT *,'agg_sgsl_data_to_target_grid finished'

   ! output to netcdf file
   undefined = -999.9
   undef_int = -999

   netcdf_filename = TRIM(sgsl_buffer_file)
   print *, 'filename: ',TRIM(netcdf_filename)

   SELECT CASE(igrid_type)
     CASE(igrid_icon) ! ICON GRID
       CALL write_netcdf_buffer_sgsl(netcdf_filename,  &
        &                                tg,            &
        &                                undefined,     &
        &                                undef_int,     &
        &                                igrid_type,    &
        &                                lon_geo,       &
        &                                lat_geo,       &
        &                                sgsl)

     CASE DEFAULT

       CALL write_netcdf_buffer_sgsl(netcdf_filename,     &
        &                                tg,              &
        &                                undefined,       &
        &                                undef_int,       &
        &                                igrid_type,      &
        &                                lon_geo,         &
        &                                lat_geo,         &
        &                                sgsl)

   END SELECT

   CALL deallocate_sgsl_fields()

  PRINT *,'============= sgsl_to_buffer done ==============='


END PROGRAM extpar_sgsl_to_buffer

