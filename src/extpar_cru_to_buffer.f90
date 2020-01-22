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

  USE info_extpar, ONLY: info_print
  USE mo_logging
  
  USE mo_kind, ONLY: wp, i4, i4

  USE mo_grid_structures, ONLY: igrid_icon,          &
       &                        igrid_cosmo


  USE mo_target_grid_data, ONLY: tg,                &
       &                         lon_geo,           &
       &                         lat_geo

  USE mo_target_grid_routines, ONLY: init_target_grid
  USE mo_icon_grid_data, ONLY: ICON_grid !< structure which contains the definition of the ICON grid

  USE  mo_cosmo_grid, ONLY: COSMO_grid

  USE mo_io_units,          ONLY: filename_max

  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_cru_data, ONLY: read_namelists_extpar_t_clim, &
       &                 allocate_cru_data, &
       &                 deallocate_cru_data, &
       &                 get_dimension_cru_data, &
       &                 get_cru_grid_and_data, &
       &                 cru_grid

  USE mo_cru_target_fields, ONLY: allocate_cru_target_fields,&
       &                          crutemp, &
       &                          cruelev,      &
       &                          i_t_cru_fine, &
       &                          i_t_cru_coarse

  USE mo_agg_cru, ONLY: agg_cru_data_to_target_grid

  USE mo_cru_output_nc, ONLY: write_netcdf_buffer_cru, &
       &                         write_netcdf_cosmo_grid_cru, &
       &                         write_netcdf_icon_grid_cru  

  IMPLICIT NONE

  CHARACTER(len=filename_max) :: filename
  CHARACTER(len=filename_max) :: netcdf_filename

  CHARACTER(len=filename_max) :: namelist_file
  CHARACTER(len=filename_max) :: namelist_grid_def

  INTEGER (i4):: raw_data_t_id !< integer switch to choose a land use raw data set
  !! 1 CRU (fine), 2 CRU (coarse)

  CHARACTER (len=filename_max) :: raw_data_t_clim_path     !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_t_clim_filename !< filename temperature climatology raw data

  CHARACTER (len=filename_max) :: t_clim_buffer_file = '' !< name for temperature climatology buffer
  CHARACTER (len=filename_max) :: t_clim_output_file = '' !< name for temperature climatology output file

  REAL (wp) :: undefined

  !--------------------------------------------------------------------------------------

  INTEGER (i4) :: nrows !< number of rows
  INTEGER (i4) :: ncolumns !< number of columns
  INTEGER (i4) :: ntime !< number of times


  INTEGER (i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO

  ! local variables
  namelist_grid_def = 'INPUT_grid_org'
  namelist_file = 'INPUT_TCLIM'
  undefined = -999.0_wp

  CALL initialize_logging("extpar_cru_to_buffer.log")
  CALL info_print ()

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*) '============= start cru_to_buffer =============='
  WRITE(logging%fileunit,*) ''
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*) '============= init grid and read namelist======='
  WRITE(logging%fileunit,*) ''

  !--------------------------------------------------------------------------------------------------------
  ! get information on target grid, allocate target fields with coordinates and determin the coordinates 
  ! for th target grid

  CALL  init_target_grid(namelist_grid_def)

  IF (verbose >= idbg_low )THEN 
    WRITE(logging%fileunit,*)'target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
  ENDIF

  igrid_type = tg%igrid_type

  !------------------------------------------------------------------------------------

  ! get information about temperature climatology data

  CALL  read_namelists_extpar_t_clim(namelist_file, &
       raw_data_t_id, &
       raw_data_t_clim_path, &
       raw_data_t_clim_filename, &
       t_clim_buffer_file, &
       t_clim_output_file)

  filename = TRIM(raw_data_t_clim_path) // TRIM(raw_data_t_clim_filename)

  IF (verbose >= idbg_low ) WRITE(logging%fileunit,*) 'CRU raw data filename ',TRIM(filename)

  ! inquire dimensions
  CALL  get_dimension_cru_data(filename,     &
       &                          nrows,        &
       &                          ncolumns,     &
       &                          ntime)
  IF (verbose >= idbg_low ) THEN
    WRITE(logging%fileunit,*) 'nrows: ',nrows
    WRITE(logging%fileunit,*) 'ncolumns: ',ncolumns
    WRITE(logging%fileunit,*) 'ntime: ',ntime
  ENDIF

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= allocate fields =================='
  WRITE(logging%fileunit,*) ''

  CALL allocate_cru_data(nrows,ncolumns,ntime)
  ! allocate target grid fields for aerosol optical thickness
  CALL allocate_cru_target_fields(tg)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*) '============= get grid and data ==============='
  WRITE(logging%fileunit,*) ''

  ! read in aot raw data
  CALL get_cru_grid_and_data(filename,     &
       &                     raw_data_t_id,   &
       &                     nrows,        &
       &                     ncolumns)

  IF (verbose >= idbg_low ) THEN 
    WRITE(logging%fileunit,*) 'cru_grid: ', cru_grid
    WRITE(logging%fileunit,*) 'target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
  ENDIF


  crutemp  = undefined  ! set target grid values to undefined
  cruelev  = undefined

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= start aggregation ================'
  WRITE(logging%fileunit,*) ''

  undefined = -999.0_wp
  CALL  agg_cru_data_to_target_grid(nrows,ncolumns,ntime,raw_data_t_id)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= write data to netcdf=============='
  WRITE(logging%fileunit,*) ''

  !write out data
  filename = TRIM(t_clim_output_file)

  SELECT CASE(igrid_type)

  CASE(igrid_icon) ! ICON GRID

    IF (t_clim_output_file == '') THEN
      CALL abort_extpar('CRU ICON output filename not defined.')
    ELSE
      WRITE(logging%fileunit,*) 'Write CRU aggregated data for ICON to '//TRIM(t_clim_output_file)
    ENDIF
    
    netcdf_filename = TRIM(t_clim_output_file)

    SELECT CASE (raw_data_t_id)

    CASE(i_t_cru_coarse)

      CALL write_netcdf_icon_grid_cru(netcdf_filename,  &
           &                          icon_grid,        &
           &                          tg,               &
           &                          undefined,        &
           &                          lon_geo,          &
           &                          lat_geo,          &
           &                          crutemp)

    CASE(i_t_cru_fine)

      CALL write_netcdf_icon_grid_cru(netcdf_filename,  &
           &                          icon_grid,        &
           &                          tg,               &
           &                          undefined,        &
           &                          lon_geo,          &
           &                          lat_geo,          &
           &                          crutemp,          &
           &                          cruelev=cruelev)

    END SELECT

  CASE(igrid_cosmo) ! COSMO grid

    netcdf_filename = TRIM(t_clim_output_file)
    IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'write out ', TRIM(netcdf_filename)

    SELECT CASE (raw_data_t_id)
    CASE(i_t_cru_coarse)
      CALL write_netcdf_cosmo_grid_cru(netcdf_filename,  &
           &                                     cosmo_grid,       &
           &                                     tg,         &
           &                                     undefined, &
           &                                     lon_geo,     &
           &                                     lat_geo, &
           &                                     crutemp)


    CASE(i_t_cru_fine)
      CALL write_netcdf_cosmo_grid_cru(netcdf_filename,  &
           &                                     cosmo_grid,       &
           &                                     tg,         &
           &                                     undefined, &
           &                                     lon_geo,     &
           &                                     lat_geo, &
           &                                     crutemp, &
           &                                     cruelev=cruelev)

    END SELECT

  END SELECT

  IF (t_clim_buffer_file == '') THEN
    CALL abort_extpar('CRU buffer output filename not defined.')
  ELSE
    WRITE(logging%fileunit,*) 'Write CRU aggregated data buffer to '//TRIM(t_clim_buffer_file)
  ENDIF

  netcdf_filename = TRIM(t_clim_buffer_file)

  SELECT CASE (raw_data_t_id)
  CASE(i_t_cru_coarse)
    CALL write_netcdf_buffer_cru(netcdf_filename,  &
         &                       tg,         &
         &                       undefined, &
         &                       lon_geo,     &
         &                       lat_geo, &
         &                       crutemp)

  CASE(i_t_cru_fine)
    CALL write_netcdf_buffer_cru(netcdf_filename,  &
         &                       tg,         &
         &                       undefined, &
         &                       lon_geo,     &
         &                       lat_geo, &
         &                       crutemp, &
         &                       cruelev=cruelev)

  END SELECT

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= deallocate fields =============='
  WRITE(logging%fileunit,*) ''

  CALL deallocate_cru_data()

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= cru_to_buffer done ============='
  
END PROGRAM extpar_cru_to_buffer


