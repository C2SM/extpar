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

  USE mo_logging
  USE info_extpar,              ONLY: info_print
  USE mo_io_units,              ONLY: filename_max
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: igrid_icon, &
       &                              igrid_cosmo

  USE mo_target_grid_routines,  ONLY: init_target_grid

  USE mo_target_grid_data,      ONLY: tg, &
       &                              lon_geo, &
       &                              lat_geo

  USE mo_icon_grid_data,        ONLY: ICON_grid

  USE mo_cosmo_grid,            ONLY: COSMO_grid

  USE mo_cru_data,              ONLY: read_namelists_extpar_t_clim, &
       &                              allocate_cru_data, &
       &                              deallocate_cru_data, &
       &                              get_dimension_cru_data, &
       &                              get_cru_grid_and_data

  USE mo_cru_target_fields,     ONLY: allocate_cru_target_fields,&
       &                              crutemp, &
       &                              cruelev,      &
       &                              i_t_cru_fine, &
       &                              i_t_cru_coarse

  USE mo_agg_cru,               ONLY: agg_cru_data_to_target_grid

  USE mo_cru_output_nc,         ONLY: write_netcdf_buffer_cru, &
       &                              write_netcdf_cosmo_grid_cru, &
       &                              write_netcdf_icon_grid_cru  

  IMPLICIT NONE

  CHARACTER(len=filename_max)  :: filename, &
       &                          netcdf_filename, &
       &                          namelist_file, &
       &                          namelist_grid_def, &
       &                          raw_data_t_clim_path, &     !< path to raw data
       &                          raw_data_t_clim_filename, & !< filename temperature climatology raw data
       &                          t_clim_buffer_file = '' , &!< name for temperature climatology buffer
       &                          t_clim_output_file = '' !< name for temperature climatology output file

  INTEGER (KIND=i4)           :: raw_data_t_id, & !< integer switch to choose a land use raw data set
       &                         nrows, & !< number of rows
       &                         ncolumns, & !< number of columns
       &                         ntime, & !< number of times
       &                         igrid_type  !< target grid type, 1 for ICON, 2 for COSMO

  REAL (KIND=wp)              :: undefined

  !--------------------------------------------------------------------------------------
  namelist_grid_def = 'INPUT_grid_org'
  namelist_file = 'INPUT_TCLIM'
  undefined = -999.0_wp

  CALL initialize_logging("extpar_cru_to_buffer.log")
  CALL info_print ()

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= start cru_to_buffer ==============')
  CALL logging%info( '')
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= init grid and read namelist=======')
  CALL logging%info( '')

  !--------------------------------------------------------------------------------------------------------
  ! get information on target grid, allocate target fields with coordinates and determin the coordinates 
  ! for th target grid

  CALL  init_target_grid(namelist_grid_def)

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

  ! inquire dimensions
  CALL  get_dimension_cru_data(filename,     &
       &                          nrows,        &
       &                          ncolumns,     &
       &                          ntime)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= allocate fields ==================')
  CALL logging%info( '')

  CALL allocate_cru_data(nrows,ncolumns,ntime)
  ! allocate target grid fields for aerosol optical thickness
  CALL allocate_cru_target_fields(tg)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info( '============= get grid and data ===============')
  CALL logging%info( '')

  ! read in aot raw data
  CALL get_cru_grid_and_data(filename,     &
       &                     raw_data_t_id,   &
       &                     nrows,        &
       &                     ncolumns)

  crutemp  = undefined  ! set target grid values to undefined
  cruelev  = undefined

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= start aggregation ================')
  CALL logging%info( '')

  CALL  agg_cru_data_to_target_grid(nrows,ncolumns,ntime,raw_data_t_id)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= write data to netcdf==============')
  CALL logging%info( '')

  !write out data
  filename = TRIM(t_clim_output_file)

  SELECT CASE(igrid_type)

    CASE(igrid_icon) ! ICON GRID

      IF (t_clim_output_file == '') THEN
        CALL logging%error('CRU ICON output filename not defined.',__FILE__,__LINE__)
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
    CALL logging%error('CRU buffer output filename not defined.',__FILE__,__LINE__)
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

  CALL logging%info( '')
  CALL logging%info('============= deallocate fields ==============')
  CALL logging%info( '')

  CALL deallocate_cru_data()

  CALL logging%info( '')
  CALL logging%info('============= cru_to_buffer done =============')
  
END PROGRAM extpar_cru_to_buffer
