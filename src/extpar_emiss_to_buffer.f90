!+  Fortran main program to read in EMISS data and aggregate to target grid
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
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran main program to read in EMISS data and aggregate to target grid
!>  
!> \author Hermann Asensio
PROGRAM extpar_emiss_to_buffer

  USE mo_logging
  USE info_extpar,              ONLY: info_print
  USE mo_kind,                  ONLY: wp,i4
  USE mo_io_units,              ONLY: filename_max

  USE mo_target_grid_data,      ONLY: lon_geo, &
       &                              lat_geo, &
       &                              tg

  USE mo_target_grid_routines,  ONLY: init_target_grid

  USE mo_grid_structures,       ONLY: igrid_icon, &
       &                              igrid_cosmo

  USE mo_icon_grid_data,        ONLY: ICON_grid !< structure which contains the definition of the ICON grid
 
  USE mo_cosmo_grid,            ONLY: COSMO_grid
    


  USE mo_emiss_data,            ONLY: emiss_raw_data_grid, &
       &                              lon_emiss, &
       &                              lat_emiss, &
       &                              ntime_emiss, &
       &                              allocate_raw_emiss_fields,&
       &                              deallocate_emiss_fields
                                  
  USE mo_emiss_tg_fields,       ONLY:  emiss_max, &
       &                               emiss_field_mom, &
       &                               emiss_ratio_mom, &
       &                               allocate_emiss_target_fields

  USE mo_emiss_routines,        ONLY: open_netcdf_EMISS_data, &
       &                              close_netcdf_EMISS_data, &
       &                              get_dimension_EMISS_data, &
       &                              read_namelists_extpar_emiss, &
       &                              get_EMISS_data_coordinates

  USE mo_agg_emiss,             ONLY: agg_emiss_data_to_target_grid

  USE mo_emiss_output_nc,       ONLY: write_netcdf_buffer_emiss, &
       &                              write_netcdf_cosmo_grid_emiss, &
       &                              write_netcdf_icon_grid_emiss

  IMPLICIT NONE

  CHARACTER(len=filename_max)     :: namelist_grid_def, &
       &                             namelist_emiss_data_input, & !< file with input namelist with EMISS data information
       &                             raw_data_emiss_filename, & !< filename emiss raw data
       &                             path_emiss_file, &      !< filename with path for EMISS raw data
       &                             netcdf_filename, &      !< filename for netcdf file with EMISS data on COSMO grid
       &                             raw_data_emiss_path, &        !< path to raw data
       &                             emiss_buffer_file, & !< name for EMISS buffer file
       &                             emiss_output_file !< name for EMISS output file

  INTEGER (KIND=i4)               :: ncid_emiss, &  !< netcdf unit file number for EMISS data netcdf file
       &                             nmonth, &  !< index for month for EMISS data
       &                             nlon_emiss, & !< number of grid elements in zonal direction for EMISS data
       &                             nlat_emiss, & !< number of grid elements in meridional direction for EMISS data
       &                             igrid_type  !< target grid type, 1 for ICON, 2 for COSMO

  INTEGER  (KIND=i4), ALLOCATABLE :: time(:) !< array with time axis values (needed for netcdf cf convention)
  
  REAL (KIND=wp)                  :: dlon_emiss, & !< grid point distance in zonal direction (in degrees) for EMISS data
       &                             dlat_emiss, & !< grid point distance in meridional direction (in degrees) for EMISS data
       &                             startlon_emiss, & !< longitude of lower left grid element for EMISS data 
       &                             startlat_emiss, & !< latitude of lower left grid element for EMISS data
       &                             undefined !< value to indicate undefined grid elements 


  namelist_emiss_data_input = 'INPUT_EMISS'
  namelist_grid_def         = 'INPUT_grid_org'

  CALL initialize_logging("extpar_emiss_to_buffer.log")
  CALL info_print ()

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= start emiss_to_buffer ============')
  CALL logging%info( '')
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= init grid and read namelist=======')
  CALL logging%info( '')

  undefined = -999.0 ! undef vlaue
      
  CALL init_target_grid(namelist_grid_def)

  igrid_type = tg%igrid_type

  ! read namelist for input EMISS data
  CALL  read_namelists_extpar_emiss(namelist_emiss_data_input, &
    &                                  raw_data_emiss_path, &
    &                                  raw_data_emiss_filename, &
    &                                  emiss_buffer_file, &
    &                                  emiss_output_file)
     
  path_emiss_file = TRIM(raw_data_emiss_path)//TRIM(raw_data_emiss_filename)
       
  ! open netcdf file with EMISS data
  CALL open_netcdf_EMISS_data(path_emiss_file, &
    &                           ncid_emiss)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= allocate fields and get dimensions')
  CALL logging%info( '')

  !> inquire dimension information for EMISS raw data 
  CALL get_dimension_EMISS_data(ncid_emiss, &
   &                                nlon_emiss, &
   &                                nlat_emiss, &
                                   ntime_emiss)

  ALLOCATE(time(ntime_emiss)) ! this array is needed for netcdf output at the end
  DO nmonth=1, ntime_emiss
    time(nmonth) = nmonth
  ENDDO

  CALL allocate_raw_emiss_fields(nlon_emiss,nlat_emiss,ntime_emiss)   
  CALL allocate_emiss_target_fields(tg,ntime_emiss)

  CALL get_EMISS_data_coordinates(ncid_emiss,      &
    &                               nlon_emiss,      &
    &                               nlat_emiss,      &
    &                               startlon_emiss,  &
    &                               startlat_emiss,  &
    &                               dlon_emiss,      &
    &                               dlat_emiss,      &
    &                               lon_emiss,       &
    &                               lat_emiss)

  ! put the values of the grid definition in the data structure emiss_raw_data_grid (type emiss_reg_lonlat_grid)
  emiss_raw_data_grid%start_lon_reg= startlon_emiss
  emiss_raw_data_grid%start_lat_reg= startlat_emiss
  emiss_raw_data_grid%dlon_reg= dlon_emiss
  emiss_raw_data_grid%dlat_reg= dlat_emiss
  emiss_raw_data_grid%nlon_reg= nlon_emiss
  emiss_raw_data_grid%nlat_reg= nlat_emiss

  emiss_raw_data_grid%end_lon_reg= lon_emiss(nlon_emiss) ! startlon_emiss + (nlon_emiss - 1) * dlon_emiss
  emiss_raw_data_grid%end_lat_reg= lat_emiss(nlat_emiss) ! startlat_emiss - (nlat_emiss - 1) * dlat_emiss 
 ! not negative increment, but EMISS latitude goes from north to south

  CALL close_netcdf_EMISS_data(ncid_emiss)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= start aggregation ================')
  CALL logging%info( '')

  CALL agg_emiss_data_to_target_grid(tg, path_emiss_file)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= write data to netcdf==============')
  CALL logging%info( '')

  SELECT CASE(igrid_type)

    CASE(igrid_icon) ! ICON GRID

      netcdf_filename = TRIM(emiss_output_file)
      undefined = -500.

      CALL write_netcdf_icon_grid_emiss(netcdf_filename,  &
           &                             icon_grid,         &
           &                             tg,         &
           &                             ntime_emiss, &
           &                             undefined, &
           &                             lon_geo,     &
           &                             lat_geo, &
           &                             emiss_max,  &
           &                             emiss_field_mom,&
           &                             emiss_ratio_mom)

    CASE(igrid_cosmo) ! COSMO grid
    
      netcdf_filename = TRIM(emiss_output_file)
      undefined = -500.


      CALL write_netcdf_cosmo_grid_emiss(netcdf_filename,  &
           &                             cosmo_grid,         &
           &                             tg,         &
           &                             ntime_emiss, &
           &                             undefined, &
           &                             emiss_max,  &
           &                             emiss_field_mom,&
           &                             emiss_ratio_mom)

  END SELECT

  netcdf_filename = TRIM(emiss_buffer_file)
  undefined = -500.

  CALL write_netcdf_buffer_emiss(netcdf_filename,  &
       &                                 tg,         &
       &                                 ntime_emiss, &
       &                                 undefined, &
       &                                 lon_geo,     &
       &                                 lat_geo, &
       &                                 emiss_max,  &
       &                                 emiss_field_mom,&
       &                                 emiss_ratio_mom)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= deallocate fields =================')
  CALL logging%info( '')

  CALL deallocate_emiss_fields()

  CALL logging%info( '')
  CALL logging%info( '============= emiss_to_buffer done =============')

END PROGRAM extpar_emiss_to_buffer
