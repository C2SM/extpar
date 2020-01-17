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

  USE info_extpar, ONLY: info_print
  USE mo_logging

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4


  USE mo_target_grid_data, ONLY: lon_geo, &
    &                            lat_geo

  USE mo_target_grid_data, ONLY: tg  !< structure with target grid description

  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo

  USE  mo_icon_grid_data, ONLY: ICON_grid !< structure which contains the definition of the ICON grid
 
  USE  mo_cosmo_grid, ONLY: COSMO_grid
    
  USE mo_io_units,          ONLY: filename_max

  USE mo_emiss_routines, ONLY: read_namelists_extpar_emiss

  USE mo_emiss_data, ONLY: emiss_raw_data_grid, &
    &                           lon_emiss, &
    &                           lat_emiss, &
    &                           ntime_emiss, &
    &                           allocate_raw_emiss_fields,&
    &                           deallocate_emiss_fields
                               
  USE mo_emiss_tg_fields, ONLY:  emiss_max, &
    &                                emiss_field_mom, &
    &                                emiss_ratio_mom, &
    &                                allocate_emiss_target_fields

  USE mo_emiss_routines, ONLY: open_netcdf_EMISS_data, &
    &                               close_netcdf_EMISS_data, &
    &                               get_dimension_EMISS_data, &
    &                               get_EMISS_data_coordinates
                                   

  USE mo_agg_emiss, ONLY: agg_emiss_data_to_target_grid

  USE mo_emiss_output_nc, ONLY: write_netcdf_buffer_emiss
  USE mo_emiss_output_nc, ONLY: write_netcdf_cosmo_grid_emiss
  USE mo_emiss_output_nc, ONLY: write_netcdf_icon_grid_emiss

  IMPLICIT NONE




  CHARACTER(len=filename_max) :: namelist_grid_def

  CHARACTER (len=filename_max) :: namelist_emiss_data_input !< file with input namelist with EMISS data information

  CHARACTER (len=filename_max) :: raw_data_emiss_filename !< filename emiss raw data
  CHARACTER (len=filename_max) :: path_emiss_file      !< filename with path for EMISS raw data
  CHARACTER (len=filename_max) :: netcdf_filename      !< filename for netcdf file with EMISS data on COSMO grid
  CHARACTER (len=filename_max) :: raw_data_emiss_path        !< path to raw data

  CHARACTER (len=filename_max) :: emiss_buffer_file !< name for EMISS buffer file
  CHARACTER (len=filename_max) :: emiss_output_file !< name for EMISS output file


  INTEGER (KIND=i4) :: ncid_emiss  !< netcdf unit file number for EMISS data netcdf file

  INTEGER  (KIND=i4) :: nlon_emiss !< number of grid elements in zonal direction for EMISS data
  INTEGER  (KIND=i4) :: nlat_emiss !< number of grid elements in meridional direction for EMISS data

  INTEGER  (KIND=i4), ALLOCATABLE :: time(:) !< array with time axis values (needed for netcdf cf convention)

  INTEGER (KIND=i4):: nmonth  !< index for month for EMISS data

  
  REAL (KIND=wp) :: dlon_emiss !< grid point distance in zonal direction (in degrees) for EMISS data
  REAL (KIND=wp) :: dlat_emiss !< grid point distance in meridional direction (in degrees) for EMISS data

  REAL (KIND=wp) :: startlon_emiss !< longitude of lower left grid element for EMISS data 

  REAL (KIND=wp) :: startlat_emiss !< latitude of lower left grid element for EMISS data

  INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO

  REAL(KIND=wp) :: undefined !< value to indicate undefined grid elements 

  CALL initialize_logging("extpar_emiss_to_buffer.log", stdout_level=debug)
  CALL info_print ()
  !--------------------------------------------------------------------------------------------------------
  undefined = -999.0 ! undef vlaue
      
  namelist_grid_def = 'INPUT_grid_org'
  CALL init_target_grid(namelist_grid_def)

  PRINT *,'target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
  !HA debug:
    print *,' MAXVAL(lat_geo): ', MAXVAL(lat_geo)
    print *,' MINVAL(lat_geo): ', MINVAL(lat_geo)


  igrid_type = tg%igrid_type
  ! get information on target grid

  ! read namelist for input EMISS data

  namelist_emiss_data_input = 'INPUT_EMISS'
  CALL  read_namelists_extpar_emiss(namelist_emiss_data_input, &
    &                                  raw_data_emiss_path, &
    &                                  raw_data_emiss_filename, &
    &                                  emiss_buffer_file, &
    &                                  emiss_output_file)
     
  path_emiss_file = TRIM(raw_data_emiss_path)//TRIM(raw_data_emiss_filename)
  !HA debug
  print *, 'after reading namelist for input EMISS data, EMISS raw data are in file:'
  print *, TRIM(path_emiss_file)

       
  ! open netcdf file with EMISS data
  CALL open_netcdf_EMISS_data(path_emiss_file, &
    &                           ncid_emiss)


   !> inquire dimension information for EMISS raw data 
   CALL get_dimension_EMISS_data(ncid_emiss, &
    &                                nlon_emiss, &
    &                                nlat_emiss, &
                                    ntime_emiss)
  !HA debug
  print *, 'after check of dimensions in EMISS raw data file'
  print *, 'nlon_emiss, nlat_emiss: ',nlon_emiss, nlat_emiss
  print *, 'ntime_emiss: ', ntime_emiss

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

  !HA debug
  print *, 'after getting EMISS data coordinates'
  print *,'startlon_emiss: ', startlon_emiss
  print *,'startlat_emiss: ', startlat_emiss
  print *,'dlon_emiss: ', dlon_emiss
  print *,'dlat_emiss: ', dlat_emiss
  print *,'lon_emiss(1) = ',lon_emiss(1) 
  print *,'lon_emiss(nlon_emiss) = ', lon_emiss(nlon_emiss) 
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
  print *,'emiss_raw_data_grid: ',emiss_raw_data_grid

  CALL close_netcdf_EMISS_data(ncid_emiss)

  ! start aggregation
  PRINT *,'aggregate EMISS data to target grid'

  CALL agg_emiss_data_to_target_grid(tg, path_emiss_file)

  PRINT *,'aggregation done'

  !write out data

  SELECT CASE(igrid_type)

    CASE(igrid_icon) ! ICON GRID

      netcdf_filename = TRIM(emiss_output_file)
      undefined = -500.

      PRINT *,'write out ', TRIM(netcdf_filename)

      CALL write_netcdf_icon_grid_emiss(netcdf_filename,  &
   &                                     icon_grid,         &
   &                                     tg,         &
   &                                     ntime_emiss, &
   &                                     undefined, &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     emiss_max,  &
   &                                     emiss_field_mom,&
   &                                     emiss_ratio_mom)





    CASE(igrid_cosmo) ! COSMO grid
    
      netcdf_filename = TRIM(emiss_output_file)
      undefined = -500.

      PRINT *,'write out ', TRIM(netcdf_filename)


      CALL write_netcdf_cosmo_grid_emiss(netcdf_filename,  &
   &                                     cosmo_grid,         &
   &                                     tg,         &
   &                                     ntime_emiss, &
   &                                     undefined, &
   &                                     emiss_max,  &
   &                                     emiss_field_mom,&
   &                                     emiss_ratio_mom)

  END SELECT

  netcdf_filename = TRIM(emiss_buffer_file)
  undefined = -500.

  PRINT *,'write out ', TRIM(netcdf_filename)


  CALL write_netcdf_buffer_emiss(netcdf_filename,  &
   &                                     tg,         &
   &                                     ntime_emiss, &
   &                                     undefined, &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     emiss_max,  &
   &                                     emiss_field_mom,&
   &                                     emiss_ratio_mom)

  CALL deallocate_emiss_fields()

  PRINT *,'============= emiss_to_buffer done ==============='

END PROGRAM extpar_emiss_to_buffer
