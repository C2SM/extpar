!+  Fortran main program to read in NDVI data and aggregate to target grid
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
!> Fortran main program to read in NDVI data and aggregate to target grid
!>  
!> \author Hermann Asensio
PROGRAM extpar_ndvi_to_buffer

  USE mo_logging
  USE info_extpar,              ONLY: info_print
  USE mo_kind,                  ONLY: wp, i4
  
  USE mo_io_units,              ONLY: filename_max

  USE mo_target_grid_data,      ONLY: lon_geo, &
       &                              tg, &
       &                              lat_geo

  USE mo_target_grid_routines,  ONLY: init_target_grid

  USE mo_grid_structures,       ONLY: igrid_icon, &
       &                              igrid_cosmo
                                
  USE mo_icon_grid_data,        ONLY: ICON_grid !< structure which contains the definition of the ICON grid
 
  USE mo_cosmo_grid,            ONLY: COSMO_grid
    
  USE mo_ndvi_data,             ONLY: ndvi_raw_data_grid, &
       &                              lon_ndvi, &
       &                              lat_ndvi, &
       &                              ntime_ndvi, &
       &                              allocate_raw_ndvi_fields,&
       &                              deallocate_ndvi_fields
                               
  USE mo_ndvi_tg_fields,        ONLY: ndvi_max, &
       &                              ndvi_field_mom, &
       &                              ndvi_ratio_mom, &
       &                              allocate_ndvi_target_fields

  USE mo_ndvi_routines,         ONLY: open_netcdf_NDVI_data, &
       &                              close_netcdf_NDVI_data, &
       &                              get_dimension_NDVI_data, &
       &                              read_namelists_extpar_ndvi, &
       &                              get_NDVI_data_coordinates

  USE mo_agg_ndvi,              ONLY: agg_ndvi_data_to_target_grid

  USE mo_ndvi_output_nc,        ONLY: write_netcdf_buffer_ndvi, &
       &                              write_netcdf_cosmo_grid_ndvi, &
       &                              write_netcdf_icon_grid_ndvi

  IMPLICIT NONE

  CHARACTER(len=filename_max)      :: namelist_grid_def, &
       &                              namelist_ndvi_data_input, & !< file with input namelist with NDVI data information
       &                              raw_data_ndvi_filename, & !< filename ndvi raw data
       &                              path_ndvi_file, &      !< filename with path for NDVI raw data
       &                              netcdf_filename, &      !< filename for netcdf file with NDVI data on COSMO grid
       &                              raw_data_ndvi_path, &        !< path to raw data
       &                              ndvi_buffer_file, & !< name for NDVI buffer file
       &                              ndvi_output_file !< name for NDVI output file
                                   
                                   
  INTEGER (KIND=i4)                :: ncid_ndvi, &  !< netcdf unit file number for NDVI data netcdf file
       &                              nlon_ndvi, & !< number of grid elements in zonal direction for NDVI data
       &                              nlat_ndvi, & !< number of grid elements in meridional direction for NDVI data
       &                              nmonth, &  !< index for month for NDVI data
       &                              igrid_type  !< target grid type, 1 for ICON, 2 for COSMO

  INTEGER  (KIND=i4), ALLOCATABLE :: time(:) !< array with time axis values (needed for netcdf cf convention)

  REAL (KIND=wp)                  :: dlon_ndvi, & !< grid point distance in zonal direction (in degrees) for NDVI data
       &                             dlat_ndvi, & !< grid point distance in meridional direction (in degrees) for NDVI data
       &                             startlon_ndvi, & !< longitude of lower left grid element for NDVI data 
       &                             startlat_ndvi, & !< latitude of lower left grid element for NDVI data
       &                             undefined !< value to indicate undefined grid elements 

  undefined = -999.0 ! undef vlaue
  namelist_grid_def =        'INPUT_grid_org'
  namelist_ndvi_data_input = 'INPUT_NDVI'

  CALL initialize_logging("extpar_ndvi_to_buffer.log")
  CALL info_print ()

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= start ndvi_to_buffer =============')
  CALL logging%info( '')

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= read namelist and init grid ======')
  CALL logging%info( '')
      
  CALL init_target_grid(namelist_grid_def)

  igrid_type = tg%igrid_type
  ! get information on target grid

  ! read namelist for input NDVI data
  CALL  read_namelists_extpar_ndvi(namelist_ndvi_data_input, &
    &                                  raw_data_ndvi_path, &
    &                                  raw_data_ndvi_filename, &
    &                                  ndvi_buffer_file, &
    &                                  ndvi_output_file)
     
  path_ndvi_file = TRIM(raw_data_ndvi_path)//TRIM(raw_data_ndvi_filename)

  ! open netcdf file with NDVI data
  CALL open_netcdf_NDVI_data(path_ndvi_file, &
    &                           ncid_ndvi)

   !> inquire dimension information for NDVI raw data 
   CALL get_dimension_NDVI_data(ncid_ndvi, &
    &                                nlon_ndvi, &
    &                                nlat_ndvi, &
                                    ntime_ndvi)

  ALLOCATE(time(ntime_ndvi)) ! this array is needed for netcdf output at the end
  DO nmonth=1, ntime_ndvi
    time(nmonth) = nmonth
  ENDDO

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info( '============= allocate fields ==================')
  CALL logging%info( '')

  CALL allocate_raw_ndvi_fields(nlon_ndvi,nlat_ndvi,ntime_ndvi)   
  CALL allocate_ndvi_target_fields(tg,ntime_ndvi)

  CALL get_NDVI_data_coordinates(ncid_ndvi,      &
    &                               nlon_ndvi,      &
    &                               nlat_ndvi,      &
    &                               startlon_ndvi,  &
    &                               startlat_ndvi,  &
    &                               dlon_ndvi,      &
    &                               dlat_ndvi,      &
    &                               lon_ndvi,       &
    &                               lat_ndvi)


  ! put the values of the grid definition in the data structure ndvi_raw_data_grid (type ndvi_reg_lonlat_grid)
  ndvi_raw_data_grid%start_lon_reg= startlon_ndvi
  ndvi_raw_data_grid%start_lat_reg= startlat_ndvi
  ndvi_raw_data_grid%dlon_reg= dlon_ndvi
  ndvi_raw_data_grid%dlat_reg= dlat_ndvi
  ndvi_raw_data_grid%nlon_reg= nlon_ndvi
  ndvi_raw_data_grid%nlat_reg= nlat_ndvi

  ndvi_raw_data_grid%end_lon_reg= lon_ndvi(nlon_ndvi) ! startlon_ndvi + (nlon_ndvi - 1) * dlon_ndvi
  ndvi_raw_data_grid%end_lat_reg= lat_ndvi(nlat_ndvi) ! startlat_ndvi - (nlat_ndvi - 1) * dlat_ndvi 
 ! not negative increment, but NDVI latitude goes from north to south

  CALL close_netcdf_NDVI_data(ncid_ndvi)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= start aggregation ================')
  CALL logging%info( '')

  CALL agg_ndvi_data_to_target_grid(tg, path_ndvi_file)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= write data to netcdf==============')
  CALL logging%info( '')

  SELECT CASE(igrid_type)

    CASE(igrid_icon) ! ICON GRID

      netcdf_filename = TRIM(ndvi_output_file)
      undefined = -500.

      CALL write_netcdf_icon_grid_ndvi(netcdf_filename,  &
           &                           icon_grid,         &
           &                           tg,         &
           &                           ntime_ndvi, &
           &                           undefined, &
           &                           lon_geo,     &
           &                           lat_geo, &
           &                           ndvi_max,  &
           &                           ndvi_field_mom,&
           &                           ndvi_ratio_mom)





    CASE(igrid_cosmo) ! COSMO grid
    
      netcdf_filename = TRIM(ndvi_output_file)
      undefined = -500.

      CALL write_netcdf_cosmo_grid_ndvi(netcdf_filename,  &
           &                            cosmo_grid,         &
           &                            tg,         &
           &                            ntime_ndvi, &
           &                            undefined, &
           &                            ndvi_max,  &
           &                            ndvi_field_mom,&
           &                            ndvi_ratio_mom)

  END SELECT

  netcdf_filename = TRIM(ndvi_buffer_file)
  undefined = -500.

  CALL write_netcdf_buffer_ndvi(netcdf_filename,  &
       &                        tg,         &
       &                        ntime_ndvi, &
       &                        undefined, &
       &                        lon_geo,     &
       &                        lat_geo, &
       &                        ndvi_max,  &
       &                        ndvi_field_mom,&
       &                        ndvi_ratio_mom)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= deallocate fields =================')
  CALL logging%info( '')

  CALL deallocate_ndvi_fields()

  CALL logging%info( '')
  CALL logging%info('============= ndvi_to_buffer done ===============')

END PROGRAM extpar_ndvi_to_buffer
