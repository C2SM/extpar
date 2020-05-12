!+ Fortran main program to aggregate lake depth data to a target grid
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
! V1_14        2014-07-18 Juergen Helmert
!  Combined COSMO Release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran main program to aggregate lake depth data to a target grid
!!
!! @par extpar_flake_to_buffer 
!!
!! 
!! @author
!!     Hermann Asensio
!!     (DWD)
!!
!!
PROGRAM extpar_flake_to_buffer

  USE mo_logging
  USE info_extpar,              ONLY: info_print
  USE mo_io_units,              ONLY: filename_max
  
  USE mo_kind,                  ONLY: wp, i4

  USE mo_target_grid_data,      ONLY: lon_geo, &
      &                               lat_geo

  USE mo_target_grid_data,      ONLY: tg
  
  USE mo_target_grid_routines,  ONLY: init_target_grid

  USE mo_agg_flake,             ONLY : agg_flake_data_to_target_grid

  USE mo_flake_routines,        ONLY: read_namelists_extpar_flake, &
      &                               get_dimension_flake_data, &
      &                               get_lonlat_flake_data
  
  USE mo_flake_data,            ONLY: flake_grid, &
      &                               lon_flake,  &
      &                               lat_flake,  &
      &                               allocate_raw_flake_fields,  &
      &                               deallocate_raw_flake_fields

  USE mo_flake_tg_fields,       ONLY: fr_lake, &
      &                               lake_depth,    &
      &                               flake_tot_npixel, &
      &                               allocate_flake_target_fields
                                

  USE mo_flake_output_nc,      ONLY: write_netcdf_buffer_flake

  
  IMPLICIT NONE

  CHARACTER(len=filename_max) :: netcdf_filename, & 
       &                         namelist_grid_def, & 
       &                         input_flake_namelist_file, &  
       &                         flake_file, & 
       &                         raw_data_flake_path, &         !< path to raw data
       &                         raw_data_flake_filename, &  !< filename flake raw data
       &                         flake_buffer_file  !< name for flake buffer file

  REAL (KIND=wp)              :: undefined

  INTEGER(KIND=i4)            :: undef_int, & 
       &                         nlon_flake, &  !< number of grid elements in zonal direction for flake data
       &                         igrid_type, &   !< target grid type, 1 for ICON, 2 for COSMO
       &                         nlat_flake !< number of grid elements in meridional direction for flake data

  !--------------------------------------------------------------------------------------

  input_flake_namelist_file = 'INPUT_FLAKE'
  namelist_grid_def = 'INPUT_grid_org'

  CALL initialize_logging("extpar_flake_to_buffer.log")
  CALL info_print ()

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= start flake_to_buffer ============')
  CALL logging%info( '')
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= read namelist and get dimension ==')
  CALL logging%info( '')

  CALL init_target_grid(namelist_grid_def)

  igrid_type = tg%igrid_type

  ! get information about FLAE data
  !---------------------------------------------------------------------------
  CALL read_namelists_extpar_flake(input_flake_namelist_file, &
    &                                      raw_data_flake_path, &
                                           raw_data_flake_filename, &
                                           flake_buffer_file)


  flake_file = TRIM(raw_data_flake_path) // TRIM(raw_data_flake_filename)

  CALL get_dimension_flake_data(flake_file, &
    &                                  nlon_flake, &
    &                                  nlat_flake)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info( '============= allocate fields ==================')
  CALL logging%info( '')

  CALL allocate_flake_target_fields(tg, l_use_array_cache=.FALSE.)

  CALL allocate_raw_flake_fields(nlat_flake,nlon_flake)

  CALL  get_lonlat_flake_data(flake_file, &
                                      nlon_flake, &
                                      nlat_flake, &
                                      lon_flake,  &
                                      lat_flake,  &
                                      flake_grid)
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= start aggregation ================')
  CALL logging%info( '')

  undefined = 0.0_wp

  CALL agg_flake_data_to_target_grid(flake_file, &
    &                                      undefined,  &
    &                                      tg,         &
    &                                      lake_depth, &
    &                                      fr_lake,    &
    &                                      flake_tot_npixel)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= write data to netcdf==============')
  CALL logging%info( '')

   undefined = -999.0_wp
   undef_int = -999


   netcdf_filename = TRIM(flake_buffer_file)

   CALL write_netcdf_buffer_flake(TRIM(netcdf_filename),  &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     lake_depth, &
    &                                     fr_lake,    &
    &                                     flake_tot_npixel)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= deallocate fields =================')
  CALL logging%info( '')

  CALL deallocate_raw_flake_fields()

  CALL logging%info( '')
  CALL logging%info('============= flake_to_buffer done ==============')

END PROGRAM extpar_flake_to_buffer
