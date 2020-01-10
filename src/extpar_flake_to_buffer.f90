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

  USE info_extpar, ONLY: info_print
  USE mo_logging
  
  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo

  USE mo_target_grid_data, ONLY: lon_geo, &
    &                            lat_geo

  USE mo_target_grid_data, ONLY: tg
  
  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_icon_grid_data, ONLY: ICON_grid  !< structure which contains the definition of the ICON grid
 
  USE  mo_cosmo_grid, ONLY: COSMO_grid

  USE mo_io_units,          ONLY: filename_max

  USE mo_flake_routines, ONLY: read_namelists_extpar_flake

  USE mo_flake_routines, ONLY:  get_dimension_flake_data, &
    &                             get_lonlat_flake_data
  

  USE mo_flake_data, ONLY: flake_grid, &
 &         lon_flake,  &
 &         lat_flake,  &
 &         allocate_raw_flake_fields,  &
 &         deallocate_raw_flake_fields

  USE mo_flake_tg_fields, ONLY: fr_lake, &
  &       lake_depth,    &
  &       flake_tot_npixel, &
  &       allocate_flake_target_fields

  USE mo_agg_flake, ONLY : agg_flake_data_to_target_grid


  USE mo_flake_output_nc, ONLY: write_netcdf_buffer_flake, &
    &                             write_netcdf_cosmo_grid_flake, &
    &                             write_netcdf_icon_grid_flake

  
  IMPLICIT NONE

  CHARACTER(len=filename_max) :: netcdf_filename

  CHARACTER(len=filename_max) :: namelist_grid_def

  CHARACTER(len=filename_max) :: input_flake_namelist_file 
  CHARACTER(len=filename_max) :: flake_file

  CHARACTER (len=filename_max) :: raw_data_flake_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_flake_filename !< filename flake raw data

  CHARACTER (len=filename_max) :: flake_buffer_file !< name for flake buffer file
  CHARACTER (len=filename_max) :: flake_output_file !< name for flake output file

  REAL (KIND=wp) :: undefined
  INTEGER :: undef_int


  INTEGER (KIND=i8) :: nlon_flake !< number of grid elements in zonal direction for flake data
  INTEGER (KIND=i8) :: nlat_flake !< number of grid elements in meridional direction for flake data

  !--------------------------------------------------------------------------------------

  INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO

  input_flake_namelist_file = 'INPUT_FLAKE'
  namelist_grid_def = 'INPUT_grid_org'

  CALL initialize_logging("extpar_flake_to_buffer.log")
  CALL info_print ()

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*) '============= read namelist and get dimension =='
  WRITE(logging%fileunit,*) ''

  CALL init_target_grid(namelist_grid_def)

  IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat

  igrid_type = tg%igrid_type

  CALL allocate_flake_target_fields(tg)

  ! get information about FLAE data
  !---------------------------------------------------------------------------
  CALL read_namelists_extpar_flake(input_flake_namelist_file, &
    &                                      raw_data_flake_path, &
                                           raw_data_flake_filename, &
                                           flake_buffer_file, &
                                           flake_output_file)


  IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'raw_data_flake_filename: ',TRIM(raw_data_flake_filename)
  flake_file = TRIM(raw_data_flake_path) // TRIM(raw_data_flake_filename)

  IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'flake file: ', TRIM(flake_file)

  CALL get_dimension_flake_data(flake_file, &
    &                                  nlon_flake, &
    &                                  nlat_flake)


  IF (verbose >= idbg_low ) THEN
    WRITE(logging%fileunit,*)'nlon_flake: ',nlon_flake
    WRITE(logging%fileunit,*)'nlat_flake: ',nlat_flake
  ENDIF

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*) '============= allocate fields =================='
  WRITE(logging%fileunit,*) ''

  CALL allocate_raw_flake_fields(nlat_flake,nlon_flake)

  CALL  get_lonlat_flake_data(flake_file, &
                                      nlon_flake, &
                                      nlat_flake, &
                                      lon_flake,  &
                                      lat_flake,  &
                                      flake_grid)

  IF (verbose >= idbg_low ) THEN
    WRITE(logging%fileunit,*)'MINVAL(lat_flake) :', MINVAL(lat_flake)
    WRITE(logging%fileunit,*)'MAXVAL(lat_flake) :', MAXVAL(lat_flake)
    WRITE(logging%fileunit,*)'flake_grid: ', flake_grid
    WRITE(logging%fileunit,*)'MINVAL(lon_flake) :', MINVAL(lon_flake)
    WRITE(logging%fileunit,*)'MAXVAL(lon_flake) :', MAXVAL(lon_flake)
    WRITE(logging%fileunit,*)'MINVAL(lat_flake) :', MINVAL(lat_flake)
    WRITE(logging%fileunit,*)'MAXVAL(lat_flake) :', MAXVAL(lat_flake)
  ENDIF


  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= start aggregation ================'
  WRITE(logging%fileunit,*) ''

  undefined = 0.0_wp

  CALL agg_flake_data_to_target_grid(flake_file, &
    &                                      undefined,  &
    &                                      tg,         &
    &                                      lake_depth, &
    &                                      fr_lake,    &
    &                                      flake_tot_npixel)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= write data to netcdf=============='
  WRITE(logging%fileunit,*) ''

   undefined = -999.0_wp
   undef_int = -999


   netcdf_filename = TRIM(flake_buffer_file)
   WRITE(logging%fileunit,*) 'INFO: Flake data buffer filename: ',TRIM(netcdf_filename)

   CALL write_netcdf_buffer_flake(TRIM(netcdf_filename),  &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     lake_depth, &
    &                                     fr_lake,    &
    &                                     flake_tot_npixel)



    SELECT CASE(igrid_type)

      CASE(igrid_icon) ! ICON GRID
        
        netcdf_filename = TRIM(flake_output_file)
        IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'write out ', TRIM(netcdf_filename)

        CALL write_netcdf_icon_grid_flake(TRIM(netcdf_filename),  &
    &                                     icon_grid,       &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     lake_depth, &
    &                                     fr_lake,    &
    &                                     flake_tot_npixel)




         
      CASE(igrid_cosmo) ! COSMO grid

         netcdf_filename = TRIM(flake_output_file)
         IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'write out ', TRIM(netcdf_filename)

         CALL write_netcdf_cosmo_grid_flake(TRIM(netcdf_filename), &
    &                                     cosmo_grid,       &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     lake_depth, &
    &                                     fr_lake,    &
    &                                     flake_tot_npixel)


    END SELECT

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= deallocate fields ================='
  WRITE(logging%fileunit,*) ''

    CALL deallocate_raw_flake_fields()

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= flake_to_buffer done =============='

END PROGRAM extpar_flake_to_buffer
