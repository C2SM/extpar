!+  Fortran main program to aggregate aerosol optical thickness raw data to target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio
!  small bug fixes accroding to Fortran compiler warnings
! V1_3         2011/04/19 Hermann Asensio
!  clean up of code
! V1_7         2013/01/25 Guenther Zaengl 
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON  
! V2_0         2013/08/08 Daniel Luethi
!   Addition of 2 alternative Aerosol Climatologies
! V4_0         2013/08/17 authors from RHM and Daniel Lthi
!   Addition of support for MACv2 aerosol fields
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!>  Fortran main program to aggregate aerosol optical thickness raw data to target grid
!>
!! @par extpar_aot_to_buffer 
!! 
!! This program interpolates monthly aerosol optical thicknesses for five different types of aerosols 
!! <ul>
!!   <li> black carbon </li>
!!   <li> dust </li>
!!   <li> organic </li>
!!   <li> SO4 </li>
!!   <li> sea salt </li>
!! </ul>
!! iaot_type = 1
!! from a global climatology from Ina Tegen (Tegen et al. 1997) to a target grid (COSMO/ICON). 
!! The raw data and the describing paper are available at NASA/GISS at the Global Aerosol Climatology Project 
!! (GACP http://gacp.giss.nasa.gov/data_sets/transport/). 
!!
!!
!! Tegen, I., P. Hollrigl, M. Chin, I. Fung, D. Jacob, and J. Penner 1997.
!!  <a href="http://pubs.giss.nasa.gov/abstracts/1997/Tegen_etal.html">
!!  Contribution of different aerosol species to the global aerosol extinction optical thickness: 
!!  Estimates from model results</a>.
!! J. Geophys. Res., <b>102</b>, 23895-23915.
!!
!! iaot_type = 2 
!! aerosol climatology from the AEROCOM project
!! (http://aerocom.met.no/aerocomhome.html)
!!
!! Kinne, S., M. Schulz, C. Textor, S. Guibert, Y. Balkanski, S.E. Bauer, 
!! T. Berntsen, T.F. Berglen, O. Boucher, M. Chin, W. Collins, F. Dentener, 
!! T. Diehl, R. Easter, J. Feichter, D. Fillmore, S. Ghan, P. Ginoux, S. Gong, 
!! A. Grini, J. Hendricks, M. Herzog, L. Horowitz, I. Isaksen, T. Iversen, 
!! A. Kirkevg, S. Kloster, D. Koch, J.E. Kristjansson, M. Krol, A. Lauer, 
!! J.F. Lamarque, G. Lesins, X. Liu, U. Lohmann, V. Montanaro, G. Myhre, 
!! J. Penner, G. Pitari, S. Reddy, . Seland, P. Stier, T. Takemura, and X. Tie:
!! An AeroCom initial assessment optical properties in aerosol component modules
!! of global models, Atmos. Chem. Phys., 6, 1815-1834, 2006.
!! 
PROGRAM extpar_aot_to_buffer

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4
  USE info_extpar,              ONLY: info_print
  
  USE mo_target_grid_data,      ONLY: lon_geo, &
    &                                 lat_geo, &
    &                                 tg

  USE mo_target_grid_routines,  ONLY: init_target_grid
 
  USE mo_io_units,              ONLY: filename_max


  USE mo_aot_data,              ONLY: allocate_aot_data, &
    &                                 deallocate_aot_data, &
    &                                 get_dimension_aot_data, &
    &                                 get_aot_grid_and_data, &
    &                                 lon_aot, &
    &                                 lat_aot, &
    &                                 aot_grid, &
    &                                 aot_data, &
    &                                 read_namelists_extpar_aerosol, &
    &                                 iaot_type

  USE mo_agg_aot,               ONLY: agg_aot_data_to_target_grid

  USE mo_aot_target_fields,     ONLY: allocate_aot_target_fields, &
    &                                 aot_tg
  
  USE mo_aot_output_nc,         ONLY: write_netcdf_buffer_aot

  USE mo_io_utilities,          ONLY: join_path 

  IMPLICIT NONE
  
  CHARACTER(len=filename_max) :: netcdf_filename, &
    &                            filename, &
    &                            namelist_grid_def, &
    &                            input_namelist_file, &
    &                            raw_data_aot_path, &        !< path to raw data
    &                            raw_data_aot_filename, & !< filename temperature climatology raw data
    &                            aot_buffer_file !< name for aerosol buffer file

  REAL (KIND=wp) :: undefined

  INTEGER (KIND=i4) :: ntype, & !< number of types of aerosols
                       nrows, & !< number of rows
                       ncolumns, & !< number of columns
                       ntime !< number of times

  !local variables
  input_namelist_file='INPUT_AOT'
  namelist_grid_def = 'INPUT_grid_org'
  undefined = -999.0_wp

  CALL initialize_logging("extpar_aot_to_buffer.log")
  CALL info_print()

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info('')
  CALL logging%info('============= start aot_to_buffer ==============')
  CALL logging%info('')

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info('')
  CALL logging%info('============= read namelist and get dimension ==')
  CALL logging%info('')

  !--------------------------------------------------------------------------------------------------------
  ! get information on target grid, allocate target fields with coordinates and determin the coordinates 
  ! for th target grid

  CALL  init_target_grid(namelist_grid_def)

  !------------------------------------------------------------------------------------

  ! get information about aerosol data
  CALL read_namelists_extpar_aerosol(input_namelist_file, &
   &                                 iaot_type,    &
   &                                 raw_data_aot_path, &
   &                                 raw_data_aot_filename, &
   &                                 aot_buffer_file)


  filename = join_path(raw_data_aot_path,raw_data_aot_filename)

  ! inquire dimensions
  CALL  get_dimension_aot_data(filename, &
                               nrows,        &
                               ncolumns,     &
                               ntime,        &
                               ntype)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  CALL logging%info('')
  CALL logging%info('============= allocate fields ==================')
  CALL logging%info('')
  
  CALL logging%info('l_use_array_cache=.FALSE. -> can only be used in consistency_check')

  ! allocate aot raw data fields
  CALL allocate_aot_data(nrows,ncolumns,ntime,ntype)

  ! allocate target grid fields for aerosol optical thickness
  CALL allocate_aot_target_fields(tg, ntime, ntype, l_use_array_cache=.FALSE.)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  CALL logging%info('')
  CALL logging%info('============= get grid and data ===============')
  CALL logging%info('')

  ! read in aot raw data
  CALL get_aot_grid_and_data(filename,     &
                             nrows,        &
                             ncolumns,     &
                             ntime,        &
                             ntype,        &
                             aot_grid,     &
                             lon_aot,      &
                             lat_aot,      &
                             aot_data)


  aot_tg  =  undefined  ! set target grid values to undefined

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info('')
  CALL logging%info('============= start aggregation ================')
  CALL logging%info('')

  CALL  agg_aot_data_to_target_grid(ntime,ntype)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= write data to netcdf==============')
  CALL logging%info( '')

  netcdf_filename = TRIM(aot_buffer_file)
  CALL write_netcdf_buffer_aot(netcdf_filename, &
      &                       tg,              &
      &                       undefined,       &
      &                       lon_geo,         &
      &                       lat_geo,         &
      &                       ntype,           &
      &                       ntime,           &
      &                       aot_tg,          &
      &                       iaot_type)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= deallocate fields =================')
  CALL logging%info( '')

  CALL deallocate_aot_data()

  CALL logging%info( '')
  CALL logging%info('============= aot_to_buffer done ================')
  
END PROGRAM extpar_aot_to_buffer
