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
!! iaot_type = 3
!! MACC-II dataset from ECMWF
!! (http://www.gmes-atmosphere.eu/)
!!
!! Morcrette, J.-J., O. Boucher, L. Jones, D. Salmond, P. Bechtold, A. Beljaars,
!! A. Benedetti, A. Bonet, J. W. Kaiser, M. Razinger, M. Schulz, S. Serrar,
!! A. J. Simmons, M. Sofiev, M. Suttie, A. M. Tompkins, and A. Untch, 2009: 
!! Aerosol analysis and forecast in the ECMWF Integrated Forecast System. 
!! Part I: Forward modelling, J. Geophys. Res., 114D, D06206, 
!! doi:10.1029/2008JD011235
!!
!! iaot_type = 4 
!! MACv2 dataset from AEROCOM project
!! ftp://ftp-projects.zmaw.de/aerocom/climatology/MACv2_2015/MACv2_13_mosk/2005/
!! 
!! Kinne, S., D. O'Donnel, P. Stier, S. Kloster, K. Zhang, H. Schmidt, S. Rast,
!! M. Giorgetta, T. F. Eck, and B. Stevens (2013), 
!! MAC-v1: A new global aerosol climatology for climate studies, 
!! J. Adv. Model. Earth Syst., 5, 704740, doi:10.1002/jame.20035
PROGRAM extpar_aot_to_buffer

  USE info_extpar, ONLY: info_print
  USE mo_logging
  USE mo_kind,              ONLY: wp, i8
  
  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo

  USE mo_target_grid_data, ONLY: lon_geo, &
    &                            lat_geo

  USE mo_target_grid_data, ONLY: tg
  USE mo_target_grid_routines, ONLY: init_target_grid
 
 
  USE  mo_icon_grid_data, ONLY: ICON_grid !< structure which contains the definition of the ICON grid
   
  USE  mo_cosmo_grid, ONLY: cosmo_grid

  USE mo_io_units,          ONLY: filename_max

  USE mo_aot_data, ONLY:  read_namelists_extpar_aerosol

  USE mo_aot_data, ONLY : allocate_aot_data, &
    &                      deallocate_aot_data, &
    &                      get_dimension_aot_data, &
    &                      get_aot_grid_and_data, &
    &                      lon_aot, &
    &                      lat_aot, &
    &                      aot_grid

  USE mo_aot_data, ONLY : iaot_type, n_spectr

  USE mo_aot_target_fields, ONLY: allocate_aot_target_fields,&
    &                              aot_tg,&
    &                              MAC_aot_tg,&
    &                              MAC_ssa_tg,&
    &                              MAC_asy_tg
  
  USE mo_agg_aot, ONLY: agg_aot_data_to_target_grid

  USE mo_aot_output_nc, ONLY: write_netcdf_buffer_aot, &
    &                         write_netcdf_cosmo_grid_aot, &
    &                         write_netcdf_icon_grid_aot


  IMPLICIT NONE
  
  CHARACTER(len=filename_max) :: netcdf_filename
  CHARACTER(len=filename_max) :: filename

  CHARACTER(len=filename_max) :: namelist_grid_def
  CHARACTER(len=filename_max) :: input_namelist_file

  CHARACTER (len=filename_max) :: raw_data_aot_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_aot_filename !< filename temperature climatology raw data

  CHARACTER (len=filename_max) :: aot_buffer_file !< name for aerosol buffer file
  CHARACTER (len=filename_max) :: aot_output_file !< name for aerosol output file

  REAL (KIND=wp) :: undefined
  INTEGER        :: undef_int

  !--------------------------------------------------------------------------------------

  INTEGER (KIND=i8) :: ntype !< number of types of aerosols
  INTEGER (KIND=i8) :: nrows !< number of rows
  INTEGER (KIND=i8) :: ncolumns !< number of columns
  INTEGER (KIND=i8) :: ntime !< number of times

  CALL initialize_logging("extpar_aot_to_buffer.log", stdout_level=debug)
  CALL info_print()
  !--------------------------------------------------------------------------------------------------------
  ! get information on target grid, allocate target fields with coordinates and determin the coordinates 
  ! for th target grid

  namelist_grid_def = 'INPUT_grid_org'
  CALL  init_target_grid(namelist_grid_def)

  PRINT *,'target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
  !------------------------------------------------------------------------------------

  ! get information about aerosol data
  input_namelist_file='INPUT_AOT'

  CALL read_namelists_extpar_aerosol(input_namelist_file, &
   &                                  iaot_type,    &
   &                                  raw_data_aot_path, &
   &                                  raw_data_aot_filename, &
   &                                  aot_buffer_file, &
   &                                  aot_output_file)


   filename = TRIM(raw_data_aot_path) // TRIM(raw_data_aot_filename)

   PRINT *,'filename ',TRIM(filename)

   ! inquire dimensions

   CALL  get_dimension_aot_data(filename, &
                                     iaot_type,    &
                                     nrows,        &
                                     ncolumns,     &
                                     ntime,        &
                                     ntype,        &
                                     n_spectr)

   PRINT *, 'nrows: ',nrows
   PRINT *, 'ncolumns: ',ncolumns
   PRINT *, 'ntime: ',ntime
   PRINT *, 'ntype: ',ntype
   PRINT *, 'n_spectr: ',n_spectr
   PRINT *, 'iaot_type: ',iaot_type
   ! allocate aot raw data fields

   CALL allocate_aot_data(iaot_type,nrows,ncolumns,ntime,ntype,n_spectr)


   ! read in aot raw data

   CALL get_aot_grid_and_data(iaot_type,           &
                                     filename,     &
                                     nrows,        &
                                     ncolumns,     &
                                     ntime,        &
                                     ntype,        &
                                     n_spectr,     &
                                     aot_grid,     &
                                     lon_aot,      &
                                     lat_aot)

    ! allocate target grid fields for aerosol optical thickness

    CALL allocate_aot_target_fields(tg, iaot_type, ntime, ntype, n_spectr)
    
    undefined = -999.0_wp
    undef_int = -999

    IF (iaot_type == 4) THEN
      MAC_aot_tg =  undefined
      MAC_ssa_tg =  undefined
      MAC_asy_tg =  undefined
    ELSE
      aot_tg  =  undefined  ! set target grid values to undefined
    ENDIF

    PRINT *,'call agg_aot_data_to_target_grid'
    CALL  agg_aot_data_to_target_grid(iaot_type,nrows,ncolumns,ntime,ntype,n_spectr)

    netcdf_filename = TRIM(aot_buffer_file)
    CALL logging%info('write BUFFER output to '//TRIM(netcdf_filename), __FILE__, __LINE__)    
    CALL write_netcdf_buffer_aot(netcdf_filename, &
         &                       tg,              &
         &                       undefined,       &
         &                       undef_int,       &
         &                       lon_geo,         &
         &                       lat_geo,         &
         &                       ntype,           &
         &                       ntime,           &
         &                       n_spectr,        &
         &                       aot_tg,          &
         &                       MAC_aot_tg,      &
         &                       MAC_ssa_tg,      &
         &                       MAC_asy_tg,      &
         &                       iaot_type)

    !write out data
    netcdf_filename =  TRIM(aot_output_file)

    SELECT CASE(tg%igrid_type)

    CASE(igrid_icon) ! ICON GRID
      
        CALL logging%info('write ICON output to '//TRIM(aot_output_file), __FILE__, __LINE__)
        CALL write_netcdf_icon_grid_aot(netcdf_filename,  &
   &                                     icon_grid,       &
   &                                     tg,              &
   &                                     undefined,       &
   &                                     undef_int,       &
   &                                     lon_geo,         &
   &                                     lat_geo,         &
   &                                     ntype,           &
   &                                     ntime,           &
   &                                     n_spectr,        &
   &                                     aot_tg,          &
   &                                     MAC_aot_tg, &
   &                                     MAC_ssa_tg, &
   &                                     MAC_asy_tg, &
   &                                     iaot_type)

      CASE(igrid_cosmo) ! COSMO grid

        CALL logging%info('write COSMO output to '//TRIM(aot_output_file), __FILE__, __LINE__)        
        CALL write_netcdf_cosmo_grid_aot(netcdf_filename, &
   &                                     cosmo_grid,      &
   &                                     tg,              &
   &                                     undefined,       &
   &                                     undef_int,       &
   &                                     lon_geo,         &
   &                                     lat_geo,         &
   &                                     ntype,           &
   &                                     ntime,           &
   &                                     n_spectr,        &
   &                                     aot_tg,          &
   &                                     MAC_aot_tg, &
   &                                     MAC_ssa_tg, &
   &                                     MAC_asy_tg, &
   &                                     iaot_type)

    END SELECT

    CALL deallocate_aot_data()

    PRINT *,'============= extpar_aot_to_buffer done ==============='
  
END PROGRAM extpar_aot_to_buffer
