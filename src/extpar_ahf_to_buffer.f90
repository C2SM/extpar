
!+ Fortran main program to read in AHF data and aggregate to target grid
!
!
! Description:
! Fortran main program to read in AHF data and aggregate to target grid
!
! Current Code Owner: DWD, Juergen Helmert
!    Frankfurter Str. 135, 63067 Offenbach, +49-69-8062-2704, juergen.helmert@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_4         2015-05-21 Hendrik Wouters
!  Initial release based on extpar_ndvi_to_buffer.f90 (V1_14)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!> Fortran main program to read in AHF data and aggregate to target grid
!>  
!> \author Hermann Asensio
PROGRAM extpar_ahf_to_buffer

  USE mo_logging
  USE info_extpar,              ONLY: info_print
  USE mo_io_units,              ONLY: filename_max

  USE mo_kind,                  ONLY: wp, i4

  USE mo_target_grid_data,      ONLY: lon_geo, &
       &                              lat_geo, &
       &                              tg


  USE mo_target_grid_routines,  ONLY: init_target_grid
                                
  USE mo_grid_structures,       ONLY: igrid_icon, &
       &                              igrid_cosmo
                                
  USE  mo_icon_grid_data,       ONLY: ICON_grid !< structure which contains the definition of the ICON grid
                                
  USE  mo_cosmo_grid,           ONLY: COSMO_grid
                                
  USE mo_ahf_routines,          ONLY: read_namelists_extpar_ahf

  USE mo_ahf_data,              ONLY: ahf_raw_data_grid, &
       &                              lon_ahf, &
       &                              lat_ahf, &
       &                              allocate_raw_ahf_fields,&
       &                              iahf_type, &
       &                              deallocate_ahf_fields

                               
  USE mo_ahf_tg_fields,         ONLY: ahf_field, &
    &                                 allocate_ahf_target_fields

  USE mo_ahf_routines,          ONLY: open_netcdf_AHF_data, &
    &                                 close_netcdf_AHF_data, &
    &                                 get_dimension_AHF_data, &
    &                                 get_AHF_data_coordinates
                                     

  USE mo_agg_ahf,               ONLY: agg_ahf_data_to_target_grid
                                
  USE mo_ahf_output_nc,         ONLY: write_netcdf_buffer_ahf, &
       &                              write_netcdf_cosmo_grid_ahf, &
       &                              write_netcdf_icon_grid_ahf

  IMPLICIT NONE

  CHARACTER(len=filename_max) :: namelist_grid_def, &
       &                         namelist_ahf_data_input, & !< file with input namelist with AHF data information
       &                         raw_data_ahf_filename, & !< filename ahf raw data
       &                         path_ahf_file, &      !< filename with path for AHF raw data
       &                         netcdf_filename, &      !< filename for netcdf file with AHF data on COSMO grid
       &                         raw_data_ahf_path, &        !< path to raw data
       &                         ahf_buffer_file, & !< name for AHF buffer file
       &                         ahf_output_file !< name for AHF output file

  INTEGER (KIND=i4)           :: ncid_ahf, &  !< netcdf unit file number for AHF data netcdf file
       &                         nlon_ahf, & !< number of grid elements in zonal direction for AHF data
       &                         igrid_type, &  !< target grid type, 1 for ICON, 2 for COSMO
       &                         nlat_ahf !< number of grid elements in meridional direction for AHF data

  REAL (KIND=wp)              :: dlon_ahf, & !< grid point distance in zonal direction (in degrees) for AHF data
       &                         dlat_ahf, & !< grid point distance in meridional direction (in degrees) for AHF data
       &                         startlon_ahf, & !< longitude of lower left grid element for AHF data 
       &                         startlat_ahf, & !< latitude of lower left grid element for AHF data
       &                         undefined !< value to indicate undefined grid elements 

  namelist_grid_def       = 'INPUT_grid_org'
  namelist_ahf_data_input = 'INPUT_AHF'

  CALL initialize_logging("extpar_ahf_to_buffer.log")
  CALL info_print ()
  !--------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= start ahf_to_buffer ==============')
  CALL logging%info( '')
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= init grid and read namelist=======')
  CALL logging%info( '')

  undefined = -999.0 ! undef vlaue
      
  CALL init_target_grid(namelist_grid_def)


  igrid_type = tg%igrid_type
  ! get information on target grid


  CALL  read_namelists_extpar_ahf(namelist_ahf_data_input, &
    &                                  iahf_type,    & !_br 14.04.16
    &                                  raw_data_ahf_path, &
    &                                  raw_data_ahf_filename, &
    &                                  ahf_buffer_file, &
    &                                  ahf_output_file)
     
  path_ahf_file = TRIM(raw_data_ahf_path)//TRIM(raw_data_ahf_filename)
       
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info( '============= allocate fields ==================')
  CALL logging%info( '')

  ! open netcdf file with AHF data
  CALL open_netcdf_AHF_data(path_ahf_file, &
    &                           ncid_ahf)


   !> inquire dimension information for AHF raw data 
  CALL get_dimension_AHF_data(ncid_ahf, &
    &                                nlon_ahf, &
    &                                nlat_ahf)

  CALL allocate_raw_ahf_fields(nlon_ahf,nlat_ahf)   
  CALL allocate_ahf_target_fields(tg)

  CALL get_AHF_data_coordinates(ncid_ahf,      &
    &                               nlon_ahf,      &
    &                               nlat_ahf,      &
    &                               startlon_ahf,  &
    &                               startlat_ahf,  &
    &                               dlon_ahf,      &
    &                               dlat_ahf,      &
    &                               lon_ahf,       &
    &                               lat_ahf)

  ! put the values of the grid definition in the data structure ahf_raw_data_grid (type ahf_reg_lonlat_grid)
  ahf_raw_data_grid%start_lon_reg= startlon_ahf
  ahf_raw_data_grid%start_lat_reg= startlat_ahf
  ahf_raw_data_grid%dlon_reg= dlon_ahf
  ahf_raw_data_grid%dlat_reg= -1. * dlat_ahf ! AHF raw data rows from North to South
  ahf_raw_data_grid%nlon_reg= nlon_ahf
  ahf_raw_data_grid%nlat_reg= nlat_ahf

  ahf_raw_data_grid%end_lon_reg= lon_ahf(nlon_ahf) ! startlon_ahf + (nlon_ahf - 1) * dlon_ahf
  ahf_raw_data_grid%end_lat_reg= lat_ahf(nlat_ahf) ! startlat_ahf - (nlat_ahf - 1) * dlat_ahf 
 ! not negative increment, but AHF latitude goes from north to south

  CALL close_netcdf_AHF_data(ncid_ahf)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= start aggregation ================')
  CALL logging%info( '')

  CALL agg_ahf_data_to_target_grid(tg,undefined, path_ahf_file)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= write data to netcdf==============')
  CALL logging%info( '')

  SELECT CASE(igrid_type)

    CASE(igrid_icon) ! ICON GRID

      netcdf_filename = TRIM(ahf_output_file)
      undefined = -500.

      CALL write_netcdf_icon_grid_ahf(netcdf_filename,  &
   &                                     icon_grid,         &
   &                                     tg,         &
   &                                     undefined, &
   &                                     lon_geo,     &
   &                                     lat_geo,   &
   &                                     ahf_field)





    CASE(igrid_cosmo) ! COSMO grid
    
      netcdf_filename = TRIM(ahf_output_file)
      undefined = -500.

      IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'write out ', TRIM(netcdf_filename)

      CALL write_netcdf_cosmo_grid_ahf(netcdf_filename,  &
   &                                     cosmo_grid,         &
   &                                     tg,         &
   &                                     undefined, &
   &                                     ahf_field)

  END SELECT

  netcdf_filename = TRIM(ahf_buffer_file)
  undefined = -500.

  CALL write_netcdf_buffer_ahf(netcdf_filename,  &
   &                                     tg,         &
   &                                     undefined, &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     ahf_field)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= deallocate fields =================')
  CALL logging%info( '')

  CALL deallocate_ahf_fields()

  CALL logging%info( '')
  CALL logging%info( '============= start ahf_to_buffer ==============')

END PROGRAM extpar_ahf_to_buffer
