
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

  USE info_extpar, ONLY: info_print
  USE mo_logging

  USE mo_kind, ONLY: wp, i4

  USE mo_target_grid_data, ONLY: lon_geo, &
    &                            lat_geo

  USE mo_target_grid_data, ONLY: tg  !< structure with target grid description

  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo

  USE  mo_icon_grid_data, ONLY: ICON_grid !< structure which contains the definition of the ICON grid
 
  USE  mo_cosmo_grid, ONLY: COSMO_grid

  USE mo_io_units,          ONLY: filename_max

  USE mo_ahf_routines, ONLY: read_namelists_extpar_ahf

  USE mo_ahf_data, ONLY: ahf_raw_data_grid, &
    &                           lon_ahf, &
    &                           lat_ahf, &
    &                           allocate_raw_ahf_fields,&
    &                           deallocate_ahf_fields

   USE mo_ahf_data, ONLY : iahf_type !_br 14.04.16

                               
  USE mo_ahf_tg_fields, ONLY: ahf_field, &
    &                                allocate_ahf_target_fields

  USE mo_ahf_routines, ONLY: open_netcdf_AHF_data, &
    &                               close_netcdf_AHF_data, &
    &                               get_dimension_AHF_data, &
    &                               get_AHF_data_coordinates
                                   

  USE mo_agg_ahf, ONLY: agg_ahf_data_to_target_grid

  USE mo_ahf_output_nc, ONLY: write_netcdf_buffer_ahf
  USE mo_ahf_output_nc, ONLY: write_netcdf_cosmo_grid_ahf
  USE mo_ahf_output_nc, ONLY: write_netcdf_icon_grid_ahf

  IMPLICIT NONE




  CHARACTER(len=filename_max) :: namelist_grid_def

  CHARACTER (len=filename_max) :: namelist_ahf_data_input !< file with input namelist with AHF data information

  CHARACTER (len=filename_max) :: raw_data_ahf_filename !< filename ahf raw data
  CHARACTER (len=filename_max) :: path_ahf_file      !< filename with path for AHF raw data
  CHARACTER (len=filename_max) :: netcdf_filename      !< filename for netcdf file with AHF data on COSMO grid
  CHARACTER (len=filename_max) :: raw_data_ahf_path        !< path to raw data

  CHARACTER (len=filename_max) :: ahf_buffer_file !< name for AHF buffer file
  CHARACTER (len=filename_max) :: ahf_output_file !< name for AHF output file


  INTEGER (KIND=i4) :: ncid_ahf  !< netcdf unit file number for AHF data netcdf file

  INTEGER  (KIND=i4) :: nlon_ahf !< number of grid elements in zonal direction for AHF data
  INTEGER  (KIND=i4) :: nlat_ahf !< number of grid elements in meridional direction for AHF data

  ! INTEGER (KIND=i4):: nmonth  !< index for month for AHF data

  
  REAL (KIND=wp) :: dlon_ahf !< grid point distance in zonal direction (in degrees) for AHF data
  REAL (KIND=wp) :: dlat_ahf !< grid point distance in meridional direction (in degrees) for AHF data

  REAL (KIND=wp) :: startlon_ahf !< longitude of lower left grid element for AHF data 

  REAL (KIND=wp) :: startlat_ahf !< latitude of lower left grid element for AHF data

  INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO

  REAL(KIND=wp) :: undefined !< value to indicate undefined grid elements 
  INTEGER (KIND=i4) :: undef_int   !< value for undefined integer

  namelist_grid_def = 'INPUT_grid_org'
  namelist_ahf_data_input = 'INPUT_AHF'
 ! Print the default information to stdout:
  CALL initialize_logging("extpar_ahf_to_buffer.log")
  CALL info_print ()
  !--------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*) '============= start ahf_to_buffer =============='
  WRITE(logging%fileunit,*) ''
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*) '============= init grid and read namelist======='
  WRITE(logging%fileunit,*) ''

  undef_int = 0 ! set undefined to zero
  undefined = -999.0 ! undef vlaue
      
  CALL init_target_grid(namelist_grid_def)


  IF (verbose >= idbg_low ) THEN
    WRITE(logging%fileunit,*)'target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
    WRITE(logging%fileunit,*)' MAXVAL(lat_geo): ', MAXVAL(lat_geo)
    WRITE(logging%fileunit,*)' MINVAL(lat_geo): ', MINVAL(lat_geo)
  ENDIF


  igrid_type = tg%igrid_type
  ! get information on target grid

  ! read namelist for input AHF data

  CALL  read_namelists_extpar_ahf(namelist_ahf_data_input, &
    &                                  iahf_type,    & !_br 14.04.16
    &                                  raw_data_ahf_path, &
    &                                  raw_data_ahf_filename, &
    &                                  ahf_buffer_file, &
    &                                  ahf_output_file)
     
  path_ahf_file = TRIM(raw_data_ahf_path)//TRIM(raw_data_ahf_filename)
  IF (verbose >= idbg_low ) THEN
    WRITE(logging%fileunit,*) 'after reading namelist for input AHF data, AHF raw data are in file:'
    WRITE(logging%fileunit,*) TRIM(path_ahf_file)
  ENDIF
       
  ! open netcdf file with AHF data
  CALL open_netcdf_AHF_data(path_ahf_file, &
    &                           ncid_ahf)


   !> inquire dimension information for AHF raw data 
   CALL get_dimension_AHF_data(ncid_ahf, &
    &                                nlon_ahf, &
    &                                nlat_ahf)
  IF (verbose >= idbg_low ) THEN
    WRITE(logging%fileunit,*) 'after check of dimensions in AHF raw data file'
    WRITE(logging%fileunit,*) 'nlon_ahf, nlat_ahf: ',nlon_ahf, nlat_ahf
  ENDIF

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

  IF (verbose >= idbg_low ) THEN
    WRITE(logging%fileunit,*) 'after getting AHF data coordinates'
    WRITE(logging%fileunit,*)'startlon_ahf: ', startlon_ahf
    WRITE(logging%fileunit,*)'startlat_ahf: ', startlat_ahf
    WRITE(logging%fileunit,*)'dlon_ahf: ', dlon_ahf
    WRITE(logging%fileunit,*)'dlat_ahf: ', dlat_ahf
    WRITE(logging%fileunit,*)'lon_ahf(1) = ',lon_ahf(1) 
    WRITE(logging%fileunit,*)'lon_ahf(nlon_ahf) = ', lon_ahf(nlon_ahf) 
  ENDIF
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
  IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'ahf_raw_data_grid: ',ahf_raw_data_grid

  CALL close_netcdf_AHF_data(ncid_ahf)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= start aggregation ================'
  WRITE(logging%fileunit,*) ''

  CALL agg_ahf_data_to_target_grid(tg,undefined, path_ahf_file)


  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= write data to netcdf=============='
  WRITE(logging%fileunit,*) ''

  SELECT CASE(igrid_type)

    CASE(igrid_icon) ! ICON GRID

      netcdf_filename = TRIM(ahf_output_file)
      undefined = -500.
      undef_int = -500

      IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'write out ', TRIM(netcdf_filename)

      CALL write_netcdf_icon_grid_ahf(netcdf_filename,  &
   &                                     icon_grid,         &
   &                                     tg,         &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     lon_geo,     &
   &                                     lat_geo,   &
   &                                     ahf_field)





    CASE(igrid_cosmo) ! COSMO grid
    
      netcdf_filename = TRIM(ahf_output_file)
      undefined = -500.
      undef_int = -500

      IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'write out ', TRIM(netcdf_filename)


      CALL write_netcdf_cosmo_grid_ahf(netcdf_filename,  &
   &                                     cosmo_grid,         &
   &                                     tg,         &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     ahf_field)

  END SELECT

  netcdf_filename = TRIM(ahf_buffer_file)
  undefined = -500.
  undef_int = -500

  IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'write out ', TRIM(netcdf_filename)


  CALL write_netcdf_buffer_ahf(netcdf_filename,  &
   &                                     tg,         &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     ahf_field)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= deallocate fields ================='
  WRITE(logging%fileunit,*) ''

  CALL deallocate_ahf_fields()

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*) '============= start ahf_to_buffer =============='

END PROGRAM extpar_ahf_to_buffer

