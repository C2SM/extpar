!+  Fortran main program to read in albedo data and aggregate to target grid
!  
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_8         2013/03/12 Frank Brenner
!  introduced MODIS albedo dataset(s) as new external parameter(s) 
! V1_12        2013-04-24 Frank Brenner
!  bug fix regarding old file paths
! V2_0         2013-10-11 Daniel Luethi
!  added support for prescribed bare soil albedo
! V2_0_3       2015-01-12 Juergen Helmert 
!  new orientation of latitudes for ialb_type=1   
! 
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran main program to read in Albedo data and aggregate to target grid
!> Following albedo files are used:
!>              month_aluvp.nc    - UV visible albedo for direct radiation
!>  
!> \author Hermann Asensio, Frank Brenner

PROGRAM extpar_albedo_to_buffer

  USE info_extpar, ONLY: info_print
  USE mo_logging
  USE mo_kind, ONLY: wp, i4

  USE mo_target_grid_data, ONLY: tg, lon_geo, lat_geo

  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_io_units, ONLY: filename_max

  USE mo_albedo_routines, ONLY: read_namelists_extpar_alb

  USE mo_albedo_data, ONLY: alb_raw_data_grid, &
       &                           lon_alb, &
       &                           lat_alb, &
       &                           ntime_alb, &
       &                           ialb_type, &
       &                           undef_alb_bs, &
       &                           allocate_raw_alb_fields, &
       &                           deallocate_raw_alb_fields

  USE mo_albedo_tg_fields, ONLY: alb_dry, &
       &                            alb_sat,     &
       &                            alb_field_mom, &
       &                            alnid_field_mom, &
       &                            aluvd_field_mom, &
       &                            allocate_alb_target_fields, &
       &                            deallocate_alb_target_fields


  USE mo_albedo_routines, ONLY: open_netcdf_ALB_data, &
       &                               close_netcdf_ALB_data, &
       &                               get_dimension_ALB_data, &
       &                               get_ALB_data_coordinates


  USE mo_agg_albedo, ONLY: agg_alb_data_to_target_grid

  USE mo_albedo_output_nc, ONLY: write_netcdf_buffer_alb

  IMPLICIT NONE

  CHARACTER(len=filename_max) :: namelist_grid_def

  CHARACTER (len=filename_max) :: namelist_alb_data_input !< file with input namelist with albedo data information

  CHARACTER (len=filename_max) :: raw_data_alb_filename !< filename alb raw data
  CHARACTER (len=filename_max) :: raw_data_alnid_filename
  CHARACTER (len=filename_max) :: raw_data_aluvd_filename
  CHARACTER (len=filename_max) :: path_alb_file      !< filename with path for Albedo raw data
  CHARACTER (len=filename_max) :: path_alnid_file
  CHARACTER (len=filename_max) :: path_aluvd_file


  CHARACTER (len=filename_max) :: netcdf_filename      !< filename for netcdf file with Albedo data on COSMO grid

  CHARACTER (len=filename_max) :: raw_data_alb_path        !< path to raw data

  CHARACTER (len=filename_max) :: alb_buffer_file !< name for aluvp buffer file
  CHARACTER (len=filename_max) :: alb_output_file !< name for aluvp output file

  CHARACTER (len=filename_max) :: alb_source, alnid_source, aluvd_source, &
                                  albdry_source, albsat_source


  INTEGER (KIND=i4) :: ncid_alb  !< netcdf unit file number for albedo data netcdf file

  INTEGER  (KIND=i4) :: nlon_alb !< number of grid elements in zonal direction for albedo data
  INTEGER  (KIND=i4) :: nlat_alb !< number of grid elements in meridional direction for albedo data



  INTEGER  (KIND=i4), ALLOCATABLE :: time(:) !< array with time axis values (needed for netcdf cf convention)


  INTEGER (KIND=i4):: nmonth  !< index for month for albedo data

  REAL (KIND=wp) :: dlon_alb !< grid point distance in zonal direction (in degrees) for albedo data
  REAL (KIND=wp) :: dlat_alb !< grid point distance in meridional direction (in degrees) for albedo data

  REAL (KIND=wp) :: startlon_alb !< longitude of lower left grid element for albedo data 

  REAL (KIND=wp) :: startlat_alb !< latitude of lower left grid element for albedo data

  REAL(KIND=wp) :: undefined !< value to indicate undefined grid elements 

  !local variables
  undefined = -999.0_wp ! undef vlaue
  namelist_grid_def = 'INPUT_grid_org'
  namelist_alb_data_input = 'INPUT_ALB'

  ! init logger and logfile
  CALL initialize_logging("extpar_alb_to_buffer.log")
  ! print compile-information
  CALL info_print ()

  
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*) '============= start alb_to_buffer =============='
  WRITE(logging%fileunit,*) ''
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*) '============= init grid and read namelist======='
  WRITE(logging%fileunit,*) ''

  CALL init_target_grid(namelist_grid_def)
 
  IF (verbose >= idbg_low ) THEN
    WRITE(logging%fileunit,*)' target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
    WRITE(logging%fileunit,*)' MAXVAL(lat_geo): ', MAXVAL(lat_geo)
    WRITE(logging%fileunit,*)' MINVAL(lat_geo): ', MINVAL(lat_geo)
  ENDIF 

  CALL  read_namelists_extpar_alb(namelist_alb_data_input, &
       &                          raw_data_alb_path,       &
       &                          raw_data_alb_filename,   &
       &                          raw_data_alnid_filename, &
       &                          raw_data_aluvd_filename, &
       &                          ialb_type,               &
       &                          alb_buffer_file,         &
       &                          alb_output_file,         &
       &                          alb_source,              &
       &                          alnid_source,            &
       &                          aluvd_source)

  IF (ialb_type == 2) THEN
    WRITE(logging%fileunit,*)'INFO: MODIS albedo data (south -> north storage direction)'
  ELSE
    WRITE(logging%fileunit,*)'INFO: albedo data (north -> south storage direction)'
  ENDIF
  
  !generate paths
  path_alb_file = TRIM(raw_data_alb_path)//TRIM(raw_data_alb_filename)
  IF (ialb_type /= 2) THEN
    path_alnid_file = TRIM(raw_data_alb_path)//TRIM(raw_data_alnid_filename)
    path_aluvd_file = TRIM(raw_data_alb_path)//TRIM(raw_data_aluvd_filename)
    alb_source = TRIM(alb_source)
    alnid_source = TRIM(alnid_source)
    aluvd_source = TRIM(aluvd_source)
    IF (verbose >= idbg_low ) THEN
      WRITE(logging%fileunit,*) 'after reading namelist for input Albedo data, Albedo raw data are in file:'
      WRITE(logging%fileunit,*)TRIM(path_alb_file)
      WRITE(logging%fileunit,*)TRIM(path_alnid_file)
      WRITE(logging%fileunit,*)TRIM(path_aluvd_file)
      WRITE(logging%fileunit,*)'name of buffer file: ', TRIM(alb_buffer_file)
      WRITE(logging%fileunit,*)'name of output file: ', TRIM(alb_output_file)
    ENDIF
  ENDIF

  ! open netcdf file with albedo data
  CALL open_netcdf_ALB_data(path_alb_file, ncid_alb)

  !> inquire dimension information for albedo raw data 
  CALL get_dimension_ALB_data(ncid_alb, nlon_alb, nlat_alb, ntime_alb)

  IF (verbose >= idbg_low ) THEN
   WRITE(logging%fileunit,*)'after check of dimensions in Albedo raw data file'
   WRITE(logging%fileunit,*)'nlon_alb, nlat_alb: ',nlon_alb, nlat_alb
   WRITE(logging%fileunit,*)'ntime_alb: ', ntime_alb
  ENDIF

  ALLOCATE(time(ntime_alb)) ! this array is needed for netcdf output at the end
  DO nmonth=1, ntime_alb
    time(nmonth) = nmonth
  ENDDO

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= allocate fields =================='
  WRITE(logging%fileunit,*) ''

  CALL allocate_raw_alb_fields(nlon_alb,nlat_alb)   
  CALL allocate_alb_target_fields(tg,ntime_alb,ialb_type)
  
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= get albedo coordinates ==========='
  WRITE(logging%fileunit,*) ''

  CALL get_ALB_data_coordinates(ncid_alb,      &
       &                        nlon_alb,      &
       &                        nlat_alb,      &
       &                        startlon_alb,  &
       &                        startlat_alb,  &
       &                        dlon_alb,      &
       &                        dlat_alb,      &
       &                        lon_alb,       &
       &                        lat_alb)

  IF (verbose >= idbg_low ) THEN
    WRITE(logging%fileunit,*) 'after getting Albedo data coordinates'
    WRITE(logging%fileunit,*)'startlon_alb: ', startlon_alb
    WRITE(logging%fileunit,*)'startlat_alb: ', startlat_alb
    WRITE(logging%fileunit,*)'dlon_alb: ', dlon_alb
    WRITE(logging%fileunit,*)'dlat_alb: ', dlat_alb
    WRITE(logging%fileunit,*)'lon_alb(1) = ',lon_alb(1) 
    WRITE(logging%fileunit,*)'lon_alb(nlon_alb) = ', lon_alb(nlon_alb) 
  ENDIF

  ! put the values of the grid definition in the data structure alb_raw_data_grid (type alb_reg_lonlat_grid)
  alb_raw_data_grid%start_lon_reg= startlon_alb
  alb_raw_data_grid%start_lat_reg= startlat_alb
  alb_raw_data_grid%dlon_reg= dlon_alb 
  alb_raw_data_grid%dlat_reg= dlat_alb
  alb_raw_data_grid%nlon_reg= nlon_alb
  alb_raw_data_grid%nlat_reg= nlat_alb
  alb_raw_data_grid%end_lon_reg= lon_alb(nlon_alb) 
  alb_raw_data_grid%end_lat_reg= lat_alb(nlat_alb) 

  IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'alb_raw_data_grid: ',alb_raw_data_grid
  CALL close_netcdf_ALB_data(ncid_alb)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= start aggregation ================'
  WRITE(logging%fileunit,*) ''

  IF (verbose >= idbg_low) WRITE(logging%fileunit,*)'aggregate Albedo data to target grid'

  IF (ialb_type == 2) THEN
    albdry_source='ALB_DRY'
    CALL agg_alb_data_to_target_grid(tg,undef_alb_bs, path_alb_file, albdry_source,alb_field_mom)
    alb_dry(:,:,:) = alb_field_mom(:,:,:,1)
    IF (verbose >= idbg_low) WRITE(logging%fileunit,*)'aggregation dry soil albedo done'

    albsat_source='ALB_SAT'
    CALL agg_alb_data_to_target_grid(tg,undef_alb_bs, path_alb_file, albsat_source,alb_field_mom)
    alb_sat(:,:,:) = alb_field_mom(:,:,:,1)
    IF (verbose >= idbg_low) WRITE(logging%fileunit,*)'aggregation saturated soil albedo done'
  ELSE
    CALL agg_alb_data_to_target_grid(tg,undefined, path_alb_file,  alb_source,alb_field_mom)
    IF (verbose >= idbg_low) WRITE(logging%fileunit,*)'aggregation month_albedo done'

    IF (ialb_type == 1) THEN
      CALL agg_alb_data_to_target_grid(tg,undefined, path_alnid_file, alnid_source,alnid_field_mom)
      IF (verbose >= idbg_low) WRITE(logging%fileunit,*)'aggregation month_alnid done'

      CALL agg_alb_data_to_target_grid(tg,undefined, path_aluvd_file, aluvd_source,aluvd_field_mom)
      IF (verbose >= idbg_low) WRITE(logging%fileunit,*)'aggregation month_aluvd done'
    ENDIF
  ENDIF

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= write data to netcdf=============='
  WRITE(logging%fileunit,*) ''

  IF (verbose >= idbg_low) THEN
    WRITE(logging%fileunit,*)'MAX albedo values:'
    IF (ialb_type == 1) THEN
      WRITE(logging%fileunit,*)MAXVAL(alb_field_mom), MAXVAL(alnid_field_mom), MAXVAL(aluvd_field_mom)
    ELSE IF (ialb_type == 2) THEN
      WRITE(logging%fileunit,*) MAXVAL(alb_sat), MAXVAL(alb_dry)
    ELSE IF (ialb_type == 3) THEN
      WRITE(logging%fileunit,*) MAXVAL(alb_field_mom)
    ENDIF
  ENDIF

  netcdf_filename = TRIM(alb_buffer_file)

  IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'write to: ', TRIM(netcdf_filename)
  IF (ialb_type == 2) THEN
    CALL write_netcdf_buffer_alb(netcdf_filename,  &
         &                            tg,         &
         &                            ntime_alb, &
         &                            undefined, &
         &                            lon_geo,     &
         &                            lat_geo, &
         &                            alb_dry=alb_dry, &
         &                            alb_sat=alb_sat)
  ELSE IF (ialb_type == 1) THEN
    CALL write_netcdf_buffer_alb(netcdf_filename,  &
         &                            tg,         &
         &                            ntime_alb, &
         &                            undefined, &
         &                            lon_geo,     &
         &                            lat_geo, &
         &                            alb_field_mom=alb_field_mom, &
         &                            alnid_field_mom=alnid_field_mom, &
         &                            aluvd_field_mom=aluvd_field_mom)
  ELSE IF (ialb_type == 3) THEN
    CALL write_netcdf_buffer_alb(netcdf_filename,  &
         &                            tg,         &
         &                            ntime_alb, &
         &                            undefined, &
         &                            lon_geo,     &
         &                            lat_geo, &
         &                            alb_field_mom=alb_field_mom)
  ENDIF

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= deallocate fields ================='
  WRITE(logging%fileunit,*) ''

  CALL deallocate_raw_alb_fields
  CALL deallocate_alb_target_fields()

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= albedo_to_buffer done ============='

END PROGRAM extpar_albedo_to_buffer
