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

  USE mo_logging
  USE info_extpar,              ONLY: info_print
  USE mo_kind,                  ONLY: wp, i4
  USE mo_io_units,              ONLY: filename_max

  USE mo_target_grid_data,      ONLY: tg, lon_geo, lat_geo

  USE mo_target_grid_routines,  ONLY: init_target_grid

  USE mo_albedo_routines,       ONLY: read_namelists_extpar_alb, &
    &                                 close_netcdf_ALB_data, &
    &                                 get_dimension_ALB_data, &
    &                                 open_netcdf_ALB_data, &
    &                                 get_ALB_data_coordinates

  USE mo_albedo_data,           ONLY: alb_raw_data_grid, &
    &                                 lon_alb, &
    &                                 lat_alb, &
    &                                 ntime_alb, &
    &                                 ialb_type, &
    &                                 undef_alb_bs, &
    &                                 allocate_raw_alb_fields, &
    &                                 deallocate_raw_alb_fields

  USE mo_albedo_tg_fields,     ONLY: alb_dry, &
    &                                alb_sat,     &
    &                                alb_field_mom, &
    &                                alnid_field_mom, &
    &                                aluvd_field_mom, &
    &                                allocate_alb_target_fields, &
    &                                deallocate_alb_target_fields


  USE mo_agg_albedo,           ONLY: agg_alb_data_to_target_grid

  USE mo_albedo_output_nc,     ONLY: write_netcdf_buffer_alb

  IMPLICIT NONE

  CHARACTER(len=filename_max)     :: namelist_grid_def, &
    &                                namelist_alb_data_input, & !< file with input namelist with albedo data information
    &                                raw_data_alb_filename, & !< filename alb raw data
    &                                raw_data_alnid_filename, &
    &                                raw_data_aluvd_filename, &
    &                                path_alb_file, &      !< filename with path for Albedo raw data
    &                                path_alnid_file, &
    &                                path_aluvd_file, &
    &                                netcdf_filename, &      !< filename for netcdf file with Albedo data on COSMO grid
    &                                raw_data_alb_path, &        !< path to raw data
    &                                alb_buffer_file, & !< name for aluvp buffer file
    &                                alb_output_file, & !< name for aluvp output file
    &                                alb_source, alnid_source, aluvd_source, &
    &                                albdry_source, albsat_source
                                  
                                  
  INTEGER (KIND=i4)               :: ncid_alb, &  !< netcdf unit file number for albedo data netcdf file
    &                                nlon_alb, & !< number of grid elements in zonal direction for albedo data
    &                                nlat_alb, & !< number of grid elements in meridional direction for albedo data
    &                                nmonth  !< index for month for albedo data

  INTEGER  (KIND=i4), ALLOCATABLE :: time(:) !< array with time axis values (needed for netcdf cf convention)

  REAL (KIND=wp)                  :: dlon_alb, & !< grid point distance in zonal direction (in degrees) for albedo data
    &                                dlat_alb, & !< grid point distance in meridional direction (in degrees) for albedo data
    &                                startlon_alb, & !< longitude of lower left grid element for albedo data 
    &                                startlat_alb, & !< latitude of lower left grid element for albedo data
    &                                undefined !< value to indicate undefined grid elements 

  !local variables
  undefined = -999.0_wp ! undef vlaue
  namelist_grid_def = 'INPUT_grid_org'
  namelist_alb_data_input = 'INPUT_ALB'

  ! init logger and logfile
  CALL initialize_logging("extpar_alb_to_buffer.log")
  CALL info_print ()

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= start alb_to_buffer ==============')
  CALL logging%info( '')
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= init grid and read namelist=======')
  CALL logging%info( '')

  CALL init_target_grid(namelist_grid_def)
 
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

  
  !generate paths
  path_alb_file = TRIM(raw_data_alb_path)//TRIM(raw_data_alb_filename)
  IF (ialb_type /= 2) THEN
    path_alnid_file = TRIM(raw_data_alb_path)//TRIM(raw_data_alnid_filename)
    path_aluvd_file = TRIM(raw_data_alb_path)//TRIM(raw_data_aluvd_filename)
    alb_source = TRIM(alb_source)
    alnid_source = TRIM(alnid_source)
    aluvd_source = TRIM(aluvd_source)
  ENDIF

  ! open netcdf file with albedo data
  CALL open_netcdf_ALB_data(path_alb_file, ncid_alb)

  !> inquire dimension information for albedo raw data 
  CALL get_dimension_ALB_data(ncid_alb, nlon_alb, nlat_alb, ntime_alb)

  ALLOCATE(time(ntime_alb)) ! this array is needed for netcdf output at the end
  DO nmonth=1, ntime_alb
    time(nmonth) = nmonth
  ENDDO

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= allocate fields ==================')
  CALL logging%info( '')

  CALL allocate_raw_alb_fields(nlon_alb,nlat_alb)   
  CALL allocate_alb_target_fields(tg,ntime_alb,ialb_type)
  
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= get albedo coordinates ===========')
  CALL logging%info( '')

  CALL get_ALB_data_coordinates(ncid_alb,      &
       &                        nlon_alb,      &
       &                        nlat_alb,      &
       &                        startlon_alb,  &
       &                        startlat_alb,  &
       &                        dlon_alb,      &
       &                        dlat_alb,      &
       &                        lon_alb,       &
       &                        lat_alb)


  ! put the values of the grid definition in the data structure alb_raw_data_grid (type alb_reg_lonlat_grid)
  alb_raw_data_grid%start_lon_reg= startlon_alb
  alb_raw_data_grid%start_lat_reg= startlat_alb
  alb_raw_data_grid%dlon_reg= dlon_alb 
  alb_raw_data_grid%dlat_reg= dlat_alb
  alb_raw_data_grid%nlon_reg= nlon_alb
  alb_raw_data_grid%nlat_reg= nlat_alb
  alb_raw_data_grid%end_lon_reg= lon_alb(nlon_alb) 
  alb_raw_data_grid%end_lat_reg= lat_alb(nlat_alb) 

  CALL close_netcdf_ALB_data(ncid_alb)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= start aggregation ================')
  CALL logging%info( '')


  IF (ialb_type == 2) THEN
    albdry_source='ALB_DRY'
    CALL agg_alb_data_to_target_grid(tg,undef_alb_bs, path_alb_file, albdry_source,alb_field_mom)
    alb_dry(:,:,:) = alb_field_mom(:,:,:,1)
    CALL logging%info('aggregation dry soil albedo done')

    albsat_source='ALB_SAT'
    CALL agg_alb_data_to_target_grid(tg,undef_alb_bs, path_alb_file, albsat_source,alb_field_mom)
    alb_sat(:,:,:) = alb_field_mom(:,:,:,1)
    CALL logging%info('aggregation saturated soil albedo done')
  ELSE
    CALL agg_alb_data_to_target_grid(tg,undefined, path_alb_file,  alb_source,alb_field_mom)
    CALL logging%info('aggregation month_albedo done')

    IF (ialb_type == 1) THEN
      CALL agg_alb_data_to_target_grid(tg,undefined, path_alnid_file, alnid_source,alnid_field_mom)
      CALL logging%info('aggregation month_alnid done')

      CALL agg_alb_data_to_target_grid(tg,undefined, path_aluvd_file, aluvd_source,aluvd_field_mom)
      CALL logging%info('aggregation month_aluvd done')
    ENDIF
  ENDIF

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= write data to netcdf==============')
  CALL logging%info( '')

  netcdf_filename = TRIM(alb_buffer_file)

  WRITE(message_text,*)'Write to: ', TRIM(netcdf_filename)
  CALL logging%info(message_text)

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

  CALL logging%info( '')
  CALL logging%info('============= deallocate fields =================')
  CALL logging%info( '')

  CALL deallocate_raw_alb_fields
  CALL deallocate_alb_target_fields()

  CALL logging%info( '')
  CALL logging%info('============= albedo_to_buffer done =============')

END PROGRAM extpar_albedo_to_buffer
