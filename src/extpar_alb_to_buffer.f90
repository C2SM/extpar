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

  USE info_extpar, ONLY: info_define, info_print

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
  CHARACTER (len=filename_max) :: filename
  CHARACTER (len=filename_max) :: raw_data_alb_path        !< path to raw data

  CHARACTER (len=filename_max) :: alb_buffer_file !< name for aluvp buffer file
  CHARACTER (len=filename_max) :: alb_output_file !< name for aluvp output file

  CHARACTER (len=filename_max) :: alb_source, alnid_source, aluvd_source


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
  INTEGER (KIND=i4) :: undef_int   !< value for undefined integer

  CALL info_define ('alb_to_buffer')      ! Pre-define the program name as binary name
  CALL info_print ()                      ! Print the information to stdout
  !--------------------------------------------------------------------------------------------------------

  undef_int = -500
  undefined = -999.0_wp ! undef vlaue

  namelist_grid_def = 'INPUT_grid_org'
  CALL init_target_grid(namelist_grid_def)

#ifdef DEBUG
  PRINT *,' target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
  PRINT *,' MAXVAL(lat_geo): ', MAXVAL(lat_geo)
  PRINT *,' MINVAL(lat_geo): ', MINVAL(lat_geo)
#endif
  
  namelist_alb_data_input = 'INPUT_ALB'

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
    PRINT *,'MODIS albedo data (south -> north storage direction)'
  ELSE
    PRINT *,'albedo data (north -> south storage direction)'
  ENDIF
  
  path_alb_file = TRIM(raw_data_alb_path)//TRIM(raw_data_alb_filename)

  IF (ialb_type /= 2) THEN
    path_alnid_file = TRIM(raw_data_alb_path)//TRIM(raw_data_alnid_filename)
    path_aluvd_file = TRIM(raw_data_alb_path)//TRIM(raw_data_aluvd_filename)
    alb_source = TRIM(alb_source)
    alnid_source = TRIM(alnid_source)
    aluvd_source = TRIM(aluvd_source)
#ifdef DEBUG  
    PRINT *, 'after reading namelist for input Albedo data, Albedo raw data are in file:'
    PRINT *, TRIM(path_alb_file)
    PRINT *, TRIM(path_alnid_file)
    PRINT *, TRIM(path_aluvd_file)
    PRINT *, ncid_alb
    PRINT *, 'name of buffer file: ', TRIM(alb_buffer_file)
    PRINT *, 'name of output file: ', TRIM(alb_output_file)
#endif
  ENDIF

  ! open netcdf file with albedo data
  CALL open_netcdf_ALB_data(path_alb_file, ncid_alb)
#ifdef DEBUG  
  PRINT *, 'open netcdf albedo file done '
#endif
  !> inquire dimension information for albedo raw data 
  CALL get_dimension_ALB_data(ncid_alb, nlon_alb, nlat_alb, ntime_alb)

#ifdef DEBUG    
  PRINT *, 'after check of dimensions in Albedo raw data file'
  PRINT *, 'nlon_alb, nlat_alb: ',nlon_alb, nlat_alb
  PRINT *, 'ntime_alb: ', ntime_alb
#endif
  ALLOCATE(time(ntime_alb)) ! this array is needed for netcdf output at the end
  DO nmonth=1, ntime_alb
    time(nmonth) = nmonth
  ENDDO

  CALL allocate_raw_alb_fields(nlon_alb,nlat_alb,ntime_alb)   
  CALL allocate_alb_target_fields(tg,ntime_alb,ialb_type)
  
  CALL get_ALB_data_coordinates(ncid_alb,      &
       &                        nlon_alb,      &
       &                        nlat_alb,      &
       &                        startlon_alb,  &
       &                        startlat_alb,  &
       &                        dlon_alb,      &
       &                        dlat_alb,      &
       &                        lon_alb,       &
       &                        lat_alb)

#ifdef DEBUG  
  PRINT *, 'after getting Albedo data coordinates'
  PRINT *,'startlon_alb: ', startlon_alb
  PRINT *,'startlat_alb: ', startlat_alb
  PRINT *,'dlon_alb: ', dlon_alb
  PRINT *,'dlat_alb: ', dlat_alb
  PRINT *,'lon_alb(1) = ',lon_alb(1) 
  PRINT *,'lon_alb(nlon_alb) = ', lon_alb(nlon_alb) 
#endif
  ! put the values of the grid definition in the data structure alb_raw_data_grid (type alb_reg_lonlat_grid)
  alb_raw_data_grid%start_lon_reg= startlon_alb
  alb_raw_data_grid%start_lat_reg= startlat_alb
  alb_raw_data_grid%dlon_reg= dlon_alb
  ! versions from 2.0_3 on use also raw MODIS albedo data running from South to North
  IF (ialb_type == 2) THEN
    alb_raw_data_grid%dlat_reg=  dlat_alb
  ELSE
    alb_raw_data_grid%dlat_reg= -dlat_alb ! albedo raw data rows from North to South
  ENDIF
  alb_raw_data_grid%nlon_reg= nlon_alb
  alb_raw_data_grid%nlat_reg= nlat_alb

  alb_raw_data_grid%end_lon_reg= lon_alb(nlon_alb) 
  alb_raw_data_grid%end_lat_reg= lat_alb(nlat_alb) 
#ifdef DEBUG  
  PRINT *,'alb_raw_data_grid: ',alb_raw_data_grid
#endif
  CALL close_netcdf_ALB_data(ncid_alb)

  ! start aggregation
#ifdef DEBUG  
  PRINT *,'aggregate Albedo data to target grid'
#endif
  IF (ialb_type == 2) THEN
    CALL agg_alb_data_to_target_grid(tg,undef_alb_bs, path_alb_file, 'ALB_DRY',alb_field_mom)
    alb_dry(:,:,:) = alb_field_mom(:,:,:,1)
#ifdef DEBUG  
    PRINT *,'aggregation dry soil albedo done'
#endif
    CALL agg_alb_data_to_target_grid(tg,undef_alb_bs, path_alb_file, 'ALB_SAT',alb_field_mom)
    alb_sat(:,:,:) = alb_field_mom(:,:,:,1)
#ifdef DEBUG  
    PRINT *,'aggregation saturated soil albedo done'
#endif
  ELSE
    CALL agg_alb_data_to_target_grid(tg,undefined, path_alb_file,  alb_source,alb_field_mom)
#ifdef DEBUG  
    PRINT *,'aggregation month_albedo done'
#endif
    IF (ialb_type == 1) THEN
      CALL agg_alb_data_to_target_grid(tg,undefined, path_alnid_file, alnid_source,alnid_field_mom)
#ifdef DEBUG  
      PRINT *,'aggregation month_alnid done'
#endif
      CALL agg_alb_data_to_target_grid(tg,undefined, path_aluvd_file, aluvd_source,aluvd_field_mom)
#ifdef DEBUG  
      PRINT *,'aggregation month_aluvd done'
#endif
    ENDIF
  ENDIF

  !write out data
  filename = TRIM(alb_output_file)
  !  filename = 'alb_extpar_icon.nc'

  PRINT *,' ======= Checking maximal albedo values  ========='
  IF (ialb_type == 1) THEN
    PRINT *, MAXVAL(alb_field_mom), MAXVAL(alnid_field_mom), MAXVAL(aluvd_field_mom)
  ELSE IF (ialb_type == 2) THEN
    PRINT *, MAXVAL(alb_sat), MAXVAL(alb_dry)
  ELSE IF (ialb_type == 3) THEN
    PRINT *, MAXVAL(alb_field_mom)
  ENDIF

  netcdf_filename = TRIM(alb_buffer_file)

#ifdef DEBUG  
  PRINT *,'write out ', TRIM(netcdf_filename)
#endif
  IF (ialb_type == 2) THEN
    CALL write_netcdf_buffer_alb(netcdf_filename,  &
         &                            tg,         &
         &                            ntime_alb, &
         &                            undefined, &
         &                            undef_int,   &
         &                            lon_geo,     &
         &                            lat_geo, &
         &                            alb_dry=alb_dry, &
         &                            alb_sat=alb_sat)
  ELSE IF (ialb_type == 1) THEN
    CALL write_netcdf_buffer_alb(netcdf_filename,  &
         &                            tg,         &
         &                            ntime_alb, &
         &                            undefined, &
         &                            undef_int,   &
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
         &                            undef_int,   &
         &                            lon_geo,     &
         &                            lat_geo, &
         &                            alb_field_mom=alb_field_mom)
  ENDIF


  CALL deallocate_raw_alb_fields
  CALL deallocate_alb_target_fields(ialb_type)


  PRINT *,'============= alb_to_buffer done ==============='


END PROGRAM extpar_albedo_to_buffer
