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

  ! Load the library information data:
  USE info_extpar, ONLY: info_define, info_readnl, info_print


  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4


  USE mo_target_grid_data, ONLY: lon_geo, &
    &                            lat_geo, &
    &                            no_raw_data_pixel

  USE mo_target_grid_data, ONLY: tg  !< structure with target grid description

  USE mo_target_grid_routines, ONLY: init_target_grid


  USE mo_grid_structures, ONLY: target_grid_def,   &
    &                            reg_lonlat_grid,   &
    &                            rotated_lonlat_grid

  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo
  USE mo_grid_structures, ONLY: igrid_gme


  USE  mo_icon_grid_data, ONLY: ICON_grid !< structure which contains the definition of the ICON grid

  USE  mo_cosmo_grid, ONLY: COSMO_grid, &
    &                       lon_rot, &
    &                       lat_rot, &
    &                       allocate_cosmo_rc, &
    &                       get_cosmo_grid_info, &
    &                       calculate_cosmo_domain_coordinates



  USE mo_base_geometry,    ONLY:  geographical_coordinates, &
    &                               cartesian_coordinates
  
  USE mo_icon_domain,          ONLY: icon_domain, &
    &                            grid_cells,               &
    &                            grid_vertices,            &
    &                            construct_icon_domain,    &
                                destruct_icon_domain
  
  USE mo_io_units,          ONLY: filename_max

  USE mo_exception,         ONLY: message_text, message, finish

  USE mo_utilities_extpar, ONLY: abort_extpar

 
  USE mo_additional_geometry,   ONLY: cc2gc,                  &
    &                            gc2cc,                  &
    &                            arc_length,             &
    &                            cos_arc_length,         &
    &                            inter_section,          &
    &                            vector_product,         &
    &                            point_in_polygon_sp

  USE mo_math_constants,  ONLY: pi, pi_2, dbl_eps,rad2deg

  USE mo_albedo_routines, ONLY: read_namelists_extpar_alb

  USE mo_albedo_data, ONLY: alb_raw_data_grid, &
    &                           alb_field_row_mom, &
    &                           alnid_field_row_mom, &
    &                           aluvd_field_row_mom, &
    &                           alb_field_row, &
    &                           alnid_field_row, &
    &                           aluvd_field_row, &
    &                           lon_alb, &
    &                           lat_alb, &
    &                           ntime_alb, &
    &                           ialb_type, &
    &                           allocate_raw_alb_fields, &
    &                           deallocate_raw_alb_fields
                               
  USE mo_albedo_tg_fields, ONLY: alb_field, alb_dry, alb_sat, &
    &                            alnid_field, &
    &                            aluvd_field, &
    &                            alb_field_mom, &
    &                            alnid_field_mom, &
    &                            aluvd_field_mom, &
    &                            allocate_alb_target_fields, &
    &                            deallocate_alb_target_fields

  USE mo_albedo_data, ONLY: undef_alb_bs

  USE mo_albedo_routines, ONLY: open_netcdf_ALB_data, &
    &                               close_netcdf_ALB_data, &
    &                               read_alb_data_input_namelist, &
    &                               get_dimension_ALB_data, &
    &                               get_ALB_data_coordinates
                                   

  USE mo_agg_albedo, ONLY: agg_alb_data_to_target_grid

  USE mo_albedo_output_nc, ONLY: write_netcdf_buffer_alb
  USE mo_albedo_output_nc, ONLY: write_netcdf_cosmo_grid_alb
  USE mo_albedo_output_nc, ONLY: write_netcdf_icon_grid_alb

  IMPLICIT NONE

  CHARACTER(len=filename_max) :: namelist_grid_def

  
  CHARACTER (len=filename_max) :: input_namelist_cosmo_grid !< file with input namelist with COSMO grid definition
  CHARACTER (len=filename_max) :: namelist_alb_data_input !< file with input namelist with albedo data information

  CHARACTER (len=filename_max) :: raw_data_path        !< path to albedo data
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


  INTEGER (KIND=i4):: nrow    !< index for number of data row of albedo data
  INTEGER (KIND=i4):: ncolumn !< index for number of data column of albedo data
  INTEGER (KIND=i4):: nmonth  !< index for month for albedo data

  REAL (KIND=wp) :: dlon_alb !< grid point distance in zonal direction (in degrees) for albedo data
  REAL (KIND=wp) :: dlat_alb !< grid point distance in meridional direction (in degrees) for albedo data

  REAL (KIND=wp) :: startlon_alb !< longitude of lower left grid element for albedo data 

  REAL (KIND=wp) :: startlat_alb !< latitude of lower left grid element for albedo data

  INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid

  REAL(KIND=wp) :: undefined !< value to indicate undefined grid elements 
  INTEGER (KIND=i4) :: undef_int   !< value for undefined integer
  INTEGER (KIND=i4) :: default_value !< default value


 ! Print the default information to stdout:
  CALL info_define ('alb_to_buffer')      ! Pre-define the program name as binary name
  CALL info_print ()                     ! Print the information to stdout
  !--------------------------------------------------------------------------------------------------------
  undef_int = 0 ! set undefined to zero
  undefined = -999.0 ! undef vlaue
  default_value = 0. ! default value
      
  namelist_grid_def = 'INPUT_grid_org'
  CALL init_target_grid(namelist_grid_def)

  PRINT *,' target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
  !HA debug:
    print *,' MAXVAL(lat_geo): ', MAXVAL(lat_geo)
    print *,' MINVAL(lat_geo): ', MINVAL(lat_geo)


  igrid_type = tg%igrid_type
  ! get information on target grid

  ! read namelist for input albedo data

  namelist_alb_data_input = 'INPUT_ALB'

  CALL  read_namelists_extpar_alb(namelist_alb_data_input,     &
    &                                  raw_data_alb_path,      &
    &                                  raw_data_alb_filename,  &
    &                                  raw_data_alnid_filename,&
    &                                  raw_data_aluvd_filename,&
    &                                  ialb_type,              &
    &                                  alb_buffer_file,        &
    &                                  alb_output_file,        &
    &                                  alb_source,             &
    &                                  alnid_source,           &
    &                                  aluvd_source)

   
  path_alb_file = TRIM(raw_data_alb_path)//TRIM(raw_data_alb_filename)

  IF (ialb_type /= 2) THEN
    path_alnid_file = TRIM(raw_data_alb_path)//TRIM(raw_data_alnid_filename)
    path_aluvd_file = TRIM(raw_data_alb_path)//TRIM(raw_data_aluvd_filename)
    alb_source = TRIM(alb_source)
    alnid_source = TRIM(alnid_source)
    aluvd_source = TRIM(aluvd_source)
  !  print *, 'after reading namelist for input Albedo data, Albedo raw data are in file:'
  !  print *, TRIM(path_alb_file)
  !  print *, TRIM(path_alnid_file)
  !  print *, TRIM(path_aluvd_file)
  !  print *, ncid_alb
  !  print *, 'name of buffer file: ', TRIM(alb_buffer_file)
  !  print *, 'name of output file: ', TRIM(alb_output_file)
  ENDIF
       
  ! open netcdf file with albedo data
  CALL open_netcdf_ALB_data(path_alb_file, &
    &                           ncid_alb)
  PRINT *, 'open netcdf albedo file done '

  !> inquire dimension information for albedo raw data 
  CALL get_dimension_ALB_data(ncid_alb,       &  
    &                                nlon_alb, &
    &                                nlat_alb, &
                                     ntime_alb)

  !HA debug
  !  print *, 'after check of dimensions in Albedo raw data file'
  !  print *, 'nlon_alb, nlat_alb: ',nlon_alb, nlat_alb
  !  print *, 'ntime_alb: ', ntime_alb

  ALLOCATE(time(ntime_alb)) ! this array is needed for netcdf output at the end
  DO nmonth=1, ntime_alb
    time(nmonth) = nmonth
  ENDDO


  CALL allocate_raw_alb_fields(nlon_alb,nlat_alb,ntime_alb)   
  CALL allocate_alb_target_fields(tg,ntime_alb,ialb_type)

  CALL get_ALB_data_coordinates(ncid_alb,      &
    &                               nlon_alb,      &
    &                               nlat_alb,      &
    &                               startlon_alb,  &
    &                               startlat_alb,  &
    &                               dlon_alb,      &
    &                               dlat_alb,      &
    &                               lon_alb,       &
    &                               lat_alb)



  !HA debug
  !  print *, 'after getting Albedo data coordinates'
  !  print *,'startlon_alb: ', startlon_alb
  !  print *,'startlat_alb: ', startlat_alb
  !  print *,'dlon_alb: ', dlon_alb
  !  print *,'dlat_alb: ', dlat_alb
  !  print *,'lon_alb(1) = ',lon_alb(1) 
  !  print *,'lon_alb(nlon_alb) = ', lon_alb(nlon_alb) 
  ! put the values of the grid definition in the data structure alb_raw_data_grid (type alb_reg_lonlat_grid)
  alb_raw_data_grid%start_lon_reg= startlon_alb
  alb_raw_data_grid%start_lat_reg= startlat_alb
  alb_raw_data_grid%dlon_reg= dlon_alb
  alb_raw_data_grid%dlat_reg= dlat_alb ! albedo raw data rows from South to North
  alb_raw_data_grid%nlon_reg= nlon_alb
  alb_raw_data_grid%nlat_reg= nlat_alb

  alb_raw_data_grid%end_lon_reg= lon_alb(nlon_alb) 
  alb_raw_data_grid%end_lat_reg= lat_alb(nlat_alb) 
  PRINT *,'alb_raw_data_grid: ',alb_raw_data_grid

  CALL close_netcdf_ALB_data(ncid_alb)

  ! start aggregation
  PRINT *,'aggregate Albedo data to target grid'

  IF (ialb_type == 2) THEN
    CALL agg_alb_data_to_target_grid(tg,undef_alb_bs, path_alb_file, &
         &                     'ALB_DRY',alb_field_mom)
    alb_dry(:,:,:) = alb_field_mom(:,:,:,1)
    PRINT *,'aggregation dry soil albedo done'
    CALL agg_alb_data_to_target_grid(tg,undef_alb_bs, path_alb_file, &
         &                     'ALB_SAT',alb_field_mom)
    alb_sat(:,:,:) = alb_field_mom(:,:,:,1)
    PRINT *,'aggregation saturated soil albedo done'
  ELSE
    CALL agg_alb_data_to_target_grid(tg,undefined, path_alb_file, &
         &                     alb_source,alb_field_mom)
    PRINT *,'aggregation month_albedo done'
    IF (ialb_type == 1) THEN
      CALL agg_alb_data_to_target_grid(tg,undefined, path_alnid_file, &
           &                     alnid_source,alnid_field_mom)
      PRINT *,'aggregation month_alnid done'
      CALL agg_alb_data_to_target_grid(tg,undefined, path_aluvd_file, &
           &                     aluvd_source,aluvd_field_mom)
      PRINT *,'aggregation month_aluvd done'
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

  SELECT CASE(igrid_type)  
  CASE(igrid_icon) ! ICON GRID

      netcdf_filename = TRIM(alb_output_file)
 !     netcdf_filename = filename
      undefined = -500.
      undef_int = -500

      PRINT *,'write out ', TRIM(netcdf_filename)

      CALL write_netcdf_icon_grid_alb(netcdf_filename,  &
   &                                     icon_grid,         &
   &                                     tg,         &
   &                                     ntime_alb, &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     alb_field_mom, &
   &                                     alnid_field_mom, &
   &                                     aluvd_field_mom)

    CASE(igrid_cosmo) ! COSMO grid
    
      netcdf_filename = TRIM(alb_output_file)
      undefined = -500.
      undef_int = -500

      PRINT *,'write out ', TRIM(netcdf_filename)


      CALL write_netcdf_cosmo_grid_alb(netcdf_filename,  &
   &                                     cosmo_grid,         &
   &                                     tg,         &
   &                                     ntime_alb, &
   &                                     undefined, &
   &                                     undef_int,   &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     alb_field_mom, &
   &                                     alnid_field_mom, &
   &                                     aluvd_field_mom)



    CASE(igrid_gme) ! GME grid   

  END SELECT

  netcdf_filename = TRIM(alb_buffer_file)
!  netcdf_filename = 'alb_buffer.nc'
  undefined = -500.
  undef_int = -500

  PRINT *,'write out ', TRIM(netcdf_filename)

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
