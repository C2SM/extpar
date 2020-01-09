!+ Fortran main program to aggregate land use data to a target grid 
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
! V1_3         2011/04/19 Hermann Asensio
! introduce Globcover 2009 land use data set for external parameters
! V1_7         2013/01/25 Guenther Zaengl 
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON 
! V2_0         2013/06/04 Martina Messmer
!   introduce a new reading routine of the Globcover data set
!   (available as 6 tiles)
!   routines are adapted from the topography
! V2_0_3       2014/09/17 Burkhardt Rockel
!  Added use of directory information to access raw data files
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran main program to aggregate land use data to a target grid
!!
!! @par extpar_landuse_to_buffer 
!!
!! 
!! This program aggregates the GLC2000 land use data and the GLCC data to a given target grid (COSMO/ICON).
!! The desired external parameters are mapped from the GLC2000 land use classes with look-up tables 
!! (see module mo_glc2000_data) and avereaged to the target grid cell (see module mo_agg_glc2000 for details),
!! this is also done for the GLCC data (which cover the whole earth, while GLC2000 data does not cover Antarctica.
!! 
!!
!! @author
!!     Hermann Asensio
!!     (DWD)
!!
!!
!!
PROGRAM extpar_landuse_to_buffer
  
  USE info_extpar, ONLY: info_print
  USE mo_logging

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  USE mo_target_grid_data, ONLY: tg,      &
       &                         lon_geo, &
       &                         lat_geo

  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_io_units,          ONLY: filename_max

  USE mo_utilities_extpar,  ONLY: abort_extpar

  USE mo_landuse_routines, ONLY: read_namelists_extpar_land_use

  USE mo_landuse_routines, ONLY:  get_dimension_glcc_data, &
    &                             get_lonlat_glcc_data, &
    &                             get_dimension_glc2000_data,       &
    &                             get_lonlat_glc2000_data
! >mes
  USE mo_landuse_routines, ONLY:  get_globcover_tiles_grid
! <mes

 USE mo_glc2000_data, ONLY: glc2000_grid, &
    &                       lon_glc2000,  &
    &                       lat_glc2000,  &
    &                       allocate_raw_glc2000_fields,  &
    &                       deallocate_glc2000_fields

  USE mo_glcc_data, ONLY: glcc_grid, &
 &                        lon_glcc,  &
 &                        lat_glcc,  &
 &                        allocate_raw_glcc_fields,&
 &                        deallocate_glcc_fields

  USE mo_ecoclimap_data, ONLY: deallocate_ecoclimap_fields

  USE mo_glc2000_lookup_tables, ONLY: nclass_glc2000

  USE mo_glcc_lookup_tables, ONLY: nclass_glcc
  USE mo_glcc_lookup_tables, ONLY: ilookup_table_glcc

  USE mo_glcc_tg_fields, ONLY:  fr_land_glcc,       &
      &                         glcc_class_fraction,&
      &                         glcc_class_npixel,  &
      &                         glcc_tot_npixel,    &
      &                         ice_glcc,           & 
      &                         z0_glcc,            &
      &                         root_glcc,          &
      &                         plcov_mn_glcc,      &
      &                         plcov_mx_glcc,      &
      &                         lai_mn_glcc,        &
      &                         lai_mx_glcc,        &
      &                         rs_min_glcc,        &
      &                         urban_glcc,         &
      &                         for_d_glcc,         &
      &                         for_e_glcc,         &
      &                         emissivity_glcc,    &
      &                         allocate_glcc_target_fields

  USE mo_agg_glc2000, ONLY : agg_glc2000_data_to_target_grid

  USE mo_agg_glcc, ONLY : agg_glcc_data_to_target_grid

  USE mo_lu_tg_fields, ONLY :  i_lu_globcover, i_lu_glc2000, i_lu_glcc
  USE mo_lu_tg_fields, ONLY :  i_lu_ecoclimap
  USE mo_lu_tg_fields, ONLY: allocate_lu_target_fields, allocate_add_lu_fields
  USE mo_lu_tg_fields, ONLY: fr_land_lu,       &
  &                          ice_lu,           &
  &                          z0_lu,            &
  &                          root_lu,          &
  &                          plcov_mn_lu,      &
  &                          plcov_mx_lu,      &
  &                          lai_mn_lu,        &
  &                          lai_mx_lu,        &
  &                          rs_min_lu,        &
  &                          urban_lu,         &
  &                          for_d_lu,         &
  &                          for_e_lu,         &
  &                          skinc_lu,         &
  &                          emissivity_lu,    &
  &                          lu_class_fraction,&
  &                          lu_class_npixel,  &
  &                          lu_tot_npixel,    &
  &                          lai12_lu,         &
  &                          plcov12_lu,       & 
  &                          z012_lu

  USE mo_landuse_output_nc, ONLY: write_netcdf_buffer_glcc

  USE mo_landuse_output_nc, ONLY: write_netcdf_buffer_ecoclimap
  USE mo_landuse_output_nc, ONLY: write_netcdf_buffer_lu
           
  USE mo_globcover_lookup_tables, ONLY: nclass_globcover

  USE mo_landuse_routines, ONLY: get_dimension_globcover_data, &
    &                            get_lonlat_globcover_data
  USE mo_ecoclimap_lookup_tables, ONLY: nclass_ecoclimap
  USE mo_landuse_routines, ONLY: get_dimension_ecoclimap_data,       &
    &                             get_lonlat_ecoclimap_data

  USE mo_globcover_data, ONLY: globcover_grid,                &
    &                          lon_globcover,                 &
    &                          lat_globcover,                 &
    &                          globcover_tiles_grid,          &
    &                          ntiles_globcover,              &
    &                          max_tiles_lu,                  &
    &                          lu_tiles_lon_min,              &
    &                          lu_tiles_lon_max,              &
    &                          lu_tiles_lat_min,              &
    &                          lu_tiles_lat_max,              &
    &                          nc_tiles_lu,                   &
    &                          allocate_raw_globcover_fields, &
    &                          allocate_globcover_data,       &
    &                          fill_globcover_data,           &
    &                          deallocate_landuse_data

 USE mo_ecoclimap_data, ONLY: ecoclimap_grid, &
    &                         lon_ecoclimap,  &
    &                         lat_ecoclimap,  &
    &                         allocate_raw_ecoclimap_fields

 USE mo_agg_globcover, ONLY : agg_globcover_data_to_target_grid
 USE mo_agg_ecoclimap, ONLY : agg_ecoclimap_data_to_target_grid

  IMPLICIT NONE
  
  CHARACTER(len=filename_max) :: netcdf_filename

  CHARACTER(len=filename_max) :: namelist_grid_def

  CHARACTER(len=filename_max) :: input_lu_namelist_file
  CHARACTER(len=filename_max), ALLOCATABLE:: lu_file(:)

  CHARACTER (len=filename_max) :: raw_data_lu_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_lu_filename(1:max_tiles_lu) !< filename glc2000 raw data
  CHARACTER(len=filename_max) :: glcc_file(1)

  CHARACTER (len=filename_max) :: lu_buffer_file !< name for glc2000 buffer file
  CHARACTER (len=filename_max) :: lu_output_file !< name for glc2000 output file

  CHARACTER (len=filename_max) :: raw_data_glcc_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_glcc_filename !< filename glcc raw data

  CHARACTER (len=filename_max) :: glcc_buffer_file !< name for glcc buffer file
  CHARACTER (len=filename_max) :: glcc_output_file !< name for glcc output file

  CHARACTER(len=filename_max) :: lu_dataset !< name of landuse data set

  INTEGER :: i,k !< counter
  INTEGER :: errorcode

  REAL (KIND=wp) :: undefined
  REAL (KIND=wp) :: tg_southern_bound

  LOGICAL :: l_use_glcc=.FALSE.

  INTEGER :: undef_int

! >mes
  INTEGER(KIND=i4)  :: ntiles_lu
! <mes

  INTEGER (KIND=i8) :: nlon_globcover !< number of grid elements in zonal direction for globcover data
  INTEGER (KIND=i8) :: nlat_globcover !< number of grid elements in meridional direction for globcover data

  INTEGER (KIND=i8) :: nlon_ecoclimap !< number of grid elements in zonal direction for ecoclimap data
  INTEGER (KIND=i8) :: nlat_ecoclimap !< number of grid elements in meridional direction for ecoclimap data

  INTEGER (KIND=i8) :: nlon_glc2000 !< number of grid elements in zonal direction for glc2000 data
  INTEGER (KIND=i8) :: nlat_glc2000 !< number of grid elements in meridional direction for glc2000 data

  INTEGER (KIND=i8) :: nlon_glcc !< number of grid elements in zonal direction for glcc data
  INTEGER (KIND=i8) :: nlat_glcc !< number of grid elements in meridional direction for glcc data

  !--------------------------------------------------------------------------------------

  INTEGER  :: i_landuse_data !<integer switch to choose a land use raw data set
  INTEGER  :: ilookup_table_lu !< integer switch to choose a lookup table
  INTEGER  :: nclass_lu !< number of land use classes 

  CALL initialize_logging("extpar_landuse_to_buffer.log")
  CALL info_print ()
  !--------------------------------------------------------------------------------------------------------

  namelist_grid_def = 'INPUT_grid_org'
  input_lu_namelist_file = 'INPUT_LU'

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*) '============= landuse_to_buffer ================'
  WRITE(logging%fileunit,*) ''
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*) '============= init grid and read namelist======='
  WRITE(logging%fileunit,*) ''

  CALL init_target_grid(namelist_grid_def)

  tg_southern_bound=MINVAL(lat_geo) ! get southern boundary of target grid

  CALL allocate_lu_target_fields(tg)

  CALL read_namelists_extpar_land_use(input_lu_namelist_file, &
    &                                 i_landuse_data,         &
    &                                 raw_data_lu_path,       &
    &                                 raw_data_lu_filename,   &
    &                                 ilookup_table_lu,       &
    &                                 lu_buffer_file,         &
    &                                 lu_output_file,         &
    &                                 raw_data_glcc_path,     &
    &                                 raw_data_glcc_filename, &
    &                                 ilookup_table_glcc,     &
    &                                 glcc_buffer_file,       &
    &                                 glcc_output_file)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= allocate fields =================='
  WRITE(logging%fileunit,*) ''

  ntiles_lu = 1
  SELECT CASE(i_landuse_data)
    CASE (i_lu_globcover)
      ntiles_lu = ntiles_globcover
       
      CALL allocate_globcover_data(ntiles_lu)                  ! allocates the data using ntiles
      CALL fill_globcover_data(raw_data_lu_path,     &
                          raw_data_lu_filename, &  ! the allocated vectors need to be filled with the respective value.
                          lu_tiles_lon_min, &
                          lu_tiles_lon_max, &    
                          lu_tiles_lat_min, &
                          lu_tiles_lat_max, &
                          nc_tiles_lu)

  IF (verbose >= idbg_high ) THEN
    DO i = 1,ntiles_globcover
      WRITE(logging%fileunit,*) 'GLOBCOVER TILES, LON, LAT (MIN,MAX): ' 
      WRITE(logging%fileunit,*)  i, lu_tiles_lon_min(i), lu_tiles_lon_max(i), &
                   lu_tiles_lat_min(i), lu_tiles_lat_max(i) 
    END DO
    WRITE(logging%fileunit,*) 'MODEL DOMAIN, LON, LAT (MIN,MAX): ' 
    WRITE(logging%fileunit,*)  MINVAL(lon_geo), MAXVAL(lon_geo), &
                MINVAL(lat_geo), MAXVAL(lat_geo)
  ENDIF

  DO i = 1,ntiles_globcover
    IF (lu_tiles_lon_min(i) < MINVAL(lon_geo).AND. &
        lu_tiles_lon_max(i) > MAXVAL(lon_geo).AND. &
        lu_tiles_lat_min(i) < MINVAL(lat_geo).AND. &
        lu_tiles_lat_max(i) > MAXVAL(lat_geo)) THEN
      WRITE(logging%fileunit,*)'MODEL DOMAIN COVERED BY GLOBCOVER TILE ',i
    END IF
  END DO

  END SELECT

  ALLOCATE(lu_file(1:ntiles_lu), STAT= errorcode)
  IF(errorcode /= 0) CALL abort_extpar('Cant allocate lu_file')
  IF (verbose >= idbg_low ) THEN
    WRITE(logging%fileunit,*)'ntiles_lu: ', ntiles_lu
    WRITE(logging%fileunit,*)'raw_data_glcc_filename: ',TRIM(raw_data_glcc_filename)
  ENDIF

  DO k = 1,ntiles_lu
    lu_file(k) = TRIM(raw_data_lu_path) // TRIM(raw_data_lu_filename(k))
    IF (verbose >= idbg_high ) WRITE(logging%fileunit,*)'lu_file: ', TRIM(lu_file(k))
  END DO

  glcc_file(1) = TRIM(raw_data_glcc_path) // TRIM(raw_data_glcc_filename)
  IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'glcc file: ', TRIM(glcc_file(1))

  SELECT CASE (i_landuse_data)
    CASE (i_lu_globcover)
      nclass_lu = nclass_globcover
      lu_dataset = 'GLOBCOVER2009'

      CALL get_dimension_globcover_data(nlon_globcover, &
        &                                  nlat_globcover)
      CALL allocate_raw_globcover_fields(nlat_globcover,nlon_globcover)
      CALL allocate_add_lu_fields(tg,nclass_globcover)
      CALL get_lonlat_globcover_data( &
        &                              nlon_globcover, &
        &                              nlat_globcover, &
        &                              lon_globcover,  &
        &                              lat_globcover,  &
        &                              globcover_grid)
        IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'globcover_grid: ',globcover_grid

        ! If southern boundary of target grid is south of southern boundary of Globcover data
        ! (Globcover 2009 does not include Antarctica) then also process GLCC data)
        IF (tg_southern_bound < globcover_grid%end_lat_reg) THEN
          l_use_glcc=.TRUE.
          CALL allocate_glcc_target_fields(tg)
        ENDIF 

      CALL get_globcover_tiles_grid(globcover_tiles_grid)
      WRITE(logging%fileunit,*)'globcover_tiles_grid(1): ', globcover_tiles_grid(1)

    CASE (i_lu_ecoclimap)
      nclass_lu = nclass_ecoclimap
      lu_dataset = 'ECOCLIMAP'

      CALL get_dimension_ecoclimap_data(nlon_ecoclimap, &
        &                                  nlat_ecoclimap)
      CALL allocate_raw_ecoclimap_fields(nlat_ecoclimap,nlon_ecoclimap)

      CALL allocate_add_lu_fields(tg,nclass_ecoclimap)

      CALL get_lonlat_ecoclimap_data( &
        &                              nlon_ecoclimap, &
        &                              nlat_ecoclimap, &
        &                              lon_ecoclimap,  &
        &                              lat_ecoclimap,  &
        &                              ecoclimap_grid)

    CASE (i_lu_glc2000)
      nclass_lu = nclass_glc2000
      lu_dataset = 'GLC2000'

      CALL get_dimension_glc2000_data(lu_file, &
        &                                  nlon_glc2000, &
        &                                  nlat_glc2000)
      CALL allocate_raw_glc2000_fields(nlat_glc2000,nlon_glc2000)
      CALL allocate_add_lu_fields(tg,nclass_glc2000)
      CALL get_lonlat_glc2000_data(lu_file, &
        &                              nlon_glc2000, &
        &                              nlat_glc2000, &
        &                              lon_glc2000,  &
        &                              lat_glc2000,  &
        &                              glc2000_grid)
      ! If southern boundary of target grid is south of southern boundary of GLC2000 data
      ! (GLC2000 does not include Antarctica) then also process GLCC data)
      IF(tg_southern_bound < glc2000_grid%end_lat_reg) THEN
        l_use_glcc=.TRUE.
        CALL allocate_glcc_target_fields(tg)
      ENDIF 
    CASE (i_lu_glcc)
      nclass_lu = nclass_glcc
      lu_dataset = 'GLCC'
    CASE DEFAULT
      WRITE(logging%fileunit,*) '***ERROR: extpar_landuse_to_buffer'
      WRITE(logging%fileunit,*) '***ERROR: Invalid land use class :',i_landuse_data
      CALL abort_extpar(' Invalid land use class')
  END SELECT

  IF (l_use_glcc.OR.(i_landuse_data==i_lu_glcc)) THEN
     
    CALL get_dimension_glcc_data(glcc_file, &
      &                           nlon_glcc, &
      &                           nlat_glcc)

    CALL allocate_raw_glcc_fields(nlat_glcc, nlon_glcc)

    CALL  get_lonlat_glcc_data(glcc_file, &
     &                                   nlon_glcc, &
     &                                   nlat_glcc, &
     &                                   lon_glcc,  &
     &                                   lat_glcc,  &
     &                                   glcc_grid)
  ENDIF

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= start aggregation ================'
  WRITE(logging%fileunit,*) ''
  
  undefined = 0.0_wp
  SELECT CASE (i_landuse_data)
    CASE(i_lu_globcover)

    CALL agg_globcover_data_to_target_grid(lu_file,                &  
    &                                        ilookup_table_lu,     &
    &                                        undefined,            &
    &                                        globcover_tiles_grid, &
    &                                        tg,                   &
    &                                        nclass_globcover,     &
    &                                        lu_class_fraction,    &
    &                                        lu_class_npixel,      &
    &                                        lu_tot_npixel,        &
    &                                        fr_land_lu ,          &
    &                                        ice_lu,               &
    &                                        z0_lu,                &
    &                                        root_lu,              &
    &                                        plcov_mn_lu,          &
    &                                        plcov_mx_lu,          &
    &                                        lai_mn_lu,            &
    &                                        lai_mx_lu,            &
    &                                        rs_min_lu,            &
    &                                        urban_lu,             &
    &                                        for_d_lu,             &
    &                                        for_e_lu,             &
    &                                        skinc_lu,             &
    &                                        emissivity_lu    )


   CASE(i_lu_ecoclimap)

!_br 17.09.14    CALL agg_ecoclimap_data_to_target_grid(lu_file,ilookup_table_lu,undefined,       &
    CALL agg_ecoclimap_data_to_target_grid(raw_data_lu_path, lu_file,ilookup_table_lu,undefined,       &
    &                                        tg,                                         &
    &                                        nclass_ecoclimap,                             &
    &                                        lu_class_fraction, &
    &                                        lu_class_npixel, &
    &                                        lu_tot_npixel,   &
    &                                        fr_land_lu ,     &
    &                                        ice_lu,          &
    &                                        z012_lu, &
    &                                        root_lu, &
    &                                        plcov12_lu, &
    &                                        lai12_lu,   &
    &                                        rs_min_lu, &
    &                                        urban_lu,  &
    &                                        for_d_lu,  &
    &                                        for_e_lu, &
    &                                        emissivity_lu )

    CASE(i_lu_glc2000)

    CALL agg_glc2000_data_to_target_grid(lu_file,ilookup_table_lu,undefined,       &
    &                                        tg,                                         &
    &                                        nclass_glc2000,                             &
    &                                        lu_class_fraction, &
    &                                        lu_class_npixel, &
    &                                        lu_tot_npixel,   &
    &                                        fr_land_lu ,     &
    &                                        ice_lu,          &
    &                                        z0_lu, &
    &                                        root_lu, &
    &                                        plcov_mn_lu, &
    &                                        plcov_mx_lu, &
    &                                        lai_mn_lu,   &
    &                                        lai_mx_lu, &
    &                                        rs_min_lu, &
    &                                        urban_lu,  &
    &                                        for_d_lu,  &
    &                                        for_e_lu, &
    &                                        emissivity_lu    )

   CASE(i_lu_glcc)

    CALL agg_glcc_data_to_target_grid(lu_file,ilookup_table_lu,undefined,       &
    &                                        tg,                                         &
    &                                        nclass_glcc,                             &
    &                                        lu_class_fraction, &
    &                                        lu_class_npixel, &
    &                                        lu_tot_npixel,   &
    &                                        fr_land_lu ,     &
    &                                        ice_lu,          &
    &                                        z0_lu, &
    &                                        root_lu, &
    &                                        plcov_mn_lu, &
    &                                        plcov_mx_lu, &
    &                                        lai_mn_lu,   &
    &                                        lai_mx_lu, &
    &                                        rs_min_lu, &
    &                                        urban_lu,  &
    &                                        for_d_lu,  &
    &                                        for_e_lu, &
    &                                        emissivity_lu    )

  END SELECT

  IF (l_use_glcc) THEN ! additionally process GLCC data
    CALL agg_glcc_data_to_target_grid(glcc_file,ilookup_table_glcc,undefined,  &
      &                                        tg,                              &
      &                                        nclass_glcc,                     &
      &                                        glcc_class_fraction, &
      &                                        glcc_class_npixel, &
      &                                        glcc_tot_npixel,   &
      &                                        fr_land_glcc ,     &
      &                                        ice_glcc,          &
      &                                        z0_glcc, &
      &                                        root_glcc, &
      &                                        plcov_mn_glcc, &
      &                                        plcov_mx_glcc, &
      &                                        lai_mn_glcc,   &
      &                                        lai_mx_glcc, &
      &                                        rs_min_glcc, &
      &                                        urban_glcc,  &
      &                                        for_d_glcc,  &
      &                                        for_e_glcc, &
      &                                        emissivity_glcc    )

  ENDIF

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= write data to netcdf=============='
  WRITE(logging%fileunit,*) ''

  ! output
  undefined = -999.0_wp
  undef_int = -999

  netcdf_filename = TRIM(lu_buffer_file)
  IF (verbose >= idbg_low ) WRITE(logging%fileunit,*) 'Land-use buffer filename: ',TRIM(netcdf_filename)

   
  SELECT CASE (i_landuse_data)
  
  CASE(i_lu_glc2000, i_lu_globcover)

  IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)' WRITE ICE_LU MAX: ', MAXVAL(ice_lu)

  CALL write_netcdf_buffer_lu(TRIM(netcdf_filename),  &
    &                          TRIM(lu_dataset), &
    &                                     tg,         &
    &                                     i_landuse_data, &
    &                                     ilookup_table_lu, &
    &                                     nclass_lu, &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     fr_land_lu, &
    &                                     lu_class_fraction,    &
    &                                     lu_class_npixel, &
    &                                     lu_tot_npixel, &
    &                                     ice_lu, &
    &                                     z0_lu, &
    &                                     root_lu, &
    &                                     plcov_mn_lu, &
    &                                     plcov_mx_lu, &
    &                                     lai_mn_lu, &
    &                                     lai_mx_lu, &
    &                                     rs_min_lu, &
    &                                     urban_lu,  &
    &                                     for_d_lu,  &
    &                                     for_e_lu, &
    &                                     skinc_lu, &
    &                                     emissivity_lu)

   IF (l_use_glcc) THEN !
     netcdf_filename = TRIM(glcc_buffer_file)
     IF (verbose >= idbg_low ) WRITE(logging%fileunit,*) 'GLCC buffer filename: ',TRIM(netcdf_filename)
     CALL write_netcdf_buffer_glcc(TRIM(netcdf_filename),  &
      &                                     tg,         &
      &                                     undefined, &
      &                                     undef_int,   &
      &                                     lon_geo,     &
      &                                     lat_geo, &
      &                                     fr_land_glcc, &
      &                                     glcc_class_fraction,    &
      &                                     glcc_class_npixel, &
      &                                     glcc_tot_npixel, &
      &                                     ice_glcc, &
      &                                     z0_glcc, &
      &                                     root_glcc, &
      &                                     plcov_mn_glcc, &
      &                                     plcov_mx_glcc, &
      &                                     lai_mn_glcc, &
      &                                     lai_mx_glcc, &
      &                                     rs_min_glcc, &
      &                                     urban_glcc,  &
      &                                     for_d_glcc,  &
      &                                     for_e_glcc, &
      &                                     emissivity_glcc)
   ENDIF

   CASE(i_lu_ecoclimap)

     netcdf_filename = TRIM(lu_buffer_file)

     CALL write_netcdf_buffer_ecoclimap(TRIM(netcdf_filename),  &
      &                                     tg,         &
      &                                     i_landuse_data, &
      &                                     ilookup_table_lu, &
      &                                     nclass_lu, &
      &                                     undefined, &
      &                                     undef_int,   &
      &                                     lon_geo,     &
      &                                     lat_geo, &
      &                                     fr_land_lu, &
      &                                     lu_class_fraction,    &
      &                                     lu_class_npixel, &
      &                                     lu_tot_npixel, &
      &                                     ice_lu, &
      &                                     z012_lu, &
      &                                     root_lu, &
      &                                     plcov12_lu, &
      &                                     lai12_lu, &
      &                                     rs_min_lu, &
      &                                     urban_lu,  &
      &                                     for_d_lu,  &
      &                                     for_e_lu, &
      &                                     emissivity_lu)
 
  END SELECT 

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= deallocate fields ================='
  WRITE(logging%fileunit,*) ''

  SELECT CASE (i_landuse_data)
    CASE(i_lu_globcover)
      CALL deallocate_landuse_data()

   CASE(i_lu_ecoclimap)
     CALL deallocate_ecoclimap_fields()

   CASE(i_lu_glc2000)
      CALL  deallocate_glc2000_fields()

   CASE(i_lu_glcc)
     CALL deallocate_glcc_fields()

  END SELECT

  IF (l_use_glcc) THEN
    CALL deallocate_glcc_fields()
  ENDIF

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= landuse_to_buffer done ============'

END PROGRAM extpar_landuse_to_buffer

