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
! ------------ ---------- ----
! V3_0         2023/02/05 Andrzej Wyszogrodzki
!  added ECOSG, modifications to GLCC output data structure
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
PROGRAM extpar_landuse_to_buffer

  USE mo_logging
  USE info_extpar,              ONLY: info_print
  USE mo_io_units,              ONLY: filename_max
  USE mo_kind,                  ONLY: wp, i4

  USE mo_target_grid_data,      ONLY: tg,      &
       &                              lon_geo, &
       &                              lat_geo

  USE mo_target_grid_routines,  ONLY: init_target_grid

  USE mo_landuse_routines,      ONLY: get_dimension_glcc_data, &
       &                              read_namelists_extpar_land_use, &
       &                              get_lonlat_glcc_data, &
       &                              get_dimension_glc2000_data,       &
       &                              get_globcover_tiles_grid, &
       &                              get_ecci_tiles_grid, &
       &                              get_dimension_ecci_data, &
       &                              get_lonlat_ecci_data, &
       &                              get_dimension_ecoclimap_data,       &
       &                              get_lonlat_ecoclimap_data, &
       &                              get_dimension_globcover_data, &
       &                              get_lonlat_globcover_data, &
       &                              get_lonlat_glc2000_data, &
       &                              get_dimension_ecosg_data, &
       &                              get_lonlat_ecosg_data

  USE mo_glc2000_data,          ONLY: glc2000_grid, &
       &                              lon_glc2000,  &
       &                              lat_glc2000,  &
       &                              allocate_raw_glc2000_fields,  &
       &                              deallocate_glc2000_fields

  USE mo_glcc_data,             ONLY: glcc_grid, &
       &                              lon_glcc,  &
       &                              lat_glcc,  &
       &                              allocate_raw_glcc_fields,&
       &                              deallocate_glcc_fields

  USE mo_ecosg_data,            ONLY: ecosg_grid, &
       &                              lon_ecosg,  &
       &                              lat_ecosg,  &
       &                              allocate_raw_ecosg_fields,&
       &                              deallocate_ecosg_fields

  USE mo_ecoclimap_data,        ONLY: deallocate_ecoclimap_fields

  USE mo_ecci_data,             ONLY: ecci_grid,                &
    &                                 lon_ecci,                 &
    &                                 lat_ecci,                 &
    &                                 ecci_tiles_grid,          &
    &                                 ntiles_ecci,              &
    &                                 lu_tiles_lon_min_ecci,              &
    &                                 lu_tiles_lon_max_ecci,              &
    &                                 lu_tiles_lat_min_ecci,              &
    &                                 lu_tiles_lat_max_ecci,              &
    &                                 nc_tiles_lu_ecci,                   &
    &                                 allocate_raw_ecci_fields, &
    &                                 allocate_ecci_data,       &
    &                                 fill_ecci_data,           &
    &                                 deallocate_landuse_data_ecci



  USE mo_glc2000_lookup_tables, ONLY: nclass_glc2000

  USE mo_glcc_lookup_tables,    ONLY: nclass_glcc, &
       &                              ilookup_table_glcc

  USE mo_ecosg_lookup_tables,    ONLY: nclass_ecosg

  USE mo_agg_glc2000,           ONLY: agg_glc2000_data_to_target_grid

  USE mo_agg_ecci,              ONLY: agg_ecci_data_to_target_grid

  USE mo_glcc_tg_fields,        ONLY: fr_land_glcc,       &
       &                              glcc_class_fraction,&
       &                              glcc_class_npixel,  &
       &                              glcc_tot_npixel,    &
       &                              ice_glcc,           &
       &                              z0_glcc,            &
       &                              root_glcc,          &
       &                              plcov_mn_glcc,      &
       &                              plcov_mx_glcc,      &
       &                              lai_mn_glcc,        &
       &                              lai_mx_glcc,        &
       &                              rs_min_glcc,        &
       &                              urban_glcc,         &
       &                              for_d_glcc,         &
       &                              for_e_glcc,         &
       &                              emissivity_glcc,    &
       &                              allocate_glcc_target_fields

  USE mo_ecosg_tg_fields,       ONLY: fr_land_ecosg,       &
       &                              ecosg_class_fraction,&
       &                              ecosg_class_npixel,  &
       &                              ecosg_tot_npixel,    &
       &                              ice_ecosg,           &
       &                              z0_ecosg,            &
       &                              root_ecosg,          &
       &                              plcov_mn_ecosg,      &
       &                              plcov_mx_ecosg,      &
       &                              lai_mn_ecosg,        &
       &                              lai_mx_ecosg,        &
       &                              rs_min_ecosg,        &
       &                              urban_ecosg,         &
       &                              for_d_ecosg,         &
       &                              for_e_ecosg,         &
       &                              skinc_ecosg,         &
       &                              emissivity_ecosg,    &
       &                              allocate_ecosg_target_fields

  USE mo_agg_glcc,              ONLY: agg_glcc_data_to_target_grid

  USE mo_agg_ecosg,             ONLY: agg_ecosg_data_to_target_grid

  USE mo_lu_tg_fields,          ONLY: i_lu_globcover, i_lu_glc2000, i_lu_glcc, i_lu_ecosg, &
       &                              i_lu_ecoclimap, i_lu_ecci, &
       &                              allocate_lu_target_fields, allocate_add_lu_fields, &
       &                              fr_land_lu,       &
       &                              ice_lu,           &
       &                              z0_lu,            &
       &                              root_lu,          &
       &                              plcov_mn_lu,      &
       &                              plcov_mx_lu,      &
       &                              lai_mn_lu,        &
       &                              lai_mx_lu,        &
       &                              rs_min_lu,        &
       &                              urban_lu,         &
       &                              for_d_lu,         &
       &                              for_e_lu,         &
       &                              skinc_lu,         &
       &                              emissivity_lu,    &
       &                              lu_class_fraction,&
       &                              lu_class_npixel,  &
       &                              lu_tot_npixel,    &
       &                              lai12_lu,         &
       &                              plcov12_lu,       &
       &                              z012_lu

  USE mo_landuse_output_nc,     ONLY: write_netcdf_buffer_glcc, &
       &                              write_netcdf_buffer_ecosg,  &
       &                              write_netcdf_buffer_ecoclimap, &
       &                              write_netcdf_buffer_lu

  USE mo_globcover_lookup_tables, ONLY: nclass_globcover

  USE mo_ecci_lookup_tables,    ONLY: nclass_ecci

  USE mo_ecoclimap_lookup_tables, ONLY: nclass_ecoclimap

  USE mo_globcover_data,        ONLY: globcover_grid,                &
       &                              lon_globcover,                 &
       &                              lat_globcover,                 &
       &                              globcover_tiles_grid,          &
       &                              ntiles_globcover,              &
       &                              max_tiles_lu,                  &
       &                              lu_tiles_lon_min,              &
       &                              lu_tiles_lon_max,              &
       &                              lu_tiles_lat_min,              &
       &                              lu_tiles_lat_max,              &
       &                              nc_tiles_lu,                   &
       &                              allocate_raw_globcover_fields, &
       &                              allocate_globcover_data,       &
       &                              fill_globcover_data,           &
       &                              deallocate_landuse_data

  USE mo_ecoclimap_data,         ONLY: ecoclimap_grid, &
       &                               lon_ecoclimap,  &
       &                               lat_ecoclimap,  &
       &                               allocate_raw_ecoclimap_fields

  USE mo_agg_globcover,          ONLY: agg_globcover_data_to_target_grid

  USE mo_agg_ecoclimap,          ONLY: agg_ecoclimap_data_to_target_grid

  USE mo_io_utilities,           ONLY: join_path

  IMPLICIT NONE

  CHARACTER(len=filename_max)             :: netcdf_filename, &
       &                                     namelist_grid_def, &
       &                                     input_lu_namelist_file, &
       &                                     raw_data_lu_path, &         !< path to raw data
       &                                     raw_data_lu_filename(1:max_tiles_lu), &  !< filename glc2000 raw data
       &                                     glcc_file(1), &
       &                                     lu_buffer_file, &  !< name for glc2000 buffer file
       &                                     raw_data_glcc_path, &         !< path to raw data
       &                                     raw_data_glcc_filename, &  !< filename glcc raw data
       &                                     glcc_buffer_file, &  !< name for glcc buffer file
       &                                     lu_dataset  !< name of landuse data set

  CHARACTER(len=filename_max), ALLOCATABLE:: lu_file(:)

  INTEGER(KIND=i4)                        :: ntiles_lu, &
       &                                     nlon_globcover, &  !< number of grid elements in zonal direction for globcover data
       &                                     nlat_globcover, &  !< number of grid elements in meridional direction for globcover data
       &                                     nlon_ecci, & !< number of grid elements in zonal direction for ecci data
       &                                     nlat_ecci, & !< number of grid elements in meridional  direction for ecci data
       &                                     nlon_ecoclimap, &  !< number of grid elements in zonal direction for ecoclimap data
       &                                     nlat_ecoclimap, &  !< number of grid elements in meridional direction for ecoclimap data
       &                                     nlon_glc2000, &  !< number of grid elements in zonal direction for glc2000 data
       &                                     nlat_glc2000, &  !< number of grid elements in meridional direction for glc2000 data
       &                                     nlon_glcc, &  !< number of grid elements in zonal direction for glcc data
       &                                     nlat_glcc, &  !< number of grid elements in meridional direction for glcc data
       &                                     nlon_ecosg, &      !< number of grid elements in zonal      direction for ecosg data
       &                                     nlat_ecosg, &      !< number of grid elements in meridional direction for ecosg data
       &                                     i,k, &  !< counter
       &                                     errorcode, &
       &                                     undef_int, &
       &                                     i_landuse_data, &  !<integer switch to choose a land use raw data set
       &                                     ilookup_table_lu, &  !< integer switch to choose a lookup table
       &                                     nclass_lu !< number of land use classes

  REAL (KIND=wp)                          :: undefined, &
       &                                     tg_southern_bound

  LOGICAL                                 :: l_use_glcc   =.FALSE., &
       &                                     l_use_corine = .FALSE., &
       &                                     l_terra_urb  = .FALSE.

  namelist_grid_def      = 'INPUT_grid_org'
  input_lu_namelist_file = 'INPUT_LU'

  CALL initialize_logging("extpar_landuse_to_buffer.log")
  CALL info_print ()

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= landuse_to_buffer ================')
  CALL logging%info( '')
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= init grid and read namelist=======')
  CALL logging%info( '')

  CALL init_target_grid(namelist_grid_def)

  tg_southern_bound=MINVAL(lat_geo) ! get southern boundary of target grid


  CALL read_namelists_extpar_land_use(input_lu_namelist_file,  &
    &                                 i_landuse_data,          &
    &                                 l_use_corine,            &
    &                                 l_terra_urb,             &
    &                                 raw_data_lu_path,        &
    &                                 raw_data_lu_filename,    &
    &                                 ilookup_table_lu,        &
    &                                 lu_buffer_file,          &
    &                                 raw_data_glcc_path,      &
    &                                 raw_data_glcc_filename,  &
    &                                 ilookup_table_glcc,      &
    &                                 glcc_buffer_file)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= allocate fields ==================')
  CALL logging%info( '')

  CALL logging%info('l_use_array_cache=.FALSE. -> can only be used in consistency_check')

  CALL allocate_lu_target_fields(tg, l_use_array_cache=.FALSE.)

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

  DO i = 1,ntiles_globcover
    CALL logging%info( 'GLOBCOVER TILES, LON, LAT (MIN,MAX): ' )
    WRITE(message_text,*)  i, lu_tiles_lon_min(i), lu_tiles_lon_max(i), &
                 lu_tiles_lat_min(i), lu_tiles_lat_max(i)
  END DO
  CALL logging%info( 'MODEL DOMAIN, LON, LAT (MIN,MAX): ' )
  WRITE(message_text,*)  MINVAL(lon_geo), MAXVAL(lon_geo), &
              MINVAL(lat_geo), MAXVAL(lat_geo)
  CALL logging%info(message_text)

  DO i = 1,ntiles_globcover
    IF (lu_tiles_lon_min(i) < MINVAL(lon_geo).AND. &
        lu_tiles_lon_max(i) > MAXVAL(lon_geo).AND. &
        lu_tiles_lat_min(i) < MINVAL(lat_geo).AND. &
        lu_tiles_lat_max(i) > MAXVAL(lat_geo)) THEN
          WRITE(message_text,*)'MODEL DOMAIN COVERED BY GLOBCOVER TILE ',i
          CALL logging%info(message_text)
    END IF
  END DO

    CASE (i_lu_ecci)
      ntiles_lu = ntiles_ecci

      CALL allocate_ecci_data(ntiles_lu)                  ! allocates the data using ntiles
      CALL fill_ecci_data(raw_data_lu_path,     &
                          raw_data_lu_filename, &  ! the allocated vectors need to be filled with the respective value.
                          lu_tiles_lon_min_ecci, &
                          lu_tiles_lon_max_ecci, &
                          lu_tiles_lat_min_ecci, &
                          lu_tiles_lat_max_ecci, &
                          nc_tiles_lu_ecci)
       DO i = 1,ntiles_ecci
        IF (lu_tiles_lon_min_ecci(i) < MINVAL(lon_geo).AND. &
            lu_tiles_lon_max_ecci(i) > MAXVAL(lon_geo).AND. &
            lu_tiles_lat_min_ecci(i) < MINVAL(lat_geo).AND. &
            lu_tiles_lat_max_ecci(i) > MAXVAL(lat_geo)) THEN
       END IF
       END DO
  END SELECT

  ALLOCATE(lu_file(1:ntiles_lu), STAT= errorcode)
  IF(errorcode /= 0) CALL logging%error('Cant allocate lu_file',__FILE__,__LINE__)

  DO k = 1,ntiles_lu
    lu_file(k) = join_path(raw_data_lu_path,raw_data_lu_filename(k))
  END DO

  glcc_file(1) = join_path(raw_data_glcc_path,raw_data_glcc_filename)

  SELECT CASE (i_landuse_data)
    CASE (i_lu_globcover)
      nclass_lu = nclass_globcover
      lu_dataset = 'GLOBCOVER2009'

      CALL get_dimension_globcover_data(nlon_globcover, &
        &                                  nlat_globcover)
      CALL allocate_raw_globcover_fields(nlat_globcover,nlon_globcover)
      CALL allocate_add_lu_fields(tg,nclass_globcover, l_use_array_cache=.FALSE.)
      CALL get_lonlat_globcover_data( &
        &                              nlon_globcover, &
        &                              nlat_globcover, &
        &                              lon_globcover,  &
        &                              lat_globcover,  &
        &                              globcover_grid)

        ! If southern boundary of target grid is south of southern boundary of Globcover data
        ! (Globcover 2009 does not include Antarctica) then also process GLCC data)
      IF (tg_southern_bound < globcover_grid%end_lat_reg) THEN
        l_use_glcc=.TRUE.
        CALL allocate_glcc_target_fields(tg, l_use_array_cache=.FALSE.)
      ENDIF

      CALL get_globcover_tiles_grid(globcover_tiles_grid)

    CASE (i_lu_ecci)
      nclass_lu = nclass_ecci
      lu_dataset = 'ECCI'

      CALL get_dimension_ecci_data(nlon_ecci, &
        &                                  nlat_ecci)
      CALL allocate_raw_ecci_fields(nlat_ecci,nlon_ecci)
      CALL allocate_add_lu_fields(tg,nclass_ecci, l_use_array_cache=.FALSE.)
      CALL get_lonlat_ecci_data( &
        &                              nlon_ecci, &
        &                              nlat_ecci, &
        &                              lon_ecci,  &
        &                              lat_ecci,  &
        &                              ecci_grid)

      CALL get_ecci_tiles_grid(ecci_tiles_grid)



    CASE (i_lu_ecoclimap)
      nclass_lu = nclass_ecoclimap
      lu_dataset = 'ECOCLIMAP'

      CALL get_dimension_ecoclimap_data(nlon_ecoclimap, &
        &                                  nlat_ecoclimap)
      CALL allocate_raw_ecoclimap_fields(nlat_ecoclimap,nlon_ecoclimap)

      CALL allocate_add_lu_fields(tg,nclass_ecoclimap, l_use_array_cache=.FALSE.)

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
      CALL allocate_add_lu_fields(tg,nclass_glc2000, l_use_array_cache=.FALSE.)
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
        CALL allocate_glcc_target_fields(tg, l_use_array_cache=.FALSE.)
      ENDIF


    CASE (i_lu_ecosg)

      nclass_lu = nclass_ecosg
      lu_dataset = 'ECOCLIMAP SG'

      CALL get_dimension_ecosg_data(lu_file, &
        &                           nlon_ecosg, &
        &                           nlat_ecosg)
      CALL allocate_raw_ecosg_fields(nlat_ecosg,nlon_ecosg)
      CALL allocate_add_lu_fields(tg,nclass_ecosg, l_use_array_cache=.FALSE.) ! CHECK ?????
      CALL get_lonlat_ecosg_data(lu_file, &
        &                        nlon_ecosg, &
        &                        nlat_ecosg, &
        &                        lon_ecosg,  &
        &                        lat_ecosg,  &
        &                        ecosg_grid)

      CALL allocate_ecosg_target_fields(tg)


    CASE (i_lu_glcc)
      nclass_lu = nclass_glcc
      lu_dataset = 'GLCC'
    CASE DEFAULT
      WRITE(message_text,*) 'Invalid land use class :',i_landuse_data
      CALL logging%error(message_text,__FILE__,__LINE__)
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

  CALL logging%info( '')
  CALL logging%info('============= start aggregation ================')
  CALL logging%info( '')

  undefined = 0.0_wp
  SELECT CASE (i_landuse_data)
    CASE(i_lu_globcover)

      CALL agg_globcover_data_to_target_grid(lu_file,                &
           &                                 ilookup_table_lu,     &
           &                                 l_use_corine,         &
           &                                 undefined,            &
           &                                 globcover_tiles_grid, &
           &                                 tg,                   &
           &                                 nclass_globcover,     &
           &                                 lu_class_fraction,    &
           &                                 lu_class_npixel,      &
           &                                 lu_tot_npixel,        &
           &                                 fr_land_lu ,          &
           &                                 ice_lu,               &
           &                                 z0_lu,                &
           &                                 root_lu,              &
           &                                 plcov_mn_lu,          &
           &                                 plcov_mx_lu,          &
           &                                 lai_mn_lu,            &
           &                                 lai_mx_lu,            &
           &                                 rs_min_lu,            &
           &                                 urban_lu,             &
           &                                 for_d_lu,             &
           &                                 for_e_lu,             &
           &                                 skinc_lu,             &
           &                                 emissivity_lu    )

    CASE(i_lu_ecci)

    CALL agg_ecci_data_to_target_grid(lu_file,                &
    &                                        ilookup_table_lu,     &
    &                                        undefined,            &
    &                                        ecci_tiles_grid, &
    &                                        tg,                   &
    &                                        nclass_ecci,     &
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

      CALL agg_ecoclimap_data_to_target_grid(raw_data_lu_path, lu_file,ilookup_table_lu,undefined,       &
      &                                       tg,                                         &
      &                                       nclass_ecoclimap,                             &
      &                                       lu_class_fraction, &
      &                                       lu_class_npixel, &
      &                                       lu_tot_npixel,   &
      &                                       fr_land_lu ,     &
      &                                       ice_lu,          &
      &                                       z012_lu, &
      &                                       root_lu, &
      &                                       plcov12_lu, &
      &                                       lai12_lu,   &
      &                                       rs_min_lu, &
      &                                       urban_lu,  &
      &                                       for_d_lu,  &
      &                                       for_e_lu, &
      &                                       emissivity_lu )

    CASE(i_lu_glc2000)

      CALL agg_glc2000_data_to_target_grid(lu_file,ilookup_table_lu,undefined,       &
      &                                       tg,                                         &
      &                                       nclass_glc2000,                             &
      &                                       lu_class_fraction, &
      &                                       lu_class_npixel, &
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

    CASE(i_lu_ecosg)

      CALL agg_ecosg_data_to_target_grid(lu_file,ilookup_table_lu,undefined,       &
      &                                        tg,                                         &
      &                                        nclass_ecosg,                             &
      &                                        ecosg_class_fraction, &
      &                                        ecosg_class_npixel, &
      &                                        ecosg_tot_npixel,   &
      &                                        fr_land_ecosg ,     &
      &                                        ice_ecosg,          &
      &                                        z0_ecosg, &
      &                                        root_ecosg, &
      &                                        plcov_mn_ecosg, &
      &                                        plcov_mx_ecosg, &
      &                                        lai_mn_ecosg,   &
      &                                        lai_mx_ecosg, &
      &                                        rs_min_ecosg, &
      &                                        urban_ecosg,  &
      &                                        for_d_ecosg,  &
      &                                        for_e_ecosg, &
      &                                        skinc_ecosg, &
      &                                        emissivity_ecosg    )

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

  CALL logging%info( '')
  CALL logging%info('============= write data to netcdf==============')
  CALL logging%info( '')

  ! output
  undefined = -999.0_wp
  undef_int = -999

  netcdf_filename = TRIM(lu_buffer_file)

  SELECT CASE (i_landuse_data)

    CASE(i_lu_glc2000, i_lu_globcover, i_lu_ecci)

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

    CASE(i_lu_ecosg)

      CALL write_netcdf_buffer_ecosg(TRIM(netcdf_filename),  &
        &                                     tg,         &
        &                                     undefined, &
        &                                     undef_int,   &
        &                                     lon_geo,     &
        &                                     lat_geo, &
        &                                     fr_land_ecosg, &
        &                                     ecosg_class_fraction,    &
        &                                     ecosg_class_npixel, &
        &                                     ecosg_tot_npixel, &
        &                                     ice_ecosg, &
        &                                     z0_ecosg, &
        &                                     root_ecosg, &
        &                                     plcov_mn_ecosg, &
        &                                     plcov_mx_ecosg, &
        &                                     lai_mn_ecosg, &
        &                                     lai_mx_ecosg, &
        &                                     rs_min_ecosg, &
        &                                     urban_ecosg,  &
        &                                     for_d_ecosg,  &
        &                                     for_e_ecosg, &
        &                                     skinc_ecosg, &
        &                                     emissivity_ecosg)

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

  CALL logging%info( '')
  CALL logging%info('============= deallocate fields =================')
  CALL logging%info( '')

  SELECT CASE (i_landuse_data)
    CASE(i_lu_globcover)
      CALL deallocate_landuse_data()

    CASE(i_lu_ecci)
      CALL deallocate_landuse_data_ecci()

   CASE(i_lu_ecoclimap)
     CALL deallocate_ecoclimap_fields()

   CASE(i_lu_glc2000)
      CALL  deallocate_glc2000_fields()

   CASE(i_lu_glcc)
     CALL deallocate_glcc_fields()

   CASE(i_lu_ecosg)
     CALL deallocate_ecosg_fields()

  END SELECT

  IF (l_use_glcc) THEN
    CALL deallocate_glcc_fields()
  ENDIF

  CALL logging%info( '')
  CALL logging%info('============= landuse_to_buffer done ============')

END PROGRAM extpar_landuse_to_buffer
