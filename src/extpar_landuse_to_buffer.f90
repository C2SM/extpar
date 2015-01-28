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
  
  ! Load the library information data:
  USE info_extpar, ONLY: info_define, info_readnl, info_print


  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  USE mo_grid_structures, ONLY: target_grid_def,   &
    &                            reg_lonlat_grid,   &
    &                            rotated_lonlat_grid
  
  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo
  USE mo_grid_structures, ONLY: igrid_gme

  USE mo_target_grid_data, ONLY: lon_geo, &
    &                            lat_geo, &
    &                            no_raw_data_pixel, &
    &                            allocate_com_target_fields
  
  USE mo_target_grid_data, ONLY: tg
  
  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_icon_grid_data, ONLY: ICON_grid  !< structure which contains the definition of the ICON grid

  USE  mo_cosmo_grid, ONLY: COSMO_grid, &
    &                       lon_rot, &
    &                       lat_rot, &
    &                       allocate_cosmo_rc, &
    &                       get_cosmo_grid_info, &
    &                       calculate_cosmo_domain_coordinates

  USE mo_base_geometry,    ONLY:  geographical_coordinates, &
    &                             cartesian_coordinates

  USE mo_icon_domain,          ONLY: icon_domain, &
    &                             grid_cells,               &
    &                             grid_vertices,            &
    &                             construct_icon_domain,    &
    &                             destruct_icon_domain

  USE mo_io_units,          ONLY: filename_max

  USE mo_exception,         ONLY: message_text, message, finish

  USE mo_utilities_extpar,  ONLY: abort_extpar

  USE mo_additional_geometry,   ONLY: cc2gc,                  &
    &                                 gc2cc,                  &
    &                                 arc_length,             &
    &                                 cos_arc_length,         &
    &                                 inter_section,          &
    &                                 vector_product,         &
    &                                 point_in_polygon_sp


  USE mo_math_constants,  ONLY: pi, pi_2, dbl_eps,rad2deg

  USE mo_landuse_routines, ONLY: read_namelists_extpar_land_use

  USE mo_landuse_routines, ONLY:  get_dimension_glcc_data, &
    &                             get_lonlat_glcc_data, &
    &                             get_dimension_glc2000_data,       &
    &                             get_lonlat_glc2000_data
! >mes
  USE mo_landuse_routines, ONLY:  get_globcover_tiles_grid,  &
    &                             det_band_globcover_data
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

  USE mo_glc2000_lookup_tables, ONLY: glc2000_legend
  USE mo_glc2000_lookup_tables, ONLY: nclass_glc2000
  USE mo_glc2000_lookup_tables, ONLY: i_gme_lookup_table,   &
    &                                 i_cosmo_lookup_table, &
    &                                 i_experimental_lookup_table
  USE mo_glc2000_lookup_tables, ONLY: init_glc2000_lookup_tables

  USE mo_glc2000_lookup_tables, ONLY:   z0_lt_glc2000, lnz0_lt_glc2000, plc_mn_lt_glc2000, plc_mx_lt_glc2000, & 
    &               lai_mn_lt_glc2000, lai_mx_lt_glc2000, rd_lt_glc2000, emiss_lt_glc2000, rs_min_lt_glc2000   

USE mo_glcc_lookup_tables, ONLY: nclass_glcc
USE mo_glcc_lookup_tables, ONLY: init_glcc_lookup_tables
USE mo_glcc_lookup_tables, ONLY: glcc_legend
USE mo_glcc_lookup_tables, ONLY: ilookup_table_glcc, &
  &                              i_gme_lookup_table_glcc, &
  &                              i_cosmo_lookup_table_glcc, &
  &                              i_experimental_lookup_table_glcc
USE mo_glcc_lookup_tables, ONLY: z0_lt_glcc, lnz0_lt_glcc, plc_mn_lt_glcc, plc_mx_lt_glcc

USE mo_glcc_lookup_tables, ONLY: lai_mn_lt_glcc, lai_mx_lt_glcc, rd_lt_glcc, emiss_lt_glcc, rs_min_lt_glcc         

  USE mo_glc2000_tg_fields, ONLY: allocate_glc2000_target_fields

  USE mo_glc2000_tg_fields, ONLY: fr_land_glc2000,       &
    &                             glc2000_class_fraction,&
    &                             glc2000_class_npixel,  &
    &                             glc2000_tot_npixel,    &
    &                             ice_glc2000,           &
    &                             z0_glc2000,            &
    &                             root_glc2000,          &
    &                             plcov_mn_glc2000,      &
    &                             plcov_mx_glc2000,      &
    &                             lai_mn_glc2000,        &
    &                             lai_mx_glc2000,        &
    &                             rs_min_glc2000,        &
    &                             urban_glc2000,         &
    &                             for_d_glc2000,         &
    &                             for_e_glc2000,         &
    &                             emissivity_glc2000

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
  &                          z0_tot,           &
  &                          root_lu,          &
  &                          plcov_mn_lu,      &
  &                          plcov_mx_lu,      &
  &                          lai_mn_lu,        &
  &                          lai_mx_lu,        &
  &                          rs_min_lu,        &
  &                          urban_lu,         &
  &                          for_d_lu,         &
  &                          for_e_lu,         &
  &                          emissivity_lu,    &
  &                          fr_ocean_lu,      &
  &                          lu_class_fraction,&
  &                          lu_class_npixel,  &
  &                          lu_tot_npixel,    &
  &                          lai12_lu,         &
  &                          plcov12_lu,       & 
  &                          z012_lu



  USE mo_landuse_output_nc, ONLY: write_netcdf_buffer_glc2000,     &
    &                             write_netcdf_cosmo_grid_glc2000, &
    &                             write_netcdf_icon_grid_glc2000

  USE mo_landuse_output_nc, ONLY: write_netcdf_buffer_glcc,     &   
    &                             write_netcdf_cosmo_grid_glcc, &
    &                             write_netcdf_icon_grid_glcc
  USE mo_landuse_output_nc, ONLY: write_netcdf_buffer_ecoclimap
  USE mo_landuse_output_nc, ONLY: write_netcdf_buffer_lu
           
  USE mo_globcover_lookup_tables, ONLY: nclass_globcover
  USE mo_globcover_lookup_tables, ONLY: init_globcover_lookup_tables

  USE mo_landuse_routines, ONLY: get_dimension_globcover_data, &
    &                            get_lonlat_globcover_data
  USE mo_ecoclimap_lookup_tables, ONLY: nclass_ecoclimap
  USE mo_ecoclimap_lookup_tables, ONLY: init_ecoclimap_lookup_tables
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
    &                         ntime_ecoclimap, &                     
    &                         allocate_raw_ecoclimap_fields

 USE mo_agg_globcover, ONLY : agg_globcover_data_to_target_grid
 USE mo_agg_ecoclimap, ONLY : agg_ecoclimap_data_to_target_grid

  IMPLICIT NONE
  
  CHARACTER(len=filename_max) :: filename
  CHARACTER(len=filename_max) :: netcdf_filename

  CHARACTER(len=filename_max) :: input_namelist_file
  CHARACTER(len=filename_max) :: input_namelist_cosmo_grid !< file with input namelist with COSMO grid definition

  CHARACTER(len=filename_max) :: namelist_grid_def

  CHARACTER (len=filename_max) :: namelist_topo_data_input !< file with input namelist with GLOBE data information
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

  INTEGER :: i,j,k !< counter
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

  INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid
  INTEGER  :: i_landuse_data !<integer switch to choose a land use raw data set
  INTEGER  :: ilookup_table_lu !< integer switch to choose a lookup table
  INTEGER  :: nclass_lu !< number of land use classes 

  
  ! Print the default information to stdout:
  CALL info_define ('extpar_landuse_to_buffer')      ! Pre-define the program name as binary name
  CALL info_print ()                     ! Print the information to stdout
  !--------------------------------------------------------------------------------------------------------

  namelist_grid_def = 'INPUT_grid_org'
  CALL init_target_grid(namelist_grid_def)

  igrid_type = tg%igrid_type
  tg_southern_bound=MINVAL(lat_geo) ! get southern boundary of target grid
  CALL allocate_lu_target_fields(tg)
  print *,'Grid defined, lu target fields allocated'

  !------------------------------------------------------------------------------------

  ! get information about landuse data

  ! get info on raw data file
  input_lu_namelist_file = 'INPUT_LU'

  !---------------------------------------------------------------------------
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

! >mes
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

      PRINT *, 'GLOBCOVER TILES, LON, LAT (MIN,MAX): ' 
      DO i = 1,ntiles_globcover
        WRITE(*,998)  i, lu_tiles_lon_min(i), lu_tiles_lon_max(i), &
                     lu_tiles_lat_min(i), lu_tiles_lat_max(i) 
998     FORMAT(I1,1X,4(F9.4,1X))      
      ENDDO

      PRINT *, 'MODEL DOMAIN, LON, LAT (MIN,MAX): ' 
      WRITE(*,999)  MINVAL(lon_geo), MAXVAL(lon_geo), &
                    MINVAL(lat_geo), MAXVAL(lat_geo)
999   FORMAT(4(F9.4,1X)) 

      DO i = 1,ntiles_globcover
        IF (lu_tiles_lon_min(i) < MINVAL(lon_geo).AND. &
            lu_tiles_lon_max(i) > MAXVAL(lon_geo).AND. &
            lu_tiles_lat_min(i) < MINVAL(lat_geo).AND. &
            lu_tiles_lat_max(i) > MAXVAL(lat_geo)) THEN
          PRINT *,'MODEL DOMAIN COVERED BY GLOBCOVER TILE ',i
        ENDIF
      ENDDO

    END SELECT

    ALLOCATE(lu_file(1:ntiles_lu), STAT= errorcode)
    IF(errorcode /= 0) CALL abort_extpar('Cant allocate lu_file')
! <mes
    print*, 'ntiles_lu: ', ntiles_lu
  PRINT *,'raw_data_glcc_filename: ',TRIM(raw_data_glcc_filename)

! >mes
  DO k = 1,ntiles_lu
    lu_file(k) = TRIM(raw_data_lu_path) // TRIM(raw_data_lu_filename(k))
    PRINT *,'lu_file: ', TRIM(lu_file(k))

  END DO
! <mes

  glcc_file(1) = TRIM(raw_data_glcc_path) // TRIM(raw_data_glcc_filename)
  PRINT *,'glcc file: ', TRIM(glcc_file(1))

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
        !HA debug
        PRINT *,'globcover_grid: ',globcover_grid
        ! If southern boundary of target grid is south of southern boundary of Globcover data
        ! (Globcover 2009 does not include Antarctica) then also process GLCC data)
        IF (tg_southern_bound < globcover_grid%end_lat_reg) THEN
          l_use_glcc=.TRUE.
          CALL allocate_glcc_target_fields(tg)
        ENDIF 

! >mes
      CALL get_globcover_tiles_grid(globcover_tiles_grid)
      print*,'globcover_tiles_grid(1): ', globcover_tiles_grid(1)
! <mes

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
!_br 17.03.14
    CASE DEFAULT
      PRINT*, 'ERROR: in extpar_landuse_to_buffer'
      PRINT*, 'ERROR: Invalid land use class ',i_landuse_data
      STOP 51
!_br 17.03.14 end
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

  undefined = 0.0_wp
  PRINT *,'aggregate land use data to target grid'
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
    PRINT *,'aggregation of land use data done'

  !--------------------------------------------------------------------------------
  ! output
   undefined = -999.0_wp
   undef_int = -999

   netcdf_filename = TRIM(lu_buffer_file)
   print *, 'Land-use buffer filename: ',TRIM(netcdf_filename)

   
 SELECT CASE (i_landuse_data)
   CASE(i_lu_glc2000, i_lu_globcover)
 
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
    &                                     emissivity_lu)

    IF (l_use_glcc) THEN !
      netcdf_filename = TRIM(glcc_buffer_file)
      print *, 'GLCC buffer filename: ',TRIM(netcdf_filename)
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


  PRINT *,'============= landuse_to_buffer done ==============='

END PROGRAM extpar_landuse_to_buffer

