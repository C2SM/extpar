!+  Fortran main program for consistency check of external parameters
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
!  introduce simple hight correction for T_CL field
!  Add DWD_min_lake_depth for minimal lake depth
! V1_3         2011/04/19 Hermann Asensio
!  introduce Globcover 2009 land use data set for external parameters
!  additional concistency check for vegetation on ice grid elements (glaciers)
! add support for GRIB1 and GRIB2
! V1_4         2011/04/21 Hermann Asensio
!  add a manual correction for the lake depth of Lake Constance
! V1_5         2011/08/11 Hermann Asensio
!  bug fix in the concistency check for FR_LAKE
!  remove erroneous simple height correction of T_CL field
! V1_7         2013/01/25 Guenther Zaengl
!   Parallel threads for ICON and COSMO using Open-MP,
!   Several bug fixes and optimizations for ICON search algorithm,
!   particularly for the special case of non-contiguous domains;
!   simplified namelist control for ICON
! V1_8         2013-03-12 Frank Brenner
!  introduced MODIS albedo dataset(s) as new external parameter(s)
! V1_11        2013/04/16 Juergen Helmert
!  Adaptions for using special points and external land-sea mask
! V1_12        2013-04-24 Frank Brenner
!  bug fix regarding old file paths
! V2_0         2013-06-04 Martina Messmer
!  introduce HWSD soil data set for external parameters (Juergen Helmert)
!  introduce lradtopo parameters for external parameters (Anne Roches)
!  add possibility to omit the SSO parameters (Martina Messmer)
!  introduce a finer temperature climatology data set (CRU) for external
!  parameters (CLM Community)
!  adjusted consistency check for soiltype on ice grid elements for
!  Globcover 2009
! V2_0_3       2015-01-12 Juergen Helmert
!  Bugfix correction covers CSCS SVN r5907-r6359
! V3_0         2015-05-18 Juergen Helmert, 2015-06-19 Daniel Luethi
!  Change tile_mode switch to integer, scale glacier soil with fr_land
!  Modification of Caspian Sea treatment (height below sea level, JH)
!  Adjust surface height (DL) and correct treatment also for ASTER topo
!  Adaptions for urban fields
! V4_0         2016/08/17 Daniel Lthi and authors from RHM
!  add support for subgrid-scale slope parameter
!  add support for MACv2 spectr. stratified monthly aerosol fields (iaot_type=4)
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran main program for consistency check of external parameters
!>
!> Purpose: read in the external paraemeters and check for consistency, especially for the land-sea mask.
!> \author Hermann Asensio
PROGRAM extpar_consistency_check

  USE info_extpar, ONLY: info_print
  USE mo_kind,     ONLY: wp, i4, i8
  USE mo_logging

  USE mo_target_grid_data, ONLY: lon_geo, lat_geo, tg

  USE mo_grid_structures, ONLY: igrid_icon, igrid_cosmo

  USE  mo_cosmo_grid, ONLY: COSMO_grid

  USE  mo_icon_grid_data, ONLY: icon_grid !< structure which contains the definition of the ICON grid
  USE  mo_icon_grid_data, ONLY: icon_grid_region

  USE mo_base_geometry,    ONLY:  geographical_coordinates

  USE mo_io_units,          ONLY: filename_max

  USE mo_read_extpar_namelists, ONLY: read_namelists_extpar_check_icon,    &
       &                              read_namelists_extpar_check_cosmo,   &
       &                              read_namelists_extpar_special_points

  USE mo_soil_routines, ONLY: read_namelists_extpar_soil

  USE mo_soil_data, ONLY: define_soiltype


  USE mo_soil_data, ONLY: undef_soiltype, default_soiltype, soiltype_ice, soiltype_water
  USE mo_soil_data, ONLY: soil_data, FAO_data, HWSD_data, HWSD_map

  USE mo_soil_tg_fields, ONLY:  fr_sand,fr_silt,fr_clay, &
       &                        fr_oc,fr_bd
  USE mo_soil_tg_fields, ONLY:  fr_sand_deep,fr_silt_deep, &
       &                        fr_clay_deep, fr_oc_deep,  &
       &                        fr_bd_deep
  USE mo_soil_tg_fields, ONLY:  fr_land_soil
  USE mo_soil_tg_fields, ONLY:  soiltype_fao, soiltype_hwsd, soiltype_deep,soiltype_hwsd_s
  USE mo_soil_tg_fields, ONLY:  allocate_soil_target_fields

  USE mo_soil_consistency, ONLY:  calculate_soiltype

  USE mo_soil_output_nc, ONLY: read_netcdf_soil_buffer

  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_glcc_lookup_tables, ONLY: nclass_glcc

  USE mo_glcc_tg_fields, ONLY:  fr_land_glcc, &
       &       glcc_class_fraction,    &
       &        glcc_class_npixel, &
       &        glcc_tot_npixel, &
       &        ice_glcc, &
       &        z0_glcc, &
       &        root_glcc, &
       &        plcov_mn_glcc, &
       &        plcov_mx_glcc, &
       &        lai_mn_glcc, &
       &        lai_mx_glcc, &
       &        rs_min_glcc, &
       &        urban_glcc,  &
       &        for_d_glcc,  &
       &        for_e_glcc, &
       &        emissivity_glcc, &
       &        allocate_glcc_target_fields

  USE mo_glc2000_lookup_tables, ONLY: get_name_glc2000_lookup_tables
  USE mo_glc2000_lookup_tables, ONLY: nclass_glc2000

  USE mo_globcover_data,   ONLY: max_tiles_lu  !_br 21.02.14 including definition of max_tiles_lu

  USE mo_globcover_lookup_tables, ONLY: nclass_globcover

  USE mo_globcover_lookup_tables, ONLY: get_name_globcover_lookup_tables
  USE mo_glcc_lookup_tables, ONLY: get_name_glcc_lookup_tables
  USE mo_ecoclimap_lookup_tables, ONLY: nclass_ecoclimap
  USE mo_ecoclimap_lookup_tables, ONLY: get_name_ecoclimap_lookup_tables


  USE mo_landuse_output_nc, ONLY: read_netcdf_buffer_glcc
  USE mo_landuse_output_nc, ONLY: read_netcdf_buffer_lu
  USE mo_landuse_output_nc, ONLY: read_netcdf_buffer_ecoclimap

  USE mo_isa_output_nc, ONLY: read_netcdf_buffer_isa

  USE mo_landuse_routines, ONLY: read_namelists_extpar_land_use

  USE mo_lu_tg_fields, ONLY :  i_lu_globcover, i_lu_glc2000, i_lu_glcc
  USE mo_lu_tg_fields, ONLY :  i_lu_ecoclimap

  USE mo_lu_tg_fields, ONLY: fr_land_lu, &
       &    fr_land_mask, &
       &        ice_lu, &
       &        z0_lu, &
       &        z0_tot, &
       &        root_lu, &
       &        plcov_mn_lu, &
       &        plcov_mx_lu, &
       &        lai_mn_lu, &
       &        lai_mx_lu, &
       &        rs_min_lu, &
       &        urban_lu,  &
       &        for_d_lu,  &
       &        for_e_lu, &
       &        emissivity_lu, &
       &        fr_ocean_lu, &
       &        lu_class_fraction,    &
       &        lu_class_npixel, &
       &        lu_tot_npixel,  &
       &        lai12_lu,     &
       &        z012_lu,      &
       &        plcov12_lu


  USE mo_lu_tg_fields, ONLY: allocate_lu_target_fields, allocate_add_lu_fields

  USE mo_albedo_tg_fields, ONLY: alb_dry, alb_sat, &
       &                        alb_field_mom, &
       &                        alnid_field_mom, &
       &                        aluvd_field_mom, &
       &                        allocate_alb_target_fields

  USE mo_albedo_data, ONLY: ntime_alb
  USE mo_albedo_data, ONLY: wso_min,wso_max,csalb,csalbw,zalso
  USE mo_albedo_data, ONLY: allocate_alb_interp_fields, &
       &                     alb_interp_data
  USE mo_albedo_data, ONLY: minimal_alb_dry, maximal_alb_dry, &
       &                      minimal_alb_sat, maximal_alb_sat, undef_alb_bs
  USE mo_albedo_data, ONLY: ialb_type

  USE mo_albedo_output_nc, ONLY: read_netcdf_buffer_alb

  USE mo_albedo_routines, ONLY: open_netcdf_ALB_data, &
       const_check_interpol_alb,&
       read_namelists_extpar_alb

  USE mo_albedo_data, ONLY: alb_raw_data_grid, &
       ntime_alb

  USE mo_isa_tg_fields, ONLY: isa_field, &
       &        isa_tot_npixel


  USE mo_isa_tg_fields, ONLY: allocate_isa_target_fields, &
       & allocate_add_isa_fields

  USE mo_isa_routines, ONLY: read_namelists_extpar_isa

  USE mo_ahf_tg_fields, ONLY: ahf_field, &
       &                                allocate_ahf_target_fields

  USE mo_ahf_data, ONLY: undef_ahf, minimal_ahf, iahf_type !_br 14.04.16
  USE mo_ahf_output_nc, ONLY: read_netcdf_buffer_ahf

  USE mo_ahf_routines, ONLY: read_namelists_extpar_ahf

  USE mo_ndvi_tg_fields, ONLY:          ndvi_max, &
       &                                ndvi_field_mom, &
       &                                ndvi_ratio_mom, &
       &                                allocate_ndvi_target_fields

  USE mo_era_tg_fields, ONLY: sst_field, &
       &                           wsnow_field, &
       &                           t2m_field, &
       &                           hsurf_field, &
       &                           allocate_era_target_fields


  USE mo_ndvi_data, ONLY: ntime_ndvi
  USE mo_ndvi_data, ONLY: undef_ndvi, minimal_ndvi

  USE mo_ndvi_output_nc, ONLY: read_netcdf_buffer_ndvi

  USE mo_era_output_nc, ONLY: read_netcdf_buffer_sst,&
       read_netcdf_buffer_t2m


  USE mo_topo_tg_fields, ONLY:  fr_land_topo,        &
       &                         hh_topo,            &
       &                         hh_topo_max,        &
       &                         hh_topo_min,        &       
       &                         stdh_topo,          &
       &                         theta_topo,         &
       &                         aniso_topo,         &
       &                         slope_topo,         &
       &                         z0_topo,            &
       &                         slope_asp_topo,     &
       &                         slope_ang_topo,     &
       &                         horizon_topo,       &
       &                         skyview_topo,       &
       &                         allocate_topo_target_fields

  USE mo_topo_tg_fields, ONLY: vertex_param, &
       &        allocate_additional_hh_param

  USE mo_topo_output_nc, ONLY: read_netcdf_buffer_topo

  USE mo_topo_routines, ONLY: read_namelists_extpar_orography, &
       &                       read_namelists_extpar_scale_sep

  USE mo_topo_data, ONLY: lradtopo, nhori, max_tiles, itopo_type

  USE mo_sgsl_tg_fields, ONLY: sgsl, allocate_sgsl_target_fields

  USE mo_sgsl_output_nc, ONLY: read_netcdf_buffer_sgsl

  USE mo_sgsl_routines, ONLY: read_namelists_extpar_sg_slope

  USE mo_sgsl_data, ONLY: idem_type

  USE mo_aot_target_fields, ONLY: allocate_aot_target_fields,&
       &                              aot_tg,&
       &                              MAC_aot_tg,&
       &                              MAC_ssa_tg,&
       &                              MAC_asy_tg


  USE mo_aot_output_nc, ONLY: read_netcdf_buffer_aot, read_netcdf_buffer_aot_MAC

  USE mo_cru_target_fields, ONLY: allocate_cru_target_fields,   &
       &                              crutemp, crutemp2, cruelev,&
       &                              i_t_cru_fine, i_t_cru_coarse

  USE mo_cru_data,  ONLY: read_namelists_extpar_t_clim

  USE mo_cru_output_nc, ONLY: read_netcdf_buffer_cru

  USE mo_aot_data, ONLY: ntype_aot, ntime_aot, iaot_type, n_spectr, nspb_aot

  USE mo_aot_data, ONLY: read_namelists_extpar_aerosol

  USE mo_flake_tg_fields, ONLY: fr_lake, &
       &       lake_depth,    &
       &       flake_tot_npixel, &
       &       allocate_flake_target_fields
  USE mo_flake_data, ONLY: flake_depth_undef !< default value for undefined lake depth
  USE mo_flake_data, ONLY: flake_depth_default !< default value for default lake depth, 10 [m]
  USE mo_flake_data, ONLY: DWD_max_lake_depth !< Maximum lake depth in [m] for FLAKE (50 m)
  USE mo_flake_data, ONLY: DWD_min_lake_depth !< Minimal lake depth in [m] for FLAKE (1 m)

  USE mo_flake_output_nc, ONLY: read_netcdf_buffer_flake

  USE mo_lsm_output_nc, ONLY: read_netcdf_buffer_lsm

  USE mo_extpar_output_nc, ONLY: write_netcdf_icon_grid_extpar, &
       &                        write_netcdf_cosmo_grid_extpar


  USE mo_lradtopo,   ONLY: read_namelists_extpar_lradtopo

  USE mo_search_target_grid, ONLY: find_nearest_target_grid_element

  USE mo_oro_filter, ONLY: read_namelists_extpar_orosmooth

  USE mo_isa_data,   ONLY: max_tiles_isa
  USE mo_isa_data, ONLY: undef_isa, minimal_isa, isa_type !_br 14.04.16

  IMPLICIT NONE

  CHARACTER (len=filename_max) :: namelist_grid_def
  CHARACTER (len=filename_max) :: netcdf_output_filename

  CHARACTER (len=filename_max) :: grib_output_filename
  CHARACTER (len=filename_max) :: grib_sample

  !-----------------------------------------------------------------------------------------------------------------------
  CHARACTER (len=filename_max) :: namelist_file !< filename with namelists for for EXTPAR settings

  ! soil
  CHARACTER (len=filename_max) :: soil_buffer_file  !< name for soil buffer file
  CHARACTER (len=filename_max) :: soil_output_file  !< name for soil output file
  CHARACTER (len=filename_max) :: soil_buffer_file_consistent !< name for soil buffer file after consistency check
  CHARACTER (len=filename_max) :: soil_output_file_consistent !< name for soil output file after consistency check

  CHARACTER (len=filename_max) :: raw_data_soil_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_soil_filename !< filename soil raw data
  CHARACTER (len=filename_max) :: raw_data_deep_soil_filename !< filename deep soil raw data

  ! orography
  CHARACTER (len=filename_max) :: orography_buffer_file !< name for orography buffer file
  CHARACTER (len=filename_max) :: orography_output_file !< name for orography output file

  CHARACTER (len=filename_max) :: raw_data_orography_path        !< path to raw data

  ! subgrid-scale slope
  CHARACTER (LEN=filename_max) :: sgsl_files(1:max_tiles)  !< filenames globe raw data
  CHARACTER (len=filename_max) :: sgsl_buffer_file !< name for orography buffer file
  CHARACTER (len=filename_max) :: raw_data_sgsl_path        !< path to raw data

  ! land use

  CHARACTER (len=filename_max) :: raw_data_lu_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_lu_filename(1:max_tiles_lu) !< filename glc2000 raw data !_br 21.02.14
  CHARACTER(len=filename_max) :: name_lookup_table_lu !< name for look up table
  CHARACTER(len=filename_max) :: lu_dataset !< name of landuse data set


  CHARACTER (len=filename_max) :: lu_buffer_file !< name for glc2000 buffer file
  CHARACTER (len=filename_max) :: lu_output_file !< name for glc2000 output file

  CHARACTER (len=filename_max) :: glcc_buffer_file    !< name for glcc buffer file

  !-----------------------------------------------------------------------------------------------------------------------
  ! albedo
  CHARACTER (len=filename_max) :: raw_data_alb_path   !< path to albedo raw input data
  CHARACTER (len=filename_max) :: raw_data_alb_filename    !< raw data filename
  CHARACTER (len=filename_max) :: raw_data_alnid_filename  !< raw data filename, NI
  CHARACTER (len=filename_max) :: raw_data_aluvd_filename  !< raw data filename, UV

  CHARACTER (len=filename_max) :: alb_buffer_file    !< name for albedo buffer file
  CHARACTER (len=filename_max) :: alb_output_file    !< name for albedo output file
  !-----------------------------------------------------------------------------------------------------------------------
  CHARACTER (len=filename_max) :: land_sea_mask_file !< name for external land-sea-mask

  ! ISA
  CHARACTER (len=filename_max) :: raw_data_isa_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_isa_filename(1:max_tiles_isa) !< filename glc2000 raw data
  INTEGER                      :: ntiles_isa
  CHARACTER (len=filename_max) :: isa_buffer_file !< name for NDVI buffer file

  !AHF
  CHARACTER (len=filename_max) :: raw_data_ahf_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_ahf_filename !< filename NDVI raw data

  CHARACTER (len=filename_max) :: ahf_buffer_file !< name for NDVI buffer file
  CHARACTER (len=filename_max) :: ahf_output_file !< name for NDVI output file

  ! NDVI
  CHARACTER (len=filename_max) :: ndvi_buffer_file !< name for NDVI buffer file

  CHARACTER (len=filename_max) :: sst_icon_file !< name for SST icon file
  CHARACTER (len=filename_max) :: t2m_icon_file !< name for SST icon file

  ! temperature climatology
  CHARACTER (len=filename_max) :: namelist_file_t_clim !< filename with namelists for for EXTPAR settings
  CHARACTER (len=filename_max) :: raw_data_t_clim_path     !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_t_clim_filename !< filename temperature climatology raw data


  CHARACTER (len=filename_max) :: t_clim_buffer_file !< name for temperature climatology buffer
  CHARACTER (len=filename_max) :: t_clim_dummy_file  !< to avoid overwriting when reading the INPUT_CHECK namelist

  CHARACTER (len=filename_max) :: t_clim_output_file !< name for temperature climatology output file

  ! aerosol optical thickness

  CHARACTER (len=filename_max) :: raw_data_aot_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_aot_filename !< filename temperature climatology raw data

  CHARACTER (len=filename_max) :: flake_buffer_file  !< name for flake buffer file

  CHARACTER (len=filename_max) :: aot_buffer_file !< name for aerosol buffer file
  CHARACTER (len=filename_max) :: aot_output_file !< name for aerosol output file

  CHARACTER (len=filename_max) :: topo_files(1:max_tiles) !< filenames globe raw data
  INTEGER (KIND=i4)  :: ntiles_column
  INTEGER (KIND=i4)  :: ntiles_row
  INTEGER (KIND=i8)  :: it_cl_type
  INTEGER (KIND=i4)  :: isoil_data
  LOGICAL            :: lsso_param,lsubtract_mean_slope
  LOGICAL            :: ldeep_soil

  REAL(KIND=wp) :: undefined !< value to indicate undefined grid elements
  INTEGER (KIND=i4) :: undef_int   !< value for undefined integer

  REAL(KIND=wp), PARAMETER :: undefined_lu = 0.0_wp !< value to indicate undefined land use grid elements
  REAL(KIND=wp) :: thr_cr !< control threshold
  REAL(KIND=wp) :: fill_value_real !< value to indicate undefined grid elements
  INTEGER (KIND=i4) :: fill_value_int   !< value for undefined integer

  INTEGER (KIND=i4), PARAMETER :: mpy=12     !< month per year

  INTEGER(i8) :: i !< counter
  INTEGER(i8) :: j !< counter
  INTEGER(i8) :: k !< counter
  INTEGER :: t !< counter

  INTEGER(i8) :: jj !< counter
  INTEGER(i8) :: ii !< counter

  INTEGER :: n
  INTEGER :: nloops

  INTEGER :: nnb !< number of neighbor grid elements with common edge
  INTEGER :: nv  !< number of vertices
  INTEGER :: n_index !< help variable
  INTEGER(i8) :: ne_ie(9) !< index for grid element neighbor
  INTEGER(i8) :: ne_je(9) !< index for grid element neighbor
  INTEGER(i8) :: ne_ke(9) !< index for grid element neighbor

  INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid
  INTEGER  :: i_landuse_data !<integer switch to choose a land use raw data set
  INTEGER  :: i_lsm_data !<integer switch to choose a land sea mask data set
  INTEGER  :: ilookup_table_lu !< integer switch to choose a lookup table
  INTEGER  :: nclass_lu !< number of land use classes

  ! variables for the T_CLIM
  LOGICAL ::  last=.FALSE. ! in TCL leave loop
  LOGICAL ::  foundtcl=.FALSE. ! in TCL

  LOGICAL :: l_use_isa=.FALSE. !< flag if additional urban data are present
  LOGICAL :: l_use_ahf=.FALSE. !< flag if additional urban data are present
  LOGICAL :: l_use_sgsl=.FALSE. !< flag if additional urban data are present
  LOGICAL :: l_use_glcc=.FALSE. !< flag if additional glcc data are present
  REAL :: lu_data_southern_boundary

  LOGICAL :: lwrite_netcdf  !< flag to enable netcdf output for COSMO
  LOGICAL :: lwrite_grib    !< flag to enable GRIB output for COSMO
  LOGICAL :: lflake_correction !< flag to correct fr_lake and depth_lake near coastlines

  !HA tests
  REAL :: timestart
  REAL :: timeend
  REAL :: timediff

  INTEGER           :: db_ice_counter,db_water_counter


  REAL (KIND=wp) :: step

  !Special Points
  INTEGER          :: number_special_points
  REAL (KIND=wp)   ::                           lon_geo_sp=-999.,           &
       lat_geo_sp=-999.,           &
       soiltype_sp=-999.,          &
       z0_sp=-999.,                &
       rootdp_sp=-999.,            &
       plcovmn_sp=-999.,           &
       plcovmx_sp=-999.,           &
       laimn_sp=-999.,             &
       laimx_sp=-999.,             &
       for_d_sp=-999.,             &
       for_e_sp=-999.,             &
       fr_land_sp=-999.

  INTEGER (KIND=i8) :: start_cell_id !< ID of starting cell for ICON search
  INTEGER (KIND=i8) :: isp,i_sp,j_sp,k_sp
  INTEGER (KIND=i4) :: ncid_alb
  INTEGER (KIND=i8) :: nlon_reg !< number of columns
  INTEGER (KIND=i8) :: nlat_reg !< number of rows
  CHARACTER (LEN=filename_max) :: path_alb_file
  CHARACTER (LEN=filename_max) :: alb_source, alnid_source, aluvd_source
  CHARACTER (LEN=filename_max) :: namelist_alb_data_input
  CHARACTER (LEN=filename_max), DIMENSION(1:23) :: glc_class
  INTEGER                      :: tile_mode
  LOGICAL                      :: tile_mask, ltcl_merge

  ! for albedo consistency check
  REAL (KIND=wp) :: albvis_min, albnir_min, albuv_min
  REAL (KIND=wp) :: hh_dead_sea, hh_cr_casp

  ! T_CL consistency check
  INTEGER (KIND=i8) :: iml, imu, ipl, ipu, jml, jmu, jpl, jpu, l
  INTEGER (KIND=i4) :: ntclct
  REAL    (KIND=wp) :: tclsum, elesum

  ! Namelist values for topography scale separation
  CHARACTER (LEN=filename_max) :: raw_data_scale_sep_orography_path !< path to raw data
  CHARACTER (LEN=filename_max) :: scale_sep_files(1:max_tiles)  !< filenames globe raw data
  LOGICAL           :: lscale_separation

  ! Namelist values for orography smoothing
  LOGICAL           :: &
       lfilter_oro,     &
       lxso_first

  INTEGER(KIND=i4)  :: &
       ilow_pass_oro,   &
       numfilt_oro,     &
       ifill_valley,    &
       ilow_pass_xso,   &
       numfilt_xso

  REAL(KIND=wp)     :: &
       eps_filter,      &
       rfill_valley,    &
       rxso_mask

  ! define string used for global attributes:
  CHARACTER (LEN=255) :: y_orofilter

  !--------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------

  CALL initialize_logging("extpar_consistency.log", stdout_level=debug)
  CALL info_print ()

  !--------------------------------------------------------------------------------------------------------
  ! Get lradtopo and nhori value from namelist

  namelist_file = 'INPUT_RADTOPO'
  CALL read_namelists_extpar_lradtopo(namelist_file,lradtopo,nhori)

  ! Get lsso_param from namelist

  namelist_file = 'INPUT_ORO'
  CALL read_namelists_extpar_orography(namelist_file,          &
       raw_data_orography_path,&
       topo_files,             &
       ntiles_column,          &
       ntiles_row,             &
       itopo_type,             &
       lsso_param,             &
       lsubtract_mean_slope,  &
       orography_buffer_file,  &
       orography_output_file)

  namelist_file = 'INPUT_SOIL'
  CALL read_namelists_extpar_soil(namelist_file,                     &
       isoil_data,                 &
       ldeep_soil,                 &
       raw_data_soil_path,         &
       raw_data_soil_filename,     &
       raw_data_deep_soil_filename,&
       soil_buffer_file,           &
       soil_output_file,           &
       soil_buffer_file_consistent,&
       soil_output_file_consistent)

  IF (ldeep_soil .AND. isoil_data /= HWSD_data) THEN  !_br 21.02.14 replaced eq by eqv
     ldeep_soil = .FALSE.
     CALL logging%error('one can use deep soil only, if HWSD data is used - ldeep_soil is set to FALSE', __FILE__, __LINE__)
  ENDIF

  WRITE(message_text,'(a,i0)') 'isoil_data: ', isoil_data
  CALL logging%info(message_text, __FILE__, __LINE__)
  WRITE(message_text,'(a,l1)') 'ldeep_soil: ', ldeep_soil
  CALL logging%info(message_text, __FILE__, __LINE__)



  !--------------------------------------------------------------
  ! get namelist for albedo fields
  !--------------------------------------------------------------
  namelist_alb_data_input = 'INPUT_ALB'
  CALL  read_namelists_extpar_alb(namelist_alb_data_input, &
       &                                  raw_data_alb_path, &
       &                                  raw_data_alb_filename, &
       &                                  raw_data_alnid_filename, &
       &                                  raw_data_aluvd_filename, &
       &                                  ialb_type,        &
       &                                  alb_buffer_file, &
       &                                  alb_output_file, &
       &                                  alb_source, &
       &                                  alnid_source, &
       &                                  aluvd_source)

  !  print*, 'ialb_type: ', ialb_type

  !--------------------------------------------------------------
  ! get namelist for aerosol fields
  !--------------------------------------------------------------
  namelist_file = 'INPUT_AOT'
  CALL read_namelists_extpar_aerosol(namelist_file, &
       &                                  iaot_type,    &
       &                                  raw_data_aot_path, &
       &                                  raw_data_aot_filename, &
       &                                  aot_buffer_file, &
       &                                  aot_output_file)

  !--------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------
  ! get namelist for topography scale separation
  !--------------------------------------------------------------
  namelist_file = 'INPUT_SCALE_SEP'
  CALL read_namelists_extpar_scale_sep(namelist_file,                     &
       &                                  raw_data_scale_sep_orography_path, &
       &                                  scale_sep_files,                   &
       &                                  lscale_separation)

  !--------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------
  ! get namelist for topography smoothing
  !--------------------------------------------------------------
  namelist_file = 'INPUT_OROSMOOTH'
  CALL read_namelists_extpar_orosmooth(namelist_file,            &
       lfilter_oro,          &
       ilow_pass_oro,        &
       numfilt_oro,          &
       eps_filter,           &
       ifill_valley,         &
       rfill_valley,         &
       ilow_pass_xso,        &
       numfilt_xso,          &
       lxso_first,           &
       rxso_mask)

  !--------------------------------------------------------------------------------------------------------
  ! get information on target grid, allocate target fields with coordinates and determin the coordinates
  ! for th target grid

  namelist_grid_def = 'INPUT_grid_org'
  CALL  init_target_grid(namelist_grid_def)
  PRINT *,' target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
  igrid_type = tg%igrid_type


  SELECT CASE(igrid_type)
  CASE(igrid_icon) ! ICON GRID
     PRINT *,'icon_grid%ncell: ',icon_grid%ncell
     PRINT *,'icon_grid%nvertex: ',icon_grid%nvertex

     CALL  allocate_additional_hh_param(icon_grid%nvertex)

  END SELECT


  !--------------------------------------------------------------------------------------------------------
  ! get info on raw data file
  namelist_file = 'INPUT_LU'

  !---------------------------------------------------------------------------
  CALL read_namelists_extpar_land_use(namelist_file, &
       &                                 i_landuse_data, &
       &                                 raw_data_lu_path, &
       &                                 raw_data_lu_filename, &
       &                                 ilookup_table_lu, &
       &                                 lu_buffer_file, &
       &                                 lu_output_file)

  lu_data_southern_boundary = -91.0  
  SELECT CASE (i_landuse_data)
  CASE (i_lu_ecoclimap)
     lu_dataset = 'ECOCLIMAP'
     CALL get_name_ecoclimap_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
     nclass_lu = nclass_ecoclimap
     lu_data_southern_boundary = -90.9
  CASE (i_lu_globcover)
     lu_dataset = 'GLOBCOVER2009'
     CALL get_name_globcover_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
     nclass_lu = nclass_globcover
     lu_data_southern_boundary = -64.99
  CASE (i_lu_glc2000)
     lu_dataset = 'GLC2000'
     CALL get_name_glc2000_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
     nclass_lu = nclass_glc2000
     lu_data_southern_boundary = -56.0
  CASE(i_lu_glcc)
     lu_dataset = 'GLCC'
     CALL get_name_glcc_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
     nclass_lu = nclass_glcc
     lu_data_southern_boundary = -90.0
  END SELECT

  CALL logging%info('Land use datatset    : '//TRIM(lu_dataset), __FILE__, __LINE__)
  CALL logging%info('Land use lookup table: '//TRIM(name_lookup_table_lu), __FILE__, __LINE__)

  !-----------------------------------------------------------------------------------------------
  ! get info on urban data file
  !-----------------------------------------------------------------------------------------------
  namelist_file = 'INPUT_ISA'
  INQUIRE(file=TRIM(namelist_file),exist=l_use_isa)
  IF (l_use_isa) THEN
     PRINT *,'URBAN DATA ISA active: Reading INPUT-File...'
     CALL read_namelists_extpar_isa(namelist_file, &
          &                                 isa_type,    & !_br 14.04.16
          &                                 raw_data_isa_path,       &
          &                                 raw_data_isa_filename,   &
          &                                 ntiles_isa,       &
          &                                 isa_buffer_file)
  END IF

  namelist_file = 'INPUT_AHF'
  INQUIRE(file=TRIM(namelist_file),exist=l_use_ahf)
  IF (l_use_ahf) THEN
     PRINT *,'URBAN DATA AHF active: Reading INPUT-File...'
     CALL  read_namelists_extpar_ahf(namelist_file, &
          &                                  iahf_type,    & !_br 14.04.16
          &                                  raw_data_ahf_path, &
          &                                  raw_data_ahf_filename, &
          &                                  ahf_buffer_file, &
          &                                  ahf_output_file)
  END IF


  !-----------------------------------------------------------------------------------------------
  ! get info on subgrid-scale slope file
  !-----------------------------------------------------------------------------------------------
  namelist_file = 'INPUT_SGSL'
  INQUIRE(file=TRIM(namelist_file),exist=l_use_sgsl)
  IF (l_use_sgsl) THEN
     CALL  read_namelists_extpar_sg_slope(namelist_file,    &
          &                                  raw_data_sgsl_path, &
          &                                  sgsl_files,         &
          &                                  ntiles_column,      &
          &                                  ntiles_row,         &
          &                                  idem_type,          &
          &                                  sgsl_buffer_file)
  END IF

  !-----------------------------------------------------------------------------------------------
  ! get filenames from namelist
  !-----------------------------------------------------------------------------------------------
  namelist_file = 'INPUT_CHECK'

  SELECT CASE(igrid_type)
  CASE(igrid_icon)
    CALL logging%info('Read INPUT_CHECK for ICON', __FILE__, __LINE__)
    CALL read_namelists_extpar_check_icon(namelist_file, &
         grib_output_filename, &
         grib_sample, &
         netcdf_output_filename, &
         orography_buffer_file, &
         soil_buffer_file, &
         lu_buffer_file, &
         glcc_buffer_file, &
         flake_buffer_file, &
         ndvi_buffer_file, &
         sst_icon_file, &
         t2m_icon_file, &
         t_clim_dummy_file, &
         aot_buffer_file, &
         alb_buffer_file, &
         i_lsm_data, &
         land_sea_mask_file,&
         lwrite_netcdf,         &
         lwrite_grib,           &
         number_special_points, &
         tile_mode,             &
         ltcl_merge )

  CASE(igrid_cosmo)
    CALL logging%info('Read INPUT_CHECK for COSMO', __FILE__, __LINE__)
    CALL read_namelists_extpar_check_cosmo(namelist_file, &
         grib_output_filename, &
         grib_sample, &
         netcdf_output_filename, &
         orography_buffer_file, &
         soil_buffer_file, &
         lu_buffer_file, &
         glcc_buffer_file, &
         flake_buffer_file, &
         ndvi_buffer_file, &
         t_clim_dummy_file, &
         aot_buffer_file, &
         alb_buffer_file, &
         i_lsm_data, &
         land_sea_mask_file,&
         lwrite_netcdf,         &
         lwrite_grib,           &
         number_special_points, &
         tile_mode,             &
         lflake_correction,     &
         ltcl_merge)

  END SELECT

  IF(ltcl_merge) THEN
 
    namelist_file_t_clim = 'INPUT_TCLIM_FINAL'
    CALL read_namelists_extpar_t_clim(namelist_file_t_clim,     &
                                      it_cl_type,               &
                                      raw_data_t_clim_path,     &
                                      raw_data_t_clim_filename, &
                                      t_clim_buffer_file , &
                                      t_clim_output_file  )

    WRITE(message_text,'(a,a)')   'T_CL Merging from CRU Coarse : ', TRIM(t_clim_output_file)
    CALL logging%info(message_text, __FILE__, __LINE__)
    WRITE(message_text,'(a,a)')  'with T_CL from CRU Fine:       ', TRIM(t_clim_buffer_file)
    CALL logging%info(message_text, __FILE__, __LINE__)
  ELSE
    namelist_file_t_clim = 'INPUT_TCLIM'
    CALL read_namelists_extpar_t_clim(namelist_file_t_clim,     &
         it_cl_type,               &
         raw_data_t_clim_path,     &
         raw_data_t_clim_filename, &
         t_clim_buffer_file      , &
         t_clim_output_file        )
  END IF

  WRITE(message_text,'(a,i0)') 'TCLIM (it_cl_type): ', it_cl_type
  CALL logging%info(message_text, __FILE__, __LINE__)

  IF (tile_mode == 1) THEN
     tile_mask=.TRUE.
     PRINT*, 'Tile mode for EXTPAR is set to tile_mode= ',tile_mode,'tile_mask= ',tile_mask
  END IF

  ! test for glcc data
  INQUIRE(file=TRIM(glcc_buffer_file),exist=l_use_glcc)
  IF (l_use_glcc) THEN
    CALL allocate_glcc_target_fields(tg)
    CALL logging%info('GLCC fields allocated', __FILE__, __LINE__)
  ENDIF
  ! allocate Land use target fields
  CALL allocate_lu_target_fields(tg)
  CALL allocate_add_lu_fields(tg,nclass_lu)
  CALL logging%info('Land Use fields allocated', __FILE__, __LINE__)

  CALL allocate_soil_target_fields(tg,ldeep_soil)
  CALL logging%info('soil fields allocated', __FILE__, __LINE__)

  IF (l_use_ahf) THEN
    CALL allocate_ahf_target_fields(tg)
    CALL logging%info('AHF fields allocated', __FILE__, __LINE__)
  END IF
  IF (l_use_isa) THEN
    CALL allocate_isa_target_fields(tg)
    CALL allocate_add_isa_fields(tg)
    CALL logging%info('ISA fields allocated', __FILE__, __LINE__)
  END IF

  CALL allocate_ndvi_target_fields(tg,ntime_ndvi)
  PRINT *,'ntime_ndvi ', ntime_ndvi
  CALL logging%info('NDVI fields allocated', __FILE__, __LINE__)

  IF (l_use_sgsl) THEN
    CALL allocate_sgsl_target_fields(tg)
    CALL logging%info('SGSL fields allocated', __FILE__, __LINE__)
  END IF

  CALL allocate_era_target_fields(tg,ntime_ndvi) ! sst clim contains also 12 monthly values as ndvi
  PRINT *,'ntime_sst ', ntime_ndvi
  CALL logging%info('ERA-I SST/W_SNOW fields allocated', __FILE__, __LINE__)
  CALL logging%info('ERA-I T2M/HSURF fields allocated', __FILE__, __LINE__)

  IF (lscale_separation .AND. (itopo_type == 2)) THEN   !_br 21.02.14 replaced eq by eqv
    lscale_separation = .FALSE.
    CALL logging%warning('Scale separation can only be used with GLOBE topography', __FILE__, __LINE__)
  ENDIF

  CALL allocate_topo_target_fields(tg,nhori)
  CALL logging%info('TOPO fields allocated', __FILE__, __LINE__)

  CALL allocate_aot_target_fields(tg, iaot_type, ntime_aot, ntype_aot, nspb_aot)
  CALL logging%info('AOT fields allocated', __FILE__, __LINE__)

  CALL allocate_cru_target_fields(tg)
  CALL logging%info('CRU temperature field allocated', __FILE__, __LINE__)

  CALL allocate_flake_target_fields(tg)
  CALL logging%info('FLAKE parameter fields allocated', __FILE__, __LINE__)

  CALL allocate_alb_target_fields(tg,ntime_alb,ialb_type)
  CALL logging%info('ALBEDO fields allocated', __FILE__, __LINE__)

  !-----------------------------------------------------------------------------------------------
  ! Start Input

  PRINT *,'Read in Land Use data'
  PRINT *,'read ', TRIM(lu_buffer_file)

  SELECT CASE (i_landuse_data)
  CASE (i_lu_ecoclimap)
     lu_dataset = 'ECOCLIMAP'
     CALL read_netcdf_buffer_ecoclimap(lu_buffer_file,  &
          &                                     tg,         &
          &                                     nclass_lu, &
          &                                     undefined, &
          &                                     undef_int,   &
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


  CASE(i_lu_globcover, i_lu_glc2000 )
     PRINT *,'read ', TRIM(lu_buffer_file)
     CALL read_netcdf_buffer_lu(lu_buffer_file,  &
          &                                     tg,         &
          &                                     nclass_lu, &
          &                                     undefined, &
          &                                     undef_int,   &
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


     PRINT *,'MAX ICE_LU GLOBCOVER: ', MAXVAL(ICE_LU)
     PRINT *,'MAX lu_class_fraction_22 (ICE) GLOBCOVER: ', MAXVAL(lu_class_fraction(:,:,:,22))

     !    ice_lu = lu_class_fraction(:,:,:,22)
     !    IF (igrid_type == igrid_icon) THEN
     !      ice_lu = 0.0_wp
     !    ENDIF


     IF (l_use_glcc) THEN
        PRINT *,'read ', TRIM(glcc_buffer_file)
        CALL read_netcdf_buffer_glcc(glcc_buffer_file,  &
             &                                     tg,         &
             &                                     undefined, &
             &                                     undef_int,   &
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
        PRINT *,'MAX ICE_LU GLCC:', MAXVAL(ICE_GLCC)
     ENDIF
  END SELECT ! GlobCover needs also GLCC!


  PRINT *,'Read in soil data'
  IF(ldeep_soil) THEN
     CALL read_netcdf_soil_buffer(soil_buffer_file,    &
          &                                     tg,          &
          &                                     isoil_data,  &
          &                                     undefined,   &
          &                                     undef_int,   &
          &                                     fr_land_soil,&
          &                                     soiltype_fao,&
          &                                     soiltype_hwsd,&
          &                                     soiltype_deep,&
          &                                     soiltype_hwsd_s)
     PRINT *,'Selected HWSD + Deep Soil'
     print*, 'MIN/MAX soiltype_FAO : ', MINVAL(soiltype_fao), MAXVAL(soiltype_fao)
     print*, 'MIN/MAX soiltype_HWSD : ', MINVAL(soiltype_hwsd), MAXVAL(soiltype_hwsd)
     print*, 'MIN/MAX soiltype_FAO_deep : ', MINVAL(soiltype_deep), MAXVAL(soiltype_deep)
     print*, 'MIN/MAX soiltype_HWSD_deep : ', MINVAL(soiltype_hwsd_s), MAXVAL(soiltype_hwsd_s)
  ELSE
     SELECT CASE(isoil_data)
     CASE(FAO_data, HWSD_map)
        PRINT *,'Selected FAO, read ', TRIM(soil_buffer_file)
        CALL read_netcdf_soil_buffer(soil_buffer_file,    &
             &                                     tg,          &
             &                                     isoil_data,  &
             &                                     undefined,   &
             &                                     undef_int,   &
             &                                     fr_land_soil,&
             &                                     soiltype_fao,&
             &                                     soiltype_hwsd)
     CASE(HWSD_data)
        print*, 'MIN/MAX soiltype_FAO : ', MINVAL(soiltype_fao), MAXVAL(soiltype_fao)
        print*, 'MIN/MAX soiltype_HWSD : ', MINVAL(soiltype_hwsd), MAXVAL(soiltype_hwsd)
        PRINT *,'Selected HWSD'
        print*, 'Read ',soil_buffer_file
        CALL read_netcdf_soil_buffer(soil_buffer_file,    &
             &                                     tg,          &
             &                                     isoil_data,  &
             &                                     undefined,   &
             &                                     undef_int,   &
             &                                     fr_land_soil,&
             &                                     soiltype_fao,&
             &                                     soiltype_hwsd )

        print*, 'MIN/MAX soiltype_FAO : ', MINVAL(soiltype_fao), MAXVAL(soiltype_fao)
        print*, 'MIN/MAX soiltype_HWSD : ', MINVAL(soiltype_hwsd), MAXVAL(soiltype_hwsd)
     END SELECT

  ENDIF

  IF (l_use_isa) THEN
     PRINT *,'Read in ISA data'
     CALL read_netcdf_buffer_isa(isa_buffer_file,  &
          &                                     tg,         &
          &                                     undefined, &
          &                                     undef_int,   &
          &                                     isa_field, &
          &                                     isa_tot_npixel )
  END IF
  IF (l_use_ahf) THEN
     PRINT *,'Read in AHF data'
     CALL read_netcdf_buffer_ahf(ahf_buffer_file,  &
          &                                     tg,         &
          &                                     undefined, &
          &                                     undef_int, &
          &                                     ahf_field )
  END IF

  IF(igrid_type == igrid_icon) THEN
     PRINT *,'Read in SST data from ', TRIM(sst_icon_file)
     CALL read_netcdf_buffer_sst(sst_icon_file,  &
          &                                     tg,         &
          &                                     ntime_ndvi, &
          &                                     undefined, &
          &                                     undef_int,   &
          &                                     sst_field,&
          &                                     wsnow_field)

     PRINT *,'Read in T2M data from ', TRIM(t2m_icon_file)
     CALL read_netcdf_buffer_t2m(t2m_icon_file,  &
          &                                     tg,         &
          &                                     ntime_ndvi, &
          &                                     undefined, &
          &                                     undef_int,   &
          &                                     t2m_field,&
          &                                     hsurf_field)

  END IF

  PRINT *,'Read in albedo data'
  IF (ialb_type == 2) THEN
     CALL read_netcdf_buffer_alb(alb_buffer_file,  &
          &                           tg, &
          &                           ntime_alb, &
          &                           undefined, &
          &                           undef_int,   &
          &                           alb_dry=alb_dry, &
          &                           alb_sat=alb_sat)
  ELSE IF (ialb_type == 1) THEN
     CALL read_netcdf_buffer_alb(alb_buffer_file,  &
          &                           tg, &
          &                           ntime_alb, &
          &                           undefined, &
          &                           undef_int,   &
          &                           alb_field_mom, &
          &                           alnid_field_mom, &
          &                           aluvd_field_mom)
  ELSE IF (ialb_type == 3) THEN
     CALL read_netcdf_buffer_alb(alb_buffer_file,  &
          &                           tg, &
          &                           ntime_alb, &
          &                           undefined, &
          &                           undef_int,   &
          &                           alb_field_mom)
  ENDIF


  PRINT *,'Read in NDVI data'
  CALL read_netcdf_buffer_ndvi(ndvi_buffer_file,  &
       &                                     tg,         &
       &                                     ntime_ndvi, &
       &                                     undefined, &
       &                                     undef_int,   &
       &                                     ndvi_max,  &
       &                                     ndvi_field_mom,&
       &                                     ndvi_ratio_mom)


  PRINT *,'MAX/MIN of ERA-I SST and W_SNOW ',MAXVAL(sst_field),MINVAL(sst_field),MAXVAL(wsnow_field),MINVAL(wsnow_field)
  PRINT *,'MAX/MIN of ERA-I T2M and HSURF ',MAXVAL(t2m_field),MINVAL(t2m_field),MAXVAL(hsurf_field),MINVAL(hsurf_field)

  PRINT *,'Read in orography data'

  SELECT CASE(igrid_type)
  CASE(igrid_icon) ! ICON GRID

    IF (lsso_param) THEN

      CALL read_netcdf_buffer_topo(orography_buffer_file,                 &
           &                                     tg,                      &
           &                                     undefined,               &
           &                                     undef_int,               &
           &                                     fr_land_topo,            &
           &                                     hh_topo,                 &
           &                                     stdh_topo,               &
           &                                     z0_topo,                 &
           &                                     hh_topo_max=hh_topo_max, &
           &                                     hh_topo_min=hh_topo_min, &           
           &                                     theta_topo=theta_topo,   &
           &                                     aniso_topo=aniso_topo,   &
           &                                     slope_topo=slope_topo,   &
           &                                     vertex_param=vertex_param)
    ELSE
      CALL read_netcdf_buffer_topo(orography_buffer_file, &
           &                                     tg,           &
           &                                     undefined,    &
           &                                     undef_int,    &
           &                                     fr_land_topo,&
           &                                     hh_topo,     &
           &                                     stdh_topo,   &
           &                                     z0_topo,      &
           &                                     vertex_param=vertex_param)
    ENDIF

  CASE DEFAULT

     IF (lradtopo) THEN
        IF (lsso_param) THEN
           CALL read_netcdf_buffer_topo(orography_buffer_file,&
                &                                     tg,           &
                &                                     undefined,    &
                &                                     undef_int,    &
                &                                     fr_land_topo,&
                &                                     hh_topo,     &
                &                                     stdh_topo,   &
                &                                     z0_topo,      &
                &                                     lrad=lradtopo,&
                &                                     nhori=nhori,  &
                &                                     theta_topo=theta_topo,&
                &                                     aniso_topo=aniso_topo,&
                &                                     slope_topo=slope_topo,&
                &                                     slope_asp_topo=slope_asp_topo,     &
                &                                     slope_ang_topo=slope_ang_topo,     &
                &                                     horizon_topo=horizon_topo,         &
                &                                     skyview_topo=skyview_topo)
        ELSE
           CALL read_netcdf_buffer_topo(orography_buffer_file,&
                &                                     tg,           &
                &                                     undefined,    &
                &                                     undef_int,    &
                &                                     fr_land_topo,&
                &                                     hh_topo,     &
                &                                     stdh_topo,   &
                &                                     z0_topo,      &
                &                                     lrad=lradtopo,&
                &                                     nhori=nhori,  &
                &                                     slope_asp_topo=slope_asp_topo,     &
                &                                     slope_ang_topo=slope_ang_topo,     &
                &                                     horizon_topo=horizon_topo,         &
                &                                     skyview_topo=skyview_topo)
        ENDIF

     ELSE
        IF (lsso_param) THEN
           CALL read_netcdf_buffer_topo(orography_buffer_file,&
                &                                     tg,           &
                &                                     undefined,    &
                &                                     undef_int,    &
                &                                     fr_land_topo,&
                &                                     hh_topo,     &
                &                                     stdh_topo,   &
                &                                     z0_topo,      &
                &                                     nhori=nhori,  &
                &                                     theta_topo=theta_topo,  &
                &                                     aniso_topo=aniso_topo,  &
                &                                     slope_topo=slope_topo)
        ELSE
           CALL read_netcdf_buffer_topo(orography_buffer_file,&
                &                                     tg,           &
                &                                     undefined,    &
                &                                     undef_int,    &
                &                                     fr_land_topo,&
                &                                     hh_topo,     &
                &                                     stdh_topo,   &
                &                                     z0_topo,      &
                &                                     nhori=nhori)
        ENDIF
     ENDIF


  END SELECT

  IF (l_use_sgsl) THEN
     PRINT *,'Read in subgrid-scale slope data: ',TRIM(sgsl_buffer_file)
     CALL read_netcdf_buffer_sgsl(sgsl_buffer_file,  &
          &                                     tg,         &
          &                                     undefined, &
          &                                     undef_int, &
          &                                     sgsl )
  END IF


  PRINT *,'Read in aot data: ', TRIM (aot_buffer_file)
  IF (iaot_type == 4) THEN
     n_spectr = 9
     CALL read_netcdf_buffer_aot_MAC (aot_buffer_file,     &
          &                                     tg,             &
          &                                     ntype_aot,      &
          &                                     ntime_aot,      &
          &                                     n_spectr,       &
          &                                     MAC_aot_tg,     &
          &                                     MAC_ssa_tg,     &
          &                                     MAC_asy_tg)
  ELSE
     CALL read_netcdf_buffer_aot(aot_buffer_file,    &
          &                                     tg,       &
          &                                     ntype_aot,&
          &                                     ntime_aot,&
          &                                     aot_tg)
  ENDIF


  PRINT *,'Read in cru data for it_cl_type:', it_cl_type
  PRINT *,'Status ltcl_merge: ', ltcl_merge

IF (ltcl_merge) THEN
   SELECT CASE(it_cl_type)
   CASE(i_t_cru_fine)
   PRINT *,'Selected CRU Fine, ltcl_merge', ltcl_merge
     CALL read_netcdf_buffer_cru(t_clim_buffer_file,&
    &                                     tg,       &
    &                                     crutemp,  &
    &                                     cruelev)
     CALL read_netcdf_buffer_cru(t_clim_output_file, &
    &                                     tg,        &
    &                                     crutemp2)
   CASE(i_t_cru_coarse)
   PRINT *,'Selected CRU Coarse, ltcl_merge', ltcl_merge
        CALL read_netcdf_buffer_cru(t_clim_buffer_file, &
             &                                     tg,        &
             &                                     crutemp)
   END SELECT
ELSE
  SELECT CASE(it_cl_type)
  CASE(i_t_cru_fine)
    PRINT *,'Selected CRU Fine'
    CALL read_netcdf_buffer_cru(t_clim_buffer_file,&
         &                                     tg,       &
         &                                     crutemp,  &
         &                                     cruelev)  
CASE(i_t_cru_coarse)
    PRINT *,'Selected CRU Coarse'
    CALL read_netcdf_buffer_cru(t_clim_buffer_file, &
         &                                     tg,        &
         &                                     crutemp)
  END SELECT

  CALL read_netcdf_buffer_cru(t_clim_buffer_file, &
       &                                     tg,      &
       &                                     crutemp)

END IF



  PRINT *,'Read in FLAKE'
  CALL read_netcdf_buffer_flake(flake_buffer_file,   &
       &                                     tg,        &
       &                                     undefined, &
       &                                     undef_int, &
       &                                     lake_depth,&
       &                                     fr_lake,   &
       &                                     flake_tot_npixel)


  IF (i_lsm_data == 2 .and. igrid_type == igrid_cosmo) THEN
     PRINT *,'Read in Land-Sea-Mask from file  ',land_sea_mask_file
     CALL read_netcdf_buffer_lsm(land_sea_mask_file,  &
          &                           tg, &
          &                           fr_land_mask)
     PRINT *,'Land-Sea-Mask file  ',land_sea_mask_file," is used for consistency tests."
  ELSE
     PRINT *,'External  Land-Sea-Mask is NOT used for consistency tests.'
     PRINT *,'External Land-Sea-Mask is only tested for the COSMO grid.'
  END IF



  !------------------------------------------------------------------------------------------
  !------------- land use data --------------------------------------------------------------
  !------------------------------------------------------------------------------------------
  PRINT *,'determine land-sea mask'
  CALL CPU_TIME(timestart)

  IF (l_use_glcc) THEN
     WHERE (lat_geo <=lu_data_southern_boundary) ! glc2000 and globcover 2009 exclude Antarctica
        fr_land_lu = fr_land_glcc !fr_land_topo
     ENDWHERE
  ENDIF

  !set land-sea mask, spread at 0.5 due to poor accuracy of values in GRIB files
  WHERE (fr_land_lu < 0.5_wp)
     fr_land_lu = MIN(0.49_wp,fr_land_lu)
  ELSEWHERE !   fr_land_glc2000 >= 0.5
     fr_land_lu = MAX(0.51_wp,fr_land_lu)
  ENDWHERE

  IF (i_lsm_data == 2) THEN
     fr_land_lu = fr_land_mask
  ENDIF

  CALL CPU_TIME(timeend)
  timediff = timeend - timestart
  PRINT *,'determine land-sea mask, WHERE, done in: ', timediff

  ! total roughness length
  z0_tot = z0_lu + z0_topo
  IF (l_use_glcc) THEN
     WHERE (lat_geo <=lu_data_southern_boundary) ! glc2000 and globcover 2009 exclude Antarctica
        ice_lu            = ice_glcc
        z0_lu             = z0_glcc
        z0_tot            = z0_glcc + z0_topo
        root_lu           = root_glcc
        plcov_mn_lu       = plcov_mn_glcc
        plcov_mx_lu       = plcov_mx_glcc
        lai_mn_lu         = lai_mn_glcc
        lai_mx_lu         = lai_mx_glcc
        rs_min_lu         = rs_min_glcc
        urban_lu          = urban_glcc
        for_d_lu          = for_d_glcc
        for_e_lu          = for_e_glcc
        emissivity_lu     = emissivity_glcc
     ENDWHERE
  ENDIF

  ! Lower Limit for roughness length
  WHERE (fr_land_lu < 0.5_wp)
     z0_tot=MAX(1.E-6_wp,z0_tot)
  ELSEWHERE !
     z0_tot=MAX(1.E-2_wp,z0_tot)
  ENDWHERE

  !Store ice_lu with glcc (antarctica) on lu_class_fraction see mo_XXX_lookup_tables.f90
  SELECT CASE (i_landuse_data)
  CASE (i_lu_glc2000)
     lu_class_fraction(:,:,:,21)=ice_lu(:,:,:)
  CASE (i_lu_globcover)
     lu_class_fraction(:,:,:,22)=ice_lu(:,:,:)
  END SELECT

  !------------------------------------------------------------------------------------------
  !------------- land use data consistency  -------------------------------------------------
  !------------------------------------------------------------------------------------------
  IF (tile_mode < 1) THEN   ! values are kept for ICON because of tile approach
     WHERE (fr_land_lu < 0.5)  ! set vegetation to undefined (0.) for water grid elements
        ! z0 and emissivity are not set to undefined_lu
        ice_lu            = undefined_lu
        root_lu           = undefined_lu
        plcov_mn_lu       = undefined_lu
        plcov_mx_lu       = undefined_lu
        lai_mn_lu         = undefined_lu
        lai_mx_lu         = undefined_lu
        rs_min_lu         = undefined_lu
        urban_lu          = undefined_lu
        for_d_lu          = undefined_lu
        for_e_lu          = undefined_lu
     ENDWHERE
  ENDIF

  !------------------------------------------------------------------------------------------
  !------------- soil data consistency ------------------------------------------------------
  !------------------------------------------------------------------------------------------

  !undef_soiltype   = 0 ! \TODO read undef_soiltype from netcdf file (_Fill_Value)
  !default_soiltype = 5 ! default soil type loam
  !soiltype_ice     = 1   !< soiltype for ice
  !soiltype_water   = 9   !< soiltype for water
  PRINT *,'Soil data consistency check'
  CALL CPU_TIME(timestart)

  CALL define_soiltype(isoil_data, ldeep_soil, &
       undef_soiltype,         &
       default_soiltype,       &
       soiltype_ice,           &
       soiltype_water,         &
       soil_data)

  SELECT CASE (isoil_data)
  CASE(FAO_data, HWSD_map)

     WHERE (fr_land_lu < MERGE(0.01,0.5,tile_mask))  ! set water soiltype for water grid elements
        !MERGE(TSOURCE, FSOURCE, MASK) is a function which joins two arrays.
        !It gives the elements in TSOURCE if the condition in MASK is .TRUE. and FSOURCE if the condition in MASK is .FALSE.
        soiltype_fao = soiltype_water
     ELSEWHERE ! fr_land_lu >= 0.5, i.e. a land grid element
        ! (soiltyp(:,:) > 8 .OR. soiltyp(:,:) < 1))
        WHERE ((soiltype_fao == undef_soiltype).OR.(soiltype_fao > 8) ) ! land grid elements must have a valid soiltype
           !  WHERE ( (lat_geo < -60.).OR.(lat_geo > 65.) ) ! Arctic and Antarctica
           !     soiltype_fao = soiltype_ice  ! set soil type to ice for Arctic or Antarctic undefined points
           WHERE ( (lat_geo < lu_data_southern_boundary) ) ! Antarctica
              soiltype_fao = soiltype_ice  ! set soil type to ice for Antarctic undefined points

           ELSEWHERE  ! rest of the World
              soiltype_fao = default_soiltype ! set default soiltype to loam
           ENDWHERE
        ENDWHERE
     ENDWHERE

     CALL CPU_TIME(timeend)
     timediff = timeend - timestart
     PRINT *,'soil data consitency check, done in: ', timediff, ' s'

     db_ice_counter = 0

     IF (i_landuse_data == i_lu_globcover) THEN
        DO k=1,tg%ke
           DO j=1,tg%je
              DO i=1,tg%ie
                 IF  ( (soiltype_fao(i,j,k) /= soiltype_ice) .AND.  &
                                !               &    (fr_land_lu(i,j,k) > 0.5).AND. (ice_lu(i,j,k) > 0.5 )) THEN
                      &    (fr_land_lu(i,j,k)*ice_lu(i,j,k) > 0.5)) THEN ! scale glacier frac with fr_land
                    soiltype_fao(i,j,k) = soiltype_ice
                    db_ice_counter = db_ice_counter +1
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        !HA debug
        PRINT *,'number of grid elements set to ice soiltype: ', db_ice_counter
        PRINT *,'Soiltype range MIN/MAX: ', MINVAL(soiltype_fao),MAXVAL(soiltype_fao)

     ELSE   ! iI_lu_glc2000 or i_lu_glcc

        DO k=1,tg%ke
           DO j=1,tg%je
              DO i=1,tg%ie
                 IF  ( (soiltype_fao(i,j,k) /= soiltype_ice) .AND.  &
                      &    (fr_land_lu(i,j,k) > 0.5).AND. (ice_lu(i,j,k) ==  fr_land_lu(i,j,k)) ) THEN
                    soiltype_fao(i,j,k) = soiltype_ice
                    db_ice_counter = db_ice_counter +1
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        !HA debug
        PRINT *,'number of grid elements set to ice soiltype: ', db_ice_counter
     ENDIF

  CASE(HWSD_data)

     IF (ldeep_soil) THEN

        CALL calculate_soiltype(tg,            &
             &                          ldeep_soil,     &
             &                          soiltype_deep,  &
             &                          soiltype_HWSD_s,  &
             &                          fr_sand,       &
             &                          fr_silt,       &
             &                          fr_clay,       &
             &                          fr_oc,         &
             &                          fr_bd,         &
             &                          fr_sand_deep,  &
             &                          fr_silt_deep,  &
             &                          fr_clay_deep,  &
             &                          fr_oc_deep,    &
             &                          fr_bd_deep     )
        !     &                          soiltype_deep = soiltype_deep)
        print*, 'MIN/MAX soiltype_FAO_deep : ', MINVAL(soiltype_deep), MAXVAL(soiltype_deep)
        print*, 'MIN/MAX soiltype_HWSD_deep : ', MINVAL(soiltype_hwsd_s), MAXVAL(soiltype_hwsd_s)
     END IF
     print*, 'calculate_soiltype for top soil: '

     CALL calculate_soiltype(tg,            &
          &                          .false.,       & ! switch off deep soil for top soil calculation
          &                          soiltype_FAO,  &
          &                          soiltype_HWSD,  &
          &                          fr_sand,       &
          &                          fr_silt,       &
          &                          fr_clay,       &
          &                          fr_oc,         &
          &                          fr_bd          )
     print*, 'MIN/MAX soiltype_FAO top: ', MINVAL(soiltype_fao), MAXVAL(soiltype_fao)
     print*, 'MIN/MAX soiltype_HWSD top : ', MINVAL(soiltype_hwsd), MAXVAL(soiltype_hwsd)

     ! Use land-use data for setting glacier points to soiltype ice

     db_ice_counter = 0

     IF (i_landuse_data == i_lu_globcover) THEN
        DO k=1,tg%ke
           DO j=1,tg%je
              DO i=1,tg%ie
                 IF  ( (soiltype_fao(i,j,k) /= soiltype_ice) .AND.  &
                      &    (fr_land_lu(i,j,k) > 0.5).AND. (ice_lu(i,j,k) > 0.5 )) THEN
                    soiltype_fao(i,j,k) = soiltype_ice
                    db_ice_counter = db_ice_counter +1
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        !HA debug
        PRINT *,'number of grid elements set to ice soiltype: ', db_ice_counter
        PRINT *,'Soiltype range MIN/MAX: ', MINVAL(soiltype_fao),MAXVAL(soiltype_fao)

     ELSE   ! iI_lu_glc2000 or i_lu_glcc

        DO k=1,tg%ke
           DO j=1,tg%je
              DO i=1,tg%ie
                 IF  ( (soiltype_fao(i,j,k) /= soiltype_ice) .AND.  &
                      &    (fr_land_lu(i,j,k) > 0.5).AND. (ice_lu(i,j,k) ==  fr_land_lu(i,j,k)) ) THEN
                    soiltype_fao(i,j,k) = soiltype_ice
                    db_ice_counter = db_ice_counter +1
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        !HA debug
        PRINT *,'number of grid elements set to ice soiltype: ', db_ice_counter
     ENDIF

     ! Final check for special cases - Dunes to sand, antarctica to ice, cities to loam

     WHERE ( (lat_geo < lu_data_southern_boundary) ) ! Antarctica
        soiltype_fao = soiltype_ice  ! set soil type to ice for Antarctic undefined points
     ENDWHERE

     IF(igrid_type == igrid_icon) THEN
        WHERE (fr_land_lu < MERGE(0.01,0.5,tile_mask))  ! set water soiltype for water grid elements
           !MERGE(TSOURCE, FSOURCE, MASK) is a function which joins two arrays.
           !It gives the elements in TSOURCE if the condition in MASK is .TRUE. and
           !FSOURCE if the condition in MASK is .FALSE.
           soiltype_fao = soiltype_water
        ENDWHERE
     ELSE
        WHERE (fr_land_lu < 0.5)  ! set water soiltype for water grid elements
           !MERGE(TSOURCE, FSOURCE, MASK) is a function which joins two arrays.
           !It gives the elements in TSOURCE if the condition in MASK is .TRUE. and FSOURCE if the condition in MASK is .FALSE.
           soiltype_fao = soiltype_water
        ENDWHERE
        IF (ldeep_soil) THEN
           WHERE (fr_land_lu < 0.5)  ! set water soiltype for water grid elements
              soiltype_deep = soiltype_water
           ENDWHERE
        END IF
     END IF

     !Consider Land-points with soiltype water
     db_water_counter=0
     DO k=1,tg%ke
        DO j=1,tg%je
           DO i=1,tg%ie
              IF  ( (soiltype_fao(i,j,k) == soiltype_water) .AND.  &
                   &    (fr_land_lu(i,j,k) > 0.5) ) THEN
                 soiltype_fao(i,j,k) = 5
                 db_water_counter = db_water_counter +1
              ENDIF
           ENDDO
        ENDDO
     ENDDO
     !HA debug
     PRINT *,'number of land points with water set to  soiltype loam: ', db_water_counter

     WHERE (soiltype_fao == 11) ! Dunes
        soiltype_fao = 3  ! set soil type to sand for dunes
     ENDWHERE
     WHERE (soiltype_fao > 12) ! undefined
        soiltype_fao = 5  ! set soil type to loam for undefined points
     ENDWHERE


     CALL CPU_TIME(timeend)
     timediff = timeend - timestart
     PRINT *,'soil data consitency check, WHERE, done in: ', timediff, ' s'


  END SELECT

  !------------------------------------------------------------------------------------------

  !    SELECT CASE (isoil_data)
  !    CASE(FAO_data)
  WHERE (soiltype_fao == soiltype_ice)  ! set vegetation to undefined (0.) for ice grid elements (e.g. glaciers)
     ! z0, rs_min and emissivity are not set to undefined_lu
     root_lu           = undefined_lu
     plcov_mn_lu       = undefined_lu
     plcov_mx_lu       = undefined_lu
     lai_mn_lu         = undefined_lu
     lai_mx_lu         = undefined_lu
     urban_lu          = undefined_lu
     for_d_lu          = undefined_lu
     for_e_lu          = undefined_lu
  ENDWHERE
  !    CASE(HWSD_data)
  WHERE (soiltype_hwsd == 1)  ! set vegetation to undefined (0.) for ice grid elements (e.g. glaciers)
     ! z0, rs_min and emissivity are not set to undefined_lu
     root_lu           = undefined_lu
     plcov_mn_lu       = undefined_lu
     plcov_mx_lu       = undefined_lu
     lai_mn_lu         = undefined_lu
     lai_mx_lu         = undefined_lu
     urban_lu          = undefined_lu
     for_d_lu          = undefined_lu
     for_e_lu          = undefined_lu
  ENDWHERE

  !    END SELECT
  !------------------------------------------------------------------------------------------
  !------------- soil data consistency ------------------------------------------------------
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  !------------- flake data consistency  ----------------------------------------------------
  !------------------------------------------------------------------------------------------
  !# Comment Merge Conflicts in flake consistency between COSMO and DWD#
  SELECT CASE(igrid_type) ! get indices for neighbour grid elements

  CASE(igrid_icon) ! ICON GRID

     PRINT *,'flake data consistency check'
     CALL CPU_TIME(timestart)

     ! determine "fraction ocean" first before considering "fraction lake"
     ! fr_ocean should be determined by ocean model if available
     ! so use (1. - lsm_ocean_model) as mask instead of fr_land_topo from the orography data
     thr_cr = 0.99
     WHERE (fr_land_topo < thr_cr)
        fr_ocean_lu = 1. - fr_land_lu
        fr_lake = 0.0
     ELSEWHERE
        fr_ocean_lu = 0.0
        fr_lake = 1. - fr_land_lu
     ENDWHERE

     ! set Death Sea to "ocean water"
     WHERE ((hh_topo < -390.).AND. &
          &     (lon_geo > 35.).AND.(lon_geo < 36.).AND. &
          &     (lat_geo > 31.).AND.(lat_geo < 32.) )
        fr_ocean_lu = 1. - fr_land_lu
        fr_lake = 0.0
     ENDWHERE

     ! set Caspian Sea to "ocean water"
     WHERE ((hh_topo < -25.).AND. &
          &     (lon_geo > 46.).AND.(lon_geo < 55.).AND. &
          &     (lat_geo > 36.).AND.(lat_geo < 48.) )
        fr_ocean_lu = 1. - fr_land_lu
        fr_lake = 0.0
     ENDWHERE
     ! here fr_ocean_lu + fr_lake +fr_land_lu = 1
     ! fr_ocean_lu + fr_lake = fr_water
     ! fr_water + fr_land_lu = 1

     ! check consistency for "lake depth"
     IF (tile_mode == 1) THEN ! subgrid lakes for ICON
        WHERE (fr_land_lu >= thr_cr ) ! 0.5)
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE
        WHERE (fr_ocean_lu >= 0.5)
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE

        WHERE ((fr_lake > 1.-thr_cr).AND.(lake_depth < 0.0)) ! fr_lake > 0.5
           lake_depth = flake_depth_default ! set lake depth to default value (10 m)
        ENDWHERE !

     ELSE
        WHERE (fr_land_lu >= 0.5 )
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE
        WHERE (fr_ocean_lu >= 0.5)
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE

        WHERE ((fr_lake > 0.5).AND.(lake_depth < 0.0))
           lake_depth = flake_depth_default ! set lake depth to default value (10 m)
        ENDWHERE !

     ENDIF
     ! restrict lake depth to maximum value (50 m)
     WHERE (lake_depth > DWD_max_lake_depth)
        lake_depth = DWD_max_lake_depth
     END WHERE

     ! restrict lake depth to minimal value (1 m)
     WHERE ( (lake_depth > 0.0).AND.(lake_depth < DWD_min_lake_depth ))
        lake_depth = DWD_min_lake_depth
     END WHERE

     DO nloops=1,3
        DO k=1,tg%ke
           DO j=1,tg%je
              DO i=1,tg%ie

                 IF (fr_lake(i,j,k)>0.05) THEN ! concistency check for neighbour ocean elements
                    ! get neighbour grid indices for ICON grid
                    ne_je(:) = 1_i8
                    ne_ke(:) = 1_i8
                    ne_ie(:) = 0_i8
                    nnb=icon_grid%nvertex_per_cell ! number of neighbours in ICON grid
                    DO nv=1, nnb
                       n_index = icon_grid_region%cells%neighbor_index(i,nv) ! get cell id of neighbour cells
                       IF (n_index > 0) THEN
                          ne_ie(nv) = n_index
                       ENDIF
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     ! manual correction for lake depth of Lake Constance
     ! the raw database of the lake depth data appears not to be correct
     ! the main part of Lake Constance has a mean depth of 98 m, so set to DWD_max_lake_depth
     WHERE ((lon_geo > 8.8).AND.(lon_geo < 9.9).AND. &
          &     (lat_geo > 47.4).AND.(lat_geo < 48.4).AND. &
          &     (lake_depth > 9.7).AND.(lake_depth < 9.9) )
        lake_depth = DWD_max_lake_depth
     ENDWHERE

  CASE(igrid_cosmo) ! COSMO grid

     PRINT *,'flake data consistency check'
     CALL CPU_TIME(timestart)

     ! determine "fraction ocean" first before considering "fraction lake"
     ! fr_ocean should be determined by ocean model if available
     ! so use (1. - lsm_ocean_model) as mask instead of fr_land_topo from the orography data

     ! set surface height of all Dead Sea points to the level given in the
     ! repective data_set, i. e. -405 m (GLOBE) and -432 m (ASTER)
     hh_dead_sea = -405._wp
     IF (itopo_type ==2) hh_dead_sea = -432._wp
     WHERE ((lon_geo > 35.).AND.(lon_geo < 36.).AND. &
          &     (lat_geo > 31.).AND.(lat_geo < 32.).AND. &
          &     (1._wp - fr_land_lu > 0.5_wp))
        hh_topo = hh_dead_sea
        fr_lake = 0.0_wp
        fr_ocean_lu = 1._wp - fr_land_lu
     ENDWHERE
     WHERE ((lon_geo > 35.).AND.(lon_geo < 36.).AND. &
          &     (lat_geo > 31.).AND.(lat_geo < 32.))
        hh_topo = MAX(hh_topo, hh_dead_sea)
        fr_ocean_lu = fr_ocean_lu + fr_lake
        fr_lake = 0.0_wp
        fr_land_lu = 1._wp - fr_ocean_lu
     ENDWHERE

     hh_cr_casp = -25._wp
     ! set surface height of all Caspian Sea points to -28. meters
     WHERE ((lon_geo > 46.).AND.(lon_geo < 55.).AND. &
          &     (lat_geo > 36.).AND.(lat_geo < 49.))
        fr_ocean_lu = 1. - fr_land_lu - fr_lake
        fr_land_topo = -1._wp
     ENDWHERE
     WHERE ((lon_geo > 45.2).AND.(lon_geo < 48.8).AND. &
          &     (lat_geo > 45.9).AND.(lat_geo < 50.0))
        fr_ocean_lu = 0.0
        fr_lake = 1._wp - fr_land_lu
        fr_land_topo = 1._wp
     ENDWHERE
     WHERE ((lon_geo > 48.8).AND.(lon_geo < 52.9).AND. &
          &     (lat_geo > 47.2).AND.(lat_geo < 50.0))
        fr_ocean_lu = 0.0
        fr_lake = 1._wp - fr_land_lu
        fr_land_topo = 1._wp
     ENDWHERE
     WHERE ((lon_geo > 52.9).AND.(lon_geo < 55.0).AND. &
          &     (lat_geo > 47.2).AND.(lat_geo < 50.0))
        fr_ocean_lu = 0.0
        fr_lake = 1._wp - fr_land_lu
        fr_land_topo = 1._wp
     ENDWHERE
     WHERE ((lon_geo > 53.4).AND.(lon_geo < 56.7).AND. &
          &     (lat_geo > 45.8).AND.(lat_geo < 47.0))
        fr_ocean_lu = 0.0
        fr_lake = 1._wp - fr_land_lu
        fr_land_topo = 1._wp
     ENDWHERE
     WHERE ((lon_geo > 45.8).AND.(lon_geo < 48.6).AND. &
          &     (lat_geo > 39.5).AND.(lat_geo < 41.4))
        fr_ocean_lu = 0.0
        fr_lake = 1._wp - fr_land_lu
        fr_land_topo = 1._wp
     ENDWHERE

     thr_cr = 0.99
     WHERE ((fr_land_topo < thr_cr).AND.(fr_land_topo >= 0._wp))
        fr_ocean_lu = 1. - fr_land_lu
        fr_lake = 0.0
     ENDWHERE
     WHERE ((fr_land_topo >= thr_cr))
        fr_ocean_lu = 0.0
        fr_lake = 1. - fr_land_lu
     ENDWHERE

     ! check consistency for "lake depth"
     IF (tile_mode == 1) THEN ! subgrid lakes for ICON
        WHERE (fr_land_lu >= thr_cr ) ! 0.5)
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE
        WHERE (fr_ocean_lu >= 0.5)
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE
        WHERE ((fr_lake > 1.-thr_cr).AND.(lake_depth < 0.0)) ! fr_lake > 0.5
           lake_depth = flake_depth_default ! set lake depth to default value (10 m)
        ENDWHERE !
     ELSE
        WHERE (fr_land_lu >= 0.5)
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE
        WHERE (fr_ocean_lu >= 0.5)
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE
        WHERE (fr_lake < 0.5)
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE
        WHERE ((fr_lake > 0.5).AND.(lake_depth < 0.0))
           lake_depth = flake_depth_default ! set lake depth to default value (10 m)
        ENDWHERE !
     ENDIF
     ! restrict lake depth to maximum value (50 m)
     WHERE (lake_depth > DWD_max_lake_depth)
        lake_depth = DWD_max_lake_depth
     END WHERE

     ! restrict lake depth to minimal value (1 m)
     WHERE ( (lake_depth > 0.0).AND.(lake_depth < DWD_min_lake_depth ))
        lake_depth = DWD_min_lake_depth
     END WHERE

     DO nloops=1,3
        DO k=1,tg%ke
           DO j=1,tg%je
              DO i=1,tg%ie

                 IF (fr_lake(i,j,k)>0.05) THEN ! concistency check for neighbour ocean elements

                    nnb = 8
                    ! northern neighbour
                    ne_ie(1) = i
                    ne_je(1) = MAX(1_i8,j-1)
                    ne_ke(1) = k
                    ! north-eastern neighbour
                    ne_ie(2) = MIN(tg%ie,INT(i+1,i8))
                    ne_je(2) = MAX(1_i8,j-1)
                    ne_ke(2) = k
                    ! eastern neighbour
                    ne_ie(3) = MIN(tg%ie,INT(i+1,i8))
                    ne_je(3) = j
                    ne_ke(3) = k
                    ! south-eastern neighbour
                    ne_ie(4) = MIN(tg%ie,INT(i+1,i8))
                    ne_je(4) = MIN(tg%je,INT(j+1,i8))
                    ne_ke(4) = k
                    ! southern neighbour
                    ne_ie(5) = i
                    ne_je(5) = MIN(tg%je,INT(j+1,i8))
                    ne_ke(5) = k
                    ! south-west neighbour
                    ne_ie(6) = MAX(1_i8,i-1)
                    ne_je(6) = MIN(tg%je,INT(j+1,i8))
                    ne_ke(6) = k
                    ! western neighbour
                    ne_ie(7) = MAX(1_i8,i-1)
                    ne_je(7) = j
                    ne_ke(7) = k
                    ! north-west neighbour
                    ne_ie(8) = MAX(1_i8,i-1)
                    ne_je(8) = MAX(1_i8,j-1)
                    ne_ke(8) = k

                    IF (lflake_correction) THEN
                       DO n=1,nnb
                          IF ((ne_ie(n)>= 1).AND.(ne_je(n)>=1).AND.(ne_ke(n)>=1)) THEN
                             IF (fr_ocean_lu(ne_ie(n),ne_je(n),ne_ke(n))>0.5) THEN ! if the direct neighbour element is ocean,
                                fr_lake(i,j,k) = 0.0                                ! set this grid element also to ocean.
                                IF ((i==391).AND.(j==267)) PRINT *,'changed: ',                    &
                                     ne_ie(n),ne_je(n),ne_ke(n),fr_ocean_lu(ne_ie(n),ne_je(n),ne_ke(n))
                                fr_ocean_lu(i,j,k) = 1.0 - fr_land_lu(i,j,k)
                                lake_depth(i,j,k) = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
                             ENDIF
                          ENDIF
                       ENDDO
                    ENDIF
                 ENDIF ! check for ocean
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     ! manual correction for lake depth of Lake Constance
     ! the raw database of the lake depth data appears not to be correct
     ! the main part of Lake Constance has a mean depth of 98 m, so set to DWD_max_lake_depth
     WHERE ((lon_geo > 8.8).AND.(lon_geo < 9.9).AND. &
          &     (lat_geo > 47.4).AND.(lat_geo < 48.4).AND. &
          &     (lake_depth > 9.7).AND.(lake_depth < 9.9) )
        lake_depth = DWD_max_lake_depth
     ENDWHERE

     ! check consistency for "lake depth" again
     WHERE (fr_lake >= fr_ocean_lu)
        fr_lake = 1.0_wp - fr_land_lu
        fr_ocean_lu = 0.0_wp
     ELSEWHERE
        fr_ocean_lu = 1.0_wp - fr_land_lu
        fr_lake = 0.0_wp
     ENDWHERE

     IF (tile_mode == 1) THEN ! subgrid lakes for ICON
        WHERE (fr_land_lu >= thr_cr ) ! 0.5)
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE
        WHERE (fr_ocean_lu >= 0.5)
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE
        WHERE ((fr_lake > 1.-thr_cr).AND.(lake_depth < 0.0)) ! fr_lake > 0.5
           lake_depth = flake_depth_default ! set lake depth to default value (10 m)
        ENDWHERE !
     ELSE
        WHERE (fr_land_lu >= 0.5)
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE
        WHERE (fr_ocean_lu >= 0.5)
           lake_depth = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
        ENDWHERE
        WHERE ((fr_lake < 0.5))
           lake_depth = flake_depth_undef !  set lake depth to flake_depth_undef (-1 m)
        ENDWHERE !
        WHERE ((fr_lake > 0.5).AND.(lake_depth < DWD_min_lake_depth))
           lake_depth = flake_depth_default ! set lake depth to default value (10 m)
        ENDWHERE !
     ENDIF
     ! restrict lake depth to maximum value (50 m)
     WHERE (lake_depth > DWD_max_lake_depth)
        lake_depth = DWD_max_lake_depth
     END WHERE

     ! restrict lake depth to minimal value (1 m)
     WHERE ( (lake_depth > 0.0).AND.(lake_depth < DWD_min_lake_depth ))
        lake_depth = DWD_min_lake_depth
     END WHERE

     ! adjust surface height of Caspian sea to -28 m
     WHERE ((lon_geo > 46.).AND.(lon_geo < 55.).AND. &
          &     (lat_geo > 36.).AND.(lat_geo < 48.).AND. &
          &     (fr_ocean_lu > 0.5))
        hh_topo = -28.0_wp
     ENDWHERE

  END SELECT




  !------------------------------------------------------------------------------------------
  !------------- Albedo data consistency ------------------------------------------------------
  !------------------------------------------------------------------------------------------
  ! set default Albedo values for land grid elements with so far undefined or unexpected values
  IF (ialb_type /= 2) THEN

     ! set default Albedo values for land grid elements with so far undefined or unexpected values
     PRINT *,'Albedo data consistency check'

     CALL CPU_TIME(timestart)


     namelist_alb_data_input = 'INPUT_ALB'

     CALL  read_namelists_extpar_alb(namelist_alb_data_input, &
          &                                  raw_data_alb_path, &
          &                                  raw_data_alb_filename, &
          &                                  raw_data_alnid_filename, &
          &                                  raw_data_aluvd_filename, &
          &                                  ialb_type, &
          &                                  alb_buffer_file, &
          &                                  alb_output_file, &
          &                                  alb_source, &
          &                                  alnid_source, &
          &                                  aluvd_source)

     path_alb_file = TRIM(raw_data_alb_path)//TRIM(raw_data_alb_filename)
     print *, TRIM(path_alb_file)

     nlon_reg = alb_raw_data_grid%nlon_reg
     nlat_reg = alb_raw_data_grid%nlat_reg

     path_alb_file = TRIM(raw_data_alb_path)//TRIM(raw_data_alb_filename)
     print *, TRIM(path_alb_file)

     CALL open_netcdf_ALB_data(path_alb_file, &
          ncid_alb)

     CALL allocate_alb_interp_fields(mpy)

     CALL alb_interp_data()

     DO i=1,9
        step = 0.1667*(wso_max(i) - wso_min(i))
        zalso(i,7) = csalb(i) - wso_min(i)*csalbw(i)
        zalso(i,1) = csalb(i) - wso_max(i)*csalbw(i)
        DO t=1,5
           zalso(i,7-t) = csalb(i) - csalbw(i)*step*t
           zalso(i,7+t) = csalb(i) - csalbw(i)*step*t
        ENDDO
     ENDDO

     albvis_min=0.07
     CALL const_check_interpol_alb(alb_field_mom,fr_land_lu,albvis_min)

     albnir_min=0.10
     CALL const_check_interpol_alb(alnid_field_mom,fr_land_lu,albnir_min)

     albuv_min=0.02
     CALL const_check_interpol_alb(aluvd_field_mom,fr_land_lu,albuv_min)


     CALL CPU_TIME(timeend)
     timediff = timeend - timestart
     PRINT *,'albedo data consistency check, WHERE, done in: ', timediff
  ENDIF

  !------------------------------------------------------------------------------------------
  !------------- soil albedo consistency check ----------------------------------------------
  !------------------------------------------------------------------------------------------
  IF (ialb_type == 2) THEN

     ! set default soil albedo values for land grid elements with so far undefined values
     PRINT *,'soil albedo data consistency check'

     CALL CPU_TIME(timestart)

     WHERE (fr_land_lu < 0.5) ! set undefined albedo value (0.0) for water grid elements
        alb_dry = undef_alb_bs
        alb_sat = undef_alb_bs
     ELSEWHERE ! fr_land_lu >= 0.5
        WHERE (alb_dry < minimal_alb_dry) ! small albedo values at land grid elements
           alb_dry = minimal_alb_dry
        ENDWHERE
        WHERE (alb_dry > maximal_alb_dry) ! large albedo values at land grid elements
           alb_dry = maximal_alb_dry
        ENDWHERE
        WHERE (alb_sat < minimal_alb_sat) ! small albedo values at land grid elements
           alb_sat = minimal_alb_sat
        ENDWHERE
        WHERE (alb_sat > maximal_alb_sat) ! large albedo values at land grid elements
           alb_sat = maximal_alb_sat
        ENDWHERE
     ENDWHERE


     CALL CPU_TIME(timeend)
     timediff = timeend - timestart
     PRINT *,'albedo data consistency check, WHERE, done in: ', timediff

  ENDIF

  !------------------------------------------------------------------------------------------
  !------------- ISA/AHF data consistency ---------------------------------------------------
  !------------------------------------------------------------------------------------------

  IF (l_use_isa.AND.l_use_ahf) THEN

     ! set default ISA/AHF values for land grid elements with so far undefined values or very small NDVI values
     PRINT *,'ISA/AHF data consistency check'

     CALL CPU_TIME(timestart)
     !minimal_ndvi = 0.09 ! bare soil value
     !undef_ndvi   = 0.0  ! no vegetation


     WHERE (fr_land_lu < MERGE(0.01,0.5,tile_mask)) ! set undefined ISA value (0.0) for water grid elements
        isa_field = undef_isa
     ELSEWHERE ! fr_land_lu >= 0.5
        WHERE (isa_field <= minimal_isa) ! small ISA values at land grid elements
           isa_field = minimal_isa
        ENDWHERE
     ENDWHERE


     WHERE (fr_land_lu < MERGE(0.01,0.5,tile_mask)) ! set undefined AHF value (0.0) for water grid elements
        ahf_field = undef_ahf
     ELSEWHERE ! fr_land_lu >= 0.5
        WHERE (ahf_field <= minimal_ahf) ! small AHF values at land grid elements
           ahf_field = minimal_ahf
        ENDWHERE
     ENDWHERE

     PRINT *,'Urban data consistency check'

     WHERE (fr_land_lu < 0.5)  ! set water soiltype for water grid elements
        isa_field=0.
     ENDWHERE
  END IF
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  !------------- NDVI data consistency ------------------------------------------------------
  !------------------------------------------------------------------------------------------

  ! set default NDVI values for land grid elements with so far undefined values or very small NDVI values
  PRINT *,'NDVI data consistency check'

  CALL CPU_TIME(timestart)
  !minimal_ndvi = 0.09 ! bare soil value
  !undef_ndvi   = 0.0  ! no vegetation

  WHERE (fr_land_lu < MERGE(0.01,0.5,tile_mask)) ! set undefined NDVI value (0.0) for water grid elements
     ndvi_max = undef_ndvi
  ELSEWHERE ! fr_land_lu >= 0.5
     WHERE (ndvi_max <= minimal_ndvi) ! small NDVI values at land grid elements
        ndvi_max = minimal_ndvi
     ENDWHERE
  ENDWHERE

  FORALL (t=1:mpy) ! mpy = month per year = 12
     WHERE (fr_land_lu < MERGE(0.01,0.5,tile_mask)) ! set undefined NDVI value (0.0) for water grid elements
        ndvi_field_mom(:,:,:,t) = undef_ndvi
        ndvi_ratio_mom(:,:,:,t) = undef_ndvi
     ELSEWHERE ! fr_land_lu >= 0.5
        WHERE (ndvi_max(:,:,:) <= minimal_ndvi) ! small NDVI values at land grid elements
           ndvi_field_mom(:,:,:,t) = minimal_ndvi
           ndvi_ratio_mom(:,:,:,t) = 1.0  ! minimal_ndvi / minimal_ndvi ! ndvi_max set to minimal_ndvi
        ENDWHERE
        WHERE (ndvi_field_mom(:,:,:,t) <= minimal_ndvi) ! small NDVI values at land grid elements
           ndvi_field_mom(:,:,:,t) = minimal_ndvi
           ndvi_ratio_mom(:,:,:,t) = minimal_ndvi / ndvi_max(:,:,:)
        ENDWHERE
     ENDWHERE
  END FORALL

  CALL CPU_TIME(timeend)
  timediff = timeend - timestart
  PRINT *,'NDVI data consitency check, WHERE, done in:  ', timediff


  !------------------------------------------------------------------------------------------
  !------------- NDVI data consistency ------------------------------------------------------
  !------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------
  !-------------TC_L Correction ------------------------------------------------------
  !------------------------------------------------------------------------------------------
  !#Comment from Merge: Check this section between COSMO and DWD!
IF (ltcl_merge) THEN
     PRINT*,'T_CL Merging of Coarse and Fine' 
    DO j=1,tg%je
      DO i=1,tg%ie
          IF ( crutemp(i,j,1) > 0.0 ) THEN  ! Fine
            crutemp(i,j,1) = crutemp(i,j,1) + 0.65 * 0.01*( cruelev(i,j,1) - hh_topo(i,j,1))
            ELSE
            crutemp(i,j,1) = crutemp2(i,j,1) ! Coarse
         END IF
      END DO
    END DO
ELSE 
  !gs
  IF (igrid_type == igrid_cosmo) THEN
    SELECT CASE(it_cl_type)
    CASE(i_t_cru_fine)

      PRINT*,'T_CL Correction'
      crutemp2 = crutemp
      DO j=1,tg%je
        DO i=1,tg%ie
          last = .FALSE.
          IF ( fr_land_lu(i,j,1) < 0.5) THEN
            crutemp(i,j,1)  = -1.E20_wp
          ELSE
            IF ( crutemp(i,j,1) > 0.0 ) THEN
              foundtcl = .TRUE.
              !   PRINT*, 'CRUTEMP', i,j, crutemp(i,j,1)
              !   PRINT*, 'ELEV DOMAIN', hh_topo(i,j,1)
              !   PRINT*, 'ELEV CRU', cruelev(i,j,1)

              crutemp(i,j,1) = crutemp(i,j,1) + 0.65 * 0.01*( cruelev(i,j,1) - hh_topo(i,j,1) )
              !    PRINT*, 'CRUTEMP', crutemp(i,j,1)

            ELSE
              ! PRINT*, 'TCL NOT DEFINED ',tg%ie, tg%je, i, j
              ! 3x3 search
              foundtcl = .FALSE.
              DO jj=-1,2
                DO ii=-1,2
                  IF (j+jj > 0 .and.  j+jj < tg%je .and. i+ii > 0 .and. i+ii < tg%ie) THEN
                    IF ( crutemp2(i+ii,j+jj,1) > 0.0 ) THEN
                      ! PRINT*, 'FOUND TCL ', i, j, ii, jj, crutemp2(i+ii,j+jj,1)
                      crutemp(i,j,1) = crutemp2(i+ii,j+jj,1) + 0.65 * 0.01*( cruelev(i+ii,j+jj,1) - hh_topo(i,j,1) )
                      foundtcl = .TRUE.
                      last = .TRUE.
                      exit
                    END IF
                  ENDIF   ! inside domain
                END DO
                IF  (last) THEN
                  exit
                ENDIF
              END DO
              ! if still missing go along longitude
              IF (.NOT. foundtcl) THEN
                tclsum = 0._wp
                elesum = 0._wp
                ntclct = 0
                l = 1
                DO WHILE (.NOT. foundtcl .AND. l .le. (tg%ie / 6))
                  iml=MAX(1_i8,i-3*l)
                  imu=i+2-3*l
                  ipl=i+3*l-2
                  ipu=MIN(tg%ie,i+3*l)
                  jml=MAX(1_i8,j-l)
                  jmu=j-l
                  jpl=j+l
                  jpu=MIN(tg%je,INT(j+l,i8))
                  IF (jml == jmu) THEN
                    DO ii=iml,ipu
                      IF ( crutemp2(ii,jml,1) > 0.0 ) THEN
                        tclsum = tclsum + crutemp2(ii,jml,1)
                        elesum = elesum + cruelev(ii,jml,1)
                        ntclct = ntclct + 1
                      ENDIF
                    ENDDO
                  ELSE
                    jml = jml - 1
                  ENDIF
                  IF (jpl == jpu) THEN
                    DO ii=iml,ipu
                      IF ( crutemp2(ii,jpu,1) > 0.0 ) THEN
                        tclsum = tclsum + crutemp2(ii,jpu,1)
                        elesum = elesum + cruelev(ii,jpu,1)
                        ntclct = ntclct + 1
                      ENDIF
                    ENDDO
                  ELSE
                    jpu = jpu + 1
                  ENDIF
                  IF (iml .LE. imu) THEN
                    DO jj = jml+1,jpu-1
                      DO ii = iml,imu
                        IF ( crutemp2(ii,jj,1) > 0.0 ) THEN
                          tclsum = tclsum + crutemp2(ii,jj,1)
                          elesum = elesum + cruelev(ii,jj,1)
                          ntclct = ntclct + 1
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDIF
                  IF (ipl .LE. ipu) THEN
                    DO jj = jml+1,jpu-1
                      DO ii = ipl,ipu
                        IF ( crutemp2(ii,jj,1) > 0.0 ) THEN
                          tclsum = tclsum + crutemp2(ii,jj,1)
                          elesum = elesum + cruelev(ii,jj,1)
                          ntclct = ntclct + 1
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDIF
                  IF (ntclct > 0) THEN
                    crutemp(i,j,1) = tclsum/REAL(ntclct,wp) + 0.0065 * (elesum/REAL(ntclct,wp) - hh_topo(i,j,1))
                    foundtcl = .TRUE.
                  ELSE
                    l = l + 1
                  ENDIF
                ENDDO  ! while
              ENDIF    ! .not. foundtcl

              IF ( .NOT. foundtcl) THEN

                PRINT*, 'ERROR NO TEMPERATURE DATA FOR T_CL CORRECTION  AT'
                PRINT *,i,j
                crutemp(i,j,1) = 288.15 - 0.0065 * hh_topo(i,j,1)

              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
   END SELECT
 ENDIF
END IF

  !------------------------------------------------------------------------------------------

  SELECT CASE(isoil_data)

  CASE(HWSD_data)
     PRINT *,'Selected HWSD - Copy HWSD data for Output'
     !   soiltype_fao=soiltype_hwsd
  END SELECT


  !------------- Special Points        ------------------------------------------------------

  !   SELECT CASE(igrid_type)
  !           CASE(igrid_cosmo)) ! ICON GRID ! COSMO grid
  IF (igrid_type == igrid_cosmo .OR. igrid_type == igrid_icon) THEN


     do isp = 1, number_special_points
        IF (number_special_points<1) THEN
           write(*,*) 'No treatment of special points: Number of special points is ',number_special_points
           EXIT
        END IF

        !------------------------------------------------------------------------------------------
        ! CAUTION! Tested only for COSMO!!!

        !--------------------------------------------------------------------------------------------------------
        ! get file

        WRITE(namelist_file,'(A9,I1)') 'INPUT_SP_',isp

        !---------------------------------------------------------------------------

        CALL read_namelists_extpar_special_points(namelist_file, &
             lon_geo_sp,           &
             lat_geo_sp,           &
             soiltype_sp,          &
             z0_sp,                &
             rootdp_sp,            &
             plcovmn_sp,           &
             plcovmx_sp,           &
             laimn_sp,             &
             laimx_sp,             &
             for_d_sp,             &
             for_e_sp,             &
             fr_land_sp            )


        ! Consider only well defined variables, default is -999.!
        IF ((lon_geo_sp < -360._wp ).OR.(lat_geo_sp < -90._wp)) THEN
           PRINT*,"CAUTION! Special points defined but not in target domain!"
        ELSE
           start_cell_id = 1

           CALL  find_nearest_target_grid_element(lon_geo_sp, &
                & lat_geo_sp, &
                & tg,            &
                & start_cell_id, &
                & i_sp,      &
                & j_sp,      &
                & k_sp)

           WRITE(*,'(A)') '-------------------------------------------------------------------------------------'
           WRITE(*,'(A26,A10,A4,2X,I1)')  "Consider special point in ",namelist_file," of ",number_special_points
           WRITE(*,'(A33,1X,2(F6.3,2X))') "         special point position (lon,lat): ",lon_geo_sp,lat_geo_sp
           WRITE(*,'(A33,1X,3(I9,2X))')   "         special point index (ie,je,ke):   ",i_sp,j_sp,k_sp
           IF ((i_sp == 0).OR.(j_sp == 0)) THEN
              PRINT*,"CAUTION! Special points out of range of target domain!"
           ELSE
              WRITE(*,'(A23,I9,2X,I7,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                   i_sp,j_sp," z0_tot old ",z0_tot (i_sp,j_sp,k_sp),"new ",z0_sp
              WRITE(*,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                   i_sp,j_sp," root_lu old ",root_lu(i_sp,j_sp,k_sp),"new ",rootdp_sp
              WRITE(*,'(A23,I9,2X,I9,A19,I5,2X,A4,I5)')"         special point: ",&
                   i_sp,j_sp," soiltype_fao old  ",soiltype_fao(i_sp,j_sp,k_sp),"new ",NINT(soiltype_sp)
              WRITE(*,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                   i_sp,j_sp," plcov_mn_lu  old  ",plcov_mn_lu (i_sp,j_sp,k_sp),"new ",plcovmn_sp
              WRITE(*,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                   i_sp,j_sp," plcov_mx_lu  old  ",plcov_mx_lu (i_sp,j_sp,k_sp),"new ",plcovmx_sp
              WRITE(*,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                   i_sp,j_sp," lai_mn_lu    old  ",lai_mn_lu (i_sp,j_sp,k_sp),"new ",laimn_sp
              WRITE(*,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                   i_sp,j_sp," lai_mx_lu    old  ",lai_mx_lu (i_sp,j_sp,k_sp),"new ",laimx_sp
              IF (for_d_sp >= 0._wp) WRITE(*,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                   i_sp,j_sp," for_d_lu    old  ", for_d_lu(i_sp,j_sp,k_sp),"new ",for_d_sp
              IF (for_e_sp >= 0._wp)WRITE(*,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                   i_sp,j_sp," for_e_lu    old  ", for_e_lu(i_sp,j_sp,k_sp),"new ",for_e_sp
              IF (fr_land_sp >= 0._wp)WRITE(*,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                   i_sp,j_sp," fr_land    old  ", fr_land_lu(i_sp,j_sp,k_sp),"new ",fr_land_sp
              SELECT CASE (i_landuse_data)
              CASE (i_lu_globcover)
                 glc_class(1)=  'Post-flooding or irrigated croplands                                              '
                 glc_class(2)=  'Rainfed croplands                                                                 '
                 glc_class(3)=  'Mosaic Cropland (50-70%) / Vegetation (grassland, shrubland, forest) (20-50%)     '
                 glc_class(4)=  'Mosaic Vegetation (grassland, shrubland, forest) (50-70%) / Cropland (20-50%)     '
                 glc_class(5)=  'Closed to open (>15%) broadleaved evergreen and/or semi-deciduous forest (>5m)    '
                 glc_class(6)=  'Closed (>40%) --broadleaved-- deciduous forest (>5m)                              '
                 glc_class(7)=  'Open (15-40%) broadleaved deciduous forest (>5m)                                  '
                 glc_class(8)=  'Closed (>40%) --needleleaved-- evergreen forest (>5m)                             '
                 glc_class(9)=  'Open (15-40%) needleleaved deciduous or evergreen forest (>5m)                    '
                 glc_class(10)= 'Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m)             '
                 glc_class(11)= 'Mosaic Forest/Shrubland (50-70%) / Grassland (20-50%)                             '
                 glc_class(12)= 'Mosaic Grassland (50-70%) / Forest/Shrubland (20-50%)                             '
                 glc_class(13)= 'Closed to open (>15%) shrubland (<5m)                                             '
                 glc_class(14)= 'Closed to open (>15%) grassland                                                   '
                 glc_class(15)= 'Sparse (>15%) vegetation (woody vegetation, shrubs, grassland)                    '
                 glc_class(16)= 'Closed (>40%) broadleaved forest regularly flooded - Fresh water                  '
                 glc_class(17)= 'Closed (>40%) broadleaved semi-deciduous and/or evergreen forest regularly flooded'
                 glc_class(18)= 'Closed to open (>15%) vegetation  on regularly flooded or waterlogged soil        '
                 glc_class(19)= 'Artificial surfaces and associated areas (urban areas >50%)                       '
                 glc_class(20)= 'Bare areas                                                                        '
                 glc_class(21)= 'Water bodies                                                                      '
                 glc_class(22)= 'Permanent snow and ice                                                            '
                 glc_class(23)= 'undefined                                                                         '

                 DO i=1,nclass_globcover
                    WRITE (*,'(A33,1X,A85,2X,F8.4)') "Land-Use Fractions for GLOBCOVER class  ", &
                         glc_class(i),lu_class_fraction(i_sp,j_sp,k_sp,i)
                 ENDDO
              END SELECT
           ENDIF

           ! Consider only well defined variables, default is -999.!
           IF(z0_sp > 0._wp)          z0_tot (i_sp,j_sp,k_sp)      = z0_sp
           IF(rootdp_sp>0._wp)        root_lu(i_sp,j_sp,k_sp)      = rootdp_sp
           IF(soiltype_sp>0._wp )     soiltype_fao(i_sp,j_sp,k_sp) = NINT(soiltype_sp)
           IF(plcovmn_sp>0._wp)       plcov_mn_lu (i_sp,j_sp,k_sp) = plcovmn_sp
           IF(plcovmx_sp>0._wp)       plcov_mx_lu (i_sp,j_sp,k_sp) = plcovmx_sp
           IF(laimn_sp>0._wp)         lai_mn_lu (i_sp,j_sp,k_sp)   = laimn_sp
           IF(laimx_sp>0._wp)         lai_mx_lu (i_sp,j_sp,k_sp)   = laimx_sp
           IF(for_d_sp>0._wp)         for_d_lu (i_sp,j_sp,k_sp)    = for_d_sp
           IF(for_e_sp>0._wp)         for_e_lu (i_sp,j_sp,k_sp)    = for_e_sp
           IF(fr_land_sp>0._wp)    THEN
              fr_land_lu (i_sp,j_sp,k_sp)    = fr_land_sp
              fr_lake(i_sp,j_sp,k_sp) = 1. - fr_land_lu(i_sp,j_sp,k_sp)
           END IF

           ! Optional parameters for special point - not considered yet
           !                                     rs_min_lu, &
           !                                     urban_lu,  &
           !                                     for_d_lu,  &
           !                                     for_e_lu, &
           !                                     emissivity_lu, &
           !                                     soiltype_fao, &
           !                                     ndvi_max,  &
           !                                     ndvi_field_mom,&
           !                                     ndvi_ratio_mom, &
           !                                     hh_globe,            &
        END IF

     END DO ! Special Points  loop

  END IF
  !              END SELECT
  !------------------------------------------------------------------------------------------

  IF (lradtopo) THEN

     WHERE (skyview_topo .lt. 1.0 .AND. fr_land_lu .lt. 0.5)
        skyview_topo = 1.0
     ENDWHERE
     WHERE (skyview_topo .gt. 1.0)
        skyview_topo = 1.0
     ENDWHERE
  END IF

  IF (lfilter_oro) THEN
     WRITE (y_orofilter,'(A,I2,A15,I2,A13,L1,A14,I2,A16,I2,A15,I2,A13,F5.1,A15,F5.1,A12,F6.1)') &
          'Orography filter setting: lfilter_oro=.TRUE., numfilt_oro=', &
          numfilt_oro,', ilow_pass_oro=',ilow_pass_oro,', lxso_first=', &
          lxso_first,', numfilt_xso=',numfilt_xso,', ilow_pass_xso=',   &
          ilow_pass_xso,', ifill_valley=',ifill_valley,', eps_filter=', &
          eps_filter,', rfill_valley=',rfill_valley,', rxso_mask=',rxso_mask
  ELSE
     y_orofilter='lfilter_oro=.FALSE.'
  ENDIF

  PRINT *,TRIM(y_orofilter)
  !------------------------------------------------------------------------------------------
  !------------- data output ----------------------------------------------------------------
  !------------------------------------------------------------------------------------------
  fill_value_real = -1.E20_wp
  fill_value_int = -999

  CALL logging%info('write out '//TRIM(netcdf_output_filename), __FILE__, __LINE__)
  
  SELECT CASE(igrid_type)
  CASE(igrid_icon) ! ICON GRID
    CALL write_netcdf_icon_grid_extpar(TRIM(netcdf_output_filename),&
         &                                     icon_grid,                     &
         &                                     tg,                            &
         &                                     isoil_data,                    &
         &                                     ldeep_soil,                    &
         &                                     itopo_type,                    &
         &                                     lsso_param,                    &
         &                                     lscale_separation,             &
         &                                     l_use_isa,                     &
         &                                     l_use_ahf,                     &
         &                                     TRIM(y_orofilter),             &
         &                                     fill_value_real,               &
         &                                     fill_value_int,                &
         &                                     TRIM(name_lookup_table_lu),    &
         &                                     TRIM(lu_dataset),              &
         &                                     nclass_lu,                     &
         &                                     lon_geo,                       &
         &                                     lat_geo,                       &
         &                                     fr_land_lu,                    &
         &                                     lu_class_fraction,             &
         &                                     ice_lu,                        &
         &                                     z0_tot,                        &
         &                                     root_lu,                       &
         &                                     plcov_mx_lu,                   &
         &                                     lai_mx_lu,                     &
         &                                     rs_min_lu,                     &
         &                                     urban_lu,                      &
         &                                     for_d_lu,                      &
         &                                     for_e_lu,                      &
         &                                     emissivity_lu,                 &
         &                                     lake_depth,                    &
         &                                     fr_lake,                       &
         &                                     soiltype_fao,                  &
         &                                     ndvi_max,                      &
         &                                     ndvi_field_mom,                &
         &                                     ndvi_ratio_mom,                &
         &                                     hh_topo,                       &
         &                                     hh_topo_max,                   &
         &                                     hh_topo_min,                   &         
         &                                     stdh_topo,                     &
         &                                     theta_topo,                    &
         &                                     aniso_topo,                    &
         &                                     slope_topo,                    &
         &                                     vertex_param,                  &
         &                                     aot_tg,                        &
         &                                     crutemp,                       &
         &                                     alb_field_mom,                 &
         &                                     alnid_field_mom,               &
         &                                     aluvd_field_mom ,             &
         &                                     fr_sand = fr_sand,           &
         &                                     fr_silt = fr_silt,           &
         &                                     fr_clay = fr_clay,           &
         &                                     fr_oc = fr_oc,               &
         &                                     fr_bd = fr_bd,               &
         &                                     soiltype_deep=soiltype_deep, &
         &                                     fr_sand_deep=fr_sand_deep,   &
         &                                     fr_silt_deep=fr_silt_deep,   &
         &                                     fr_clay_deep=fr_clay_deep,   &
         &                                     fr_oc_deep=fr_oc_deep,       &
         &                                     fr_bd_deep=fr_bd_deep,       &
         &                                     isa_field=isa_field,         &
         &                                     ahf_field=ahf_field,         &
         &                                     sst_field=sst_field,         &
         &                                     wsnow_field=wsnow_field,     &
         &                                     t2m_field=t2m_field,         &
         &                                     hsurf_field=hsurf_field      )

  CASE(igrid_cosmo) ! COSMO grid

    IF(ldeep_soil) THEN
      CALL  write_netcdf_cosmo_grid_extpar(TRIM(netcdf_output_filename),&
           &                                     cosmo_grid,                  &
           &                                     tg,                          &
           &                                     isoil_data,                  &
           &                                     ldeep_soil,                  &
           &                                     itopo_type,                  &
           &                                     lsso_param,                  &
           &                                     l_use_isa,                   &
           &                                     l_use_ahf,                   &
           &                                     l_use_sgsl,                  &
           &                                     lscale_separation,           &
           &                                     TRIM(y_orofilter),           &
           &                                     lradtopo,                    &
           &                                     nhori,                       &
           &                                     fill_value_real,             &
           &                                     fill_value_int,              &
           &                                     TRIM(name_lookup_table_lu),  &
           &                                     TRIM(lu_dataset),            &
           &                                     i_landuse_data,              &
           &                                     nclass_lu,                   &
           &                                     lon_geo,                     &
           &                                     lat_geo,                     &
           &                                     fr_land_lu,                  &
           &                                     lu_class_fraction,           &
           &                                     ice_lu,                      &
           &                                     z0_tot,                      &
           &                                     z0_lu,                       &
           &                                     z0_topo,                     &
           &                                     z012_lu,                     &
           &                                     root_lu,                     &
           &                                     plcov_mn_lu,                 &
           &                                     plcov_mx_lu,                 &
           &                                     plcov12_lu,                  &
           &                                     lai_mn_lu,                   &
           &                                     lai_mx_lu,                   &
           &                                     lai12_lu,                    &
           &                                     rs_min_lu,                   &
           &                                     urban_lu,                    &
           &                                     for_d_lu,                    &
           &                                     for_e_lu,                    &
           &                                     emissivity_lu,               &
           &                                     lake_depth,                  &
           &                                     fr_lake,                     &
           &                                     soiltype_fao,                &
           &                                     ndvi_max,                    &
           &                                     ndvi_field_mom,              &
           &                                     ndvi_ratio_mom,              &
           &                                     hh_topo,                     &
           &                                     stdh_topo,                   &
           &                                     aot_tg,                      &
           &                                     MAC_aot_tg,                  &
           &                                     MAC_ssa_tg,                  &
           &                                     MAC_asy_tg,                  &
           &                                     crutemp,                     &
           &                                     alb_field_mom,               &
           &                                     alnid_field_mom,             &
           &                                     aluvd_field_mom,             &
           &                                     alb_dry = alb_dry,           &
           &                                     alb_sat = alb_sat,           &
           &                                     fr_sand = fr_sand,           &
           &                                     fr_silt = fr_silt,           &
           &                                     fr_clay = fr_clay,           &
           &                                     fr_oc = fr_oc,               &
           &                                     fr_bd = fr_bd,               &
           &                                     soiltype_deep=soiltype_deep, &
           &                                     fr_sand_deep=fr_sand_deep,   &
           &                                     fr_silt_deep=fr_silt_deep,   &
           &                                     fr_clay_deep=fr_clay_deep,   &
           &                                     fr_oc_deep=fr_oc_deep,       &
           &                                     fr_bd_deep=fr_bd_deep,       &
           &                                     theta_topo=theta_topo,       &
           &                                     aniso_topo=aniso_topo,       &
           &                                     slope_topo=slope_topo,       &
           &                                     slope_asp_topo=slope_asp_topo,&
           &                                     slope_ang_topo=slope_ang_topo,&
           &                                     horizon_topo=horizon_topo,    &
           &                                     skyview_topo=skyview_topo,    &
           &                                     isa_field=isa_field,          &
           &                                     ahf_field=ahf_field,          &
           &                                     sgsl = sgsl                   )

    ELSE
      CALL  write_netcdf_cosmo_grid_extpar(TRIM(netcdf_output_filename),&
           &                                     cosmo_grid,                  &
           &                                     tg,                          &
           &                                     isoil_data,                  &
           &                                     ldeep_soil,                  &
           &                                     itopo_type,                  &
           &                                     lsso_param,                  &
           &                                     l_use_isa,                   &
           &                                     l_use_ahf,                   &
           &                                     l_use_sgsl,                  &
           &                                     lscale_separation,           &
           &                                     y_orofilter,                 &
           &                                     lradtopo,                    &
           &                                     nhori,                       &
           &                                     fill_value_real,             &
           &                                     fill_value_int,              &
           &                                     TRIM(name_lookup_table_lu),  &
           &                                     TRIM(lu_dataset),            &
           &                                     i_landuse_data,              &
           &                                     nclass_lu,                   &
           &                                     lon_geo,                     &
           &                                     lat_geo,                     &
           &                                     fr_land_lu,                  &
           &                                     lu_class_fraction,           &
           &                                     ice_lu,                      &
           &                                     z0_tot,                      &
           &                                     z0_lu,                       &
           &                                     z0_topo,                     &
           &                                     z012_lu,                     &
           &                                     root_lu,                     &
           &                                     plcov_mn_lu,                 &
           &                                     plcov_mx_lu,                 &
           &                                     plcov12_lu,                  &
           &                                     lai_mn_lu,                   &
           &                                     lai_mx_lu,                   &
           &                                     lai12_lu,                    &
           &                                     rs_min_lu,                   &
           &                                     urban_lu,                    &
           &                                     for_d_lu,                    &
           &                                     for_e_lu,                    &
           &                                     emissivity_lu,               &
           &                                     lake_depth,                  &
           &                                     fr_lake,                     &
           &                                     soiltype_fao,                &
           &                                     ndvi_max,                    &
           &                                     ndvi_field_mom,              &
           &                                     ndvi_ratio_mom,              &
           &                                     hh_topo,                     &
           &                                     stdh_topo,                   &
           &                                     aot_tg,                      &
           &                                     MAC_aot_tg,                  &
           &                                     MAC_ssa_tg,                  &
           &                                     MAC_asy_tg,                  &
           &                                     crutemp,                     &
           &                                     alb_field_mom,               &
           &                                     alnid_field_mom,             &
           &                                     aluvd_field_mom,             &
           &                                     alb_dry = alb_dry,           &
           &                                     alb_sat = alb_sat,           &
           &                                     fr_sand = fr_sand,           &
           &                                     fr_silt = fr_silt,           &
           &                                     fr_clay = fr_clay,           &
           &                                     fr_oc = fr_oc,               &
           &                                     fr_bd = fr_bd,               &
           &                                     theta_topo=theta_topo,       &
           &                                     aniso_topo=aniso_topo,       &
           &                                     slope_topo=slope_topo,       &
           &                                     slope_asp_topo=slope_asp_topo,&
           &                                     slope_ang_topo=slope_ang_topo,&
           &                                     horizon_topo=horizon_topo,    &
           &                                     skyview_topo=skyview_topo,   &
           &                                     isa_field=isa_field,         &
           &                                     ahf_field=ahf_field,          &
           &                                     sgsl = sgsl                   )
    ENDIF
  END SELECT

END PROGRAM
