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

  USE mo_logging
  USE info_extpar,              ONLY: info_print
  USE mo_kind,                  ONLY: wp, i4
  USE mo_io_units,              ONLY: filename_max

  USE mo_target_grid_data,      ONLY: lon_geo, lat_geo, tg

  USE mo_grid_structures,       ONLY: igrid_icon, igrid_cosmo

  USE mo_cosmo_grid,            ONLY: COSMO_grid

  USE mo_icon_grid_data,        ONLY: icon_grid, &
       &                              icon_grid_region

  USE mo_read_extpar_namelists, ONLY: read_namelists_extpar_check_icon,    &
       &                              read_namelists_extpar_check_cosmo,   &
       &                              read_namelists_extpar_special_points

  USE mo_soil_routines,         ONLY: read_namelists_extpar_soil

  USE mo_soil_data,             ONLY: define_soiltype, &
       &                              undef_soiltype, default_soiltype,&
       &                              soiltype_ice, soiltype_water, &
       &                              soil_data, FAO_data, HWSD_data, HWSD_map

  USE mo_soil_tg_fields,        ONLY: fr_sand,fr_silt,fr_clay, &
       &                              fr_oc, fr_bd, &
       &                              fr_sand_deep,fr_silt_deep, &
       &                              fr_clay_deep, fr_oc_deep,  &
       &                              fr_bd_deep, &
       &                              fr_land_soil, &
       &                              soiltype_fao, soiltype_hwsd, soiltype_deep,soiltype_hwsd_s, &
       &                              allocate_soil_target_fields

  USE mo_soil_consistency,      ONLY: calculate_soiltype
                                
  USE mo_soil_output_nc,        ONLY: read_netcdf_soil_buffer
                                
  USE mo_target_grid_routines,  ONLY: init_target_grid
                                
  USE mo_glcc_lookup_tables,    ONLY: nclass_glcc

  USE mo_glcc_tg_fields,        ONLY: fr_land_glcc, &
       &                              glcc_class_fraction,    &
       &                              glcc_class_npixel, &
       &                              glcc_tot_npixel, &
       &                              ice_glcc, &
       &                              z0_glcc, &
       &                              root_glcc, &
       &                              plcov_mn_glcc, &
       &                              plcov_mx_glcc, &
       &                              lai_mn_glcc, &
       &                              lai_mx_glcc, &
       &                              rs_min_glcc, &
       &                              urban_glcc,  &
       &                              for_d_glcc,  &
       &                              for_e_glcc, &
       &                              emissivity_glcc, &
       &                              allocate_glcc_target_fields

  USE mo_glc2000_lookup_tables, ONLY: get_name_glc2000_lookup_tables, &
       &                              nclass_glc2000

  USE mo_globcover_data,        ONLY: max_tiles_lu

  USE mo_globcover_lookup_tables, ONLY: nclass_globcover, &
       &                                get_name_globcover_lookup_tables

  USE mo_ecci_lookup_tables,      ONLY: nclass_ecci, &
      &                                 get_name_ecci_lookup_tables

  USE mo_ecoclimap_lookup_tables, ONLY: nclass_ecoclimap, &
       &                                get_name_ecoclimap_lookup_tables

  USE mo_glcc_lookup_tables,    ONLY: get_name_glcc_lookup_tables


  USE mo_landuse_output_nc,     ONLY: read_netcdf_buffer_glcc, &
       &                              read_netcdf_buffer_lu, &
       &                              read_netcdf_buffer_ecoclimap

  USE mo_isa_output_nc,         ONLY: read_netcdf_buffer_isa

  USE mo_landuse_routines,      ONLY: read_namelists_extpar_land_use

  USE mo_lu_tg_fields,          ONLY: i_lu_globcover, i_lu_glc2000, i_lu_glcc, &
       &                              i_lu_ecoclimap, i_lu_ecci, & 
       &                              fr_land_lu, &
       &                              fr_land_mask, &
       &                              ice_lu, &
       &                              z0_lu, &
       &                              z0_tot, &
       &                              root_lu, &
       &                              plcov_mn_lu, &
       &                              plcov_mx_lu, &
       &                              lai_mn_lu, &
       &                              lai_mx_lu, &
       &                              rs_min_lu, &
       &                              urban_lu,  &
       &                              for_d_lu,  &
       &                              for_e_lu, &
       &                              skinc_lu, &
       &                              emissivity_lu, &
       &                              fr_ocean_lu, &
       &                              lu_class_fraction,    &
       &                              lu_class_npixel, &
       &                              lu_tot_npixel,  &
       &                              lai12_lu,     &
       &                              z012_lu,      &
       &                              plcov12_lu, &
       &                              allocate_lu_target_fields, allocate_add_lu_fields

  USE mo_isa_tg_fields,         ONLY: isa_field, &
       &                              isa_tot_npixel, &
       &                              allocate_isa_target_fields, &
       &                              allocate_add_isa_fields

  USE mo_isa_routines,          ONLY: read_namelists_extpar_isa

  USE mo_ahf_tg_fields,         ONLY: ahf_field, &
       &                              allocate_ahf_target_fields

  USE mo_ahf_data,              ONLY: undef_ahf, minimal_ahf, iahf_type

  USE mo_ahf_output_nc,         ONLY: read_netcdf_buffer_ahf

  USE mo_ahf_routines,          ONLY: read_namelists_extpar_ahf


  USE mo_era_output_nc,         ONLY: read_netcdf_buffer_sst,&
       &                              read_netcdf_buffer_t2m

  USE mo_topo_tg_fields,        ONLY: fr_land_topo,       &
       &                              hh_topo,            &
       &                              hh_topo_max,        &
       &                              hh_topo_min,        &       
       &                              stdh_topo,          &
       &                              theta_topo,         &
       &                              aniso_topo,         &
       &                              slope_topo,         &
       &                              z0_topo,            &
       &                              slope_asp_topo,     &
       &                              slope_ang_topo,     &
       &                              horizon_topo,       &
       &                              skyview_topo,       &
       &                              sgsl,               &
       &                              vertex_param,       &
       &                              allocate_additional_param, &
       &                              allocate_topo_target_fields

  USE mo_topo_output_nc,        ONLY: read_netcdf_buffer_topo

  USE mo_topo_routines,         ONLY: read_namelists_extpar_orography, &
       &                              read_namelists_extpar_scale_sep

  USE mo_topo_data,             ONLY: lradtopo, nhori, max_tiles, itopo_type, &
       &                              radius, min_circ_cov, max_missing, itype_scaling

  USE mo_aot_target_fields,     ONLY: allocate_aot_target_fields,&
       &                              aot_tg,&
       &                              MAC_aot_tg,&
       &                              MAC_ssa_tg,&
       &                              MAC_asy_tg

  USE mo_aot_output_nc,         ONLY: read_netcdf_buffer_aot, read_netcdf_buffer_aot_MAC

  USE mo_aot_data,              ONLY: ntype_aot, ntime_aot, iaot_type, n_spectr, nspb_aot

  USE mo_aot_data,              ONLY: read_namelists_extpar_aerosol

  USE mo_flake_routines,        ONLY: read_namelists_extpar_flake

  USE mo_flake_tg_fields,       ONLY: fr_lake, &
       &                              lake_depth,    &
       &                              flake_tot_npixel, &
       &                              allocate_flake_target_fields

  USE mo_flake_data,            ONLY: flake_depth_undef, & !< default value for undefined lake depth
       &                              flake_depth_default, & !< default value for default lake depth, 10 [m]
       &                              DWD_max_lake_depth, & !< Maximum lake depth in [m] for FLAKE (50 m)
       &                              DWD_min_lake_depth !< Minimal lake depth in [m] for FLAKE (1 m)

  USE mo_flake_output_nc,       ONLY: read_netcdf_buffer_flake
                                
  USE mo_lsm_output_nc,         ONLY: read_netcdf_buffer_lsm
                                
  USE mo_extpar_output_nc,      ONLY: write_netcdf_cosmo_grid_extpar, &
       &                              write_cdi_icon_grid_extpar

  USE mo_lradtopo,              ONLY: read_namelists_extpar_lradtopo

  USE mo_search_target_grid,    ONLY: find_nearest_target_grid_element

  USE mo_oro_filter,            ONLY: read_namelists_extpar_orosmooth

  USE mo_isa_data,              ONLY: max_tiles_isa, &
       &                              undef_isa, minimal_isa, isa_type

  USE mo_python_data,           ONLY: &
  ! emiss                                      
       &                              ntime_emiss, &
       &                              minimal_emiss, &
  ! ndvi                                      
       &                              ntime_ndvi, &
       &                              undef_ndvi, minimal_ndvi, &
  ! albedo                                      
       &                              ntime_alb, &
       &                              wso_min,wso_max,csalb,csalbw,zalso, &
       &                              allocate_alb_interp_fields, &
       &                              alb_interp_data, &
       &                              minimal_alb_dry, maximal_alb_dry, &
       &                              ntime_alb, &
       &                              minimal_alb_sat, maximal_alb_sat, undef_alb_bs, &
       &                              ialb_type, &
  ! era
       &                              ntime_era

  USE mo_python_routines,       ONLY: read_namelists_extpar_emiss, &
       &                              read_namelists_extpar_t_clim, &
       &                              read_namelists_extpar_ndvi, &
       &                              read_namelists_extpar_alb, &
       &                              open_netcdf_ALB_data, &
       &                              const_check_interpol_alb, &
       &                              read_namelists_extpar_era

  USE mo_python_tg_fields,      ONLY: &
  ! emiss                                      
       &                              emiss_max, &
       &                              emiss_field_mom, &
       &                              emiss_ratio_mom, &
       &                              allocate_emiss_target_fields, &
  ! ndvi                                      
       &                              ndvi_max, &
       &                              ndvi_field_mom, &
       &                              ndvi_ratio_mom, &
       &                              allocate_ndvi_target_fields, &
       &                              allocate_cru_target_fields,   &
  ! cru                                      
       &                              ndvi_max, &
       &                              crutemp, crutemp2, cruelev, &
  ! albedo                                      
       &                              alb_dry, alb_sat, &
       &                              alb_field_mom, &
       &                              alnid_field_mom, &
       &                              aluvd_field_mom, &
       &                              allocate_alb_target_fields, &
  ! era
       &                              sst_field, &
       &                              wsnow_field, &
       &                              t2m_field, &
       &                              hsurf_field, &
       &                              allocate_era_target_fields

  USE mo_python_output_nc,      ONLY: read_netcdf_buffer_emiss, &
       &                              read_netcdf_buffer_ndvi, &
       &                              read_netcdf_buffer_cru, &
       &                              read_netcdf_buffer_alb, &
       &                              read_netcdf_buffer_era

  IMPLICIT NONE

  CHARACTER (len=filename_max)                  :: namelist_grid_def, &
       &                                           netcdf_output_filename, &
       &                                           grib_output_filename, &
       &                                           grib_sample, &
       &                                           namelist_file, & !< filename with namelists for for EXTPAR settings
  ! soil                                        
       &                                           soil_buffer_file, &  !< name for soil buffer file
       &                                           soil_buffer_file_consistent, & !< name for soil buffer file after consistency check
       &                                           soil_output_file_consistent, & !< name for soil output file after consistency check
       &                                           raw_data_soil_path, &        !< path to raw data
       &                                           raw_data_soil_filename, & !< filename soil raw data
       &                                           raw_data_deep_soil_filename, & !< filename deep soil raw data
  ! orography                                   
       &                                           orography_output_file,  &
       &                                           orography_buffer_file, & !< name for orography buffer file
       &                                           raw_data_orography_path, &        !< path to raw data
  ! subgrid-scale slope                         
       &                                           sgsl_files(1:max_tiles), &  !< filenames globe raw data
       &                                           sgsl_buffer_file, & !< name for orography buffer file
  ! land use                                    
       &                                           raw_data_lu_path, &        !< path to raw data
       &                                           raw_data_lu_filename(1:max_tiles_lu), & !< filename glc2000 raw data !_br 21.02.14
       &                                           name_lookup_table_lu, & !< name for look up table
       &                                           lu_dataset, & !< name of landuse data set
       &                                           lu_buffer_file, & !< name for glc2000 buffer file
       &                                           raw_data_path, &
       &                                           glcc_buffer_file, &    !< name for glcc buffer file
  ! albedo                                      
       &                                           raw_data_alb_path, &   !< path to albedo raw input data
       &                                           raw_data_alb_filename, &    !< raw data filename
       &                                           raw_data_alnid_filename, &  !< raw data filename, NI
       &                                           raw_data_aluvd_filename, &  !< raw data filename, UV
       &                                           alb_buffer_file, &    !< name for albedo buffer file
       &                                           alb_output_file, &    !< name for albedo output file
       &                                           land_sea_mask_file, & !< name for external land-sea-mask
  ! ISA                                         
       &                                           raw_data_isa_path, &        !< path to raw data
       &                                           raw_data_isa_filename(1:max_tiles_isa), & !< filename glc2000 raw data
       &                                           isa_buffer_file, & !< name for NDVI buffer file
  !AHF                                          
       &                                           raw_data_ahf_path, &        !< path to raw data
       &                                           raw_data_ahf_filename, & !< filename NDVI raw data
       &                                           ahf_buffer_file, & !< name for NDVI buffer file
  ! NDVI                                        
       &                                           raw_data_ndvi_path, &
       &                                           raw_data_ndvi_filename, &
       &                                           ndvi_buffer_file, & !< name for NDVI buffer file
       &                                           ndvi_output_file, &
 ! EMISS                                        
       &                                           emiss_buffer_file, & !< name for EMISS buffer file
       &                                           raw_data_emiss_path, & !< dummy var for routine read_namelist_extpar_emiss
       &                                           raw_data_emiss_filename, &!< dummy var for routine read_namelist_extpar_emiss
       &                                           emiss_output_file, &!< dummy var for routine read_namelist_extpar_emiss
  ! ERA
       &                                           era_buffer_file, & !< name for SST buffer-file
       &                                           sst_file_legacy, & !< name of SST-file (legacy)
       &                                           t2m_file_legacy, & !< name of T2M-file (legacy)
  ! temperature climatology                     
       &                                           namelist_file_t_clim, & !< filename with namelists for for EXTPAR settings
       &                                           raw_data_t_clim_path, &     !< path to raw data
       &                                           raw_data_t_clim_filename, & !< filename temperature climatology raw data
       &                                           t_clim_buffer_file, & !< name for temperature climatology buffer
       &                                           t_clim_output_file, & !< name for temperature climatology output file
  ! aerosol optical thickness                   
       &                                           raw_data_aot_path, &        !< path to raw data
       &                                           raw_data_aot_filename, & !< filename temperature climatology raw data
       &                                           aot_buffer_file, & !< name for aerosol buffer file
       &                                           topo_files(1:max_tiles), & !< filenames globe raw data
  ! flake
       &                                           raw_data_flake_path, &
       &                                           raw_data_flake_filename, &
       &                                           flake_buffer_file, &  !< name for flake buffer file
  !special points                               
       &                                           path_alb_file, &
       &                                           namelist_alb_data_input, &
  ! Namelist values for topography scale separation
       &                                           raw_data_scale_sep_orography_path, & !< path to raw data
       &                                           scale_sep_files(1:max_tiles)  !< filenames globe raw data

  CHARACTER (LEN=filename_max), DIMENSION(1:23) :: glc_class
  CHARACTER (LEN=filename_max), DIMENSION(1:38) :: ecci_class
  ! define string used for global attributes:
  CHARACTER (LEN=255)                           :: y_orofilter

  INTEGER (KIND=i4)                             :: fill_value_int, &   !< value for undefined integer
       &                                           ntiles_column, &
       &                                           ntiles_row, &
       &                                           it_cl_type, &
       &                                           isoil_data, &
       &                                           ntiles_isa, &
       &                                           i, j ,k, t, & !< counter
       &                                           jj, ii, & !< counter
       &                                           n, &
       &                                           nloops, &
       &                                           nnb, & !< number of neighbor grid elements with common edge
       &                                           nv, &  !< number of vertices
       &                                           n_index, & !< help variable
       &                                           ne_ie(9), & !< index for grid element neighbor
       &                                           ne_je(9), & !< index for grid element neighbor
       &                                           ne_ke(9), & !< index for grid element neighbor
       &                                           igrid_type, &  !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid
       &                                           i_landuse_data, & !<integer switch to choose a land use raw data set
       &                                           i_lsm_data, & !<integer switch to choose a land sea mask data set
       &                                           ilookup_table_lu, & !< integer switch to choose a lookup table
       &                                           nclass_lu, & !< number of land use classes
       &                                           count_ice2tclim,count_ice2tclim_tile, &
       &                                           start_cell_id, & !< ID of starting cell for ICON search
       &                                           isp,i_sp,j_sp,k_sp, &
       &                                           ncid_alb, &
       &                                           tile_mode, &
       &                                           db_ice_counter,db_water_counter, &
       &                                           number_special_points, &
       &                                           ilow_pass_oro,   &
       &                                           numfilt_oro,     &
       &                                           ifill_valley,    &
       &                                           ilow_pass_xso,   &
       &                                           numfilt_xso, &
  ! T_CL consistency check
       &                                           iml, imu, ipl, ipu, jml, jmu, jpl, jpu, l, &
       &                                           ntclct

  INTEGER (KIND=i4), PARAMETER                  :: mpy=12, &     !< month per year
       &                                           i_gcv__snow_ice = 22, & ! GlobCover land-use class for glaciers
       &                                           i_gcv_bare_soil = 20, & ! GlobCover land-use class for bare-soil
       &                                           i_ecci__snow_ice = 38, & ! ESA CCI LU CLASS for glaciers
       &                                           i_ecci_bare_soil = 34    ! ESA CCI LU CLASS for bare-soil


  LOGICAL                                       :: last=.FALSE., & ! in TCL leave loop
       &                                           foundtcl=.FALSE., & ! in TCL
       &                                           lsso_param,lsubtract_mean_slope, &
       &                                           ldeep_soil, &
       &                                           l_use_isa =.FALSE., & !< flag if additional urban data are present
       &                                           l_use_ahf =.FALSE., & !< flag if additional urban data are present
       &                                           l_use_sgsl=.FALSE., & !< flag if sgsl is used in topo
       &                                           l_preproc_oro=.FALSE., & 
       &                                           l_use_glcc=.FALSE., & !< flag if additional glcc data are present
       &                                           l_use_emiss=.FALSE., &!< flag if additional CAMEL emissivity data are present
       &                                           l_unified_era_buffer=.FALSE., &!< flag if ERA-data from extpar_era_to_buffer.py is used
       &                                           lwrite_netcdf, &  !< flag to enable netcdf output for COSMO
       &                                           lwrite_grib, &    !< flag to enable GRIB output for COSMO
       &                                           lflake_correction, & !< flag to correct fr_lake and depth_lake near coastlines
       &                                           tile_mask, &
  ! Namelist values for topography scale separation
       &                                           lscale_separation, &
  ! Namelist values for orography smoothing
       &                                           lfilter_oro,     &
       &                                           lscale_file=.FALSE., &
       &                                           lxso_first, & 
       &                                           l_use_corine


  REAL (KIND=wp)                                :: t2mclim_hc, &
       &                                           step, &
       &                                           thr_cr, & !< control threshold
       &                                           fill_value_real, & !< value to indicate undefined grid elements
  ! for albedo consistency check
       &                                           albvis_min, albnir_min, albuv_min, &
       &                                           hh_dead_sea, &
       &                                           lu_data_southern_boundary, &
  !Special Points
       &                                           lon_geo_sp=-999.,           &
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
                                                   fr_land_sp=-999., &
  ! T_CL consistency check
       &                                           tclsum, elesum, &
       &                                           eps_filter,      &
       &                                           rfill_valley,    &
       &                                           rxso_mask

  REAL(KIND=wp), PARAMETER                      :: dtdz_clim = -5.e-3_wp, & !-5 K/km indicate undefined land use grid elements
       &                                           tmelt = 273.15_wp, &
       &                                           frlndtile_thrhld=0.05_wp, &
       &                                           undefined_lu = 0.0_wp !< value to indicate undefined land use grid elements

  LOGICAL :: l_use_array_cache
  
  !---------------------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------------------

  CALL initialize_logging("extpar_consistency.log")
  CALL info_print ()

  !--------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= start consistency_check ==========')
  CALL logging%info( '')

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= read namelist ====================')
  CALL logging%info( '')

  !--------------------------------------------------------------------------------------------------------
  ! Get flake buffer file name from namelist
  namelist_file = 'INPUT_FLAKE'
  CALL read_namelists_extpar_flake(namelist_file, &
       &                           raw_data_flake_path, &
       &                           raw_data_flake_filename, &
       &                           flake_buffer_file)

  ! Get ndvi buffer file name from namelist
  namelist_file = 'INPUT_NDVI'

  CALL read_namelists_extpar_ndvi(namelist_file, &
       &                          raw_data_ndvi_path, &
       &                          raw_data_ndvi_filename, &
       &                          ndvi_buffer_file, &
       &                          ndvi_output_file)
  

  ! Get lradtopo and nhori value from namelist

  namelist_file = 'INPUT_RADTOPO'
  CALL read_namelists_extpar_lradtopo(namelist_file,lradtopo,nhori, radius,min_circ_cov,max_missing, itype_scaling)

  ! Get lsso_param from namelist

  namelist_file = 'INPUT_ORO'
  CALL read_namelists_extpar_orography(namelist_file,          &
       &                               raw_data_orography_path,&
       &                               topo_files,             &
       &                               sgsl_files,             &
       &                               ntiles_column,          &
       &                               ntiles_row,             &
       &                               itopo_type,             &
       &                               l_use_sgsl,             &
       &                               l_preproc_oro,          &
       &                               lsso_param,             &
       &                               lsubtract_mean_slope,   &
       &                               orography_buffer_file,  &
       &                               orography_output_file,  &
       &                               sgsl_buffer_file)

  IF (l_use_sgsl) THEN
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    CALL logging%info( '')
    CALL logging%warning( 'Subgrid-slope (SGSL) active')
    CALL logging%info( '')
  ENDIF

  namelist_file = 'INPUT_SOIL'
  CALL read_namelists_extpar_soil(namelist_file,                     &
       isoil_data,                 &
       ldeep_soil,                 &
       raw_data_soil_path,         &
       raw_data_soil_filename,     &
       raw_data_deep_soil_filename,&
       soil_buffer_file,           &
       soil_buffer_file_consistent,&
       soil_output_file_consistent)

  IF (ldeep_soil .AND. isoil_data /= HWSD_data) THEN
     ldeep_soil = .FALSE.
     CALL logging%warning('One can use deep soil only, if HWSD data is used => ldeep_soil is set to .FALSE.')
   ENDIF

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
       &                                  alb_output_file)

  !--------------------------------------------------------------
  ! get namelist for aerosol fields
  !--------------------------------------------------------------
  namelist_file = 'INPUT_AOT'
  CALL read_namelists_extpar_aerosol(namelist_file, &
       &                                  iaot_type,    &
       &                                  raw_data_aot_path, &
       &                                  raw_data_aot_filename, &
       &                                  aot_buffer_file)

  !--------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------
  ! get namelist for topography scale separation
  !--------------------------------------------------------------
  namelist_file = 'INPUT_SCALE_SEP'
  INQUIRE(file=TRIM(namelist_file),exist=lscale_file)
  IF (lscale_file) THEN
    CALL logging%info('Scale separation active')
    CALL read_namelists_extpar_scale_sep(namelist_file,                     &
         &                                  raw_data_scale_sep_orography_path, &
         &                                  scale_sep_files,                   &
         &                                  lscale_separation)
  ENDIF

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
  igrid_type = tg%igrid_type

  SELECT CASE(igrid_type)
    CASE(igrid_icon) ! ICON GRID
      CALL  allocate_additional_param(icon_grid%nvertex,l_use_sgsl, l_use_array_cache)
  END SELECT

  ! get info on raw data file
  namelist_file = 'INPUT_LU'
  CALL read_namelists_extpar_land_use(namelist_file, &
       &                                 i_landuse_data, &
       &                                 l_use_corine,   &
       &                                 raw_data_lu_path, &
       &                                 raw_data_lu_filename, &
       &                                 ilookup_table_lu, &
       &                                 lu_buffer_file, &
       &                                 raw_data_glcc_path_opt = raw_data_path, &
       &                                 glcc_buffer_file_opt = glcc_buffer_file)

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
       lu_data_southern_boundary = -56.0 ! Needed to capture the Antarctic peninsula
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
    CASE (i_lu_ecci)
     lu_dataset = 'ESA CCI'
     CALL get_name_ecci_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
     nclass_lu = nclass_ecci
     lu_data_southern_boundary = -90.0 ! No need to capture the Antarctic peninsula
  END SELECT

  WRITE(message_text,*)'Land use datatset    : '//TRIM(lu_dataset)
  CALL logging%info(message_text)
  WRITE(message_text,*)'Land use lookup table: '//TRIM(name_lookup_table_lu)
  CALL logging%info(message_text)

  !-----------------------------------------------------------------------------------------------
  ! get info on urban data file
  !-----------------------------------------------------------------------------------------------
  namelist_file = 'INPUT_ISA'
  INQUIRE(file=TRIM(namelist_file),exist=l_use_isa)
  IF (l_use_isa) THEN
     CALL logging%info('Urban data ISA active')
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
     CALL logging%info('Urban data AHF active')
     CALL  read_namelists_extpar_ahf(namelist_file, &
          &                                  iahf_type,    & !_br 14.04.16
          &                                  raw_data_ahf_path, &
          &                                  raw_data_ahf_filename, &
          &                                  ahf_buffer_file)
  END IF


  !-----------------------------------------------------------------------------------------------
  ! get filenames from namelist
  !-----------------------------------------------------------------------------------------------
  namelist_file = 'INPUT_CHECK'

  SELECT CASE(igrid_type)
  CASE(igrid_icon)
    CALL read_namelists_extpar_check_icon(namelist_file, &
         grib_output_filename, &
         grib_sample, &
         netcdf_output_filename, &
         sst_file_legacy, &
         t2m_file_legacy, &
         i_lsm_data, &
         land_sea_mask_file,&
         lwrite_netcdf,         &
         lwrite_grib,           &
         number_special_points, &
         tile_mode,             &
         l_use_glcc,            &
         l_use_array_cache)

    INQUIRE(file=TRIM(glcc_buffer_file),exist=l_use_glcc)

  CASE(igrid_cosmo)
    CALL read_namelists_extpar_check_cosmo( &
         namelist_file, &
         grib_output_filename, &
         grib_sample, &
         netcdf_output_filename, &
         i_lsm_data, &
         land_sea_mask_file,&
         lwrite_netcdf,         &
         lwrite_grib,           &
         number_special_points, &
         tile_mode,             &
         lflake_correction,     &
         l_use_array_cache)

    INQUIRE(file=TRIM(glcc_buffer_file),exist=l_use_glcc)
  END SELECT

  namelist_file_t_clim = 'INPUT_TCLIM'
  CALL read_namelists_extpar_t_clim(namelist_file_t_clim,     &
       it_cl_type,               &
       raw_data_t_clim_path,     &
       raw_data_t_clim_filename, &
       t_clim_buffer_file      , &
       t_clim_output_file        )

  ! read namelist for input EMISS data
  namelist_file = 'INPUT_EMISS'
  INQUIRE(FILE=TRIM(namelist_file), EXIST=l_use_emiss)
  IF (l_use_emiss) THEN
    CALL  read_namelists_extpar_emiss(namelist_file, &
      &                               raw_data_emiss_path, &
      &                               raw_data_emiss_filename, &
      &                               emiss_buffer_file, &
      &                               emiss_output_file)
  ENDIF

  ! determine type of ERA-climatologies to use
  IF (igrid_type == igrid_icon) THEN

    ! read namelist INPUT_ERA
    namelist_file = 'INPUT_ERA'
    INQUIRE(FILE=TRIM(namelist_file), EXIST=l_unified_era_buffer)
    IF (l_unified_era_buffer) THEN
      CALL logging%info('Use ERA-climatologies provided by extpar_era_to_buffer.py')
      CALL  read_namelists_extpar_era(namelist_file, era_buffer_file)
    ELSE
      CALL logging%info('Use legacy ERA-climatologies generated by ICON-REMAP')
    ENDIF
  ENDIF

  IF (tile_mode == 1) THEN
     tile_mask=.TRUE.
     WRITE(message_text,*)'Tile mode for EXTPAR is set to tile_mode= ',tile_mode,'tile_mask= ',tile_mask
     CALL logging%info(message_text)
  ELSE
    tile_mask=.FALSE.
     WRITE(message_text,*)'Tile mode for EXTPAR is set to tile_mode= ',tile_mode,'tile_mask= ',tile_mask
     CALL logging%info(message_text)
  END IF

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info( '============= allocate fields ==================')
  CALL logging%info( '')

  ! test for glcc data
  IF (l_use_glcc) THEN
    CALL allocate_glcc_target_fields(tg, l_use_array_cache)
  ENDIF
  ! allocate Land use target fields
  CALL allocate_lu_target_fields(tg, l_use_array_cache)
  CALL allocate_add_lu_fields(tg,nclass_lu, l_use_array_cache)

  CALL allocate_soil_target_fields(tg,ldeep_soil, l_use_array_cache)

  IF (l_use_ahf) THEN
    CALL allocate_ahf_target_fields(tg, l_use_array_cache)
  END IF
  IF (l_use_isa) THEN
    CALL allocate_isa_target_fields(tg, l_use_array_cache)
    CALL allocate_add_isa_fields(tg, l_use_array_cache)
  END IF

  CALL allocate_ndvi_target_fields(tg,ntime_ndvi, l_use_array_cache)

  CALL allocate_emiss_target_fields(tg,ntime_emiss, l_use_array_cache)

  CALL allocate_era_target_fields(tg,ntime_era, l_use_array_cache) ! sst clim contains also 12 monthly values as ndvi

  IF (lscale_separation .AND. (itopo_type == 2)) THEN
    lscale_separation = .FALSE.
    CALL logging%warning('Scale separation can only be used with GLOBE topography => lscale_separation set to .FALSE.')
  ENDIF

  CALL allocate_topo_target_fields(tg,nhori,l_use_sgsl, l_use_array_cache)

  CALL allocate_aot_target_fields(tg, iaot_type, ntime_aot, ntype_aot, nspb_aot, l_use_array_cache)

  CALL allocate_cru_target_fields(tg, l_use_array_cache)

  CALL allocate_flake_target_fields(tg, l_use_array_cache)

  CALL allocate_alb_target_fields(tg,ntime_alb,ialb_type, l_use_array_cache)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info( '============= read input ======================')
  CALL logging%info( '')
 
  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Landuse')

  SELECT CASE (i_landuse_data)
    CASE (i_lu_ecoclimap)
       lu_dataset = 'ECOCLIMAP'
       CALL read_netcdf_buffer_ecoclimap(lu_buffer_file,  &
            &                                     tg,         &
            &                                     nclass_lu, &
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


    CASE(i_lu_globcover, i_lu_glc2000, i_lu_ecci)
       CALL read_netcdf_buffer_lu(lu_buffer_file,  &
            &                                     tg,         &
            &                                     nclass_lu, &
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

       IF (l_use_glcc) THEN
          CALL read_netcdf_buffer_glcc(glcc_buffer_file,  &
               &                                     tg,         &
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
  END SELECT ! GlobCover needs also GLCC!

  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Soil')

  IF(ldeep_soil) THEN
     CALL read_netcdf_soil_buffer(soil_buffer_file,    &
          &                       tg,          &
          &                       isoil_data,  &
          &                       fr_land_soil,&
          &                       soiltype_fao,&
          &                       soiltype_hwsd,&
          &                       soiltype_deep,&
          &                       soiltype_hwsd_s)
  ELSE
     SELECT CASE(isoil_data)
       CASE(FAO_data, HWSD_map)
          CALL read_netcdf_soil_buffer(soil_buffer_file,    &
               &                       tg,          &
               &                       isoil_data,  &
               &                       fr_land_soil,&
               &                       soiltype_fao,&
               &                       soiltype_hwsd)
       CASE(HWSD_data)
          CALL read_netcdf_soil_buffer(soil_buffer_file,    &
               &                       tg,          &
               &                       isoil_data,  &
               &                       fr_land_soil,&
               &                       soiltype_fao,&
               &                       soiltype_hwsd )
     END SELECT
  ENDIF

  !-------------------------------------------------------------------------
  IF (l_use_isa) THEN
    CALL logging%info( '')
    CALL logging%info('ISA')
    CALL read_netcdf_buffer_isa(isa_buffer_file,  &
          &                      tg,         &
          &                      isa_field, &
          &                      isa_tot_npixel )
  END IF

  !-------------------------------------------------------------------------
  IF (l_use_ahf) THEN
    CALL logging%info( '')
    CALL logging%info('AHF')
    CALL read_netcdf_buffer_ahf(ahf_buffer_file,  &
          &                     tg,         &
          &                     ahf_field )
  END IF

  !-------------------------------------------------------------------------
  IF(igrid_type == igrid_icon) THEN
    CALL logging%info( '')
    CALL logging%info('ERA')
    CALL read_netcdf_buffer_era(era_buffer_file, &
          &                     sst_file_legacy, &
          &                     t2m_file_legacy, &
          &                     tg,         &
          &                     ntime_era, &
          &                     l_unified_era_buffer, &
          &                     sst_field,&
          &                     wsnow_field,&
          &                     t2m_field,&
          &                     hsurf_field)
  END IF

  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Albedo')
  IF (ialb_type == 2) THEN
     CALL read_netcdf_buffer_alb(alb_buffer_file,  &
          &                           tg, &
          &                           ntime_alb, &
          &                           alb_dry=alb_dry, &
          &                           alb_sat=alb_sat)
  ELSE IF (ialb_type == 1) THEN
     CALL read_netcdf_buffer_alb(alb_buffer_file,  &
          &                           tg, &
          &                           ntime_alb, &
          &                           alb_field_mom, &
          &                           alnid_field_mom, &
          &                           aluvd_field_mom)
  ELSE IF (ialb_type == 3) THEN
     CALL read_netcdf_buffer_alb(alb_buffer_file,  &
          &                           tg, &
          &                           ntime_alb, &
          &                           alb_field_mom)
  ENDIF


  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('NDVI')

  CALL read_netcdf_buffer_ndvi(ndvi_buffer_file,  &
       &                                     tg,         &
       &                                     ntime_ndvi, &
       &                                     ndvi_max,  &
       &                                     ndvi_field_mom,&
       &                                     ndvi_ratio_mom)

  !-------------------------------------------------------------------------
  IF (l_use_emiss) THEN
    CALL logging%info( '')
    CALL logging%info('Emiss')
    CALL read_netcdf_buffer_emiss(TRIM(emiss_buffer_file),  &
         &                                     tg,         &
         &                                     ntime_emiss, &
         &                                     emiss_max,  &
         &                                     emiss_field_mom,&
         &                                     emiss_ratio_mom)
  END IF

  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Orography')

  CALL read_netcdf_buffer_topo(orography_buffer_file,                 &
       &                                     tg,                      &
       &                                     igrid_type,              &
       &                                     fr_land_topo,            &
       &                                     hh_topo,                 &
       &                                     stdh_topo,               &
       &                                     z0_topo,                 &
       &                                     lradtopo,                &
       &                                     lsso_param,              &
       &                                     l_use_sgsl,              &
       &                                     nhori,                   &
       &                                     vertex_param,            &
       &                                     hh_topo_max,             &
       &                                     hh_topo_min,             &           
       &                                     theta_topo,              &
       &                                     aniso_topo,              &
       &                                     slope_topo,              &
       &                                     slope_asp_topo,          &
       &                                     slope_ang_topo,          &
       &                                     horizon_topo,            &
       &                                     skyview_topo, &
       &                                     sgsl)

   IF ( (igrid_type == igrid_icon) .AND. (.NOT. lsso_param) ) THEN
     ! Provide also SSO fields, filled with zero
     theta_topo = 0._wp
     aniso_topo = 0._wp
     slope_topo = 0._wp
     CALL logging%warning('Fields theta_topo, aniso_topo and slope_topo: --> Set to 0._wp!')
   ENDIF

  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('AOT')

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

  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('CRU')

  SELECT CASE(it_cl_type)
    CASE(1)
      CALL read_netcdf_buffer_cru(t_clim_buffer_file,&
       &                                     tg,       &
       &                                     crutemp,  &
       &                                     cruelev)  
    CASE(2)
      CALL read_netcdf_buffer_cru(t_clim_buffer_file, &
       &                                     tg,        &
       &                                     crutemp)
  END SELECT

  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Flake')

  CALL read_netcdf_buffer_flake(flake_buffer_file,   &
       &                        lake_depth,&
       &                        fr_lake,   &
       &                        flake_tot_npixel)

  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Land-Sea mask')

  IF (i_lsm_data == 2 .and. igrid_type == igrid_cosmo) THEN
    CALL read_netcdf_buffer_lsm(land_sea_mask_file,  &
          &                           fr_land_mask)

    CALL logging%info('External  Land-Sea-Mask is NOT used for consistency tests.')
  ELSE
    CALL logging%info('External  Land-Sea-Mask is NOT used for consistency tests.')
    CALL logging%warning('External Land-Sea-Mask is only tested for the COSMO grid.')
  END IF
  
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info( '============= consistency check ===============')
  CALL logging%info( '')
  
  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Landuse')

  !determine land-sea mask
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
      skinc_lu          = 200.0_wp
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
  CASE (i_lu_ecci)
     lu_class_fraction(:,:,:,38)=ice_lu(:,:,:)
  END SELECT

  !----------------------- FIX for wrong classified Glacier points - Land-Use  --------------------------

  ! Calculate height-corrected annual maximum of T2M climatology, including contribution from SSO standard 
  ! deviation. This is used below to reset misclassified glacier points (e.g. salt lakes) to bare soil
  !
  ! This correction requires a monthly T2M climatology
  !

  ! climatological temperature gradient used for height correction of T2M climatology
  IF (i_landuse_data == i_lu_globcover) THEN
    count_ice2tclim = 0
    count_ice2tclim_tile = 0
    DO k=1,tg%ke
      DO j=1,tg%je
        DO i=1,tg%ie
          t2mclim_hc = MAXVAL(t2m_field(i,j,k,:)) + dtdz_clim *                &
                     ( hh_topo(i,j,k) + 1.5_wp*stdh_topo(i,j,k) - hsurf_field(i,j,k))

          ! consistency corrections for glaciered points
          ! a) set soiltype to ice if landuse = ice (already done in extpar for dominant glacier points)
          !
          ! a) plausibility check for glacier points based on T2M climatology (if available):
          !    if the warmest month exceeds 10 deg C, then it is unlikely for glaciers to exist

          IF ( t2mclim_hc > (tmelt + 10.0_wp) ) THEN
            IF (lu_class_fraction(i,j,k,i_gcv__snow_ice)> 0._wp) count_ice2tclim=count_ice2tclim + 1
            IF (lu_class_fraction(i,j,k,i_gcv__snow_ice)>= frlndtile_thrhld) THEN
               count_ice2tclim_tile = count_ice2tclim_tile + 1  ! Statistics >= frlndtile_thrhld in ICON           
            ENDIF
            lu_class_fraction(i,j,k,i_gcv_bare_soil) = lu_class_fraction(i,j,k,i_gcv_bare_soil) + &
            lu_class_fraction(i,j,k,i_gcv__snow_ice) ! add always wrong ice frac to bare soil
            lu_class_fraction(i,j,k,i_gcv__snow_ice) = 0._wp ! remove always wrong ice frac
            ice_lu(i,j,k)  = 0._wp
          ENDIF ! t2mclim_hc(i,j,k) > tmelt + 10.0_wp
        ENDDO
      ENDDO
    ENDDO
    
    WRITE(message_text,*)"Number of corrected false glacier points in EXTPAR: ", &
         &               count_ice2tclim, " with fraction >= 0.05 (TILE): ",count_ice2tclim_tile
    CALL logging%info(message_text)
  END IF
    
      IF (i_landuse_data == i_lu_ecci) THEN
         count_ice2tclim = 0
         count_ice2tclim_tile = 0
        DO k=1,tg%ke
          DO j=1,tg%je
            DO i=1,tg%ie
              t2mclim_hc = MAXVAL(t2m_field(i,j,k,:)) + dtdz_clim *                &
                         ( hh_topo(i,j,k) + 1.5_wp*stdh_topo(i,j,k) - hsurf_field(i,j,k))

              ! consistency corrections for glaciered points
              ! a) set soiltype to ice if landuse = ice (already done in extpar for dominant glacier points)
              !
              ! a) plausibility check for glacier points based on T2M climatology (if available):
              !    if the warmest month exceeds 10 deg C, then it is unlikely for glaciers to exist

              IF ( t2mclim_hc > (tmelt + 10.0_wp) ) THEN
                IF (lu_class_fraction(i,j,k,i_ecci__snow_ice)> 0._wp) count_ice2tclim=count_ice2tclim + 1
                IF (lu_class_fraction(i,j,k,i_ecci__snow_ice)>= frlndtile_thrhld) THEN
                   count_ice2tclim_tile = count_ice2tclim_tile + 1  ! Statistics >= frlndtile_thrhld in ICON           
                ENDIF
                lu_class_fraction(i,j,k,i_ecci_bare_soil) = lu_class_fraction(i,j,k,i_ecci_bare_soil) + &
                lu_class_fraction(i,j,k,i_ecci__snow_ice) ! add always wrong ice frac to bare soil
                lu_class_fraction(i,j,k,i_ecci__snow_ice) = 0._wp ! remove always wrong ice frac
                ice_lu(i,j,k)  = 0._wp
              ENDIF ! t2mclim_hc(i,j,k) > tmelt + 10.0_wp
            ENDDO
          ENDDO
        ENDDO
        
        WRITE(message_text,*)"Number of corrected false glacier points in EXTPAR: ", &
                      count_ice2tclim, " with fraction >= 0.05 (TILE): ",count_ice2tclim_tile

        CALL logging%info(message_text)

     END IF

 
  IF (tile_mode < 1) THEN   ! values are kept for ICON because of tile approach
    WHERE (fr_land_lu < 0.5)  ! set vegetation to undefined (0.) for water grid elements
      ! z0, emissivity, and skin_conductivity are not set to undefined_lu
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

  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Soil')

  !------------------------------------------------------------------------------------------
  !------------- soil data consistency ------------------------------------------------------
  !------------------------------------------------------------------------------------------

  !undef_soiltype   = 0 ! \TODO read undef_soiltype from netcdf file (_Fill_Value)
  !default_soiltype = 5 ! default soil type loam
  !soiltype_ice     = 1   !< soiltype for ice
  !soiltype_water   = 9   !< soiltype for water


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
        WHERE ((soiltype_fao == undef_soiltype).OR.(soiltype_fao > 8) ) ! land grid elements must have a valid soiltype
           !     soiltype_fao = soiltype_ice  ! set soil type to ice for Arctic or Antarctic undefined points
          WHERE ( (lat_geo < lu_data_southern_boundary) ) ! Antarctica
             soiltype_fao = soiltype_ice  ! set soil type to ice for Antarctic undefined points
      
          ELSEWHERE  ! rest of the World
            soiltype_fao = default_soiltype ! set default soiltype to loam
          ENDWHERE
        ENDWHERE
      ENDWHERE

      db_ice_counter = 0

      IF (i_landuse_data == i_lu_globcover) THEN
        DO k=1,tg%ke
          DO j=1,tg%je
            DO i=1,tg%ie
              IF  ( (soiltype_fao(i,j,k) /= soiltype_ice) .AND.  &
                  &    (fr_land_lu(i,j,k)*ice_lu(i,j,k) > 0.5)) THEN ! scale glacier frac with fr_land
                soiltype_fao(i,j,k) = soiltype_ice
                db_ice_counter = db_ice_counter +1
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        WRITE(message_text,*)'Number of grid elements set to ice soiltype: ', db_ice_counter
        CALL logging%info(message_text)

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
        WRITE(message_text,*)'Number of grid elements set to ice soiltype: ', db_ice_counter
        CALL logging%info(message_text)
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
      END IF

      CALL calculate_soiltype(tg,            &
          &                          .false.,       & ! switch off deep soil for top soil calculation
          &                          soiltype_FAO,  &
          &                          soiltype_HWSD,  &
          &                          fr_sand,       &
          &                          fr_silt,       &
          &                          fr_clay,       &
          &                          fr_oc,         &
          &                          fr_bd          )

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
        WRITE(message_text,*)'Number of grid elements set to ice soiltype: ', db_ice_counter
        CALL logging%info(message_text)

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
        WRITE(message_text,*)'Number of grid elements set to ice soiltype: ', db_ice_counter
        CALL logging%info(message_text)
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
     WRITE(message_text,*)'Number of land points with water set to  soiltype loam: ', db_water_counter
     CALL logging%info(message_text)

     WHERE (soiltype_fao == 11) ! Dunes
       soiltype_fao = 3  ! set soil type to sand for dunes
     ENDWHERE
     WHERE (soiltype_fao > 12) ! undefined
       soiltype_fao = 5  ! set soil type to loam for undefined points
     ENDWHERE

  END SELECT

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

  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Flake')

  !------------------------------------------------------------------------------------------
  !------------- flake data consistency  ----------------------------------------------------
  !------------------------------------------------------------------------------------------

  SELECT CASE(igrid_type) ! get indices for neighbour grid elements

    CASE(igrid_icon) ! ICON GRID
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
                ne_je(:) = 1_i4
                ne_ke(:) = 1_i4
                ne_ie(:) = 0_i4
                nnb=icon_grid%nvertex_per_cell ! number of neighbours in ICON grid
                DO nv=1, nnb
                  n_index = icon_grid_region%cells%neighbor_index(i,nv) ! get cell id of neighbour cells
                  IF (n_index > 0) THEN
                    ne_ie(nv) = n_index
                  ENDIF
                ENDDO
                DO n=1,nnb
                  IF ((ne_ie(n)>= 1).AND.(ne_je(n)>=1).AND.(ne_ke(n)>=1)) THEN
                    IF (fr_ocean_lu(ne_ie(n),ne_je(n),ne_ke(n))>0.5) THEN ! if the direct neighbour element is ocean,
                      fr_lake(i,j,k) = 0.0                                ! set this grid element also to ocean.
                      fr_ocean_lu(i,j,k) = 1.0 - fr_land_lu(i,j,k)
                      lake_depth(i,j,k) = flake_depth_undef ! set lake depth to flake_depth_undef (-1 m)
                    ENDIF  
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
                ne_je(1) = MAX(1_i4,j-1)
                ne_ke(1) = k
                ! north-eastern neighbour
                ne_ie(2) = MIN(tg%ie,INT(i+1,i4))
                ne_je(2) = MAX(1_i4,j-1)
                ne_ke(2) = k
                ! eastern neighbour
                ne_ie(3) = MIN(tg%ie,INT(i+1,i4))
                ne_je(3) = j
                ne_ke(3) = k
                ! south-eastern neighbour
                ne_ie(4) = MIN(tg%ie,INT(i+1,i4))
                ne_je(4) = MIN(tg%je,INT(j+1,i4))
                ne_ke(4) = k
                ! southern neighbour
                ne_ie(5) = i
                ne_je(5) = MIN(tg%je,INT(j+1,i4))
                ne_ke(5) = k
                ! south-west neighbour
                ne_ie(6) = MAX(1_i4,i-1)
                ne_je(6) = MIN(tg%je,INT(j+1,i4))
                ne_ke(6) = k
                ! western neighbour
                ne_ie(7) = MAX(1_i4,i-1)
                ne_je(7) = j
                ne_ke(7) = k
                ! north-west neighbour
                ne_ie(8) = MAX(1_i4,i-1)
                ne_je(8) = MAX(1_i4,j-1)
                ne_ke(8) = k
      
                IF (lflake_correction) THEN
                  DO n=1,nnb
                    IF ((ne_ie(n)>= 1).AND.(ne_je(n)>=1).AND.(ne_ke(n)>=1)) THEN
                      IF (fr_ocean_lu(ne_ie(n),ne_je(n),ne_ke(n))>0.5) THEN ! if the direct neighbour element is ocean,
                        fr_lake(i,j,k) = 0.0                                ! set this grid element also to ocean.
                        IF ((i==391).AND.(j==267))THEN
                          WRITE(message_text,*)'changed: ',                    &
                          ne_ie(n),ne_je(n),ne_ke(n),fr_ocean_lu(ne_ie(n),ne_je(n),ne_ke(n))
                          CALL logging%info(message_text)
                        ENDIF
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

  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Albedo')

  !------------------------------------------------------------------------------------------
  !------------- Albedo data consistency ------------------------------------------------------
  !------------------------------------------------------------------------------------------
  ! set default Albedo values for land grid elements with so far undefined or unexpected values
  IF (ialb_type /= 2) THEN

     ! set default Albedo values for land grid elements with so far undefined or unexpected values

    namelist_alb_data_input = 'INPUT_ALB'

    CALL  read_namelists_extpar_alb(namelist_alb_data_input, &
         &                                  raw_data_alb_path, &
         &                                  raw_data_alb_filename, &
         &                                  raw_data_alnid_filename, &
         &                                  raw_data_aluvd_filename, &
         &                                  ialb_type, &
         &                                  alb_buffer_file, &
         &                                  alb_output_file)

    path_alb_file = TRIM(raw_data_alb_path)//TRIM(raw_data_alb_filename)

    path_alb_file = TRIM(raw_data_alb_path)//TRIM(raw_data_alb_filename)

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
  ENDIF

  !------------------------------------------------------------------------------------------
  !------------- soil albedo consistency check ----------------------------------------------
  !------------------------------------------------------------------------------------------
  IF (ialb_type == 2) THEN

     ! set default soil albedo values for land grid elements with so far undefined values
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

  ENDIF

  !------------------------------------------------------------------------------------------
  !------------- ISA/AHF data consistency ---------------------------------------------------
  !------------------------------------------------------------------------------------------

  IF (l_use_isa.AND.l_use_ahf) THEN

    !-------------------------------------------------------------------------
    CALL logging%info( '')
    CALL logging%info('AHF/ISA')

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

     WHERE (fr_land_lu < 0.5)  ! set water soiltype for water grid elements
       isa_field=0.
     ENDWHERE
  END IF

  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('NDVI')

  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  !------------- NDVI data consistency ------------------------------------------------------
  !------------------------------------------------------------------------------------------

  ! set default NDVI values for land grid elements with so far undefined values or very small NDVI values

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

  IF (l_use_emiss) THEN
  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Emiss')
  
  !------------------------------------------------------------------------------------------
  !------------- EMISS data consistency ------------------------------------------------------
  !------------------------------------------------------------------------------------------

  ! set default EMISS values for land grid elements with so far undefined values or very small EMISS values
  minimal_emiss = 0.8 

  FORALL (t=1:mpy) ! mpy = month per year = 12
    WHERE (fr_land_lu < MERGE(0.01,0.5,tile_mask)) ! set undefined EMISS value (0.0) for water grid elements
      emiss_field_mom(:,:,:,t) = 0.991
    ELSEWHERE ! fr_land_lu >= 0.5
      WHERE (emiss_max(:,:,:) <= minimal_emiss) ! small EMISS values at land grid elements
        emiss_field_mom(:,:,:,t) = 0.95 ! bare areas 
      ENDWHERE
      WHERE (emiss_field_mom(:,:,:,t) <= minimal_emiss) ! small EMISS values at land grid elements
        emiss_field_mom(:,:,:,t) = 0.95 ! bare areas 
      ENDWHERE
    ENDWHERE
  END FORALL

  END IF

  !------------------------------------------------------------------------------------------
  !-------------TC_L Correction ------------------------------------------------------
  !------------------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('T_CL correction')

  IF (igrid_type == igrid_cosmo) THEN
    SELECT CASE(it_cl_type)
      CASE(1)

        crutemp2 = crutemp
        DO j=1,tg%je
          DO i=1,tg%ie
            last = .FALSE.
            IF ( fr_land_lu(i,j,1) < 0.5) THEN
              crutemp(i,j,1)  = -1.E20_wp
            ELSE
              IF ( crutemp(i,j,1) > 0.0 ) THEN
                foundtcl = .TRUE.
                crutemp(i,j,1) = crutemp(i,j,1) + 0.65 * 0.01*( cruelev(i,j,1) - hh_topo(i,j,1) )
              ELSE
                ! 3x3 search
                foundtcl = .FALSE.
                DO jj=-1,2
                  DO ii=-1,2
                    IF (j+jj > 0 .and.  j+jj < tg%je .and. i+ii > 0 .and. i+ii < tg%ie) THEN
                      IF ( crutemp2(i+ii,j+jj,1) > 0.0 ) THEN
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
                    iml=MAX(1_i4,i-3*l)
                    imu=i+2-3*l
                    ipl=i+3*l-2
                    ipu=MIN(tg%ie,i+3*l)
                    jml=MAX(1_i4,j-l)
                    jmu=j-l
                    jpl=j+l
                    jpu=MIN(tg%je,INT(j+l,i4))
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
                  WRITE(message_text,*) 'No temperature data for T_CL correction  at point: ', i,j
                  CALL logging%warning(message_text)
                  crutemp(i,j,1) = 288.15 - 0.0065 * hh_topo(i,j,1)
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
    END SELECT
  ENDIF

  SELECT CASE(isoil_data)
    CASE(HWSD_data)
      CALL logging%info('Selected HWSD => Copy HWSD data for Output')
  END SELECT


  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Special points')

  IF (igrid_type == igrid_cosmo .OR. igrid_type == igrid_icon) THEN
    DO isp = 1, number_special_points
      IF (number_special_points<1) THEN
        WRITE(message_text,*)'No treatment of special points: Number of special points is ',number_special_points
        CALL logging%error(message_text,__FILE__,__LINE__)
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
         CALL logging%warning(' Special points defined but not in target domain!')
       ELSE
         start_cell_id = 1

         CALL  find_nearest_target_grid_element(lon_geo_sp, &
               & lat_geo_sp, &
               & tg,            &
               & start_cell_id, &
               & i_sp,      &
               & j_sp,      &
               & k_sp)

          WRITE(message_text,'(A)') '-------------------------------------------------------------------------------------'
          CALL logging%info(message_text)
          WRITE(message_text,'(A26,A10,A4,2X,I1)')  "Consider special point in ",namelist_file," of ",number_special_points
          CALL logging%info(message_text)
          WRITE(message_text,'(A33,1X,2(F6.3,2X))') "         special point position (lon,lat): ",lon_geo_sp,lat_geo_sp
          CALL logging%info(message_text)
          WRITE(message_text,'(A33,1X,3(I9,2X))')   "         special point index (ie,je,ke):   ",i_sp,j_sp,k_sp
          CALL logging%info(message_text)
          IF ((i_sp == 0).OR.(j_sp == 0)) THEN
            WRITE(message_text,*)"WARNING: Special points out of range of target domain!"
            CALL logging%info(message_text)
          ELSE
            WRITE(message_text,'(A23,I9,2X,I7,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                  i_sp,j_sp," z0_tot old ",z0_tot (i_sp,j_sp,k_sp),"new ",z0_sp
            CALL logging%info(message_text)
            WRITE(message_text,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                  i_sp,j_sp," root_lu old ",root_lu(i_sp,j_sp,k_sp),"new ",rootdp_sp
            CALL logging%info(message_text)
            WRITE(message_text,'(A23,I9,2X,I9,A19,I5,2X,A4,I5)')"         special point: ",&
                  i_sp,j_sp," soiltype_fao old  ",soiltype_fao(i_sp,j_sp,k_sp),"new ",NINT(soiltype_sp)
            CALL logging%info(message_text)
            WRITE(message_text,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                  i_sp,j_sp," plcov_mn_lu  old  ",plcov_mn_lu (i_sp,j_sp,k_sp),"new ",plcovmn_sp
            CALL logging%info(message_text)
            WRITE(message_text,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                  i_sp,j_sp," plcov_mx_lu  old  ",plcov_mx_lu (i_sp,j_sp,k_sp),"new ",plcovmx_sp
            CALL logging%info(message_text)
            WRITE(message_text,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                  i_sp,j_sp," lai_mn_lu    old  ",lai_mn_lu (i_sp,j_sp,k_sp),"new ",laimn_sp
            CALL logging%info(message_text)
            WRITE(message_text,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                  i_sp,j_sp," lai_mx_lu    old  ",lai_mx_lu (i_sp,j_sp,k_sp),"new ",laimx_sp
            CALL logging%info(message_text)
            IF (for_d_sp >= 0._wp) THEN
              WRITE(message_text,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                  i_sp,j_sp," for_d_lu    old  ", for_d_lu(i_sp,j_sp,k_sp),"new ",for_d_sp
              CALL logging%info(message_text)
            ENDIF
            IF (for_e_sp >= 0._wp)THEN
              WRITE(message_text,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                  i_sp,j_sp," for_e_lu    old  ", for_e_lu(i_sp,j_sp,k_sp),"new ",for_e_sp
              CALL logging%info(message_text)
            ENDIF
            IF (fr_land_sp >= 0._wp)THEN
              WRITE(message_text,'(A23,I9,2X,I9,A19,F6.4,2X,A4,F6.4)')"         special point: ",&
                  i_sp,j_sp," fr_land    old  ", fr_land_lu(i_sp,j_sp,k_sp),"new ",fr_land_sp
              CALL logging%info(message_text)
            ENDIF
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
                  WRITE(message_text,'(A33,1X,A85,2X,F8.4)') "Land-Use Fractions for GLOBCOVER class  ", &
                         glc_class(i),lu_class_fraction(i_sp,j_sp,k_sp,i)
                  CALL logging%info(message_text)
                ENDDO

              CASE (i_lu_ecci)

                ecci_class(1)=   'No data                                        '     ! 1.
                ecci_class(2)=   'Cropland, rainfed                              '     ! 2.
                ecci_class(3)=   'Herbaceous cover                               '     ! 3.
                ecci_class(4)=   'Tree or shrub cover                            '     ! 4.
                ecci_class(5)=   'Cropland, irrigated or post-flooding           '     ! 5.
                ecci_class(6)=   'Mosaic cropland(>50%)/natural vegetation(<50%) '     ! 6.
                ecci_class(7)=   'Mosaic natural vegetation(>50%)/cropland(<50%) '     ! 7.
                ecci_class(8)=   'TC, BL, evergreen, closed to open (>15%)       '     ! 8.
                ecci_class(9)=   'TC, BL, deciduous, closed to open (>15%)       '     ! 9.
                ecci_class(10)=   'TC, BL, deciduous, closed (>40%)               '     ! 10.
                ecci_class(11)=   'TC, BL, deciduous, open (15-40%)               '     ! 11.
                ecci_class(12)=   'TC, NL, evergreen, closed to open (>15%)       '     ! 12.
                ecci_class(13)=   'TC, NL, evergreen, closed (>40%)               '     ! 13.
                ecci_class(14)=   'TC, NL, evergreen, open (15-40%)               '     ! 14.
                ecci_class(15)=   'TC, NL, deciduous, closed to open (>15%)       '     ! 15.
                ecci_class(16)=   'TC, NL, deciduous, closed (>40%)               '     ! 16.
                ecci_class(17)=   'TC, NL, deciduous, open (15-40%)               '     ! 17.
                ecci_class(18)=   'TC, mixed leaf type (BL and NL)                '     ! 18.
                ecci_class(19)=   'Mosaic tree and shrub (>50%)/HC (<50%)         '     ! 19.
                ecci_class(20)=   'Mosaic HC (>50%) / tree and shrub (<50%)       '     ! 20.
                ecci_class(21)=   'Shrubland                                      '     ! 21.
                ecci_class(22)=   'Shrubland evergreen                            '     ! 22.
                ecci_class(23)=   'Shrubland deciduous                            '     ! 23.
                ecci_class(24)=   'Grassland                                      '     ! 24.
                ecci_class(25)=   'Lichens and mosses                             '     ! 25.
                ecci_class(26)=   'Sparse vegetation (tree, shrub, HC) (<15%)     '     ! 26.
                ecci_class(27)=   'Sparse tree (<15%)                             '     ! 27.
                ecci_class(28)=   'Sparse shrub (<15%)                            '     ! 28.
                ecci_class(29)=   'Sparse herbaceous cover (<15%)                 '     ! 29.
                ecci_class(30)=   'TC, flooded, fresh or brakish water            '     ! 30.
                ecci_class(31)=   'TC, flooded, saline water                      '     ! 31.
                ecci_class(32)=   'Shrub or HC,flooded,fresh/saline/brakish water '     ! 32.
                ecci_class(33)=   'Urban areas                                    '     ! 33.
                ecci_class(34)=   'Bare areas                                     '     ! 34.
                ecci_class(35)=   'Consolidated bare areas                        '     ! 35.
                ecci_class(36)=   'Unconsolidated bare areas                      '     ! 36.
                ecci_class(37)=   'Water bodies                                   '     ! 37.
                ecci_class(38)=   'Permanent snow and ice                         '     ! 38.

                DO i=1,nclass_ecci
                  WRITE(message_text,'(A33,1X,A85,2X,F8.4)') "Land-Use Fractions for ESA CCI class  ", &
                        ecci_class(i),lu_class_fraction(i_sp,j_sp,k_sp,i)
                  CALL logging%info(message_text)
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
       END IF
    END DO ! Special Points  loop
  END IF

  !-------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info('Orography')

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

  WRITE(message_text,*)TRIM(y_orofilter)
  CALL logging%info(message_text)

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= write data to netcdf==============')
  CALL logging%info( '')

  fill_value_real = -1.E20_wp
  fill_value_int = -999

  SELECT CASE(igrid_type)
  CASE(igrid_icon) ! ICON GRID
    CALL write_cdi_icon_grid_extpar(TRIM(netcdf_output_filename),&
         &                                     icon_grid,                     &
         &                                     tg,                            &
         &                                     isoil_data,                    &
         &                                     ldeep_soil,                    &
         &                                     itopo_type,                    &
         &                                     lsso_param,                    &
         &                                     l_use_isa,                     &
         &                                     l_use_ahf,                     &
         &                                     l_use_emiss,                   &
         &                                     lradtopo,                      &
         &                                     nhori,                         &
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
         &                                     skinc_lu,                      &
         &                                     emissivity_lu,                 &
         &                                     lake_depth,                    &
         &                                     fr_lake,                       &
         &                                     soiltype_fao,                  &
         &                                     ndvi_max,                      &
         &                                     ndvi_field_mom,                &
         &                                     ndvi_ratio_mom,                &
         &                                     emiss_field_mom,               &
         &                                     hh_topo,                       &
         &                                     hh_topo_max,                   &
         &                                     hh_topo_min,                   &         
         &                                     stdh_topo,                     &
         &                                     theta_topo,                    &
         &                                     aniso_topo,                    &
         &                                     slope_topo,                    &
         &                                     aot_tg,                        &
         &                                     crutemp,                       &
         &                                     alb_field_mom,                 &
         &                                     alnid_field_mom,               &
         &                                     aluvd_field_mom,               &
         &                                     fr_sand = fr_sand,             &
         &                                     fr_silt = fr_silt,             &
         &                                     fr_clay = fr_clay,             &
         &                                     fr_oc = fr_oc,                 &
         &                                     fr_bd = fr_bd,                 &
         &                                     soiltype_deep=soiltype_deep,   &
         &                                     fr_sand_deep=fr_sand_deep,     &
         &                                     fr_silt_deep=fr_silt_deep,     &
         &                                     fr_clay_deep=fr_clay_deep,     &
         &                                     fr_oc_deep=fr_oc_deep,         &
         &                                     fr_bd_deep=fr_bd_deep,         &
         &                                     isa_field=isa_field,           &
         &                                     ahf_field=ahf_field,           &
         &                                     sst_field=sst_field,           &
         &                                     wsnow_field=wsnow_field,       &
         &                                     t2m_field=t2m_field,           &
         &                                     hsurf_field=hsurf_field,       &
         &                                     horizon_topo=horizon_topo,     &
         &                                     skyview_topo=skyview_topo      )

    CASE(igrid_cosmo) ! COSMO grid

    IF(ldeep_soil) THEN
      CALL  write_netcdf_cosmo_grid_extpar(TRIM(netcdf_output_filename),      &
         &                                     cosmo_grid,                    &
         &                                     tg,                            &
         &                                     isoil_data,                    &
         &                                     ldeep_soil,                    &
         &                                     itopo_type,                    &
         &                                     lsso_param,                    &
         &                                     l_use_isa,                     &
         &                                     l_use_ahf,                     &
         &                                     l_use_sgsl,                    &
         &                                     lscale_separation,             &
         &                                     TRIM(y_orofilter),             &
         &                                     lradtopo,                      &
         &                                     nhori,                         &
         &                                     fill_value_real,               &
         &                                     TRIM(name_lookup_table_lu),    &
         &                                     TRIM(lu_dataset),              &
         &                                     i_landuse_data,                &
         &                                     nclass_lu,                     &
         &                                     lon_geo,                       &
         &                                     lat_geo,                       &
         &                                     fr_land_lu,                    &
         &                                     lu_class_fraction,             &
         &                                     ice_lu,                        &
         &                                     z0_tot,                        &
         &                                     z0_topo,                       &
         &                                     z012_lu,                       &
         &                                     root_lu,                       &
         &                                     plcov_mn_lu,                   &
         &                                     plcov_mx_lu,                   &
         &                                     plcov12_lu,                    &
         &                                     lai_mn_lu,                     &
         &                                     lai_mx_lu,                     &
         &                                     lai12_lu,                      &
         &                                     rs_min_lu,                     &
         &                                     urban_lu,                      &
         &                                     for_d_lu,                      &
         &                                     for_e_lu,                      &
         &                                     skinc_lu,                      &
         &                                     emissivity_lu,                 &
         &                                     lake_depth,                    &
         &                                     fr_lake,                       &
         &                                     soiltype_fao,                  &
         &                                     ndvi_max,                      &
         &                                     ndvi_field_mom,                &
         &                                     ndvi_ratio_mom,                &
         &                                     emiss_field_mom,               &
         &                                     hh_topo,                       &
         &                                     stdh_topo,                     &
         &                                     aot_tg,                        &
         &                                     MAC_aot_tg,                    &
         &                                     MAC_ssa_tg,                    &
         &                                     MAC_asy_tg,                    &
         &                                     crutemp,                       &
         &                                     alb_field_mom,                 &
         &                                     alnid_field_mom,               &
         &                                     aluvd_field_mom,               &
         &                                     alb_dry = alb_dry,             &
         &                                     alb_sat = alb_sat,             &
         &                                     fr_sand = fr_sand,             &
         &                                     fr_silt = fr_silt,             &
         &                                     fr_clay = fr_clay,             &
         &                                     fr_oc = fr_oc,                 &
         &                                     fr_bd = fr_bd,                 &
         &                                     soiltype_deep=soiltype_deep,   &
         &                                     fr_sand_deep=fr_sand_deep,     &
         &                                     fr_silt_deep=fr_silt_deep,     &
         &                                     fr_clay_deep=fr_clay_deep,     &
         &                                     fr_oc_deep=fr_oc_deep,         &
         &                                     fr_bd_deep=fr_bd_deep,         &
         &                                     theta_topo=theta_topo,         &
         &                                     aniso_topo=aniso_topo,         &
         &                                     slope_topo=slope_topo,         &
         &                                     slope_asp_topo=slope_asp_topo, &
         &                                     slope_ang_topo=slope_ang_topo, &
         &                                     horizon_topo=horizon_topo,     &
         &                                     skyview_topo=skyview_topo,     &
         &                                     isa_field=isa_field,           &
         &                                     ahf_field=ahf_field,           &
         &                                     sgsl = sgsl                    )

      ELSE
        CALL  write_netcdf_cosmo_grid_extpar(TRIM(netcdf_output_filename),    &
         &                                     cosmo_grid,                    &
         &                                     tg,                            &
         &                                     isoil_data,                    &
         &                                     ldeep_soil,                    &
         &                                     itopo_type,                    &
         &                                     lsso_param,                    &
         &                                     l_use_isa,                     &
         &                                     l_use_ahf,                     &
         &                                     l_use_sgsl,                    &
         &                                     lscale_separation,             &
         &                                     y_orofilter,                   &
         &                                     lradtopo,                      &
         &                                     nhori,                         &
         &                                     fill_value_real,               &
         &                                     TRIM(name_lookup_table_lu),    &
         &                                     TRIM(lu_dataset),              &
         &                                     i_landuse_data,                &
         &                                     nclass_lu,                     &
         &                                     lon_geo,                       &
         &                                     lat_geo,                       &
         &                                     fr_land_lu,                    &
         &                                     lu_class_fraction,             &
         &                                     ice_lu,                        &
         &                                     z0_tot,                        &
         &                                     z0_topo,                       &
         &                                     z012_lu,                       &
         &                                     root_lu,                       &
         &                                     plcov_mn_lu,                   &
         &                                     plcov_mx_lu,                   &
         &                                     plcov12_lu,                    &
         &                                     lai_mn_lu,                     &
         &                                     lai_mx_lu,                     &
         &                                     lai12_lu,                      &
         &                                     rs_min_lu,                     &
         &                                     urban_lu,                      &
         &                                     for_d_lu,                      &
         &                                     for_e_lu,                      &
         &                                     skinc_lu,                      &
         &                                     emissivity_lu,                 &
         &                                     lake_depth,                    &
         &                                     fr_lake,                       &
         &                                     soiltype_fao,                  &
         &                                     ndvi_max,                      &
         &                                     ndvi_field_mom,                &
         &                                     ndvi_ratio_mom,                &
         &                                     emiss_field_mom,               &
         &                                     hh_topo,                       &
         &                                     stdh_topo,                     &
         &                                     aot_tg,                        &
         &                                     MAC_aot_tg,                    &
         &                                     MAC_ssa_tg,                    &
         &                                     MAC_asy_tg,                    &
         &                                     crutemp,                       &
         &                                     alb_field_mom,                 &
         &                                     alnid_field_mom,               &
         &                                     aluvd_field_mom,               &
         &                                     alb_dry = alb_dry,             &
         &                                     alb_sat = alb_sat,             &
         &                                     fr_sand = fr_sand,             &
         &                                     fr_silt = fr_silt,             &
         &                                     fr_clay = fr_clay,             &
         &                                     fr_oc = fr_oc,                 &
         &                                     fr_bd = fr_bd,                 &
         &                                     theta_topo=theta_topo,         &
         &                                     aniso_topo=aniso_topo,         &
         &                                     slope_topo=slope_topo,         &
         &                                     slope_asp_topo=slope_asp_topo, &
         &                                     slope_ang_topo=slope_ang_topo, &
         &                                     horizon_topo=horizon_topo,     &
         &                                     skyview_topo=skyview_topo,     &
         &                                     isa_field=isa_field,           &
         &                                     ahf_field=ahf_field,           &
         &                                     sgsl = sgsl                    )
      ENDIF
  END SELECT

  WRITE(logging%fileunit,*) ''
  WRITE(logging%fileunit,*)'============= consistency_check done ============'

END PROGRAM
