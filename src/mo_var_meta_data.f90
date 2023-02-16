!+ Fortran module with definitions of meta information for variables and dimension for the output
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio
!  Add meta information for aerosol types
! V1_2         2011/03/25 Hermann Asensio
!  Use capital letters in variable name for output variables, but minuscule letters for coordinates
!  Adapt variable names to COSMO-CLM conventions and "standard name" to cf convention
! V1_3         2011/04/19 Hermann Asensio
!  introduce Globcover 2009 land use data set for external parameters
! V1_4         2011/04/21 Anne Roches
!  implementation of orography smoothing
! V1_8         2013-03-12 Frank Brenner
!  introduced MODIS albedo dataset(s) as new external parameter(s)
! V1_11        2013/04/16 Juergen Helmert
!  Adaptions for external land-sea mask
! V2_0         2013-06-04 Anne Roches
!  introduced topographical corrected radiation parameters
! V2_0         2013-06-04 Martina Messmer/Daniel Luethi
!  introduced HWSD soil data set as new external parameters (Juergen Helmert)
!  contains the topsoil and the subsoil
!  allow for different AOT climatologies
! V2_0_3       2015-01-12 Juergen Helmert
!  Bugfix correction covers CSCS SVN r5907-r6359
! V3_0         2015-05-21 Juergen Helmert
!  Add information for urban fields ISA and AHF
! V3_1         2016-04-14 Burkhardt Rockel
!  Add meta information for SA and AHF fields
! V4_0         2016-08-05 Daniel Luethi
!  Add information for subgrid scale slope fields
!              2016-08-23 authors from RHM and Daniel Lthi
!  Add information for MACv2 aerosol fields (iaot_type == 4)
! ------------ ---------- ----
! V5_0         2023/02/05 Andrzej Wyszogrodzki
!  added ECOSG options
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with definitions of meta information for variables and dimension for the output
!> \author Hermann Asensio
MODULE mo_var_meta_data

  USE mo_logging
  USE mo_kind,                  ONLY: i4

  USE mo_io_utilities,          ONLY: dim_meta_info, var_meta_info, &
       &                              vartype_int, vartype_real,    &
       &                              netcdf_grid_mapping

  USE mo_grid_structures,       ONLY: target_grid_def, &
       &                              rotated_lonlat_grid, &
       &                              icosahedral_triangular_grid

  USE mo_topo_data,             ONLY: itype_scaling
  USE mo_python_data,           ONLY: iera_type, isa_type, iahf_type

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: &

            ! grid specification
       &    dim_2d_tg, dim_3d_tg, dim_4d_tg, def_dimension_info_buffer, &
       &    dim_rlon_cosmo, dim_rlat_cosmo, dim_nhori_cosmo, dim_2d_cosmo, dim_3d_cosmo, def_dimension_info_cosmo, &
       &    rlon_meta, rlat_meta, &
       &    dim_icon, def_dimension_info_icon, &
       &    clon_meta, clat_meta, &
       &    clon_vertices_meta, clat_vertices_meta, &
       &    dim_buffer_cell, dim_buffer_vertex, &
       &    lon_geo_meta, lat_geo_meta, no_raw_data_pixel_meta, def_com_target_fields_meta, &
       &    crutemp_meta, def_crutemp_meta, &
       &    cruelev_meta, def_cruelev_meta, &
       &    nc_grid_def_cosmo, set_nc_grid_def_cosmo, &
       &    nc_grid_def_icon, set_nc_grid_def_icon, &

            ! landuse
       &    dim_lu_tg, dim_ecoclimap_tg, dim_ecoclimap_tg2, &
       &    fr_land_lu_meta, fr_land_mask_meta,lu_tot_npixel_meta, &
       &    lu_class_fraction_meta, lu_class_npixel_meta, &
       &    ice_lu_meta, z0_lu_meta, &
       &    plcov_mx_lu_meta, plcov_mn_lu_meta, &
       &    lai_mx_lu_meta, lai_mn_lu_meta, &
       &    rs_min_lu_meta, urban_lu_meta, &
       &    for_d_lu_meta, for_e_lu_meta, &
       &    skinc_lu_meta, &
       &    emissivity_lu_meta, root_lu_meta, &
       &    fr_ocean_lu_meta, &
       &    plcov12_lu_meta, lai12_lu_meta, &
       &    z012_lu_meta, z012_tot_meta, &
       &    def_lu_fields_meta, &
       &    def_ecoclimap_fields_meta, &
       &    def_glc2000_fields_meta, &
       &    fr_land_glc2000_meta, glc2000_tot_npixel_meta, &
       &    glc2000_class_fraction_meta, glc2000_class_npixel_meta, &
       &    ice_glc2000_meta, z0_glc2000_meta, &
       &    plcov_mx_glc2000_meta, plcov_mn_glc2000_meta, &
       &    lai_mx_glc2000_meta, lai_mn_glc2000_meta, &
       &    rs_min_glc2000_meta, urban_glc2000_meta, &
       &    for_d_glc2000_meta, for_e_glc2000_meta, &
       &    emissivity_glc2000_meta, root_glc2000_meta, &
       &    dim_glcc_tg, &
       &    fr_land_glcc_meta, glcc_tot_npixel_meta, &
       &    glcc_class_fraction_meta, glcc_class_npixel_meta, &
       &    ice_glcc_meta, z0_glcc_meta, &
       &    plcov_mx_glcc_meta, plcov_mn_glcc_meta, &
       &    lai_mx_glcc_meta, lai_mn_glcc_meta, &
       &    rs_min_glcc_meta, urban_glcc_meta, &
       &    for_d_glcc_meta, for_e_glcc_meta, &
       &    emissivity_glcc_meta, root_glcc_meta, &
       &    def_glcc_fields_meta, &
       &    dim_glc2000_tg, &
       &    dim_ecosg_tg, &
       &    def_ecosg_fields_meta, &
       &    fr_land_ecosg_meta, &
       &    ecosg_tot_npixel_meta, &
       &    ecosg_class_fraction_meta, &
       &    ecosg_class_npixel_meta, &
       &    ice_ecosg_meta, &
       &    z0_ecosg_meta, &
       &    plcov_mx_ecosg_meta, &
       &    plcov_mn_ecosg_meta, &
       &    lai_mx_ecosg_meta, &
       &    lai_mn_ecosg_meta, &
       &    rs_min_ecosg_meta, &
       &    urban_ecosg_meta, &
       &    for_d_ecosg_meta, &
       &    for_e_ecosg_meta, &
       &    skinc_ecosg_meta, &
       &    emissivity_ecosg_meta, &
       &    root_ecosg_meta, &

            ! topography
       &    hh_topo_meta, fr_land_topo_meta,          &
       &    hh_topo_max_meta, hh_topo_min_meta,       &
       &    stdh_topo_meta, theta_topo_meta,          &
       &    aniso_topo_meta, slope_topo_meta,         &
       &    hh_vert_meta, npixel_vert_meta,           &
       &    hh_fis_meta, z0_topo_meta,                &
       &    slope_asp_topo_meta, slope_ang_topo_meta, &
       &    horizon_topo_meta, skyview_topo_meta, &
       &    def_topo_meta, def_topo_vertex_meta, &
       &    sgsl_meta, &

            ! soil
       &    def_soil_meta, &
       &    fr_land_soil_meta, soiltype_fao_meta, soiltype_hwsd_meta, &
       &    soiltype_FAO_deep_meta,soiltype_HWSD_deep_meta, &
       &    HWSD_SAND_meta, HWSD_SILT_meta, HWSD_CLAY_meta, &
       &    HWSD_OC_meta, HWSD_BD_meta,HWSD_DM_meta, &
       &    HWSD_SAND_DEEP_meta, HWSD_SILT_DEEP_meta, HWSD_CLAY_DEEP_meta, &
       &    HWSD_OC_DEEP_meta, HWSD_BD_DEEP_meta,HWSD_DM_DEEP_meta, &
       &    lake_depth_meta, fr_lake_meta, flake_tot_npixel_meta, &
       &    def_flake_fields_meta, &
       &    def_lsm_fields_meta, &

            ! ahf/isa
       &    dim_isa_tg, &
       &    isa_field_meta, &
       &    def_isa_fields_meta, &
       &    dim_ahf_tg, def_ahf_meta, &
       &    ahf_field_meta, &

            ! ndvi
       &    dim_ndvi_tg, def_ndvi_meta, &
       &    ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta, &
       &    dim_emiss_tg, def_emiss_meta, &
       &    emiss_max_meta, emiss_field_mom_meta, emiss_ratio_mom_meta, &
       &    dim_era_tg, def_era_meta, &
       &    sst_field_meta, wsnow_field_meta, t2m_field_meta, hsurf_field_meta, &
       &    dim_alb_tg, def_alb_meta, &

            ! albedo
       &    alb_field_mom_meta, &
       &    dim_aot_tg, dim_aot_ty, &
       &    aot_tg_meta, aer_bc_meta, aer_dust_meta, aer_org_meta, aer_so4_meta, aer_ss_meta, &
       &    aot_tg_MAC_meta, ssa_tg_MAC_meta, asy_tg_MAC_meta, &
       &    CAMS_SS1_tg_meta,CAMS_SS2_tg_meta,CAMS_SS3_tg_meta, &
       &    CAMS_DUST1_tg_meta,CAMS_DUST2_tg_meta,CAMS_DUST3_tg_meta, &
       &    CAMS_OCphilic_tg_meta,CAMS_OCphobic_tg_meta,&
       &    CAMS_BCphilic_tg_meta,CAMS_BCphobic_tg_meta,&
       &    CAMS_SU_tg_meta,CAMS_plev_tg_meta,&
       &    def_aot_tg_meta, &
       &    aot_type_shortname, &
       &    alnid_field_mom_meta, &
       &    aluvd_field_mom_meta, &
       &    alb_interpol_meta, &
       &    alb_dry_meta, alb_sat_meta

  TYPE(dim_meta_info), TARGET              :: dim_2d_tg(1:2), &
       &                                      dim_3d_tg(1:3), &
       &                                      dim_4d_tg(1:4), &
       &                                      dim_rlon_cosmo(1:1), &
       &                                      dim_rlat_cosmo(1:1), &
       &                                      dim_nhori_cosmo(1:1), &
       &                                      dim_2d_cosmo(1:2), &
       &                                      dim_3d_cosmo(1:3), &
       &                                      dim_icon(1:7), &
       &                                      dim_cells_icon(1:1), &
       &                                      dim_2d_icon(1:2), &
       &                                      dim_buffer_cell(1:3), &
       &                                      dim_buffer_vertex(1:3)

  TYPE(dim_meta_info), TARGET, ALLOCATABLE :: dim_aot_tg(:), & !< dimensions for field with all aerosol types
       &                                      dim_aot_ty(:), & !< dimensions for fields with single aerosol types
       &                                      dim_glc2000_tg(:), &
       &                                      dim_glcc_tg(:), &
       &                                      dim_ecosg_tg(:), &
       &                                      dim_lu_tg(:), &
       &                                      dim_isa_tg(:), &
       &                                      dim_ahf_tg(:), &
       &                                      dim_ecoclimap_tg(:), &
       &                                      dim_ecoclimap_tg2(:), &
       &                                      dim_ndvi_tg(:), &
       &                                      dim_emiss_tg(:), &
       &                                      dim_era_tg(:), &
       &                                      dim_alb_tg(:)

  TYPE(var_meta_info)                      :: aot_tg_meta, & !< variable aot_tg with all aerosol fields
       &                                      aer_bc_meta, & !< variable with aerosol optical thickness of black carbon
       &                                      aer_dust_meta, & !< variable with aerosol optical thickness of dust
       &                                      aer_org_meta, & !< variable with aerosol optical thickness of organic matter
       &                                      aer_so4_meta, & !< variable with aerosol optical thickness of sulfate
       &                                      aer_ss_meta, & !< avariable with aerosol optical thickness of sea salt
       &                                      aot_tg_MAC_meta, & !< meta data for MACv2 AOT field
       &                                      ssa_tg_MAC_meta, & !< meta data for MACv2 SSA field
       &                                      asy_tg_MAC_meta, & !< meta data for MACv2 ASY field
       &                                      CAMS_SS1_tg_meta, & !< meta data for CAMS aerosols
       &                                      CAMS_SS2_tg_meta, & !< meta data for CAMS aerosols
       &                                      CAMS_SS3_tg_meta, & !< meta data for CAMS aerosols
       &                                      CAMS_DUST1_tg_meta, & !< meta data for CAMS aerosols
       &                                      CAMS_DUST2_tg_meta, & !< meta data for CAMS aerosols
       &                                      CAMS_DUST3_tg_meta, & !< meta data for CAMS aerosols
       &                                      CAMS_OCphilic_tg_meta, & !< meta data for CAMS aerosols
       &                                      CAMS_OCphobic_tg_meta, & !< meta data for CAMS aerosols
       &                                      CAMS_BCphilic_tg_meta, & !< meta data for CAMS aerosols
       &                                      CAMS_BCphobic_tg_meta, & !< meta data for CAMS aerosols
       &                                      CAMS_SU_tg_meta, & !< meta data for CAMS aerosols
       &                                      CAMS_plev_tg_meta, & !< meta data for CAMS aerosols
       &                                      ahf_field_meta, & !< additional information for variable
       &                                      sst_field_meta, & !< additional information for variable
       &                                      wsnow_field_meta, & !< additional information for variable
       &                                      t2m_field_meta, & !< additional information for variable
       &                                      hsurf_field_meta, & !< additional information for variable
       &                                      ndvi_max_meta, & !< additional information for variable
       &                                      ndvi_field_mom_meta, & !< additional information for variable
       &                                      ndvi_ratio_mom_meta, & !< additional information for variable
       &                                      emiss_max_meta, & !< additional information for variable
       &                                      emiss_field_mom_meta, & !< additional information for variable
       &                                      emiss_ratio_mom_meta, & !< additional information for variable
       &                                      alb_field_mom_meta, & !< additional information for variable
       &                                      alnid_field_mom_meta, & !< additional information for variable
       &                                      aluvd_field_mom_meta, & !< additional information for variable
       &                                      alb_interpol_meta, & !< additional information for variable
       &                                      alb_dry_meta, &
       &                                      alb_sat_meta, &
       &                                      crutemp_meta, & !< additional information for variable crutemp
       &                                      cruelev_meta, & !< additional information for variable cruelev
       &                                      lon_geo_meta, & !< additional information for variable lon_geo_meta
       &                                      lat_geo_meta, & !< additional information for variable lat_geo_meta
       &                                      no_raw_data_pixel_meta, & !< additional information for variable no_raw_data_pixel_meta
       &                                      rlon_meta , &
       &                                      rlat_meta , &
       &                                      clon_meta , &
       &                                      clat_meta , &
       &                                      clon_vertices_meta , &
       &                                      clat_vertices_meta , &
       &                                      fr_land_glc2000_meta  , &
       &                                      glc2000_tot_npixel_meta , &
       &                                      glc2000_class_npixel_meta , &
       &                                      glc2000_class_fraction_meta , &
       &                                      ice_glc2000_meta , &
       &                                      z0_glc2000_meta , &
       &                                      root_glc2000_meta , &
       &                                      plcov_mx_glc2000_meta , &
       &                                      plcov_mn_glc2000_meta , &
       &                                      lai_mx_glc2000_meta , &
       &                                      lai_mn_glc2000_meta , &
       &                                      rs_min_glc2000_meta , &
       &                                      urban_glc2000_meta , &
       &                                      for_d_glc2000_meta , &
       &                                      for_e_glc2000_meta , &
       &                                      emissivity_glc2000_meta , &
       &                                      isa_field_meta  , &
       &                                      fr_land_glcc_meta  , &
       &                                      glcc_tot_npixel_meta , &
       &                                      glcc_class_npixel_meta , &
       &                                      glcc_class_fraction_meta , &
       &                                      ice_glcc_meta , &
       &                                      z0_glcc_meta , &
       &                                      root_glcc_meta , &
       &                                      plcov_mx_glcc_meta , &
       &                                      plcov_mn_glcc_meta , &
       &                                      lai_mx_glcc_meta , &
       &                                      lai_mn_glcc_meta , &
       &                                      rs_min_glcc_meta , &
       &                                      urban_glcc_meta , &
       &                                      for_d_glcc_meta , &
       &                                      for_e_glcc_meta , &
       &                                      emissivity_glcc_meta , &
       &                                      fr_land_ecosg_meta  , &
       &                                      ecosg_tot_npixel_meta , &
       &                                      ecosg_class_npixel_meta , &
       &                                      ecosg_class_fraction_meta , &
       &                                      ice_ecosg_meta , &
       &                                      z0_ecosg_meta , &
       &                                      root_ecosg_meta , &
       &                                      plcov_mx_ecosg_meta , &
       &                                      plcov_mn_ecosg_meta , &
       &                                      lai_mx_ecosg_meta , &
       &                                      lai_mn_ecosg_meta , &
       &                                      rs_min_ecosg_meta , &
       &                                      urban_ecosg_meta , &
       &                                      for_d_ecosg_meta , &
       &                                      for_e_ecosg_meta , &
       &                                      skinc_ecosg_meta , &
       &                                      emissivity_ecosg_meta , &
       &                                      fr_land_lu_meta  , &
       &                                      fr_land_mask_meta  , &
       &                                      lu_tot_npixel_meta , &
       &                                      lu_class_npixel_meta , &
       &                                      lu_class_fraction_meta , &
       &                                      ice_lu_meta , &
       &                                      z0_lu_meta , &
       &                                      root_lu_meta , &
       &                                      plcov_mx_lu_meta , &
       &                                      plcov_mn_lu_meta , &
       &                                      lai_mx_lu_meta , &
       &                                      lai_mn_lu_meta , &
       &                                      z012_lu_meta , &
       &                                      z012_tot_meta , &
       &                                      plcov12_lu_meta , &
       &                                      lai12_lu_meta , &
       &                                      rs_min_lu_meta , &
       &                                      urban_lu_meta , &
       &                                      for_d_lu_meta , &
       &                                      for_e_lu_meta , &
       &                                      skinc_lu_meta , &
       &                                      emissivity_lu_meta , &
       &                                      fr_ocean_lu_meta  , &
       &                                      hh_topo_meta      , &
       &                                      hh_topo_max_meta  , &
       &                                      hh_topo_min_meta  , &
       &                                      hh_fis_meta       , &
       &                                      fr_land_topo_meta , &
       &                                      stdh_topo_meta    , &
       &                                      theta_topo_meta   , &
       &                                      aniso_topo_meta   , &
       &                                      slope_topo_meta   , &
       &                                      hh_vert_meta      , &
       &                                      npixel_vert_meta  , &
       &                                      z0_topo_meta      , &
       &                                      slope_asp_topo_meta  , &
       &                                      slope_ang_topo_meta  , &
       &                                      horizon_topo_meta  , &
       &                                      skyview_topo_meta  , &
       &                                      sgsl_meta , &
       &                                      fr_land_soil_meta , &
       &                                      soiltype_fao_meta , &
       &                                      soiltype_hwsd_meta , &
       &                                      soiltype_FAO_deep_meta , &
       &                                      soiltype_HWSD_deep_meta , &
       &                                      HWSD_SAND_meta , &
       &                                      HWSD_SILT_meta , &
       &                                      HWSD_CLAY_meta , &
       &                                      HWSD_OC_meta , &
       &                                      HWSD_BD_meta , &
       &                                      HWSD_DM_meta , &
       &                                      HWSD_SAND_DEEP_meta , &
       &                                      HWSD_SILT_DEEP_meta , &
       &                                      HWSD_CLAY_DEEP_meta , &
       &                                      HWSD_OC_DEEP_meta , &
       &                                      HWSD_BD_DEEP_meta , &
       &                                      HWSD_DM_DEEP_meta , &
       &                                      lake_depth_meta , &
       &                                      fr_lake_meta , &
       &                                      flake_tot_npixel_meta

  TYPE(netcdf_grid_mapping)               :: nc_grid_def_cosmo, &
       &                                     nc_grid_def_icon

  CHARACTER (len=1), PARAMETER            :: c_undef = "-" !< default character for undefined string
  CHARACTER (len=80)                      :: aot_type_shortname(1:5) !< short names for optical thickness of aerosol types

  CONTAINS

  !> define buffer dimensions for netcdf output
  SUBROUTINE def_dimension_info_buffer(tg,nhori)

    TYPE(target_grid_def), INTENT(IN)      :: tg !< structure with target grid description
    INTEGER(KIND=i4), INTENT(IN), OPTIONAL :: nhori

    ! set meta information for strucutre dim_2d_tg
    dim_2d_tg(1)%dimname = 'ie'
    dim_2d_tg(1)%dimsize = tg%ie
    dim_2d_tg(2)%dimname = 'je'
    dim_2d_tg(2)%dimsize = tg%je

    ! set meta information for strucutre dim_3d_tg
    dim_3d_tg(1)%dimname = 'ie'
    dim_3d_tg(1)%dimsize = tg%ie
    dim_3d_tg(2)%dimname = 'je'
    dim_3d_tg(2)%dimsize = tg%je
    dim_3d_tg(3)%dimname = 'ke'
    dim_3d_tg(3)%dimsize = tg%ke

    IF(PRESENT(nhori)) THEN
      ! set meta information for strucutre dim_3d_tg
      dim_4d_tg(1)%dimname = 'ie'
      dim_4d_tg(1)%dimsize = tg%ie
      dim_4d_tg(2)%dimname = 'je'
      dim_4d_tg(2)%dimsize = tg%je
      dim_4d_tg(3)%dimname = 'ke'
      dim_4d_tg(3)%dimsize = tg%ke
      dim_4d_tg(4)%dimname = 'nhori'
      dim_4d_tg(4)%dimsize = nhori
    ENDIF

  END SUBROUTINE def_dimension_info_buffer


  !> define COSMO grid dimensions for netcdf output
  SUBROUTINE def_dimension_info_cosmo(cosmo_grid,nhori)
    TYPE(rotated_lonlat_grid), INTENT(IN) :: cosmo_grid !< structure which contains the definition of the COSMO grid
    INTEGER(KIND=i4), INTENT(IN), OPTIONAL :: nhori

    ! set meta information for strucutre dim_rlon_cosmo
    dim_rlon_cosmo(1)%dimname = 'rlon'
    dim_rlon_cosmo(1)%dimsize = cosmo_grid%nlon_rot

    ! set meta information for strucutre dim_rlat_cosmo
    dim_rlat_cosmo(1)%dimname = 'rlat'
    dim_rlat_cosmo(1)%dimsize = cosmo_grid%nlat_rot

    ! set meta information for strucutre dim_2d_cosmo
    dim_2d_cosmo(1) = dim_rlon_cosmo(1)
    dim_2d_cosmo(2) = dim_rlat_cosmo(1)

    ! set meta information for variable rlon

    rlon_meta%varname = 'rlon'
    rlon_meta%n_dim = 1
    rlon_meta%diminfo => dim_rlon_cosmo
    rlon_meta%vartype = vartype_real !REAL variable
    rlon_meta%standard_name = "grid_longitude" !_br 08.04.14
    rlon_meta%long_name =  "longitude in rotated pole grid"
    rlon_meta%shortName = c_undef
    rlon_meta%stepType = 'instant'
    rlon_meta%units =  "degrees"
    rlon_meta%grid_mapping = c_undef
    rlon_meta%coordinates = c_undef
    rlon_meta%data_set = c_undef

    ! set meta information for variable rlat

    rlat_meta%varname = 'rlat'
    rlat_meta%n_dim = 1
    rlat_meta%diminfo => dim_rlat_cosmo
    rlat_meta%vartype = vartype_real !REAL variable
    rlat_meta%standard_name = "grid_latitude" !_br 08.04.14
    rlat_meta%long_name =  "latitude in rotated pole grid"
    rlat_meta%shortName = c_undef
    rlat_meta%stepType = 'instant'
    rlat_meta%units =  "degrees"
    rlat_meta%grid_mapping = c_undef
    rlat_meta%coordinates = c_undef
    rlat_meta%data_set = c_undef

    IF (PRESENT(nhori)) THEN
      ! set meta information for strucutre dim_nhori_cosmo
      dim_nhori_cosmo(1)%dimname = 'nhori'
      dim_nhori_cosmo(1)%dimsize = nhori

      ! set meta information for strucutre dim_2d_cosmo
      dim_3d_cosmo(1) = dim_rlon_cosmo(1)
      dim_3d_cosmo(2) = dim_rlat_cosmo(1)
      dim_3d_cosmo(3) = dim_nhori_cosmo(1)

    ENDIF


  END SUBROUTINE def_dimension_info_cosmo


  !> define ICON grid dimensions for netcdf output
  SUBROUTINE def_dimension_info_icon(icon_grid)
    TYPE(icosahedral_triangular_grid), INTENT(IN) :: icon_grid !< structure which contains the definition of the ICON grid

    ! set meta information for strucutre dim_icon
    dim_icon(1)%dimname = 'cell'
    dim_icon(1)%dimsize = icon_grid%ncell ! total number of cells in ICON grid
    dim_icon(2)%dimname = 'vertex'
    dim_icon(2)%dimsize = icon_grid%nvertex ! total number of vertices in ICON grid
    dim_icon(3)%dimname = 'edge'
    dim_icon(3)%dimsize = icon_grid%nedge ! total number of edges in ICON grid

    dim_icon(4)%dimname = 'nc'
    dim_icon(4)%dimsize = icon_grid%ncells_per_edge ! number of cell per edge in ICON grid
    dim_icon(5)%dimname = 'nv'
    dim_icon(5)%dimsize = icon_grid%nvertex_per_cell ! number of vertices per cell
    dim_icon(6)%dimname = 'ne'
    dim_icon(6)%dimsize = icon_grid%nedges_per_vertex ! number of edges per vertex
    dim_icon(7)%dimname = 'no'
    dim_icon(7)%dimsize = icon_grid%nchilds_per_cell ! number of childs per cell

    dim_cells_icon(1)= dim_icon(1)

    dim_2d_icon(1) = dim_icon(5)
    dim_2d_icon(2) = dim_icon(1)



   ! set meta information for variable clon
    clon_meta%varname = 'clon'
    clon_meta%n_dim = 1
    clon_meta%diminfo => dim_cells_icon
    clon_meta%vartype = vartype_real !REAL variable
    clon_meta%standard_name = "grid_longitude" !_br 08.04.14
    clon_meta%long_name =  "longitude of icon grid cell centre"
    clon_meta%shortName = c_undef
    clon_meta%stepType = 'instant'
    clon_meta%units =  "radians"
    clon_meta%grid_mapping = c_undef
    clon_meta%coordinates = c_undef
    clon_meta%data_set = c_undef

    ! set meta information for variable clat

    clat_meta%varname = 'clat'
    clat_meta%n_dim = 1
    clat_meta%diminfo =>  dim_cells_icon
    clat_meta%vartype = vartype_real !REAL variable
    clat_meta%standard_name = "grid_latitude" !_br 08.04.14
    clat_meta%long_name =  "latitude of icon grid cell centre"
    clat_meta%shortName = c_undef
    clat_meta%stepType = 'instant'
    clat_meta%units =  "radians"
    clat_meta%grid_mapping = c_undef
    clat_meta%coordinates = c_undef
    clat_meta%data_set = c_undef

 ! set meta information for variable clon_vertices_meta

    clon_vertices_meta%varname = 'clon_vertices'
    clon_vertices_meta%n_dim = 2
    clon_vertices_meta%diminfo => dim_2d_icon
    clon_vertices_meta%vartype = vartype_real !REAL variable
    clon_vertices_meta%standard_name = c_undef !_br 08.04.14
    clon_vertices_meta%long_name =  "longitude of icon grid cell vertices"
    clon_vertices_meta%shortName = c_undef
    clon_vertices_meta%stepType = 'instant'
    clon_vertices_meta%units =  "radians"
    clon_vertices_meta%grid_mapping = c_undef
    clon_vertices_meta%coordinates = c_undef
    clon_vertices_meta%data_set = c_undef

    ! set meta information for variable clat_vertices_meta

    clat_vertices_meta%varname = 'clat_vertices'
    clat_vertices_meta%n_dim = 2
    clat_vertices_meta%diminfo => dim_2d_icon
    clat_vertices_meta%vartype = vartype_real !REAL variable
    clat_vertices_meta%standard_name = c_undef !_br 08.04.14
    clat_vertices_meta%long_name =  "latitude of icon grid cell vertices"
    clat_vertices_meta%shortName = c_undef
    clat_vertices_meta%stepType = 'instant'
    clat_vertices_meta%units =  "radians"
    clat_vertices_meta%grid_mapping = c_undef
    clat_vertices_meta%coordinates = c_undef
    clat_vertices_meta%data_set = c_undef


  END SUBROUTINE def_dimension_info_icon


  !> define meta information for variable crutemp for netcdf output
  SUBROUTINE def_crutemp_meta(diminfo,coordinates,grid_mapping)
    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)
    !ha change varname to CAPITAL letters
    crutemp_meta%varname = 'T_CL'
    crutemp_meta%n_dim = n_dim
    crutemp_meta%diminfo => diminfo
    crutemp_meta%vartype = vartype_real !REAL variable
    crutemp_meta%standard_name = 'soil_temperature'
    crutemp_meta%long_name = 'CRU near surface temperature climatology'
    crutemp_meta%shortName = 'T_2M_CL'
    crutemp_meta%stepType = 'avg'
    crutemp_meta%units = 'K'
    crutemp_meta%grid_mapping = gridmp
    crutemp_meta%coordinates = coord
    crutemp_meta%data_set = 'CRU'

  END SUBROUTINE def_crutemp_meta

 !> define meta information for variable crutemp for netcdf output
  SUBROUTINE def_cruelev_meta(diminfo,coordinates,grid_mapping)
    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)
    !ha change varname to CAPITAL letters
    cruelev_meta%varname = 'HSURF'
    cruelev_meta%n_dim = n_dim
    cruelev_meta%diminfo => diminfo
    cruelev_meta%vartype = vartype_real !REAL variable
    cruelev_meta%standard_name = 'surface_altitude'
    cruelev_meta%long_name = 'CRU grid elevation'
    cruelev_meta%shortName = 'T_2M_CL'
    cruelev_meta%stepType = 'instant'
    cruelev_meta%units = 'm'
    cruelev_meta%grid_mapping = gridmp
    cruelev_meta%coordinates = coord
    cruelev_meta%data_set = 'CRU'

  END SUBROUTINE def_cruelev_meta

  !> define meta information for soil data for netcdf output
  SUBROUTINE def_soil_meta(diminfo,isoil_data,coordinates,grid_mapping)

  USE mo_soil_data,       ONLY: FAO_data, HWSD_data, HWSD_map

    TYPE(dim_meta_info), TARGET  :: diminfo(:)     !< pointer to dimensions of variable
    INTEGER (KIND=i4), INTENT(IN):: isoil_data
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)


    fr_land_soil_meta%varname = 'FR_LAND_SOIL'
    fr_land_soil_meta%n_dim = n_dim
    fr_land_soil_meta%diminfo => diminfo
    fr_land_soil_meta%vartype = vartype_real !REAL variable
    fr_land_soil_meta%standard_name = 'land_area_fraction' !_br 08.04.14
    SELECT CASE (isoil_data)
    CASE(FAO_data)
      fr_land_soil_meta%long_name = 'Fraction Land due to FAO Digital Soil Map of the World'
      fr_land_soil_meta%data_set  = 'FAO Digital Soil Map of the World'
    CASE(HWSD_data,HWSD_map)
      fr_land_soil_meta%long_name = 'Fraction Land due to HWSD Digital Soil Map of the World'
      fr_land_soil_meta%data_set = 'HWSD Digital Soil Map of the World'
    END SELECT
    fr_land_soil_meta%shortName = 'FR_LAND'
    fr_land_soil_meta%stepType = 'instant'
    fr_land_soil_meta%units = c_undef
    fr_land_soil_meta%grid_mapping = gridmp
    fr_land_soil_meta%coordinates = coord


    soiltype_fao_meta%varname = 'SOILTYP'
    soiltype_fao_meta%n_dim = n_dim
    soiltype_fao_meta%diminfo => diminfo
    soiltype_fao_meta%vartype = vartype_int !REAL variable
    soiltype_fao_meta%standard_name = 'soil_type'
    SELECT CASE (isoil_data)
    CASE(FAO_data)
      soiltype_fao_meta%long_name = 'soil type derived from FAO Digital Soil Map of the World'
      soiltype_fao_meta%data_set = 'FAO Digital Soil Map of the World'
    CASE(HWSD_data,HWSD_map)
      soiltype_fao_meta%long_name = 'soil type derived from HWSD Digital Soil Map of the World'
      soiltype_fao_meta%data_set = 'HWSD Digital Soil Map of the World'
    END SELECT
    soiltype_fao_meta%shortName = 'SOILTYP'
    soiltype_fao_meta%stepType = 'instant'
    soiltype_fao_meta%units = c_undef
    soiltype_fao_meta%grid_mapping = gridmp
    soiltype_fao_meta%coordinates = coord

    soiltype_hwsd_meta%varname = 'HWSDTYP'
    soiltype_hwsd_meta%n_dim = n_dim
    soiltype_hwsd_meta%diminfo => diminfo
    soiltype_hwsd_meta%vartype = vartype_int !REAL variable
    soiltype_hwsd_meta%standard_name = 'HWSDTYP'
    SELECT CASE (isoil_data)
    CASE(FAO_data)
      soiltype_hwsd_meta%long_name = 'soil id derived from HWSD'
    CASE(HWSD_data)
      soiltype_hwsd_meta%long_name = 'soil id derived from HWSD '
    END SELECT
    soiltype_hwsd_meta%shortName = 'HWSDTYP'
    soiltype_hwsd_meta%units = c_undef
    soiltype_hwsd_meta%grid_mapping = gridmp
    soiltype_hwsd_meta%coordinates = coord

    soiltype_FAO_deep_meta%varname = 'SUBSOILTYP_FAO'
    soiltype_FAO_deep_meta%n_dim = n_dim
    soiltype_FAO_deep_meta%diminfo => diminfo
    soiltype_FAO_deep_meta%vartype = vartype_int !REAL variable
    soiltype_FAO_deep_meta%standard_name = 'subsoil_type'
    soiltype_FAO_deep_meta%long_name = 'deep soil type from HWSD similar to FAO'
    soiltype_FAO_deep_meta%shortName = 'SUBSOILTYP_FAO'
    soiltype_FAO_deep_meta%units = c_undef
    soiltype_FAO_deep_meta%grid_mapping = gridmp
    soiltype_FAO_deep_meta%coordinates = coord

    soiltype_HWSD_deep_meta%varname = 'SUBSOILTYP_HWSD'
    soiltype_HWSD_deep_meta%n_dim = n_dim
    soiltype_HWSD_deep_meta%diminfo => diminfo
    soiltype_HWSD_deep_meta%vartype = vartype_int !REAL variable
    soiltype_HWSD_deep_meta%standard_name = 'subsoil_type'
    soiltype_HWSD_deep_meta%long_name = 'deep soil index derived from HWSD'
    soiltype_HWSD_deep_meta%shortName = 'SUBSOILTYP_HWSD'
    soiltype_HWSD_deep_meta%units = c_undef
    soiltype_HWSD_deep_meta%grid_mapping = gridmp
    soiltype_HWSD_deep_meta%coordinates = coord

    IF (isoil_data == HWSD_data) THEN
      HWSD_SAND_meta%varname = 'FR_SAND'
      HWSD_SAND_meta%n_dim = n_dim
      HWSD_SAND_meta%diminfo => diminfo
      HWSD_SAND_meta%vartype = vartype_real !REAL variable
      HWSD_SAND_meta%standard_name = c_undef !_br 08.04.14
      HWSD_SAND_meta%long_name = 'fraction of sand for soil index'
      HWSD_SAND_meta%shortName = 'ru-103d'
      HWSD_SAND_meta%stepType = 'instant'
      HWSD_SAND_meta%units = c_undef
      HWSD_SAND_meta%grid_mapping = gridmp
      HWSD_SAND_meta%coordinates = coord
      HWSD_SAND_meta%data_set = 'HWSD Digital Soil Map of the World'

      HWSD_SILT_meta%varname = 'FR_SILT'
      HWSD_SILT_meta%n_dim = n_dim
      HWSD_SILT_meta%diminfo => diminfo
      HWSD_SILT_meta%vartype = vartype_real !REAL variable
      HWSD_SILT_meta%standard_name = c_undef !_br 08.04.14
      HWSD_SILT_meta%long_name = 'fraction of silt for soil index'
      HWSD_SILT_meta%shortName = 'ru-103w'
      HWSD_SILT_meta%stepType = 'instant'
      HWSD_SILT_meta%units = c_undef
      HWSD_SILT_meta%grid_mapping = gridmp
      HWSD_SILT_meta%coordinates = coord
      HWSD_SILT_meta%data_set = 'HWSD Digital Soil Map of the World'

      HWSD_CLAY_meta%varname = 'FR_CLAY'
      HWSD_CLAY_meta%n_dim = n_dim
      HWSD_CLAY_meta%diminfo => diminfo
      HWSD_CLAY_meta%vartype = vartype_real !REAL variable
      HWSD_CLAY_meta%standard_name = c_undef !_br 08.04.14
      HWSD_CLAY_meta%long_name = 'fraction of clay for soil index'
      HWSD_CLAY_meta%shortName = 'sr-90d'
      HWSD_CLAY_meta%stepType = 'instant'
      HWSD_CLAY_meta%units = c_undef
      HWSD_CLAY_meta%grid_mapping = gridmp
      HWSD_CLAY_meta%coordinates = coord
      HWSD_CLAY_meta%data_set = 'HWSD Digital Soil Map of the World'

      HWSD_OC_meta%varname = 'FR_OC'
      HWSD_OC_meta%n_dim = n_dim
      HWSD_OC_meta%diminfo => diminfo
      HWSD_OC_meta%vartype = vartype_real !REAL variable
      HWSD_OC_meta%standard_name = c_undef !_br 08.04.14
      HWSD_OC_meta%long_name = 'fraction of oc for soil index'
      HWSD_OC_meta%shortName = 'sr-90w'
      HWSD_OC_meta%stepType = 'instant'
      HWSD_OC_meta%units = c_undef
      HWSD_OC_meta%grid_mapping = gridmp
      HWSD_OC_meta%coordinates = coord
      HWSD_OC_meta%data_set = 'HWSD Digital Soil Map of the World'

      HWSD_BD_meta%varname = 'BULK_DENS'
      HWSD_BD_meta%n_dim = n_dim
      HWSD_BD_meta%diminfo => diminfo
      HWSD_BD_meta%vartype = vartype_real !REAL variable
      HWSD_BD_meta%standard_name = c_undef !_br 08.04.14
      HWSD_BD_meta%long_name = 'bulk density for soil index'
      HWSD_BD_meta%shortName = 'i-131ad'
      HWSD_BD_meta%stepType = 'instant'
      HWSD_BD_meta%units = c_undef
      HWSD_BD_meta%grid_mapping = gridmp
      HWSD_BD_meta%coordinates = coord
      HWSD_BD_meta%data_set = 'HWSD Digital Soil Map of the World'

      HWSD_DM_meta%varname = 'DUMMY'
      HWSD_DM_meta%n_dim = n_dim
      HWSD_DM_meta%diminfo => diminfo
      HWSD_DM_meta%vartype = vartype_real !REAL variable
      HWSD_DM_meta%standard_name = c_undef !_br 08.04.14
      HWSD_DM_meta%long_name = 'bulk density for soil index'
      HWSD_DM_meta%shortName = 'i-131aw'
      HWSD_DM_meta%stepType = 'instant'
      HWSD_DM_meta%units = c_undef
      HWSD_DM_meta%grid_mapping = gridmp
      HWSD_DM_meta%coordinates = coord
      HWSD_DM_meta%data_set = 'HWSD Digital Soil Map of the World'

      HWSD_SAND_DEEP_meta%varname = 'SUB_FR_SAND'
      HWSD_SAND_DEEP_meta%n_dim = n_dim
      HWSD_SAND_DEEP_meta%diminfo => diminfo
      HWSD_SAND_DEEP_meta%vartype = vartype_real !REAL variable
      HWSD_SAND_DEEP_meta%standard_name = c_undef !_br 08.04.14
      HWSD_SAND_DEEP_meta%long_name = 'fraction of sand for deep soil index'
      HWSD_SAND_DEEP_meta%shortName = 'ru-103d'
      HWSD_SAND_DEEP_meta%stepType = 'instant'
      HWSD_SAND_DEEP_meta%units = c_undef
      HWSD_SAND_DEEP_meta%grid_mapping = gridmp
      HWSD_SAND_DEEP_meta%coordinates = coord
      HWSD_SAND_DEEP_meta%data_set = 'HWSD Digital Soil Map of the World'

      HWSD_SILT_DEEP_meta%varname = 'SUB_FR_SILT'
      HWSD_SILT_DEEP_meta%n_dim = n_dim
      HWSD_SILT_DEEP_meta%diminfo => diminfo
      HWSD_SILT_DEEP_meta%vartype = vartype_real !REAL variable
      HWSD_SILT_DEEP_meta%standard_name = c_undef !_br 08.04.14
      HWSD_SILT_DEEP_meta%long_name = 'fraction of silt for deep soil index'
      HWSD_SILT_DEEP_meta%shortName = 'ru-103w'
      HWSD_SILT_DEEP_meta%stepType = 'instant'
      HWSD_SILT_DEEP_meta%units = c_undef
      HWSD_SILT_DEEP_meta%grid_mapping = gridmp
      HWSD_SILT_DEEP_meta%coordinates = coord
      HWSD_SILT_DEEP_meta%data_set = 'HWSD Digital Soil Map of the World'

      HWSD_CLAY_DEEP_meta%varname = 'SUB_FR_CLAY'
      HWSD_CLAY_DEEP_meta%n_dim = n_dim
      HWSD_CLAY_DEEP_meta%diminfo => diminfo
      HWSD_CLAY_DEEP_meta%vartype = vartype_real !REAL variable
      HWSD_CLAY_DEEP_meta%standard_name = c_undef !_br 08.04.14
      HWSD_CLAY_DEEP_meta%long_name = 'fraction of clay for deep soil index'
      HWSD_CLAY_DEEP_meta%shortName = 'sr-90d'
      HWSD_CLAY_DEEP_meta%stepType = 'instant'
      HWSD_CLAY_DEEP_meta%units = c_undef
      HWSD_CLAY_DEEP_meta%grid_mapping = gridmp
      HWSD_CLAY_DEEP_meta%coordinates = coord
      HWSD_CLAY_DEEP_meta%data_set = 'HWSD Digital Soil Map of the World'

      HWSD_OC_DEEP_meta%varname = 'SUB_FR_OC'
      HWSD_OC_DEEP_meta%n_dim = n_dim
      HWSD_OC_DEEP_meta%diminfo => diminfo
      HWSD_OC_DEEP_meta%vartype = vartype_real !REAL variable
      HWSD_OC_DEEP_meta%standard_name = c_undef !_br 08.04.14
      HWSD_OC_DEEP_meta%long_name = 'fraction of oc for deep soil index'
      HWSD_OC_DEEP_meta%shortName = 'sr-90w'
      HWSD_OC_DEEP_meta%stepType = 'instant'
      HWSD_OC_DEEP_meta%units = c_undef
      HWSD_OC_DEEP_meta%grid_mapping = gridmp
      HWSD_OC_DEEP_meta%coordinates = coord
      HWSD_OC_DEEP_meta%data_set = 'HWSD Digital Soil Map of the World'

      HWSD_BD_DEEP_meta%varname = 'SUB_BULK_DENS'
      HWSD_BD_DEEP_meta%n_dim = n_dim
      HWSD_BD_DEEP_meta%diminfo => diminfo
      HWSD_BD_DEEP_meta%vartype = vartype_real !REAL variable
      HWSD_BD_DEEP_meta%standard_name = c_undef !_br 08.04.14
      HWSD_BD_DEEP_meta%long_name = 'bulk density for deep soil index'
      HWSD_BD_DEEP_meta%shortName = 'i-131ad'
      HWSD_BD_DEEP_meta%stepType = 'instant'
      HWSD_BD_DEEP_meta%units = c_undef
      HWSD_BD_DEEP_meta%grid_mapping = gridmp
      HWSD_BD_DEEP_meta%coordinates = coord
      HWSD_BD_DEEP_meta%data_set = 'HWSD Digital Soil Map of the World'

      HWSD_DM_DEEP_meta%varname = 'SUB_DUMMY'
      HWSD_DM_DEEP_meta%n_dim = n_dim
      HWSD_DM_DEEP_meta%diminfo => diminfo
      HWSD_DM_DEEP_meta%vartype = vartype_real !REAL variable
      HWSD_DM_DEEP_meta%standard_name = c_undef !_br 08.04.14
      HWSD_DM_DEEP_meta%long_name = 'bulk density for deep soil index'
      HWSD_DM_DEEP_meta%shortName = 'i-131aw'
      HWSD_DM_DEEP_meta%stepType = 'instant'
      HWSD_DM_DEEP_meta%units = c_undef
      HWSD_DM_DEEP_meta%grid_mapping = gridmp
      HWSD_DM_DEEP_meta%coordinates = coord
      HWSD_DM_DEEP_meta%data_set = 'HWSD Digital Soil Map of the World'
    END IF


  END SUBROUTINE def_soil_meta

  SUBROUTINE def_alb_meta(ntime,diminfo,coordinates,grid_mapping)

    INTEGER (KIND=i4), INTENT(IN) :: ntime !< number of times
    TYPE(dim_meta_info),TARGET    :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL  :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL  :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER                       :: n_dim      !< number of dimensions
    CHARACTER (len=80)            :: gridmp
    CHARACTER (len=80)            :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)

    IF (ALLOCATED(dim_alb_tg)) DEALLOCATE(dim_alb_tg)
    ALLOCATE(dim_alb_tg(1:n_dim+1))
    SELECT CASE(n_dim)

    CASE (1)
      dim_alb_tg(1)%dimname = diminfo(1)%dimname
      dim_alb_tg(1)%dimsize = diminfo(1)%dimsize
      dim_alb_tg(2)%dimname = 'time'
      dim_alb_tg(2)%dimsize = ntime
    CASE (2)
      dim_alb_tg(1)%dimname = diminfo(1)%dimname
      dim_alb_tg(1)%dimsize = diminfo(1)%dimsize
      dim_alb_tg(2)%dimname = diminfo(2)%dimname
      dim_alb_tg(2)%dimsize = diminfo(2)%dimsize
      dim_alb_tg(3)%dimname = 'time'
      dim_alb_tg(3)%dimsize = ntime
    CASE (3)
      dim_alb_tg(1)%dimname = diminfo(1)%dimname
      dim_alb_tg(1)%dimsize = diminfo(1)%dimsize
      dim_alb_tg(2)%dimname = diminfo(2)%dimname
      dim_alb_tg(2)%dimsize = diminfo(2)%dimsize
      dim_alb_tg(3)%dimname = diminfo(3)%dimname
      dim_alb_tg(3)%dimsize = diminfo(3)%dimsize
      dim_alb_tg(4)%dimname = 'time'
      dim_alb_tg(4)%dimsize = ntime
    END SELECT

    alb_field_mom_meta%varname = 'ALB_DIF12'
    alb_field_mom_meta%n_dim = n_dim + 1
    alb_field_mom_meta%diminfo => dim_alb_tg
    alb_field_mom_meta%vartype = vartype_real !REAL variable
    alb_field_mom_meta%standard_name = 'surface_albedo' !_br 08.04.14
    alb_field_mom_meta%long_name = 'Albedo'
    alb_field_mom_meta%shortName = 'ALB_DIF12'
    alb_field_mom_meta%units = '%'
    alb_field_mom_meta%grid_mapping = gridmp
    alb_field_mom_meta%coordinates = coord
    alb_field_mom_meta%data_set = 'MODIS'
    alb_field_mom_meta%stepType = 'avg'

    alnid_field_mom_meta%varname = 'ALNID12'
    alnid_field_mom_meta%n_dim = n_dim + 1
    alnid_field_mom_meta%diminfo => dim_alb_tg
    alnid_field_mom_meta%vartype = vartype_real !REAL variable
    alnid_field_mom_meta%standard_name = c_undef !_br 08.04.14
    alnid_field_mom_meta%long_name = 'NI_Albedo'
    alnid_field_mom_meta%shortName = 'ALB_DIF12'
    alnid_field_mom_meta%units = '%'
    alnid_field_mom_meta%grid_mapping = gridmp
    alnid_field_mom_meta%coordinates = coord
    alnid_field_mom_meta%data_set = 'MODIS'
    alnid_field_mom_meta%stepType = 'avg'

    aluvd_field_mom_meta%varname = 'ALUVD12'
    aluvd_field_mom_meta%n_dim = n_dim + 1
    aluvd_field_mom_meta%diminfo => dim_alb_tg
    aluvd_field_mom_meta%vartype = vartype_real !REAL variable
    aluvd_field_mom_meta%standard_name = c_undef !_br 08.04.14
    aluvd_field_mom_meta%long_name = 'UV_Albedo'
    aluvd_field_mom_meta%shortName = 'ALB_DIF12'
    aluvd_field_mom_meta%units = '%'
    aluvd_field_mom_meta%grid_mapping = gridmp
    aluvd_field_mom_meta%coordinates = coord
    aluvd_field_mom_meta%data_set = 'MODIS'
    aluvd_field_mom_meta%stepType = 'avg'

    alb_interpol_meta%varname = 'ALB_I'
    alb_interpol_meta%n_dim = n_dim + 1
    alb_interpol_meta%diminfo => dim_alb_tg
    alb_interpol_meta%vartype = vartype_real !REAL variable
    alb_interpol_meta%standard_name = c_undef !_br 08.04.14
    alb_interpol_meta%long_name = 'Interpolated Albedo'
    alb_interpol_meta%shortName = 'ALB_RAD'
    alb_interpol_meta%units = c_undef
    alb_interpol_meta%grid_mapping = gridmp
    alb_interpol_meta%coordinates = coord
    alb_interpol_meta%data_set = 'MODIS'

    alb_dry_meta%varname = 'ALB_DRY'
    alb_dry_meta%n_dim = n_dim
    alb_dry_meta%diminfo => diminfo
    alb_dry_meta%vartype = vartype_real !REAL variable
    alb_dry_meta%standard_name = 'surface_albedo' !_br 08.04.14
    alb_dry_meta%long_name = 'soil albedo for dry soil'
    alb_dry_meta%shortName = 'ALB_DRY'
    alb_dry_meta%units = '1'
    alb_dry_meta%grid_mapping = gridmp
    alb_dry_meta%coordinates = coord
    alb_dry_meta%data_set = 'MODIS soil color derived soil albedo'

    alb_sat_meta%varname = 'ALB_SAT'
    alb_sat_meta%n_dim = n_dim
    alb_sat_meta%diminfo => diminfo
    alb_sat_meta%vartype = vartype_real !REAL variable
    alb_sat_meta%standard_name = 'surface_albedo' !_br 08.04.14
    alb_sat_meta%long_name = 'soil albedo for saturated soil'
    alb_sat_meta%shortName = 'ALB_SAT'
    alb_sat_meta%units = '1'
    alb_sat_meta%grid_mapping = gridmp
    alb_sat_meta%coordinates = coord
    alb_sat_meta%data_set = 'MODIS soil color derived soil albedo'

  END SUBROUTINE def_alb_meta

  !> define meta information for AHF data for netcdf output
  SUBROUTINE def_ahf_meta(diminfo,coordinates,grid_mapping)

    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord
    CHARACTER (len=80) :: dataset     ! info dataset !_br 14.04.16

    gridmp = c_undef
    coord = c_undef
    dataset = c_undef  !_br 14.04.16

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)

    IF (iahf_type == 1 ) THEN
      dataset = "For 2006 after Flanner(2009) 2,5'"
    ELSE IF  (iahf_type == 2 ) THEN
      dataset = 'For 2006 after Flanner(2009) 30"'
    ENDIF

    ahf_field_meta%varname = 'AHF'
    ahf_field_meta%n_dim = n_dim
    ahf_field_meta%diminfo => diminfo
    ahf_field_meta%vartype = vartype_real !REAL variable
    ahf_field_meta%standard_name = c_undef !_br 14.04.16
    ahf_field_meta%long_name = 'Anthropogenic heat flux' !_br 14.04.16
    ahf_field_meta%shortName = 'GFLUX'
    ahf_field_meta%units = 'W m-2' !_br 14.04.16
    ahf_field_meta%grid_mapping = gridmp
    ahf_field_meta%coordinates = coord
    ahf_field_meta%data_set = dataset !_br 14.04.16

  END SUBROUTINE def_ahf_meta

  !> define meta information for  landuse target fields
  SUBROUTINE def_isa_fields_meta(diminfo,coordinates,grid_mapping)


    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping


    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord
    CHARACTER (len=80) :: dataset     ! info dataset !_br 14.04.16

    gridmp = c_undef
    coord = c_undef
    dataset = c_undef  !_br 14.04.16

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)

    n_dim = SIZE(diminfo)

    IF (isa_type == 1 ) THEN
      dataset = 'NOAA 30"'
    ELSE IF  (isa_type == 2 ) THEN
      dataset = 'European Environmental Agency 10"'
    ENDIF

    ! urban_isa_meta
    !isa_field_meta%varname = 'FR_PAVED'
    isa_field_meta%varname = 'ISA'
    isa_field_meta%n_dim = n_dim
    isa_field_meta%diminfo => diminfo
    isa_field_meta%vartype = vartype_real !REAL variable
    isa_field_meta%standard_name = c_undef !_br 14.04.16
    !isa_field_meta%long_name = 'Fraction of impervious surface area'
    isa_field_meta%long_name = 'impervious surface area'
    !isa_field_meta%shortName = 'FR_PAVED' ! dummy for GRIB2
    isa_field_meta%shortName = 'ISA' ! dummy for GRIB2
    isa_field_meta%units =  c_undef
    isa_field_meta%grid_mapping = gridmp
    isa_field_meta%coordinates = coord
    isa_field_meta%data_set = dataset !_br 14.04.16

  END SUBROUTINE def_isa_fields_meta

  !> define meta information for NDVI data for netcdf output
  SUBROUTINE def_ndvi_meta(ntime,diminfo,coordinates,grid_mapping)
    INTEGER (KIND=i4), INTENT(IN) :: ntime !< number of times
    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)

    ! set meta information for strucutre dim_ndvi_tg
    IF (ALLOCATED(dim_ndvi_tg)) DEALLOCATE(dim_ndvi_tg)
    ALLOCATE(dim_ndvi_tg(1:n_dim+1))
    SELECT CASE(n_dim)
      CASE (1)
      dim_ndvi_tg(1)%dimname = diminfo(1)%dimname
      dim_ndvi_tg(1)%dimsize = diminfo(1)%dimsize
      dim_ndvi_tg(2)%dimname = 'time'
      dim_ndvi_tg(2)%dimsize = ntime
    CASE (2)
      dim_ndvi_tg(1)%dimname = diminfo(1)%dimname
      dim_ndvi_tg(1)%dimsize = diminfo(1)%dimsize
      dim_ndvi_tg(2)%dimname = diminfo(2)%dimname
      dim_ndvi_tg(2)%dimsize = diminfo(2)%dimsize
      dim_ndvi_tg(3)%dimname = 'time'
      dim_ndvi_tg(3)%dimsize = ntime
    CASE (3)
      dim_ndvi_tg(1)%dimname = diminfo(1)%dimname
      dim_ndvi_tg(1)%dimsize = diminfo(1)%dimsize
      dim_ndvi_tg(2)%dimname = diminfo(2)%dimname
      dim_ndvi_tg(2)%dimsize = diminfo(2)%dimsize
      dim_ndvi_tg(3)%dimname = diminfo(3)%dimname
      dim_ndvi_tg(3)%dimsize = diminfo(3)%dimsize
      dim_ndvi_tg(4)%dimname = 'time'
      dim_ndvi_tg(4)%dimsize = ntime
    END SELECT


    ndvi_max_meta%varname = 'NDVI_MAX'
    ndvi_max_meta%n_dim = n_dim
    ndvi_max_meta%diminfo => diminfo
    ndvi_max_meta%vartype = vartype_real !REAL variable
    ndvi_max_meta%standard_name = c_undef !_br 08.04.14
    ndvi_max_meta%long_name = 'NDVI yearly maximum for climatology 1998-2003'
    ndvi_max_meta%shortName = 'NDVI_MAX'
    ndvi_max_meta%stepType = 'max'
    ndvi_max_meta%units = c_undef
    ndvi_max_meta%grid_mapping = gridmp
    ndvi_max_meta%coordinates = coord
    ndvi_max_meta%data_set = 'NASA/GSFS climatology 1998-2003'

    ndvi_field_mom_meta%varname = 'NDVI'
    ndvi_field_mom_meta%n_dim = n_dim + 1
    ndvi_field_mom_meta%diminfo => dim_ndvi_tg
    ndvi_field_mom_meta%vartype = vartype_real !REAL variable
    ndvi_field_mom_meta%standard_name = c_undef !_br 08.04.14
    ndvi_field_mom_meta%long_name = 'monthly mean NDVI climatology 1998-2003'
    ndvi_field_mom_meta%shortName = 'NDVI'
    ndvi_field_mom_meta%stepType = 'avg'
    ndvi_field_mom_meta%units = c_undef
    ndvi_field_mom_meta%grid_mapping = gridmp
    ndvi_field_mom_meta%coordinates = coord
    ndvi_field_mom_meta%data_set = ' NASA/GSFS climatology 1998-2003'

    ndvi_ratio_mom_meta%varname = 'NDVI_MRAT'
    ndvi_ratio_mom_meta%n_dim = n_dim + 1
    ndvi_ratio_mom_meta%diminfo => dim_ndvi_tg
    ndvi_ratio_mom_meta%vartype = vartype_real !REAL variable
    ndvi_ratio_mom_meta%standard_name = c_undef !_br 08.04.14
    ndvi_ratio_mom_meta%long_name = '(monthly) proportion of actual value/maximum normalized differential vegetation index'
    ndvi_ratio_mom_meta%shortName = 'NDVI_MRAT'
    ndvi_ratio_mom_meta%stepType = 'avg'
    ndvi_ratio_mom_meta%units = c_undef
    ndvi_ratio_mom_meta%grid_mapping = gridmp
    ndvi_ratio_mom_meta%coordinates = coord
    ndvi_ratio_mom_meta%data_set = ' NASA/GSFS climatology 1998-2003'


  END SUBROUTINE def_ndvi_meta


  !> define meta information for EMISS data for netcdf output
  SUBROUTINE def_emiss_meta(ntime,diminfo,coordinates,grid_mapping)
    INTEGER (KIND=i4), INTENT(IN) :: ntime !< number of times
    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)

    ! set meta information for strucutre dim_emiss_tg
    IF (ALLOCATED(dim_emiss_tg)) DEALLOCATE(dim_emiss_tg)
    ALLOCATE(dim_emiss_tg(1:n_dim+1))
    SELECT CASE(n_dim)
      CASE (1)
      dim_emiss_tg(1)%dimname = diminfo(1)%dimname
      dim_emiss_tg(1)%dimsize = diminfo(1)%dimsize
      dim_emiss_tg(2)%dimname = 'time'
      dim_emiss_tg(2)%dimsize = ntime
    CASE (2)
      dim_emiss_tg(1)%dimname = diminfo(1)%dimname
      dim_emiss_tg(1)%dimsize = diminfo(1)%dimsize
      dim_emiss_tg(2)%dimname = diminfo(2)%dimname
      dim_emiss_tg(2)%dimsize = diminfo(2)%dimsize
      dim_emiss_tg(3)%dimname = 'time'
      dim_emiss_tg(3)%dimsize = ntime
    CASE (3)
      dim_emiss_tg(1)%dimname = diminfo(1)%dimname
      dim_emiss_tg(1)%dimsize = diminfo(1)%dimsize
      dim_emiss_tg(2)%dimname = diminfo(2)%dimname
      dim_emiss_tg(2)%dimsize = diminfo(2)%dimsize
      dim_emiss_tg(3)%dimname = diminfo(3)%dimname
      dim_emiss_tg(3)%dimsize = diminfo(3)%dimsize
      dim_emiss_tg(4)%dimname = 'time'
      dim_emiss_tg(4)%dimsize = ntime
    END SELECT


    emiss_max_meta%varname = 'EMISS_MAX'
    emiss_max_meta%n_dim = n_dim
    emiss_max_meta%diminfo => diminfo
    emiss_max_meta%vartype = vartype_real !REAL variable
    emiss_max_meta%standard_name = c_undef !_br 08.04.14
    emiss_max_meta%long_name = 'EMISS yearly maximum for climatology 1998-2003'
    emiss_max_meta%shortName = 'EMISS_MAX'
    emiss_max_meta%stepType = 'max'
    emiss_max_meta%units = c_undef
    emiss_max_meta%grid_mapping = gridmp
    emiss_max_meta%coordinates = coord
    emiss_max_meta%data_set = 'NASA/GSFS climatology 1998-2003'

    emiss_field_mom_meta%varname = 'EMISS'
    emiss_field_mom_meta%n_dim = n_dim + 1
    emiss_field_mom_meta%diminfo => dim_emiss_tg
    emiss_field_mom_meta%vartype = vartype_real !REAL variable
    emiss_field_mom_meta%standard_name = c_undef !_br 08.04.14
    emiss_field_mom_meta%long_name = 'monthly mean EMISS climatology 1998-2003'
    emiss_field_mom_meta%shortName = 'EMISS'
    emiss_field_mom_meta%stepType = 'avg'
    emiss_field_mom_meta%units = c_undef
    emiss_field_mom_meta%grid_mapping = gridmp
    emiss_field_mom_meta%coordinates = coord
    emiss_field_mom_meta%data_set = ' NASA/GSFS climatology 1998-2003'

    emiss_ratio_mom_meta%varname = 'EMISS_MRAT'
    emiss_ratio_mom_meta%n_dim = n_dim + 1
    emiss_ratio_mom_meta%diminfo => dim_emiss_tg
    emiss_ratio_mom_meta%vartype = vartype_real !REAL variable
    emiss_ratio_mom_meta%standard_name = c_undef !_br 08.04.14
    emiss_ratio_mom_meta%long_name = '(monthly) proportion of actual value/maximum normalized differential vegetation index'
    emiss_ratio_mom_meta%shortName = 'EMISS_MRAT'
    emiss_ratio_mom_meta%stepType = 'avg'
    emiss_ratio_mom_meta%units = c_undef
    emiss_ratio_mom_meta%grid_mapping = gridmp
    emiss_ratio_mom_meta%coordinates = coord
    emiss_ratio_mom_meta%data_set = ' NASA/GSFS climatology 1998-2003'


  END SUBROUTINE def_emiss_meta


  !> define meta information for SST data for netcdf output
  SUBROUTINE def_era_meta(ntime,diminfo,coordinates,grid_mapping)
    INTEGER (KIND=i4), INTENT(IN) :: ntime !< number of times
    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)

    ! set meta information for strucutre dim_era_tg
    IF (ALLOCATED(dim_era_tg)) DEALLOCATE(dim_era_tg)
    ALLOCATE(dim_era_tg(1:n_dim+1))
    SELECT CASE(n_dim)
      CASE (1)
      dim_era_tg(1)%dimname = diminfo(1)%dimname
      dim_era_tg(1)%dimsize = diminfo(1)%dimsize
      dim_era_tg(2)%dimname = 'time'
      dim_era_tg(2)%dimsize = ntime
    CASE (2)
      dim_era_tg(1)%dimname = diminfo(1)%dimname
      dim_era_tg(1)%dimsize = diminfo(1)%dimsize
      dim_era_tg(2)%dimname = diminfo(2)%dimname
      dim_era_tg(2)%dimsize = diminfo(2)%dimsize
      dim_era_tg(3)%dimname = 'time'
      dim_era_tg(3)%dimsize = ntime
    CASE (3)
      dim_era_tg(1)%dimname = diminfo(1)%dimname
      dim_era_tg(1)%dimsize = diminfo(1)%dimsize
      dim_era_tg(2)%dimname = diminfo(2)%dimname
      dim_era_tg(2)%dimsize = diminfo(2)%dimsize
      dim_era_tg(3)%dimname = diminfo(3)%dimname
      dim_era_tg(3)%dimsize = diminfo(3)%dimsize
      dim_era_tg(4)%dimname = 'time'
      dim_era_tg(4)%dimsize = ntime
    END SELECT



    sst_field_meta%varname = 'T_SEA'
    sst_field_meta%n_dim = n_dim + 1
    sst_field_meta%diminfo => dim_era_tg
    sst_field_meta%vartype = vartype_real !REAL variable
    sst_field_meta%standard_name = 'T_SEA'
    IF (iera_type == 1) THEN
      sst_field_meta%long_name = 'monthly mean SST climatology 1990-2019'
    ELSEIF (iera_type == 2) THEN
      sst_field_meta%long_name = 'monthly mean SST climatology 1986-2015'
    ENDIF
    sst_field_meta%shortName = 'T_SEA'
    sst_field_meta%units = c_undef
    sst_field_meta%grid_mapping = gridmp
    sst_field_meta%coordinates = coord

    wsnow_field_meta%varname = 'W_SNOW'
    wsnow_field_meta%n_dim = n_dim + 1
    wsnow_field_meta%diminfo => dim_era_tg
    wsnow_field_meta%vartype = vartype_real !REAL variable
    wsnow_field_meta%standard_name = 'W_SNOW'
    IF (iera_type == 1) THEN
      wsnow_field_meta%long_name = 'monthly mean WSNOW climatology 1990-2019'
    ELSEIF (iera_type == 2) THEN
      wsnow_field_meta%long_name = 'monthly mean WSNOW climatology 1986-2015'
    ENDIF
    wsnow_field_meta%shortName = 'W_SNOW'
    wsnow_field_meta%units = c_undef
    wsnow_field_meta%grid_mapping = gridmp
    wsnow_field_meta%coordinates = coord

    t2m_field_meta%varname = 'T_2M_CLIM'
    t2m_field_meta%n_dim = n_dim + 1
    t2m_field_meta%diminfo => dim_era_tg
    t2m_field_meta%vartype = vartype_real !REAL variable
    t2m_field_meta%standard_name = 'T_2M_CLIM'
    IF (iera_type == 1) THEN
      t2m_field_meta%long_name = 'monthly mean T2M climatology 1990-2019'
    ELSEIF (iera_type == 2) THEN
      t2m_field_meta%long_name = 'monthly mean T2M climatology 1986-2015'
    ENDIF
    t2m_field_meta%shortName = 'T_2M_S'
    t2m_field_meta%units = c_undef
    t2m_field_meta%grid_mapping = gridmp
    t2m_field_meta%coordinates = coord

    hsurf_field_meta%varname = 'TOPO_CLIM'
    hsurf_field_meta%n_dim = n_dim
    hsurf_field_meta%diminfo => diminfo
    hsurf_field_meta%vartype = vartype_real !REAL variable
    hsurf_field_meta%standard_name = 'TOPO_CLIM'
    IF (iera_type == 1) THEN
      hsurf_field_meta%long_name = 'TOPO_CLIM for climatology 1990-2019'
    ELSEIF (iera_type == 2) THEN
      hsurf_field_meta%long_name = 'TOPO_CLIM for climatology 1986-2015'
    ENDIF
    hsurf_field_meta%shortName = 'FIS'
    hsurf_field_meta%units = c_undef
    hsurf_field_meta%grid_mapping = gridmp
    hsurf_field_meta%coordinates = coord

  END SUBROUTINE def_era_meta



  !> define dimensions and meta information for variable aot_tg for netcdf output
  SUBROUTINE def_aot_tg_meta(ntime,ntype,diminfo,coordinates,grid_mapping,n_spectr)

    USE mo_aot_data, ONLY : iaot_type, nlevel_cams

    INTEGER (KIND=i4), INTENT(IN) :: ntime !< number of times
    INTEGER (KIND=i4), INTENT(IN) :: ntype !< number of types of aerosols
    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid map
    INTEGER (KIND=i4),  OPTIONAL :: n_spectr !< number of spectral new

    ! local variables
    INTEGER  :: n_dim, &      !< number of dimensions
                nspb          !< number of spectral bands

    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord
    CHARACTER (len=80) :: dataset     ! info dataset

    gridmp = c_undef
    coord = c_undef

    IF (iaot_type == 1 ) THEN
      dataset = 'Tegen JGR 1997 (NASA/GISS)'
    ELSEIF(iaot_type == 2 ) THEN
      dataset = 'AeroCom1 (MPI_MET)'
    ELSEIF(iaot_type == 3 ) THEN
      dataset = 'MACC (ECMWF/KIT)'
    ELSEIF (iaot_type == 4 ) THEN
      dataset = 'MACv2'
      IF (PRESENT(n_spectr)) THEN
        nspb = n_spectr
      ELSE
        nspb = 9
      ENDIF
    ELSEIF (iaot_type == 5) THEN
      dataset = 'CAMS'
    ELSE
      CALL logging%error('Unknown AOT data option', __FILE__, __LINE__)
    ENDIF

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)

    ! set meta information for strucutre dim_aot_tg
    IF (ALLOCATED(dim_aot_tg)) DEALLOCATE(dim_aot_tg)
    ALLOCATE(dim_aot_tg(1:n_dim+2))

    ! set meta information for strucutre dim_aot_ty
    IF (ALLOCATED(dim_aot_ty)) DEALLOCATE(dim_aot_ty)
    ALLOCATE(dim_aot_ty(1:n_dim+1))

    IF (iaot_type == 4) THEN
      dim_aot_tg(1)%dimname = diminfo(1)%dimname
      dim_aot_tg(1)%dimsize = diminfo(1)%dimsize
      dim_aot_tg(2)%dimname = diminfo(2)%dimname
      dim_aot_tg(2)%dimsize = diminfo(2)%dimsize
      dim_aot_tg(3)%dimname = 'spectr'
      dim_aot_tg(3)%dimsize = nspb
      dim_aot_tg(4)%dimname = 'time'
      dim_aot_tg(4)%dimsize = ntime

    ELSE

      SELECT CASE(n_dim)
      CASE (1)
        dim_aot_tg(1)%dimname = diminfo(1)%dimname
        dim_aot_tg(1)%dimsize = diminfo(1)%dimsize

        IF(iaot_type == 5) THEN
          dim_aot_tg(2)%dimname = 'level'
          dim_aot_tg(2)%dimsize = nlevel_cams
        ELSE
          dim_aot_tg(2)%dimname = 'ntype'
          dim_aot_tg(2)%dimsize = ntype
        ENDIF

        dim_aot_tg(3)%dimname = 'time'
        dim_aot_tg(3)%dimsize = ntime

        dim_aot_ty(1) = dim_aot_tg(1)
        dim_aot_ty(2) = dim_aot_tg(3)

      CASE (2)
        dim_aot_tg(1)%dimname = diminfo(1)%dimname
        dim_aot_tg(1)%dimsize = diminfo(1)%dimsize
        dim_aot_tg(2)%dimname = diminfo(2)%dimname
        dim_aot_tg(2)%dimsize = diminfo(2)%dimsize

        IF(iaot_type == 5) THEN
          dim_aot_tg(3)%dimname = 'level'
          dim_aot_tg(3)%dimsize = nlevel_cams
        ELSE
          dim_aot_tg(3)%dimname = 'ntype'
          dim_aot_tg(3)%dimsize = ntype
        ENDIF

        dim_aot_tg(4)%dimname = 'time'
        dim_aot_tg(4)%dimsize = ntime

        dim_aot_ty(1) = dim_aot_tg(1)
        dim_aot_ty(2) = dim_aot_tg(2)
        dim_aot_ty(3) = dim_aot_tg(4)

      CASE (3)
        dim_aot_tg(1)%dimname = diminfo(1)%dimname
        dim_aot_tg(1)%dimsize = diminfo(1)%dimsize
        dim_aot_tg(2)%dimname = diminfo(2)%dimname
        dim_aot_tg(2)%dimsize = diminfo(2)%dimsize

        IF(iaot_type == 5) THEN
          dim_aot_tg(3)%dimname = 'level'
          dim_aot_tg(3)%dimsize = nlevel_cams
          dim_aot_tg(4)%dimname = 'time'
          dim_aot_tg(4)%dimsize = ntime
        ELSE
          dim_aot_tg(3)%dimname = diminfo(3)%dimname
          dim_aot_tg(3)%dimsize = diminfo(3)%dimsize
          dim_aot_tg(4)%dimname = 'ntype'
          dim_aot_tg(4)%dimsize = ntype
          dim_aot_tg(5)%dimname = 'time'
          dim_aot_tg(5)%dimsize = ntime

          dim_aot_ty(1) = dim_aot_tg(1)
          dim_aot_ty(2) = dim_aot_tg(2)
          dim_aot_ty(3) = dim_aot_tg(3)
          dim_aot_ty(4) = dim_aot_tg(5)
        ENDIF


      END SELECT
        ! set meta information for strucutre dim_aot_tg
    ENDIF

    IF (iaot_type == 4) THEN
    !-------------MACv2---------------
      aot_tg_MAC_meta%varname = 'AOT12'

      aot_tg_MAC_meta%n_dim = n_dim + 2
      aot_tg_MAC_meta%diminfo => dim_aot_tg
      aot_tg_MAC_meta%vartype = vartype_real !REAL variable
      aot_tg_MAC_meta%standard_name = c_undef
      aot_tg_MAC_meta%long_name = 'Aerosol optical thickness from Climatology MACv2'
      aot_tg_MAC_meta%shortName = 'AOT'
      aot_tg_MAC_meta%stepType = 'avg'
      aot_tg_MAC_meta%units = c_undef
      aot_tg_MAC_meta%grid_mapping = gridmp
      aot_tg_MAC_meta%coordinates = coord
      aot_tg_MAC_meta%data_set = dataset

      ssa_tg_MAC_meta%varname = 'SSA12'

      ssa_tg_MAC_meta%n_dim = n_dim + 2
      ssa_tg_MAC_meta%diminfo => dim_aot_tg
      ssa_tg_MAC_meta%vartype = vartype_real !REAL variable
      ssa_tg_MAC_meta%standard_name = c_undef
      ssa_tg_MAC_meta%long_name = 'Single scattering albedo from Climatology MACv2'
      ssa_tg_MAC_meta%shortName = 'SSA'
      ssa_tg_MAC_meta%stepType = 'avg'
      ssa_tg_MAC_meta%units = c_undef
      ssa_tg_MAC_meta%grid_mapping = gridmp
      ssa_tg_MAC_meta%coordinates = coord
      ssa_tg_MAC_meta%data_set = dataset

      asy_tg_MAC_meta%varname = 'ASY12'
      asy_tg_MAC_meta%n_dim = n_dim + 2
      asy_tg_MAC_meta%diminfo => dim_aot_tg
      asy_tg_MAC_meta%vartype = vartype_real !REAL variable
      asy_tg_MAC_meta%standard_name = c_undef
      asy_tg_MAC_meta%long_name = 'Factor asymmetry from Climatology MACv2'
      asy_tg_MAC_meta%shortName = 'ASY'
      asy_tg_MAC_meta%stepType = 'avg'
      asy_tg_MAC_meta%units = c_undef
      asy_tg_MAC_meta%grid_mapping = gridmp
      asy_tg_MAC_meta%coordinates = coord
      asy_tg_MAC_meta%data_set = dataset
    ELSEIF (iaot_type == 5) THEN

      CAMS_SS1_tg_meta%varname = 'Sea_Salt_bin1'
      CAMS_SS1_tg_meta%n_dim = n_dim + 2
      CAMS_SS1_tg_meta%diminfo => dim_aot_tg
      CAMS_SS1_tg_meta%vartype = vartype_real !REAL variable
      CAMS_SS1_tg_meta%standard_name = c_undef
      CAMS_SS1_tg_meta%long_name = 'aerosol layer-integrated mass of Sea Salt bin 1 (kg/m**2)'
      CAMS_SS1_tg_meta%shortName = 'AOT_SS1'
      CAMS_SS1_tg_meta%stepType = 'avg'
      CAMS_SS1_tg_meta%units = c_undef
      CAMS_SS1_tg_meta%grid_mapping = gridmp
      CAMS_SS1_tg_meta%coordinates = coord
      CAMS_SS1_tg_meta%data_set = dataset

      CAMS_SS2_tg_meta%varname = 'Sea_Salt_bin2'
      CAMS_SS2_tg_meta%n_dim = n_dim + 2
      CAMS_SS2_tg_meta%diminfo => dim_aot_tg
      CAMS_SS2_tg_meta%vartype = vartype_real !REAL variable
      CAMS_SS2_tg_meta%standard_name = c_undef
      CAMS_SS2_tg_meta%long_name = 'aerosol layer-integrated mass of Sea Salt bin 2 (kg/m**2)'
      CAMS_SS2_tg_meta%shortName = 'AOT_SS2'
      CAMS_SS2_tg_meta%stepType = 'avg'
      CAMS_SS2_tg_meta%units = c_undef
      CAMS_SS2_tg_meta%grid_mapping = gridmp
      CAMS_SS2_tg_meta%coordinates = coord
      CAMS_SS2_tg_meta%data_set = dataset

      CAMS_SS3_tg_meta%varname = 'Sea_Salt_bin3'
      CAMS_SS3_tg_meta%n_dim = n_dim + 2
      CAMS_SS3_tg_meta%diminfo => dim_aot_tg
      CAMS_SS3_tg_meta%vartype = vartype_real !REAL variable
      CAMS_SS3_tg_meta%standard_name = c_undef
      CAMS_SS3_tg_meta%long_name = 'aerosol layer-integrated mass of Sea Salt bin 3 (kg/m**2)'
      CAMS_SS3_tg_meta%shortName = 'AOT_SS3'
      CAMS_SS3_tg_meta%stepType = 'avg'
      CAMS_SS3_tg_meta%units = c_undef
      CAMS_SS3_tg_meta%grid_mapping = gridmp
      CAMS_SS3_tg_meta%coordinates = coord
      CAMS_SS3_tg_meta%data_set = dataset

      CAMS_DUST1_tg_meta%varname = 'Mineral_Dust_bin1'
      CAMS_DUST1_tg_meta%n_dim = n_dim + 2
      CAMS_DUST1_tg_meta%diminfo => dim_aot_tg
      CAMS_DUST1_tg_meta%vartype = vartype_real !REAL variable
      CAMS_DUST1_tg_meta%standard_name = c_undef
      CAMS_DUST1_tg_meta%long_name = 'aerosol layer-integrated mass of Mineral Dust bin 1 (kg/m**2)'
      CAMS_DUST1_tg_meta%shortName = 'AOT_DUST1'
      CAMS_DUST1_tg_meta%stepType = 'avg'
      CAMS_DUST1_tg_meta%units = c_undef
      CAMS_DUST1_tg_meta%grid_mapping = gridmp
      CAMS_DUST1_tg_meta%coordinates = coord
      CAMS_DUST1_tg_meta%data_set = dataset

      CAMS_DUST2_tg_meta%varname = 'Mineral_Dust_bin2'
      CAMS_DUST2_tg_meta%n_dim = n_dim + 2
      CAMS_DUST2_tg_meta%diminfo => dim_aot_tg
      CAMS_DUST2_tg_meta%vartype = vartype_real !REAL variable
      CAMS_DUST2_tg_meta%standard_name = c_undef
      CAMS_DUST2_tg_meta%long_name = 'aerosol layer-integrated mass of Mineral Dust bin 2 (kg/m**2)'
      CAMS_DUST2_tg_meta%shortName = 'AOT_DUST2'
      CAMS_DUST2_tg_meta%stepType = 'avg'
      CAMS_DUST2_tg_meta%units = c_undef
      CAMS_DUST2_tg_meta%grid_mapping = gridmp
      CAMS_DUST2_tg_meta%coordinates = coord
      CAMS_DUST2_tg_meta%data_set = dataset

      CAMS_DUST3_tg_meta%varname = 'Mineral_Dust_bin3'
      CAMS_DUST3_tg_meta%n_dim = n_dim + 2
      CAMS_DUST3_tg_meta%diminfo => dim_aot_tg
      CAMS_DUST3_tg_meta%vartype = vartype_real !REAL variable
      CAMS_DUST3_tg_meta%standard_name = c_undef
      CAMS_DUST3_tg_meta%long_name = 'aerosol layer-integrated mass of Mineral Dust bin 3 (kg/m**2)'
      CAMS_DUST3_tg_meta%shortName = 'AOT_DUST3'
      CAMS_DUST3_tg_meta%stepType = 'avg'
      CAMS_DUST3_tg_meta%units = c_undef
      CAMS_DUST3_tg_meta%grid_mapping = gridmp
      CAMS_DUST3_tg_meta%coordinates = coord
      CAMS_DUST3_tg_meta%data_set = dataset

      CAMS_OCphilic_tg_meta%varname = 'OC_hydrophilic'
      CAMS_OCphilic_tg_meta%n_dim = n_dim + 2
      CAMS_OCphilic_tg_meta%diminfo => dim_aot_tg
      CAMS_OCphilic_tg_meta%vartype = vartype_real !REAL variable
      CAMS_OCphilic_tg_meta%standard_name = c_undef
      CAMS_OCphilic_tg_meta%long_name = 'aerosol layer-integrated mass of Organic Matter hydrophilic(kg/m**2)'
      CAMS_OCphilic_tg_meta%shortName = 'AOT_OCphilic'
      CAMS_OCphilic_tg_meta%stepType = 'avg'
      CAMS_OCphilic_tg_meta%units = c_undef
      CAMS_OCphilic_tg_meta%grid_mapping = gridmp
      CAMS_OCphilic_tg_meta%coordinates = coord
      CAMS_OCphilic_tg_meta%data_set = dataset

      CAMS_OCphobic_tg_meta%varname = 'OC_hydrophobic'
      CAMS_OCphobic_tg_meta%n_dim = n_dim + 2
      CAMS_OCphobic_tg_meta%diminfo => dim_aot_tg
      CAMS_OCphobic_tg_meta%vartype = vartype_real !REAL variable
      CAMS_OCphobic_tg_meta%standard_name = c_undef
      CAMS_OCphobic_tg_meta%long_name = 'aerosol layer-integrated mass of Organic Matter hydrophobic (kg/m**2)'
      CAMS_OCphobic_tg_meta%shortName = 'AOT_OCphobic'
      CAMS_OCphobic_tg_meta%stepType = 'avg'
      CAMS_OCphobic_tg_meta%units = c_undef
      CAMS_OCphobic_tg_meta%grid_mapping = gridmp
      CAMS_OCphobic_tg_meta%coordinates = coord
      CAMS_OCphobic_tg_meta%data_set = dataset

      CAMS_BCphilic_tg_meta%varname = 'BC_hydrophilic'
      CAMS_BCphilic_tg_meta%n_dim = n_dim + 2
      CAMS_BCphilic_tg_meta%diminfo => dim_aot_tg
      CAMS_BCphilic_tg_meta%vartype = vartype_real !REAL variable
      CAMS_BCphilic_tg_meta%standard_name = c_undef
      CAMS_BCphilic_tg_meta%long_name = 'aerosol layer-integrated mass of Black Carbon hydrophilic (kg/m**2)'
      CAMS_BCphilic_tg_meta%shortName = 'AOT_BCphilic'
      CAMS_BCphilic_tg_meta%stepType = 'avg'
      CAMS_BCphilic_tg_meta%units = c_undef
      CAMS_BCphilic_tg_meta%grid_mapping = gridmp
      CAMS_BCphilic_tg_meta%coordinates = coord
      CAMS_BCphilic_tg_meta%data_set = dataset

      CAMS_BCphobic_tg_meta%varname = 'BC_hydrophobic'
      CAMS_BCphobic_tg_meta%n_dim = n_dim + 2
      CAMS_BCphobic_tg_meta%diminfo => dim_aot_tg
      CAMS_BCphobic_tg_meta%vartype = vartype_real !REAL variable
      CAMS_BCphobic_tg_meta%standard_name = c_undef
      CAMS_BCphobic_tg_meta%long_name = 'aerosol layer-integrated mass of Black Carbon hydrophobic (kg/m**2)'
      CAMS_BCphobic_tg_meta%shortName = 'AOT_BCphobic'
      CAMS_BCphobic_tg_meta%stepType = 'avg'
      CAMS_BCphobic_tg_meta%units = c_undef
      CAMS_BCphobic_tg_meta%grid_mapping = gridmp
      CAMS_BCphobic_tg_meta%coordinates = coord
      CAMS_BCphobic_tg_meta%data_set = dataset

      CAMS_SU_tg_meta%varname = 'Sulfates'
      CAMS_SU_tg_meta%n_dim = n_dim + 2
      CAMS_SU_tg_meta%diminfo => dim_aot_tg
      CAMS_SU_tg_meta%vartype = vartype_real !REAL variable
      CAMS_SU_tg_meta%standard_name = c_undef
      CAMS_SU_tg_meta%long_name = 'aerosol layer-integrated mass of Sulfates (kg/m**2)'
      CAMS_SU_tg_meta%shortName = 'AOT_SU'
      CAMS_SU_tg_meta%stepType = 'avg'
      CAMS_SU_tg_meta%units = c_undef
      CAMS_SU_tg_meta%grid_mapping = gridmp
      CAMS_SU_tg_meta%coordinates = coord
      CAMS_SU_tg_meta%data_set = dataset

      CAMS_plev_tg_meta%varname = 'half_level_pressure'
      CAMS_plev_tg_meta%n_dim = n_dim + 2
      CAMS_plev_tg_meta%diminfo => dim_aot_tg
      CAMS_plev_tg_meta%vartype = vartype_real !REAL variable
      CAMS_plev_tg_meta%standard_name = c_undef
      CAMS_plev_tg_meta%long_name = 'half level pressure (Pa)'
      CAMS_plev_tg_meta%shortName = 'p_lev_CAMS'
      CAMS_plev_tg_meta%stepType = 'avg'
      CAMS_plev_tg_meta%units = c_undef
      CAMS_plev_tg_meta%grid_mapping = gridmp
      CAMS_plev_tg_meta%coordinates = coord
      CAMS_plev_tg_meta%data_set = dataset
    ELSE
      ! set meta information for variable aot_tg
      aot_tg_meta%varname = 'AOT_TG'
      aot_tg_meta%n_dim = n_dim + 2
      aot_tg_meta%diminfo => dim_aot_tg
      aot_tg_meta%vartype = vartype_real !REAL variable
      aot_tg_meta%standard_name = c_undef !_br 08.04.14
      aot_tg_meta%long_name = 'aerosol optical thickness'
      aot_tg_meta%shortName = 'AOT'
      aot_tg_meta%stepType = 'avg'
      aot_tg_meta%units = c_undef
      aot_tg_meta%grid_mapping = gridmp
      aot_tg_meta%coordinates = coord
      aot_tg_meta%data_set = dataset

      aot_type_shortname(1) = 'AER_BC12'
      aot_type_shortname(2) = 'AER_DUST12'
      aot_type_shortname(3) = 'AER_ORG12'
      aot_type_shortname(4) = 'AER_SO412'
      aot_type_shortname(5) = 'AER_SS12'

      aer_bc_meta%varname = 'AER_BC12'
      aer_bc_meta%n_dim = n_dim + 1
      aer_bc_meta%diminfo => dim_aot_ty
      aer_bc_meta%vartype = vartype_real !REAL variable
      aer_bc_meta%standard_name = 'atmosphere_absorption_optical_thickness_due_to_black_carbon_ambient_aerosol'
      aer_bc_meta%long_name = 'aerosol optical thickness of black carbon'
      aer_bc_meta%shortName = 'AER_BC12'
      aer_bc_meta%stepType = 'avg'
      aer_bc_meta%units = c_undef
      aer_bc_meta%grid_mapping = gridmp
      aer_bc_meta%coordinates = coord
      aer_bc_meta%data_set = dataset

      aer_dust_meta%varname = 'AER_DUST12'
      aer_dust_meta%n_dim = n_dim + 1
      aer_dust_meta%diminfo => dim_aot_ty
      aer_dust_meta%vartype = vartype_real !REAL variable
      aer_dust_meta%standard_name = c_undef !_br 08.04.14
      aer_dust_meta%long_name = 'atmosphere_absorption_optical_thickness_due_to_dust_ambient_aerosol'
      aer_dust_meta%shortName = 'AER_DUST12'
      aer_dust_meta%stepType = 'avg'
      aer_dust_meta%units = c_undef
      aer_dust_meta%grid_mapping = gridmp
      aer_dust_meta%coordinates = coord
      aer_dust_meta%data_set = dataset

      aer_org_meta%varname = 'AER_ORG12'
      aer_org_meta%n_dim = n_dim + 1
      aer_org_meta%diminfo => dim_aot_ty
      aer_org_meta%vartype = vartype_real !REAL variable
      aer_org_meta%standard_name = c_undef !_br 08.04.14
      aer_org_meta%long_name = 'atmosphere_absorption_optical_thickness_due_to_particulate_organic_matter_ambient_aerosol'
      aer_org_meta%shortName = 'AER_ORG12'
      aer_org_meta%stepType = 'avg'
      aer_org_meta%units = c_undef
      aer_org_meta%grid_mapping = gridmp
      aer_org_meta%coordinates = coord
      aer_org_meta%data_set = dataset

      aer_so4_meta%varname = 'AER_SO412'
      aer_so4_meta%n_dim = n_dim + 1
      aer_so4_meta%diminfo => dim_aot_ty
      aer_so4_meta%vartype = vartype_real !REAL variable
      aer_so4_meta%standard_name = c_undef !_br 08.04.14
      aer_so4_meta%long_name = 'atmosphere_absorption_optical_thickness_due_to_sulfate_ambient_aerosol'
      aer_so4_meta%shortName = 'AER_SO412'
      aer_so4_meta%stepType = 'avg'
      aer_so4_meta%units = c_undef
      aer_so4_meta%grid_mapping = gridmp
      aer_so4_meta%coordinates = coord
      aer_so4_meta%data_set = dataset

      aer_ss_meta%varname = 'AER_SS12'
      aer_ss_meta%n_dim = n_dim + 1
      aer_ss_meta%diminfo => dim_aot_ty
      aer_ss_meta%vartype = vartype_real !REAL variable
      aer_ss_meta%standard_name = c_undef !_br 08.04.14
      aer_ss_meta%long_name = 'atmosphere_absorption_optical_thickness_due_to_seasalt_ambient_aerosol'
      aer_ss_meta%shortName = 'AER_SS12'
      aer_ss_meta%stepType = 'avg'
      aer_ss_meta%units = c_undef
      aer_ss_meta%grid_mapping = gridmp
      aer_ss_meta%coordinates = coord
      aer_ss_meta%data_set = dataset

    ENDIF

  END SUBROUTINE def_aot_tg_meta


  ! define meta information for target field variables lon_geo, lat_geo and no_raw_data_pixel
  SUBROUTINE def_com_target_fields_meta(diminfo,coordinates,grid_mapping)

    TYPE(dim_meta_info),TARGET   :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER                      :: n_dim      !< number of dimensions
    CHARACTER (len=80)           :: gridmp
    CHARACTER (len=80)           :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)

    lon_geo_meta%varname = 'lon'
    lon_geo_meta%n_dim = n_dim
    lon_geo_meta%diminfo => diminfo
    lon_geo_meta%vartype = vartype_real !REAL variable
    lon_geo_meta%standard_name = 'longitude'
    lon_geo_meta%long_name = 'geographical longitude'
    lon_geo_meta%shortName = 'rlon'
    lon_geo_meta%stepType = 'instant'
    lon_geo_meta%units =  'degrees_east'
    lon_geo_meta%grid_mapping = c_undef
    lon_geo_meta%coordinates = c_undef
    lon_geo_meta%data_set = c_undef

    lat_geo_meta%varname = 'lat'
    lat_geo_meta%n_dim = n_dim
    lat_geo_meta%diminfo => diminfo
    lat_geo_meta%vartype = vartype_real !REAL variable
    lat_geo_meta%standard_name = 'latitude'
    lat_geo_meta%long_name = 'geographical latitude'
    lat_geo_meta%shortName = 'rlat'
    lat_geo_meta%stepType = 'instant'
    lat_geo_meta%units =  'degrees_north'
    lat_geo_meta%grid_mapping = c_undef
    lat_geo_meta%coordinates = c_undef
    lat_geo_meta%data_set = c_undef

    no_raw_data_pixel_meta%varname = 'NO_RAW_DATA_PIXEL'
    no_raw_data_pixel_meta%n_dim = n_dim
    no_raw_data_pixel_meta%diminfo => diminfo
    no_raw_data_pixel_meta%vartype = vartype_int !REAL variable
    no_raw_data_pixel_meta%standard_name = c_undef !_br 08.04.14
    no_raw_data_pixel_meta%long_name = 'number of raw data pixel in target grid element'
    no_raw_data_pixel_meta%shortName = c_undef
    no_raw_data_pixel_meta%stepType = 'instant'
    no_raw_data_pixel_meta%units = c_undef
    no_raw_data_pixel_meta%grid_mapping = gridmp
    no_raw_data_pixel_meta%coordinates = coord
    no_raw_data_pixel_meta%data_set = c_undef


  END SUBROUTINE def_com_target_fields_meta

  !> define meta information for GLC2000 target fields
  SUBROUTINE def_glc2000_fields_meta(nclass_glc2000,diminfo,coordinates,grid_mapping)
    INTEGER (KIND=i4) :: nclass_glc2000 !< GLC2000 has 23 classes for the land use description
    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)

    n_dim = SIZE(diminfo)

    ! set meta information for strucutre dim_glc2000_tg
    IF (ALLOCATED(dim_glc2000_tg)) DEALLOCATE(dim_glc2000_tg)
    ALLOCATE(dim_glc2000_tg(1:n_dim+1))
    SELECT CASE(n_dim)
    CASE (1)
      dim_glc2000_tg(1)%dimname = diminfo(1)%dimname
      dim_glc2000_tg(1)%dimsize = diminfo(1)%dimsize
      dim_glc2000_tg(2)%dimname = 'nclass'
      dim_glc2000_tg(2)%dimsize = nclass_glc2000
    CASE (2)
      dim_glc2000_tg(1)%dimname = diminfo(1)%dimname
      dim_glc2000_tg(1)%dimsize = diminfo(1)%dimsize
      dim_glc2000_tg(2)%dimname = diminfo(2)%dimname
      dim_glc2000_tg(2)%dimsize = diminfo(2)%dimsize
      dim_glc2000_tg(3)%dimname = 'nclass'
      dim_glc2000_tg(3)%dimsize = nclass_glc2000
    CASE (3)
      dim_glc2000_tg(1)%dimname = diminfo(1)%dimname
      dim_glc2000_tg(1)%dimsize = diminfo(1)%dimsize
      dim_glc2000_tg(2)%dimname = diminfo(2)%dimname
      dim_glc2000_tg(2)%dimsize = diminfo(2)%dimsize
      dim_glc2000_tg(3)%dimname = diminfo(3)%dimname
      dim_glc2000_tg(3)%dimsize = diminfo(3)%dimsize
      dim_glc2000_tg(4)%dimname = 'nclass'
      dim_glc2000_tg(4)%dimsize = nclass_glc2000
    END SELECT




    ! fr_land_glc2000_meta
    fr_land_glc2000_meta%varname = 'FR_LAND_GLC2000'
    fr_land_glc2000_meta%n_dim = n_dim
    fr_land_glc2000_meta%diminfo => diminfo
    fr_land_glc2000_meta%vartype = vartype_real !REAL variable
    fr_land_glc2000_meta%standard_name = 'land_area_fraction' !_br 08.04.14
    fr_land_glc2000_meta%long_name = 'Fraction land due to GLC2000 Data'
    fr_land_glc2000_meta%shortName = 'FR_LAND'
    fr_land_glc2000_meta%stepType = 'instant'
    fr_land_glc2000_meta%units =  c_undef
    fr_land_glc2000_meta%grid_mapping = gridmp
    fr_land_glc2000_meta%coordinates = coord
    fr_land_glc2000_meta%data_set = 'GLC2000'


   ! glc2000_tot_npixel_meta
    glc2000_tot_npixel_meta%varname = 'GLC2000_TOT_NPIXEL'
    glc2000_tot_npixel_meta%n_dim = n_dim
    glc2000_tot_npixel_meta%diminfo => diminfo
    glc2000_tot_npixel_meta%vartype = vartype_int !INTEGER variable
    glc2000_tot_npixel_meta%standard_name = c_undef !_br 08.04.14
    glc2000_tot_npixel_meta%long_name = 'number of raw data pixel in target grid element'
    glc2000_tot_npixel_meta%shortName = c_undef
    glc2000_tot_npixel_meta%stepType = 'instant'
    glc2000_tot_npixel_meta%units = c_undef
    glc2000_tot_npixel_meta%grid_mapping = gridmp
    glc2000_tot_npixel_meta%coordinates = coord
    glc2000_tot_npixel_meta%data_set = 'GLC2000'


    ! glc2000_class_fraction_meta
    glc2000_class_fraction_meta%varname = 'GLC2000_CLASS_FRACTION'
    glc2000_class_fraction_meta%n_dim = n_dim + 1
    glc2000_class_fraction_meta%diminfo => dim_glc2000_tg
    glc2000_class_fraction_meta%vartype = vartype_real !REAL variable
    glc2000_class_fraction_meta%standard_name = c_undef !_br 08.04.14
    glc2000_class_fraction_meta%long_name = 'Fraction of GLC2000 land use classes in target grid element'
    glc2000_class_fraction_meta%shortName = c_undef
    glc2000_class_fraction_meta%stepType = 'instant'
    glc2000_class_fraction_meta%units =  c_undef
    glc2000_class_fraction_meta%grid_mapping = gridmp
    glc2000_class_fraction_meta%coordinates = coord
    glc2000_class_fraction_meta%data_set = 'GLC2000'

    ! glc2000_class_npixel_meta
    glc2000_class_npixel_meta%varname = 'GLC2000_CLASS_NPIXEL'
    glc2000_class_npixel_meta%n_dim = n_dim + 1
    glc2000_class_npixel_meta%diminfo => dim_glc2000_tg
    glc2000_class_npixel_meta%vartype = vartype_int !INTEGER variable
    glc2000_class_npixel_meta%standard_name = c_undef !_br 08.04.14
    glc2000_class_npixel_meta%long_name = 'number of pixels of GLC2000 land use classes in target grid element'
    glc2000_class_npixel_meta%shortName = c_undef
    glc2000_class_npixel_meta%stepType = 'instant'
    glc2000_class_npixel_meta%units = c_undef
    glc2000_class_npixel_meta%grid_mapping = gridmp
    glc2000_class_npixel_meta%coordinates = coord
    glc2000_class_npixel_meta%data_set = 'GLC2000'

    ! ice_glc2000_meta
    ice_glc2000_meta%varname = 'ICE_GLC2000'
    ice_glc2000_meta%n_dim = n_dim
    ice_glc2000_meta%diminfo => diminfo
    ice_glc2000_meta%vartype = vartype_real !REAL variable
    ice_glc2000_meta%standard_name = c_undef !_br 08.04.14
    ice_glc2000_meta%long_name = 'Ice fraction due to GLC2000 Data'
    ice_glc2000_meta%shortName = c_undef
    ice_glc2000_meta%stepType = 'instant'
    ice_glc2000_meta%units =  c_undef
    ice_glc2000_meta%grid_mapping = gridmp
    ice_glc2000_meta%coordinates = coord
    ice_glc2000_meta%data_set = 'GLC2000'

    ! z0_glc2000_meta
    z0_glc2000_meta%varname = 'Z0_GLC2000'
    z0_glc2000_meta%n_dim = n_dim
    z0_glc2000_meta%diminfo => diminfo
    z0_glc2000_meta%vartype = vartype_real !REAL variable
    z0_glc2000_meta%standard_name = 'surface_roughness_length' !_br 08.04.14
    z0_glc2000_meta%long_name =  'Roughness length z0 due to GLC2000 land use data'
    z0_glc2000_meta%shortName = 'Z0_LU'
    z0_glc2000_meta%stepType = 'instant'
    z0_glc2000_meta%units = 'm'
    z0_glc2000_meta%grid_mapping = gridmp
    z0_glc2000_meta%coordinates = coord
    z0_glc2000_meta%data_set = 'GLC2000'

    ! root_glc2000_meta
    root_glc2000_meta%varname = 'ROOT_GLC2000'
    root_glc2000_meta%n_dim = n_dim
    root_glc2000_meta%diminfo => diminfo
    root_glc2000_meta%vartype = vartype_real !REAL variable
    root_glc2000_meta%standard_name = 'root_depth' !_br 08.04.14
    root_glc2000_meta%long_name = 'Root depth due to GLC2000 land use data'
    root_glc2000_meta%shortName = 'ROOTDP'
    root_glc2000_meta%stepType = 'instant'
    root_glc2000_meta%units =  'm'
    root_glc2000_meta%grid_mapping = gridmp
    root_glc2000_meta%coordinates = coord
    root_glc2000_meta%data_set = 'GLC2000'

    ! plcov_mx_glc2000_meta
    plcov_mx_glc2000_meta%varname = 'PLCOV_MX_GLC2000'
    plcov_mx_glc2000_meta%n_dim = n_dim
    plcov_mx_glc2000_meta%diminfo => diminfo
    plcov_mx_glc2000_meta%vartype = vartype_real !REAL variable
    plcov_mx_glc2000_meta%standard_name = 'vegetation_area_fraction' !_br 08.04.14
    plcov_mx_glc2000_meta%long_name = 'Plant cover maximum due to GLC2000 land use data'
    plcov_mx_glc2000_meta%shortName = 'PLCOV_MX'
    plcov_mx_glc2000_meta%stepType = 'max'
    plcov_mx_glc2000_meta%units =  c_undef
    plcov_mx_glc2000_meta%grid_mapping = gridmp
    plcov_mx_glc2000_meta%coordinates = coord
    plcov_mx_glc2000_meta%data_set = 'GLC2000'


    ! plcov_mn_glc2000_meta
    plcov_mn_glc2000_meta%varname = 'PLCOV_MN_GLC2000'
    plcov_mn_glc2000_meta%n_dim = n_dim
    plcov_mn_glc2000_meta%diminfo => diminfo
    plcov_mn_glc2000_meta%vartype = vartype_real !REAL variable
    plcov_mn_glc2000_meta%standard_name = 'vegetation_area_fraction' !_br 08.04.14
    plcov_mn_glc2000_meta%long_name = 'Plant cover minimum due to GLC2000 land use data'
    plcov_mn_glc2000_meta%shortName = 'PLCOV_MN'
    plcov_mn_glc2000_meta%stepType = 'min'
    plcov_mn_glc2000_meta%units =  c_undef
    plcov_mn_glc2000_meta%grid_mapping = gridmp
    plcov_mn_glc2000_meta%coordinates = coord
    plcov_mn_glc2000_meta%data_set = 'GLC2000'

    ! lai_mx_glc2000_meta
    lai_mx_glc2000_meta%varname = 'LAI_MX_GLC2000'
    lai_mx_glc2000_meta%n_dim = n_dim
    lai_mx_glc2000_meta%diminfo => diminfo
    lai_mx_glc2000_meta%vartype = vartype_real !REAL variable
    lai_mx_glc2000_meta%standard_name = 'leaf_area_index' !_br 08.04.14
    lai_mx_glc2000_meta%long_name = 'Leaf Area Index Maximum'
    lai_mx_glc2000_meta%shortName = 'LAI_MX'
    lai_mx_glc2000_meta%stepType = 'max'
    lai_mx_glc2000_meta%units =  c_undef
    lai_mx_glc2000_meta%grid_mapping = gridmp
    lai_mx_glc2000_meta%coordinates = coord
    lai_mx_glc2000_meta%data_set = 'GLC2000'


    ! lai_mn_glc2000_meta
    lai_mn_glc2000_meta%varname = 'LAI_MN_GLC2000'
    lai_mn_glc2000_meta%n_dim = n_dim
    lai_mn_glc2000_meta%diminfo => diminfo
    lai_mn_glc2000_meta%vartype = vartype_real !REAL variable
    lai_mn_glc2000_meta%standard_name = 'leaf_area_index' !_br 08.04.14
    lai_mn_glc2000_meta%long_name = 'Leaf Area Minimum'
    lai_mn_glc2000_meta%shortName = 'LAI_MN'
    lai_mn_glc2000_meta%stepType = 'min'
    lai_mn_glc2000_meta%units =  c_undef
    lai_mn_glc2000_meta%grid_mapping = gridmp
    lai_mn_glc2000_meta%coordinates = coord
    lai_mn_glc2000_meta%data_set = 'GLC2000'


    ! rs_min_glc2000_meta
    rs_min_glc2000_meta%varname = 'RSMIN_GLC2000'
    rs_min_glc2000_meta%n_dim = n_dim
    rs_min_glc2000_meta%diminfo => diminfo
    rs_min_glc2000_meta%vartype = vartype_real !REAL variable
    rs_min_glc2000_meta%standard_name = c_undef !_br 08.04.14
    rs_min_glc2000_meta%long_name = 'Minimal stomata resistence'
    rs_min_glc2000_meta%shortName = 'RSMIN'
    rs_min_glc2000_meta%stepType = 'instant'
    rs_min_glc2000_meta%units =  's/m'
    rs_min_glc2000_meta%grid_mapping =gridmp
    rs_min_glc2000_meta%coordinates = coord
    rs_min_glc2000_meta%data_set = 'GLC2000'

    ! urban_glc2000_meta
    urban_glc2000_meta%varname = 'URBAN_GLC2000'
    urban_glc2000_meta%n_dim = n_dim
    urban_glc2000_meta%diminfo => diminfo
    urban_glc2000_meta%vartype = vartype_real !REAL variable
    urban_glc2000_meta%standard_name = c_undef !_br 08.04.14
    urban_glc2000_meta%long_name = 'Urban land use fraction'
    urban_glc2000_meta%shortName = 'URBAN'
    urban_glc2000_meta%stepType = 'instant'
    urban_glc2000_meta%units =  c_undef
    urban_glc2000_meta%grid_mapping = gridmp
    urban_glc2000_meta%coordinates = coord
    urban_glc2000_meta%data_set = 'GLC2000'


    ! for_d_glc2000_meta
    for_d_glc2000_meta%varname = 'FOR_D_GLC2000'
    for_d_glc2000_meta%n_dim = n_dim
    for_d_glc2000_meta%diminfo => diminfo
    for_d_glc2000_meta%vartype = vartype_real !REAL variable
    for_d_glc2000_meta%standard_name = c_undef !_br 08.04.14
    for_d_glc2000_meta%long_name = 'Fraction of deciduous forest'
    for_d_glc2000_meta%shortName = 'FOR_D'
    for_d_glc2000_meta%stepType = 'instant'
    for_d_glc2000_meta%units = '1'
    for_d_glc2000_meta%grid_mapping = gridmp
    for_d_glc2000_meta%coordinates = coord
    for_d_glc2000_meta%data_set = 'GLC2000'


    ! for_e_glc2000_meta
    for_e_glc2000_meta%varname = 'FOR_E_GLC2000'
    for_e_glc2000_meta%n_dim = n_dim
    for_e_glc2000_meta%diminfo => diminfo
    for_e_glc2000_meta%vartype = vartype_real !REAL variable
    for_e_glc2000_meta%standard_name = c_undef !_br 08.04.14
    for_e_glc2000_meta%long_name = 'Fraction of evergreen forest'
    for_e_glc2000_meta%shortName =  'FOR_E'
    for_e_glc2000_meta%stepType = 'instant'
    for_e_glc2000_meta%units =  '1'
    for_e_glc2000_meta%grid_mapping = gridmp
    for_e_glc2000_meta%coordinates = coord
    for_e_glc2000_meta%data_set = 'GLC2000'


    ! emissivity_glc2000_meta
    emissivity_glc2000_meta%varname = 'EMISSIVITY_GLC2000'
    emissivity_glc2000_meta%n_dim = n_dim
    emissivity_glc2000_meta%diminfo => diminfo
    emissivity_glc2000_meta%vartype = vartype_real !REAL variable
    emissivity_glc2000_meta%standard_name = c_undef !_br 08.04.14
    emissivity_glc2000_meta%long_name = 'longwave surface emissivity'
    emissivity_glc2000_meta%shortName = 'EMIS_RAD'
    emissivity_glc2000_meta%stepType = 'instant'
    emissivity_glc2000_meta%units =  '1'
    emissivity_glc2000_meta%grid_mapping = gridmp
    emissivity_glc2000_meta%coordinates = coord
    emissivity_glc2000_meta%data_set = 'GLC2000'



  END SUBROUTINE def_glc2000_fields_meta

    !> define meta information for GLCC target fields
  SUBROUTINE def_glcc_fields_meta(nclass_glcc,diminfo,coordinates,grid_mapping)
    INTEGER (KIND=i4) :: nclass_glcc !< GLCC has 23 classes for the land use description
    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)

    n_dim = SIZE(diminfo)


    ! set meta information for strucutre dim_glcc_tg
    IF (ALLOCATED(dim_glcc_tg)) DEALLOCATE(dim_glcc_tg)
    ALLOCATE(dim_glcc_tg(1:n_dim+1))
    SELECT CASE(n_dim)
    CASE (1)
      dim_glcc_tg(1)%dimname = diminfo(1)%dimname
      dim_glcc_tg(1)%dimsize = diminfo(1)%dimsize
      dim_glcc_tg(2)%dimname = 'nclass'
      dim_glcc_tg(2)%dimsize = nclass_glcc
    CASE (2)
      dim_glcc_tg(1)%dimname = diminfo(1)%dimname
      dim_glcc_tg(1)%dimsize = diminfo(1)%dimsize
      dim_glcc_tg(2)%dimname = diminfo(2)%dimname
      dim_glcc_tg(2)%dimsize = diminfo(2)%dimsize
      dim_glcc_tg(3)%dimname = 'nclass'
      dim_glcc_tg(3)%dimsize = nclass_glcc
    CASE (3)
      dim_glcc_tg(1)%dimname = diminfo(1)%dimname
      dim_glcc_tg(1)%dimsize = diminfo(1)%dimsize
      dim_glcc_tg(2)%dimname = diminfo(2)%dimname
      dim_glcc_tg(2)%dimsize = diminfo(2)%dimsize
      dim_glcc_tg(3)%dimname = diminfo(3)%dimname
      dim_glcc_tg(3)%dimsize = diminfo(3)%dimsize
      dim_glcc_tg(4)%dimname = 'nclass'
      dim_glcc_tg(4)%dimsize = nclass_glcc
    END SELECT



    ! fr_land_glcc_meta
    fr_land_glcc_meta%varname = 'FR_LAND_GLCC'
    fr_land_glcc_meta%n_dim = n_dim
    fr_land_glcc_meta%diminfo => diminfo
    fr_land_glcc_meta%vartype = vartype_real !REAL variable
    fr_land_glcc_meta%standard_name = 'land_area_fraction' !_br 08.04.14
    fr_land_glcc_meta%long_name = 'Fraction land due to GLCC Data'
    fr_land_glcc_meta%shortName = 'FR_LAND'
    fr_land_glcc_meta%stepType = 'instant'
    fr_land_glcc_meta%units =  c_undef
    fr_land_glcc_meta%grid_mapping = gridmp
    fr_land_glcc_meta%coordinates = coord
    fr_land_glcc_meta%data_set = 'GLCC'


   ! glcc_tot_npixel_meta
    glcc_tot_npixel_meta%varname = 'GLCC_TOT_NPIXEL'
    glcc_tot_npixel_meta%n_dim = n_dim
    glcc_tot_npixel_meta%diminfo => dim_glcc_tg
    glcc_tot_npixel_meta%vartype = vartype_int !INTEGER variable
    glcc_tot_npixel_meta%standard_name = c_undef !_br 08.04.14
    glcc_tot_npixel_meta%long_name = 'number of raw data pixel in target grid element'
    glcc_tot_npixel_meta%shortName = c_undef
    glcc_tot_npixel_meta%stepType = 'instant'
    glcc_tot_npixel_meta%units = c_undef
    glcc_tot_npixel_meta%grid_mapping = gridmp
    glcc_tot_npixel_meta%coordinates = coord
    glcc_tot_npixel_meta%data_set = 'GLCC'


    ! glcc_class_fraction_meta
    glcc_class_fraction_meta%varname = 'GLCC_CLASS_FRACTION'
    glcc_class_fraction_meta%n_dim = n_dim + 1
    glcc_class_fraction_meta%diminfo => dim_glcc_tg
    glcc_class_fraction_meta%vartype = vartype_real !REAL variable
    glcc_class_fraction_meta%standard_name = c_undef !_br 08.04.14
    glcc_class_fraction_meta%long_name = 'Fraction of GLCC land use classes in target grid element'
    glcc_class_fraction_meta%shortName = c_undef
    glcc_class_fraction_meta%stepType = 'instant'
    glcc_class_fraction_meta%units =  c_undef
    glcc_class_fraction_meta%grid_mapping = gridmp
    glcc_class_fraction_meta%coordinates = coord
    glcc_class_fraction_meta%data_set = 'GLCC'

    ! glcc_class_npixel_meta
    glcc_class_npixel_meta%varname = 'GLCC_CLASS_NPIXEL'
    glcc_class_npixel_meta%n_dim = n_dim + 1
    glcc_class_npixel_meta%diminfo => dim_glcc_tg
    glcc_class_npixel_meta%vartype = vartype_int !INTEGER variable
    glcc_class_npixel_meta%standard_name = c_undef !_br 08.04.14
    glcc_class_npixel_meta%long_name = 'number of pixels of GLCC land use classes in target grid element'
    glcc_class_npixel_meta%shortName = c_undef
    glcc_class_npixel_meta%stepType = 'instant'
    glcc_class_npixel_meta%units = c_undef
    glcc_class_npixel_meta%grid_mapping = gridmp
    glcc_class_npixel_meta%coordinates = coord
    glcc_class_npixel_meta%data_set = 'GLCC'

    ! ice_glcc_meta
    ice_glcc_meta%varname = 'ICE_GLCC'
    ice_glcc_meta%n_dim = n_dim
    ice_glcc_meta%diminfo => diminfo
    ice_glcc_meta%vartype = vartype_real !REAL variable
    ice_glcc_meta%standard_name = c_undef !_br 08.04.14
    ice_glcc_meta%long_name = 'Ice fraction due to GLCC Data'
    ice_glcc_meta%shortName = c_undef
    ice_glcc_meta%stepType = 'instant'
    ice_glcc_meta%units =  '1'
    ice_glcc_meta%grid_mapping = gridmp
    ice_glcc_meta%coordinates = coord
    ice_glcc_meta%data_set = 'GLCC'

    ! z0_glcc_meta
    z0_glcc_meta%varname = 'Z0_GLCC'
    z0_glcc_meta%n_dim = n_dim
    z0_glcc_meta%diminfo => diminfo
    z0_glcc_meta%vartype = vartype_real !REAL variable
    z0_glcc_meta%standard_name = 'surface_roughness_length' !_br 08.04.14
    z0_glcc_meta%long_name =  'Roughness length z0 due to GLCC land use data'
    z0_glcc_meta%shortName = 'Z0_LU'
    z0_glcc_meta%stepType = 'instant'
    z0_glcc_meta%units = 'm'
    z0_glcc_meta%grid_mapping = gridmp
    z0_glcc_meta%coordinates = coord
    z0_glcc_meta%data_set = 'GLCC'

    ! root_glcc_meta
    root_glcc_meta%varname = 'ROOT_GLCC'
    root_glcc_meta%n_dim = n_dim
    root_glcc_meta%diminfo => diminfo
    root_glcc_meta%vartype = vartype_real !REAL variable
    root_glcc_meta%standard_name = c_undef !_br 08.04.14
    root_glcc_meta%long_name = 'Root depth due to GLCC land use data'
    root_glcc_meta%shortName = 'ROOTDP'
    root_glcc_meta%stepType = 'instant'
    root_glcc_meta%units =  'm'
    root_glcc_meta%grid_mapping = gridmp
    root_glcc_meta%coordinates = coord
    root_glcc_meta%data_set = 'GLCC'

    ! plcov_mx_glcc_meta
    plcov_mx_glcc_meta%varname = 'PLCOV_MX_GLCC'
    plcov_mx_glcc_meta%n_dim = n_dim
    plcov_mx_glcc_meta%diminfo => diminfo
    plcov_mx_glcc_meta%vartype = vartype_real !REAL variable
    plcov_mx_glcc_meta%standard_name = 'vegetation_area_fraction' !_br 08.04.14
    plcov_mx_glcc_meta%long_name = 'Plant cover maximum due to GLCC land use data'
    plcov_mx_glcc_meta%shortName = 'PLCOV_MX'
    plcov_mx_glcc_meta%stepType = 'max'
    plcov_mx_glcc_meta%units =  '1'
    plcov_mx_glcc_meta%grid_mapping = gridmp
    plcov_mx_glcc_meta%coordinates = coord
    plcov_mx_glcc_meta%data_set = 'GLCC'


    ! plcov_mn_glcc_meta
    plcov_mn_glcc_meta%varname = 'PLCOV_MN_GLCC'
    plcov_mn_glcc_meta%n_dim = n_dim
    plcov_mn_glcc_meta%diminfo => diminfo
    plcov_mn_glcc_meta%vartype = vartype_real !REAL variable
    plcov_mn_glcc_meta%standard_name = 'vegetation_area_fraction' !_br 08.04.14
    plcov_mn_glcc_meta%long_name = 'Plant cover minimum due to GLCC land use data'
    plcov_mn_glcc_meta%shortName = 'PLCOV_MN'
    plcov_mn_glcc_meta%stepType = 'min'
    plcov_mn_glcc_meta%units =  '1'
    plcov_mn_glcc_meta%grid_mapping = gridmp
    plcov_mn_glcc_meta%coordinates = coord
    plcov_mn_glcc_meta%data_set = 'GLCC'

    ! lai_mx_glcc_meta
    lai_mx_glcc_meta%varname = 'LAI_MX_GLCC'
    lai_mx_glcc_meta%n_dim = n_dim
    lai_mx_glcc_meta%diminfo => diminfo
    lai_mx_glcc_meta%vartype = vartype_real !REAL variable
    lai_mx_glcc_meta%standard_name = 'leaf_area_index' !_br 08.04.14
    lai_mx_glcc_meta%long_name = 'Leaf Area Index Maximum'
    lai_mx_glcc_meta%shortName = 'LAI_MX'
    lai_mx_glcc_meta%stepType = 'max'
    lai_mx_glcc_meta%units =  c_undef
    lai_mx_glcc_meta%grid_mapping = gridmp
    lai_mx_glcc_meta%coordinates = coord
    lai_mx_glcc_meta%data_set = 'GLCC'


    ! lai_mn_glcc_meta
    lai_mn_glcc_meta%varname = 'LAI_MN_GLCC'
    lai_mn_glcc_meta%n_dim = n_dim
    lai_mn_glcc_meta%diminfo => diminfo
    lai_mn_glcc_meta%vartype = vartype_real !REAL variable
    lai_mn_glcc_meta%standard_name = 'leaf_area_index' !_br 08.04.14
    lai_mn_glcc_meta%long_name = 'Leaf Area Minimum'
    lai_mn_glcc_meta%shortName = 'LAI_MN'
    lai_mn_glcc_meta%stepType = 'min'
    lai_mn_glcc_meta%units =  c_undef
    lai_mn_glcc_meta%grid_mapping = gridmp
    lai_mn_glcc_meta%coordinates = coord
    lai_mn_glcc_meta%data_set = 'GLCC'


    ! rs_min_glcc_meta
    rs_min_glcc_meta%varname = 'RSMIN_GLCC'
    rs_min_glcc_meta%n_dim = n_dim
    rs_min_glcc_meta%diminfo => diminfo
    rs_min_glcc_meta%vartype = vartype_real !REAL variable
    rs_min_glcc_meta%standard_name = c_undef !_br 08.04.14
    rs_min_glcc_meta%long_name = 'Minimal stomata resistence'
    rs_min_glcc_meta%shortName = 'RSMIN'
    rs_min_glcc_meta%stepType = 'instant'
    rs_min_glcc_meta%units =  's/m'
    rs_min_glcc_meta%grid_mapping = gridmp
    rs_min_glcc_meta%coordinates = coord
    rs_min_glcc_meta%data_set = 'GLCC'

    ! urban_glcc_meta
    urban_glcc_meta%varname = 'URBAN_GLCC'
    urban_glcc_meta%n_dim = n_dim
    urban_glcc_meta%diminfo => diminfo
    urban_glcc_meta%vartype = vartype_real !REAL variable
    urban_glcc_meta%standard_name = c_undef !_br 08.04.14
    urban_glcc_meta%long_name = 'Urban land use fraction'
    urban_glcc_meta%shortName = 'URBAN'
    urban_glcc_meta%stepType = 'instant'
    urban_glcc_meta%units =  '1'
    urban_glcc_meta%grid_mapping = gridmp
    urban_glcc_meta%coordinates = coord
    urban_glcc_meta%data_set = 'GLCC'


    ! for_d_glcc_meta
    for_d_glcc_meta%varname = 'FOR_D_GLCC'
    for_d_glcc_meta%n_dim = n_dim
    for_d_glcc_meta%diminfo => diminfo
    for_d_glcc_meta%vartype = vartype_real !REAL variable
    for_d_glcc_meta%standard_name = c_undef !_br 08.04.14
    for_d_glcc_meta%long_name = 'Fraction of deciduous forest'
    for_d_glcc_meta%shortName = 'FOR_D'
    for_d_glcc_meta%stepType = 'instant'
    for_d_glcc_meta%units =  '1'
    for_d_glcc_meta%grid_mapping = gridmp
    for_d_glcc_meta%coordinates = coord
    for_d_glcc_meta%data_set = 'GLCC'


    ! for_e_glcc_meta
    for_e_glcc_meta%varname = 'FOR_E_GLCC'
    for_e_glcc_meta%n_dim = n_dim
    for_e_glcc_meta%diminfo => diminfo
    for_e_glcc_meta%vartype = vartype_real !REAL variable
    for_e_glcc_meta%standard_name = c_undef !_br 08.04.14
    for_e_glcc_meta%long_name = 'Fraction of evergreen forest'
    for_e_glcc_meta%shortName =  'FOR_E'
    for_e_glcc_meta%stepType = 'instant'
    for_e_glcc_meta%units =  '1'
    for_e_glcc_meta%grid_mapping = gridmp
    for_e_glcc_meta%coordinates = coord
    for_e_glcc_meta%data_set = 'GLCC'


    ! emissivity_glcc_meta
    emissivity_glcc_meta%varname = 'EMISSIVITY_GLCC'
    emissivity_glcc_meta%n_dim = n_dim
    emissivity_glcc_meta%diminfo => diminfo
    emissivity_glcc_meta%vartype = vartype_real !REAL variable
    emissivity_glcc_meta%standard_name = c_undef !_br 08.04.14
    emissivity_glcc_meta%long_name = 'longwave surface emissivity'
    emissivity_glcc_meta%shortName = 'EMIS_RAD'
    emissivity_glcc_meta%stepType = 'instant'
    emissivity_glcc_meta%units =  '1'
    emissivity_glcc_meta%grid_mapping = gridmp
    emissivity_glcc_meta%coordinates = coord
    emissivity_glcc_meta%data_set = 'GLCC'



  END SUBROUTINE def_glcc_fields_meta

  !> define meta information for ECOSG target fields
  SUBROUTINE def_ecosg_fields_meta(nclass_ecosg,diminfo,coordinates,grid_mapping)
    INTEGER (KIND=i4), INTENT(IN) :: nclass_ecosg !< ECOSG has 23 classes for the land use description
    TYPE(dim_meta_info),TARGET :: diminfo(:)      !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates   !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping  !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)

    n_dim = SIZE(diminfo)


    ! set meta information for strucutre dim_ecosg_tg
    IF (ALLOCATED(dim_ecosg_tg)) DEALLOCATE(dim_ecosg_tg)
    ALLOCATE(dim_ecosg_tg(1:n_dim+1))
    SELECT CASE(n_dim)
    CASE (1)
      dim_ecosg_tg(1)%dimname = diminfo(1)%dimname
      dim_ecosg_tg(1)%dimsize = diminfo(1)%dimsize
      dim_ecosg_tg(2)%dimname = 'nclass_lu'
      dim_ecosg_tg(2)%dimsize = nclass_ecosg
    CASE (2)
      dim_ecosg_tg(1)%dimname = diminfo(1)%dimname
      dim_ecosg_tg(1)%dimsize = diminfo(1)%dimsize
      dim_ecosg_tg(2)%dimname = diminfo(2)%dimname
      dim_ecosg_tg(2)%dimsize = diminfo(2)%dimsize
      dim_ecosg_tg(3)%dimname = 'nclass_lu'
      dim_ecosg_tg(3)%dimsize = nclass_ecosg
    CASE (3)
      dim_ecosg_tg(1)%dimname = diminfo(1)%dimname
      dim_ecosg_tg(1)%dimsize = diminfo(1)%dimsize
      dim_ecosg_tg(2)%dimname = diminfo(2)%dimname
      dim_ecosg_tg(2)%dimsize = diminfo(2)%dimsize
      dim_ecosg_tg(3)%dimname = diminfo(3)%dimname
      dim_ecosg_tg(3)%dimsize = diminfo(3)%dimsize
      dim_ecosg_tg(4)%dimname = 'nclass_lu'
      dim_ecosg_tg(4)%dimsize = nclass_ecosg
    END SELECT

    ! fr_land_ecosg_meta
    fr_land_ecosg_meta%varname = 'FR_LAND'
    fr_land_ecosg_meta%n_dim = n_dim
    fr_land_ecosg_meta%diminfo => diminfo
    fr_land_ecosg_meta%vartype = vartype_real !REAL variable
    fr_land_ecosg_meta%standard_name = 'land_area_fraction' !_br 08.04.14
    fr_land_ecosg_meta%long_name = 'Fraction land ECOSG Data'
    fr_land_ecosg_meta%shortName = 'FR_LAND'
    fr_land_ecosg_meta%stepType = 'instant'
    fr_land_ecosg_meta%units =  c_undef
    fr_land_ecosg_meta%grid_mapping = gridmp
    fr_land_ecosg_meta%coordinates = coord
    fr_land_ecosg_meta%data_set = 'ECOSG'

   ! ecosg_tot_npixel_meta
    ecosg_tot_npixel_meta%varname = 'LU_TOT_NPIXEL'
    ecosg_tot_npixel_meta%n_dim = n_dim
    ecosg_tot_npixel_meta%diminfo => dim_ecosg_tg
    ecosg_tot_npixel_meta%vartype = vartype_int !INTEGER variable
    ecosg_tot_npixel_meta%standard_name = c_undef !_br 08.04.14
    ecosg_tot_npixel_meta%long_name = 'number of raw data pixel in target grid element'
    ecosg_tot_npixel_meta%shortName = c_undef
    ecosg_tot_npixel_meta%stepType = 'instant'
    ecosg_tot_npixel_meta%units = c_undef
    ecosg_tot_npixel_meta%grid_mapping = gridmp
    ecosg_tot_npixel_meta%coordinates = coord
    ecosg_tot_npixel_meta%data_set = 'ECOSG'

    ! ecosg_class_fraction_meta
    ecosg_class_fraction_meta%varname = 'LU_CLASS_FRACTION'
    ecosg_class_fraction_meta%n_dim = n_dim + 1
    ecosg_class_fraction_meta%diminfo => dim_ecosg_tg
    ecosg_class_fraction_meta%vartype = vartype_real !REAL variable
    ecosg_class_fraction_meta%standard_name = c_undef !_br 08.04.14
    ecosg_class_fraction_meta%long_name = 'Fraction of ECOSG land use classes in target grid element'
    ecosg_class_fraction_meta%shortName = c_undef
    ecosg_class_fraction_meta%stepType = 'instant'
    ecosg_class_fraction_meta%units =  c_undef
    ecosg_class_fraction_meta%grid_mapping = gridmp
    ecosg_class_fraction_meta%coordinates = coord
    ecosg_class_fraction_meta%data_set = 'ECOSG'

    ! ecosg_class_npixel_meta
    ecosg_class_npixel_meta%varname = 'LU_CLASS_NPIXEL'
    ecosg_class_npixel_meta%n_dim = n_dim + 1
    ecosg_class_npixel_meta%diminfo => dim_ecosg_tg
    ecosg_class_npixel_meta%vartype = vartype_int !INTEGER variable
    ecosg_class_npixel_meta%standard_name = c_undef !_br 08.04.14
    ecosg_class_npixel_meta%long_name = 'number of pixels of ECOSG land use classes in target grid element'
    ecosg_class_npixel_meta%shortName = c_undef
    ecosg_class_npixel_meta%stepType = 'instant'
    ecosg_class_npixel_meta%units = c_undef
    ecosg_class_npixel_meta%grid_mapping = gridmp
    ecosg_class_npixel_meta%coordinates = coord
    ecosg_class_npixel_meta%data_set = 'ECOSG'

    ! ice_ecosg_meta
    ice_ecosg_meta%varname = 'ICE'
    ice_ecosg_meta%n_dim = n_dim
    ice_ecosg_meta%diminfo => diminfo
    ice_ecosg_meta%vartype = vartype_real !REAL variable
    ice_ecosg_meta%standard_name = c_undef !_br 08.04.14
    ice_ecosg_meta%long_name = 'Ice fraction due to ECOSG Data'
    ice_ecosg_meta%shortName = c_undef
    ice_ecosg_meta%stepType = 'instant'
    ice_ecosg_meta%units =  '1'
    ice_ecosg_meta%grid_mapping = gridmp
    ice_ecosg_meta%coordinates = coord
    ice_ecosg_meta%data_set = 'ECOSG'

    ! z0_ecosg_meta
    z0_ecosg_meta%varname = 'Z0'
    z0_ecosg_meta%n_dim = n_dim
    z0_ecosg_meta%diminfo => diminfo
    z0_ecosg_meta%vartype = vartype_real !REAL variable
    z0_ecosg_meta%standard_name = 'surface_roughness_length' !_br 08.04.14
    z0_ecosg_meta%long_name =  'Roughness length z0 due to ECOSG land use data'
    z0_ecosg_meta%shortName = 'Z0_LU'
    z0_ecosg_meta%stepType = 'instant'
    z0_ecosg_meta%units = 'm'
    z0_ecosg_meta%grid_mapping = gridmp
    z0_ecosg_meta%coordinates = coord
    z0_ecosg_meta%data_set = 'ECOSG'

    ! root_ecosg_meta
    root_ecosg_meta%varname = 'ROOTDP'
    root_ecosg_meta%n_dim = n_dim
    root_ecosg_meta%diminfo => diminfo
    root_ecosg_meta%vartype = vartype_real !REAL variable
    root_ecosg_meta%standard_name = c_undef !_br 08.04.14
    root_ecosg_meta%long_name = 'Root depth due to ECOSG land use data'
    root_ecosg_meta%shortName = 'ROOTDP'
    root_ecosg_meta%stepType = 'instant'
    root_ecosg_meta%units =  'm'
    root_ecosg_meta%grid_mapping = gridmp
    root_ecosg_meta%coordinates = coord
    root_ecosg_meta%data_set = 'ECOSG'

    ! plcov_mx_ecosg_meta
    plcov_mx_ecosg_meta%varname = 'PLCOV_MX'
    plcov_mx_ecosg_meta%n_dim = n_dim
    plcov_mx_ecosg_meta%diminfo => diminfo
    plcov_mx_ecosg_meta%vartype = vartype_real !REAL variable
    plcov_mx_ecosg_meta%standard_name = 'vegetation_area_fraction' !_br 08.04.14
    plcov_mx_ecosg_meta%long_name = 'Plant cover maximum due to ECOSG land use data'
    plcov_mx_ecosg_meta%shortName = 'PLCOV_MX'
    plcov_mx_ecosg_meta%stepType = 'max'
    plcov_mx_ecosg_meta%units =  '1'
    plcov_mx_ecosg_meta%grid_mapping = gridmp
    plcov_mx_ecosg_meta%coordinates = coord
    plcov_mx_ecosg_meta%data_set = 'ECOSG'


    ! plcov_mn_ecosg_meta
    plcov_mn_ecosg_meta%varname = 'PLCOV_MN'
    plcov_mn_ecosg_meta%n_dim = n_dim
    plcov_mn_ecosg_meta%diminfo => diminfo
    plcov_mn_ecosg_meta%vartype = vartype_real !REAL variable
    plcov_mn_ecosg_meta%standard_name = 'vegetation_area_fraction' !_br 08.04.14
    plcov_mn_ecosg_meta%long_name = 'Plant cover minimum due to ECOSG land use data'
    plcov_mn_ecosg_meta%shortName = 'PLCOV_MN'
    plcov_mn_ecosg_meta%stepType = 'min'
    plcov_mn_ecosg_meta%units =  '1'
    plcov_mn_ecosg_meta%grid_mapping = gridmp
    plcov_mn_ecosg_meta%coordinates = coord
    plcov_mn_ecosg_meta%data_set = 'ECOSG'

    ! lai_mx_ecosg_meta
    lai_mx_ecosg_meta%varname = 'LAI_MX'
    lai_mx_ecosg_meta%n_dim = n_dim
    lai_mx_ecosg_meta%diminfo => diminfo
    lai_mx_ecosg_meta%vartype = vartype_real !REAL variable
    lai_mx_ecosg_meta%standard_name = 'leaf_area_index' !_br 08.04.14
    lai_mx_ecosg_meta%long_name = 'Leaf Area Index Maximum'
    lai_mx_ecosg_meta%shortName = 'LAI_MX'
    lai_mx_ecosg_meta%stepType = 'max'
    lai_mx_ecosg_meta%units =  c_undef
    lai_mx_ecosg_meta%grid_mapping = gridmp
    lai_mx_ecosg_meta%coordinates = coord
    lai_mx_ecosg_meta%data_set = 'ECOSG'


    ! lai_mn_ecosg_meta
    lai_mn_ecosg_meta%varname = 'LAI_MN'
    lai_mn_ecosg_meta%n_dim = n_dim
    lai_mn_ecosg_meta%diminfo => diminfo
    lai_mn_ecosg_meta%vartype = vartype_real !REAL variable
    lai_mn_ecosg_meta%standard_name = 'leaf_area_index' !_br 08.04.14
    lai_mn_ecosg_meta%long_name = 'Leaf Area Minimum'
    lai_mn_ecosg_meta%shortName = 'LAI_MN'
    lai_mn_ecosg_meta%stepType = 'min'
    lai_mn_ecosg_meta%units =  c_undef
    lai_mn_ecosg_meta%grid_mapping = gridmp
    lai_mn_ecosg_meta%coordinates = coord
    lai_mn_ecosg_meta%data_set = 'ECOSG'


    ! rs_min_ecosg_meta
    rs_min_ecosg_meta%varname = 'RSMIN'
    rs_min_ecosg_meta%n_dim = n_dim
    rs_min_ecosg_meta%diminfo => diminfo
    rs_min_ecosg_meta%vartype = vartype_real !REAL variable
    rs_min_ecosg_meta%standard_name = c_undef !_br 08.04.14
    rs_min_ecosg_meta%long_name = 'Minimal stomata resistence'
    rs_min_ecosg_meta%shortName = 'RSMIN'
    rs_min_ecosg_meta%stepType = 'instant'
    rs_min_ecosg_meta%units =  's/m'
    rs_min_ecosg_meta%grid_mapping = gridmp
    rs_min_ecosg_meta%coordinates = coord
    rs_min_ecosg_meta%data_set = 'ECOSG'

    ! urban_ecosg_meta
    urban_ecosg_meta%varname = 'URBAN'
    urban_ecosg_meta%n_dim = n_dim
    urban_ecosg_meta%diminfo => diminfo
    urban_ecosg_meta%vartype = vartype_real !REAL variable
    urban_ecosg_meta%standard_name = c_undef !_br 08.04.14
    urban_ecosg_meta%long_name = 'Urban land use fraction'
    urban_ecosg_meta%shortName = 'URBAN'
    urban_ecosg_meta%stepType = 'instant'
    urban_ecosg_meta%units =  '1'
    urban_ecosg_meta%grid_mapping = gridmp
    urban_ecosg_meta%coordinates = coord
    urban_ecosg_meta%data_set = 'ECOSG'


    ! for_d_ecosg_meta
    for_d_ecosg_meta%varname = 'FOR_D'
    for_d_ecosg_meta%n_dim = n_dim
    for_d_ecosg_meta%diminfo => diminfo
    for_d_ecosg_meta%vartype = vartype_real !REAL variable
    for_d_ecosg_meta%standard_name = c_undef !_br 08.04.14
    for_d_ecosg_meta%long_name = 'Fraction of deciduous forest'
    for_d_ecosg_meta%shortName = 'FOR_D'
    for_d_ecosg_meta%stepType = 'instant'
    for_d_ecosg_meta%units =  '1'
    for_d_ecosg_meta%grid_mapping = gridmp
    for_d_ecosg_meta%coordinates = coord
    for_d_ecosg_meta%data_set = 'ECOSG'


    ! for_e_ecosg_meta
    for_e_ecosg_meta%varname = 'FOR_E'
    for_e_ecosg_meta%n_dim = n_dim
    for_e_ecosg_meta%diminfo => diminfo
    for_e_ecosg_meta%vartype = vartype_real !REAL variable
    for_e_ecosg_meta%standard_name = c_undef !_br 08.04.14
    for_e_ecosg_meta%long_name = 'Fraction of evergreen forest'
    for_e_ecosg_meta%shortName =  'FOR_E'
    for_e_ecosg_meta%stepType = 'instant'
    for_e_ecosg_meta%units =  '1'
    for_e_ecosg_meta%grid_mapping = gridmp
    for_e_ecosg_meta%coordinates = coord
    for_e_ecosg_meta%data_set = 'ECOSG'

! skinc_lu_meta
    skinc_ecosg_meta%varname = 'SKC'
    skinc_ecosg_meta%n_dim = n_dim
    skinc_ecosg_meta%diminfo => diminfo
    skinc_ecosg_meta%vartype = vartype_real !REAL variable
    skinc_ecosg_meta%standard_name = 'skin_conductivity'
    skinc_ecosg_meta%long_name = 'Skin conductivity'
    skinc_ecosg_meta%shortName = 'SKC'
    skinc_ecosg_meta%units =  c_undef
    skinc_ecosg_meta%grid_mapping = gridmp
    skinc_ecosg_meta%coordinates = coord
    skinc_ecosg_meta%data_set = 'ECOSG'

    ! emissivity_ecosg_meta
    emissivity_ecosg_meta%varname = 'EMIS_RAD'
    emissivity_ecosg_meta%n_dim = n_dim
    emissivity_ecosg_meta%diminfo => diminfo
    emissivity_ecosg_meta%vartype = vartype_real !REAL variable
    emissivity_ecosg_meta%standard_name = c_undef !_br 08.04.14
    emissivity_ecosg_meta%long_name = 'longwave surface emissivity'
    emissivity_ecosg_meta%shortName = 'EMIS_RAD'
    emissivity_ecosg_meta%stepType = 'instant'
    emissivity_ecosg_meta%units =  c_undef
    emissivity_ecosg_meta%grid_mapping = gridmp
    emissivity_ecosg_meta%coordinates = coord
    emissivity_ecosg_meta%data_set = 'ECOSG'

  END SUBROUTINE def_ecosg_fields_meta



  !> define meta information for  landuse target fields
  SUBROUTINE def_lu_fields_meta(nclass_lu,diminfo,lu_dataset,coordinates,grid_mapping)
    INTEGER (KIND=i4), INTENT(IN) :: nclass_lu !< Number of classes for the land use description
    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (LEN=*), OPTIONAL :: lu_dataset !< name of landuse data set
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping


    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord
    CHARACTER (len=80) :: dataset

    gridmp = c_undef
    coord = c_undef
    dataset = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)

    n_dim = SIZE(diminfo)

    ! set meta information for strucutre dim_ndvi_tg
    IF (ALLOCATED(dim_lu_tg)) DEALLOCATE(dim_lu_tg)
    ALLOCATE(dim_lu_tg(1:n_dim+1))
    SELECT CASE(n_dim)
    CASE (1)
      dim_lu_tg(1)%dimname = diminfo(1)%dimname
      dim_lu_tg(1)%dimsize = diminfo(1)%dimsize
      dim_lu_tg(2)%dimname = 'nclass_lu'
      dim_lu_tg(2)%dimsize = nclass_lu
    CASE (2)
      dim_lu_tg(1)%dimname = diminfo(1)%dimname
      dim_lu_tg(1)%dimsize = diminfo(1)%dimsize
      dim_lu_tg(2)%dimname = diminfo(2)%dimname
      dim_lu_tg(2)%dimsize = diminfo(2)%dimsize
      dim_lu_tg(3)%dimname = 'nclass_lu'
      dim_lu_tg(3)%dimsize = nclass_lu
    CASE (3)
      dim_lu_tg(1)%dimname = diminfo(1)%dimname
      dim_lu_tg(1)%dimsize = diminfo(1)%dimsize
      dim_lu_tg(2)%dimname = diminfo(2)%dimname
      dim_lu_tg(2)%dimsize = diminfo(2)%dimsize
      dim_lu_tg(3)%dimname = diminfo(3)%dimname
      dim_lu_tg(3)%dimsize = diminfo(3)%dimsize
      dim_lu_tg(4)%dimname = 'nclass_lu'
      dim_lu_tg(4)%dimsize = nclass_lu
    END SELECT

    ! fr_land_lu_meta
    fr_land_lu_meta%varname = 'FR_LAND'
    fr_land_lu_meta%n_dim = n_dim
    fr_land_lu_meta%diminfo => diminfo
    fr_land_lu_meta%vartype = vartype_real !REAL variable
    fr_land_lu_meta%standard_name = 'land_area_fraction'
    fr_land_lu_meta%long_name = 'Fraction land'
    fr_land_lu_meta%shortName = 'FR_LAND'
    fr_land_lu_meta%stepType = 'instant'
    fr_land_lu_meta%units =  c_undef
    fr_land_lu_meta%grid_mapping = gridmp
    fr_land_lu_meta%coordinates = coord
    fr_land_lu_meta%data_set = dataset



   ! lu_tot_npixel_meta
    lu_tot_npixel_meta%varname = 'LU_TOT_NPIXEL'
    lu_tot_npixel_meta%n_dim = n_dim
    lu_tot_npixel_meta%diminfo => diminfo
    lu_tot_npixel_meta%vartype = vartype_int !INTEGER variable
    lu_tot_npixel_meta%standard_name = c_undef !_br 08.04.14
    lu_tot_npixel_meta%long_name = 'number of raw data pixel in target grid element'
    lu_tot_npixel_meta%shortName = c_undef
    lu_tot_npixel_meta%stepType = 'instant'
    lu_tot_npixel_meta%units = c_undef
    lu_tot_npixel_meta%grid_mapping = gridmp
    lu_tot_npixel_meta%coordinates = coord
    lu_tot_npixel_meta%data_set = dataset


    ! lu_class_fraction_meta
    lu_class_fraction_meta%varname = 'LU_CLASS_FRACTION'
    lu_class_fraction_meta%n_dim = n_dim + 1
    lu_class_fraction_meta%diminfo => dim_lu_tg
    lu_class_fraction_meta%vartype = vartype_real !REAL variable
    lu_class_fraction_meta%standard_name = c_undef !_br 08.04.14
    IF ( PRESENT(lu_dataset) ) THEN
      lu_class_fraction_meta%long_name = 'Fraction of '//TRIM(lu_dataset)// &
                                         ' land use classes in target grid element'
    ELSE
      lu_class_fraction_meta%long_name = ''
    ENDIF
    lu_class_fraction_meta%shortName = c_undef
    lu_class_fraction_meta%stepType = 'instant'
    lu_class_fraction_meta%units =  c_undef
    lu_class_fraction_meta%grid_mapping = gridmp
    lu_class_fraction_meta%coordinates = coord
    lu_class_fraction_meta%data_set = dataset

    ! lu_class_npixel_meta
    lu_class_npixel_meta%varname = 'LU_CLASS_NPIXEL'
    lu_class_npixel_meta%n_dim = n_dim + 1
    lu_class_npixel_meta%diminfo => dim_lu_tg
    lu_class_npixel_meta%vartype = vartype_int !INTEGER variable
    lu_class_npixel_meta%standard_name = c_undef !_br 08.04.14
    IF ( PRESENT(lu_dataset) ) THEN
      lu_class_npixel_meta%long_name = 'number of pixels of '//TRIM(lu_dataset)// &
                                       ' land use classes in target grid element'
    ELSE
      lu_class_fraction_meta%long_name = ''
    ENDIF
    lu_class_npixel_meta%shortName = c_undef
    lu_class_npixel_meta%stepType = 'instant'
    lu_class_npixel_meta%units = c_undef
    lu_class_npixel_meta%grid_mapping = gridmp
    lu_class_npixel_meta%coordinates = coord
    lu_class_npixel_meta%data_set = dataset

    ! ice_lu_meta
    ice_lu_meta%varname = 'ICE'
    ice_lu_meta%n_dim = n_dim
    ice_lu_meta%diminfo => diminfo
    ice_lu_meta%vartype = vartype_real !REAL variable
    ice_lu_meta%standard_name = c_undef !_br 08.04.14
    IF ( PRESENT(lu_dataset) ) THEN
      ice_lu_meta%long_name = 'Ice fraction due to '//TRIM(lu_dataset)//' Data'
    ELSE
      ice_lu_meta%long_name = ''
    ENDIF
    ice_lu_meta%shortName = c_undef
    ice_lu_meta%stepType = 'instant'
    ice_lu_meta%units =  c_undef
    ice_lu_meta%grid_mapping = gridmp
    ice_lu_meta%coordinates = coord
    ice_lu_meta%data_set = dataset

    ! z0_lu_meta
    z0_lu_meta%varname = 'Z0'
    z0_lu_meta%n_dim = n_dim
    z0_lu_meta%diminfo => diminfo
    z0_lu_meta%vartype = vartype_real !REAL variable
    z0_lu_meta%standard_name = 'surface_roughness_length'
    z0_lu_meta%long_name = 'Roughness length'
    z0_lu_meta%shortName = 'Z0'
    z0_lu_meta%stepType = 'instant'
    z0_lu_meta%units =  c_undef
    z0_lu_meta%grid_mapping = gridmp
    z0_lu_meta%coordinates = coord
    z0_lu_meta%data_set = dataset

    ! root_lu_meta
    root_lu_meta%varname = 'ROOTDP'
    root_lu_meta%n_dim = n_dim
    root_lu_meta%diminfo => diminfo
    root_lu_meta%vartype = vartype_real !REAL variable
    root_lu_meta%standard_name = 'root_depth'
    root_lu_meta%long_name = 'Root depth'
    root_lu_meta%shortName = 'ROOTDP'
    root_lu_meta%stepType = 'instant'
    root_lu_meta%units =  'm'
    root_lu_meta%grid_mapping = gridmp
    root_lu_meta%coordinates = coord
    root_lu_meta%data_set = dataset

    ! plcov_mx_lu_meta
    plcov_mx_lu_meta%varname = 'PLCOV_MX'
    plcov_mx_lu_meta%n_dim = n_dim
    plcov_mx_lu_meta%diminfo => diminfo
    plcov_mx_lu_meta%vartype = vartype_real !REAL variable
    plcov_mx_lu_meta%standard_name = 'vegetation_area_fraction' !_br 08.04.14
    IF ( PRESENT(lu_dataset) ) THEN
      plcov_mx_lu_meta%long_name = 'Plant cover maximum due to '//TRIM(lu_dataset)//' land use data'
    ELSE
      plcov_mx_lu_meta%long_name = ''
    ENDIF
    plcov_mx_lu_meta%shortName = 'PLCOV_MX'
    plcov_mx_lu_meta%stepType = 'max'
    plcov_mx_lu_meta%units =  c_undef
    plcov_mx_lu_meta%grid_mapping = gridmp
    plcov_mx_lu_meta%coordinates = coord
    plcov_mx_lu_meta%data_set = dataset


    ! plcov_mn_lu_meta
    plcov_mn_lu_meta%varname = 'PLCOV_MN'
    plcov_mn_lu_meta%n_dim = n_dim
    plcov_mn_lu_meta%diminfo => diminfo
    plcov_mn_lu_meta%vartype = vartype_real !REAL variable
    plcov_mn_lu_meta%standard_name = 'vegetation_area_fraction' !_br 08.04.14
    IF ( PRESENT(lu_dataset) ) THEN
      plcov_mn_lu_meta%long_name = 'Plant cover minimum due to '//TRIM(lu_dataset)//' land use data'
    ELSE
      plcov_mn_lu_meta%long_name = ''
    ENDIF
    plcov_mn_lu_meta%shortName = 'PLCOV_MN'
    plcov_mn_lu_meta%stepType = 'min'
    plcov_mn_lu_meta%units =  c_undef
    plcov_mn_lu_meta%grid_mapping = gridmp
    plcov_mn_lu_meta%coordinates = coord
    plcov_mn_lu_meta%data_set = dataset

    ! lai_mx_lu_meta
    lai_mx_lu_meta%varname = 'LAI_MX'
    lai_mx_lu_meta%n_dim = n_dim
    lai_mx_lu_meta%diminfo => diminfo
    lai_mx_lu_meta%vartype = vartype_real !REAL variable
    lai_mx_lu_meta%standard_name = 'leaf_area_index' !_br 08.04.14
    lai_mx_lu_meta%long_name = 'Leaf Area Index Maximum'
    lai_mx_lu_meta%shortName = 'LAI_MX'
    lai_mx_lu_meta%stepType = 'max'
    lai_mx_lu_meta%units =  c_undef
    lai_mx_lu_meta%grid_mapping = gridmp
    lai_mx_lu_meta%coordinates = coord
    lai_mx_lu_meta%data_set = dataset


    ! lai_mn_lu_meta
    lai_mn_lu_meta%varname = 'LAI_MN'
    lai_mn_lu_meta%n_dim = n_dim
    lai_mn_lu_meta%diminfo => diminfo
    lai_mn_lu_meta%vartype = vartype_real !REAL variable
    lai_mn_lu_meta%standard_name = 'leaf_area_index' !_br 08.04.14
    lai_mn_lu_meta%long_name = 'Leaf Area Minimum'
    lai_mn_lu_meta%shortName = 'LAI_MN'
    lai_mn_lu_meta%stepType = 'min'
    lai_mn_lu_meta%units =  c_undef
    lai_mn_lu_meta%grid_mapping =gridmp
    lai_mn_lu_meta%coordinates = coord
    lai_mn_lu_meta%data_set = dataset


    ! rs_min_lu_meta
    rs_min_lu_meta%varname = 'RSMIN'
    rs_min_lu_meta%n_dim = n_dim
    rs_min_lu_meta%diminfo => diminfo
    rs_min_lu_meta%vartype = vartype_real !REAL variable
    rs_min_lu_meta%standard_name = c_undef !_br 08.04.14
    rs_min_lu_meta%long_name = 'Minimal stomata resistence'
    rs_min_lu_meta%shortName = 'RSMIN'
    rs_min_lu_meta%stepType = 'instant'
    rs_min_lu_meta%units =  's/m'
    rs_min_lu_meta%grid_mapping = gridmp
    rs_min_lu_meta%coordinates = coord
    rs_min_lu_meta%data_set = dataset

    ! urban_lu_meta
    urban_lu_meta%varname = 'URBAN'
    urban_lu_meta%n_dim = n_dim
    urban_lu_meta%diminfo => diminfo
    urban_lu_meta%vartype = vartype_real !REAL variable
    urban_lu_meta%standard_name = c_undef !_br 08.04.14
    urban_lu_meta%long_name = 'urban area fraction'
    urban_lu_meta%shortName = 'URBAN'
    urban_lu_meta%stepType = 'instant'
    urban_lu_meta%units =  c_undef
    urban_lu_meta%grid_mapping = gridmp
    urban_lu_meta%coordinates = coord
    urban_lu_meta%data_set = dataset


    ! for_d_lu_meta
    for_d_lu_meta%varname = 'FOR_D'
    for_d_lu_meta%n_dim = n_dim
    for_d_lu_meta%diminfo => diminfo
    for_d_lu_meta%vartype = vartype_real !REAL variable
    for_d_lu_meta%standard_name = c_undef !_br 08.04.14
    for_d_lu_meta%long_name = 'Fraction of deciduous forest'
    for_d_lu_meta%shortName = 'FOR_D'
    for_d_lu_meta%stepType = 'instant'
    for_d_lu_meta%units =  c_undef
    for_d_lu_meta%grid_mapping = gridmp
    for_d_lu_meta%coordinates = coord
    for_d_lu_meta%data_set = dataset


    ! for_e_lu_meta
    for_e_lu_meta%varname = 'FOR_E'
    for_e_lu_meta%n_dim = n_dim
    for_e_lu_meta%diminfo => diminfo
    for_e_lu_meta%vartype = vartype_real !REAL variable
    for_e_lu_meta%standard_name = c_undef !_br 08.04.14
    for_e_lu_meta%long_name = 'Fraction of evergreen forest'
    for_e_lu_meta%shortName = 'FOR_E'
    for_e_lu_meta%stepType = 'instant'
    for_e_lu_meta%units =  c_undef
    for_e_lu_meta%grid_mapping = gridmp
    for_e_lu_meta%coordinates = coord
    for_e_lu_meta%data_set = dataset

! skinc_lu_meta
    skinc_lu_meta%varname = 'SKC'
    skinc_lu_meta%n_dim = n_dim
    skinc_lu_meta%diminfo => diminfo
    skinc_lu_meta%vartype = vartype_real !REAL variable
    skinc_lu_meta%standard_name = 'skin_conductivity'
    skinc_lu_meta%long_name = 'Skin conductivity'
    skinc_lu_meta%shortName = 'SKC'
    skinc_lu_meta%units =  c_undef
    skinc_lu_meta%grid_mapping = gridmp
    skinc_lu_meta%coordinates = coord
    skinc_lu_meta%data_set = dataset

    ! emissivity_lu_meta
    emissivity_lu_meta%varname = 'EMIS_RAD'
    emissivity_lu_meta%n_dim = n_dim
    emissivity_lu_meta%diminfo => diminfo
    emissivity_lu_meta%vartype = vartype_real !REAL variable
    emissivity_lu_meta%standard_name = c_undef !_br 08.04.14'
    emissivity_lu_meta%long_name = 'longwave surface emissivity'
    emissivity_lu_meta%shortName = 'EMIS_RAD'
    emissivity_lu_meta%stepType = 'instant'
    emissivity_lu_meta%units =  c_undef
    emissivity_lu_meta%grid_mapping = gridmp
    emissivity_lu_meta%coordinates = coord
    emissivity_lu_meta%data_set = dataset

    ! fr_ocean_lu_meta
    fr_ocean_lu_meta%varname = 'FR_OCEAN'
    fr_ocean_lu_meta%n_dim = n_dim
    fr_ocean_lu_meta%diminfo => diminfo
    fr_ocean_lu_meta%vartype = vartype_real !REAL variable
    fr_ocean_lu_meta%standard_name = c_undef !_br 08.04.14
    fr_ocean_lu_meta%long_name = 'Fraction ocean'
    fr_ocean_lu_meta%shortName = 'FR_OCEAN'
    fr_ocean_lu_meta%stepType = 'instant'
    fr_ocean_lu_meta%units =  c_undef
    fr_ocean_lu_meta%grid_mapping = gridmp
    fr_ocean_lu_meta%coordinates = coord
    fr_ocean_lu_meta%data_set = dataset

  END SUBROUTINE def_lu_fields_meta



  SUBROUTINE def_ecoclimap_fields_meta(ntime,nclass_lu,diminfo,coordinates,grid_mapping)
    INTEGER (KIND=i4), INTENT(IN) :: ntime !< number of times
    INTEGER (KIND=i4), INTENT(IN) :: nclass_lu !<
    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)

    n_dim = SIZE(diminfo)

    ! set meta information for strucutre dim_ndvi_tg
    IF (ALLOCATED(dim_ecoclimap_tg)) DEALLOCATE(dim_ecoclimap_tg)
    ALLOCATE(dim_ecoclimap_tg(1:n_dim+2))
    SELECT CASE(n_dim)
    CASE (1)
      dim_ecoclimap_tg(1)%dimname = diminfo(1)%dimname
      dim_ecoclimap_tg(1)%dimsize = diminfo(1)%dimsize
      dim_ecoclimap_tg(2)%dimname = 'nclass_lu'
      dim_ecoclimap_tg(2)%dimsize = nclass_lu
      dim_ecoclimap_tg(3)%dimname = 'time'
      dim_ecoclimap_tg(3)%dimsize = ntime

    CASE (2)
      dim_ecoclimap_tg(1)%dimname = diminfo(1)%dimname
      dim_ecoclimap_tg(1)%dimsize = diminfo(1)%dimsize
      dim_ecoclimap_tg(2)%dimname = diminfo(2)%dimname
      dim_ecoclimap_tg(2)%dimsize = diminfo(2)%dimsize
      dim_ecoclimap_tg(3)%dimname = 'nclass_lu'
      dim_ecoclimap_tg(3)%dimsize = nclass_lu
      dim_ecoclimap_tg(4)%dimname = 'time'
      dim_ecoclimap_tg(4)%dimsize = ntime

    CASE (3)
      dim_ecoclimap_tg(1)%dimname = diminfo(1)%dimname
      dim_ecoclimap_tg(1)%dimsize = diminfo(1)%dimsize
      dim_ecoclimap_tg(2)%dimname = diminfo(2)%dimname
      dim_ecoclimap_tg(2)%dimsize = diminfo(2)%dimsize
      dim_ecoclimap_tg(3)%dimname = diminfo(3)%dimname
      dim_ecoclimap_tg(3)%dimsize = diminfo(3)%dimsize
      dim_ecoclimap_tg(4)%dimname = 'nclass_lu'
      dim_ecoclimap_tg(4)%dimsize = nclass_lu
      dim_ecoclimap_tg(5)%dimname = 'time'
      dim_ecoclimap_tg(5)%dimsize = ntime

    END SELECT

       ! set meta information for strucutre
    IF (ALLOCATED(dim_ecoclimap_tg2)) DEALLOCATE(dim_ecoclimap_tg2)
    ALLOCATE(dim_ecoclimap_tg2(1:n_dim+1))
    SELECT CASE(n_dim)
    CASE (1)
      dim_ecoclimap_tg2(1)%dimname = diminfo(1)%dimname
      dim_ecoclimap_tg2(1)%dimsize = diminfo(1)%dimsize
      dim_ecoclimap_tg2(2)%dimname = 'time'
      dim_ecoclimap_tg2(2)%dimsize = ntime
    CASE (2)
      dim_ecoclimap_tg2(1)%dimname = diminfo(1)%dimname
      dim_ecoclimap_tg2(1)%dimsize = diminfo(1)%dimsize
      dim_ecoclimap_tg2(2)%dimname = diminfo(2)%dimname
      dim_ecoclimap_tg2(2)%dimsize = diminfo(2)%dimsize
      dim_ecoclimap_tg2(3)%dimname = 'time'
      dim_ecoclimap_tg2(3)%dimsize = ntime
    CASE (3)
      dim_ecoclimap_tg2(1)%dimname = diminfo(1)%dimname
      dim_ecoclimap_tg2(1)%dimsize = diminfo(1)%dimsize
      dim_ecoclimap_tg2(2)%dimname = diminfo(2)%dimname
      dim_ecoclimap_tg2(2)%dimsize = diminfo(2)%dimsize
      dim_ecoclimap_tg2(3)%dimname = diminfo(3)%dimname
      dim_ecoclimap_tg2(3)%dimsize = diminfo(3)%dimsize
      dim_ecoclimap_tg2(4)%dimname = 'time'
      dim_ecoclimap_tg2(4)%dimsize = ntime
    END SELECT

    ! fr_land_lu_meta
    fr_land_lu_meta%varname = 'FR_LAND'
    fr_land_lu_meta%n_dim = n_dim
    fr_land_lu_meta%diminfo => diminfo
    fr_land_lu_meta%vartype = vartype_real !REAL variable
    fr_land_lu_meta%standard_name = 'land_area_fraction'
    fr_land_lu_meta%long_name = 'Fraction land'
    fr_land_lu_meta%shortName = 'FR_LAND'
    fr_land_lu_meta%stepType = 'instant'
    fr_land_lu_meta%units =  '1'
    fr_land_lu_meta%grid_mapping = gridmp
    fr_land_lu_meta%coordinates = coord
    fr_land_lu_meta%data_set = 'ECOCLIMAP'

   ! lu_tot_npixel_meta
    lu_tot_npixel_meta%varname = 'LU_TOT_NPIXEL'
    lu_tot_npixel_meta%n_dim = n_dim
    lu_tot_npixel_meta%diminfo => diminfo
    lu_tot_npixel_meta%vartype = vartype_int !INTEGER variable
    lu_tot_npixel_meta%standard_name = c_undef
    lu_tot_npixel_meta%long_name = 'number of raw data pixel in target grid element'
    lu_tot_npixel_meta%shortName = c_undef
    lu_tot_npixel_meta%stepType = 'instant'
    lu_tot_npixel_meta%units = '1'
    lu_tot_npixel_meta%grid_mapping = gridmp
    lu_tot_npixel_meta%coordinates = coord
    lu_tot_npixel_meta%data_set = 'ECOCLIMAP'

    ! lu_class_fraction_meta
    lu_class_fraction_meta%varname = 'LU_CLASS_FRACTION'
    lu_class_fraction_meta%n_dim = n_dim + 1
    lu_class_fraction_meta%diminfo => dim_ecoclimap_tg
    lu_class_fraction_meta%vartype = vartype_real !REAL variable
    lu_class_fraction_meta%standard_name = c_undef
    lu_class_fraction_meta%long_name = 'Fraction of land use classes in target grid element'
    lu_class_fraction_meta%shortName = c_undef
    lu_class_fraction_meta%stepType = 'instant'
    lu_class_fraction_meta%units =  '1'
    lu_class_fraction_meta%grid_mapping = gridmp
    lu_class_fraction_meta%coordinates = coord
    lu_class_fraction_meta%data_set = 'ECOCLIMAP'

    ! lu_class_npixel_meta
    lu_class_npixel_meta%varname = 'LU_CLASS_NPIXEL'
    lu_class_npixel_meta%n_dim = n_dim + 1
    lu_class_npixel_meta%diminfo => dim_ecoclimap_tg
    lu_class_npixel_meta%vartype = vartype_int !INTEGER variable
    lu_class_npixel_meta%standard_name = c_undef
    lu_class_npixel_meta%long_name = 'number of pixels of land use classes in target grid element'
    lu_class_npixel_meta%shortName = c_undef
    lu_class_npixel_meta%stepType = 'instant'
    lu_class_npixel_meta%units = '1'
    lu_class_npixel_meta%grid_mapping = gridmp
    lu_class_npixel_meta%coordinates = coord
    lu_class_npixel_meta%data_set = 'ECOCLIMAP'

    ! ice_lu_meta
    ice_lu_meta%varname = 'ICE'
    ice_lu_meta%n_dim = n_dim
    ice_lu_meta%diminfo => diminfo
    ice_lu_meta%vartype = vartype_real !REAL variable
    ice_lu_meta%standard_name = c_undef !_br 08.04.14
    ice_lu_meta%long_name = 'Ice fraction'
    ice_lu_meta%shortName = c_undef
    ice_lu_meta%stepType = 'instant'
    ice_lu_meta%units =  '1'
    ice_lu_meta%grid_mapping = gridmp
    ice_lu_meta%coordinates = coord
    ice_lu_meta%data_set = 'ECOCLIMAP'

    ! z0_lu_meta
    z012_tot_meta%varname = 'Z012'
    z012_tot_meta%n_dim = n_dim + 1
    z012_tot_meta%diminfo => dim_ecoclimap_tg2
    z012_tot_meta%vartype = vartype_real !REAL variable
    z012_tot_meta%standard_name = 'surface_roughness_length ' !_br 08.04.14
    z012_tot_meta%long_name = 'Monthly roughness length'
    z012_tot_meta%shortName = 'Z0_VEG'
    z012_tot_meta%stepType = 'avg'
    z012_tot_meta%units =  'm'
    z012_tot_meta%grid_mapping = gridmp
    z012_tot_meta%coordinates = coord
    z012_tot_meta%data_set = 'ECOCLIMAP'

    ! z0_lu_meta
    z012_lu_meta%varname = 'Z0_VEG'
    z012_lu_meta%n_dim = n_dim + 1
    z012_lu_meta%diminfo => dim_ecoclimap_tg2
    z012_lu_meta%vartype = vartype_real !REAL variable
    z012_lu_meta%standard_name = 'surface_roughness_length due to vegetation' !_br 08.04.14
    z012_lu_meta%long_name = 'Monthly roughness length due to vegetation'
    z012_lu_meta%shortName = 'Z0_VEG'
    z012_lu_meta%stepType = 'avg'
    z012_lu_meta%units =  'm'
    z012_lu_meta%grid_mapping = gridmp
    z012_lu_meta%coordinates = coord
    z012_lu_meta%data_set = 'ECOCLIMAP'

    ! root_lu_meta
    root_lu_meta%varname = 'ROOTDP'
    root_lu_meta%n_dim = n_dim
    root_lu_meta%diminfo => diminfo
    root_lu_meta%vartype = vartype_real !REAL variable
    root_lu_meta%standard_name = 'root_depth'
    root_lu_meta%long_name = 'Root depth'
    root_lu_meta%shortName = 'ROOTDP'
    root_lu_meta%stepType = 'instant'
    root_lu_meta%units =  'm'
    root_lu_meta%grid_mapping = gridmp
    root_lu_meta%coordinates = coord
    root_lu_meta%data_set = 'ECOCLIMAP'

   ! plcov12_lu_meta
    plcov12_lu_meta%varname = 'PLCOV12'
    plcov12_lu_meta%n_dim = n_dim + 1
    plcov12_lu_meta%diminfo =>dim_ecoclimap_tg2
    plcov12_lu_meta%vartype = vartype_real !REAL variable
    plcov12_lu_meta%standard_name = 'vegetation_area_fraction' !_br 08.04.14
    plcov12_lu_meta%long_name = 'monthly plant cover'
    plcov12_lu_meta%shortName = 'PLCOV12'
    plcov12_lu_meta%stepType = 'avg'
    plcov12_lu_meta%units =  '1'
    plcov12_lu_meta%grid_mapping = gridmp
    plcov12_lu_meta%coordinates = coord
    plcov12_lu_meta%data_set = 'ECOCLIMAP'

    ! lai12_lu_meta
    lai12_lu_meta%varname = 'LAI12'
    lai12_lu_meta%n_dim = n_dim + 1
    lai12_lu_meta%diminfo => dim_ecoclimap_tg2
    lai12_lu_meta%vartype = vartype_real !REAL variable
    lai12_lu_meta%standard_name = 'leaf_area_index' !_br 08.04.14
    lai12_lu_meta%long_name = 'Monthly Leaf Area Index'
    lai12_lu_meta%shortName = 'LAI12'
    lai12_lu_meta%stepType = 'avg'
    lai12_lu_meta%units =  '1'
    lai12_lu_meta%grid_mapping =gridmp
    lai12_lu_meta%coordinates = coord
    lai12_lu_meta%data_set = 'ECOCLIMAP'

    ! rs_min_lu_meta
    rs_min_lu_meta%varname = 'RSMIN'
    rs_min_lu_meta%n_dim = n_dim
    rs_min_lu_meta%diminfo => diminfo
    rs_min_lu_meta%vartype = vartype_real !REAL variable
    rs_min_lu_meta%standard_name = c_undef
    rs_min_lu_meta%long_name = c_undef
    rs_min_lu_meta%shortName = 'prs_min'
    rs_min_lu_meta%stepType = 'instant'
    rs_min_lu_meta%units =  's/m'
    rs_min_lu_meta%grid_mapping = gridmp
    rs_min_lu_meta%coordinates = coord
    rs_min_lu_meta%data_set = 'ECOCLIMAP'

    ! urban_lu_meta
    urban_lu_meta%varname = 'URBAN'
    urban_lu_meta%n_dim = n_dim
    urban_lu_meta%diminfo => diminfo
    urban_lu_meta%vartype = vartype_real !REAL variable
    urban_lu_meta%standard_name = c_undef
    urban_lu_meta%long_name = 'urban area fraction'
    urban_lu_meta%shortName = 'URBAN'
    urban_lu_meta%stepType = 'instant'
    urban_lu_meta%units =  '1'
    urban_lu_meta%grid_mapping = gridmp
    urban_lu_meta%coordinates = coord
    urban_lu_meta%data_set = 'ECOCLIMAP'

    ! for_d_lu_meta
    for_d_lu_meta%varname = 'FOR_D'
    for_d_lu_meta%n_dim = n_dim
    for_d_lu_meta%diminfo => diminfo
    for_d_lu_meta%vartype = vartype_real !REAL variable
    for_d_lu_meta%standard_name = c_undef
    for_d_lu_meta%long_name = 'Fraction of deciduous forest'
    for_d_lu_meta%shortName = 'FOR_D'
    for_d_lu_meta%stepType = 'instant'
    for_d_lu_meta%units =  '1'
    for_d_lu_meta%grid_mapping = gridmp
    for_d_lu_meta%coordinates = coord
    for_d_lu_meta%data_set = 'ECOCLIMAP'

    ! for_e_lu_meta
    for_e_lu_meta%varname = 'FOR_E'
    for_e_lu_meta%n_dim = n_dim
    for_e_lu_meta%diminfo => diminfo
    for_e_lu_meta%vartype = vartype_real !REAL variable
    for_e_lu_meta%standard_name = c_undef
    for_e_lu_meta%long_name = 'Fraction of evergreen forest'
    for_e_lu_meta%shortName = 'FOR_E'
    for_e_lu_meta%stepType = 'instant'
    for_e_lu_meta%units =  '1'
    for_e_lu_meta%grid_mapping = gridmp
    for_e_lu_meta%coordinates = coord
    for_e_lu_meta%data_set = 'ECOCLIMAP'

    ! skinc_lu_meta
    skinc_lu_meta%varname = 'SKC'
    skinc_lu_meta%n_dim = n_dim
    skinc_lu_meta%diminfo => diminfo
    skinc_lu_meta%vartype = vartype_real !REAL variable
    skinc_lu_meta%standard_name = 'skin_conductivity'
    skinc_lu_meta%long_name = 'Skin conductivity'
    skinc_lu_meta%shortName = 'SKC'
    skinc_lu_meta%stepType = 'instant'
    skinc_lu_meta%units =  c_undef
    skinc_lu_meta%grid_mapping = gridmp
    skinc_lu_meta%coordinates = coord
    skinc_lu_meta%data_set = 'ECOCLIMAP'

    ! emissivity_lu_meta
    emissivity_lu_meta%varname = 'EMIS_RAD'
    emissivity_lu_meta%n_dim = n_dim
    emissivity_lu_meta%diminfo => diminfo
    emissivity_lu_meta%vartype = vartype_real !REAL variable
    emissivity_lu_meta%standard_name = c_undef !_br 08.04.14
    emissivity_lu_meta%long_name = 'longwave surface emissivity'
    emissivity_lu_meta%shortName = 'EMIS_RAD'
    emissivity_lu_meta%stepType = 'instant'
    emissivity_lu_meta%units =  c_undef
    emissivity_lu_meta%grid_mapping = gridmp
    emissivity_lu_meta%coordinates = coord
    emissivity_lu_meta%data_set = 'ECOCLIMAP'

    ! fr_ocean_lu_meta
    fr_ocean_lu_meta%varname = 'FR_OCEAN'
    fr_ocean_lu_meta%n_dim = n_dim
    fr_ocean_lu_meta%diminfo => diminfo
    fr_ocean_lu_meta%vartype = vartype_real !REAL variable
    fr_ocean_lu_meta%standard_name = c_undef !_br 08.04.14
    fr_ocean_lu_meta%long_name = 'Fraction ocean'
    fr_ocean_lu_meta%shortName = 'FR_OCEAN'
    fr_ocean_lu_meta%stepType = 'instant'
    fr_ocean_lu_meta%units =  '1'
    fr_ocean_lu_meta%grid_mapping = gridmp
    fr_ocean_lu_meta%coordinates = coord
    fr_ocean_lu_meta%data_set = 'ECOCLIMAP'

    ! plcov_mx_lu_meta
    plcov_mx_lu_meta%varname = 'PLCOV_MX'
    plcov_mx_lu_meta%n_dim = n_dim
    plcov_mx_lu_meta%diminfo => diminfo
    plcov_mx_lu_meta%vartype = vartype_real !REAL variable
    plcov_mx_lu_meta%long_name = 'vegetation_area_fraction'
    plcov_mx_lu_meta%standard_name = 'Plant cover maximum'
    plcov_mx_lu_meta%shortName = 'PLCOV_MX'
    plcov_mx_lu_meta%stepType = 'max'
    plcov_mx_lu_meta%units =  c_undef
    plcov_mx_lu_meta%grid_mapping = gridmp
    plcov_mx_lu_meta%coordinates = coord
    plcov_mx_lu_meta%data_set = 'ECOCLIMAP'

    ! plcov_mn_lu_meta
    plcov_mn_lu_meta%varname = 'PLCOV_MN'
    plcov_mn_lu_meta%n_dim = n_dim
    plcov_mn_lu_meta%diminfo => diminfo
    plcov_mn_lu_meta%vartype = vartype_real !REAL variable
    plcov_mn_lu_meta%standard_name = 'vegetation_area_fraction'
    plcov_mn_lu_meta%long_name = 'Plant cover minimum'
    plcov_mn_lu_meta%shortName = 'PLCOV_MN'
    plcov_mn_lu_meta%stepType = 'min'
    plcov_mn_lu_meta%units =  c_undef
    plcov_mn_lu_meta%grid_mapping = gridmp
    plcov_mn_lu_meta%coordinates = coord
    plcov_mn_lu_meta%data_set = 'ECOCLIMAP'

    ! lai_mx_lu_meta
    lai_mx_lu_meta%varname = 'LAI_MX'
    lai_mx_lu_meta%n_dim = n_dim
    lai_mx_lu_meta%diminfo => diminfo
    lai_mx_lu_meta%vartype = vartype_real !REAL variable
    lai_mx_lu_meta%standard_name = 'leaf_area_index'
    lai_mx_lu_meta%long_name = 'Leaf Area Index Maximum'
    lai_mx_lu_meta%shortName = 'LAI_MX'
    lai_mx_lu_meta%stepType = 'max'
    lai_mx_lu_meta%units =  c_undef
    lai_mx_lu_meta%grid_mapping = gridmp
    lai_mx_lu_meta%coordinates = coord
    lai_mx_lu_meta%data_set = 'ECOCLIMAP'

    ! lai_mn_lu_meta
    lai_mn_lu_meta%varname = 'LAI_MN'
    lai_mn_lu_meta%n_dim = n_dim
    lai_mn_lu_meta%diminfo => diminfo
    lai_mn_lu_meta%vartype = vartype_real !REAL variable
    lai_mn_lu_meta%standard_name = 'leaf_area_index'
    lai_mn_lu_meta%long_name = 'Leaf Area Minimum'
    lai_mn_lu_meta%shortName = 'LAI_MN'
    lai_mn_lu_meta%stepType = 'min'
    lai_mn_lu_meta%units =  c_undef
    lai_mn_lu_meta%grid_mapping =gridmp
    lai_mn_lu_meta%coordinates = coord
    lai_mn_lu_meta%data_set = 'ECOCLIMAP'

    ! z0_lu_meta
    z0_lu_meta%varname = 'Z0'
    z0_lu_meta%n_dim = n_dim
    z0_lu_meta%diminfo => diminfo
    z0_lu_meta%vartype = vartype_real !REAL variable
    z0_lu_meta%standard_name = 'surface_roughness_length'
    z0_lu_meta%long_name = 'Roughness length'
    z0_lu_meta%shortName = 'Z0'
    z0_lu_meta%stepType = 'instant'
    z0_lu_meta%units =  c_undef
    z0_lu_meta%grid_mapping = gridmp
    z0_lu_meta%coordinates = coord
    z0_lu_meta%data_set = 'ECOCLIMAP'

  END SUBROUTINE def_ecoclimap_fields_meta


  !> define meta information for flake data for netcdf output
  SUBROUTINE def_flake_fields_meta(diminfo,coordinates,grid_mapping)
    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord

    gridmp = c_undef
    coord = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)




    lake_depth_meta%varname = 'DEPTH_LK'
    lake_depth_meta%n_dim = n_dim
    lake_depth_meta%diminfo => diminfo
    lake_depth_meta%vartype = vartype_real !REAL variable
    lake_depth_meta%standard_name = 'sea_floor_depth_below_sea_level' !_br 08.04.14
    lake_depth_meta%long_name = 'Lake depth'
    lake_depth_meta%shortName = 'DEPTH_LK'
    lake_depth_meta%stepType = 'instant'
    lake_depth_meta%units = 'm'
    lake_depth_meta%grid_mapping = gridmp
    lake_depth_meta%coordinates = coord
    lake_depth_meta%data_set = 'DWD/RSHU/MeteoFrance'


    fr_lake_meta%varname = 'FR_LAKE'
    fr_lake_meta%n_dim = n_dim
    fr_lake_meta%diminfo => diminfo
    fr_lake_meta%vartype = vartype_real !REAL variable
    fr_lake_meta%standard_name = c_undef !_br 08.04.14
    fr_lake_meta%long_name = 'fraction lake'
    fr_lake_meta%shortName = 'FR_LAKE'
    fr_lake_meta%stepType = 'instant'
    fr_lake_meta%units = '1'
    fr_lake_meta%grid_mapping = gridmp
    fr_lake_meta%coordinates = coord
    fr_lake_meta%data_set = 'DWD/RSHU/MeteoFrance'

    ! flake_tot_npixel_meta
    flake_tot_npixel_meta%varname = 'FLAKE_TOT_NPIXEL'
    flake_tot_npixel_meta%n_dim = n_dim
    flake_tot_npixel_meta%diminfo => diminfo
    flake_tot_npixel_meta%vartype = vartype_int !INTEGER variable
    flake_tot_npixel_meta%standard_name = c_undef !_br 08.04.14
    flake_tot_npixel_meta%long_name = 'number of raw data pixel in target grid element'
    flake_tot_npixel_meta%shortName = c_undef
    flake_tot_npixel_meta%stepType = 'instant'
    flake_tot_npixel_meta%units = c_undef
    flake_tot_npixel_meta%grid_mapping = gridmp
    flake_tot_npixel_meta%coordinates = coord
    flake_tot_npixel_meta%data_set = 'DWD/RSHU/MeteoFrance'



  END SUBROUTINE def_flake_fields_meta


  !> define meta information for flake data for netcdf output
  SUBROUTINE def_lsm_fields_meta(diminfo,coordinates,grid_mapping)
    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord
    CHARACTER (len=80) :: dataset

    gridmp = c_undef
    coord = c_undef
    dataset = c_undef

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)



      ! fr_land_mask_meta
    fr_land_mask_meta%varname = 'FR_LAND'
    fr_land_mask_meta%n_dim = n_dim
    fr_land_mask_meta%diminfo => diminfo
    fr_land_mask_meta%vartype = vartype_real !REAL variable
    fr_land_mask_meta%standard_name = 'land_area_fraction'
    fr_land_mask_meta%long_name = 'Fraction land'
    fr_land_mask_meta%shortName = 'FR_LAND'
    fr_land_mask_meta%stepType = 'instant'
    fr_land_mask_meta%units =  c_undef
    fr_land_mask_meta%grid_mapping = gridmp
    fr_land_mask_meta%coordinates = coord
    fr_land_mask_meta%data_set = dataset



  END SUBROUTINE def_lsm_fields_meta



  !> define meta information for target fields derived from GLOBE data
  SUBROUTINE def_topo_meta(diminfo,itopo_type,igrid_type,coordinates,grid_mapping,diminfohor)

    TYPE(dim_meta_info),TARGET :: diminfo(:)     !< pointer to dimensions of variable
    INTEGER (KIND=i4), INTENT(IN):: itopo_type, &!< defines the desired topography (ASTER or GLOBE)
         &                          igrid_type   !< COSMO or ICON grid
    CHARACTER (len=80), OPTIONAL :: coordinates  !< netcdf attribute coordinates
    CHARACTER (len=80), OPTIONAL :: grid_mapping !< netcdf attribute grid mapping
    TYPE(dim_meta_info),TARGET, OPTIONAL :: diminfohor(:)     !< pointer to dimensions of variable

    ! local variables
    INTEGER  :: n_dim      !< number of dimensions
    CHARACTER (len=80) :: gridmp
    CHARACTER (len=80) :: coord, coordhor, dataset
    INTEGER  :: n_dimhor   !< number of dimensions
    INTEGER (KIND=i4), PARAMETER  :: topo_aster = 2, &
         &                           topo_gl = 1, &
         &                           topo_merit = 3, &
         &                           igrid_icon = 1, &
         &                           igrid_cosmo = 2


    gridmp = c_undef
    coord = c_undef
    coordhor = c_undef
    dataset = c_undef
    SELECT CASE(itopo_type)
      CASE(topo_aster)
        dataset = 'ASTER'
      CASE(topo_gl)
        dataset = 'GLOBE'
      CASE(topo_merit)
        dataset = 'MERIT'
    END SELECT

    IF (PRESENT(grid_mapping)) gridmp = TRIM(grid_mapping)
    IF (PRESENT(coordinates)) coord = TRIM(coordinates)
    n_dim = SIZE(diminfo)
    IF (PRESENT(diminfohor)) n_dimhor = SIZE(diminfohor)
    IF (PRESENT(diminfohor) .AND. PRESENT(coordinates)) coordhor =  TRIM(coordinates)// ' nhori'

    ! set meta information for strucutre dim_buffer_cell
    dim_buffer_cell = dim_3d_tg

    hh_topo_meta%varname = 'HSURF'
    hh_topo_meta%n_dim = n_dim
    hh_topo_meta%diminfo => diminfo
    hh_topo_meta%vartype = vartype_real !REAL variable
    hh_topo_meta%standard_name = 'surface_altitude' !_br 08.04.14
    hh_topo_meta%long_name = 'geometric height of the earths surface above sea level'
    hh_topo_meta%shortName = 'HSURF'
    hh_topo_meta%stepType = 'instant'
    hh_topo_meta%units = 'm'
    hh_topo_meta%grid_mapping = gridmp
    hh_topo_meta%coordinates = coord
    hh_topo_meta%data_set = dataset

    hh_topo_min_meta%varname = 'SSO_OROMIN'
    hh_topo_min_meta%n_dim = n_dim
    hh_topo_min_meta%diminfo => diminfo
    hh_topo_min_meta%vartype = vartype_real !REAL variable
    hh_topo_min_meta%standard_name = 'minimum_contributing_surface_height'
    hh_topo_min_meta%long_name = 'minimum geometric height of contributing raw data height points'
    hh_topo_min_meta%shortName = 'SSO_OROMIN'
    hh_topo_min_meta%units = 'm'
    hh_topo_min_meta%grid_mapping = gridmp
    hh_topo_min_meta%coordinates = coord

    hh_topo_max_meta%varname = 'SSO_OROMAX'
    hh_topo_max_meta%n_dim = n_dim
    hh_topo_max_meta%diminfo => diminfo
    hh_topo_max_meta%vartype = vartype_real !REAL variable
    hh_topo_max_meta%standard_name = 'maximum_contributing_surface_height'
    hh_topo_max_meta%long_name = 'maximum geometric height of contributing raw data height points'
    hh_topo_max_meta%shortName = 'SSO_OROMAX'
    hh_topo_max_meta%units = 'm'
    hh_topo_max_meta%grid_mapping = gridmp
    hh_topo_max_meta%coordinates = coord

    hh_fis_meta%varname = 'FIS'
    hh_fis_meta%n_dim = n_dim
    hh_fis_meta%diminfo => diminfo
    hh_fis_meta%vartype = vartype_real !REAL variable
    hh_fis_meta%standard_name = 'surface_geopotential_height' !_br 08.04.14
    hh_fis_meta%long_name = 'Geopotential (S)'
    hh_fis_meta%shortName = 'FIS'
    hh_fis_meta%stepType = 'instant'
    hh_fis_meta%units = 'm**2 s**-2'
    hh_fis_meta%grid_mapping = gridmp
    hh_fis_meta%coordinates = coord
    hh_fis_meta%data_set = dataset

    stdh_topo_meta%varname = 'SSO_STDH'
    stdh_topo_meta%n_dim = n_dim
    stdh_topo_meta%diminfo => diminfo
    stdh_topo_meta%vartype = vartype_real !REAL variable
    stdh_topo_meta%standard_name = c_undef !_br 08.04.14
    stdh_topo_meta%long_name = 'standard deviation of subgrid scale orography'
    stdh_topo_meta%shortName = 'SSO_STDH'
    stdh_topo_meta%stepType = 'instant'
    stdh_topo_meta%units = 'm'
    stdh_topo_meta%grid_mapping = gridmp
    stdh_topo_meta%coordinates = coord
    stdh_topo_meta%data_set = dataset

    theta_topo_meta%varname = 'SSO_THETA'
    theta_topo_meta%n_dim = n_dim
    theta_topo_meta%diminfo => diminfo
    theta_topo_meta%vartype = vartype_real !REAL variable
    theta_topo_meta%standard_name = c_undef !_br 08.04.14
    theta_topo_meta%long_name = 'Angle of sub-gridscale orography'
    theta_topo_meta%shortName = 'SSO_THETA'
    theta_topo_meta%stepType = 'instant'
    theta_topo_meta%units = 'rad'
    theta_topo_meta%grid_mapping = gridmp
    theta_topo_meta%coordinates = coord
    theta_topo_meta%data_set = dataset

    aniso_topo_meta%varname = 'SSO_GAMMA'
    aniso_topo_meta%n_dim = n_dim
    aniso_topo_meta%diminfo => diminfo
    aniso_topo_meta%vartype = vartype_real !REAL variable
    aniso_topo_meta%standard_name = c_undef !_br 08.04.14
    aniso_topo_meta%long_name = 'anisotropy of sub-gridscale orography'
    aniso_topo_meta%shortName = 'SSO_GAMMA'
    aniso_topo_meta%stepType = 'instant'
    aniso_topo_meta%units =  c_undef
    aniso_topo_meta%grid_mapping = gridmp
    aniso_topo_meta%coordinates = coord
    aniso_topo_meta%data_set = dataset

    slope_topo_meta%varname = 'SSO_SIGMA'
    slope_topo_meta%n_dim = n_dim
    slope_topo_meta%diminfo => diminfo
    slope_topo_meta%vartype = vartype_real !REAL variable
    slope_topo_meta%standard_name = c_undef !_br 08.04.14
    slope_topo_meta%long_name = 'Slope of sub-gridscale orography'
    slope_topo_meta%shortName = 'SSO_SIGMA'
    slope_topo_meta%stepType = 'instant'
    slope_topo_meta%units = c_undef
    slope_topo_meta%grid_mapping = gridmp
    slope_topo_meta%coordinates = coord
    slope_topo_meta%data_set = dataset

    fr_land_topo_meta%varname = 'FR_LAND_TOPO'
    fr_land_topo_meta%n_dim = n_dim
    fr_land_topo_meta%diminfo => diminfo
    fr_land_topo_meta%vartype = vartype_real !REAL variable
    fr_land_topo_meta%standard_name = 'land_area_fraction' !_br 08.04.14
    SELECT CASE(itopo_type)
      CASE(topo_aster)
        fr_land_topo_meta%long_name = 'fraction land due to ASTER data'
      CASE(topo_gl)
        fr_land_topo_meta%long_name = 'fraction land due to GLOBE data'
      CASE(topo_merit)
        fr_land_topo_meta%long_name = 'fraction land due to MERIT data'
      END SELECT
    fr_land_topo_meta%shortName = 'FR_LAND'
    fr_land_topo_meta%stepType = 'instant'
    fr_land_topo_meta%units =  c_undef
    fr_land_topo_meta%grid_mapping = gridmp
    fr_land_topo_meta%coordinates = coord
    fr_land_topo_meta%data_set = dataset

    ! z0_topo_meta
    z0_topo_meta%varname = 'Z0_TOPO'
    z0_topo_meta%n_dim = n_dim
    z0_topo_meta%diminfo => diminfo
    z0_topo_meta%vartype = vartype_real !REAL variable
    z0_topo_meta%standard_name = 'surface_roughness_length' !_br 08.04.14
    z0_topo_meta%long_name = 'Roughness length'
    z0_topo_meta%shortName = 'Z0'
    z0_topo_meta%stepType = 'instant'
    z0_topo_meta%units =  c_undef
    z0_topo_meta%grid_mapping = gridmp
    z0_topo_meta%coordinates = coord
    z0_topo_meta%data_set = dataset

    !lradtopo parameters
    slope_asp_topo_meta%varname = 'SLO_ASP'
    slope_asp_topo_meta%n_dim = n_dim
    slope_asp_topo_meta%diminfo => diminfo
    slope_asp_topo_meta%vartype = vartype_real !REAL variable
    slope_asp_topo_meta%standard_name = c_undef !_br 08.04.14
    slope_asp_topo_meta%long_name = 'slope aspect - topography'
    slope_asp_topo_meta%shortName = 'SLO_ASP'
    slope_asp_topo_meta%stepType = 'instant'
    slope_asp_topo_meta%units = 'rad'
    slope_asp_topo_meta%grid_mapping = gridmp
    slope_asp_topo_meta%coordinates = coord
    slope_asp_topo_meta%data_set = dataset

    slope_ang_topo_meta%varname = 'SLO_ANG'
    slope_ang_topo_meta%n_dim = n_dim
    slope_ang_topo_meta%diminfo => diminfo
    slope_ang_topo_meta%vartype = vartype_real !REAL variable
    slope_ang_topo_meta%standard_name = c_undef !_br 08.04.14
    slope_ang_topo_meta%long_name = 'slope angle - topography'
    slope_ang_topo_meta%shortName = 'SLO_ANG'
    slope_ang_topo_meta%stepType = 'instant'
    slope_ang_topo_meta%units = 'rad'
    slope_ang_topo_meta%grid_mapping = gridmp
    slope_ang_topo_meta%coordinates = coord
    slope_ang_topo_meta%data_set = dataset


    IF (PRESENT(diminfohor)) THEN
      horizon_topo_meta%varname = 'HORIZON'
      horizon_topo_meta%n_dim = n_dimhor
      horizon_topo_meta%diminfo => diminfohor
      horizon_topo_meta%vartype = vartype_real !REAL variable
      horizon_topo_meta%standard_name = c_undef !_br 08.04.14
      horizon_topo_meta%long_name = 'horizon angle - topography'
      horizon_topo_meta%shortName = 'HORIZON'
      horizon_topo_meta%stepType = 'instant'
      horizon_topo_meta%units = 'deg'
      horizon_topo_meta%grid_mapping = gridmp
      horizon_topo_meta%coordinates = coordhor
      horizon_topo_meta%data_set = dataset
    ENDIF

    skyview_topo_meta%varname = 'SKYVIEW'
    skyview_topo_meta%n_dim = n_dim
    skyview_topo_meta%diminfo => diminfo
    skyview_topo_meta%vartype = vartype_real !REAL variable
    skyview_topo_meta%standard_name = c_undef !_br 08.04.14
    SELECT CASE(igrid_type)
      CASE(igrid_cosmo)
        skyview_topo_meta%long_name = 'sky-view factor'
      CASE(igrid_icon)
        IF (itype_scaling == 0) THEN
          skyview_topo_meta%long_name = 'geometric sky-view factor'
        ELSEIF(itype_scaling == 1) THEN
          skyview_topo_meta%long_name = 'geometric sky-view factor scaled with sinus(horizon)'
        ELSEIF(itype_scaling == 2) THEN
          skyview_topo_meta%long_name = 'geometric sky-view factor scaled with sinus(horizon)**2'
        ENDIF
      END SELECT

    skyview_topo_meta%shortName = 'SKYVIEW'
    skyview_topo_meta%stepType = 'instant'
    skyview_topo_meta%units = '-'
    skyview_topo_meta%grid_mapping = gridmp
    skyview_topo_meta%coordinates = coord
    skyview_topo_meta%data_set = dataset

    sgsl_meta%varname = 'S_ORO'
    sgsl_meta%n_dim = n_dim
    sgsl_meta%diminfo => diminfo
    sgsl_meta%vartype = vartype_real !REAL variable
    sgsl_meta%standard_name = c_undef
    sgsl_meta%long_name = 'subgrid-scale slope based on raw DEM data'
    sgsl_meta%shortName = 'S_ORO'
    sgsl_meta%stepType = 'instant'
    sgsl_meta%units = c_undef
    sgsl_meta%grid_mapping = gridmp
    sgsl_meta%coordinates = coord
    sgsl_meta%data_set = dataset

  END SUBROUTINE def_topo_meta

  !> define meta information for target fields defined on vertices derived from GLOBE data
  SUBROUTINE def_topo_vertex_meta(nvertex)
  INTEGER, INTENT(IN) :: nvertex !< total number of vertices
  CHARACTER (len=80)  :: dataset

  dataset = c_undef

    ! set meta information for strucutre dim_buffer_vertex
    dim_buffer_vertex(1)%dimname = 'vertex'
    dim_buffer_vertex(1)%dimsize = nvertex
    dim_buffer_vertex(2)%dimname = 'je_v'
    dim_buffer_vertex(2)%dimsize = 1
    dim_buffer_vertex(3)%dimname = 'ke_v'
    dim_buffer_vertex(3)%dimsize = 1

    hh_vert_meta%varname = 'topography_v'
    hh_vert_meta%n_dim = 3
    hh_vert_meta%diminfo => dim_buffer_vertex
    hh_vert_meta%vartype = vartype_real !REAL variable
    hh_vert_meta%standard_name = 'surface_altitude' !_br 08.04.14
    hh_vert_meta%long_name = 'topographic height at cell vertices'
    hh_vert_meta%shortName = 'HSURF'
    hh_vert_meta%stepType = 'instant'
    hh_vert_meta%units = 'm'
    hh_vert_meta%grid_mapping = c_undef
    hh_vert_meta%coordinates = c_undef
    hh_vert_meta%data_set = dataset

    npixel_vert_meta%varname = 'npixel_vert'
    npixel_vert_meta%n_dim = 3
    npixel_vert_meta%diminfo => dim_buffer_vertex
    npixel_vert_meta%vartype = vartype_int !INTEGER variable
    npixel_vert_meta%standard_name = c_undef !_br 08.04.14
    npixel_vert_meta%long_name = 'number of raw data pixel'
    npixel_vert_meta%shortName = c_undef
    npixel_vert_meta%stepType = 'instant'
    npixel_vert_meta%units = c_undef
    npixel_vert_meta%grid_mapping = c_undef
    npixel_vert_meta%coordinates = c_undef
    npixel_vert_meta%data_set = dataset
  END SUBROUTINE def_topo_vertex_meta


  !> set projection information for netcdf output for the COSMO grid
  SUBROUTINE set_nc_grid_def_cosmo(cosmo_grid,grid_mapping)
    TYPE(rotated_lonlat_grid), INTENT(IN) :: cosmo_grid !< structure which contains the definition of the COSMO grid
    CHARACTER (len=80), INTENT(IN):: grid_mapping !< netcdf attribute grid mapping
    INTEGER :: errorcode

    IF (.NOT.ALLOCATED(nc_grid_def_cosmo%map_param)) THEN
      ALLOCATE(nc_grid_def_cosmo%map_param(1:3),STAT=errorcode)
      IF (errorcode /= 0 ) CALL logging%error('Cant nc_grid_def_cosmo%map_param',__FILE__,__LINE__)
    ENDIF
    nc_grid_def_cosmo%grid_mapping_varname =  TRIM(grid_mapping)
    nc_grid_def_cosmo%grid_mapping_name%attname='grid_mapping_name'
    nc_grid_def_cosmo%grid_mapping_name%attributetext = 'rotated_latitude_longitude'
! netcdf attribute with grid mapping name according to cf, e.g. http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.4/apf.html
    nc_grid_def_cosmo%n_r_att = 3 ! number of projection parameters
    nc_grid_def_cosmo%map_param(1)%attname     = 'grid_north_pole_longitude'
    nc_grid_def_cosmo%map_param(1)%att_value_r = REAL(cosmo_grid%pollon) ! this is for type conversion
    nc_grid_def_cosmo%map_param(2)%attname     = 'grid_north_pole_latitude'
    nc_grid_def_cosmo%map_param(2)%att_value_r = REAL(cosmo_grid%pollat) ! this is for type conversion
    nc_grid_def_cosmo%map_param(3)%attname     = 'north_pole_grid_longitude'
    nc_grid_def_cosmo%map_param(3)%att_value_r = REAL(cosmo_grid%polgam) ! this is for type conversion

  END SUBROUTINE set_nc_grid_def_cosmo


  !> set projection information for netcdf output for the ICON grid
  SUBROUTINE set_nc_grid_def_icon(grid_mapping)
    USE mo_physical_constants, ONLY: re !< av. radius of the earth [m]
    CHARACTER (len=80), INTENT(IN):: grid_mapping !< netcdf attribute grid mapping
    INTEGER :: errorcode

    ALLOCATE(nc_grid_def_icon%map_param(1:2),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant nc_grid_def_icon%map_param',__FILE__,__LINE__)
    nc_grid_def_icon%grid_mapping_varname = TRIM(grid_mapping)
    nc_grid_def_icon%grid_mapping_name%attname='grid_mapping_name'
    nc_grid_def_icon%grid_mapping_name%attributetext = 'latitude_longitude'
! netcdf attribute with grid mapping name according to cf, e.g. http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.4/apf.html
    nc_grid_def_icon%n_r_att = 2 ! number of projection parameters
    nc_grid_def_icon%map_param(1)%attname     = 'semi_major_axis'
    nc_grid_def_icon%map_param(1)%att_value_r = REAL(re)  ! type conversion to standard real
    nc_grid_def_icon%map_param(2)%attname     = 'inverse_flattening'
    nc_grid_def_icon%map_param(2)%att_value_r = 0.

 END SUBROUTINE set_nc_grid_def_icon

END MODULE mo_var_meta_data
