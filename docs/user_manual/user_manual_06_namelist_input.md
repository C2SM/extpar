# Namelist Input for Extpar {#namelist_input_for_extpar}
Extpar uses 3 types of namelists in order to determine in which way data
is processed.

-   Fortran namelists (INPUT_)

-   Python dictionaries (input_in namelist.py)

-   Fortran namelists written by Python scripts

Whereas for the Fortran namelists and the Python dictionaries the user
can specify parameters and filenames, the Fortran namelists generated
during runtime by the Python scripts do not allow any user interaction.

## Namelist files {#namelist_input_for_extpar_namelist_files}

  **Namelist file**    **Purpose**                                                      **Made by script**             **Used by program**
  -------------------- ---------------------------------------------------------------- ------------------------------ -----------------------------
  INPUT_grid_org     define target grid type                                          runscript                      extpar_consistency_check
                                                                                                                       extpar_aot_to_buffer
                                                                                                                       extpar_landuse_to_buffer
                                                                                                                       extpar_topo_to_buffer
                                                                                                                       extpar_cru_to_buffer
                                                                                                                       extpar_ndvi_to_buffer
                                                                                                                       extpar_soil_to_buffer
                                                                                                                       extpar_flake_to_buffer
                                                                                                                       extpar_isa_to_buffer
                                                                                                                       extpar_ahf_to_buffer
                                                                                                                       extpar_emiss_to_buffer
                                                                                                                       extpar_hwsdART_to_buffer
  INPUT_COSMO_GRID   define target domain for COSMO grid                              runscript                      extpar_consistency_check
                                                                                                                       extpar_aot_to_buffer
                                                                                                                       extpar_landuse_to_buffer
                                                                                                                       extpar_topo_to_buffer
                                                                                                                       extpar_cru_to_buffer
                                                                                                                       extpar_ndvi_to_buffer
                                                                                                                       extpar_soil_to_buffer
                                                                                                                       extpar_flake_to_buffer
                                                                                                                       extpar_isa_to_buffer
                                                                                                                       extpar_ahf_to_buffer
                                                                                                                       extpar_emiss_to_buffer
                                                                                                                       extpar_hwsdART_to_buffer
  INPUT_ICON_GRID    define target domain for ICON grid                               runscript                      extpar_consistency_check
                                                                                                                       extpar_aot_to_buffer
                                                                                                                       extpar_landuse_to_buffer
                                                                                                                       extpar_topo_to_buffer
                                                                                                                       extpar_cru_to_buffer
                                                                                                                       extpar_ndvi_to_buffer
                                                                                                                       extpar_soil_to_buffer
                                                                                                                       extpar_flake_to_buffer
                                                                                                                       extpar_isa_to_buffer
                                                                                                                       extpar_ahf_to_buffer
                                                                                                                       extpar_emiss_to_buffer
  INPUT_ORO           settings for orography data                                      runscript                      extpar_topo_to_buffer
  INPUT_OROSMOOTH     settings for orography smoothing                                 runscript                      extpar_topo_to_buffer
  INPUT_RADTOPO       settings for generating topographical shading fields             runscript                      extpar_topo_to_buffer
  INPUT_SCALE_SEP    settings to control scale separation for SSO an Z0 calculation   runscript                      extpar_topo_to_buffer
  INPUT_LU            settings for landuse data                                        runscript                      extpar_landuse_to_buffer
  INPUT_AOT           settings for aerosol data                                        runscript                      extpar_aot_to_buffer
  INPUT_TCLIM         settings for temperature data                                    extpar_cru_-to_buffer       extpar_consistency_check
  INPUT_NDVI          settings for NDVI data                                           extpar_ndvi_-to_buffer      extpar_consistency_check
  INPUT_SOIL          settings for soil data                                           runscript                      extpar_soil_to_buffer
  INPUT_FLAKE         settings for lake data                                           runscript                      extpar_flake_to_buffer
  INPUT_ALB           settings for albedo data                                         extpar_albedo_-to_buffer    extpar_consistency_check
  INPUT_ISA           settings for fraction of impervious surface area data            extpar_isa_to_buffer        extpar_consistency_check
  INPUT_AHF           settings for anthropogenic heat flux data                        extpar_ahf_to_buffer        extpar_consistency_check
  INPUT_EMISS         settings for emissivity data                                     extpar_emiss_-to_buffer     extpar_consistency_check
  INPUT_hwsdART       settings for HWSD USDA data                                      extpar_hwsdART_-to_buffer   
  INPUT_edgar         settings for EDGAR data                                          extpar_edgar_-to_buffer     extpar_consistency_check
  INPUT_CDNC          settings for cdnc data                                           extpar_cdnc_-to_buffer      extpar_consistency_check
  INPUT_ERA           settings for ERA data                                            extpar_era_-to_buffer       extpar_consistency_check
  INPUT_CHECK         settings for the consistency check                               runscript                      extpar_consistency_check

## Grid definition {#namelist_input_for_extpar_grid_def}

The specification of the model type (COSMO or ICON) is done in the
namelist file INPUT_grid_org, the detailed target grid description for
the model domain has to be provided in the namelists files
INPUT_COSMO_GRID or INPUT_ICON_GRID.
[\[namelist_target\]]{#namelist_target label="namelist_target"}

### general {#namelist_input_for_extpar_grid_def_general}

### NAMELIST /grid_def/ (INPUT_grid_org) {#namelist-grid_def-input_grid_org .unnumbered}

The namelist /grid_def/ defines the target grid type and the filenames
with the namelists of the detailed target grid definition.

  **Parameter**           **Type**    **Default**   **Unit**   **Description**
  ----------------------- ----------- ------------- ---------- --------------------------------------------------------------------------
  igrid_type             integer                              target grid type, 1 for ICON, 2 for COSMO
  domain_def_namelist   character                            namelist file with domain definition
  domain_refinement      character                            namelist file with domain refinement definition (e.g. for the ICON grid)

### Icon {#namelist_input_for_extpar_grid_def_icon}

### NAMELIST /icon_grid_info/ (INPUT_ICON_GRID)  {#namelist-icon_grid_info-input_icon_grid .unnumbered}

The namelist /icon_grid_info/ specifies the filenames and the
directory of the Icon grid files with the coordinates of the Icon grid.

  **Parameter**          **Type**               **Default**   **Unit**   **Description**
  ---------------------- ---------------------- ------------- ---------- --------------------------------------------------------------------------
  icon_grid_dir        character                                       path to directory which contains the ICON grid file with the coordinates
  icon_grid_nc_file   character (max_dom)                            filename of the ICON grid file with the coordinates

### COSMO {#namelist_input_for_extpar_grid_def_cosmo}

### NAMELIST /lmgrid/ (INPUT_COSMO_GRID)  {#namelist-lmgrid-input_cosmo_grid .unnumbered}

The COSMO grid is defined by a rotated latlon-grid.

  **Parameter**   **Type**   **Default**   **Unit**   **Description**
  --------------- ---------- ------------- ---------- --------------------------------------------------------------------------------------------
  pollon          real       -170.         deg        longitude of the rotated north pole (in degrees, $E>0$)
  pollat          real       32.5          deg        latitude of the rotated north pole (in degrees, $N>0$)
  polgam          real       0\.           deg        longitude (in the rotated system) of the geographical north pole
  dlon            real       0.08          deg        grid point distance in zonal direction (in degrees)
  dlat            real       0.08          deg        grid point distance in meridional direction (in degrees)
  startlon_tot   real       -1.252        deg        transformed longitude of the lower left grid point of the total domain (in degrees, $E>0$)
  startlat_tot   real       -7.972        deg        transformed latitude of the lower left grid point of the total domain (in degrees, $N>0$)
  ie_tot         integer    51                       number of grid points in zonal direction
  je_tot         integer    51                       number of grid points in meridional direction
  ke_tot         integer    0                        number of grid points in vertical direction

## Orography {#namelist_input_for_extpar_orography}

### NAMELIST /oro_runcontrol/ (INPUT_ORO) {#namelist-oro_runcontrol-input_oro .unnumbered}

  **Parameter**    **Type**   **Default**   **Unit**   **Description**
  ---------------- ---------- ------------- ---------- ----------------------------------------------
  lcompute_sgsl   logical    .false.                  switch to activate subgrid-slope calculation

### NAMELIST /orography_raw_data/ (INPUT_ORO) {#namelist-orography_raw_data-input_oro .unnumbered}

  **Parameter**                **Type**    **Default**                        **Unit**   **Description**
  ---------------------------- ----------- ---------------------------------- ---------- ----------------------------------------------------------------------------------------
  itopo_type                  integer                                                   switch to choose an orography raw data set; 1 GLOBE, 2 ASTER, 3 MERIT/REMA
  lsso_param                  logical                                                   switch to choose if SSO parameters should be generated or not
  raw_data_orography_path   character                                                 path to orography raw data
  ntiles_column               integer     GLOBE: 4 ASTER, MERIT/REMA: x                 number of tile columns of desired region
  ntiles_row                  integer     GLOBE: 4 ASTER, MERIT/REMA: x                 number of tile rows of desired region
  topo_files                  character                                                 filenames of GLOBE (16 tiles) / ASTER (240 tiles)/ MERIT/REMA (72 tiles) raw data sets
  lsubtract_mean_slope       logical     .FALSE. for operational NWP-ICON              treatment of mean slope in computation of SSO parameters for ICON

### NAMELIST /orography_io_extpar/ (INPUT_ORO) {#namelist-orography_io_extpar-input_oro .unnumbered}

  **Parameter**             **Type**    **Default**   **Unit**   **Description**
  ------------------------- ----------- ------------- ---------- --------------------------------
  orography_buffer_file   character                            name for orography buffer file
  orography_output_file   character                            name for orography output file

### NAMELIST /sgsl_io_extpar/ (INPUT_ORO) {#namelist-sgsl_io_extpar-input_oro .unnumbered}

  **Parameter**   **Type**    **Default**   **Unit**   **Description**
  --------------- ----------- ------------- ---------- -----------------------------------------------------------------------------------------------------------
  lpreproc_oro   logical     .false.                  read S_ORO from existing NetCDF (.false.) or preprocess from raw topography datasets (.true.)
  sgsl_files     character                            filenames of raw data tiles to be used S_ORO_A10 to S_ORO_P10 (GLOBE) or S_ORO_T001 to S_ORO_T240

### NAMELIST /orography_smoothing/ (INPUT_OROSMOOTH) {#namelist-orography_smoothing-input_orosmooth .unnumbered}

  **Parameter**     **Type**   **Default**   **Unit**   **Description**
  ----------------- ---------- ------------- ---------- --------------------------------------------------------------------------
  lfilter_oro      logical    FALSE                    **Cosmo-only:** switch for orogaphy smoothing
  ilow_pass_oro   integer    0                        type of orogaphy smoothing and stencil width
  numfilt_oro      integer    1                        number of filter applications
  eps_filter       real       10                       smoothing parameter (\"strength\" of the filtering)
  ifill_valley     integer    1                        fill valleys before or after oro smoothing (1: before, 2: after)
  rfill_valley     real       0             m          mask for valley filling (threshold value)
  ilow_pass_xso   integer    1                        type of orogaphy eXtra SmOothing and stencil width (for steep orography)
  numfilt_xso      integer    1                        number of applications of the eXtra filter
  lxso_first       logical    FALSE                    eXtra SmOothing before or after orography smoothing (TRUE/FALSE)
  rxso_mask        real       0             m          mask for eXtra SmOothing (threshold value)

### NAMELIST /radtopo/ (INPUT_RADTOPO) {#namelist-radtopo-input_radtopo .unnumbered}

  **Parameter**    **Type**   **Default**   **Unit**   **Description**
  ---------------- ---------- ------------- ---------- ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  lradtopo         logical                             Switch for radiation corrected topography parameters. Not recommended to use if orographical smoothing is false and smoothing is performed in Int2lm later, because of resulting inconsistencies.
  nhori            integer    24                       Number of horizon angles
  radius           integer    40000         m          **Icon-only:** Radial distance considered for computation of horizon
  min_circ_cov   integer    1             \-         **Icon-only:** Number of gridcells to be skipped at circumference of circle. A value of 1 considers all points, whereas a value of 5 only consider every fifth point at the circumference. Note that the effect of this switch is dependent on the resolution of the grid as well on the radius choosen.
  max_missing     real       0.9           \-         **Icon-only:** Upper limit for fraction of missingness for the horizon parameter. Grid-cells with values above will be set to 0.
  itype_scaling   integer    2             \-         **Icon-only:** Power of the caling factor *SIN(horizon-angle)* applied to the geometric skyview factor to account for the anisotropic nature of longwave radiation.

### NAMELIST /scale_separated_raw_data/ (INPUT_SCALE_SEP) {#namelist-scale_separated_raw_data-input_scale_sep .unnumbered}

  **Parameter**                 **Type**    **Default**   **Unit**   **Description**
  ----------------------------- ----------- ------------- ---------- ---------------------------------------------------------------------------------------------
  lscale_separation            logical                              Switch for scale separation. It can only be used in combination with GLOBE as raw data set.
  raw_data_scale_sep_path   character                            path to 3 km filtered topography
  scale_sep_files             character                            filename of 3 km filtered topography

## Land Use data {#namelist_input_for_extpar_lu}

### NAMELIST /lu_raw_data/ (INPUT_LU) {#namelist-lu_raw_data-input_lu .unnumbered}

  **Parameter**             **Type**    **Default**   **Unit**   **Description**
  ------------------------- ----------- ------------- ---------- --------------------------------------------------------------------------------------------------------------------------------------
  raw_data_lu_path       character                            path to landuse data
  raw_data_lu_filename   character                            filename of landuse raw data
  i_landuse_data          integer                              switch to choose a land use raw data set; 1 Globcover2009, 2 GLC2000, 3 GLCC, 5 ESA CCI-LC, 6 Ecoclimap-SG
  l_use_corine            logical     .false.                  switch to use corine landuse dataset; only possible if i_landuse_data = 1
  ilookup_table_lu        integer                              switch to choose a lookup table
                                                                 GLC2000 and GLCC:
                                                                 1: operational settings of GME (Ritter, 2007)
                                                                 2: operational settings of COSMO (Heise, 2005)
                                                                 3: experimental setting, analog to look-up tables of ECOCLIMAP (Asensio 2010)
                                                                 GLOBCOVER 2009:
                                                                 1: operational settings (Asensio, 2011)
                                                                 2: experimental settings, analog to look-up tables of ECOCLIMAP (Asensio 2010)
                                                                 ESA CCI-LC:
                                                                 1: experimental settings (Helmert, 2019)
                                                                 Ecoclimap-SG:
                                                                 1: Globcover analogue with added LCZs from Oke
  ntiles_globcover         integer     6                        number of tiles for GLOBCOVER data
  ncolumn_tiles            integer                              number of columns in tile matrix
  l_terra_urb             logical     .false.                  switch to use TERRA-URB (see [3.2.2](#terra_urb){reference-type="ref" reference="terra_urb"}); only possible if i_landuse_data = 6

### NAMELIST /glcc_raw_data/ (INPUT_LU) {#namelist-glcc_raw_data-input_lu .unnumbered}

  **Parameter**               **Type**    **Default**   **Unit**   **Description**
  --------------------------- ----------- ------------- ---------- -------------------------------------------------------------------------------
  raw_data_glcc_path       character                            path to glcc data
  raw_data_glcc_filename   character                            filename of glcc raw data
  ilookup_table_glcc        integer                              switch to choose a lookup table
                                                                   1: operational settings of GME (Ritter, 2007)
                                                                   2: operational settings of COSMO (Heise, 2005)
                                                                   3: experimental setting, analog to look-up tables of ECOCLIMAP (Asensio 2010)

### NAMELIST /glcc_io_extpar/ (INPUT_LU) {#namelist-glcc_io_extpar-input_lu .unnumbered}

  **Parameter**        **Type**    **Default**   **Unit**   **Description**
  -------------------- ----------- ------------- ---------- ---------------------------
  glcc_buffer_file   character                            name for glcc buffer file

## Aerosol optical depth {#namelist_input_for_extpar_aot}

### NAMELIST /aerosol_raw_data/ (INPUT_AOT) {#namelist-aerosol_raw_data-input_aot .unnumbered}

  **Parameter**              **Type**    **Default**   **Unit**   **Description**
  -------------------------- ----------- ------------- ---------- -------------------------------------------------------------------------------
  raw_data_aot_path       character                            path to aerosol raw data
  raw_data_aot_filename   character                            filename of aerosol raw data
  iaot_type                 integer     1                        index to specify AOD raw data set: 1:Tegen 2:AeroCom 3:MACC-II 4:MACv2 5:CAMS

### NAMELIST /aerosol_io_extpar/ (INPUT_AOT) {#namelist-aerosol_io_extpar-input_aot .unnumbered}

  **Parameter**       **Type**    **Default**   **Unit**   **Descriptio**n
  ------------------- ----------- ------------- ---------- ------------------------------
  aot_buffer_file   character                            name for aerosol buffer file

## Climatological 2m temperature {#namelist_input_for_extpar_cru}

### DICT input_tclim (namelist.py) {#dict-input_tclim-namelist.py .unnumbered}

  **Parameter**                **Type**    **Default**   **Unit**   **Description**
  ---------------------------- ----------- ------------- ---------- --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  raw_data_t_clim_path     character                            path to T2m climatology data
  raw_data_t_clim_coarse   character                            filename of coarse T2m climatology data
  raw_data_t_clim_fine     character                            filename of fine T2m climatology data
  it_cl_type                 integer                              switch to choose between the new and fine (1) and the old and coarse over sea and the fine over land (2) raw data set. Note that the fine data set (1) is topographically corrected.
  t_clim_buffer_file        character                            name for t_clim buffer file

## NDVI data {#namelist_input_for_extpar_ndvi}

### DICT input_ndvi (namelist.py) {#dict-input_ndvi-namelist.py .unnumbered}

  **Parameter**               **Type**    **Default**   **Unit**   **Description**
  --------------------------- ----------- ------------- ---------- ---------------------------
  raw_data_ndvi_path       character                            path to ndvi raw data
  raw_data_ndvi_filename   character                            filename of ndvi raw data
  ndvi_buffer_file          character                            name for ndvi buffer file

## EDGAR data {#namelist_input_for_extpar_edgar}

### DICT input_edgar (namelist.py) {#dict-input_edgar-namelist.py .unnumbered}

  **Parameter**                     **Type**    **Default**   **Unit**   **Description**
  --------------------------------- ----------- ------------- ---------- --------------------------------------------
  raw_data_edgar_path            character                            path to EDGAR raw data
  raw_data_edgar_filename_bc    character                            filename of EDGAR black carbon raw data
  raw_data_edgar_filename_oc    character                            filename of EDGAR organic carbon raw data
  raw_data_edgar_filename_so2   character                            filename of EDGAR sulfur dioxide raw data
  raw_data_edgar_filename_nh3   character                            filename of EDGAR ammonia raw data
  raw_data_edgar_filename_nox   character                            filename of EDGAR nitrogen oxides raw data

## CDNC data {#namelist_input_for_cdnc}

### DICT input_cdnc (namelist.py) {#dict-input_cdnc-namelist.py .unnumbered}

  **Parameter**               **Type**    **Default**   **Unit**   **Description**
  --------------------------- ----------- ------------- ---------- ---------------------------
  raw_data_cdnc_path       character                            path to cdnc raw data
  raw_data_cdnc_filename   character                            filename of cdnc raw data

## hwsdART data {#namelist_input_for_hwsdART}

### NAMELIST /hwsdART_nml/ (INPUT_hwsdART) {#namelist-hwsdart_nml-input_hwsdart .unnumbered}

  **Parameter**                  **Type**    **Default**   **Unit**   **Description**
  ------------------------------ ----------- ------------- ---------- ------------------------------
  raw_data_hwsdART_path       character                            path to hwsdART raw data
  raw_data_hwsdART_filename   character                            filename of hwsdART raw data
  hwsdART_output_file          character                            name for hwsdART output file

## Soil data {#namelist_input_for_extpar_soil}

### NAMELIST /soil_raw_data/ (INPUT_SOIL) {#namelist-soil_raw_data-input_soil .unnumbered}

  **Parameter**                     **Type**    **Default**   **Unit**   **Description**
  --------------------------------- ----------- ------------- ---------- ------------------------------------------------------------------------------------------------------------
  isoil_data                       integer                              switch to choose between the raw soil data, 1: FAO, 2: HWSD, 3:HWSD with terra mapping
  ldeep_soil                       logical                              switch for the deep soil, can only be set to .TRUE. if the top and sub soil are deduced from the HWSD data
  raw_data_soil_path             character                            path to soil raw data
  raw_data_soil_filename         character                            filename of soil raw data
  raw_data_deep_soil_filename   character                            filename of deep soil raw data

### NAMELIST /soil_io_extpar/ (INPUT_SOIL) {#namelist-soil_io_extpar-input_soil .unnumbered}

  **Parameter**                    **Type**    **Default**   **Unit**   **Description**
  -------------------------------- ----------- ------------- ---------- ---------------------------------------------------
  soil_buffer_file               character                            name for soil buffer file
  soil_buffer_file_consistent   character                            name for soil buffer file after consistency check
  soil_output_file_consistent   character                            name for soil output file after consistency check

### NAMELIST /HWSD_index_files/ (INPUT_SOIL) {#namelist-hwsd_index_files-input_soil .unnumbered}

  **Parameter**              **Type**    **Default**   **Unit**   **Description**
  -------------------------- ----------- ------------- ---------- -------------------------------------------------------------------------------------------------------------------------------------
  path_HWSD_index_files   character                            path to HWSD lookup tables
  lookup_table_HWSD        character                            lookup table to convert soiltype index from global to TERRA soiltype
  HWSD_data                 character                            lookup table to convert the global soiltype in fractions of sand, silt, clay and organic carbon and in bulk density for the topsoil
  HWSD_data_deep           character                            lookup table to convert the global soiltype in fractions of sand, silt, clay and organic carbon and in bulk density for the subsoil
  HWSD_data_extpar         character                            parameter used only for development purposes

## Freshwater Lake data {#namelist_input_for_extpar_flake}

### NAMELIST /flake_raw_data/ (INPUT_FLAKE) {#namelist-flake_raw_data-input_flake .unnumbered}

  **Parameter**                **Type**    **Default**   **Unit**   **Description**
  ---------------------------- ----------- ------------- ---------- ----------------------------
  raw_data_flake_path       character                            path to flake raw data
  raw_data_flake_filename   character                            filename of flake raw data

### NAMELIST /flake_io_extpar/ (INPUT_FLAKE) {#namelist-flake_io_extpar-input_flake .unnumbered}

  **Parameter**         **Type**    **Default**   **Unit**   **Description**
  --------------------- ----------- ------------- ---------- ----------------------------
  flake_buffer_file   character                            name for flake buffer file

## Albedo data {#namelist_input_for_extpar_albedo}

### DICT input_alb (namelist.py) {#dict-input_alb-namelist.py .unnumbered}

  **Parameter**                **Type**    **Default**   **Unit**   **Description**
  ---------------------------- ----------- ------------- ---------- ---------------------------------------------------------------------------------------------
  raw_data_alb_path         character                            path to raw albedo data
  raw_data_alb_filename     character                            filename of the raw albedo data
  raw_data_alnid_filename   character                            filename of the raw NIR-albedo data
  raw_data_aluvd_filename   character                            filename of the raw UV-albedo data
  ialb_type                   integer                              switch to indicate albedo type 1:total albedo 2:soil albedo 3:as 1 without NI and UV fields
  alb_buffer_file            character                            name for the albedo buffer file

## ISA data {#namelist_input_for_extpar_isa}

### DICT input_isa (namelist.py) {#dict-input_isa-namelist.py .unnumbered}

  **Parameter**              **Type**    **Default**   **Unit**   **Description**
  -------------------------- ----------- ------------- ---------- ------------------------------
  raw_data_isa_path       character                            path to ISA raw data
  raw_data_isa_filename   character                            filename of ISA raw data
  isa_type                  integer                              type of used ISA data source
  isa_buffer_file          character                            name for ISA buffer file

## AHF data {#namelist_input_for_extpar_ahf}

### DICT input_ahf (namelist.py) {#dict-input_ahf-namelist.py .unnumbered}

  **Parameter**              **Type**    **Default**   **Unit**   **Description**
  -------------------------- ----------- ------------- ---------- ------------------------------
  raw_data_ahf_path       character                            path to AHF raw data
  raw_data_ahf_filename   character                            filename of AHF raw data
  iahf_type                 integer                              type of used AHF data source
  ahf_buffer_file          character                            name for AHF buffer file

## Emissivity parameter {#namelist_input_for_extpar_emissivity}

### DICT input_emiss (namelist.py) {#dict-input_emiss-namelist.py .unnumbered}

  **Parameter**                **Type**    **Default**   **Unit**   **Description**
  ---------------------------- ----------- ------------- ---------- ---------------------------------------------------------------------------
  iemiss_type                 integer                              switch to choose between full-range (1) and long-wave (2) emissivity data
  raw_data_emiss_path       character                            path to emissivity parameter raw data
  raw_data_emiss_filename   character                            filenames of emissivity raw data
  emiss_buffer_file          character                            name for emissivity parameter buffer file

## ERA parameter {#namelist_input_for_extpar_era}

### DICT input_era (namelist.py) {#dict-input_era-namelist.py .unnumbered}

  **Parameter**          **Type**    **Default**   **Unit**   **Description**
  ---------------------- ----------- ------------- ---------- -------------------------------------------------
  iera_type             integer                              type of ERA climatology: ERA5 (1) and ERA-I (2)
  raw_data_era_path   character                            path to ERA raw data
  raw_data_era_ORO    character                            filenames of ERA ORO raw data
  raw_data_era_SD     character                            filenames of ERA SD raw data
  raw_data_era_T2M    character                            filenames of ERA T2M raw data
  raw_data_era_SST    character                            filenames of ERA SST raw data
  era_buffer_file      character                            name for era parameter buffer file

## Consistency check {#namelist_input_for_extpar_consistency_check}

### NAMELIST /extpar_consistency_check_io/ (INPUT_CHECK) {#namelist-extpar_consistency_check_io-input_check .unnumbered}

  **Parameter**              **Type**    **Default**   **Unit**   **Description**
  -------------------------- ----------- ------------- ---------- --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  l_use_array_cache       flag        F                        flag indicating whether mmap-caching is used (reduced memory consumption, but longer runtime)
  netcdf_output_filename   character                            filename for netcdf output filename
  i_lsm_data               integer                              integer switch to choose if an external land sea mask is desired or not. 0: no external land sea mask, 1: use external land sea mask
  land_sea_mask_file      character                            name of the file which can be used as the external land sea mask.
  number_special_points    integer                              number of points that should be treated specially (maximal value: 3, if no special treatment is desired choose 0).
  lwrite_netcdf             logical     T                        flag indicating whether NetCDF output for COSMO grid is desired.
  tile_mode                 integer     0                        if activated tile_mode=1 process output for ICON tile structure
  lflake_correction         logical     T                        if activated then fr_lake values of grid points next to ocean are set to ocean values the lake depth value is set to undefined. This is the default behavior in EXTPAR version 4.0, but not in DWD EXTPAR version 2.10.

### NAMELIST /special_points/ (INPUT_SP_1) {#namelist-special_points-input_sp_1 .unnumbered}

Modifications for Falkenberg.\

  **Parameter**   **Type**   **Default**   **Unit**    **Description**
  --------------- ---------- ------------- ----------- ------------------------------------------------------------------
  lon_geo_sp    real       14.115        deg east    longitude coordinate of the special point
  lat_geo_sp    real       52.156        deg north   latitude coordinate of the special point
  soiltype_sp    real       3.0           \-          soiltype of the special point
  z0_sp          real       0.03          m           roughness lenght of the special point
  rootdp_sp      real       0.6           m           rooting depth of the special point
  plcovmn_sp     real       0.55          1           plant cover minimum of the special point
  plcovmx_sp     real       0.8           1           plant cover maximum of the special point
  laimn_sp       real       0.5           1           leaf area index minimum of the special point
  laimx_sp       real       2.5           1           leaf area index maximum of the special point
  for_d_sp      real                     1           ground fraction covered by deciduous forest of the special point
  for_e_sp      real                     1           ground fraction covered by evergreen forest of the special point
  fr_land_sp    real                     1           fraction land cover of the special point

### NAMELIST /special_points/ (INPUT_SP_2) {#namelist-special_points-input_sp_2 .unnumbered}

Modifications for Waldstation.\

  **Parameter**   **Type**   **Default**   **Unit**    **Description**
  --------------- ---------- ------------- ----------- ------------------------------------------------------------------
  lon_geo_sp    real       13.954        deg east    longitude coordinate of the special point
  lat_geo_sp    real       52.186        deg north   latitude coordinate of the special point
  soiltype_sp    real       3\.           \-          soiltype of the special point
  z0_sp          real       0.78          m           roughness lenght of the special point
  rootdp_sp      real       0.6           m           rooting depth of the special point
  plcovmn_sp     real       0.79          1           plant cover minimum of the special point
  plcovmx_sp     real       0.81          1           plant cover maximum of the special point
  laimn_sp       real       3.0           1           leaf area index minimum of the special point
  laimx_sp       real       4.0           1           leaf area index maximum of the special point
  for_d_sp      real                     1           ground fraction covered by deciduous forest of the special point
  for_e_sp      real                     1           ground fraction covered by evergreen forest of the special point
  fr_land_sp    real                     1           fraction land cover of the special point

### NAMELIST /special_points/ (INPUT_SP_3) {#namelist-special_points-input_sp_3 .unnumbered}

Modifications for Lindenberg.\

  **Parameter**   **Type**   **Default**   **Unit**    **Description**
  --------------- ---------- ------------- ----------- ------------------------------------------------------------------
  lon_geo_sp    real       14.119        deg east    longitude coordinate of the special point
  lat_geo_sp    real       52.205        deg north   latitude coordinate of the special point
  soiltype_sp    real       5\.           \-          soiltype of the special point
  z0_sp          real       0.07          m           roughness lenght of the special point
  rootdp_sp      real       0.6           m           rooting depth of the special point
  plcovmn_sp     real       0.5           1           plant cover minimum of the special point
  plcovmx_sp     real       0.9           1           plant cover maximum of the special point
  laimn_sp       real       0.7           1           leaf area index minimum of the special point
  laimx_sp       real       3.3           1           leaf area index maximum of the special point
  for_d_sp      real                     1           ground fraction covered by deciduous forest of the special point
  for_e_sp      real                     1           ground fraction covered by evergreen forest of the special point
  fr_land_sp    real                     1           fraction land cover of the special point

[^1]: Tegen, I., P. Hollrigl, M. Chin, I. Fung, D. Jacob, and J. Penner
    1997. Contribution of different aerosol species to the global
    aerosol extinction optical thickness: Estimates from model results.
    J. Geophys. Res., 102, 23895-23915.
    <http://pubs.giss.nasa.gov/abstracts/1997/Tegen_etal.html>

[^2]: Kinne, S., M. Schulz, C. Textor, S. Guibert, Y. Balkanski, S.E.
    Bauer, T. Berntsen, T.F. Berglen, O. Boucher, M. Chin, W. Collins,
    F. Dentener, T. Diehl, R. Easter, J. Feichter, D. Fillmore, S. Ghan,
    P. Ginoux, S. Gong, A. Grini, J. Hendricks, M. Herzog, L. Horowitz,
    I. Isaksen, T. Iversen, A. Kirkevåg, S. Kloster, D. Koch, J.E.
    Kristjansson, M. Krol, A. Lauer, J.F. Lamarque, G. Lesins, X. Liu,
    U. Lohmann, V. Montanaro, G. Myhre, J. Penner, G. Pitari, S. Reddy,
    Ø. Seland, P. Stier, T. Takemura, and X. Tie: An AeroCom initial
    assessment optical properties in aerosol component modules of global
    models, Atmos. Chem. Phys., 6, 1815-1834, 2006.
    <http://aerocom.met.no/aerocomhome.html>

[^3]: Morcrette, J.-J., O. Boucher, L. Jones, D. Salmond, P. Bechtold,
    A. Beljaars, A. Benedetti, A. Bonet, J. W. Kaiser, M. Razinger, M.
    Schulz, S. Serrar, A. J. Simmons, M. Sofiev, M. Suttie, A. M.
    Tompkins, and A. Untch, 2009: Aerosol analysis and forecast in the
    ECMWF Integrated Forecast System. Part I: Forward modelling, J.
    Geophys. Res., 114D, D06206,doi:10.1029/2008JD011235
    <http://www.gmes-atmosphere.eu/>

[^4]: Kinne, S., D. O'Donnel, P. Stier, S. Kloster, K. Zhang, H.
    Schmidt, S. Rast, M. Giorgetta, T. F. Eck, and B. Stevens (2013),
    MAC-v1: A new global aerosol climatology for climate studies, J.
    Adv. Model. Earth Syst., 5, 704--740, doi:10.1002/jame.20035.

[^5]: Bozzo, A., Benedetti, A., Flemming, J., Kipling, Z., Rémy, S.
    (2020). An aerosol climatology for global models based on the
    tropospheric aerosol scheme in the Integrated Forecasting System of
    ECMWF. Geoscientific Model Development, 13(3), 1007-1034.

[^6]: Release 5.0 of COSMO and 2.0 of INT2LM do not support HWSD data,
    as the representation of the soil associated with this new data set
    has changed and is based on the use of pedotransfer functions and on
    fraction of soil components (e.g. clay, silt, \...)

[^7]: Soiltypes written in bold indicate a special soiltype.

[^8]: J.M. Hollisa, J. Hannamb and P.H. Bellamyb, February 2012,
    Empirically-derived pedotransfer functions for predicting bulk
    density in European soils, European Journal of Soil Science, 63,
    96109 doi: 10.1111/j.1365-2389.2011.01412.x

[^9]: J.H.M. Woesten, A. Lilly, A. Nemes and C. Le Bas, Development and
    use of a database of hydraulic properties of European soils,
    Geoderma, Volume 90, Issues 34, July 1999, Pages 169-185, doi:
    10.1016/S0016-7061(98)00132-3.

[^10]: derived from the land use data

[^11]: derived from the land use data

[^12]: derived from the lake data set

[^13]: https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/rawdata/mksrf_soilcol.081008.nc\
    Lawrence, P. J. and T. N. Chase (2007). \"Representing a new MODIS
    consistent land surface in the Community Land Model (CLM 3.0).\"
    Journal of Geophysical Research-Biogeosciences 112(G1).\
    Table 3.3 in: Oleson, K.W., D.M. Lawrence, G.B. Bonan, M.G. Flanner,
    E. Kluzek, P.J. Lawrence, S. Levis, S.C. Swenson, P.E. Thornton, A.
    Dai, M. Decker, R. Dickinson, J. Feddema, C.L. Heald, F. Hoffman,
    J.-F. Lamarque, N. Mahowald, G.-Y. Niu, T. Qian, J. Randerson, S.
    Running, K. Sakaguchi, A. Slater, R. Stockli, A. Wang, Z.-L. Yang,
    Xi. Zeng, and Xu. Zeng, 2010: Technical Description of version 4.0
    of the Community Land Model (CLM). NCAR Technical Note
    NCAR/TN-478+STR, National Center for Atmospheric Research, Boulder,
    CO, 257 pp.
