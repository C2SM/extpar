# Namelist Input {#namelist_input_for_extpar}

EXTPAR uses 3 types of namelists in order to determine in which way data is processed.

- Fortran namelists (INPUT_)
- Python dictionaries (input_in namelist.py)
- Fortran namelists written by Python scripts

Whereas for the Fortran namelists and the Python dictionaries the user can specify parameters and filenames, the Fortran namelists generated during runtime by the Python scripts do not allow any user interaction.

## Namelist files {#namelist_input_for_extpar_namelist_files}

| **Namelist file**    | **Purpose**                                                      | **Made by script**             | **Used by program**             |
| -------------------- | ---------------------------------------------------------------- | ------------------------------ | -----------------------------   |
| INPUT_grid_org       | define target grid type                                          | runscript                      | extpar_consistency_check        |
|                      |                                                                  |                                | extpar_aot_to_buffer            |
|                      |                                                                  |                                | extpar_landuse_to_buffer        |
|                      |                                                                  |                                | extpar_topo_to_buffer           |
|                      |                                                                  |                                | extpar_cru_to_buffer            |
|                      |                                                                  |                                | extpar_ndvi_to_buffer           |
|                      |                                                                  |                                | extpar_soil_to_buffer           |
|                      |                                                                  |                                | extpar_flake_to_buffer          |
|                      |                                                                  |                                | extpar_isa_to_buffer            |
|                      |                                                                  |                                | extpar_ahf_to_buffer            |
|                      |                                                                  |                                | extpar_emiss_to_buffer          |
|                      |                                                                  |                                | extpar_hwsdART_to_buffer        |
| INPUT_COSMO_GRID     | define target domain for COSMO grid                              | runscript                      | extpar_consistency_check        |
|                      |                                                                  |                                | extpar_aot_to_buffer            |
|                      |                                                                  |                                | extpar_landuse_to_buffer        |
|                      |                                                                  |                                | extpar_topo_to_buffer           |
|                      |                                                                  |                                | extpar_cru_to_buffer            |
|                      |                                                                  |                                | extpar_ndvi_to_buffer           |
|                      |                                                                  |                                | extpar_soil_to_buffer           |
|                      |                                                                  |                                | extpar_flake_to_buffer          |
|                      |                                                                  |                                | extpar_isa_to_buffer            |
|                      |                                                                  |                                | extpar_ahf_to_buffer            |
|                      |                                                                  |                                | extpar_emiss_to_buffer          |
|                      |                                                                  |                                | extpar_hwsdART_to_buffer        |
| INPUT_ICON_GRID      | define target domain for ICON grid                               | runscript                      | extpar_consistency_check        |
|                      |                                                                  |                                | extpar_aot_to_buffer            |
|                      |                                                                  |                                | extpar_landuse_to_buffer        |
|                      |                                                                  |                                | extpar_topo_to_buffer           |
|                      |                                                                  |                                | extpar_cru_to_buffer            |
|                      |                                                                  |                                | extpar_ndvi_to_buffer           |
|                      |                                                                  |                                | extpar_soil_to_buffer           |
|                      |                                                                  |                                | extpar_flake_to_buffer          |
|                      |                                                                  |                                | extpar_isa_to_buffer            |
|                      |                                                                  |                                | extpar_ahf_to_buffer            |
|                      |                                                                  |                                | extpar_emiss_to_buffer          |
| INPUT_ORO            | settings for orography data                                      | runscript                      | extpar_topo_to_buffer           |
| INPUT_OROSMOOTH      | settings for orography smoothing                                 | runscript                      | extpar_topo_to_buffer           |
| INPUT_RADTOPO        | settings for generating topographical shading fields             | runscript                      | extpar_topo_to_buffer           |
| INPUT_SCALE_SEP      | settings to control scale separation for SSO an Z0 calculation   | runscript                      | extpar_topo_to_buffer           |
| INPUT_LU             | settings for landuse data                                        | runscript                      | extpar_landuse_to_buffer        |
| INPUT_AOT            | settings for aerosol data                                        | runscript                      | extpar_aot_to_buffer            |
| INPUT_TCLIM          | settings for temperature data                                    | extpar_cru_-to_buffer          | extpar_consistency_check        |
| INPUT_NDVI           | settings for NDVI data                                           | extpar_ndvi_-to_buffer         | extpar_consistency_check        |
| INPUT_SOIL           | settings for soil data                                           | runscript                      | extpar_soil_to_buffer           |
| INPUT_FLAKE          | settings for lake data                                           | runscript                      | extpar_flake_to_buffer          |
| INPUT_ALB            | settings for albedo data                                         | extpar_albedo_-to_buffer       | extpar_consistency_check        |
| INPUT_ISA            | settings for fraction of impervious surface area data            | extpar_isa_to_buffer           | extpar_consistency_check        |
| INPUT_AHF            | settings for anthropogenic heat flux data                        | extpar_ahf_to_buffer           | extpar_consistency_check        |
| INPUT_EMISS          | settings for emissivity data                                     | extpar_emiss_-to_buffer        | extpar_consistency_check        |
| INPUT_hwsdART        | settings for HWSD USDA data                                      | extpar_hwsdART_-to_buffer      |                                 |
| INPUT_edgar          | settings for EDGAR data                                          | extpar_edgar_-to_buffer        | extpar_consistency_check        |
| INPUT_CDNC           | settings for cdnc data                                           | extpar_cdnc_-to_buffer         | extpar_consistency_check        |
| INPUT_ERA            | settings for ERA data                                            | extpar_era_-to_buffer          | extpar_consistency_check        |
| INPUT_CHECK          | settings for the consistency check                               | runscript                      | extpar_consistency_check        |

## Grid definition {#namelist_input_for_extpar_grid_def}

The specification of the model type (COSMO or ICON) is done in the namelist file INPUT_grid_org, the detailed target grid description for the model domain has to be provided in the namelists files INPUT_COSMO_GRID or INPUT_ICON_GRID.
[\[namelist_target\]]{#namelist_target label="namelist_target"}

### general {#namelist_input_for_extpar_grid_def_general}

### NAMELIST /grid_def/ (INPUT_grid_org) {#namelist-grid_def-input_grid_org .unnumbered}

The namelist /grid_def/ defines the target grid type and the filenames with the namelists of the detailed target grid definition.

| **Parameter**           | **Type**    | **Default**   | **Unit**   | **Description**                                                                 |
| ----------------------- | ----------- | ------------- | ---------- | ------------------------------------------------------------------------------- |
| igrid_type              | integer     |               |            | target grid type, 1 for ICON, 2 for COSMO                                       |
| domain_def_namelist     | character   |               |            | namelist file with domain definition                                            |
| domain_refinement       | character   |               |            | namelist file with domain refinement definition (e.g. for the ICON grid)        |

### Icon {#namelist_input_for_extpar_grid_def_icon}

### NAMELIST /icon_grid_info/ (INPUT_ICON_GRID)  {#namelist-icon_grid_info-input_icon_grid .unnumbered}

The namelist /icon_grid_info/ specifies the filenames and the directory of the Icon grid files with the coordinates of the Icon grid.

| **Parameter**          | **Type**               | **Default**   | **Unit**   | **Description**                                                                 |
| ---------------------- | ---------------------- | ------------- | ---------- | ------------------------------------------------------------------------------- |
| icon_grid_dir          | character              |               |            | path to directory which contains the ICON grid file with the coordinates        |
| icon_grid_nc_file      | character (max_dom)    |               |            | filename of the ICON grid file with the coordinates                             |

### COSMO {#namelist_input_for_extpar_grid_def_cosmo}

### NAMELIST /lmgrid/ (INPUT_COSMO_GRID)  {#namelist-lmgrid-input_cosmo_grid .unnumbered}

The COSMO grid is defined by a rotated latlon-grid.

| **Parameter**   | **Type**   | **Default**   | **Unit**   | **Description**                                                                 |
| --------------- | ---------- | ------------- | ---------- | --------------------------------------------------------------------------------|
| pollon          | real       | -170.         | deg        | longitude of the rotated north pole (in degrees, $E>0$)                         |
| pollat          | real       | 32.5          | deg        | latitude of the rotated north pole (in degrees, $N>0$)                          |
| polgam          | real       | 0\.           | deg        | longitude (in the rotated system) of the geographical north pole                |
| dlon            | real       | 0.08          | deg        | grid point distance in zonal direction (in degrees)                             |
| dlat            | real       | 0.08          | deg        | grid point distance in meridional direction (in degrees)                        |
| startlon_tot    | real       | -1.252        | deg        | transformed longitude of the lower left grid point of the total domain (in degrees, $E>0$) |
| startlat_tot    | real       | -7.972        | deg        | transformed latitude of the lower left grid point of the total domain (in degrees, $N>0$) |
| ie_tot          | integer    | 51            |            | number of grid points in zonal direction                                        |
| je_tot          | integer    | 51            |            | number of grid points in meridional direction                                   |
| ke_tot          | integer    | 0             |            | number of grid points in vertical direction                                     |

## Orography {#namelist_input_for_extpar_orography}

### NAMELIST /oro_runcontrol/ (INPUT_ORO) {#namelist-oro_runcontrol-input_oro .unnumbered}

| **Parameter**    | **Type**   | **Default**   | **Unit**   | **Description**                                                                 |
| ---------------- | ---------- | ------------- | ---------- | ------------------------------------------------------------------------------- |
| lcompute_sgsl    | logical    | .false.       |            | switch to activate subgrid-slope calculation                                    |

### NAMELIST /orography_raw_data/ (INPUT_ORO) {#namelist-orography_raw_data-input_oro .unnumbered}

| **Parameter**                | **Type**    | **Default**                        | **Unit**   | **Description**                                                                 |
| ---------------------------- | ----------- | ---------------------------------- | ---------- | --------------------------------------------------------------------------------|
| itopo_type                   | integer     |                                    |            | switch to choose an orography raw data set; 1 GLOBE, 2 ASTER, 3 MERIT/REMA      |
| lsso_param                   | logical     |                                    |            | switch to choose if SSO parameters should be generated or not                   |
| raw_data_orography_path      | character   |                                    |            | path to orography raw data                                                      |
| ntiles_column                | integer     | GLOBE: 4 ASTER, MERIT/REMA: x      |            | number of tile columns of desired region                                        |
| ntiles_row                   | integer     | GLOBE: 4 ASTER, MERIT/REMA: x      |            | number of tile rows of desired region                                           |
| topo_files                   | character   |                                    |            | filenames of GLOBE (16 tiles) / ASTER (240 tiles)/ MERIT/REMA (72 tiles) raw data sets |
| lsubtract_mean_slope         | logical     | .FALSE. for operational NWP-ICON   |            | treatment of mean slope in computation of SSO parameters for ICON               |

### NAMELIST /orography_io_extpar/ (INPUT_ORO) {#namelist-orography_io_extpar-input_oro .unnumbered}

| **Parameter**             | **Type**    | **Default**   | **Unit**   | **Description**                                                                 |
| ------------------------- | ----------- | ------------- | ---------- | --------------------------------------------------------------------------------|
| orography_buffer_file     | character   |               |            | name for orography buffer file                                                  |
| orography_output_file     | character   |               |            | name for orography output file                                                  |

### NAMELIST /sgsl_io_extpar/ (INPUT_ORO) {#namelist-sgsl_io_extpar-input_oro .unnumbered}

| **Parameter**   | **Type**    | **Default**   | **Unit**   | **Description**                                                                 |
| --------------- | ----------- | ------------- | ---------- | --------------------------------------------------------------------------------|
| lpreproc_oro    | logical     | .false.       |            | read S_ORO from existing NetCDF (.false.) or preprocess from raw topography datasets (.true.) |
| sgsl_files      | character   |               |            | filenames of raw data tiles to be used S_ORO_A10 to S_ORO_P10 (GLOBE) or S_ORO_T001 to S_ORO_T240 |

### NAMELIST /orography_smoothing/ (INPUT_OROSMOOTH) {#namelist-orography_smoothing-input_orosmooth .unnumbered}

| **Parameter**     | **Type**   | **Default**   | **Unit**   | **Description**                                                                 |
| ----------------- | ---------- | ------------- | ---------- | --------------------------------------------------------------------------------|
| lfilter_oro       | logical    | FALSE         |            | **Cosmo-only:** switch for orogaphy smoothing                                   |
| ilow_pass_oro     | integer    | 0             |            | type of orogaphy smoothing and stencil width                                    |
| numfilt_oro       | integer    | 1             |            | number of filter applications                                                   |
| eps_filter        | real       | 10            |            | smoothing parameter ("strength" of the filtering)                               |
| ifill_valley      | integer    | 1             |            | fill valleys before or after oro smoothing (1: before, 2: after)                |
| rfill_valley      | real       | 0             | m          | mask for valley filling (threshold value)                                       |
| ilow_pass_xso     | integer    | 1             |            | type of orogaphy eXtra SmOothing and stencil width (for steep orography)        |
| numfilt_xso       | integer    | 1             |            | number of applications of the eXtra filter                                      |
| lxso_first        | logical    | FALSE         |            | eXtra SmOothing before or after orography smoothing (TRUE/FALSE)                |
| rxso_mask         | real       | 0             | m          | mask for eXtra SmOothing (threshold value)                                      |

### NAMELIST /radtopo/ (INPUT_RADTOPO) {#namelist-radtopo-input_radtopo .unnumbered}

| **Parameter**    | **Type**   | **Default**   | **Unit**   | **Description**                                                                 |
| ---------------- | ---------- | ------------- | ---------- | ------------------------------------------------------------------------------- |
| lradtopo         | logical    |               |            | Switch for radiation corrected topography parameters. Not recommended to use if orographical smoothing is false and smoothing is performed in Int2lm later, because of resulting inconsistencies. |
| nhori            | integer    | 24            |            | Number of horizon angles                                                        |
| radius           | integer    | 40000         | m          | **Icon-only:** Radial distance considered for computation of horizon            |
| min_circ_cov     | integer    | 1             | -          | **Icon-only:** Number of gridcells to be skipped at circumference of circle. A value of 1 considers all points, whereas a value of 5 only consider every fifth point at the circumference. Note that the effect of this switch is dependent on the resolution of the grid as well on the radius choosen. |
| max_missing      | real       | 0.9           | -          | **Icon-only:** Upper limit for fraction of missingness for the horizon parameter. Grid-cells with values above will be set to 0. |
| itype_scaling    | integer    | 2             | -          | **Icon-only:** Power of the caling factor *SIN(horizon-angle)* applied to the geometric skyview factor to account for the anisotropic nature of longwave radiation. |

### NAMELIST /scale_separated_raw_data/ (INPUT_SCALE_SEP) {#namelist-scale_separated_raw_data-input_scale_sep .unnumbered}

| **Parameter**                 | **Type**    | **Default**   | **Unit**   | **Description**                                                                 |
| ----------------------------- | ----------- | ------------- | ---------- | --------------------------------------------------------------------------------|
| lscale_separation             | logical     |               |            | Switch for scale separation. It can only be used in combination with GLOBE as raw data set. |
| raw_data_scale_sep_path       | character   |               |            | path to 3 km filtered topography                                                |
| scale_sep_files               | character   |               |            | filename of 3 km filtered topography                                            |

## Land Use data {#namelist_input_for_extpar_lu}

### NAMELIST /lu_raw_data/ (INPUT_LU) {#namelist-lu_raw_data-input_lu .unnumbered}

| **Parameter**             | **Type**    | **Default**   | **Unit**   | **Description**                                                                 |
| ------------------------- | ----------- | ------------- | ---------- | --------------------------------------------------------------------------------|
| raw_data_lu_path          | character   |               |            | path to landuse data                                                            |
| raw_data_lu_filename      | character   |               |            | filename of landuse raw data                                                    |
| i_landuse_data            | integer     |               |            | switch to choose a land use raw data set; 1 Globcover2009, 2 GLC2000, 3 GLCC, 5 ESA CCI-LC, 6 Ecoclimap-SG |
| l_use_corine              | logical     | .false.       |            | switch to use corine landuse dataset; only possible if i_landuse_data = 1       |
| ilookup_table_lu          | integer     |               |            | switch to choose a lookup table                                                 |
|                           |             |               |            | GLC2000 and GLCC:                                                               |
|                           |             |               |            | 1: operational settings of GME (Ritter, 2007)                                   |
|                           |             |               |            | 2: operational settings of COSMO (Heise, 2005)                                  |
|                           |             |               |            | 3: experimental setting, analog to look-up tables of ECOCLIMAP (Asensio 2010)   |
|                           |             |               |            | GLOBCOVER 2009:                                                                 |
|                           |             |              # Namelist Input for EXTPAR {#namelist_input_for_extpar}

EXTPAR uses 3 types of namelists in order to determine in which way data is processed.

- Fortran namelists (INPUT_)
- Python dictionaries (input_in namelist.py)
- Fortran namelists written by Python scripts

Whereas for the Fortran namelists and the Python dictionaries the user can specify parameters and filenames, the Fortran namelists generated during runtime by the Python scripts do not allow any user interaction.

## Namelist files {#namelist_input_for_extpar_namelist_files}

| **Namelist file**    | **Purpose**                                                      | **Made by script**             | **Used by program**             |
| -------------------- | ---------------------------------------------------------------- | ------------------------------ | -----------------------------   |
| INPUT_grid_org       | define target grid type                                          | runscript                      | extpar_consistency_check        |
|                      |                                                                  |                                | extpar_aot_to_buffer            |
|                      |                                                                  |                                | extpar_landuse_to_buffer        |
|                      |                                                                  |                                | extpar_topo_to_buffer           |
|                      |                                                                  |                                | extpar_cru_to_buffer            |
|                      |                                                                  |                                | extpar_ndvi_to_buffer           |
|                      |                                                                  |                                | extpar_soil_to_buffer           |
|                      |                                                                  |                                | extpar_flake_to_buffer          |
|                      |                                                                  |                                | extpar_isa_to_buffer            |
|                      |                                                                  |                                | extpar_ahf_to_buffer            |
|                      |                                                                  |                                | extpar_emiss_to_buffer          |
|                      |                                                                  |                                | extpar_hwsdART_to_buffer        |
| INPUT_COSMO_GRID     | define target domain for COSMO grid                              | runscript                      | extpar_consistency_check        |
|                      |                                                                  |                                | extpar_aot_to_buffer            |
|                      |                                                                  |                                | extpar_landuse_to_buffer        |
|                      |                                                                  |                                | extpar_topo_to_buffer           |
|                      |                                                                  |                                | extpar_cru_to_buffer            |
|                      |                                                                  |                                | extpar_ndvi_to_buffer           |
|                      |                                                                  |                                | extpar_soil_to_buffer           |
|                      |                                                                  |                                | extpar_flake_to_buffer          |
|                      |                                                                  |                                | extpar_isa_to_buffer            |
|                      |                                                                  |                                | extpar_ahf_to_buffer            |
|                      |                                                                  |                                | extpar_emiss_to_buffer          |
|                      |                                                                  |                                | extpar_hwsdART_to_buffer        |
| INPUT_ICON_GRID      | define target domain for ICON grid                               | runscript                      | extpar_consistency_check        |
|                      |                                                                  |                                | extpar_aot_to_buffer            |
|                      |                                                                  |                                | extpar_landuse_to_buffer        |
|                      |                                                                  |                                | extpar_topo_to_buffer           |
|                      |                                                                  |                                | extpar_cru_to_buffer            |
|                      |                                                                  |                                | extpar_ndvi_to_buffer           |
|                      |                                                                  |                                | extpar_soil_to_buffer           |
|                      |                                                                  |                                | extpar_flake_to_buffer          |
|                      |                                                                  |                                | extpar_isa_to_buffer            |
|                      |                                                                  |                                | extpar_ahf_to_buffer            |
|                      |                                                                  |                                | extpar_emiss_to_buffer          |
| INPUT_ORO            | settings for orography data                                      | runscript                      | extpar_topo_to_buffer           |
| INPUT_OROSMOOTH      | settings for orography smoothing                                 | runscript                      | extpar_topo_to_buffer           |
| INPUT_RADTOPO        | settings for generating topographical shading fields             | runscript                      | extpar_topo_to_buffer           |
| INPUT_SCALE_SEP      | settings to control scale separation for SSO an Z0 calculation   | runscript                      | extpar_topo_to_buffer           |
| INPUT_LU             | settings for landuse data                                        | runscript                      | extpar_landuse_to_buffer        |
| INPUT_AOT            | settings for aerosol data                                        | runscript                      | extpar_aot_to_buffer            |
| INPUT_TCLIM          | settings for temperature data                                    | extpar_cru_-to_buffer          | extpar_consistency_check        |
| INPUT_NDVI           | settings for NDVI data                                           | extpar_ndvi_-to_buffer         | extpar_consistency_check        |
| INPUT_SOIL           | settings for soil data                                           | runscript                      | extpar_soil_to_buffer           |
| INPUT_FLAKE          | settings for lake data                                           | runscript                      | extpar_flake_to_buffer          |
| INPUT_ALB            | settings for albedo data                                         | extpar_albedo_-to_buffer       | extpar_consistency_check        |
| INPUT_ISA            | settings for fraction of impervious surface area data            | extpar_isa_to_buffer           | extpar_consistency_check        |
| INPUT_AHF            | settings for anthropogenic heat flux data                        | extpar_ahf_to_buffer           | extpar_consistency_check        |
| INPUT_EMISS          | settings for emissivity data                                     | extpar_emiss_-to_buffer        | extpar_consistency_check        |
| INPUT_hwsdART        | settings for HWSD USDA data                                      | extpar_hwsdART_-to_buffer      |                                 |
| INPUT_edgar          | settings for EDGAR data                                          | extpar_edgar_-to_buffer        | extpar_consistency_check        |
| INPUT_CDNC           | settings for cdnc data                                           | extpar_cdnc_-to_buffer         | extpar_consistency_check        |
| INPUT_ERA            | settings for ERA data                                            | extpar_era_-to_buffer          | extpar_consistency_check        |
| INPUT_CHECK          | settings for the consistency check                               | runscript                      | extpar_consistency_check        |

## Grid definition {#namelist_input_for_extpar_grid_def}

The specification of the model type (COSMO or ICON) is done in the namelist file INPUT_grid_org, the detailed target grid description for the model domain has to be provided in the namelists files INPUT_COSMO_GRID or INPUT_ICON_GRID.
[\[namelist_target\]]{#namelist_target label="namelist_target"}

### general {#namelist_input_for_extpar_grid_def_general}

### NAMELIST /grid_def/ (INPUT_grid_org) {#namelist-grid_def-input_grid_org .unnumbered}

The namelist /grid_def/ defines the target grid type and the filenames with the namelists of the detailed target grid definition.

| **Parameter**           | **Type**    | **Default**   | **Unit**   | **Description**                                                                 |
| ----------------------- | ----------- | ------------- | ---------- | ------------------------------------------------------------------------------- |
| igrid_type              | integer     |               |            | target grid type, 1 for ICON, 2 for COSMO                                       |
| domain_def_namelist     | character   |               |            | namelist file with domain definition                                            |
| domain_refinement       | character   |               |            | namelist file with domain refinement definition (e.g. for the ICON grid)        |

### Icon {#namelist_input_for_extpar_grid_def_icon}

### NAMELIST /icon_grid_info/ (INPUT_ICON_GRID)  {#namelist-icon_grid_info-input_icon_grid .unnumbered}

The namelist /icon_grid_info/ specifies the filenames and the directory of the Icon grid files with the coordinates of the Icon grid.

| **Parameter**          | **Type**               | **Default**   | **Unit**   | **Description**                                                                 |
| ---------------------- | ---------------------- | ------------- | ---------- | ------------------------------------------------------------------------------- |
| icon_grid_dir          | character              |               |            | path to directory which contains the ICON grid file with the coordinates        |
| icon_grid_nc_file      | character (max_dom)    |               |            | filename of the ICON grid file with the coordinates                             |

### COSMO {#namelist_input_for_extpar_grid_def_cosmo}

### NAMELIST /lmgrid/ (INPUT_COSMO_GRID)  {#namelist-lmgrid-input_cosmo_grid .unnumbered}

The COSMO grid is defined by a rotated latlon-grid.

| **Parameter**   | **Type**   | **Default**   | **Unit**   | **Description**                                                                 |
| --------------- | ---------- | ------------- | ---------- | --------------------------------------------------------------------------------|
| pollon          | real       | -170.         | deg        | longitude of the rotated north pole (in degrees, $E>0$)                         |
| pollat          | real       | 32.5          | deg        | latitude of the rotated north pole (in degrees, $N>0$)                          |
| polgam          | real       | 0\.           | deg        | longitude (in the rotated system) of the geographical north pole                |
| dlon            | real       | 0.08          | deg        | grid point distance in zonal direction (in degrees)                             |
| dlat            | real       | 0.08          | deg        | grid point distance in meridional direction (in degrees)                        |
| startlon_tot    | real       | -1.252        | deg        | transformed longitude of the lower left grid point of the total domain (in degrees, $E>0$) |
| startlat_tot    | real       | -7.972        | deg        | transformed latitude of the lower left grid point of the total domain (in degrees, $N>0$) |
| ie_tot          | integer    | 51            |            | number of grid points in zonal direction                                        |
| je_tot          | integer    | 51            |            | number of grid points in meridional direction                                   |
| ke_tot          | integer    | 0             |            | number of grid points in vertical direction                                     |

## Orography {#namelist_input_for_extpar_orography}

### NAMELIST /oro_runcontrol/ (INPUT_ORO) {#namelist-oro_runcontrol-input_oro .unnumbered}

| **Parameter**    | **Type**   | **Default**   | **Unit**   | **Description**                                                                 |
| ---------------- | ---------- | ------------- | ---------- | ------------------------------------------------------------------------------- |
| lcompute_sgsl    | logical    | .false.       |            | switch to activate subgrid-slope calculation                                    |

### NAMELIST /orography_raw_data/ (INPUT_ORO) {#namelist-orography_raw_data-input_oro .unnumbered}

| **Parameter**                | **Type**    | **Default**                        | **Unit**   | **Description**                                                                 |
| ---------------------------- | ----------- | ---------------------------------- | ---------- | --------------------------------------------------------------------------------|
| itopo_type                   | integer     |                                    |            | switch to choose an orography raw data set; 1 GLOBE, 2 ASTER, 3 MERIT/REMA      |
| lsso_param                   | logical     |                                    |            | switch to choose if SSO parameters should be generated or not                   |
| raw_data_orography_path      | character   |                                    |            | path to orography raw data                                                      |
| ntiles_column                | integer     | GLOBE: 4 ASTER, MERIT/REMA: x      |            | number of tile columns of desired region                                        |
| ntiles_row                   | integer     | GLOBE: 4 ASTER, MERIT/REMA: x      |            | number of tile rows of desired region                                           |
| topo_files                   | character   |                                    |            | filenames of GLOBE (16 tiles) / ASTER (240 tiles)/ MERIT/REMA (72 tiles) raw data sets |
| lsubtract_mean_slope         | logical     | .FALSE. for operational NWP-ICON   |            | treatment of mean slope in computation of SSO parameters for ICON               |

### NAMELIST /orography_io_extpar/ (INPUT_ORO) {#namelist-orography_io_extpar-input_oro .unnumbered}

| **Parameter**             | **Type**    | **Default**   | **Unit**   | **Description**                                                                 |
| ------------------------- | ----------- | ------------- | ---------- | --------------------------------------------------------------------------------|
| orography_buffer_file     | character   |               |            | name for orography buffer file                                                  |
| orography_output_file     | character   |               |            | name for orography output file                                                  |

### NAMELIST /sgsl_io_extpar/ (INPUT_ORO) {#namelist-sgsl_io_extpar-input_oro .unnumbered}

| **Parameter**   | **Type**    | **Default**   | **Unit**   | **Description**                                                                 |
| --------------- | ----------- | ------------- | ---------- | --------------------------------------------------------------------------------|
| lpreproc_oro    | logical     | .false.       |            | read S_ORO from existing NetCDF (.false.) or preprocess from raw topography datasets (.true.) |
| sgsl_files      | character   |               |            | filenames of raw data tiles to be used S_ORO_A10 to S_ORO_P10 (GLOBE) or S_ORO_T001 to S_ORO_T240 |

### NAMELIST /orography_smoothing/ (INPUT_OROSMOOTH) {#namelist-orography_smoothing-input_orosmooth .unnumbered}

| **Parameter**     | **Type**   | **Default**   | **Unit**   | **Description**                                                                 |
| ----------------- | ---------- | ------------- | ---------- | --------------------------------------------------------------------------------|
| lfilter_oro       | logical    | FALSE         |            | **Cosmo-only:** switch for orogaphy smoothing                                   |
| ilow_pass_oro     | integer    | 0             |            | type of orogaphy smoothing and stencil width                                    |
| numfilt_oro       | integer    | 1             |            | number of filter applications                                                   |
| eps_filter        | real       | 10            |            | smoothing parameter ("strength" of the filtering)                               |
| ifill_valley      | integer    | 1             |            | fill valleys before or after oro smoothing (1: before, 2: after)                |
| rfill_valley      | real       | 0             | m          | mask for valley filling (threshold value)                                       |
| ilow_pass_xso     | integer    | 1             |            | type of orogaphy eXtra SmOothing and stencil width (for steep orography)        |
| numfilt_xso       | integer    | 1             |            | number of applications of the eXtra filter                                      |
| lxso_first        | logical    | FALSE         |            | eXtra SmOothing before or after orography smoothing (TRUE/FALSE)                |
| rxso_mask         | real       | 0             | m          | mask for eXtra SmOothing (threshold value)                                      |

### NAMELIST /radtopo/ (INPUT_RADTOPO) {#namelist-radtopo-input_radtopo .unnumbered}

| **Parameter**    | **Type**   | **Default**   | **Unit**   | **Description**                                                                 |
| ---------------- | ---------- | ------------- | ---------- | ------------------------------------------------------------------------------- |
| lradtopo         | logical    |               |            | Switch for radiation corrected topography parameters. Not recommended to use if orographical smoothing is false and smoothing is performed in Int2lm later, because of resulting inconsistencies. |
| nhori            | integer    | 24            |            | Number of horizon angles                                                        |
| radius           | integer    | 40000         | m          | **Icon-only:** Radial distance considered for computation of horizon            |
| min_circ_cov     | integer    | 1             | -          | **Icon-only:** Number of gridcells to be skipped at circumference of circle. A value of 1 considers all points, whereas a value of 5 only consider every fifth point at the circumference. Note that the effect of this switch is dependent on the resolution of the grid as well on the radius choosen. |
| max_missing      | real       | 0.9           | -          | **Icon-only:** Upper limit for fraction of missingness for the horizon parameter. Grid-cells with values above will be set to 0. |
| itype_scaling    | integer    | 2             | -          | **Icon-only:** Power of the caling factor *SIN(horizon-angle)* applied to the geometric skyview factor to account for the anisotropic nature of longwave radiation. |

### NAMELIST /scale_separated_raw_data/ (INPUT_SCALE_SEP) {#namelist-scale_separated_raw_data-input_scale_sep .unnumbered}

| **Parameter**                 | **Type**    | **Default**   | **Unit**   | **Description**                                                                 |
| ----------------------------- | ----------- | ------------- | ---------- | --------------------------------------------------------------------------------|
| lscale_separation             | logical     |               |            | Switch for scale separation. It can only be used in combination with GLOBE as raw data set. |
| raw_data_scale_sep_path       | character   |               |            | path to 3 km filtered topography                                                |
| scale_sep_files               | character   |               |            | filename of 3 km filtered topography                                            |

## Land Use data {#namelist_input_for_extpar_lu}

### NAMELIST /lu_raw_data/ (INPUT_LU) {#namelist-lu_raw_data-input_lu .unnumbered}

| **Parameter**             | **Type**    | **Default**   | **Unit**   | **Description**                                                                 |
| ------------------------- | ----------- | ------------- | ---------- | --------------------------------------------------------------------------------|
| raw_data_lu_path          | character   |               |            | path to landuse data                                                            |
| raw_data_lu_filename      | character   |               |            | filename of landuse raw data                                                    |
| i_landuse_data            | integer     |               |            | switch to choose a land use raw data set; 1 Globcover2009, 2 GLC2000, 3 GLCC, 5 ESA CCI-LC, 6 Ecoclimap-SG |
| l_use_corine              | logical     | .false.       |            | switch to use corine landuse dataset; only possible if i_landuse_data = 1       |
| ilookup_table_lu          | integer     |               |            | switch to choose a lookup table                                                 |
|                           |             |               |            | GLC2000 and GLCC:                                                               |
|                           |             |               |            | 1: operational settings of GME (Ritter, 2007)                                   |
|                           |             |               |            | 2: operational settings of COSMO (Heise, 2005)                                  |
|                           |             |               |            | 3: experimental setting, analog to look-up tables of ECOCLIMAP (Asensio 2010)   |
|                           |             |               |            | GLOBCOVER 2009:                                                                 |
|                           |             |              

## Land Use Data {#namelist_input_for_extpar_lu}

### NAMELIST /lu_raw_data/ (INPUT_LU) {#namelist-lu_raw_data-input_lu .unnumbered}

| **Parameter**          | **Type**   | **Default** | **Unit** | **Description** |
|------------------------|-----------|-------------|----------|----------------------------------------------------------------------------------------------------------------------------------------|
| raw_data_lu_path      | character |             |          | Path to land use data |
| raw_data_lu_filename  | character |             |          | Filename of land use raw data |
| i_landuse_data        | integer   |             |          | Switch to choose a land use raw data set; 1 Globcover2009, 2 GLC2000, 3 GLCC, 5 ESA CCI-LC, 6 Ecoclimap-SG |
| l_use_corine          | logical   | .false.     |          | Switch to use Corine land use dataset; only possible if i_landuse_data = 1 |
| ilookup_table_lu      | integer   |             |          | Switch to choose a lookup table:<br>GLC2000 and GLCC:<br> 1: operational settings of GME (Ritter, 2007)<br> 2: operational settings of COSMO (Heise, 2005)<br> 3: experimental setting, analog to lookup tables of ECOCLIMAP (Asensio 2010)<br>GLOBCOVER 2009:<br> 1: operational settings (Asensio, 2011)<br> 2: experimental settings, analog to lookup tables of ECOCLIMAP (Asensio 2010)<br>ESA CCI-LC:<br> 1: experimental settings (Helmert, 2019)<br>Ecoclimap-SG:<br> 1: Globcover analogue with added LCZs from Oke |
| ntiles_globcover      | integer   | 6           |          | Number of tiles for GLOBCOVER data |
| ncolumn_tiles         | integer   |             |          | Number of columns in tile matrix |
| l_terra_urb          | logical   | .false.     |          | Switch to use TERRA-URB (see [3.2.2](#terra_urb)); only possible if i_landuse_data = 6 |

### NAMELIST /glcc_raw_data/ (INPUT_LU) {#namelist-glcc_raw_data-input_lu .unnumbered}

| **Parameter**            | **Type**   | **Default** | **Unit** | **Description** |
|--------------------------|-----------|-------------|----------|-------------------------------------------------------------------------|
| raw_data_glcc_path      | character |             |          | Path to GLCC data |
| raw_data_glcc_filename  | character |             |          | Filename of GLCC raw data |
| ilookup_table_glcc      | integer   |             |          | Switch to choose a lookup table:<br> 1: operational settings of GME (Ritter, 2007)<br> 2: operational settings of COSMO (Heise, 2005)<br> 3: experimental setting, analog to lookup tables of ECOCLIMAP (Asensio 2010) |

### NAMELIST /glcc_io_extpar/ (INPUT_LU) {#namelist-glcc_io_extpar-input_lu .unnumbered}

| **Parameter**      | **Type**   | **Default** | **Unit** | **Description** |
|--------------------|-----------|-------------|----------|----------------|
| glcc_buffer_file | character |             |          | Name for GLCC buffer file |

## Aerosol Optical Depth {#namelist_input_for_extpar_aot}

### NAMELIST /aerosol_raw_data/ (INPUT_AOT) {#namelist-aerosol_raw_data-input_aot .unnumbered}

| **Parameter**            | **Type**   | **Default** | **Unit** | **Description** |
|--------------------------|-----------|-------------|----------|-----------------------------------------------------------------|
| raw_data_aot_path      | character |             |          | Path to aerosol raw data |
| raw_data_aot_filename  | character |             |          | Filename of aerosol raw data |
| iaot_type              | integer   | 1           |          | Index to specify AOD raw data set: 1: Tegen, 2: AeroCom, 3: MACC-II, 4: MACv2, 5: CAMS |

### NAMELIST /aerosol_io_extpar/ (INPUT_AOT) {#namelist-aerosol_io_extpar-input_aot .unnumbered}

| **Parameter**     | **Type**   | **Default** | **Unit** | **Description** |
|-------------------|-----------|-------------|----------|----------------|
| aot_buffer_file | character |             |          | Name for aerosol buffer file |

## Climatological 2m Temperature {#namelist_input_for_extpar_cru}

### DICT input_tclim (namelist.py) {#dict-input_tclim-namelist.py .unnumbered}

| **Parameter**               | **Type**   | **Default** | **Unit** | **Description** |
|-----------------------------|-----------|-------------|----------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| raw_data_t_clim_path       | character |             |          | Path to T2m climatology data |
| raw_data_t_clim_coarse     | character |             |          | Filename of coarse T2m climatology data |
| raw_data_t_clim_fine       | character |             |          | Filename of fine T2m climatology data |
| it_cl_type                 | integer   |             |          | Switch to choose between the new and fine (1) and the old and coarse over sea and the fine over land (2) raw data set. Note that the fine data set (1) is topographically corrected. |
| t_clim_buffer_file         | character |             |          | Name for T_clim buffer file |
