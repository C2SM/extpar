# Overall Description

Numerical Weather Prediction (NWP) models and Climate models require
geographical localized datasets like the topographic height of the earth
surface, the plant cover, the distribution of land and sea and,
dependent on the schemes used, a variety of other external parameters.

The EXTPAR software system (EXTPAR - External Parameter for Numerical
Weather Prediction and Climate Application) is able to generate external
parameters for the different models COSMO and ICON. The software can run
on a UNIX or Linux system where the raw data is stored. It allows
operators (experienced users) running the scripts to create new external
parameters controlled by user specifications like the model domain.

The following steps are performed for the generation of external
parameters:

1.  The target grid has to be specified. The supported target grids are
    -   Rotated and non-rotated longitude-latitude grid (COSMO)
    -   Icosahedral Triangular grids (ICON) with optionally higher
        resolution in selected regions ('local zooming')

2.  The different raw data sets are aggregated to the target grid
    considering all raw data elements which are within the target grid
    element. If the target grid has a higher resolution than the input
    grid on which the raw data is available either an interpolation is
    performed or the target grid is filled with the nearest neighbor,
    but sub-grid scale statistical calculations (e.g. subgrid scale
    variance of orograhic height) are dismissed.

3.  All the different external parameter sets have to be checked for
    consistency against each other. In case of conflicts default values
    are set automatically. In the NetCDF output, information on the
    input data and the processing software is given.

## Input Raw Datasets {#main_input}

The information for the external parameters is aggregated from various
raw datasets for land use data, orography or soil data, see table
below for a detailed list of the raw datasets.

The input data for EXTPAR is stored in a git-LFS repository at
[https://gitlab.dkrz.de/extpar-data/extpar-input-data :material-open-in-new:](https://gitlab.dkrz.de/extpar-data/extpar-input-data){:target="_blank"}. Instructions for
downloading the whole repository or updating with new datasets can be
found in the git-LFS repository. For access to the input data
repository, contact the current EXTPAR source code administrator.

<a name="input_raw_data"></a>

| **Dataset**                                                                    | **Source**                  | **Resolution** |
|---------------------------------------------------------------------------------------|--------------------------------------|------------------------------------------|
| GLOBE orography                                                                       | NOAA/NGDC                            | 30''                |
| ASTER orography (limited domain: 60°N - 60°S)        | METI/NASA                            | 1''                 |
| MERIT/REMA orography                                                                  | Composite DEM                        | 3'' (90m)           |
| Globcover 2009                                                                        | ESA                                  | 10''                |
| GLC2000 land use                                                                      | JRC Ispra                            | 30''                |
| GLCC land use                                                                         | USGS                                 | 30''                |
| Ecoclimap-SG land use                                                                 | CNRS and Meteo France                | 300m                |
| ESA CCI-LC                                                                            | ESA                                  | 10''                |
| DSMW Digital Soil Map of the World                                                    | FAO                                  | 5'                  |
| HWSD Harmonized World Soil Database                                                   | FAO/IIASA/ISRIC/ISSCAS/JRC           | 30''                                     |
| HWSD Harmonized World Soil Database USDA                                              | KIT                                  | 30''                |
| NDVI Climatology, SEAWiFS                                                             | NASA/GSFC                            | 2.5'                |
| CRU near surface climatology                                                          | CRU University of East Anglia        | 0.5 degree          |
| Aerosol Optical thickness                                                             | NASA/GISS                            | 4x5 degree                |
|                                                                                       | (Global Aerosol Climatology Project) |                     |
| AeroCom Global AOD data                                                               | AeroCom Project                      | 1 degree            |
| MACC-II climatological AOD (2003-2012)                                                | ECMWF                                | 1.125 degree        |
| MACv2 monthly AOD, SSA and ASY data                                                   | MPI, RHM                             | 1 degree            |
| CAMS monthly 3D-climatology                                                           |                                      |                     |
| 11 types of aerosols                                                                  | ECMWF, RHM                           | 3 degree            |
| Global lake database (GLDB)                                                           | DWD/RSHU/MeteoFrance                 | 30''                |
| MODIS albedo                                                                          | NASA                                 | 5'                  |
| MODIS derived soil albedo values                                                      | Community Land Model 3.5             | 30'                 |
| CAMEL Emissivity                                                                      | NASA                                 | 5km                 |
| EDGAR Emissions                                                                       | European Commission /JRC/PBL         | 0.1 degree          |
| MODIS cloud droplet number climatology Q06                                            | NASA                                 | 1 degree            |

<center>*Table 1: Input raw datasets*</center>



## Output External Parameters {#main_output}

The output fields with the external parameters are shown here:

<a name="extpar_output"></a>

| **External parameter**                                                         | **Short name**      | **Unit**        | **Raw dataset**          |
|---------------------------------------------------------------------------------------|------------------------------|--------------------------|----------------------------------------------------|
| geometrical height                                                                    | HSURF                        | $m$                      | GLOBE/ASTER/MERIT/REMA                        |
| geopotential of earth surface                                                         | FIS                          | $m^{2} s^{-1}$          | GLOBE/ASTER/MERIT/REMA                        |
| standard deviation of subgrid scale orographic height                                 | SSO_STDH                  | $m$                      | GLOBE/ASTER/MERIT/REMA                        |
| anisotropy of topography                                                              | SSO_GAMMA                 | 1                        | GLOBE/ASTER/MERIT/REMA                        |
| angle between principal axis of orography and global E                                | SSO_THETA                 | 1                        | GLOBE/ASTER/MERIT/REMA                                       |
| mean slope of subgrid scale orography                                                 | SSO_SIGMA                 | 1                        | GLOBE/ASTER/MERIT/REMA                                       |
| surface roughness                                                                     | Z0                           | $m$                      | GLC2000, GLOBE/ASTER/MERIT/REMA                              |
| Slope aspect                                                                          | SLOPE_ASP                 | deg                      | GLOBE/ASTER/MERIT/REMA                                       |
| Slope angle                                                                           | SLOPE_ANG                 | deg                      | GLOBE/ASTER/MERIT/REMA                                       |
| Horizon angles (resolution from 15deg)                                                | HORIZON                      | deg                      | GLOBE/ASTER/MERIT/REMA                                       |
| Skyview factor                                                                        | SKYVIEW                      | -                        | GLOBE/ASTER/MERIT/REMA                                       |
| soil texture                                                                          | SOILTYP                      | -                        | DSMW/HWSD                                          |
| fraction of sand                                                                      | FR_SAND                   | %                       | HWSD                                               |
| fraction of silt                                                                      | FR_SILT                   | %                       | HWSD                                               |
| fraction of clay                                                                      | FR_CLAY                   | %                       | HWSD                                               |
| fraction of organic carbon                                                            | FR_OC                     | %                       | HWSD                                               |
| bulk density                                                                          | BULK_DENS                 | $g cm^{-3}$              | HWSD                                               |
| deep soil texture                                                                     | SUBSOILTYP                   | -                        | HWSD                                               |
| deep soil fraction of sand                                                            | SUB_FR_SAND            | %                       | HWSD                                               |
| deep soil fraction of silt                                                            | SUB_FR_SILT            | %                       | HWSD                                               |
| deep soil fraction of clay                                                            | SUB_FR_CLAY            | %                       | HWSD                                               |
| deep soil fraction of organic carbon                                                  | SUB_FR_OC              | %                       | HWSD                                               |
| deep soil bulk density                                                                | SUB_BULK_DENS          | $g cm^{-3}$              | HWSD                                               |
| Fraction of Heavy Clay                                                                | fr_hcla                     | 1                        | HWSD_USDA                    |
| Fraction of Silty Clay                                                                | fr_silc                     | 1                        | HWSD_USDA                    |
| Fraction of Light Clay                                                                | fr_lcla                     | 1                        | HWSD_USDA                    |
| Fraction of Silty Clay Loam                                                           | fr_sicl                     | 1                        | HWSD_USDA                    |
| Fraction of Clay Loam                                                                 | fr_cloa                     | 1                        | HWSD_USDA                    |
| Fraction of Silt                                                                      | fr_silt                     | 1                        | HWSD_USDA                    |
| Fraction of Silty Loam                                                                | fr_silo                     | 1                        | HWSD_USDA                    |
| Fraction of Sandy Clay                                                                | fr_scla                     | 1                        | HWSD_USDA                    |
| Fraction of Loam                                                                      | fr_loam                     | 1                        | HWSD_USDA                    |
| Fraction of Sandy Clay Loam                                                           | fr_sclo                     | 1                        | HWSD_USDA                    |
| Fraction of Sandy Loam                                                                | fr_sloa                     | 1                        | HWSD_USDA                    |
| Fraction of Loamy Sand                                                                | fr_lsan                     | 1                        | HWSD_USDA                    |
| Fraction of Sand                                                                      | fr_sand                     | 1                        | HWSD_USDA                    |
| Fraction of Undefined or Water                                                        | fr_udef                     | 1                        | HWSD_USDA                    |
| ground fraction covered by plants max (vegetation period)                             | PLCOV_MX                  | 1                        | GLC2000/Globcover/ ESA CCI-LC                      |
| ground fraction covered by plants min (vegetation period)                             | PLCOV_MN                  | 1                        | GLC2000/Globcover/ ESA CCI-LC                      |
| ground fraction covered by artificial (urban) areas                                   | URBAN                        | 1                        | GLC2000/Globcover/ ESA CCI-LC                      |
| ground fraction covered by artificial (urban) areas                                   | URBAN                        | 1                        | GLC2000/Globcover/ESA CCI-LC/ LCZs with TERRA_URB |
| ground fraction covered by deciduous forest                                           | FOR_D                     | 1                        | GLC2000/Globcover/ ESA CCI-LC                      |
| skin conductivity                                                                     | SKC                          | $W m^{-1} K^{-1}$        | Globcover/ESA CCI-LC                               |
| root depth                                                                            | ROOTDP                       | $m$                      | GLC2000/Globcover/ ESA CCI-LC                      |
| leaf area index max(vegetation period)                                                | LAI_MX                    | 1                        | GLC2000/Globcover/ESA CCI-LC                       |
| leaf area index min (vegetation period)                                               | LAI_MN                    | 1                        | GLC2000/Globcover/ ESA CCI-LC                      |
| plant resistance                                                                      | PRS_MIN                   | $s m^{-1}$               | GLC2000/Globcover/ ESA CCI-LC                      |
| long wave surface emissivity                                                          | EMISS_RAD                 | 1                        | GLC2000/Globcover/ ESA CCI-LC                      |
| (monthly) normalized differential vegetation index                                    | NDVI                         | 1                        | SEAWIFS                                            |
| Annual maximum of normalized differential vegetation index                            | NDVI_MAX                  | 1                        | SEAWIFS                                            |
| (monthly) proportion of actual value/ maximum normalized                              |                              |                          |                                                    |
| differential vegetation index                                                         | NDVI_RATIO                | 1                        | SEAWIFS                                            |
| (monthly) optical thickness from black carbon aerosol                                 | AER_BC                    | 1                        | GACP                                               |
| (monthly) optical thickness from dust aerosol                                         | AER_DUST                  | 1                        | GACP                                               |
| (monthly) optical thickness from organic aerosol                                      | AER_ORG                   | 1                        | GACP                                               |
| (monthly) optical thickness from SO4 aerosol                                          | AER_SO4                   | 1                        | GACP                                               |
| (monthly) optical thickness from sea salt aerosol                                     | AER_SS                    | 1                        | GACP                                               |
| (monthly) aerosol optical thickness for RG92 spectral bands                           | AOT12                        | 1                        | MACv2                                              |
| (monthly) single scattering albedo for RG92 spectral bands                            | SSA12                        | 1                        | MACv2                                              |
| (monthly) asymmetry factor for RG92 spectral bands                                    | ASY12                        | 1                        | MACv2                                              |
| (monthly) layer-integrated mass of Sea Salt with                                      |                              |                          |                                                    |
| dry radius in the range 0.03-0.5 microns                                              | Sea_Salt_bin1              | $kg m^{-2}$              | CAMS                                               |
| (monthly) layer-integrated mass of Sea Salt with                                      |                              |                          |                                                    |
| dry radius in the range 0.5-5.0  microns                                              | Sea_Salt_bin2              | $kg m^{-2}$              | CAMS                                               |
| (monthly) layer-integrated mass of Sea Salt with                                      |                              |                          |                                                    |
| dry radius in the range 5.0-20.0  microns                                             | Sea_Salt_bin3              | $kg m^{-2}$              | CAMS                                               |
| (monthly) layer-integrated mass of Mineral Dust with                                  |                              |                          |                                                    |
| dry radius in the range 0.03-0.55 microns                                             | Mineral_Dust_bin1          | $kg m^{-2}$              | CAMS                                               |
| (monthly) layer-integrated mass of Mineral Dust with                                  |                              |                          |                                                    |
| dry radius in the range 0.55-0.9  microns                                             | Mineral_Dust_bin2          | $kg m^{-2}$              | CAMS                                               |
| (monthly) layer-integrated mass of Mineral Dust with                                  |                              |                          |                                                    |
| dry radius in the range 0.9-20.0  microns                                             | Mineral_Dust_bin3          | $kg m^{-2}$              | CAMS                                               |
| (monthly) layer-integrated mass of hydrophilic                                        |                              |                          |                                                    |
| Organic Matter                                                                        | Organic_Matter_hydrophilic | $kg m^{-2}$              | CAMS                                               |
| (monthly) layer-integrated mass of hydrophobic                                        |                              |                          |                                                    |
| Organic Matter                                                                        | Organic_Matter_hydrophobic | $kg m^{-2}$              | CAMS                                               |
| (monthly) layer-integrated mass of hydrophilic                                        |                              |                          |                                                    |
| Black Carbon                                                                          | Black_Carbon_hydrophilic   | $kg m^{-2}$              | CAMS                                               |
| (monthly) layer-integrated mass of hydrophobic                                        |                              |                          |                                                    |
| Black Carbon                                                                          | Black_Carbon_hydrophobic   | $kg m^{-2}$              | CAMS                                               |
| (monthly) layer-integrated mass of Sulfates                                           | Sulfates                     | $kg m^{-2}$              | CAMS                                               |
| (monthly) Pressure at base of layer                                                   | half_level_pressure        | Pa                       | CAMS                                               |
| Near surface temperature (climatological mean)                                        | T_2M_CL                    | $K$                    | CRU                                                |
| Lake Depth                                                                            | DEPTH_LK                  | $m$                    | GLDB                                               |
| Lake Fraction                                                                         | FR_LAKE                   | 1                        | GLDB                                               |
| (monthly) albedo                                                                      | ALB_DIF12                   | %                       | MODIS                                              |
| (monthly) Near Infrared Albedo                                                        | ALNID                        | %                       | MODIS                                              |
| (monthly) Ultra Violet Albedo                                                         | ALUVD                        | %                       | MODIS                                              |
| soil albedo for dry soils                                                             | ALB_DRY                   | %                       | Community Land Model 3.5                           |
| soil albedo for saturated soils                                                       | ALB_SAT                   | %                       | Community Land Model 3.5                           |
| fraction of impervious surface area                                                   | ISA                          | 1                        | NOAA, EEA or LCZs with TERRA_URB                  |
| anthropogenic heat flux                                                               | AHF                          | $W m^{-2}$               | NOAA or LCZs with TERRA_URB                       |
| subgrid-scale slope parameter                                                         | S_ORO                       | 1                        | GLOBE, ASTER,                                      |
|                                                                                       |                              |                          | MERIT/REMA                    |
| EMISS yearly maximum for climatology 1998-2003                                        | EMISS_MAX                   | 1                        | CAMEL                                              |
| monthly mean EMISS climatology 1998-2003                                              | EMISS                        | 1                        | CAMEL                                              |
| (monthly) proportion of actual value/maximum normalized differential vegetation index | EMISS_MRAT                  | 1                        | CAMEL                                              |
| Urban paved fraction                                                                  | FR_PAVED                    | 1                        | LCZs with TERRA_URB          |
| Urban building fraction                                                               | URB_BLDFR                   | 1                        | LCZs with TERRA_URB          |
| Urban building height                                                                 | URB_BLDH                    | $m$                      | LCZs with TERRA_URB          |
| Urban canyon height-to-width ratio                                                    | URB_H2W                     | 1                        | LCZs with TERRA_URB          |
| Urban shortwave albedo                                                                | URB_SALB                    | 1                        | LCZs with TERRA_URB          |
| Urban thermal albedo                                                                  | URB_TALB                    | 1                        | LCZs with TERRA_URB          |
| Urban emissivity                                                                      | URB_EMIS                    | 1                        | LCZs with TERRA_URB          |
| Urban heat conductivity                                                               | URB_HCON                    | 1                        | LCZs with TERRA_URB          |
| Urban heat capacity                                                                   | URB_HCAP                    | $J/K$                    | LCZs with TERRA_URB          |
| Annual black carbon emissions                                                         | emi_bc                      | $kg\,m^{-2}\,s^{-1}$ | EDGAR                         |
| Annual organic carbon emissions                                                       | emi_oc                      | $kg\,m^{-2}\,s^{-1}$ | EDGAR                         |
| Annual sulfur dioxide carbon emissions                                                | emi_so2                     | $kg\,m^{-2}\,s^{-1}$ | EDGAR                         |
| Annual ammonia emissions                                                              | emi_nh3                     | $kg\,m^{-2}\,s^{-1}$ | EDGAR                         |
| Annual nitrogen oxides emissions                                                      | emi_nox                     | $kg\,m^{-2}\,s^{-1}$ | EDGAR                         |
| Monthly cloud droplet number climatology                                              | cdnc                         | $cm^{-3}$            | MODIS                         |

<center>*Table 2: Output external parameters*</center>