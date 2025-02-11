
# Python modules {#Python modules}

## General workflow {#General workflow}

The general workflow of all Python modules is the same. An exemplary
workflow of the Python modules is described below:

At the beginning of the program information about the environment on
which EXTPAR is running is written to the logfile and all left-overs
from prior executions of the same EXTPAR module are deleted. In a next
step each parameter from the namelist `namelist.py` is checked for
correctness as well as assigned to an internal variable for later use in
the program. The specifaction about the target grid is directly read
from the Fortan namelist-files `INPUT_grid_org`. The next step in the
modules involves the generation of all necessary meta-data for the
buffer files and the write of a namelist files in the style of a Fortran
namelist, containing all information needed for the consistency_check
at the end. In case of a COSMO target grid, a grid specification file is
written, that is later used by CDO for the interpolation.

After this is all setup, the most compute-intensive parts like the
remapping to the target grid or data modifications are done using CDO.
The Python program launches a subshell executing the CDO in it. The
output from CDO is reshaped in order to fit the dimensions of the buffer
files. After the reshape, the fields and its corresponding metadata is
written to a netCDF buffer file. The last step of the programme again
deletes all intermediate netCDF or other files written during runtime,
that do not serve any purpose.

Module-specific modifications or additional computations are described
in the paragraph *Data processing* of each Python module.

### Namelist

The namelist `namelist.py` contains the Python dictionaries
`input_alb`, `input_tclim`, `input_emiss`, `input_ndvi`,
`input_ahf`, `input_isa` and `input_edgar`. These dictionaries
replace their corresponding Fortran namelist files `INPUT_`.

`input_alb` provides information about the albedo data type and the
paths and filenames of the input/output data.

`input_tclim` contains a switch to determine the type of data (coarse
or fine) as well as the paths and filenames of the input/output data.

`input_emiss` contains a switch to determine the type of emissivity
data (full range or long-wave) and the filename and paths of the
input/output data.

`input_ndvi` only provides information about the the path and the
filenames of the input/output data.

`input_era` only provides information about the the path and the
filenames of the input/output data.

`input_isa` contains a switch determine the type of ISA data and
provides information about the the path and the filenames of the
input/output data.

`input_ahf` contains a switch determine the type of AHF data and
provides information about the the path and the filenames of the
input/output data.

`input_edgar` only provides information about the the path and the
filenames of the input/output data.

## extpar_alb_to_buffer
-----------------------

### Short description of the subprogram *extpar_alb_to_buffer*

The executable *extpar_alb_to_buffer* allows the aggregation of two
different kinds of albedo data to the target grid. The first kind is a
climatology of monthly values of total albedo derived from MODIS
satellite data for the 3 spectral bands visible, near infrared and
ultraviolet. The second kind contains information for soil albedo only
in dry and saturated conditions. It originates from the Community Land
Model[^13].

#### Data processing

The data is remapped to the target grid using the *distance-weighted
average* interpolation. CDO first generates the weights for the
interpolation from one of the input files and then applies these weights
to all input files. After the interpolation took place, all values in
the range of -100000 - 0.02 are set to 0.02 to prevent unrealistic
data-points. All other steps in extpar_alb_to_buffer are following
the general workflow of the Python scrips.

### Used namelist files and data in-/output

-   namelist files: namelist.py (dict: input_alb), INPUT_grid_org,
    INPUT_COSMO_GRID,
    INPUT_ICON_GRID

-   generate namelist: INPUT_ALB

-   data input: month_alb.nc, month_alnid.nc, month_aluvd.nc,
    global_soil_albedo.nc

-   data output: buffer file with albedo data (input_alb:
    alb_buffer_file)

## extpar_cru_to_buffer

### Short description of the subprogram *extpar_cru_to_buffer*

The executable *extpar_cru_to_buffer* aggregates the temperature
climatology of the Climate Research Unit (CRU) to the target grid. The
namelist `namelist.py` gives the information of the path and the name of
the raw temperature climatology data file. Additionally, the filename
for the buffer is provided. There is an integer switch (it_cl_type),
which allows the choice between a newer higher resolved data set for
land surfaces only (1) and an older coarser raw data set for sea
surfaces in combination with the higher resolved data set over land (2).
Aggregation of the coarse data set over land surfaces is no longer
supported since EXTPAR Version 5.4.

#### Data processing

The data processing with CDO for it_cl_type = 1 involves 4 steps:

1.  Set seapoints in the data to missing value.

2.  Extract the fields `HSURF` from the fine data set.

3.  Merge the fields from step 1 and step 2.

4.  Remap data from step 3 to the target grid using *distance-weighted
    average* interpolation.

The data processing with CDO for it_cl_type = 2 involves 5 steps:

1.  Convert coarse data from Celsius to Kelvin, calculate yearly mean
    values and remap coarse data to the grid of the higher resolved data
    set.

2.  Take landpoints from the fine data set and the seapoints from the
    data processed in step 1.

3.  Extract surface height from the buffer file of
    *extpar_topo_to_buffer*

4.  Smooth data processed in step 2 and remap to target grid using
    *distance-weighted average* interpolation.

5.  Correct data processed in step 4 with the surface height extracted
    in step 3.

All subsequent processing on the data follows the general workflow of
the Python scripts.

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input_tclim), INPUT_grid_org,
    INPUT_COSMO_GRID,
    INPUT_ICON_GRID

-   generate namelist: INPUT_TCLIM

-   data input: absolute_hadcrut3.nc, CRU_T2M_SURF_clim.nc,
    CRU_T_SOIL_clim.nc,
    orography_buffer_file (it_cl_type = 2 only)

-   Output: buffer file with CRU data (input_tclim:
    t_clim_buffer_file)

## extpar_emiss_to_buffer

### Short description of the subprogram *extpar_emiss_to_buffer*

The executable *extpar_emiss_to_buffer* aggregates CAMEL (Combined
ASTER and MODIS Emissivity for Land) data to the target grid. For the
aggregation of the emissivity the namelist `namelist.py` provides the
path and the file name of the input data. The buffer file name is
defined as well. There exists the integer switch (iemiss_type) to
determine whether one wants to use the broad band emissivity for the 3.6
and 14.3 micron spectral range (1) or the broad band emissivity between
8.0 and 13.5 micron spectral range (2).

#### Data processing

After the generation of the interpolation weights artificial low values
below 0.5 are set to -999. In a subsequent processing step -999 is set
to the value for missing data. In order to not have missing data in the
field to interpolate, all missing values are set to the values of its
nearest neighbour. The last step involves the *first order conservative*
interpolation to the target grid. After the remapping with CDO two
additional fields are computed:

-   EMISS_MAX, the maximum EMISS value over 12 months

-   EMISS_MRAT, the monthly ratio with respect to the maximum EMISS

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input_emiss) INPUT_grid_org,
    INPUT_COSMO_GRID,
    INPUT_ICON_GRID

-   generate namelist: INPUT_EMISS

-   data input: CAM_bbe_full_2010-2015ṅc or CAM_bbe_lw_2010-2015ṅc

-   Output: buffer file with CAMEL data (input_emiss:
    emiss_buffer_file)

## extpar_ndvi_to_buffer

### Short description of the subprogram *extpar_ndvi_to_buffer*

The executable *extpar_ndvi_to_buffer* aggregates NDVI data
(Normalized Differential Vegetation Index) to the target grid. The
namelist `namelist.py` only contains the path and the file name of the
raw NDVI data. No other parameters can be set.

For the aggregation of the normalized differential vegetation index the
namelist `namelist.py` is simple. It contains the path and the filename
of the raw data set, as well as the names of the buffer. No other
parameters can be set.

#### Data processing

The remapping to the target grid uses the *first order conservative*
interpolation. After the remapping with CDO two additional fields are
computed:

-   NDVI_MAX, the maximum NDVI value over 12 months

-   NDVI_MRAT, the monthly ratio with respect to the maximum NDVI

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input_ndvi), INPUT_grid_org,
    INPUT_COSMO_GRID,
    INPUT_ICON_GRID

-   generate namelist: INPUT_NDVI

-   data input: NDVI_1998_2003.nc

-   Output: buffer file with NDVI data (input_ndvi: ndvi_buffer_file)

## extpar_era_to_buffer

### Short description of the subprogram *extpar_era_to_buffer*

The executable *extpar_era_to_buffer* aggregates ERA data (T2M, SST,
W_SNOW and ORO) to the target grid. It replaces the two NetCDF-Files
generated by ICON-REMAP at DWD. Note that this executable is for
Icon-grids only.

For the aggregation of the ERA climatologies the namelist `namelist.py`
is simple again. It contains the type of ERA climatology (ERA5 (1) or
ERA-I (2)) the path and the filenames of the raw data sets, as well as
the names of the buffer. No other parameters can be set.

#### Data processing

The remapping to the target grid uses the *first order conservative*
interpolation. After the remapping with CDO the field *W_SNOW* is
scaled by a factor 1000. No other processing steps take place.

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input_era), INPUT_grid_org,
    INPUT_COSMO_GRID,
    INPUT_ICON_GRID

-   generate namelist: INPUT_ERA

-   data input: ERA5_ORO_1990.nc, ERA5_SD_1990_2019.nc,
    ERA5_SST_1990_2019.nc and ERA5_T2M_1990_2019.nc
    ERA-I_ORO_1986.nc, ERA-I_SD_1986_2015.nc,
    ERA-I_SST_1986_2015.nc and ERA-I_T2M_1986_2015

-   Output: buffer file with ERA data (input_era: era_buffer_file)

## extpar_isa_to_buffer

### Short description of the subprogram *extpar_isa_to_buffer*

The executable *extpar_isa_to_buffer* allows the aggregation or
interpolation of data on the fraction of impervious surface area needed
by TERRA_URB to the target grid.

For the aggregation of the ISA the namelist `namelist.py` is simple
again. It contains the type of ISA (NOAA (1) or EEA (2)) the path and
the filenames of the raw data sets, as well as the names of the buffer.
No other parameters can be set. Note that the underlying processing does
not differ between different types of ISA

The remapping to the target grid uses the *bilinear* interpolation. No
other processing steps take place.

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input_isa), INPUT_grid_org,
    INPUT_COSMO_GRID,
    INPUT_ICON_GRID

-   generate namelist: INPUT_ISA

-   data input: EEA_ISA_16bit_lonlat.nc( isa_type=2),
    NOAA_ISA_16bit_lonlat.nc (isa_type=1)

-   Output: buffer file with ISA data (input_isa: isa_buffer_file)

## extpar_ahf_to_buffer

### Short description of the subprogram *extpar_ahf_to_buffer*

The executable *extpar_ahf_to_buffer* allows the aggregation or
interpolation of data on the anthropogenic heat flux needed by
TERRA_URB to the target grid.

For the aggregation of the AHF the namelist `namelist.py` is simple
again. It contains the type of AHF (2.5min (1) or 30sec (2)) the path
and the filenames of the raw data sets, as well as the names of the
buffer. No other parameters can be set. Note that the underlying
processing does not differ between different types of AHF

The remapping to the target grid uses the *bilinear* interpolation. No
other processing steps take place.

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input_ahf), INPUT_grid_org,
    INPUT_COSMO_GRID,
    INPUT_ICON_GRID

-   generate namelist: INPUT_AHF

-   data input: AHF_2006_2.5min_lonlat.nc
    (ahf_type=1),AHF_2006_NOAA_30sec_lonlat.nc (ahf_type=2)

-   Output: buffer file with AHF data (input_ahf: ahf_buffer_file)

## extpar_edgar_to_buffer

### Short description of the subprogram *extpar_edgar_to_buffer*

The executable *extpar_edgar_to_buffer* allows the interpolation of
global emission data for black carbon, organic carbon and sulfur dioxide
needed for the 2D-Aerosol in ICON to the target grid.

The namelist contains only the path to the raw data, the raw data file
names and the name of the buffer file.

The remapping to the target grid uses the *first order conservative*
interpolation. No other processing steps take place.

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input_edgar), INPUT_grid_org,
    INPUT_ICON_GRID

-   generate namelist: INPUT_edgar

-   data input: v8.1_FT2022_AP_NH3_2022_TOTALS_flx.nc,
    v8.1_FT2022_AP_OC_2022_TOTALS_flx.nc,
    v8.1_FT2022_AP_BC_2022_TOTALS_flx.nc,
    v8.1_FT2022_AP_NOx_2022_TOTALS_flx.nc,
    v8.1_FT2022_AP_SO2_2022_TOTALS_flx.nc

-   Output: buffer file with EDGAR data (input_edgar:
    edgar_buffer_file)

## extpar_cdnc_to_buffer

### Short description of the subprogram *extpar_cdnc_to_buffer*

The executable *extpar_cdnc_to_buffer* allows the interpolation of
climatology data for cloud droplet number needed for the Cloud-Aerosol
in ICON to the target grid.

The namelist contains only the path to the raw data, the raw data file
names and the name of the buffer file.

The remapping to the target grid uses the *first order conservative*
interpolation. No other processing steps take place.

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input_cdnc), INPUT_grid_org,
    INPUT_ICON_GRID

-   generate namelist: INPUT_cdnc

-   data input: cdnc_climatology_Q06.nc

-   Output: buffer file with cloud droplet number data (input_cdnc:
    cdnc_buffer_file)
