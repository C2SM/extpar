
# Python modules {#Python modules}

## General workflow {#General workflow}

The general workflow of all Python modules is the same. An exemplary
workflow of the Python modules is described below:

At the beginning of the program information about the environment on
which Extpar is running is written to the logfile and all left-overs
from prior executions of the same Extpar module are deleted. In a next
step each parameter from the namelist 'namelist.py' is checked for
correctness as well as assigned to an internal variable for later use in
the program. The specifaction about the target grid is directly read
from the Fortan namelist-files 'INPUT\_grid\_org. The next step in the
modules involves the generation of all necessary meta-data for the
buffer files and the write of a namelist files in the style of a Fortran
namelist, containing all information needed for the consistency\_check
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

The namelist 'namelist.py' contains the Python dictionaries
'input\_alb', 'input\_tclim', 'input\_emiss', 'input\_ndvi',
'input\_ahf', 'input\_isa' and 'input\_edgar'. These dictionaries
replace their corresponding Fortran namelist files 'INPUT\_'.

'input\_alb' provides information about the albedo data type and the
paths and filenames of the input/output data.

'input\_tclim' contains a switch to determine the type of data (coarse
or fine) as well as the paths and filenames of the input/output data.

'input\_emiss' contains a switch to determine the type of emissivity
data (full range or long-wave) and the filename and paths of the
input/output data.

'input\_ndvi only provides information about the the path and the
filenames of the input/output data.

'input\_era only provides information about the the path and the
filenames of the input/output data.

'input\_isa contains a switch determine the type of ISA data and
provides information about the the path and the filenames of the
input/output data.

'input\_ahf contains a switch determine the type of AHF data and
provides information about the the path and the filenames of the
input/output data.

'input\_edgar only provides information about the the path and the
filenames of the input/output data.

## extpar\_alb\_to\_buffer
-----------------------

### Short description of the subprogram *extpar\_alb\_to\_buffer*

The executable *extpar\_alb\_to\_buffer* allows the aggregation of two
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
data-points. All other steps in extpar\_alb\_to\_buffer are following
the general workflow of the Python scrips.

### Used namelist files and data in-/output

-   namelist files: namelist.py (dict: input\_alb), INPUT\_grid\_org,
    INPUT\_COSMO\_GRID,\
    INPUT\_ICON\_GRID

-   generate namelist: INPUT\_ALB

-   data input: month\_alb.nc, month\_alnid.nc, month\_aluvd.nc,
    global\_soil\_albedo.nc

-   data output: buffer file with albedo data (input\_alb:
    alb\_buffer\_file)

## extpar\_cru\_to\_buffer

### Short description of the subprogram *extpar\_cru\_to\_buffer*

The executable *extpar\_cru\_to\_buffer* aggregates the temperature
climatology of the Climate Research Unit (CRU) to the target grid. The
namelist 'namelist.py' gives the information of the path and the name of
the raw temperature climatology data file. Additionally, the filename
for the buffer is provided. There is an integer switch (it\_cl\_type),
which allows the choice between a newer higher resolved data set for
land surfaces only (1) and an older coarser raw data set for sea
surfaces in combination with the higher resolved data set over land (2).
Aggregation of the coarse data set over land surfaces is no longer
supported since Extpar Version 5.4.

#### Data processing

The data processing with CDO for it\_cl\_type = 1 involves 4 steps:

1.  Set seapoints in the data to missing value.

2.  Extract the fields 'HSURF' from the fine data set.

3.  Merge the fields from step 1 and step 2.

4.  Remap data from step 3 to the target grid using *distance-weighted
    average* interpolation.

The data processing with CDO for it\_cl\_type = 2 involves 5 steps:

1.  Convert coarse data from Celsius to Kelvin, calculate yearly mean
    values and remap coarse data to the grid of the higher resolved data
    set.

2.  Take landpoints from the fine data set and the seapoints from the
    data processed in step 1.

3.  Extract surface height from the buffer file of
    *extpar\_topo\_to\_buffer*

4.  Smooth data processed in step 2 and remap to target grid using
    *distance-weighted average* interpolation.

5.  Correct data processed in step 4 with the surface height extracted
    in step 3.

All subsequent processing on the data follows the general workflow of
the Python scripts.

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input\_tclim), INPUT\_grid\_org,
    INPUT\_COSMO\_GRID,\
    INPUT\_ICON\_GRID

-   generate namelist: INPUT\_TCLIM

-   data input: absolute\_hadcrut3.nc, CRU\_T2M\_SURF\_clim.nc,
    CRU\_T\_SOIL\_clim.nc,\
    orography\_buffer\_file (it\_cl\_type = 2 only)

-   Output: buffer file with CRU data (input\_tclim:
    t\_clim\_buffer\_file)

## extpar\_emiss\_to\_buffer

### Short description of the subprogram *extpar\_emiss\_to\_buffer*

The executable *extpar\_emiss\_to\_buffer* aggregates CAMEL (Combined
ASTER and MODIS Emissivity for Land) data to the target grid. For the
aggregation of the emissivity the namelist 'namelist.py' provides the
path and the file name of the input data. The buffer file name is
defined as well. There exists the integer switch (iemiss\_type) to
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

-   EMISS\_MAX, the maximum EMISS value over 12 months

-   EMISS\_MRAT, the monthly ratio with respect to the maximum EMISS

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input\_emiss) INPUT\_grid\_org,
    INPUT\_COSMO\_GRID,\
    INPUT\_ICON\_GRID

-   generate namelist: INPUT\_EMISS

-   data input: CAM\_bbe\_full\_2010-2015ṅc or CAM\_bbe\_lw\_2010-2015ṅc

-   Output: buffer file with CAMEL data (input\_emiss:
    emiss\_buffer\_file)

## extpar\_ndvi\_to\_buffer

### Short description of the subprogram *extpar\_ndvi\_to\_buffer*

The executable *extpar\_ndvi\_to\_buffer* aggregates NDVI data
(Normalized Differential Vegetation Index) to the target grid. The
namelist 'namelist.py' only contains the path and the file name of the
raw NDVI data. No other parameters can be set.

For the aggregation of the normalized differential vegetation index the
namelist 'namelist.py' is simple. It contains the path and the filename
of the raw data set, as well as the names of the buffer. No other
parameters can be set.

#### Data processing

The remapping to the target grid uses the *first order conservative*
interpolation. After the remapping with CDO two additional fields are
computed:

-   NDVI\_MAX, the maximum NDVI value over 12 months

-   NDVI\_MRAT, the monthly ratio with respect to the maximum NDVI

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input\_ndvi), INPUT\_grid\_org,
    INPUT\_COSMO\_GRID,\
    INPUT\_ICON\_GRID

-   generate namelist: INPUT\_NDVI

-   data input: NDVI\_1998\_2003.nc

-   Output: buffer file with NDVI data (input\_ndvi: ndvi\_buffer\_file)

## extpar\_era\_to\_buffer

### Short description of the subprogram *extpar\_era\_to\_buffer*

The executable *extpar\_era\_to\_buffer* aggregates ERA data (T2M, SST,
W\_SNOW and ORO) to the target grid. It replaces the two NetCDF-Files
generated by ICON-REMAP at DWD. Note that this executable is for
Icon-grids only.

For the aggregation of the ERA climatologies the namelist 'namelist.py'
is simple again. It contains the type of ERA climatology (ERA5 (1) or
ERA-I (2)) the path and the filenames of the raw data sets, as well as
the names of the buffer. No other parameters can be set.

#### Data processing

The remapping to the target grid uses the *first order conservative*
interpolation. After the remapping with CDO the field *W\_SNOW* is
scaled by a factor 1000. No other processing steps take place.

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input\_era), INPUT\_grid\_org,
    INPUT\_COSMO\_GRID,\
    INPUT\_ICON\_GRID

-   generate namelist: INPUT\_ERA

-   data input: ERA5\_ORO\_1990.nc, ERA5\_SD\_1990\_2019.nc,
    ERA5\_SST\_1990\_2019.nc and ERA5\_T2M\_1990\_2019.nc\
    ERA-I\_ORO\_1986.nc, ERA-I\_SD\_1986\_2015.nc,
    ERA-I\_SST\_1986\_2015.nc and ERA-I\_T2M\_1986\_2015

-   Output: buffer file with ERA data (input\_era: era\_buffer\_file)

## extpar\_isa\_to\_buffer

### Short description of the subprogram *extpar\_isa\_to\_buffer*

The executable *extpar\_isa\_to\_buffer* allows the aggregation or
interpolation of data on the fraction of impervious surface area needed
by TERRA\_URB to the target grid.

For the aggregation of the ISA the namelist 'namelist.py' is simple
again. It contains the type of ISA (NOAA (1) or EEA (2)) the path and
the filenames of the raw data sets, as well as the names of the buffer.
No other parameters can be set. Note that the underlying processing does
not differ between different types of ISA

The remapping to the target grid uses the *bilinear* interpolation. No
other processing steps take place.

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input\_isa), INPUT\_grid\_org,
    INPUT\_COSMO\_GRID,\
    INPUT\_ICON\_GRID

-   generate namelist: INPUT\_ISA

-   data input: EEA\_ISA\_16bit\_lonlat.nc( isa\_type=2),
    NOAA\_ISA\_16bit\_lonlat.nc (isa\_type=1)

-   Output: buffer file with ISA data (input\_isa: isa\_buffer\_file)

## extpar\_ahf\_to\_buffer

### Short description of the subprogram *extpar\_ahf\_to\_buffer*

The executable *extpar\_ahf\_to\_buffer* allows the aggregation or
interpolation of data on the anthropogenic heat flux needed by
TERRA\_URB to the target grid.

For the aggregation of the AHF the namelist 'namelist.py' is simple
again. It contains the type of AHF (2.5min (1) or 30sec (2)) the path
and the filenames of the raw data sets, as well as the names of the
buffer. No other parameters can be set. Note that the underlying
processing does not differ between different types of AHF

The remapping to the target grid uses the *bilinear* interpolation. No
other processing steps take place.

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input\_ahf), INPUT\_grid\_org,
    INPUT\_COSMO\_GRID,\
    INPUT\_ICON\_GRID

-   generate namelist: INPUT\_AHF

-   data input: AHF\_2006\_2.5min\_lonlat.nc
    (ahf\_type=1),AHF\_2006\_NOAA\_30sec\_lonlat.nc (ahf\_type=2)

-   Output: buffer file with AHF data (input\_ahf: ahf\_buffer\_file)

## extpar\_edgar\_to\_buffer

### Short description of the subprogram *extpar\_edgar\_to\_buffer*

The executable *extpar\_edgar\_to\_buffer* allows the interpolation of
global emission data for black carbon, organic carbon and sulfur dioxide
needed for the 2D-Aerosol in ICON to the target grid.

The namelist contains only the path to the raw data, the raw data file
names and the name of the buffer file.

The remapping to the target grid uses the *first order conservative*
interpolation. No other processing steps take place.

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input\_edgar), INPUT\_grid\_org,
    INPUT\_ICON\_GRID

-   generate namelist: INPUT\_edgar

-   data input: v8.1\_FT2022\_AP\_NH3\_2022\_TOTALS\_flx.nc,
    v8.1\_FT2022\_AP\_OC\_2022\_TOTALS\_flx.nc,
    v8.1\_FT2022\_AP\_BC\_2022\_TOTALS\_flx.nc,
    v8.1\_FT2022\_AP\_NOx\_2022\_TOTALS\_flx.nc,
    v8.1\_FT2022\_AP\_SO2\_2022\_TOTALS\_flx.nc

-   Output: buffer file with EDGAR data (input\_edgar:
    edgar\_buffer\_file)

## extpar\_cdnc\_to\_buffer

### Short description of the subprogram *extpar\_cdnc\_to\_buffer*

The executable *extpar\_cdnc\_to\_buffer* allows the interpolation of
climatology data for cloud droplet number needed for the Cloud-Aerosol
in ICON to the target grid.

The namelist contains only the path to the raw data, the raw data file
names and the name of the buffer file.

The remapping to the target grid uses the *first order conservative*
interpolation. No other processing steps take place.

### Used namelist files and data in-/output:

-   namelists files: namelist.py (dict: input\_cdnc), INPUT\_grid\_org,
    INPUT\_ICON\_GRID

-   generate namelist: INPUT\_cdnc

-   data input: cdnc\_climatology\_Q06.nc

-   Output: buffer file with cloud droplet number data (input\_cdnc:
    cdnc\_buffer\_file)
