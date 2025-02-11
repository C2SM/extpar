# Dependencies

EXTPAR contains Fortran Code as well as hybrid Python-CDO scripts.
Both code bases need external libraries and installations.

## Libraries

### Fortran

EXTPAR needs the following libraries for the Fortran-Code:

- NetCDF
- JASPER
- PNG
- Z LIB

### Python-CDO

EXTPAR needs the following Python packages and installations:

- at least Python 3.6
- module netCDF4

The module netCDF is the Python interface to the netCDF C library. It allows the user create and manipulate netCDF files with Python.
For more detailed information please visit [netCDF Python :material-open-in-new:](https://unidata.github.io/netcdf4-python/){:target="_blank"}.

Additionally, an installation of CDO (Climate Data Operators) is required. All necessary information about this tool can be found at [CDO-MPI]( https://code.mpimet.mpg.de/projects/cdo/).

Be sure that these libraries are installed on your system
or install them yourself by following the installation
instructions provided with the libraries.

#### On the CSCS machines

All the required libraries are already installed on the CSCS machines.  

#### On the DKRZ machine Levante

All the required libraries are already installed on the DKRZ machine.

