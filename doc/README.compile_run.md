# How to clone, compile and run
## Code structure
EXTPAR contains FORTAN-Code as well as hybrid Python-CDO scripts.
Both code bases need external libraries and installations.
### Libraries
#### Fortran
EXTPAR needs the following libraries for the Fortran-Code:
- NetCDF
- JASPER
- PNG
- Z LIB
#### Python-CDO
EXTPAR needs the following Python-packages and installations:
- at least Python 3.6
- module netCDF4

The module netCDF is the Python interface to the netCDF C library. It allows the user create and manipulate netCDF files with Python.
For more detailed information please visit [netCDF Python](https://unidata.github.io/netcdf4-python/netCDF4/index.html).

Additionally an installation of CDO (Climate Data Operator) is required. All necessary information about this tool can be found at [CDO-MPI]( https://code.mpimet.mpg.de/projects/cdo/).

Be sure that these libraries are installed on your system
or install them yourself by following the installation
instructions provided with the libraries.

#### On the CSCS machines

All the required libraries are already installed on the CSCS machines.  
Ready-to-use Python environments for Tsa and Daint can be found at */project/g110/extpar_envs/venv_jenkins_your_machine*.  
Just *source /project/g110/extpar_envs/venv_jenkins_your_machine/bin/activate* in your Extpar runscript to activate these Python environments.

#### On the DKRZ machine mistral

Please use the provided script contrib/install-extpar-swstack.sh. This
script requires just an decent installed GCC with gcc, g++, and
gfortran (minimum requirement gcc-6.2).

This script can also be changed to use gcc/intel-18.0.2(Fortran only)
or gcc/NAGfor-6.2(Fortran only). But the main production environment
is supposed to be GCC.

This script fetches all required libraries and tools, compiles and
installs in $HOME/local.gcc.

### Cloning from GitHub
Because of the embedded submodule CDI in Extpar, please clone Extpar from GitHub using the following commands:

* git clone --recursive git@github.com:C2SM-RCM/extpar.git
* git submodule update

### Compilation
Since Version 5.4, Extpar is built with an autotool based build-system. This has been necessary to accomodate for the additional C source code files and newly required libraries.

There are three options to compile Extpar: 

#### In-source build
./configure.*hostname.compiler*  
source modules.env  
make   

#### Out-of-source build
mkdir build-*my_self_defined_note*  
cd build-*my_self_defined_note*  
*path_to_the_extpar_installation*/configure.*hostname.compiler*  
source modules.env  
make  

#### Install binaries only in external directory
./configure.*hostname.compiler* --prefix=*my_external_directory*  
source modules.env  
make install  
The binaries will be installed in  *my_external_directory/bin*

#### Restart build from scratch
make distclean 

#### Configure for new machines
The first step in creating a new machine setup is to take one of
the existing configure wrapper scripts and adapt it to your local
environment. The scripts are called configure.*hostname.compiler*.

### Run
Some runscript examples are available under /run_scripts.
Just adapt them to your needs!

### Check
1. that you haven't got any error message at runtime,
2. that you got a file in NetCDF format (.nc), and
3. in case of use with cosmo that Int2lm is able to read these 2 files 
