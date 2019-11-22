# EXTPAR: how to compile and run

### Libraries

EXTPAR needs the following libraries:
- NetCDF
- JASPER
- PNG
- Z LIB
Be sure that these libraries are installed on your system
or install them yourself by following the installation
instructions provided with the libraries.

#### On the CSCS machines

All the required libraries are already installed on the CSCS machines. 

#### On the DKRZ machine mistral

Please use the provided script contrib/install-extpar-swstack.sh. This
script requires just an decent installed GCC with gcc, g++, and
gfortran (minimum requirement gcc-6.2).

This script can also be changed to use gcc/intel-18.0.2(Fortran only)
or gcc/NAGfor-6.2(Fortran only). But the main production environment
is supposed to be GCC.

This script fetches all required libraries and tools, compiles and
installs in $HOME/local.gcc.

### Compilation

Simply create a new Options.xxx file (where xxx is the name of your machine) and adapt the compiler name, the flags to use and the library paths in this file.

In case you have a new source file, simply add it to the file Objfiles.
The dependencies between the files are generated automatically.

The build procedure of extpar does not provide a installation target
yet. So binaries are installed in the local extpar directory bin/

#### On the CSCS machines

##### Piz Daint

You have to load the following modules:
    module load cray-netcdf
    module swap PrgEnv-cray PrgEnv-gnu

Then simply type "make" to compile with optimization options or
"make debug" to compile with debug options.
In case you change the target (from debug to opt or the other way around), don't 
forget to issue a "make clean" in between.

##### Kesch

You have to load the following modules:
    module load netcdf
    module load PrgEnv-gnu

Then simply type "make" to compile with optimization options or
"make debug" to compile with debug options.
In case you change the target (from debug to opt or the other way around), don't 
forget to issue a "make clean" in between.

##### Arolla

First, you have to add to your module path:
    export MODULEPATH=$MODULEPATH:/oprusers/osm/opr.arolla/modules/modulefiles
Then load the following modules:
    module load PrgEnv-gnu
    module load netcdf/4.2.1.1-gnu-7.4.0

Then simply type "make" to compile with optimization options or
"make debug" to compile with debug options.
In case you change the target (from debug to opt or the other way around), don't 
forget to issue a "make clean" in between.

#### On the DKRZ machine mistral

Finally, following the libraries and tools installation, extpar can be
installed using the command

MACH=mistral.gcc make

### Run

Some runscript examples are available under /run_scripts.
Just adapt them to your needs!


### Check

1. that you haven't got any error message at runtime,
2. that you got a file in NetCDF format (.nc), and
3. in case of use with cosmo that Int2lm is able to read these 2 files 
