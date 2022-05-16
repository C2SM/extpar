# EXTPAR

# General Information
EXTPAR (External Parameters for Numerical Weather Prediction and Climate Application) is an official software of the [COSMO Consortium](www.cosmo-model.org).  It is used to prepare the external parameter data files that are used as input for the COSMO and the ICON model.

The code is written in Fortran 90 and in Python. The Python scripts use CDO for the most compute-intensive parts.  Currently the code is tested regularly using the gcc, NAG, and Intel Fortran compilers.  The code is also accelerated in some places with OpenMP parallelization.  

The code once compiled generates 6 Fortran executables and 7 Python scripts, which can be run simultaneously except for the final extpar_consistency_check.exe, which is used to tie together all the external parameter results into one output file.  


Information about the latest changes can be found in the [Release Notes](ReleaseNotes.md)

A full documentation of code can be found in the [manual](doc/user_and_implementation_manual.pdf)

# Quick Start
### Daint

```
git clone --recursive git@github.com:C2SM-RCM/extpar.git
cd extpar
git submodule update
./configure.daint.gcc
source modules.env
make -j 4
module load daint-gpu
module load CDO
source /project/g110/extpar/venv_daint/bin/activate
export PYTHONPATH=$PYTHONPATH:$(pwd)/python/lib
```

### Tsa

```
git clone --recursive git@github.com:C2SM-RCM/extpar.git
cd extpar
git submodule update
./configure.tsa.gcc
source modules.env
make -j 4
source /oprusers/osm/.opr_setup_dir
export MODULEPATH=$MODULEPATH\:$OPR_SETUP_DIR/modules/modulefiles
source /project/g110/extpar/venv_tsa/bin/activate
module load cdo
export PYTHONPATH=$PYTHONPATH:$(pwd)/python/lib
```

### Mistral

```
git clone --recursive git@github.com:C2SM-RCM/extpar.git
cd extpar
git submodule update
./configure.mistral.gcc # or ./configure.mistral.intel
source modules.env
make -j 4
export PYTHONPATH=$PYTHONPATH:$(pwd)/python/lib
```

### RCL@DWD

```
git clone --recursive git@github.com:C2SM-RCM/extpar.git
cd extpar
git submodule update --init --recursive
mkdir build
cd build
.././configure.rcl.gcc # or ./configure.rcl.intel
gmake -j 8
```

The final step includes for all target machines to copy 
the .exe and .py files from [bin](bin) to the directory in which the namelist and all required input-data is present.

For more detailed compilation instructions see: [compile_run](doc/compile_run.md)

# Input Data

#### Data Location
In order to run Extpar, input data files for the external parameter variables are needed. The data is provided on all supported machines:
*  Tsa: _/store/c2sm/extpar_raw_data/linked_data_
*  Daint: _/store/c2sm/extpar_raw_data/linked_data_
*  Mistral: _/work/pd1167/extpar-input-data/linked_data_

The input data files are also stored in a git-LFS data repository found at: https://gitlab.dkrz.de/extpar-data/extpar-input-data.  
Instructions to download or update the input data files can be found in this repository.  
To gain access to the git-LFS input data repository, contact the Extpar source code administrator.

#### Fast Data Access
On CSCS-machines, access from comupte-nodes to _/store_ is generally slow. Instead copy all required input-data prior to your Exptar runs to $SCRATCH.
The script [extract_inputfiles_from_namelist.py](test/testsuite/bin/extract_inputfiles_from_namelist.py) is able to to extract all input-data needed from any given set of Extpar namelists (INPUT_* and namelist.py).

For Mistral no such performance penalty has been observed.

# Testing
The extpar code comes with a technical testsuite to ensure the accuracy of the results. Weekly tests run for compilers
* GCC
* Daint
* Nag

For more information about how the testsuite can be run or new test added see [testsuite-documentation](doc/testing.md)

# Information for developers
In case you want to contribute to Extpar please have a look at our [coding rules and development workflow](doc/development.md).

# Support 
In the case of issues or questions, please contact the current source code administrator (Jonas Jucker) at jonas.jucker@c2sm.ethz.ch.  


