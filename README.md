# EXTPAR

# General Information
EXTPAR (External Parameters for Numerical Weather Prediction and Climate Application) is an official software of the [COSMO Consortium](www.cosmo-model.org).  It is used to prepare the external parameter data files that are used as input for the COSMO and the ICON model.

The code is written in Fortran 90 and in Python. The Python scripts use CDO for the most compute-intensive parts.  Currently the code is tested regularly using the gcc, NAG, and Intel Fortran compilers.  The code is also accelerated in some places with OpenMP parallelization.

The code once compiled generates 6 Fortran executables and 7 Python scripts, which can be run simultaneously except for the final extpar_consistency_check.exe, which is used to tie together all the external parameter results into one output file.


Information about the latest changes can be found in the [Release Notes](ReleaseNotes.md).

A full documentation of the code can be found as an [assets of each release](https://github.com/C2SM-RCM/extpar/releases).

# Quick Start


We support the setup for extpar on three different HPC-infrastructures.
*Daint*, *Tsa* and *Levante*. The installation steps are

1. clone the source repository
2. run the configuration script for the corresponding HPC-infrastructure
3. copy all the .exe and .py files from [bin](bin) to the directory in which
   the namelist and all required input-data is present.

You do then have two choices to run extpar:

1. configure the `PYTHONPATH` variable such that it includes to the `python/lib`
   folder of the source repository
2. build and install a python package for your user account

## Preparing extpar

### Daint

```bash
git clone --recursive git@github.com:C2SM-RCM/extpar.git
cd extpar
git submodule update
./configure.daint.gcc
source modules.env
make -j 4
module load daint-gpu
module load CDO
```

### Tsa

```bash
git clone --recursive git@github.com:C2SM-RCM/extpar.git
cd extpar
git submodule update
./configure.tsa.gcc
source modules.env
make -j 4
source /oprusers/osm/.opr_setup_dir
export MODULEPATH=$MODULEPATH\:$OPR_SETUP_DIR/modules/modulefiles
source /project/g110/extpar/venv_tsa/bin/activate
```

### Levante

```bash
git clone --recursive git@github.com:C2SM-RCM/extpar.git
cd extpar
git submodule update
./configure.levante.gcc
source modules.env
make -j 4
```

## Installing extpar

After you prepared extpar (see above), you have to options to install and run
the software.

### Option 1: PYTHONPATH

If you like to run the extpar scripts without installing a package, make sure
to have the `python/lib` folder in your `PYTHONPATH` variable. You can do this
via

```bash
export PYTHONPATH=$PYTHONPATH:$(pwd)/python/lib
```

Afterwards you can `cd` into the [`bin/`](bin) directory and run the
corresponding executables, e.g.

```bash
cd bin
./extpar_aot_to_buffer.exe
```

For more detailed compilation instructions see: [compile_run](doc/compile_run.md)

### Option 2: Build and install a python package

Alternatively you can build a python package and install it to your libraries.
This has the advantages that the executables can be ran from anywhere in the
system without the need to copy the executables themselves.

Before building the package, we need to resolve the symbolic links that
have been added during the `make` command and turn them into regular files:

```bash
for f in $(find bin -type l);do cp --remove-destination $(readlink $f) $f;done
```

To build the package, now run

```bash
python setup.py bdist_wheel
```

You can then install it via

```bash
pip install dist/extpar_*.whl
```

---
**Note:** If you do not have the permissions to install it into the system-wide python
library, it will be installed for your user account only (you can also add the
`--user` flag to `pip` to force this behaviour).

If you did not install `extpar` into the system-libraries, make sure
that the `bin` folder of your local user is on your `PATH` variable to be able
to run the extpar scripts. This is usually done via

```bash
export PATH="$HOME/.local/bin:$PATH"
```

You can add this line to your `.bashrc` file if you want.

---

You can then call the functionalities of `WrapExtpar.py` via

```bash
python -m extpar.WrapExtpar
```

or import the script in python via

```python
from extpar.WrapExtpar import generate_external_parameters
```

Or you call the executable scripts in your run directory, e.g.

```bash
extpar_aot_to_buffer.exe
```


# Input Data

#### Data Location
In order to run Extpar, input data files for the external parameter variables are needed. The data is provided on all supported machines:
*  Tsa: _/store/c2sm/extpar_raw_data/linked_data_
*  Daint: _/store/c2sm/extpar_raw_data/linked_data_
*  Levante: _/work/pd1167/extpar-input-data/linked_data_

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

For more information about how the testsuite can be run or new test added see [testsuite-documentation](doc/testing.md)

# Information for developers
In case you want to contribute to Extpar please have a look at our [coding rules and development workflow](doc/development.md).

# Support
In the case of issues or questions, please contact the current source code administrator (Jonas Jucker) at jonas.jucker@c2sm.ethz.ch.


