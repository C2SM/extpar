# Home

## General Information
EXTPAR (**Ext**ernal **Par**ameters for Numerical Weather Prediction and Climate Application) is an official software of the [COSMO Consortium](http://www.cosmo-model.org/content/default.htm).  It is used to prepare the external parameter data files that are used as input for the COSMO and the ICON model.

The code is written in Fortran 90 and in Python. The Python scripts use CDO for the most compute-intensive parts. The code is also accelerated in some places with OpenMP parallelization.

The code once compiled generates 6 Fortran executables and 9 Python scripts, which can be run simultaneously except for the final extpar_consistency_check.exe, which is used to tie together all the external parameter results into one output file.

Information about the latest changes can be found in the [Release Notes](https://github.com/C2SM-RCM/extpar/releases).

A full documentation of the code can be found as an [assets of each release](https://github.com/C2SM-RCM/extpar/releases).

## Quick Start

### Container
The easiest way to use extpar is through the container provided with [Dockerfile](Dockerfile). 
A ready-to-use image can be downloaded from [C2SM docker hub](https://hub.docker.com/repository/docker/c2sm/extpar/general) 
or even simpler via CLI `docker pull c2sm/extpar:tagname`.

Alternatively an image is provided as an [asset of each release](https://github.com/C2SM-RCM/extpar/releases)

#### WrapExtpar
The image provides a wrapper that only requires to set basic options, all other details are handled by the wrapper.

The wrapper needs two different kinds of input:


_1. Extpar settings as JSON, see official docs_

```json
{
  "extpar": {
    "igrid_type": 1,
    "iaot_type": 1,
    "ilu_type": 1,
    "ialb_type": 1,
    "isoil_type": 1,
    "itopo_type": 1,
    "lsgls": false,
    "lfilter_oro": false,
    "lurban": false
  }
}
  ```

_2. Execution options_
```console
  --input-grid INPUT_GRID
                        COSMO: Fortran Namelist "INPUT_COSMO_GRID", ICON: Icon
                        grid file
  --raw-data-path RAW_DATA_PATH
                        Path to folder "linked_data" of exptar-input-data
                        repository
  --run-dir RUN_DIR     Folder for running Extpar
  --account ACCOUNT     Account for slurm job
  --host HOST           Host
  --no-batch-job        Run jobscript not as batch job
  ```

An example call could look like
```bash
docker run -v /c2sm-data/extpar-input-data:/data \
           -v /icon-grids:/grid \
           -v /my_local_dir:/work \
           extpar \ 
           python3 -m extpar.WrapExtpar \
           --run-dir /work \
           --raw-data-path /data/linked_data \
           --account none \
           --no-batch-job \
           --host docker \
           --input-grid /grid/icon_grid.nc \
           --extpar-config /work/config.json
```
Below is a more detailed explanation about the mounted volumes:

* `-v /c2sm-data/extpar-input-data:/data`: Mounts the input data at `/data` inside the container. This should be aligned with the `--raw-data-path` argument.
* `-v /icon-grids:/grid`: Mounts a local folder with icon grids under `/grid` inside the container. This should be aligned with the `--input-grid` argument.
* `-v /my_local_dir:/work`: Mounts a local folder for extpar output at `/work` inside the container. This should be aligned with the `--run-dir` argument.

### Individual executables
For those who require a more custom setup of Extpar or need settings that are not possible to specify through the wrapper, you can run each executable within the image too. For example:

```bash
docker run extpar bash -c "extpar_topo_to_buffer"

## Bare metal build on Levante

The installation steps are

```bash
git clone --recursive git@github.com:C2SM-RCM/extpar.git
cd extpar
git submodule update
./configure.levante.gcc
source modules.env
make -j 4
```

Furthermore copy all the .exe and .py files from [bin](bin) to the directory 
in which the namelist and all required input-data is present.

You do then have two choices to run extpar:

1. configure the `PYTHONPATH` variable such that it includes to the `python/lib`
   folder of the source repository
2. build and install a python package for your user account

#### Installing extpar

After you prepared extpar (see above), you have to options to install and run
the software.

##### Option 1: PYTHONPATH

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

##### Option 2: Build and install a python package

Alternatively you can build a python package and install it to your libraries.
This has the advantages that the executables can be ran from anywhere in the
system without the need to copy the executables themselves.

To build the package, now run

```bash
python setup.py sdist
```

You can then install it via

```bash
pip install dist/extpar-*.tar.gz
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


## Input Data

### Data Location
In order to run Extpar, input data files for the external parameter variables are needed. The data is provided on all supported machines:
*  Levante: _/work/pd1167/extpar-input-data/linked_data_
*  CO2 (ETHZ): _/c2sm-data/extpar-input-data_

The input data files are also stored in a git-LFS data repository found at: https://gitlab.dkrz.de/extpar-data/extpar-input-data.
Instructions to download or update the input data files can be found in this repository.
To gain access to the git-LFS input data repository, contact the Extpar source code administrator.

## Testing
The extpar code comes with a technical testsuite to ensure the accuracy of the results. Weekly tests run for compilers
* GCC

For more information about how the testsuite can be run or new test added see [testsuite-documentation](doc/testing.md)

## Information for developers
In case you want to contribute to Extpar please have a look at our [coding rules and development workflow](doc/development.md).

## Support
In the case of issues or questions, please contact the current source code administrator (Jonas Jucker) at jonas.jucker@c2sm.ethz.ch.
