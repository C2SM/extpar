# Home

## General Information

EXTPAR (**Ext**ernal **Par**ameters for Numerical Weather Prediction and Climate Application) is an official software of the [COSMO Consortium :material-open-in-new:](http://www.cosmo-model.org/content/default.htm){:target="_blank"}.  It is used to prepare the external parameter data files that are used as input for the COSMO and the ICON model.

The code is written in Fortran 90 and in Python. The Python scripts use [CDO :material-open-in-new:](https://code.mpimet.mpg.de/projects/cdo){:target="_blank"} (Climate Data Operators) 
for the most compute-intensive parts. The code is also accelerated in some places with OpenMP parallelization.

Once compiled, the code generates 6 Fortran executables and 9 Python scripts, which can be run simultaneously except for the final extpar_consistency_check.exe, which is used to tie together all the external parameter results into one output file.

Information about the latest changes can be found in the [Release Notes on GitHub :material-open-in-new:](https://github.com/C2SM/extpar/releases){:target="_blank"}.

The technical and scientific documentation can be found in the [User and Implementation Guide](user_manual/index.md).

## Quick Start

### Container

The easiest way to use EXTPAR is through the container provided with [Dockerfile :material-open-in-new:](https://github.com/C2SM/extpar/blob/master/Dockerfile){:target="_blank"}. 
A ready-to-use image can be downloaded from [C2SM docker hub :material-open-in-new:](https://hub.docker.com/repository/docker/c2sm/extpar/general){:target="_blank"} 
or even simpler via CLI:

```shell 
docker pull c2sm/extpar:tagname
```

Alternatively, an image is provided as an [asset of each release :material-open-in-new:](https://github.com/C2SM/extpar/releases){:target="_blank"}

#### WrapExtpar

The image provides a wrapper that only requires to set basic options, all other details are handled by the wrapper.

The wrapper needs two different kinds of input:

_1. EXTPAR settings as JSON, see official docs_

```json
{
  "extpar": {
    "igrid_type": 1,
    "iaot_type": 1,
    "ilu_type": 1,
    "ialb_type": 1,
    "isoil_type": 1,
    "itopo_type": 1,
    "it_cl_type": 1,
    "iera_type": 1,
    "iemiss_type": 1,
    "enable_cdnc": false,
    "enable_edgar": false,
    "enable_art": false,
    "use_array_cache": false,
    "lsgsl": false,
    "lfilter_oro": false,
    "lurban": false,
    "lradtopo": true,
    "radtopo_radius": 40000
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
  --run-dir RUN_DIR     Folder for running EXTPAR
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
* `-v /my_local_dir:/work`: Mounts a local folder for EXTPAR output at `/work` inside the container. This should be aligned with the `--run-dir` argument.

### Individual Executables

For those who require a more custom setup of EXTPAR or need settings that are not possible to specify through the wrapper, you can run each executable within the image too. For example:

```bash
docker run extpar bash -c "extpar_topo_to_buffer"
```

### Bare Metal Build on Levante

The installation steps are

```bash
git clone --recursive git@github.com:C2SM/extpar.git
cd extpar
git submodule update
./configure.levante.gcc
source modules.env
make -j 4
```

Furthermore copy all the `.exe` and `.py` files from `bin` to the directory 
in which the namelist and all required input-data is present.

You do then have two choices to run EXTPAR:

1. configure the `PYTHONPATH` variable such that it includes to the `python/lib`
   folder of the source repository
2. build and install a python package for your user account

#### Installing EXTPAR

After you prepared EXTPAR (see above), you have two options to install and run
the software.

##### Option 1: PYTHONPATH

If you like to run the EXTPAR scripts without installing a package, make sure
to have the `python/lib` folder in your `PYTHONPATH` variable. You can do this
via

```bash
export PYTHONPATH=$PYTHONPATH:$(pwd)/python/lib
```

Afterwards you can `cd` into the `bin/` directory and run the
corresponding executables, e.g.

```bash
cd bin
./extpar_aot_to_buffer.exe
```

For more detailed compilation instructions see the [Compile and Run](compile_run.md) section.

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

!!! note
    If you do not have the permissions to install it into the system-wide python
    library, it will be installed for your user account only (you can also add the
    `--user` flag to `pip` to force this behaviour).

    If you did not install `extpar` into the system libraries, make sure
    that the `bin` folder of your local user is on your `PATH` variable to be able
    to run the EXTPAR scripts. This is usually done via

    ```bash
    export PATH="$HOME/.local/bin:$PATH"
    ```

You can then call the functionalities of `WrapExtpar.py` via

```bash
python -m extpar.WrapExtpar
```

or import the script in Python via

```python
from extpar.WrapExtpar import generate_external_parameters
```

Or you call the executable scripts in your run directory, e.g.

```bash
extpar_aot_to_buffer.exe
```

## Input Data

### Data Location

In order to run EXTPAR, input data files for the external parameter variables are needed. The data is provided on all supported machines:
 
=== "Levante (DKRZ)"
    ```console
    /work/pd1167/extpar-input-data/linked_data
    ```

=== "co2 (ETHZ)"
    ```console
    /c2sm-data/extpar-input-data
    ```

The input data files are also stored in a git-LFS data repository found at: [https://gitlab.dkrz.de/extpar-data/extpar-input-data :material-open-in-new:](https://gitlab.dkrz.de/extpar-data/extpar-input-data){:target="_blank"}.
Instructions to download or update the input data files can be found in this repository.
To gain access to the git-LFS input data repository, contact the EXTPAR source code administrator.

## Testing

The EXTPAR code comes with a technical testsuite to ensure the accuracy of the results. Weekly tests run for compilers:

* GCC

For more information about how the testsuite can be run or new test added see the [Testing](testing.md) section.

## Information for Developers

In case you want to contribute to EXTPAR please have a look at our [coding rules and development workflow](development.md).

## Support

In the case of issues or questions, please create an [issue on GitHub :material-open-in-new:](https://github.com/C2SM/extpar/issues){:target="_blank"}.
