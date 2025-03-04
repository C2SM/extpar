# Compile and Run

## Compilation

Since Version 5.4, EXTPAR is built with an autotools based build-system. This has been necessary to accomodate for the additional C source code files and newly required libraries.

There are three options to compile Extpar: 

### In-source build

```bash
./configure.<hostname>.<compiler>
source modules.env  
make   
```

### Out-of-source build

```bash
mkdir build-<my_self_defined_note>
cd build-<my_self_defined_note>
path/to/extpar/installation/configure.<hostname>.<compiler>
source modules.env  
make  
```

### Install binaries only in external directory

```bash
./configure.<hostname>.<compiler> --prefix=<my_external_directory>
source modules.env  
make install  
```

The binaries will be installed in `my_external_directory/bin`.

### Restart build from scratch

```bash
make distclean 
```

### Configure for new machines

The first step in creating a new machine setup is to take one of
the existing configure wrapper scripts and adapt it to your local
environment. The scripts are called `configure.<hostname>.<compiler>`.

## Running

The Fortran executables 

*   *extpar_aot_to_buffer.exe*
*   *extpar_landuse_to_buffer.exe*
*   *exptar_soil_to_buffer.exe*
*   *extpar_topo_to_buffer.exe*
*   *extpar_flake_to_buffer.exe* 
*   *extpar_consistency_check.exe* 

can simply be copied to the `run/` directory.
 
The main python scripts 

*   *extpar_alb_to_buffer.py*
*   *extpar_cru_to_buffer.py*
*   *extpar_ndvi_to_buffer.py*
*   *extpar_emiss_to_buffer.py*
*   *exptar_era_to_buffer.py*
*   *extpar_ahf_to_buffer.py*
*   *extpar_isa_to_buffer.py*  
*   *extpar_cdnc_to_buffer.py*
*   *extpar_edgar_to_buffer.py*  


can be treated like the Fortran binaries and copied to the `run/` directory. Make sure the *namelist.py* is also present in the `run/` directory.  

All self-written Python modules are stored in [`python/lib` :material-open-in-new:](https://github.com/C2SM/extpar/tree/master/python/lib){:target="_blank"}
and do not need to be copied to the respective run directory, but the environment variable `PYTHONPATH` needs to be set to the following:  
 
```bash
export PYTHONPATH=$PYTHONPATH:<absolute_path_to_python>/lib
```

Some runscript examples are available under [`run_scripts` :material-open-in-new:](https://github.com/C2SM/extpar/tree/master/run_scripts){:target="_blank"}.
Just adapt them to your needs!
