#!/bin/bash
      
# import functions to launch Extpar executables
. ./runcontrol_functions.sh

ulimit -s unlimited
ulimit -c unlimited

# get hostname
hostname="`echo $HOSTNAME`"
logfile="runscript.log"

rm ${logfile}

#--------------------------------------------------------------------------------
# define host-dependent paths and variables

# Daint
if [[ $hostname == daint* || $hostname == nid* ]]; then

    data_dir="$PWD/../../../input-data"

# Tsa
elif [[ $hostname == tsa* || $hostname == arolla* ]]; then

    # NetCDF raw data for external parameter
    data_dir="$PWD/../../../input-data"

# mistral
elif [[ $hostname == m* ]]; then

    # NetCDF raw data for external parameter
    data_dir=/work/pd1167/extpar-input-data/linked_data

# unkown host
else

    # exit script in case of unknown host
    echo ERROR: Unkown host: $hostname >> ${logfile}
    exit 1
fi

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define paths and variables independent from host

# directories
currentdir=$(pwd)
rootdir=${currentdir}/../../../../..
src_python=${rootdir}/python/lib

# change dir to src_python to get absolute path
cd $src_python
export PYTHONPATH=$PYTHONPATH:$(pwd)
cd - >> ${logfile}

echo PYTHONPATH: ${PYTHONPATH} >> ${logfile}

# Names of executables

# python executables
binary_alb=extpar_alb_to_buffer.py
binary_ndvi=extpar_ndvi_to_buffer.py
binary_tclim=extpar_cru_to_buffer.py
binary_ahf=extpar_ahf_to_buffer.py
binary_isa=extpar_isa_to_buffer.py

# fortran executables
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_consistency_check=extpar_consistency_check.exe

# link raw data files to local workdir
ln -s -f ${data_dir}/*.nc .
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define test-specific paths and variables 

type_of_test=`echo $currentdir | rev | cut -d"/" -f2 | rev`
name_of_test=`echo $currentdir | rev | cut -d"/" -f1 | rev`

echo Current test is $type_of_test/$name_of_test  >> ${logfile}
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# launch extpar executables

echo ">>>> Data will be processed and produced in `pwd` <<<<"

if [[ $type_of_test == clm ]]; then

    # remove S_ORO fields
    rm S_ORO_*
fi

run_sequential ${binary_topo}
run_sequential ${binary_alb}
run_sequential ${binary_aot}
run_sequential ${binary_tclim}
run_sequential ${binary_lu}
run_sequential ${binary_ndvi} 
run_sequential ${binary_soil} 
run_sequential ${binary_flake}

if [[ $name_of_test == c7_globe ]]; then
    run_sequential ${binary_ahf}
    run_sequential ${binary_isa}
fi

run_sequential ${binary_consistency_check}

echo ">>>> External parameters for COSMO model generated <<<<"
