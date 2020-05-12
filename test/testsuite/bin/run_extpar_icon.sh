#!/bin/ksh

# import functions to launch Extpar executables
. ./runcontrol_functions.sh

ulimit -s unlimited
ulimit -c 0

# get hostname
hostname="`echo $HOSTNAME`"
logfile="extpar_runscript.log"

rm ${logfile}

#--------------------------------------------------------------------------------
# define host-dependent paths and variables

# CSCS
if [[ $hostname == daint* || $hostname == nid* ]]; then

    # NetCDF raw data for external parameter
    data_dir=/store/c2sm/extpar_raw_data/linked_data

# mistral
elif [[ $hostname == m* ]]; then

    export OMP_NUM_THREADS=8
    
    # directories
    data_dir=/scratch/b/b381001/extpar-input-data/linked_data

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
binary_emiss=extpar_emiss_to_buffer.py
binary_tclim=extpar_cru_to_buffer.py

# fortran executables
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_isa=extpar_isa_to_buffer.exe
binary_consistency_check=extpar_consistency_check.exe

# link raw data files to local workdir
ln -s -f ${data_dir}/*.nc .
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define test-specific paths and variables 

type_of_test=`echo $currentdir | rev | cut -d"/" -f2 | rev`
name_of_test=`echo $currentdir | rev | cut -d"/" -f1 | rev`
icon_grid_dir=$rootdir/test/testsuite/data/$type_of_test/$name_of_test/
icon_grid_file=icon_grid*

# mpim
if [[ $type_of_test == mpim ]]; then

    ln -sf ${icon_grid_dir}/ei_sst_an1986-2015_0013_R02B04_G_BUFFER.nc .
    ln -sf ${icon_grid_dir}/ei_t2m_an1986-2015_0013_R02B04_G_BUFFER.nc .

    # jj_tmp: hack to replace missing Fortran namelist to launch emiss_to_buffer.py
    touch INPUT_EMISS

# dwd
elif [[ $type_of_test == dwd ]]; then

    ln -sf ${icon_grid_dir}/ei_2t_an1986-2015_domain2_DOM01_BUFFER.nc .
    ln -sf ${icon_grid_dir}/ei_an1986-2015_domain2_DOM01_BUFFER.nc .

# ecmwf
elif [[ $type_of_test == ecmwf ]]; then

    ln -sf ${icon_grid_dir}/ei_2t_an1986-2015_0099_R19B10_BUFFER.nc
    ln -sf ${icon_grid_dir}/ei_an1986-2015_0099_R19B10_BUFFER.nc

#unknown test

else

    # exit script in case of unknown host
    echo ERROR: Unkown test: $type_of_test >> ${logfile}
    exit 1
fi

ln -sf ${icon_grid_dir}/${icon_grid_file} .

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# launch extpar executables

echo ">>>> Data will be processed and produced in `pwd` <<<<"

# 1) topography needs to be processed first - result is input fro the
#    CRU data processing

run_sequential ${binary_topo}

#________________________________________________________________________________
# 2) all other executables

run_sequential ${binary_alb}

run_sequential ${binary_ndvi}

run_sequential ${binary_tclim}

run_sequential ${binary_aot}

run_sequential ${binary_lu}

run_sequential ${binary_soil}

run_sequential ${binary_flake}

if [ -f INPUT_EMISS ] ; then
    run_sequential ${binary_emiss}
fi

run_sequential ${binary_consistency_check}
#________________________________________________________________________________

echo ">>>> External parameters for ICON model generated <<<<"
