#!/bin/ksh
      
# import functions to launch Extpar executables
. ./runcontrol_functions.sh

ulimit -s unlimited
ulimit -c unlimited

# get hostname
hostname="`echo $HOSTNAME`"
logfile="extpar_runscript.log"

rm ${logfile}

#--------------------------------------------------------------------------------
# define host-dependent paths and variables

# CSCS
if [[ $hostname == kesch* || $hostname == daint* || $hostname == tsa* || $hostname == arolla* || $hostname == nid* ]]; then

    # NetCDF raw data for external parameter
    data_dir=/store/c2sm/extpar_raw_data/linked_data

# mistral
elif [[ $hostname == m* ]]; then

    # NetCDF raw data for external parameter
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
echo PYTHONPATH: ${PYTHONPATH} >> ${logfile}
cd $src_python
export PYTHONPATH=$PYTHONPATH:$(pwd)
cd - >> ${logfile}

echo PYTHONPATH: ${PYTHONPATH} >> ${logfile}

# Names of executables

# python executables
binary_alb=extpar_alb_to_buffer.py
binary_ndvi=extpar_ndvi_to_buffer.py
binary_tclim=extpar_cru_to_buffer.py

# fortran executables
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_ahf=extpar_ahf_to_buffer.exe
binary_isa=extpar_isa_to_buffer.exe
binary_consistency_check=extpar_consistency_check.exe

# link raw data files to local workdir
ln -s -f ${data_dir}/*.nc .
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define test-specific paths and variables 

type_of_test=`echo $currentdir | rev | cut -d"/" -f2 | rev`
name_of_test=`echo $currentdir | rev | cut -d"/" -f1 | rev`

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# launch extpar executables

echo ">>>> Data will be processed and produced in `pwd` <<<<"

if [[ $type_of_test == clm ]]; then

    # remove S_ORO fields
    rm S_ORO_*
fi

# dwd
if [[ $type_of_test == dwd ]]; then

    # binary_tclim needs output of binary_topo
    run_sequential ${binary_topo}
else
    run_parallel ${binary_topo} 
fi

run_parallel ${binary_alb}
run_parallel ${binary_aot}
run_parallel ${binary_tclim}
run_parallel ${binary_lu}
run_parallel ${binary_ndvi} 
run_parallel ${binary_soil} 
run_parallel ${binary_flake}

if [ -f INPUT_AHF ] ; then
    run_parallel ${binary_ahf}
fi
if [ -f INPUT_ISA ] ; then
    run_parallel ${binary_isa}
fi

#--------------------------------------------------------------------------------
# IMPORTANT WAIT FOR ALL PARALLEL EXECUTABLES TO END
wait
#--------------------------------------------------------------------------------

# count non-zero exit status
error_count=0

if [[ $type_of_test != dwd ]]; then
    check_exit_status ${binary_topo}  error_count
fi

check_exit_status ${binary_alb} error_count
check_exit_status ${binary_aot} error_count
check_exit_status ${binary_tclim} error_count
check_exit_status ${binary_lu} error_count
check_exit_status ${binary_ndvi}  error_count
check_exit_status ${binary_soil}  error_count
check_exit_status ${binary_flake} error_count

if [ -f INPUT_AHF ] ; then
    check_exit_status ${binary_ahf} error_count
fi
if [ -f INPUT_ISA ] ; then
    check_exit_status ${binary_isa} error_count
fi

# if execution of some Extpar executables failed exit script
if [[ $error_count > 0 ]]; then

    echo "*****************************************"
    echo ""
    echo "Some Extpar executables did not terminate correctly!"
    echo "See ${logfile} for more information"
    echo ""
    echo "*****************************************"
    exit 1 

fi

run_sequential ${binary_consistency_check}

#--------------------------------------------------------------------------------
# clean-up
rm exit_status_*
rm time_*

echo ">>>> External parameters for COSMO model generated <<<<"
