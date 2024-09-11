#!/bin/bash

# import functions to launch Extpar executables
. ./runcontrol_functions.sh

ulimit -s unlimited
ulimit -c 0

# get hostname
hostname="`echo $HOSTNAME`"
logfile="runscript.log"

rm ${logfile}

#--------------------------------------------------------------------------------
# define host-dependent paths and variables

# Levante
if [[ $hostname == l* ]]; then

    # directories
    data_dir=/work/pd1167/extpar-input-data/linked_data
    wrapper_host='levante'
    account=mh0287

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

cd $rootdir/python
./WrapExtpar.py --input_cosmo_grid ../test/testsuite/data/mch/c1_aster/INPUT_COSMO_GRID   --iaot_type 1 --ilu_type 1 --ialb_type 1 --isoil_type 1 --itopo_type 1 --raw_data_path $data_dir --run_dir $currentdir/sandbox_wrapper --account $account --host $wrapper_host --no_batch_job >> ${logfile}

cd - >> ${logfile}
