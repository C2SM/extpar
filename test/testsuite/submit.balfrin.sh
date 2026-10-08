#!/bin/bash

#SBATCH --job-name="extpar"
#SBATCH --nodes=1
#SBATCH --output="job.out"
#SBATCH --time=00:57:00
#SBATCH --partition=pp-short
#SBATCH --exclusive
#SBATCH --cpus-per-task=64

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export USER_ENV_ROOT=/mch-environment/v6

source ../../modules.env
source ../../.venv/bin/activate

export NETCDF_OUTPUT_FILETYPE=NETCDF4
./bin/runcontrol_functions.sh

#./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 --testlist=testlist_icon.xml --mpicmd='sleep 1 &&'  
#./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 --testlist=testlist_landuse.xml --mpicmd='sleep 1 &&'  
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 --testlist=testlist_icon.xml --mpicmd='sleep 1 &&'

