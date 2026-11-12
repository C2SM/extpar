#!/bin/bash

#SBATCH --job-name="extpar"
#SBATCH --nodes=1
#SBATCH --output="job.out"
#SBATCH --time=01:00:00
#SBATCH --partition=pp-short
#SBATCH --exclusive
#SBATCH --cpus-per-task=64

export USER_ENV_ROOT=/mch-environment/v8
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NETCDF_OUTPUT_FILETYPE=NETCDF4

source ../../modules.env
source ../../.venv/bin/activate

./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 --testlist=testlist_icon.xml --mpicmd='sleep 1 &&'
