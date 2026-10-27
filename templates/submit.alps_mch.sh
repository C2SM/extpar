#!/bin/bash
#SBATCH --job-name="extpar"
#SBATCH --nodes=1
#SBATCH --output="job.out"
#SBATCH --time=02:00:00
#SBATCH --partition=postproc
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --cpus-per-task=64
#SBATCH --account=@ACCOUNT@

export USER_ENV_ROOT=/mch-environment/v6
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NETCDF_OUTPUT_FILETYPE=NETCDF4
export PYTHONPATH=@PYTHONPATH@:$PYTHONPATH

ulimit -s unlimited
ulimit -c unlimited

logfile=extpar.log
rm -f $logfile

if [ -e modules.env ]; then source modules.env; fi
if [ -n "${VIRTUAL_ENV:-}" ]; then source "${VIRTUAL_ENV}/bin/activate"; fi
source runcontrol_functions.sh

for exe in @EXTPAR_EXECUTABLES@
do
    run_sequential $exe
done
