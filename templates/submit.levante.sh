#!/bin/bash
#SBATCH --job-name="extpar"
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12 
#SBATCH --output="job.out"
#SBATCH --time=02:00:00
#SBATCH --partition=shared,compute
#SBATCH --account=@ACCOUNT@

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONPATH=@PYTHONPATH@:$PYTHONPATH

logfile=extpar.log
rm $logfile

if [ -e modules.env ]; then source modules.env; fi
source runcontrol_functions.sh

for exe in @EXTPAR_EXECUTABLES@
do
    run_sequential $exe
done

