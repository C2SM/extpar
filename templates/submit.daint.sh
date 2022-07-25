#!/bin/bash
#SBATCH --constraint=gpu
#SBATCH --output="job.out"
#SBATCH --account=@ACCOUNT@
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=12
module load daint-gpu
module load CDO
source modules.env
source /project/g110/extpar/venv_daint/bin/activate

export PYTHONPATH=@PYTHONPATH@:$PYTHONPATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

logfile=extpar.log
rm $logfile
source runcontrol_functions.sh

# this is needed to not have problems with the emissivity input files
export HDF5_USE_FILE_LOCKING=FALSE

# avoid messages for CDO calls when submitted from login node
export PMI_NO_PREINITIALIZE=1


for exe in @EXTPAR_EXECUTABLES@
do
    run_sequential $exe
done

