#!/bin/bash -l
#
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --partition=xfer

module load cray-python
source /project/g110/extpar_envs/venv_jenkins_daint/bin/activate

python parse_namelist.py

command="cp"
echo -e "$SLURM_JOB_NAME started on $(date):\n $command $1 $2\n"
for file in $(cat files_to_copy.txt); do
    echo $file
    srun -n $SLURM_NTASKS $command $file input-data/.
done
echo -e "$SLURM_JOB_NAME finished on $(date)\n"


