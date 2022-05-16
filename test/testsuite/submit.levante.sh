#!/bin/bash
#SBATCH --job-name="extpar"
#SBATCH --mem=36G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12 
#SBATCH --output="job.out"
#SBATCH --time=02:00:00
#SBATCH --partition=gpu,shared,compute
#SBATCH --account=mh0287

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

source ../../modules.env

./src/testsuite.py --exe=run_extpar_cosmo.sh -v 2 -o testsuite.out --testlist=testlist_cosmo.xml --mpicmd='srun -u -n' 
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_icon.xml --mpicmd='srun -u -n'  
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_landuse.xml --mpicmd='srun -u -n'  
