#!/bin/bash
#SBATCH --output="job.out"
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=128
source ../../modules.env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_icon.xml --mpicmd='srun -u -n'  
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_landuse.xml --mpicmd='srun -u -n'
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_art.xml --mpicmd='srun -u -n'  
