#!/bin/bash
#SBATCH --job-name="extpar"
#SBATCH --nodes=1
#SBATCH --output="job.out"
#SBATCH --time=03:00:00
#SBATCH --partition=compute
#SBATCH --account=mh0287

module load cdo
export OMP_NUM_THREADS=1
./src/testsuite.py --exe=run_extpar_mistral.sh -v 1 -o testsuite.out --testlist=testlist_cosmo.xml --mpicmd='srun -u -n'  
