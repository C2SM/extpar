#!/bin/bash
#SBATCH --job-name="extpar"
#SBATCH --mem-per-cpu=10G
#SBATCH --ntasks=12
#SBATCH --output="job.out"
#SBATCH --time=02:00:00
#SBATCH --partition=prepost,shared,gpu
#SBATCH --account=mh0287

export OMP_NUM_THREADS=12

source ../../modules.env
./src/testsuite.py --exe=run_extpar_cosmo.sh -v 1 -o testsuite.out --testlist=testlist_cosmo_intel.xml --mpicmd='srun -u -n' 
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 1 -o testsuite.out --testlist=testlist_icon_intel.xml --mpicmd='srun -u -n'  
