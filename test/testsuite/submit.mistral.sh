#!/bin/bash
#SBATCH --job-name="extpar"
#SBATCH --mem-per-cpu=10G
#SBATCH --ntasks=1
#SBATCH --output="job.out"
#SBATCH --time=03:00:00
#SBATCH --partition=prepost,shared,gpu
#SBATCH --account=mh0287

source ../../modules.env
export OMP_NUM_THREADS=2 

./src/testsuite.py --exe=run_extpar_cosmo.sh -v 1 -o testsuite.out --testlist=testlist_cosmo.xml --mpicmd='srun -u -n' 
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 1 -o testsuite.out --testlist=testlist_icon.xml --mpicmd='srun -u -n'  
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 1 -o testsuite.out --testlist=testlist_landuse.xml --mpicmd='srun -u -n'  
