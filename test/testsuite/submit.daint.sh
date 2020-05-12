#!/bin/bash
#SBATCH --constraint=gpu
#SBATCH --output="job.out"
#SBATCH --account=g110
#SBATCH --time=01:00:00
module load daint-gpu
module load CDO
source ../../modules.env
module load cray-python
module load netcdf-python

./src/testsuite.py --exe=run_extpar_cosmo.sh -v 1 -o testsuite.out --testlist=testlist_cosmo.xml --mpicmd='srun -u -n'  
