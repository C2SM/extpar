#!/bin/bash
#SBATCH --constraint=gpu
#SBATCH --account=c15
#SBATCH --time=01:00:00
module load daint-gpu
module load CDO
./src/testsuite.py --exe=run_extpar_mch.sh -v 1 -o testsuite.out --testlist=testlist_cosmo.xml --mpicmd='srun -u -n'  
