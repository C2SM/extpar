#!/bin/bash
#SBATCH --job-name="extpar"
#SBATCH --nodes=1
#SBATCH --output="job.out"
#SBATCH --time=00:20:00
#SBATCH --partition=compute
#SBATCH --account=mh0287

module load python
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 3 -o testsuite.out --testlist=testlist_corine.xml --mpicmd='srun -u -n'  
