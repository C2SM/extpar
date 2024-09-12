#!/bin/bash
#SBATCH --output="job.out"
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=128
source ../../modules.env

python -m venv .venv
source .venv/bin/activate
pip install -r ../../requirements.txt

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./src/testsuite.py --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_icon.xml --mpicmd='srun -u -n'  
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_landuse.xml --mpicmd='srun -u -n'
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_art.xml --mpicmd='srun -u -n'  
