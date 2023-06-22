#!/bin/bash
#SBATCH --constraint=gpu
#SBATCH --output="job.out"
#SBATCH --account=g110
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=12
module load daint-gpu
module load CDO
source ../../modules.env
source /project/g110/extpar/venv_daint/bin/activate

# this is needed to not have problems with the emissivity input files
export HDF5_USE_FILE_LOCKING=FALSE

# avoid messages for CDO calls when submitted from login node
export PMI_NO_PREINITIALIZE=1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


./src/testsuite.py --exe=run_extpar_cosmo.sh -v 2 -o testsuite.out --testlist=testlist_cosmo.xml --mpicmd='srun -u -n'  
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_icon.xml --mpicmd='srun -u -n'  
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_landuse.xml --mpicmd='srun -u -n'
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_art.xml --mpicmd='srun -u -n'  
./src/testsuite.py -a --exe=run_extpar_wrapper.sh -v 2 -o testsuite.out --testlist=testlist_tools.xml --mpicmd='srun -u -n'  
