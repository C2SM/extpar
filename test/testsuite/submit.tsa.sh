#!/bin/bash
#SBATCH --job-name="extpar"
#SBATCH --nodes=1
#SBATCH --output="job.out"
#SBATCH --time=01:30:00
#SBATCH --partition=postproc
#SBATCH --cpus-per-task=12

source /oprusers/osm/.opr_setup_dir
export MODULEPATH=$MODULEPATH\:$OPR_SETUP_DIR/modules/modulefiles
source /project/g110/extpar_envs/venv_jenkins_tsa/bin/activate
source ../../modules.env
module load cdo
export HDF5_DISABLE_VERSION_CHECK=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./src/testsuite.py  --exe=run_extpar_cosmo.sh -v 2 -o testsuite.out --testlist=testlist_cosmo.xml --mpicmd='srun -u -n'
./src/testsuite.py  -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_icon.xml --mpicmd='srun -u -n'
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_landuse.xml --mpicmd='srun -u -n'
