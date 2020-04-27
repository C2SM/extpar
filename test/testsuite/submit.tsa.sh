#!/bin/bash
#SBATCH --job-name="extpar"
#SBATCH --nodes=1
#SBATCH --output="job.out"
#SBATCH --time=01:00:00
#SBATCH --partition=postproc

source /oprusers/osm/.opr_setup_dir
export MODULEPATH=$MODULEPATH\:$OPR_SETUP_DIR/modules/modulefiles
source ../../modules.env
module load cdo
export HDF5_DISABLE_VERSION_CHECK=1

./src/testsuite.py --exe=run_extpar_cosmo.sh -v 1 -o testsuite.out --testlist=testlist_cosmo.xml --mpicmd='srun -u -n'
