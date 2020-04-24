#!/bin/bash
#SBATCH --job-name="extpar"
#SBATCH --nodes=1
#SBATCH --output="job.out"
#SBATCH --time=01:00:00
#SBATCH --partition=postproc

source /oprusers/osm/.opr_setup_dir
export MODULEPATH=$MODULEPATH\:$OPR_SETUP_DIR/modules/modulefiles
module load gnu_PE/17.02 
module load PrgEnv-gnu 
module load netcdf/4.2.1.1-gnu-5.4.0
source /users/juckerj/venv_jenkins/bin/activate
module load cdo
export HDF5_DISABLE_VERSION_CHECK=1

./src/testsuite.py --exe=run_extpar_cosmo.sh -v 1 -o testsuite.out --testlist=testlist_cosmo.xml --mpicmd='srun -u -n'  

