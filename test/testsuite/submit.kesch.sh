#!/bin/bash
#SBATCH --job-name="extpar"
#SBATCH --nodes=1
#SBATCH --output="job.out"
#SBATCH --time=01:00:00
#SBATCH --partition=postproc

export MODULEPATH=$MODULEPATH:/oprusers/owm/modules/RH7.5/modulefiles
module load PE/17.06
module load gcc
module load netcdf/4.4.1.1-gmvolf-17.02
module load cdo
export HDF5_DISABLE_VERSION_CHECK=1
./src/testsuite.py --exe=run_extpar_kesch.sh -v 1 -o testsuite.out --testlist=testlist_cosmo.xml --mpicmd='srun -u -n'  
./src/testsuite.py --exe=run_extpar_kesch_dwd.sh -a -v 1 -o testsuite.out --testlist=testlist_cosmo_dwd.xml --mpicmd='srun -u -n'  

