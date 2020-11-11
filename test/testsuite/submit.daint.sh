#!/bin/bash
#SBATCH --constraint=gpu
#SBATCH --output="job.out"
#SBATCH --account=g110
#SBATCH --time=02:00:00
module load daint-gpu
module load CDO
source ../../modules.env
source /project/g110/extpar_envs/venv_jenkins_daint/bin/activate


./src/testsuite.py --exe=run_extpar_cosmo.sh -v 1 -o testsuite.out --testlist=testlist_cosmo.xml --mpicmd='srun -u -n'  

# only execute dwd icon test, because of problems
# with the emissivity input data on Daint compute-nodes
#icon_tests=("dwd,icon_d2" "dwd,icon_d2_caching")
#
#for icon_test in ${icon_tests[@]}; do
#    ./src/testsuite.py -a --exe=run_extpar_icon.sh -v 1 -o testsuite.out --testlist=testlist_icon.xml --only=$icon_test --mpicmd='srun -u -n'  
#done
