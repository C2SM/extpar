                                         #### BATCH_SYSTEM=CRAY ####
########################################################################
#PBS -S /bin/bash
#PBS -q np
#PBS -N diurnal
#PBS -o end_LL_det
#PBS -j oe
#PBS -v STHOST=sc2
#PBS -m n
#PBS -r n
#PBS -l EC_nodes=5
#PBS -l EC_threads_per_task=4
#PBS -l EC_hyperthreads=2
#PBS -l EC_tasks_per_node=18
#PBS -l EC_total_tasks=60
#PBS -l walltime=05:00:10
########################################################################



module load gcc/6.3.0
#ICON_GRID="0026_R03B07_G"
ICON_GRID="0044_R19B07" # 0099_R19B08
ICON_GRID="0099_R19B08"
#for ICON_GRID in 0026_R03B07_G; do
module load cdo
module load nco

GRIDDIR=$PERM/grid/
WORKDIR=$SCRATCH/extpar_data/${ICON_GRID}

cd $SCRATCH/extpar/run_scripts

# 
#./dkrz-monmean-t2m.sh ${ICON_GRID}
#./dkrz-monmean-sst-seaice-snow.sh ${ICON_GRID}
#
./ecmwf-monmean-sst4icon.sh ${ICON_GRID} ${GRIDDIR} ${WORKDIR} | tee log_${ICON_GRID}_${TODAY}_ERA-I_T_SEA
./ecmwf-monmean-t2m4icon.sh ${ICON_GRID} ${GRIDDIR} ${WORKDIR} | tee log_${ICON_GRID}_${TODAY}_ERA-I_T_2M

./ecmwf-at-dwd-extpar.sh ${ICON_GRID}
