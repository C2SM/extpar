                                         #### BATCH_SYSTEM=CRAY ####
########################################################################
#PBS -S /bin/bash
#PBS -q np
#PBS -N extpar
#PBS -o end_extpar
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
ICON_GRID="0099_R19B10"
module load cdo
module load nco

GRIDDIR=$PERM/grid/
WORKDIR=$SCRATCH/extpar_data/${ICON_GRID}
mkdir -p $WORKDIR
cp $SCRATCH/extpar_data/* $WORKDIR/
cd $SCRATCH/extpar_ecmwf/run_scripts

./ecmwf-monmean-sst4icon.sh ${ICON_GRID} ${GRIDDIR} ${WORKDIR} | tee log_${ICON_GRID}_${TODAY}_ERA-I_T_SEA
./ecmwf-monmean-t2m4icon.sh ${ICON_GRID} ${GRIDDIR} ${WORKDIR} | tee log_${ICON_GRID}_${TODAY}_ERA-I_T_2M

./ecmwf-at-dwd-extpar.sh ${ICON_GRID}


cdo -f grb2 copy $WORKDIR/external_parameter_icon_${ICON_GRID}_tiles.nc $WORKDIR/external_parameter_icon_${ICON_GRID}_tiles.g2

cp $WORKDIR/external_parameter_icon_${ICON_GRID}_tiles.nc $PERM/grid/icon_extpar_${ICON_GRID}_L_20180625_tiles.nc
cp $WORKDIR/external_parameter_icon_${ICON_GRID}_tiles.g2 $PERM/grid/icon_extpar_${ICON_GRID}_L_20180625_tiles.g2





