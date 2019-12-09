module load gcc/6.3.0
#ICON_GRID="0026_R03B07_G"
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
qsub ecmwf-monmean-sst4icon.sh ${ICON_GRID} ${GRIDDIR} ${WORKDIR}
#./ecmwf-monmean-t2m4icon.sh ${ICON_GRID} ${GRIDDIR} ${WORKDIR} | tee log_${ICON_GRID}_${TODAY}_ERA-I_T_2M
#./ecmwf-at-dwd-extpar.sh ${ICON_GRID}

