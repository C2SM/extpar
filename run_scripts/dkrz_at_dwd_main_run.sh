#!/bin/ksh -l
#-----------------------------------------------------------------------------
#PBS -q lc_big
#PBS -l pmem=64gb
#PBS -l pvmem=64gb
#PBS -l select=1:ncpus=8:ompthreads=8:mem=64gb:vmem=64gb
#PBS -l walltime=04:30:00
#PBS -l file=128gb
#PBS -j oe
# ----------------------------------------------------------------------

#ICON_GRID="0026_R03B07_G"
for ICON_GRID in 0027_R03B08_N02 0026_R03B07_G; do
#for ICON_GRID in 0026_R03B07_G; do

GRIDDIR=/e/rhome/routfor/routfox/icon/grids/public/edzw
WORKDIR=$TMP/${ICON_GRID}

cd /e/uhome/jhelmert/TOOLS/runscript/ICON

# 
#./dkrz-monmean-t2m.sh ${ICON_GRID}
#./dkrz-monmean-sst-seaice-snow.sh ${ICON_GRID}
#
./dwd-monmean-sst4icon.sh ${ICON_GRID} ${GRIDDIR} ${WORKDIR} | tee log_${ICON_GRID}_${TODAY}_ERA-I_T_SEA
./dwd-monmean-t2m4icon.sh ${ICON_GRID} ${GRIDDIR} ${WORKDIR} | tee log_${ICON_GRID}_${TODAY}_ERA-I_T_2M
./dkrz-at-dwd-extpar.sh ${ICON_GRID}

done
