#! /bin/bash
#_________________________________________________________________________________________________
#SBATCH --job-name=prepare-sst-seaice
#SBATCH --partition=gpu
#SBATCH --exclusive
##SBATCH --nodelist=mg204
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=18
##SBATCH --mem=864G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0287
#SBATCH --output=LOG.prepare-sst-seaice.run.%j
#SBATCH --error=LOG.prepare-sst-seaice.run.%j
#_________________________________________________________________________________________________
#
set -eu
#_________________________________________________________________________________________________
#
declare -A grid_ids=( ["0019_R02B05"]="80km"
                      ["0021_R02B06"]="40km"
                      ["0023_R02B07"]="20km"
                      ["0025_R02B08"]="10km" )

#                      ["0015_R02B09"]="5km" )

for grid_id in ${!grid_ids[@]}
do

work_dir=/scratch/m/m214089/preprocessing/dyamond/${grid_ids[$grid_id]}
output_dir=/home/mpim/m214089/experiments/input/${grid_ids[$grid_id]}
#_________________________________________________________________________________________________
#
cd $work_dir

grid_file=icon_grid_${grid_id}_G.nc
grid_dir=${grid_id%%_*}
output_file=sst-sic-runmean_${grid_id}_G.nc

raw_data_dir=/home/mpim/m214089/icon_preprocessing/source/sst_and_seaice/dyamond/ecmwf-sst-sic-dyamond

if [[ ! -f "$grid_file"  ]]
then
    cp /pool/data/ICON/grids/public/mpim/$grid_dir/$grid_file .
    #wget http://icon-downloads.mpimet.mpg.de/grids/public/mpim/$grid_dir/$grid_file
fi

intermediate_file=sst-sic-runmean.grb
preliminary_file=sst_sic_missing.nc

for file in $(ls -1 $raw_data_dir)
do
    ln -sf $raw_data_dir/$file .
done

rm -f tmp.?? $intermediate_file

for hour in 00 06 12 18; do
  cdo runmean,7 -select,hour=$hour 'sst-sic_*' tmp.$hour
done

cdo mergetime tmp.?? $intermediate_file

rm -f tmp.??

edate=$(date -d "2016-09-14" +%s)
cdate=$(date -d "2016-07-28" +%s)
dates=$(while [[ "$cdate" -le "$edate" ]]; do date -d "@$cdate" +%F; let cdate+=86400; done)

cat > grib2cf.tab << EOF
&parameter name = var34 out_name = SST longname = sea_surface_temperature units = K /
&parameter name = var31 out_name = SIC longname = sea_ice_are_fraction    units = 1 /
EOF

cdo -f nc4 \
    remapnn,$grid_file \
    -setgridtype,regular \
    -setpartabn,grib2cf.tab \
    $intermediate_file $preliminary_file

for d in $dates
do
    cdo seldate,$d $preliminary_file $preliminary_file.$d
done

for d in $dates
do
    cdo setmisstonn $preliminary_file.$d $output_file$d
done

cdo mergetime  ${output_file}?* $output_file

rm -f grib2cf.tab $output_file$d

mv $output_file $output_dir/.

done
#_________________________________________________________________________________________________
