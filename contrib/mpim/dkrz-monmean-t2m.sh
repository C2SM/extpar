#! /bin/bash
#_________________________________________________________________________________________________
#SBATCH --job-name=prepare-t2m
#SBATCH --partition=gpu
#SBATCH --exclusive
##SBATCH --nodelist=mg204
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=18
##SBATCH --mem=864G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0287
#SBATCH --output=LOG.prepare-t2m.run.%j
#SBATCH --error=LOG.prepare-t2m.run.%j
#_________________________________________________________________________________________________
#
set -eu
#_________________________________________________________________________________________________
#
declare -A grid_ids=( ["0019_R02B05"]="80km"
                      ["0021_R02B06"]="40km"
                      ["0023_R02B07"]="20km"
                      ["0025_R02B08"]="10km"
                      ["0015_R02B09"]="5km" )

for grid_id in ${!grid_ids[@]}
do

work_dir=/scratch/m/m214089/preprocessing/dyamond/${grid_ids[$grid_id]}
output_dir=/home/mpim/m214089/experiments/input/${grid_ids[$grid_id]}
#_________________________________________________________________________________________________
#
cd $work_dir

grid_file=icon_grid_${grid_id}_G.nc
grid_dir=${grid_id%%_*}

raw_data_dir=/work/mh0287/m214089/extpar_T_input

if [[ ! -f "$grid_file"  ]]
then
    cp /pool/data/ICON/grids/public/mpim/$grid_dir/$grid_file .
    #wget http://icon-downloads.mpimet.mpg.de/grids/public/mpim/$grid_dir/$grid_file
fi

first_year=1986   # first year of climatology
last_year=2015    # last year of climatology

ifs_file="ei_2t_an${first_year}-${last_year}"
ifs_oro_file="ei_oro_${first_year}"
output_file="ei_t2m_${first_year}-${last_year}_${grid_id}_G_BUFFER.nc"

#____________________________________________________________________________________
#

ln -sf ${raw_data_dir}/${ifs_file} ${ifs_file}
ln -sf ${raw_data_dir}/${ifs_oro_file} ${ifs_oro_file}.tmp

# calculate monthly mean over all years and adapt to DWD convention on static boundary conditions

cdo -setyear,1111 -setday,11 -ymonmean ${ifs_file} ${ifs_file}.mean

# convert from geopotential to height (multiply by 1/9.80665)

cdo -setparam,8.2 -mulc,0.102 ${ifs_oro_file}.tmp ${ifs_oro_file}
rm ${ifs_oro_file}.tmp

# interpolate from IFS to ICON grid, on two threads using  nearest-neighbor scalar interpolation

number_of_threads=2
interpolation_method=nn

# using the following parameter table

cat > parameterTable <<EOF
&parameter
  param=167.128
  name=2T
  out_name='T_2M_CLIM'
  long_name="2 metre temperature"
  units="K"
/
&parameter
  param=8.2
  name=HSURF
  out_name='TOPO_CLIM'
  long_name="Geometrical height"
  units="m"
/
EOF

cdo -r -f nc4 -P ${number_of_threads} setpartabp,parameterTable -remap${interpolation_method},${grid_file} ${ifs_oro_file}  ei-oro.nc
cdo -r -f nc4 -P ${number_of_threads} setpartabp,parameterTable -remap${interpolation_method},${grid_file} ${ifs_file}.mean ei.nc

# convert the resulting netcdf files to extpar buffer format for consistency check

$HOME/dyamond/cdo2monmean-t2m.py

mv ei_BUFFER.nc $output_dir/$output_file

rm -f ei-oro.nc ei.nc  ${ifs_file}.mean ${ifs_file} ${ifs_oro_file}.tmp ${ifs_oro_file}

done
