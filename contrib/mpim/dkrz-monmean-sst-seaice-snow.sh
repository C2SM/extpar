#! /bin/bash
#_________________________________________________________________________________________________
#SBATCH --job-name=prepare-sst-seaice-snow
#SBATCH --partition=gpu
#SBATCH --exclusive
##SBATCH --nodelist=mg204
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=18
##SBATCH --mem=864G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0287
#SBATCH --output=LOG.prepare-sst-seaice-snow.run.%j
#SBATCH --error=LOG.prepare-sst-seaice-snow.run.%j
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

ifs_file="ei_an${first_year}-${last_year}"
output_file="ei_sst_${first_year}-${last_year}_${grid_id}_G_BUFFER.nc"

#____________________________________________________________________________________
#

ln -sf ${raw_data_dir}/${ifs_file} ${ifs_file}

# calculate monthly mean over all years and adapt to DWD convention on static boundary conditions

cdo -setyear,1111 -setday,11 -ymonmean ${ifs_file} ${ifs_file}.mean

# interpolate from IFS to ICON grid, on two threads using  nearest-neighbor scalar interpolation

number_of_threads=2
interpolation_method=nn               

# using the following parameter table

cat > parameterTable <<EOF
&parameter
  param=235.128
  name=SKT
  out_name=T_S
  long_name="Skin temperature"
  units="K" 
/
&parameter
  param=34.128
  name=SST
  out_name=T_SEA
  long_name="Sea surface temperature"
  units="K"
/
&parameter
  param=141.128
  name=SD
  out_name=W_SNOW
  long_name="Snow depth"
  units="m of water equivalent"
  standard_name=lwe_thickness_of_surface_snow_amount
/
EOF

cdo -r -f nc4 -P ${number_of_threads} setpartabp,parameterTable -remap${interpolation_method},${grid_file} ${ifs_file}.mean ei.nc

# convert the resulting netcdf file to extpar buffer format for consistency check

$HOME/dyamond/cdo2monmean-sst-seaice-snow.py

mv ei_BUFFER.nc $output_dir/$output_file

rm -f ei.nc  ${ifs_file}.mean ${ifs_file} 

done

