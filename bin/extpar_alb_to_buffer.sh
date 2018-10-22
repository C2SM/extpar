#! /bin/bash
#_______________________________________________________________________________
#
# ALBEDO processing
# EXAMPLE:
# module load anaconda3/bleeding_edge
# ./extpar_alb_to_buffer.sh -r month_alb_new.nc -u month_aluvd_new.nc -i month_alnid_new.nc  -g icon_grid_0026_R03B07_G.nc -p /home/mpim/m214089/icon_preprocessing/source/extpar_input.2016 
#_______________________________________________________________________________
# disable core dumps
ulimit -c 0
# limit stacksize
ulimit -s unlimited

set -eux
#-----------------------------------------------------------------------------

#
usage() {
    echo "Usage: $0 <options>"
    echo 
    echo "Required Options:"
    echo 
    echo "   -r <ALB raw data file>"
    echo "   -u <ALUVD raw data file>"
    echo "   -i <ALNID raw data file>"
    echo "   -g <ICON grid file>"
    echo 
    echo "Optional options:"
    echo 
    echo "  [-p <raw data path>]"
}

raw_data_path=""
raw_data_alb=""
raw_data_aluvd=""
raw_data_alnid=""
buffer_alb=""
icon_grid_file=""
icon_grid_dir=/pool/data/ICON/grids/public/edzw

while getopts ":r:u:i:g:p:" opt; do
  case $opt in
    r)
      echo "-c was triggered, Parameter: $OPTARG" >&2
      raw_data_alb="$OPTARG"
      ;;
    u)
      echo "-c was triggered, Parameter: $OPTARG" >&2
      raw_data_aluvd="$OPTARG"
      ;;
    i)
      echo "-c was triggered, Parameter: $OPTARG" >&2
      raw_data_alnid="$OPTARG"
      ;;
    g)
      echo "-g was triggered, Parameter: $OPTARG" >&2
      icon_grid_file="$OPTARG"
      ;;
    
    p)
      echo "-p was triggered, Parameter: $OPTARG" >&2        
      raw_data_path="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires arguments" >&2
      usage
      exit 1
      ;;
  esac
done

raw_data_dir="${raw_data_path%%/}"

test -f "${raw_data_dir}/${raw_data_alb}" || echo "ERROR: ALB raw data could not be found"
test -f "${raw_data_dir}/${raw_data_aluvd}" || echo "ERROR: ALUVD raw data could not be found"
test -f "${raw_data_dir}/${raw_data_alnid}" || echo "ERROR: ALNID raw data could not be found"
test -f "${icon_grid_dir}/${icon_grid_file}"                         || echo "ERROR: ICON grid file could not be found" 

data_dir=/home/mpim/m214089/icon_preprocessing/source/extpar_input.2016/
progdir=~/EXTPAR_GIT/bin

filename=0013_R02B04_G
gridID=1
icon_grid_dir=/pool/data/ICON/grids/public/edzw
#icon_grid_file=icon_grid_${filename}.nc

export OMP_NUM_THREADS=8

echo Current working directory: $(pwd)

# raw data: AVAILABLE
#raw_data_alb='month_alb_new.nc'
#raw_data_alnid='month_alnid_new.nc'
#raw_data_aluvd='month_aluvd_new.nc'

# buffer output file
buffer_alb='month_alb_BUFFER.nc'

# proper CF convention files for checking 
output_alb='month_alb_extpar_ICON.nc'
output_alnid='month_alnid_extpar_ICON.nc'
output_aluvd='month_aluvd_extpar_ICON.nc'

echo "Process albedo data ..."

rm -f weights.nc

cdo -f nc4 -P ${OMP_NUM_THREADS} \
    gendis,$icon_grid_dir/$icon_grid_file:$gridID ${data_dir}${raw_data_alb} weights.nc
cdo -f nc4 -P ${OMP_NUM_THREADS} \
    setrtoc,-1000000,0.02,0.02 \
    -remap,$icon_grid_dir/$icon_grid_file:$gridID,weights.nc ${data_dir}${raw_data_alb} alb-dis.nc
cdo -f nc4 -P ${OMP_NUM_THREADS} \
    setrtoc,-1000000,0.02,0.02 \
    -remap,$icon_grid_dir/$icon_grid_file:$gridID,weights.nc ${data_dir}${raw_data_alnid} alnid-dis.nc
cdo -f nc4 -P ${OMP_NUM_THREADS} \
    setrtoc,-1000000,0.02,0.02 \
    -remap,$icon_grid_dir/$icon_grid_file:$gridID,weights.nc ${data_dir}${raw_data_aluvd} aluvd-dis.nc

${progdir}/cdo2alb-buffer.py
mv alb-dis.nc ${output_alb}
mv alnid-dis.nc ${output_alnid}
mv aluvd-dis.nc ${output_aluvd}
mv alb-dis_BUFFER.nc ${buffer_alb}

echo done

