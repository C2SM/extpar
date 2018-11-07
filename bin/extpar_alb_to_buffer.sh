#! /bin/bash
#_______________________________________________________________________________
#
# ALBEDO processing
#
# EXAMPLE:
#
# ./extpar_alb_to_buffer.sh -r month_alb_new.nc \
#                           -u month_aluvd_new.nc \
#                           -i month_alnid_new.nc \
#                           -b month_alb_buffer.nc \
#                           -g icon_grid_0026_R03B07_G.nc \
#                           -p /home/mpim/m214089/icon_preprocessing/source/extpar_input.2016 
#_______________________________________________________________________________
# disable core dumps
ulimit -c 0
# limit stacksize
ulimit -s unlimited

set -eu
#_______________________________________________________________________________

usage() {
    echo "Usage: $0 <options>"
    echo 
    echo "Required Options:"
    echo 
    echo "   -r <alb raw data file>"
    echo "   -u <aluvd raw data file>"
    echo "   -i <alnid raw data file>"
    echo "   -b <output buffer file>"
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

while getopts ":r:u:i:b:g:p:" opt; do
  case $opt in
    r)
      raw_data_alb="$OPTARG"
      ;;
    u)
      raw_data_aluvd="$OPTARG"
      ;;
    i)
      raw_data_alnid="$OPTARG"
      ;;
    b)
      buffer_alb="$OPTARG"
      ;;
    g)
      icon_grid_file="$OPTARG"
      ;;
    
    p)
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

test -f "${raw_data_dir}/${raw_data_alb}"   || echo "ERROR: alb raw data could not be found"
test -f "${raw_data_dir}/${raw_data_aluvd}" || echo "ERROR: aluvd raw data could not be found"
test -f "${raw_data_dir}/${raw_data_alnid}" || echo "ERROR: alnid raw data could not be found"
test -f "${icon_grid_file}"                 || echo "ERROR: ICON grid file could not be found" 

export OMP_NUM_THREADS=8

# proper CF convention files for checking 
output_alb=${buffer_alb/_BUFFER.nc/_ICON.nc}
output_alnid=${output_alb/_ICON.nc/_nid_ICON.nc}
output_aluvd=${output_alb/_ICON.nc/_uvd_ICON.nc}

rm -f weights.nc

cdo -f nc4 -P ${OMP_NUM_THREADS} \
    gendis,$icon_grid_file ${raw_data_dir}/${raw_data_alb} weights.nc
cdo -f nc4 -P ${OMP_NUM_THREADS} \
    setrtoc,-1000000,0.02,0.02 \
    -remap,$icon_grid_file,weights.nc ${raw_data_dir}/${raw_data_alb} alb-dis.nc
cdo -f nc4 -P ${OMP_NUM_THREADS} \
    setrtoc,-1000000,0.02,0.02 \
    -remap,$icon_grid_file,weights.nc ${raw_data_dir}/${raw_data_alnid} alnid-dis.nc
cdo -f nc4 -P ${OMP_NUM_THREADS} \
    setrtoc,-1000000,0.02,0.02 \
    -remap,$icon_grid_file,weights.nc ${raw_data_dir}/${raw_data_aluvd} aluvd-dis.nc

cdo2alb-buffer.py

mv alb-dis.nc ${output_alb}
mv alnid-dis.nc ${output_alnid}
mv aluvd-dis.nc ${output_aluvd}
mv alb-dis_BUFFER.nc ${buffer_alb}

exit 0
