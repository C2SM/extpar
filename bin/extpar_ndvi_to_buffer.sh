#! /bin/ksh -l
#_______________________________________________________________________________
#
# NDVI processing
#
# EXAMPLE:
#
# ./extpar_ndvi_to_buffer_new.sh -r NDVI_1998_2003.nc \
#                                -b ndvi_buffer.nc \    
#                                -g icon_grid_0026_R03B07_G.nc \
#                                -p /home/mpim/m214089/icon_preprocessing/source/extpar_input.2016/
#_______________________________________________________________________________
# disable core dumps
ulimit -c 0
# limit stacksize
ulimit -s unlimited

module load unsupported
module unload python
module unload cdo
module load netcdf
module load cdo/1.9.5-prerelease 
module load intel/2018.2.199
module load python/3.5.2


set -eu
#_______________________________________________________________________________

usage() {
    echo "Usage: $0 <options>"
    echo 
    echo "Required Options:"
    echo 
    echo "   -r <ndvi raw data file>"
    echo "   -b <output buffer file>"    
    echo "   -g <ICON grid file>"
    echo 
    echo "Optional options:"
    echo 
    echo "  [-p <raw data path>]"
}

raw_data_path=""
raw_data_ndvi=""
buffer_ndvi=""
icon_grid_file=""

while getopts ":r:b:g:p:" opt; do
  case $opt in
    r)
      raw_data_ndvi="$OPTARG"
      ;;
    b)
      buffer_ndvi="$OPTARG"
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

test -f "${raw_data_dir}/${raw_data_ndvi}" || echo "ERROR: ndvi raw data could not be found"
test -f "${icon_grid_file}"                || echo "ERROR: ICON grid file could not be found" 

export OMP_NUM_THREADS=8

# proper CF convention files for checking 
output_ndvi=${buffer_ndvi/_BUFFER.nc/_ICON.nc}

rm -rf weights.nc

cdo -f nc4 -P ${OMP_NUM_THREADS} \
    genycon,$icon_grid_file ${raw_data_dir}/${raw_data_ndvi} weights.nc
cdo -f nc4 -P ${OMP_NUM_THREADS} \
    settaxis,1111-01-01,0,1mo -remap,$icon_grid_file,weights.nc ${raw_data_dir}/${raw_data_ndvi} ndvi-ycon.nc

cdo2ndvi-buffer.py

mv ndvi-ycon.nc ${output_ndvi}
mv ndvi-ycon_BUFFER.nc ${buffer_ndvi}

exit 0
