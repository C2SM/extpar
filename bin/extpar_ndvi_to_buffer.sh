#! /bin/bash
#_______________________________________________________________________________
#
# NDVI processing
# EXAMPLE:
# module load anaconda3/bleeding_edge
# ./extpar_ndvi_to_buffer_new.sh -r NDVI_1998_2003.nc -g icon_grid_0026_R03B07_G.nc -p /home/mpim/m214089/icon_preprocessing/source/extpar_input.2016/
#_______________________________________________________________________________
# disable core dumps
ulimit -c 0
# limit stacksize
ulimit -s unlimited

set -eux
#_______________________________________________________________________________

#
usage() {
    echo "Usage: $0 <options>"
    echo 
    echo "Required Options:"
    echo 
    echo "   -r <NDVI raw data file>"
    echo "   -g <ICON grid file>"
    echo 
    echo "Optional options:"
    echo 
    echo "  [-p <raw data path>]"
}

raw_data_path=""
raw_data_ndvi=""
raw_data_tclim_fine=""
buffer_ndvi=""
icon_grid_file=""
icon_grid_dir=/pool/data/ICON/grids/public/edzw

while getopts ":r:g:p:" opt; do
  case $opt in
    r)
      echo "-c was triggered, Parameter: $OPTARG" >&2
      raw_data_ndvi="$OPTARG"
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

test -f "${raw_data_dir}/${raw_data_ndvi}" || echo "ERROR: NDVI raw data could not be found"
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
#raw_data_ndvi='NDVI_1998_2003.nc'

# buffer output file
buffer_ndvi='NDVI_BUFFER.nc'

# proper CF convention files for checking 
output_ndvi='ndvi_extpar_ICON.nc'

echo "Process ndvi data ... "

cdo -f nc4 -P ${OMP_NUM_THREADS} \
    genycon,$icon_grid_dir/$icon_grid_file:$gridID ${data_dir}${raw_data_ndvi} weights.nc
cdo -f nc4 -P ${OMP_NUM_THREADS} \
    settaxis,1111-01-01,0,1mo -remap,$icon_grid_dir/$icon_grid_file:$gridID,weights.nc ${data_dir}${raw_data_ndvi} ndvi-ycon.nc

${progdir}/cdo2ndvi-buffer.py
mv ndvi-ycon.nc ${output_ndvi}
mv ndvi-ycon_BUFFER.nc ${buffer_ndvi}

echo done
