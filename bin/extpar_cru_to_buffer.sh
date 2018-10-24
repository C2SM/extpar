#! /bin/bash
#_______________________________________________________________________________
#
# T_CL processing
#
# EXAMPLE: 
#
# ./extpar_cru_to_buffer.sh -c absolute_hadcrut3.nc \
#                           -f CRU_T2M_SURF_clim.nc \
#                           -b crutemp_climF_extpar_BUFFER.nc \
#                           -g icon_grid_0026_R03B07_G.nc \
#                           -p /home/mpim/m214089/icon_preprocessing/source/extpar_input.2016/
#_______________________________________________________________________________
# disable core dumps
ulimit -c 0
# limit stacksize
ulimit -s unlimited

set -eu
#_______________________________________________________________________________
#
usage() {
    echo "Usage: $0 <options>"
    echo 
    echo "Required Options:"
    echo 
    echo "   -c <coarse CRU raw data file>"
    echo "   -f <fine CRU raw data file>"
    echo "   -b <output buffer file>"
    echo "   -g <ICON grid file>"
    echo 
    echo "Optional options:"
    echo 
    echo "  [-p <raw data path>]"
}

raw_data_path=""
raw_data_tclim_coarse=""
raw_data_tclim_fine=""
buffer_tclim=""
icon_grid_file=""

while getopts ":c:f:b:g:p:" opt; do
  case $opt in
    c)
      raw_data_tclim_coarse="$OPTARG"
      ;;
    f)
      raw_data_tclim_fine="$OPTARG"
      ;;
    b)
      buffer_tclim="$OPTARG"
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

test -f "${raw_data_dir}/${raw_data_tclim_coarse}" || echo "ERROR: Coarse CRU raw data could not be found"
test -f "${raw_data_dir}/${raw_data_tclim_fine}"   || echo "ERROR: Fine CRU raw data could not be found"
test -f "${icon_grid_file}"                         || echo "ERROR: ICON grid file could not be found" 

export OMP_NUM_THREADS=8

# proper CF convention files for checking 
output_tclim=${buffer_tclim/_BUFFER.nc/_ICON.nc}
output_topo='cru_topography_ICON.nc'

# raw_data_tclim_fine:
# Gridsize    Miss :     Minimum        Mean     Maximum : Parameter name
#   259200       0 :     -78.000      159.91      5734.0 : HSURF         
#   259200       0 :      0.0000     0.26011      1.0000 : FR_LAND       
#   259200       0 :      0.0000      73.159      303.23 : T_CL        
# raw_data_tclim_coarse:
# Gridsize    Miss :     Minimum        Mean     Maximum : Parameter name
#     2592       0 :     -48.800      3.1647      32.000 : tem           
#     2592       0 :     -48.700      2.4530      30.800 : tem           
#     2592       0 :     -62.300      2.2681      30.900 : tem           
#     2592       0 :     -68.500      3.5152      32.700 : tem           
#     2592       0 :     -69.200      5.4965      34.100 : tem           
#     2592       0 :     -69.000      6.9265      35.200 : tem           
#     2592       0 :     -70.300      7.2876      35.700 : tem           
#     2592       0 :     -70.200      6.9544      35.000 : tem           
#     2592       0 :     -68.900      5.8637      32.900 : tem           
#     2592       0 :     -60.800      4.6208      30.700 : tem           
#     2592       0 :     -47.000      3.8041      31.000 : tem           
#     2592       0 :     -46.300      3.5452      31.800 : tem      
# crut_coarse_prepared.nc:
# Gridsize    Miss :     Minimum        Mean     Maximum : Parameter name
#   259200       0 :      214.51      277.77      302.07 : tem           
# t_cl-prepared.nc:
# Gridsize    Miss :     Minimum        Mean     Maximum : Parameter name
#   259200       0 :      214.51      277.70      303.23 : T_CL          
# icon_topo-icon_grid.nc: 
# Gridsize    Miss :     Minimum        Mean     Maximum : Parameter name
# 83886080       0 :     -405.00      234.62      8296.2 : HH_TOPO       

cdo -f nc4 -P ${OMP_NUM_THREADS} \
    addc,273.15 \
    -yearmonmean \
    -remapdis,${raw_data_dir}/${raw_data_tclim_fine} \
    ${raw_data_dir}/${raw_data_tclim_coarse} crut_coarse_prepared.nc

cdo expr,'T_CL = ((FR_LAND != 0.0)) ? T_CL : tem; HSURF; FR_LAND;' \
    -merge ${raw_data_dir}/${raw_data_tclim_fine} crut_coarse_prepared.nc \
    t_cl-fine-prepared.nc

cdo -f nc4 -P ${OMP_NUM_THREADS} \
    -chname,topography_c,HH_TOPO \
    -selname,topography_c topography_ICON.nc \
    icon_topo-icon_grid.nc    

cdo -f nc4 -P ${OMP_NUM_THREADS} \
    smooth,maxpoints=16 \
    -setmisstonn \
    -remapdis,$icon_grid_file \
    t_cl-fine-prepared.nc crut_fine-icon_grid.nc

cdo expr,'T_CL = ((FR_LAND != 0.0)) ? T_CL+0.0065*(HSURF-HH_TOPO) : T_CL; HSURF;' \
    -merge crut_fine-icon_grid.nc icon_topo-icon_grid.nc \
    t_cl-dis.nc

cdo2t_cl-buffer.py

mv t_cl-dis_BUFFER.nc ${buffer_tclim}
mv t_cl-dis.nc ${output_tclim}

exit 0
