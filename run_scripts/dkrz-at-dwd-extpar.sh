#!/bin/ksh -l

set -eu
module load unsupported
module unload intel/14.0.0.080
module load intel/15.0.3.187
module unload python
module load netcdf
module load cdo

module load python/3.5.2

ulimit -s unlimited
ulimit -c 0

#________________________________________________________________________________
error_count=0
run_command()
{
    set +e
    echo ">> Run ${1%% *} ..."    
    start=$(date +%s.%N)
    eval $1 >> ${logfile} 2>&1 
    rc=$?
    printf "   Return code: %i\n" $rc
    end=$(date +%s.%N)
    runtime=$(bc -l <<< "$end - $start")
    if (( rc > 0 ))
    then
        (( error_count++ ))
    fi 
    case $rc in
        0)
            echo "   SUCCESS ${1%% *}"            
            ;;
        127)
            echo "   ERROR ${1%% *}: command not found"
            ;;
        130)
            echo "   ERROR ${1%% *}: script terminated by Ctrl-C"
            ;;             
        *)
            echo "   ERROR ${1%% *}: fatal error - return code $rc"
            ;;
    esac
    echo "   execution time: $runtime s"
    set -e
}

ap_command='time '

#________________________________________________________________________________

for grid_id in $1
do

work_dir=$TMP/${grid_id}
output_dir=${work_dir}
raw_data_dir=$SCR/ep_netcdf/
grid_dir=/e/rhome/routfor/routfox/icon/grids/public/edzw
extpar_dir=~/EXTPAR_GIT/release5.0

if [[ ! -d ${work_dir} ]] ; then
  mkdir -p ${work_dir} 
fi



scriptpath=$0
scriptname=${scriptpath##*/}
logfile=${scriptname%.*}_$(date +%Y%m%d%H%M%S).log

currentdir=$(pwd)
rootdir=${extpar_dir}
progdir=${rootdir}/bin

export PATH=${progdir}:$PATH

cd $work_dir

icon_grid_file=icon_grid_${grid_id}.nc
grid_dir=${grid_id%%_*}

if [[ ! -f "$icon_grid_file"  ]]
then
    if [[ -d /e/rhome/routfor/routfox/icon/grids/public/edzw ]]
    then
        cp /e/rhome/routfor/routfox/icon/grids/public/edzw/$icon_grid_file .
    else
        wget http://icon-downloads.mpimet.mpg.de/grids/public/mpim/$grid_dir/$icon_grid_file
    fi
fi

icon_grid_dir=${work_dir}

cp $rootdir/test/testsuite/data/mpim/icon_r3b7/INPUT* .
# Use version with Fortran software
mv INPUT_CHECK_FT INPUT_CHECK
# Use version with CDO scripts
#mv INPUT_CHECK_SH INPUT_CHECK
# TCLIM_ files needed by consistency check for merging of TCLIM
cp INPUT_TCLIM INPUT_TCLIM_COARSE
cp INPUT_TCLIM INPUT_TCLIM_FINE
cp INPUT_TCLIM INPUT_TCLIM_FINAL

sed -i 's#@icon_grid_dirname@#'"${icon_grid_dir}"'#' INPUT_ICON_GRID
sed -i 's#@icon_grid_filename@#'"${icon_grid_file}"'#' INPUT_ICON_GRID

#________________________________________________________________________________
# NetCDF raw data for external parameter; adjust the path setting!
# (this case is for mistral.dkrz.de)
#datadir=/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/extpar_input.2016/
datadir=${raw_data_dir}

# Output file format and names; adjust!
netcdf_output_filename='external_parameter.nc'
#________________________________________________________________________________
# Names of executables

export OMP_NUM_THREADS=8
# Scripts
#binary_alb=extpar_alb_to_buffer.sh
binary_alb=extpar_alb_to_buffer.exe
#binary_ndvi=extpar_ndvi_to_buffer.sh
binary_ndvi=extpar_ndvi_to_buffer.exe
#binary_tclim=extpar_cru_to_buffer.sh
binary_tclim=extpar_cru_to_buffer.exe
binary_tclim_ref=/lustre2/uwork/jhelmert/EXTPAR_V2_10/tstextpar_cru_to_buffer
#
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe

binary_consistency_check=extpar_consistency_check.exe
#________________________________________________________________________________
if [[ -e ${logfile} ]] ; then
  rm -f ${logfile}
fi
echo ">>>> Data will be processed and produced in $(pwd) <<<<"
#________________________________________________________________________________
# 1) topography needs to be processed first - result is input fro the
#    CRU data processing
raw_data_globe_A10='GLOBE_A10.nc'
raw_data_globe_B10='GLOBE_B10.nc'
raw_data_globe_C10='GLOBE_C10.nc'
raw_data_globe_D10='GLOBE_D10.nc'
raw_data_globe_E10='GLOBE_E10.nc'
raw_data_globe_F10='GLOBE_F10.nc'
raw_data_globe_G10='GLOBE_G10.nc'
raw_data_globe_H10='GLOBE_H10.nc'
raw_data_globe_I10='GLOBE_I10.nc'
raw_data_globe_J10='GLOBE_J10.nc'
raw_data_globe_K10='GLOBE_K10.nc'
raw_data_globe_L10='GLOBE_L10.nc'
raw_data_globe_M10='GLOBE_M10.nc'
raw_data_globe_N10='GLOBE_N10.nc'
raw_data_globe_O10='GLOBE_O10.nc'
raw_data_globe_P10='GLOBE_P10.nc'

buffer_topo='topography_BUFFER.nc'
output_topo='topography_ICON.nc'

sed -i 's#@orography_buffer_filename@#'"${buffer_topo}"'#' INPUT_ORO
sed -i 's#@orography_output_filename@#'"${output_topo}"'#' INPUT_ORO
sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_ORO

run_command ${binary_topo}

#________________________________________________________________________________
# 2) drive the cdo repacement scripts of the failing extpar routines
# because of algorithmic problems for high res output with respect to
# low res source data

raw_data_alb='month_alb_new.nc'
raw_data_alnid='month_alnid_new.nc'
raw_data_aluvd='month_aluvd_new.nc'
buffer_alb='month_alb_BUFFER.nc'
output_alb='month_alb_ICON.nc'

sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_ALB
sed -i 's#@alb_buffer_filename@#'"${buffer_alb}"'#' INPUT_ALB
sed -i 's#@alb_output_filename@#'"${output_alb}"'#' INPUT_ALB

#run_command "${binary_alb} -r ${raw_data_alb} -u ${raw_data_aluvd} -i ${raw_data_alnid} -g ${icon_grid_file} -b ${buffer_alb} -p ${datadir}"
run_command ${binary_alb}

raw_data_ndvi='NDVI_1998_2003.nc'
buffer_ndvi='ndvi_BUFFER.nc'
output_ndvi='ndvi_ICON.nc'

sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_NDVI
sed -i 's#@ndvi_buffer_filename@#'"${buffer_ndvi}"'#' INPUT_NDVI
sed -i 's#@ndvi_output_filename@#'"${output_ndvi}"'#' INPUT_NDVI

#run_command "${binary_ndvi} -r ${raw_data_ndvi} -g ${icon_grid_file} -b ${buffer_ndvi} -p ${datadir}"
run_command "${binary_ndvi}"


raw_data_tclim_coarse='absolute_hadcrut3.nc'
raw_data_tclim_fine='CRU_T2M_SURF_clim.nc'
buffer_tclim_coarse='crutemp_climC_extpar_BUFFER.nc'
output_tclim_coarse='crutemp_climC_extpar_ICON.nc'
buffer_tclim_fine='crutemp_climF_extpar_BUFFER.nc'
output_tclim_fine='crutemp_climF_extpar_ICON.nc'

# for later use in extpar_consistency
sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_TCLIM_COARSE
sed -i 's#@raw_data_tclim@#'"${raw_data_tclim_coarse}"'#' INPUT_TCLIM_COARSE
sed -i 's#@it_cl_type@#'"2"'#' INPUT_TCLIM_COARSE
#ACTIVATE ADDITIONALLY ONLY FOR V2_10, i.e. replace "it_cl_type" with "raw_data_t_id"
#sed -i 's#it_cl_type#'"raw_data_t_id"'#' INPUT_TCLIM_COARSE
#
sed -i 's#@tclim_buffer_filename@#'"${buffer_tclim_coarse}"'#' INPUT_TCLIM_COARSE
sed -i 's#@tclim_output_filename@#'"${output_tclim_coarse}"'#' INPUT_TCLIM_COARSE

cp INPUT_TCLIM_COARSE INPUT_TCLIM

#run_command "${binary_tclim} -c ${raw_data_tclim_coarse} -f ${raw_data_tclim_fine} -g ${icon_grid_file} -b ${buffer_tclim} -p ${datadir}"
run_command "${binary_tclim}"
#run_command "${binary_tclim_ref}"


sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_TCLIM_FINE
sed -i 's#@raw_data_tclim@#'"${raw_data_tclim_fine}"'#' INPUT_TCLIM_FINE
sed -i 's#@it_cl_type@#'"1"'#' INPUT_TCLIM_FINE
sed -i 's#@tclim_buffer_filename@#'"${buffer_tclim_fine}"'#' INPUT_TCLIM_FINE
sed -i 's#@tclim_output_filename@#'"${output_tclim_fine}"'#' INPUT_TCLIM_FINE

cp INPUT_TCLIM_FINE INPUT_TCLIM

run_command "${binary_tclim}"

#________________________________________________________________________________
# 3) handle all the remaining files

raw_data_aot='aerosol_optical_thickness.nc'
buffer_aot='extpar_aot_BUFFER.nc'
output_aot='aot_extpar_ICON.nc'

sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_AOT
sed -i 's#@aot_buffer_filename@#'"${buffer_aot}"'#' INPUT_AOT
sed -i 's#@aot_output_filename@#'"${output_aot}"'#' INPUT_AOT

run_command ${binary_aot}

raw_data_globcover_0='GLOBCOVER_0_16bit.nc'
raw_data_globcover_1='GLOBCOVER_1_16bit.nc'
raw_data_globcover_2='GLOBCOVER_2_16bit.nc'
raw_data_globcover_3='GLOBCOVER_3_16bit.nc'
raw_data_globcover_4='GLOBCOVER_4_16bit.nc'
raw_data_globcover_5='GLOBCOVER_5_16bit.nc'

buffer_lu='extpar_landuse_BUFFER.nc'
output_lu='extpar_landuse_ICON.nc'

sed -i 's#@raw_data_pathname@#'"${datadir}"'#g' INPUT_LU
sed -i 's#@lu_buffer_filename@#'"${buffer_lu}"'#' INPUT_LU
sed -i 's#@lu_output_filename@#'"${output_lu}"'#' INPUT_LU

#raw_data_glc2000='glc2000_byte.nc'
#buffer_glc2000='extpar_landuse_BUFFER.nc'
#output_glc2000='extpar_landuse_ICON.nc'

raw_data_glcc='glcc_usgs_class_byte.nc'
buffer_glcc='glcc_landuse_BUFFER.nc'
output_glcc='glcc_landuse_ICON.nc'

sed -i 's#@glcc_buffer_filename@#'"${buffer_glcc}"'#' INPUT_LU
sed -i 's#@glcc_output_filename@#'"${output_glcc}"'#' INPUT_LU

run_command ${binary_lu}

raw_data_soil_FAO='FAO_DSMW_DP.nc'
raw_data_soil_HWSD='HWSD0_30_texture_2.nc'
raw_data_deep_soil='HWSD30_100_texture_2.nc'
buffer_soil='soil_BUFFER.nc'
output_soil='soil_ICON.nc'

raw_lookup_table_HWSD='LU_TAB_HWSD_UF.data'
raw_HWSD_data='HWSD_DATA_COSMO.data'
raw_HWSD_data_deep='HWSD_DATA_COSMO_S.data'
raw_HWSD_data_extpar='HWSD_DATA_COSMO_EXTPAR.asc'

sed -i 's#@raw_data_pathname@#'"${datadir}"'#g' INPUT_SOIL
sed -i 's#@soil_buffer_filename@#'"${buffer_soil}"'#' INPUT_SOIL
sed -i 's#@soil_output_filename@#'"${output_soil}"'#' INPUT_SOIL

run_command ${binary_soil}

raw_data_flake='lakedepth.nc'
buffer_flake='flake_BUFFER.nc'
output_flake='flake_ICON.nc'

sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_FLAKE
sed -i 's#@flake_buffer_filename@#'"${buffer_flake}"'#' INPUT_FLAKE
sed -i 's#@flake_output_filename@#'"${output_flake}"'#' INPUT_FLAKE

run_command ${binary_flake}

#________________________________________________________________________________
# the consistency check requires the output of 
# ${binary_aot}, ${binary_tclim}, ${binary_lu}, ${binary_globe}, 
# ${binary_ndvi}, ${binary_soil} and ${binary_flake}

output_extpar="external_parameter_icon_${grid_id}_tiles.nc"
output_extpar_grib="external_parameter_icon_${grid_id}_tiles.g2"

sed -i 's#@final_extpar_output_grib@#'"${output_extpar_grib}"'#' INPUT_CHECK
sed -i 's#@final_extpar_output@#'"${output_extpar}"'#' INPUT_CHECK

sed -i 's#@orography_buffer_filename@#'"${buffer_topo}"'#' INPUT_CHECK
sed -i 's#@soil_buffer_filename@#'"${buffer_soil}"'#' INPUT_CHECK
sed -i 's#@lu_buffer_filename@#'"${buffer_lu}"'#' INPUT_CHECK
sed -i 's#@glcc_buffer_filename@#'"${buffer_glcc}"'#' INPUT_CHECK
sed -i 's#@flake_buffer_filename@#'"${buffer_flake}"'#' INPUT_CHECK
sed -i 's#@ndvi_buffer_filename@#'"${buffer_ndvi}"'#' INPUT_CHECK
sed -i 's#@tclim_buffer_filename@#'"${buffer_tclim_fine}"'#' INPUT_CHECK
sed -i 's#@aot_buffer_filename@#'"${buffer_aot}"'#' INPUT_CHECK
sed -i 's#@alb_buffer_filename@#'"${buffer_alb}"'#' INPUT_CHECK

#sed -i 's#@sst_icon_filename@#'"ei_sst_1986-2015_${grid_id}_BUFFER.nc"'#' INPUT_CHECK
#sed -i 's#@t2m_icon_filename@#'"ei_t2m_1986-2015_${grid_id}_BUFFER.nc"'#' INPUT_CHECK

# Alternative Files from ICON-REMAP
sed -i 's#@sst_icon_filename@#'"ei_an1986-2015_${grid_id}_BUFFER.nc"'#' INPUT_CHECK
sed -i 's#@t2m_icon_filename@#'"ei_2t_an1986-2015_${grid_id}_BUFFER.nc"'#' INPUT_CHECK

echo ${icon_grid_file} 

case ${icon_grid_file} in
  *_G*) sed -i 's#@l_use_glcc@#'".true."'#' INPUT_CHECK
        ;;
  *_N02*) sed -i 's#@l_use_glcc@#'".false."'#' INPUT_CHECK
        ;;
  *_L_*) sed -i 's#@l_use_glcc@#'".false."'#' INPUT_CHECK
        ;;
  *_LN02*) sed -i 's#@l_use_glcc@#'".false."'#' INPUT_CHECK
        ;;
esac

#ln -sf ${output_dir}/ei_sst_1986-2015_${grid_id}_G_BUFFER.nc .
#ln -sf ${output_dir}/ei_t2m_1986-2015_${grid_id}_G_BUFFER.nc .

sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_TCLIM_FINAL
sed -i 's#@raw_data_tclim@#'"${raw_data_tclim_fine}"'#' INPUT_TCLIM_FINAL
sed -i 's#@it_cl_type@#'"1"'#' INPUT_TCLIM_FINAL
sed -i 's#@tclim_buffer_filename@#'"${buffer_tclim_fine}"'#' INPUT_TCLIM_FINAL
sed -i 's#@tclim_output_filename@#'"${buffer_tclim_coarse}"'#' INPUT_TCLIM_FINAL







run_command ${binary_consistency_check}

#mv ${output_extpar} ${output_dir}/${output_extpar}

done

#________________________________________________________________________________
#
if (( error_count > 0 ))
then
    echo "CRITICAL: External parameter generation for ICON failed!"
else
    echo ">>>> External parameters for ICON generated <<<<"
fi
#________________________________________________________________________________


