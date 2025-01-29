#!/bin/ksh -l

#PBS -w NAME=extpar_prep
#PBS -N ${NAME}
#PBS -S /bin/ksh
#PBS -q %queue_version%
#PBS -m n
#PBS -r n
#PBS -o %loggingfile%
#PBS -j o
#PBS -l elapstim_req=1800

set -x

set -e

set +e
type module 2>&1 >/dev/null
rcmod=$?
set -e
if (( rcmod != 0 )) ; then
  . /usr/share/Modules/init/ksh
fi

module list
module clear <<< 'y'
module list

ulimit -s unlimited
ulimit -c 0

eval $(routine_config)

TMP=${TMPDIR}
cd $TMP

# for grid_id in $(<igrid)
for icon_grid_file in %igrid%
do

work_dir=%iwork_dir%
output_dir=${work_dir}
# @ckoziar, 20190731:
# raw_data_dir=/lustre2/uscratch/jhelmert/ep_netcdf/
raw_data_dir=~routfor/routfox/extpar/rawdata/
grid_dir=%igrid_dir%
progdir=%progdir_lc%

if [[ ! -d ${work_dir} ]] ; then
  mkdir -p ${work_dir} 
fi



scriptpath=$0
scriptname=${scriptpath##*/}
logfile=${scriptname%.*}_$(date +%Y%m%d%H%M%S).log

currentdir=$(pwd)

export PATH=${progdir}:$PATH

cd $work_dir

# icon_grid_file=icon_grid_${grid_id}.nc
grid_id=${icon_grid_file%.nc}
grid_id=${grid_id#icon_grid_}
cp ${grid_dir}/$icon_grid_file .
grid_dir=${grid_id%%_*}

#if [[ ! -f "$icon_grid_file"  ]]
#then
#    if [[ -d /e/rhome/routfor/routfox/icon/grids/public/edzw ]]
#    then
#        cp /e/rhome/routfor/routfox/icon/grids/public/edzw/$icon_grid_file .
#    else
#        wget http://icon-downloads.mpimet.mpg.de/grids/public/mpim/$grid_dir/$icon_grid_file
#    fi
#fi




icon_grid_dir=${work_dir}

#cp $rootdir/test/testsuite/data/mpim/icon_r3b7/INPUT* .
#cp /hpc/rhome/routfor/routfox/extpar/const/INP* .

#GITLAB
cp /hpc/rhome/for0adm/nwp/x86/util/run_scripts/routfor/routfox/const/extpar/INP* .
cp /hpc/rhome/for0adm/nwp/x86/util/run_scripts/routfor/routfox/const/extpar/namelist.py .
#OLD
#cp ${ROUTINE_CONST}/extpar/INP* .

#cp /hpc/uhome/jhelmert/const/extpar/*.py .
# Copy Python Modules from the GIT into the WORKDIR 
##cp /hpc/rhome/routfor/routfox/extpar/job/extpar_*.py .
##cp /hpc/rhome/routfor/routfox/extpar/lib/*.py .

cp /hpc/rhome/for0adm/nwp/x86/util/bin/*.py .

#cp /hpc/rhome/routfor/routfox/extpar/const/namelist.py .
#OLD
#cp ${ROUTINE_CONST}/extpar/namelist.py .

# Copy some ERA FILES
#cp /hpc/rhome/routfor/routfox/extpar/const/*.nc .
cp ${raw_data_dir%/}/ea_oro_1986.nc .
cp ${raw_data_dir%/}/ei_an1986-2015_mean_*.nc .
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Check for UUID to identify operational grids 0026, 0027, 0047 for adaption 
# of special points INPUT_SP_{1,2,3}

uuidOfHGrid_iglo_0026='a27b8de6-18c4-11e4-820a-b5b098c6a5c0'
uuidOfHGrid_ieu_0027='ec13b8bc-b82d-11e4-b13f-4d55411d42e6'
uuidOfHGrid_ilam_0047='c6b12daa-91ad-6404-5b26-c1b6452a2a20'
#-----------------------------------------------------------------------------------------
# Coordinates for grid points Falkenberg (1), Waldstation (2), and Lindenberg MOL/RAO (3) 
#  in ICON-GLOBAL (iglo - 0026), ICON-NEST-EU (ieu - 0027), and ICON-LAM-D2 (ilam - 0047) 
#-----------------------------------------------------------------------------------------
lon_geo_iglo_1=14.1077
lat_geo_iglo_1=52.1215

lon_geo_iglo_2=13.9894
lat_geo_iglo_2=52.1873

lon_geo_iglo_3=14.3130
lat_geo_iglo_3=52.2600
#-----------------------------------------------------------------------------------------
lon_geo_ieu_1=14.1222
lat_geo_ieu_1=52.1665

lon_geo_ieu_2=13.9525
lat_geo_ieu_2=52.1817

lon_geo_ieu_3=14.2320
lat_geo_ieu_3=52.2418
#-----------------------------------------------------------------------------------------
lon_geo_ilam_1=14.1222
lat_geo_ilam_1=52.1665

lon_geo_ilam_2=13.9525
lat_geo_ilam_2=52.1817

lon_geo_ilam_3=14.1199
lat_geo_ilam_3=52.2096

# reads grid uuid from netcdf attributes
uuid=$(ncks -M $icon_grid_file | awk -F\" '/:uuidOfHGrid = / { att_uuid=$(NF-1) } END { print( att_uuid ) }' -)

if [[ "${uuid}" == "${uuidOfHGrid_iglo_0026}" ]] ; then
  print -- "Used grid file is $icon_grid_file with identical UUID as grid number 0026 ICON-GLOBAL \n";

  print -- "Processing Falkenberg - ${lat_geo_iglo_1}, ${lon_geo_iglo_1}";

  cat > INPUT_SP_1 << EOF_SP
&special_points
    lon_geo_sp=${lon_geo_iglo_1},
    lat_geo_sp=${lat_geo_iglo_1},
    soiltype_sp=3.,
    z0_sp=0.03,
    rootdp_sp=0.6,
    plcovmn_sp=0.55,
    plcovmx_sp=0.8,
    laimn_sp=0.5,
    laimx_sp=2.5,
    for_d_sp=-1.,
    for_e_sp=-1.,
    fr_land_sp=-1.,
/
EOF_SP

  print -- "Processing Waldstation Kehrigk - ${lat_geo_iglo_2}, ${lon_geo_iglo_2}";

  cat > INPUT_SP_2 << EOF_SP
&special_points
    lon_geo_sp=${lon_geo_iglo_2},
    lat_geo_sp=${lat_geo_iglo_2},
    soiltype_sp=3.,
    z0_sp=0.91,
    rootdp_sp=0.6,
    plcovmn_sp=0.79,
    plcovmx_sp=0.81,
    laimn_sp=3.0,
    laimx_sp=4.0,
    for_d_sp=-1.,
    for_e_sp=-1.,
    fr_land_sp=-1.,
/
EOF_SP

  print -- "Processing Lindenberg MOL/RAO - ${lat_geo_iglo_3}, ${lon_geo_iglo_3}";

  cat > INPUT_SP_3 << EOF_SP
&special_points
    lon_geo_sp=${lon_geo_iglo_3},
    lat_geo_sp=${lat_geo_iglo_3},
    soiltype_sp=5.,
    z0_sp=0.13,
    rootdp_sp=0.6,
    plcovmn_sp=0.55,
    plcovmx_sp=0.88,
    laimn_sp=0.78,
    laimx_sp=3.18,
    for_d_sp=-1.,
    for_e_sp=-1.,
    fr_land_sp=-1.,
/
EOF_SP

#-----------------------------------------------------------------------------------------

elif [[ "${uuid}" == "${uuidOfHGrid_ieu_0027}" ]] ; then
  print -- "Used grid file is $icon_grid_file with identical UUID as grid number 0027 ICON-EU \n";

  print -- "Processing Falkenberg - ${lat_geo_ieu_1}, ${lon_geo_ieu_1}";

  cat > INPUT_SP_1 << EOF_SP
&special_points
    lon_geo_sp=${lon_geo_ieu_1},
    lat_geo_sp=${lat_geo_ieu_1},
    soiltype_sp=3.,
    z0_sp=0.03,
    rootdp_sp=0.6,
    plcovmn_sp=0.55,
    plcovmx_sp=0.8,
    laimn_sp=0.5,
    laimx_sp=2.5,
    for_d_sp=-1.,
    for_e_sp=-1.,
    fr_land_sp=-1.,
/
EOF_SP

  print -- "Processing Waldstation Kehrigk - ${lat_geo_ieu_2}, ${lon_geo_ieu_2}";

  cat > INPUT_SP_2 << EOF_SP
&special_points
    lon_geo_sp=${lon_geo_ieu_2},
    lat_geo_sp=${lat_geo_ieu_2},
    soiltype_sp=3.,
    z0_sp=0.91,
    rootdp_sp=0.6,
    plcovmn_sp=0.79,
    plcovmx_sp=0.81,
    laimn_sp=3.0,
    laimx_sp=4.0,
    for_d_sp=-1.,
    for_e_sp=-1.,
    fr_land_sp=-1.,
/
EOF_SP

  print -- "Processing Lindenberg MOL/RAO - ${lat_geo_ieu_3}, ${lon_geo_ieu_3}";

  cat > INPUT_SP_3 << EOF_SP
&special_points
    lon_geo_sp=${lon_geo_ieu_3},
    lat_geo_sp=${lat_geo_ieu_3},
    soiltype_sp=5.,
    z0_sp=0.13,
    rootdp_sp=0.6,
    plcovmn_sp=0.55,
    plcovmx_sp=0.88,
    laimn_sp=0.78,
    laimx_sp=3.18,
    for_d_sp=-1.,
    for_e_sp=-1.,
    fr_land_sp=-1.,
/
EOF_SP

#-----------------------------------------------------------------------------------------

elif [[ "${uuid}" == "${uuidOfHGrid_ilam_0047}" ]] ; then
  print -- "Used grid file is $icon_grid_file with identical UUID as grid number 0047 ICON-LAM-D2 \n";

  print -- "Processing Falkenberg - ${lat_geo_ilam_1}, ${lon_geo_ilam_1}";

  cat > INPUT_SP_1 << EOF_SP
&special_points
    lon_geo_sp=${lon_geo_ilam_1},
    lat_geo_sp=${lat_geo_ilam_1},
    soiltype_sp=3.,
    z0_sp=0.03,
    rootdp_sp=0.6,
    plcovmn_sp=0.55,
    plcovmx_sp=0.8,
    laimn_sp=0.5,
    laimx_sp=2.5,
    for_d_sp=-1.,
    for_e_sp=-1.,
    fr_land_sp=-1.,
/
EOF_SP

  print -- "Processing Waldstation Kehrigk - ${lat_geo_ilam_2}, ${lon_geo_ilam_2}";

  cat > INPUT_SP_2 << EOF_SP
&special_points
    lon_geo_sp=${lon_geo_ilam_2},
    lat_geo_sp=${lat_geo_ilam_2},
    soiltype_sp=3.,
    z0_sp=0.91,
    rootdp_sp=0.6,
    plcovmn_sp=0.79,
    plcovmx_sp=0.81,
    laimn_sp=3.0,
    laimx_sp=4.0,
    for_d_sp=-1.,
    for_e_sp=-1.,
    fr_land_sp=-1.,
/
EOF_SP

  print -- "Processing Lindenberg MOL/RAO - ${lat_geo_ilam_3}, ${lon_geo_ilam_3}";

  cat > INPUT_SP_3 << EOF_SP
&special_points
    lon_geo_sp=${lon_geo_ilam_3},
    lat_geo_sp=${lat_geo_ilam_3},
    soiltype_sp=5.,
    z0_sp=0.13,
    rootdp_sp=0.6,
    plcovmn_sp=0.55,
    plcovmx_sp=0.88,
    laimn_sp=0.78,
    laimx_sp=3.18,
    for_d_sp=-1.,
    for_e_sp=-1.,
    fr_land_sp=-1.,
/
EOF_SP
fi


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

# Use version with Fortran software
#mv INPUT_CHECK_FT INPUT_CHECK
# Use version with CDO scripts
#mv INPUT_CHECK_SH INPUT_CHECK
# TCLIM_ files needed by consistency check for merging of TCLIM
#cp INPUT_TCLIM INPUT_TCLIM_COARSE
#cp INPUT_TCLIM INPUT_TCLIM_FINE
#cp INPUT_TCLIM INPUT_TCLIM_FINAL

sed -i 's#@icon_grid_dirname@#'"${icon_grid_dir}"'#' INPUT_ICON_GRID
sed -i 's#@icon_grid_filename@#'"${icon_grid_file}"'#' INPUT_ICON_GRID

#________________________________________________________________________________
# NetCDF raw data for external parameter; adjust the path setting!
# (this case is for mistral.dkrz.de)
#datadir=/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/extpar_input.2016/
datadir=${raw_data_dir}

#________________________________________________________________________________
# NetCDF raw data for external parameter; adjust the path setting!
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

sed -i 's#@orography_buffer_filename@#'"${buffer_topo}"'#' INPUT_ORO_GLOBE
sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_ORO_GLOBE

sed -i 's#@orography_buffer_filename@#'"${buffer_topo}"'#' INPUT_ORO_ASTER
sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_ORO_ASTER

sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_ORO_MERIT

#________________________________________________________________________________
# 3) handle all the remaining files

raw_data_aot='aerosol_optical_thickness.nc'
buffer_aot='extpar_aot_BUFFER.nc'
output_aot='aot_extpar_ICON.nc'

sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_AOT
sed -i 's#@aot_buffer_filename@#'"${buffer_aot}"'#' INPUT_AOT



raw_data_globcover_0='GLOBCOVER_0_16bit.nc'
raw_data_globcover_1='GLOBCOVER_1_16bit.nc'
raw_data_globcover_2='GLOBCOVER_2_16bit.nc'
raw_data_globcover_3='GLOBCOVER_3_16bit.nc'
raw_data_globcover_4='GLOBCOVER_4_16bit.nc'
raw_data_globcover_5='GLOBCOVER_5_16bit.nc'

buffer_lu='extpar_landuse_BUFFER.nc'

sed -i 's#@raw_data_pathname@#'"${datadir}"'#g' INPUT_LU
sed -i 's#@lu_buffer_filename@#'"${buffer_lu}"'#' INPUT_LU

#raw_data_glc2000='glc2000_byte.nc'
#buffer_glc2000='extpar_landuse_BUFFER.nc'
#output_glc2000='extpar_landuse_ICON.nc'

raw_data_glcc='glcc_usgs_class_byte.nc'
buffer_glcc='glcc_landuse_BUFFER.nc'

sed -i 's#@glcc_buffer_filename@#'"${buffer_glcc}"'#' INPUT_LU


raw_data_soil_FAO='FAO_DSMW_DP.nc'
raw_data_soil_HWSD='HWSD0_30_texture_2.nc'
raw_data_deep_soil='HWSD30_100_texture_2.nc'
buffer_soil='soil_BUFFER.nc'

raw_lookup_table_HWSD='LU_TAB_HWSD_UF.data'
raw_HWSD_data='HWSD_DATA_COSMO.data'
raw_HWSD_data_deep='HWSD_DATA_COSMO_S.data'
raw_HWSD_data_extpar='HWSD_DATA_COSMO_EXTPAR.asc'

sed -i 's#@raw_data_pathname@#'"${datadir}"'#g' INPUT_SOIL
sed -i 's#@soil_buffer_filename@#'"${buffer_soil}"'#' INPUT_SOIL



raw_data_flake='lakedepth.nc'
buffer_flake='flake_BUFFER.nc'

sed -i 's#@raw_data_pathname@#'"${datadir}"'#' INPUT_FLAKE
sed -i 's#@flake_buffer_filename@#'"${buffer_flake}"'#' INPUT_FLAKE



#________________________________________________________________________________
# the consistency check requires the output of 
# ${binary_aot}, ${binary_tclim}, ${binary_lu}, ${binary_globe}, 
# ${binary_ndvi}, ${binary_soil} and ${binary_flake}

# output_extpar="external_parameter_icon_${grid_id}_tiles.nc"
# output_extpar_grib="external_parameter_icon_${grid_id}_tiles.g2"
output_extpar="external_parameter_icon_${icon_grid_file%.nc}_tiles.nc"
output_extpar_grib="external_parameter_icon_${icon_grid_file%.nc}_tiles.g2"

sed -i 's#@final_extpar_output_grib@#'"${output_extpar_grib}"'#' INPUT_CHECK
sed -i 's#@final_extpar_output@#'"${output_extpar}"'#' INPUT_CHECK

sed -i 's#@orography_buffer_filename@#'"${buffer_topo}"'#' INPUT_CHECK
sed -i 's#@soil_buffer_filename@#'"${buffer_soil}"'#' INPUT_CHECK
sed -i 's#@lu_buffer_filename@#'"${buffer_lu}"'#' INPUT_CHECK
sed -i 's#@glcc_buffer_filename@#'"${buffer_glcc}"'#' INPUT_CHECK
sed -i 's#@flake_buffer_filename@#'"${buffer_flake}"'#' INPUT_CHECK
sed -i 's#@ndvi_buffer_filename@#'"${buffer_ndvi}"'#' INPUT_CHECK
sed -i 's#@emiss_buffer_filename@#'"${buffer_emiss}"'#' INPUT_CHECK
sed -i 's#@tclim_buffer_filename@#'"${buffer_tclim_fine}"'#' INPUT_CHECK
sed -i 's#@aot_buffer_filename@#'"${buffer_aot}"'#' INPUT_CHECK
sed -i 's#@alb_buffer_filename@#'"${buffer_alb}"'#' INPUT_CHECK

# Alternative Files from ICON-REMAP
sed -i 's#@sst_icon_filename@#'"era_buffer.nc"'#' INPUT_CHECK
sed -i 's#@t2m_icon_filename@#'"era_buffer.nc"'#' INPUT_CHECK

echo ${icon_grid_file} 

case ${icon_grid_file} in
  *_G.nc) sed -i 's#@l_use_glcc@#'".true."'#' INPUT_CHECK
        ;;
  *_N02.nc) sed -i 's#@l_use_glcc@#'".false."'#' INPUT_CHECK
        ;;
  *_L.nc) sed -i 's#@l_use_glcc@#'".false."'#' INPUT_CHECK
        ;;
  *_LN02.nc) sed -i 's#@l_use_glcc@#'".false."'#' INPUT_CHECK
        ;;
  *) sed -i 's#@l_use_glcc@#'".false."'#' INPUT_CHECK
     # @ckoziar, 20190805: add default: l_use_glcc=.false.
     ;;
esac

#ln -sf ${output_dir}/ei_sst_1986-2015_${grid_id}_G_BUFFER.nc .
#ln -sf ${output_dir}/ei_t2m_1986-2015_${grid_id}_G_BUFFER.nc .


done

exit 0
