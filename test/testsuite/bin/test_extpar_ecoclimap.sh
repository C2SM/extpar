#!/usr/bin/ksh

# import functions to launch Extpar executables
. /store/g142/jcanton/repos/extpar/test/testsuite/bin/runcontrol_functions.sh

ulimit -s unlimited
ulimit -c unlimited

export OMP_NUM_THREADS=4


# get hostname
hostname="`echo $HOSTNAME`"
logfile="extpar_runscript.log"

################################################

# paths to define by user

# Output file format and names
netcdf_output_filename='test_extpar_ecoclimap.nc'

grid_type='COSMO'

# Sandbox (make sure you have enough disk place at that location)!
sandboxdir=/scratch/snx3000/jcanton/extpar_files/test_extpar_ecoclimap

###############################################

# auto-set paths

# directory and name of this file
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
thisScript=`basename "$0"`

# directory of runscripts => need to be as in original repository
scriptdir="/store/g142/jcanton/repos/extpar/run_scripts"
src_python=${scriptdir}/../python/lib

# change dir to src_python to get absolute path
cd $src_python
export PYTHONPATH=$PYTHONPATH:$(pwd)
cd - > /dev/null 2>&1

# directory of compiled extpar executables
exedir=$scriptdir/../bin

source ${scriptdir}/../modules.env
module load daint-gpu
module load CDO
source /project/g110/extpar/venv_daint/bin/activate


#--------------------------------------------------------------------------------
# define host-dependent paths and variables

# CSCS-machines
if [[ $hostname == tsa* || $hostname == arolla* || $hostname == nid* \
    || $hostname == daint* ]]; then
    # NetCDF raw data for external parameter
    data_dir=/store/c2sm/extpar_raw_data/linked_data

else
    # Exit script in case of unknown host
    echo ERROR: Unknown host: $hostname >> ${logfile}
    exit 1

fi

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# define paths and variables independent from host or model

# Names of executables

# python executables
binary_alb=extpar_alb_to_buffer.py
binary_ndvi=extpar_ndvi_to_buffer.py
binary_tclim=extpar_cru_to_buffer.py
binary_ahf=extpar_ahf_to_buffer.py
binary_isa=extpar_isa_to_buffer.py
binary_era=extpar_era_to_buffer.py

# fortran executables
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_consistency_check=extpar_consistency_check.exe

# Output file format and names; adjust!
grib_sample='rotated_ll_pl_grib1'

# -------------------------------------------------------------------------------
# Names of input fields
# -------------------------------------------------------------------------------

# Orography raw data (from cosmo c1)
lsso_param=".TRUE."
ntiles_column=2
ntiles_row=4
topo_files="'ASTER_orig_T006.nc' 'ASTER_orig_T007.nc' 'ASTER_orig_T018.nc' 'ASTER_orig_T019.nc' 'ASTER_orig_T030.nc' 'ASTER_orig_T031.nc' 'ASTER_orig_T042.nc' 'ASTER_orig_T043.nc'"

# Orography smoothing (from cosmo c1)
lsmooth_oro=".TRUE." # .TRUE. not supported for ICON
ilow_pass_oro=1
numfilt_oro=2
eps_filter=1.7

# Albedo
raw_data_alb='alb_new.nc'
raw_data_alnid='alnid_new.nc'
raw_data_aluvd='aluvd_new.nc'
buffer_alb='month_alb_buffer.nc'

# Aerosols
raw_data_aot='aod_AeroCom1.nc'
buffer_aot='extpar_buffer_aot.nc'

# 2 m air temperature climatology
raw_data_tclim_fine='CRU_T_SOIL_clim.nc'
buffer_tclim='crutemp_clim_extpar_buffer.nc'

# GLCC land use
raw_data_glcc='GLCC_usgs_class_byte.nc'
buffer_glcc='glcc_landuse_buffer.nc'

# Globcover 2009 land use
raw_data_globcover_0='GLOBCOVER_0_16bit.nc'
raw_data_globcover_1='GLOBCOVER_1_16bit.nc'
raw_data_globcover_2='GLOBCOVER_2_16bit.nc'
raw_data_globcover_3='GLOBCOVER_3_16bit.nc'
raw_data_globcover_4='GLOBCOVER_4_16bit.nc'
raw_data_globcover_5='GLOBCOVER_5_16bit.nc'

# ESA CCI land use
raw_data_ESACCI_0='ECCI_300m_0.nc'
raw_data_ESACCI_1='ECCI_300m_1.nc'
raw_data_ESACCI_2='ECCI_300m_2.nc'
raw_data_ESACCI_3='ECCI_300m_3.nc'
raw_data_ESACCI_4='ECCI_300m_4.nc'
raw_data_ESACCI_5='ECCI_300m_5.nc'

# ECOCLIMAP-SG land use
raw_data_ecosg_path='/store/g142/jcanton/repos/ECOCLIMAP_SG/'
raw_data_ecosg='ECOCLIMAP_SG.nc'
ln -s -f ${raw_data_ecosg_path}${raw_data_ecosg} ${sandboxdir}/${raw_data_ecosg}

buffer_lu='extpar_landuse_buffer.nc'

# lanczos filter is recommended when activating scale separation
raw_filt_globe_A10='GLOBE_A_filt_lanczos_window.nc'
raw_filt_globe_B10='GLOBE_B_filt_lanczos_window.nc'
raw_filt_globe_C10='GLOBE_C_filt_lanczos_window.nc'
raw_filt_globe_D10='GLOBE_D_filt_lanczos_window.nc'
raw_filt_globe_E10='GLOBE_E_filt_lanczos_window.nc'
raw_filt_globe_F10='GLOBE_F_filt_lanczos_window.nc'
raw_filt_globe_G10='GLOBE_G_filt_lanczos_window.nc'
raw_filt_globe_H10='GLOBE_H_filt_lanczos_window.nc'
raw_filt_globe_I10='GLOBE_I_filt_lanczos_window.nc'
raw_filt_globe_J10='GLOBE_J_filt_lanczos_window.nc'
raw_filt_globe_K10='GLOBE_K_filt_lanczos_window.nc'
raw_filt_globe_L10='GLOBE_L_filt_lanczos_window.nc'
raw_filt_globe_M10='GLOBE_M_filt_lanczos_window.nc'
raw_filt_globe_N10='GLOBE_N_filt_lanczos_window.nc'
raw_filt_globe_O10='GLOBE_O_filt_lanczos_window.nc'
raw_filt_globe_P10='GLOBE_P_filt_lanczos_window.nc'

buffer_topo='topography_buffer.nc'
output_topo='topography_COSMO.nc'

raw_data_sgsl_A10='S_ORO_A10.nc'
raw_data_sgsl_B10='S_ORO_B10.nc'
raw_data_sgsl_C10='S_ORO_C10.nc'
raw_data_sgsl_D10='S_ORO_D10.nc'
raw_data_sgsl_E10='S_ORO_E10.nc'
raw_data_sgsl_F10='S_ORO_F10.nc'
raw_data_sgsl_G10='S_ORO_G10.nc'
raw_data_sgsl_H10='S_ORO_H10.nc'
raw_data_sgsl_I10='S_ORO_I10.nc'
raw_data_sgsl_J10='S_ORO_J10.nc'
raw_data_sgsl_K10='S_ORO_K10.nc'
raw_data_sgsl_L10='S_ORO_L10.nc'
raw_data_sgsl_M10='S_ORO_M10.nc'
raw_data_sgsl_N10='S_ORO_N10.nc'
raw_data_sgsl_O10='S_ORO_O10.nc'
raw_data_sgsl_P10='S_ORO_P10.nc'

buffer_sgsl='sgsl_buffer.nc'

# NDVI
raw_data_ndvi='NDVI_1998_2003.nc'
buffer_ndvi='ndvi_buffer.nc'

# Soil
raw_data_soil_HWSD='HWSD0_30_topsoil.nc'
raw_data_deep_soil='HWSD30_100_subsoil.nc'
buffer_soil='soil_buffer.nc'
raw_lookup_table_HWSD='LU_TAB_HWSD_UF.data'
raw_HWSD_data='HWSD_DATA_COSMO.data'
raw_HWSD_data_deep='HWSD_DATA_COSMO_S.data'
raw_HWSD_data_extpar='HWSD_DATA_COSMO_EXTPAR.asc'

# Lakes
raw_data_flake='GLDB_lakedepth.nc'
buffer_flake='flake_buffer.nc'

# Impervious surface area
raw_data_isa='NOAA_ISA_16bit_lonlat.nc'
buffer_isa='isa_buffer.nc'

# Anthropogenic heat flux
raw_data_ahf='AHF_2006_NOAA_30sec_lonlat.nc' # AHF_2006_NOAA.nc <- this one does not work
buffer_ahf='ahf_buffer.nc'

# ERA climatology (needed by ICON)
raw_data_era_oro='ERA5_ORO_1990.nc'
raw_data_era_sd='ERA5_SD_1990_2019.nc'
raw_data_era_t2m='ERA5_T2M_1990_2019.nc'
raw_data_era_sst='ERA5_SST_1990_2019.nc'
buffer_era='era_buffer.nc'

#--------------------------------------------------------------------------------
# Prepare working directory and create namelists

if [[ ! -d ${sandboxdir} ]] ; then
  mkdir -p ${sandboxdir}
fi

# link raw data files to local workdir
ln -s -f ${data_dir}/*.nc ${sandboxdir}

cp $thisDir/$thisScript ${sandboxdir}/.
cp $exedir/* ${sandboxdir}/.
cd ${sandboxdir}
if [[ -e ${logfile} ]] ; then
  rm -f ${logfile}
fi

cd ${sandboxdir}

echo "\n>>>> Data will be processed and produced in `pwd` <<<<\n"

echo PYTHONPATH: ${PYTHONPATH} >> ${logfile}

# create input namelists
#---
cat > namelist.py << EOF_namelist_python
input_alb = {
        'ialb_type': 1,
        'raw_data_alb_path': '',
        'raw_data_alb_filename': '${raw_data_alb}',
        'raw_data_alnid_filename': '${raw_data_alnid}',
        'raw_data_aluvd_filename': '${raw_data_aluvd}',
        'alb_buffer_file': '${buffer_alb}',
        }

input_tclim = {
        'raw_data_t_clim_path': '',
        'raw_data_tclim_coarse': '',
        'raw_data_tclim_fine': '${raw_data_tclim_fine}',
        't_clim_buffer_file': '${buffer_tclim}',
        'it_cl_type': 1
        }

input_ndvi = {
        'raw_data_ndvi_path': '',
        'raw_data_ndvi_filename': '${raw_data_ndvi}',
        'ndvi_buffer_file': '${buffer_ndvi}',
        }

# input_isa = {
#         'isa_type': 1,
#         'raw_data_isa_path': '',
#         'raw_data_isa_filename': '${raw_data_isa}',
#         'isa_buffer_file': '${buffer_isa}',
#         }
#
# input_ahf = {
#         'iahf_type': 2,
#         'raw_data_ahf_path': '',
#         'raw_data_ahf_filename': '${raw_data_ahf}',
#         'ahf_buffer_file': '${buffer_ahf}',
#         }

input_era = {
        'iera_type': 1,
        'raw_data_era_path': '',
        'raw_data_era_ORO': '${raw_data_era_oro}',
        'raw_data_era_SD': '${raw_data_era_sd}',
        'raw_data_era_T2M': '${raw_data_era_t2m}',
        'raw_data_era_SST': '${raw_data_era_sst}',
        'era_buffer_file': '${buffer_era}',
        }
EOF_namelist_python

# set target grid definition
if [[ $grid_type == 'ICON_10km' ]]; then
    icon_grid_dir='/store/g142/icon_input/grids/CHplus0.2deg_10km'
    icon_grid_nc_file='child_grid_DOM01.nc'
    tile_mode=1
    lradtopo='.FALSE.'
    lsmooth_oro='.FALSE.'
    cat > INPUT_grid_org << EOF_go
&GRID_DEF
 igrid_type = 1,
 domain_def_namelist='INPUT_ICON_GRID'
/
EOF_go
elif [[ $grid_type == 'ICON_1km' ]]; then
    icon_grid_dir='/store/g142/icon_input/grids/CHplus0.2deg'
    icon_grid_nc_file='child_grid_DOM01.nc'
    tile_mode=1
    lradtopo='.FALSE.'
    lsmooth_oro='.FALSE.'
    cat > INPUT_grid_org << EOF_go
&GRID_DEF
 igrid_type = 1,
 domain_def_namelist='INPUT_ICON_GRID'
/
EOF_go
elif [[ $grid_type == 'COSMO' ]]; then
    tile_mode=0
    lradtopo='.TRUE.'
    cat > INPUT_grid_org << EOF_go
&GRID_DEF
 igrid_type = 2,
 domain_def_namelist='INPUT_COSMO_GRID'
/
EOF_go
else
    echo "ERROR: unsupported grid_type"
    exit 1
fi

cat > INPUT_ICON_GRID << EOF_grid
&icon_grid_info
 icon_grid_dir='${icon_grid_dir}'
 icon_grid_nc_file='${icon_grid_nc_file}'
/
EOF_grid
cat > INPUT_COSMO_GRID << EOF_grid
&lmgrid
 pollon=-170.0,
 pollat=43.0,
 startlon_tot=-3.25,
 startlat_tot=-1.70,
 dlon=0.01,
 dlat=0.01,
 ie_tot=400,
 je_tot=315,
/
EOF_grid

cat > INPUT_AOT << EOF_aot
&aerosol_raw_data
  raw_data_aot_path='',
  raw_data_aot_filename='${raw_data_aot}'
  iaot_type=2
/
&aerosol_io_extpar
  aot_buffer_file='${buffer_aot}',
/
EOF_aot

# MCH:
#    raw_data_lu_filename = \
# '${raw_data_globcover_0}' '${raw_data_globcover_1}' \
# '${raw_data_globcover_2}' '${raw_data_globcover_3}' \
# '${raw_data_globcover_4}' '${raw_data_globcover_5}',
# or:
#    raw_data_lu_filename = \
# '${raw_data_ESACCI_0}' '${raw_data_ESACCI_1}' \
# '${raw_data_ESACCI_2}' '${raw_data_ESACCI_3}' \
# '${raw_data_ESACCI_4}' '${raw_data_ESACCI_5}',
# or:
#    raw_data_lu_filename='${raw_data_ecosg}',

cat > INPUT_LU << EOF_lu
&lu_raw_data
   raw_data_lu_path='',
   raw_data_lu_filename='${raw_data_ecosg}',
   i_landuse_data=6,
   ilookup_table_lu=2,
   l_terra_urb=.true.,
/
&lu_io_extpar
   lu_buffer_file='${buffer_lu}',
/
&glcc_raw_data
   raw_data_glcc_path='',
   raw_data_glcc_filename='${raw_data_glcc}'
/
&glcc_io_extpar
   glcc_buffer_file='${buffer_glcc}',
/
EOF_lu

cat > INPUT_ORO << EOF_oro
&oro_runcontrol
  lcompute_sgsl=.FALSE. ,
/
&orography_io_extpar
  orography_buffer_file='${buffer_topo}',
  orography_output_file='${output_topo}'
/
&orography_raw_data
 itopo_type = 2,
 lsso_param = ${lsso_param},
 raw_data_orography_path='',
 ntiles_column = ${ntiles_column},
 ntiles_row = ${ntiles_row},
 topo_files = ${topo_files}
/
EOF_oro
# &sgsl_io_extpar
#  lpreproc_oro=.FALSE.
#  sgsl_buffer_file='${buffer_sgsl}',
#  sgsl_files = \
# '${raw_data_sgsl_A10}' '${raw_data_sgsl_B10}' \
# '${raw_data_sgsl_C10}' '${raw_data_sgsl_D10}' \
# '${raw_data_sgsl_E10}' '${raw_data_sgsl_F10}' \
# '${raw_data_sgsl_G10}' '${raw_data_sgsl_H10}' \
# '${raw_data_sgsl_I10}' '${raw_data_sgsl_J10}' \
# '${raw_data_sgsl_K10}' '${raw_data_sgsl_L10}' \
# '${raw_data_sgsl_M10}' '${raw_data_sgsl_N10}' \
# '${raw_data_sgsl_O10}' '${raw_data_sgsl_P10}'
# /

cat > INPUT_OROSMOOTH << EOF_orosm
&orography_smoothing
  lfilter_oro=${lsmooth_oro},
  ilow_pass_oro=${ilow_pass_oro},
  numfilt_oro=${numfilt_oro},
  ilow_pass_xso=5,
  lxso_first=.FALSE.,
  numfilt_xso=1,
  rxso_mask=750.0,
  eps_filter=${eps_filter},
  rfill_valley=0.0,
  ifill_valley=2
/
EOF_orosm

cat > INPUT_RADTOPO << EOF_radtopo
&radtopo
  lradtopo=${lradtopo}
  itype_scaling=2
  nhori=24
  radius=40000
  min_circ_cov=2
  max_missing=0.95
/
EOF_radtopo

cat > INPUT_SCALE_SEP << EOF_scale_sep
&scale_separated_raw_data
  lscale_separation = .FALSE.,
  raw_data_scale_sep_path = '',
  scale_sep_files = '${raw_filt_globe_A10}' '${raw_filt_globe_B10}'  '${raw_filt_globe_C10}'  '${raw_filt_globe_D10}'  '${raw_filt_globe_E10}'  '${raw_filt_globe_F10}'  '${raw_filt_globe_G10}'  '${raw_filt_globe_H10}'  '${raw_filt_globe_I10}'  '${raw_filt_globe_J10}'  '${raw_filt_globe_K10}'  '${raw_filt_globe_L10}'  '${raw_filt_globe_M10}'  '${raw_filt_globe_N10}'  '${raw_filt_globe_O10}'  '${raw_filt_globe_P10}'
/
EOF_scale_sep

cat > INPUT_SOIL << EOF_soil
&soil_raw_data
 isoil_data = 3,
 ldeep_soil = .FALSE.,
 raw_data_soil_path='',
 raw_data_soil_filename='${raw_data_soil_HWSD}',
 raw_data_deep_soil_filename='${raw_data_deep_soil}'
/
&soil_io_extpar
  soil_buffer_file='${buffer_soil}',
/
&HWSD_index_files
 path_HWSD_index_files='',
 lookup_table_HWSD='${raw_lookup_table_HWSD}',
 HWSD_data='${raw_HWSD_data}',
 HWSD_data_deep='${raw_HWSD_data_deep}',
 HWSD_data_extpar='${raw_HWSD_data_extpar}'
/
EOF_soil

cat > INPUT_FLAKE << EOF_flake
&flake_raw_data
   raw_data_flake_path='',
   raw_data_flake_filename='${raw_data_flake}'
/
&flake_io_extpar
   flake_buffer_file='${buffer_flake}'
/
EOF_flake

# consistency check
cat > INPUT_CHECK << EOF_check
&extpar_consistency_check_io
  grib_output_filename="${grib_output_filename}",
  grib_sample="${grib_sample}",
  netcdf_output_filename="${netcdf_output_filename}",
  i_lsm_data=1,
  land_sea_mask_file="",
  number_special_points=0,
  lwrite_grib=.FALSE.,
  tile_mode=${tile_mode},
/
EOF_check

#--------------------------------------------------------------------------------
# launch extpar executables

run_sequential ${binary_alb}
run_sequential ${binary_aot}
run_sequential ${binary_tclim}
run_sequential ${binary_lu}
run_sequential ${binary_soil}
run_sequential ${binary_flake}
run_sequential ${binary_ndvi}
run_sequential ${binary_topo}
#run_sequential ${binary_isa}
#run_sequential ${binary_ahf}
run_sequential ${binary_era}

# the consistency check requires the output of the previous executables
run_sequential ${binary_consistency_check}

#--------------------------------------------------------------------------------
# clean-up

if [[ -e "exit_status_*" ]] ; then
    rm exit_status_*
fi
if [[ -e "time_*" ]] ; then
    rm time_*
fi

echo ">>>> External parameters for $grid_type model generated <<<<"
