#!/bin/bash

#SBATCH --job-name="extpar"
#SBATCH --nodes=1
#SBATCH --output="job.out"
#SBATCH --time=01:00:00
#SBATCH --partition=postproc
#SBATCH --cpus-per-task=12


#--------------------------------------------------------------------------------
# variables to define by user
#--------------------------------------------------------------------------------

# define model for which Extpar should run: c1, c2, i1, i2, i1_dev, i2_dev
model="i1"

# Sandbox (make sure you have enough disk place at that location)!
sandboxdir=$SCRATCH/output_extpar/i1


#--------------------------------------------------------------------------------
# environment
#--------------------------------------------------------------------------------

source /project/g110/extpar/venv_tsa/bin/activate
source ../modules.env
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK # manually set to 1 if run as ./script.sh

export NETCDF_OUTPUT_FILETYPE=NETCDF4

# import functions to launch Extpar executables
. ../test/testsuite/bin/runcontrol_functions.sh

ulimit -s unlimited
ulimit -c unlimited

# get hostname
hostname="`echo $HOSTNAME`"
logfile="extpar_runscript.log"


# auto-set paths

# directory of runscripts => need to be as in original repository
scriptdir=`pwd`
src_python=${scriptdir}/../python/lib

# change dir to src_python to get absolute path
cd $src_python
export PYTHONPATH=$PYTHONPATH:$(pwd)
cd - > /dev/null 2>&1

# directory of compiled extpar executables
exedir=$scriptdir/../bin

# define host-dependent paths and variables
# CSCS-machines
#if [[ $hostname == tsa* || $hostname == arolla* || $hostname == daint*]]; then
if [[ $hostname == tsa* || $hostname == arolla* ]]; then

    # NetCDF raw data for external parameter
    data_dir=/store/c2sm/extpar_raw_data/linked_data

# unkown host
else

    # exit script in case of unknown host
    echo ERROR: Unkown host: $hostname >> ${logfile}
    exit 1
fi

#--------------------------------------------------------------------------------
# define model (COSMO-1, COSMO-2, ICON-1, ICON-2) dependent variables
#--------------------------------------------------------------------------------

if [[ $model == "c1" ]]; then

    #output file names
    netcdf_output_filename='external_parameter_mch_cosmo1.nc'

    # grid definition
    rot_pol_lon=-170.0
    rot_pol_lat=43.0
    startlon_tot=-9.0
    startlat_tot=-9.0
    dlon=0.01
    dlat=0.01
    ie_tot=1801
    je_tot=1801

    lsso_param=".TRUE."
    lsubtract_mean_slope=".FALSE."

    # orography raw data
    ntiles_column=2
    ntiles_row=4
    topo_files="'ASTER_orig_T006.nc' 'ASTER_orig_T007.nc' 'ASTER_orig_T018.nc' 'ASTER_orig_T019.nc' 'ASTER_orig_T030.nc' 'ASTER_orig_T031.nc' 'ASTER_orig_T042.nc' 'ASTER_orig_T043.nc'"

    #orography smoothing
    lsmooth_oro=".TRUE."
    ilow_pass_oro=1
    numfilt_oro=2
    eps_filter=1.7

    # soil: tiles
    itile_mode=0

    # model grid
    model_grid_type=2
    name_model_grid="INPUT_COSMO_GRID"

elif [[ $model == "c2" ]]; then

    #output file names
    netcdf_output_filename='external_parameter_mch_cosmo2.nc'

    # grid definition
    rot_pol_lon=-170.0
    rot_pol_lat=43.0
    startlon_tot=-9.0
    startlat_tot=-9.0
    dlon=0.02
    dlat=0.02
    ie_tot=901
    je_tot=901

    # orography raw data
    lsso_param='.TRUE.'
    ntiles_column=2
    ntiles_row=4
    topo_files="'ASTER_orig_T006.nc' 'ASTER_orig_T007.nc' 'ASTER_orig_T018.nc' 'ASTER_orig_T019.nc' 'ASTER_orig_T030.nc' 'ASTER_orig_T031.nc' 'ASTER_orig_T042.nc' 'ASTER_orig_T043.nc'"

    #orography smoothing
    lsmooth_oro='.TRUE.'
    ilow_pass_oro=1
    numfilt_oro=2
    eps_filter=1.7

    # soil: tiles
    itile_mode=0

    # model grid
    model_grid_type=2
    name_model_grid='INPUT_COSMO_GRID'

elif [[ $model == "i1_dev" ]]; then

    #output file names
    netcdf_output_filename="external_parameter_mch_ICON_1E_R19B08_DOM1.nc"

    # grid definition
    grid_dir="/store/s83/tsm/ICON_INPUT/icon-1e_dev/"
    grid_nc="ICON-1E_DOM01.nc"

    lsso_param='.TRUE.'
    lsubtract_mean_slope='.FALSE.'

    # orography raw data
    ntiles_column=2
    ntiles_row=4
    topo_files="'ASTER_orig_T006.nc' 'ASTER_orig_T007.nc' 'ASTER_orig_T018.nc' 'ASTER_orig_T019.nc' 'ASTER_orig_T030.nc' 'ASTER_orig_T031.nc' 'ASTER_orig_T042.nc' 'ASTER_orig_T043.nc'"

    #orography smoothing
    lsmooth_oro='.FALSE.'
    ilow_pass_oro=1
    numfilt_oro=2
    eps_filter=1.7

    # soil: tiles
    itile_mode=1

    # model grid type
    model_grid_type=1
    name_model_grid='INPUT_ICON_GRID'

elif [[ $model == "i2_dev" ]]; then

    #output file names
    netcdf_output_filename="external_parameter_mch_ICON_2E_R19B07_DOM1.nc"

    # grid definition
    grid_dir="/store/s83/tsm/ICON_INPUT/icon-2e_dev/"
    grid_nc="ICON-2E_DOM01.nc"

    lsso_param='.TRUE.'
    lsubtract_mean_slope='.FALSE.'

    # orography raw data
    ntiles_column=2
    ntiles_row=4
    topo_files="'ASTER_orig_T006.nc' 'ASTER_orig_T007.nc' 'ASTER_orig_T018.nc' 'ASTER_orig_T019.nc' 'ASTER_orig_T030.nc' 'ASTER_orig_T031.nc' 'ASTER_orig_T042.nc' 'ASTER_orig_T043.nc'"

    #orography smoothing
    lsmooth_oro='.FALSE.'
    ilow_pass_oro=1
    numfilt_oro=2
    eps_filter=1.7

    # soil: tiles
    itile_mode=1

    # model grid type
    model_grid_type=1
    name_model_grid='INPUT_ICON_GRID'

elif [[ $model == "i1" ]]; then

    #output file names
    netcdf_output_filename="external_parameter_icon_grid_0001_R19B08_mch.nc"

    # grid definition
    grid_dir="/store/s83/tsm/ICON_INPUT/icon-1e/"
    grid_nc="icon_grid_0001_R19B08_mch.nc"

    lsso_param='.TRUE.'
    lsubtract_mean_slope='.FALSE.'

    # orography raw data
    ntiles_column=2
    ntiles_row=4
    topo_files="'ASTER_orig_T006.nc' 'ASTER_orig_T007.nc' 'ASTER_orig_T018.nc' 'ASTER_orig_T019.nc' 'ASTER_orig_T030.nc' 'ASTER_orig_T031.nc' 'ASTER_orig_T042.nc' 'ASTER_orig_T043.nc'"

    #orography smoothing
    lsmooth_oro='.FALSE.'
    ilow_pass_oro=1
    numfilt_oro=2
    eps_filter=1.7

    # soil: tiles
    itile_mode=1

    # model grid type
    model_grid_type=1
    name_model_grid='INPUT_ICON_GRID'

elif [[ $model == "i2" ]]; then

    #output file names
    netcdf_output_filename="external_parameter_icon_grid_0002_R19B07_mch.nc"

    # grid definition
    grid_dir="/store/s83/tsm/ICON_INPUT/icon-2e/"
    grid_nc="icon_grid_0002_R19B07_mch.nc"

    lsso_param='.TRUE.'
    lsubtract_mean_slope='.FALSE.'

    # orography raw data
    ntiles_column=2
    ntiles_row=4
    topo_files="'ASTER_orig_T006.nc' 'ASTER_orig_T007.nc' 'ASTER_orig_T018.nc' 'ASTER_orig_T019.nc' 'ASTER_orig_T030.nc' 'ASTER_orig_T031.nc' 'ASTER_orig_T042.nc' 'ASTER_orig_T043.nc'"

    #orography smoothing
    lsmooth_oro='.FALSE.'
    ilow_pass_oro=1
    numfilt_oro=2
    eps_filter=1.7

    # soil: tiles
    itile_mode=1

    # model grid type
    model_grid_type=1
    name_model_grid='INPUT_ICON_GRID'


else
    echo $model
    echo " Please specify one of the following models: c1, c2, i1, i2, i1_dev, i2_dev"
    exit

fi

#--------------------------------------------------------------------------------
# define paths and variables independent from host or model
#--------------------------------------------------------------------------------

# Names of executables

# python executables
binary_alb=extpar_alb_to_buffer.py
binary_ndvi=extpar_ndvi_to_buffer.py
binary_emiss=extpar_emiss_to_buffer.py
binary_era=extpar_era_to_buffer.py
binary_ahf=extpar_ahf_to_buffer.py
binary_isa=extpar_isa_to_buffer.py
binary_tclim=extpar_cru_to_buffer.py

# fortran executables
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_consistency_check=extpar_consistency_check.exe

# Names of raw data for INPUT_* namelists
raw_data_alb='alb_new.nc'
raw_data_alnid='alnid_new.nc'
raw_data_aluvd='aluvd_new.nc'
buffer_alb='month_alb_buffer.nc'

raw_data_aot='aot_GACP.nc'
buffer_aot='extpar_buffer_aot.nc'

raw_data_tclim_coarse='CRU_T2M_SURF_clim.nc'
raw_data_tclim_fine='CRU_T_SOIL_clim.nc'
buffer_tclim='crutemp_clim_extpar_buffer.nc'

raw_data_glc2000='GLC2000_byte.nc'
buffer_glc2000='extpar_landuse_buffer.nc'
raw_data_glcc='GLCC_usgs_class_byte.nc'
buffer_glcc='glcc_landuse_buffer.nc'

raw_data_globcover_0='GLOBCOVER_0_16bit.nc'
raw_data_globcover_1='GLOBCOVER_1_16bit.nc'
raw_data_globcover_2='GLOBCOVER_2_16bit.nc'
raw_data_globcover_3='GLOBCOVER_3_16bit.nc'
raw_data_globcover_4='GLOBCOVER_4_16bit.nc'
raw_data_globcover_5='GLOBCOVER_5_16bit.nc'

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
output_topo="topography_${model}.nc"

raw_data_ndvi='NDVI_1998_2003.nc'
buffer_ndvi='ndvi_buffer.nc'

raw_data_soil_FAO='FAO_DSMW_double.nc'
raw_data_soil_HWSD='HWSD0_30_topsoil.nc'
raw_data_deep_soil='HWSD30_100_subsoil.nc'
buffer_soil='soil_buffer.nc'

raw_lookup_table_HWSD='LU_TAB_HWSD_UF.data'
raw_HWSD_data='HWSD_DATA_ICON.data'
raw_HWSD_data_deep='HWSD_DATA_ICON_S.data'
raw_HWSD_data_extpar='HWSD_DATA_ICON_EXTPAR.asc'

raw_data_flake='GLDB_lakedepth.nc'
buffer_flake='flake_buffer.nc'

#--------------------------------------------------------------------------------
# Prepare working directory and create namelists
#--------------------------------------------------------------------------------

if [[ ! -d ${sandboxdir} ]] ; then
  mkdir -p ${sandboxdir}
fi

# link raw data files to local workdir
ln -s -f ${data_dir}/*.nc ${sandboxdir}

cp $scriptdir/* ${sandboxdir}/.
cp $exedir/* ${sandboxdir}/.
cd ${sandboxdir}
if [[ -e ${logfile} ]] ; then
  rm -f ${logfile}
fi

cd ${sandboxdir}

echo "\n>>>> Data will be processed and produced in `pwd` <<<<\n"

echo PYTHONPATH: ${PYTHONPATH} >> ${logfile}

# create input namelists 

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
                
input_emiss = { 
        'iemiss_type': 2,
        'raw_data_emiss_path': '',
        'raw_data_emiss_filename': 'CAMEL_bbe_lw_2010-2015.nc',
        'emiss_buffer_file': 'emiss_buffer.nc'
        }

input_era = {
        'iera_type': 1,
        'raw_data_era_path': '',
        'raw_data_era_ORO': 'ERA5_ORO_1990.nc',
        'raw_data_era_T2M': 'ERA5_T2M_1990_2019.nc',
        'raw_data_era_SST': 'ERA5_SST_1990_2019.nc',
        'raw_data_era_SD': 'ERA5_SD_1990_2019.nc',
        'era_buffer_file': 'era_buffer.nc',
        }

input_ahf = {
        'iahf_type': 1,
        'raw_data_ahf_path': '',
        'raw_data_ahf_filename': 'AHF_2006_CDO.nc',
        'ahf_buffer_file': 'ahf_buffer.nc',
        }

input_isa = {
        'raw_data_isa_path': '',
        'raw_data_isa_filename': 'NOAA_ISA_CDO.nc',
        'isa_buffer_file': 'isa_buffer.nc',
        }
EOF_namelist_python

# set target grid definition 
cat > INPUT_grid_org << EOF_go
&GRID_DEF 
 igrid_type=${model_grid_type},
 domain_def_namelist='${name_model_grid}'
/ 
EOF_go

# COSMO grid
cat > INPUT_COSMO_GRID << EOF_grid_cosmo
&lmgrid
 pollon=${rot_pol_lon}, 
 pollat=${rot_pol_lat},
 startlon_tot=${startlon_tot}, 
 startlat_tot=${startlat_tot},
 dlon=${dlon},
 dlat=${dlat},
 ie_tot=${ie_tot},
 je_tot=${je_tot},
/
EOF_grid_cosmo

# ICON grid
cat > INPUT_ICON_GRID << EOF_grid_icon
&icon_grid_info
  icon_grid_dir="${grid_dir}",
  icon_grid_nc_file="${grid_nc}",
/
EOF_grid_icon

cat > INPUT_AOT << EOF_aot
&aerosol_raw_data
  raw_data_aot_path='./',
  raw_data_aot_filename='${raw_data_aot}'
/  
&aerosol_io_extpar
  aot_buffer_file='${buffer_aot}',
/
EOF_aot

cat > INPUT_LU << EOF_lu
&lu_raw_data
  raw_data_lu_path='',
  raw_data_lu_filename='${raw_data_globcover_0}' '${raw_data_globcover_1}' '${raw_data_globcover_2}' '${raw_data_globcover_3}' '${raw_data_globcover_4}' '${raw_data_globcover_5}',
  i_landuse_data=1,
  ilookup_table_lu=1 
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
itopo_type=2,
lsso_param=${lsso_param},
raw_data_orography_path='./',
ntiles_column=${ntiles_column},
ntiles_row=${ntiles_row},
topo_files=${topo_files}
/
EOF_oro

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
  lradtopo=.TRUE.
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
  netcdf_output_filename="${netcdf_output_filename}",
  tile_mode=${itile_mode},
  lwrite_netcdf=.TRUE.
  i_lsm_data=1,
  land_sea_mask_file="",
  number_special_points=0,
/  
EOF_check

#--------------------------------------------------------------------------------
# launch extpar executables
#--------------------------------------------------------------------------------

# 1. run topography first (is needed for CRU data processing)
#--------------------------------
run_sequential ${binary_topo}    # topography

# 2. run rest which is needed for ICON
#--------------------------------
run_sequential ${binary_alb}      # albedo
run_sequential ${binary_aot}      # aerosol optical thickness
run_sequential ${binary_emiss}    # emissivity
run_sequential ${binary_era}      # era climatology
run_sequential ${binary_flake}    # fraction lake
run_sequential ${binary_lu}       # land use
run_sequential ${binary_ndvi}     # normalized difference vegetation index
run_sequential ${binary_tclim}    # cru: temperature climatology
run_sequential ${binary_soil}     # soil


# 3. binaries which are normally not needed for icon-nwp
#--------------------------------
# run_sequential ${binary_ahf}    # anthropogenic heat flux
# run_sequential ${binary_isa}    # impervious surface area


# 4. finally: run consistency check
#--------------------------------
run_sequential ${binary_consistency_check}

#--------------------------------------------------------------------------------
# clean-up
#--------------------------------------------------------------------------------
rm exit_status_*
rm time_*

echo ">>>> External parameters for COSMO/ICON model generated <<<<"
