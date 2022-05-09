#!/bin/ksh

# import functions to launch Extpar executables
. ../test/testsuite/bin/runcontrol_functions.sh

ulimit -s unlimited
ulimit -c unlimited


# get hostname
hostname="`echo $HOSTNAME`"
logfile="extpar_runscript.log"

################################################

# variables to define by user

# Sandbox (make sure you have enough disk place at that location)!
sandboxdir=/scratch/juckerj/sandbox_c7_aster/

# define model for which Extpar should run
model="c7"

###############################################

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


#--------------------------------------------------------------------------------
# define host-dependent paths and variables

# CSCS-machines
if [[ $hostname == tsa* || $hostname == arolla* ]]; then

    # NetCDF raw data for external parameter
    data_dir=/store/c2sm/extpar_raw_data/linked_data

# unkown host
else

    # exit script in case of unknown host
    echo ERROR: Unkown host: $hostname >> ${logfile}
    exit 1
fi

#---------------------------------------------------------------------------------------------------------
# define model (Cosmo1, Cosmo2 or Cosmo7) dependent variables
if [[ $model == "c1" ]]; then

    #output file names
    grib_output_filename='external_parameter_mch_cosmo1.g1'
    netcdf_output_filename='external_parameter_mch_cosmo1.nc'

    # grid definition
    startlon_tot=-9.0
    startlat_tot=-9.0
    dlon=0.01
    dlat=0.01
    ie_tot=1801
    je_tot=1801

    # orography raw data
    lsso_param=".FALSE."
    ntiles_column=2
    ntiles_row=4
    topo_files="'ASTER_orig_T006.nc' 'ASTER_orig_T007.nc' 'ASTER_orig_T018.nc' 'ASTER_orig_T019.nc' 'ASTER_orig_T030.nc' 'ASTER_orig_T031.nc' 'ASTER_orig_T042.nc' 'ASTER_orig_T043.nc'"

    #orography smoothing
    ilow_pass_oro=1
    numfilt_oro=2
    eps_filter=1.7

elif [[ $model == "c2" ]]; then

    #output file names
    grib_output_filename='external_parameter_mch_cosmo2.g1'
    netcdf_output_filename='external_parameter_mch_cosmo2.nc'

    # grid definition
    startlon_tot=-9.0
    startlat_tot=-9.0
    dlon=0.02
    dlat=0.02
    ie_tot=901
    je_tot=901

    # orography raw data
    lsso_param=".TRUE."
    ntiles_column=2
    ntiles_row=4
    topo_files="'ASTER_orig_T006.nc' 'ASTER_orig_T007.nc' 'ASTER_orig_T018.nc' 'ASTER_orig_T019.nc' 'ASTER_orig_T030.nc' 'ASTER_orig_T031.nc' 'ASTER_orig_T042.nc' 'ASTER_orig_T043.nc'"

    #orography smoothing
    ilow_pass_oro=1
    numfilt_oro=2
    eps_filter=1.7

elif [[ $model == "c7" ]]; then

    #output file names
    grib_output_filename='external_parameter_mch_cosmo7.g1'
    netcdf_output_filename='external_parameter_mch_cosmo7.nc'

    # grid definition
    startlon_tot=-18.0
    startlat_tot=-12.9
    dlon=0.06
    dlat=0.06
    ie_tot=601
    je_tot=421

    # orography raw data
    lsso_param=".TRUE."
    ntiles_column=3
    ntiles_row=5
    topo_files="'ASTER_orig_T006.nc' 'ASTER_orig_T007.nc' 'ASTER_orig_T008.nc' 'ASTER_orig_T018.nc' 'ASTER_orig_T019.nc' 'ASTER_orig_T020.nc' 'ASTER_orig_T030.nc' 'ASTER_orig_T031.nc' 'ASTER_orig_T032.nc' 'ASTER_orig_T042.nc' 'ASTER_orig_T043.nc' 'ASTER_orig_T044.nc' 'ASTER_orig_T054.nc' 'ASTER_orig_T055.nc' 'ASTER_orig_T056.nc'"

    #orography smoothing
    ilow_pass_oro=4
    numfilt_oro=1
    eps_filter=0.1

else
    echo $model
    echo " Please specify one of the following models: Cosmo1 => c1, Cosmo2 => c2 or Cosmo7 => c7"
    exit

fi
#---------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define paths and variables independent from host or model

# Names of executables

# python executables
binary_alb=extpar_alb_to_buffer.py
binary_ndvi=extpar_ndvi_to_buffer.py
binary_tclim=extpar_cru_to_buffer.py

# fortran executables
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_ahf=extpar_ahf_to_buffer.exe
binary_isa=extpar_isa_to_buffer.exe
binary_consistency_check=extpar_consistency_check.exe

# Output file format and names; adjust!
grib_sample='rotated_ll_pl_grib1'

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
output_topo='topography_COSMO.nc'

raw_data_ndvi='NDVI_1998_2003.nc'
buffer_ndvi='ndvi_buffer.nc'

raw_data_soil_FAO='FAO_DSMW_double.nc'
raw_data_soil_HWSD='HWSD0_30_topsoil.nc'
raw_data_deep_soil='HWSD30_100_subsoil.nc'
buffer_soil='soil_buffer.nc'

raw_lookup_table_HWSD='LU_TAB_HWSD_UF.data'
raw_HWSD_data='HWSD_DATA_COSMO.data'
raw_HWSD_data_deep='HWSD_DATA_COSMO_S.data'
raw_HWSD_data_extpar='HWSD_DATA_COSMO_EXTPAR.asc'

raw_data_flake='GLDB_lakedepth.nc'
buffer_flake='flake_buffer.nc'

#--------------------------------------------------------------------------------
# Prepare working directory and create namelists

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
EOF_namelist_python

# set target grid definition 
cat > INPUT_grid_org << EOF_go
&GRID_DEF 
 igrid_type = 2,
 domain_def_namelist='INPUT_COSMO_GRID'
/ 
EOF_go
cat > INPUT_COSMO_GRID << EOF_grid
&lmgrid
 pollon=-170.0, 
 pollat=43.0, 
 startlon_tot=${startlon_tot}, 
 startlat_tot=${startlat_tot},
 dlon=${dlon},
 dlat=${dlat},
 ie_tot=${ie_tot},
 je_tot=${je_tot},
/
EOF_grid
#---

cat > INPUT_AOT << EOF_aot
&aerosol_raw_data
  raw_data_aot_path='',
  raw_data_aot_filename='${raw_data_aot}'
/  
&aerosol_io_extpar
  aot_buffer_file='${buffer_aot}',
/
EOF_aot
#---

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
#---
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

cat > INPUT_OROSMOOTH << EOF_orosm
&orography_smoothing
  lfilter_oro=.TRUE.,
  ilow_pass_oro=${ilow_pass_oro},
  numfilt_oro=${numfilt_oro},
  ilow_pass_xso=5,
  lxso_first=.FALSE.,
  numfilt_xso=1,
  rxso_mask=750.0,
  eps_filter=${eps_filter},
  rfill_valley=0.0,
  ifill_valley=7
/
EOF_orosm
#---
cat > INPUT_RADTOPO << EOF_radtopo
&radtopo
  lradtopo=.TRUE.,
  nhori=24
/
EOF_radtopo
#---
cat > INPUT_SCALE_SEP << EOF_scale_sep
&scale_separated_raw_data
  lscale_separation = .FALSE.,
  raw_data_scale_sep_path = '',
  scale_sep_files = '${raw_filt_globe_A10}' '${raw_filt_globe_B10}'  '${raw_filt_globe_C10}'  '${raw_filt_globe_D10}'  '${raw_filt_globe_E10}'  '${raw_filt_globe_F10}'  '${raw_filt_globe_G10}'  '${raw_filt_globe_H10}'  '${raw_filt_globe_I10}'  '${raw_filt_globe_J10}'  '${raw_filt_globe_K10}'  '${raw_filt_globe_L10}'  '${raw_filt_globe_M10}'  '${raw_filt_globe_N10}'  '${raw_filt_globe_O10}'  '${raw_filt_globe_P10}'
/
EOF_scale_sep
#---
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
#---
cat > INPUT_FLAKE << EOF_flake
&flake_raw_data
   raw_data_flake_path='',
   raw_data_flake_filename='${raw_data_flake}'
/
&flake_io_extpar
   flake_buffer_file='${buffer_flake}'
/
EOF_flake
#---
# consistency check
cat > INPUT_CHECK << EOF_check
&extpar_consistency_check_io
  grib_output_filename="${grib_output_filename}",
  grib_sample="${grib_sample}",
  netcdf_output_filename="${netcdf_output_filename}",
  i_lsm_data=1,
  land_sea_mask_file="",
  number_special_points=0,
/  
EOF_check

#--------------------------------------------------------------------------------
# launch extpar executables

run_parallel ${binary_alb}
run_parallel ${binary_aot}
run_parallel ${binary_tclim}
run_parallel ${binary_lu}
run_parallel ${binary_soil} 
run_parallel ${binary_flake}
run_parallel ${binary_ndvi} 
run_parallel ${binary_topo} 

#--------------------------------------------------------------------------------
# IMPORTANT WAIT FOR ALL PARALLEL EXECUTABLES TO END
wait
#--------------------------------------------------------------------------------

# count non-zero exit status
error_count=0

check_exit_status ${binary_alb} error_count
check_exit_status ${binary_aot} error_count
check_exit_status ${binary_tclim} error_count
check_exit_status ${binary_lu} error_count
check_exit_status ${binary_topo}  error_count
check_exit_status ${binary_ndvi}  error_count
check_exit_status ${binary_soil}  error_count
check_exit_status ${binary_flake} error_count

# if execution of some Extpar executables failed exit script
if [[ $error_count > 0 ]]; then

    echo "*****************************************"
    echo ""
    echo "Some Extpar executables did not terminate correctly!"
    echo "See ${logfile} for more information"
    echo ""
    echo "*****************************************"
    exit 

fi

# the consistency check requires the output of 
# ${binary_aot}, ${binary_tclim}, ${binary_lu}, ${binary_globe}, 
# ${binary_ndvi}, ${binary_soil} and ${binary_flake}
run_sequential ${binary_consistency_check}

#--------------------------------------------------------------------------------
# clean-up
rm exit_status_*
rm time_*

echo ">>>> External parameters for COSMO model generated <<<<"
