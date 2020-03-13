#!/bin/ksh
      
# import functions to launch Extpar executables
. ./runcontrol_functions.sh
      
ulimit -s unlimited
ulimit -c unlimited

# get hostname
hostname="`echo $HOSTNAME`"
logfile="extpar_runscript.log"

rm ${logfile}

#--------------------------------------------------------------------------------
# define host-dependent paths and variables

# kesch
if [[ $hostname == kesch* || $hostname == daint* ]]; then


    # NetCDF raw data for external parameter
    data_dir=/store/s83/tsm/extpar/raw_data_nc/

    # NetCDF raw data file names
    raw_data_alb='MODIS_month_alb.nc'
    raw_data_alnid='MODIS_month_alnid.nc'
    raw_data_aluvd='MODIS_month_aluvd.nc'

    raw_data_tclim_fine='CRU_T2M_SURF_clim_fine.nc'
    raw_data_glc2000='GLC2000_byte.nc'
    raw_data_glcc='GLCC_usgs_class_byte.nc'
    aster_prefix='ASTER_orig'
    raw_data_flake='GLDB_lakedepth.nc'

# mistral
elif [[ $hostname == m* ]]; then

    # NetCDF raw data for external parameter
    data_dir=/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/extpar_input.2016/

    # NetCDF raw data file names
    raw_data_alb='month_alb.nc'
    raw_data_alnid='month_alnid.nc'
    raw_data_aluvd='month_aluvd.nc'

    raw_data_tclim_fine='CRU_T2M_SURF_clim.nc'
    raw_data_glc2000='glc2000_byte.nc'
    raw_data_glcc='glcc_usgs_class_byte.nc'
    aster_prefix='topo.ASTER_orig'
    raw_data_flake='lakedepth.nc'

# unkown host
else

    # exit script in case of unknown host
    echo ERROR: Unkown host: $hostname >> ${logfile}
    exit 1
fi

# substitute host-dependent namelist parameters
sed -i 's#@raw_data_alb_filename@#'"$raw_data_alb"'#' INPUT_ALB
sed -i 's#@raw_data_alnid_filename@#'"$raw_data_alnid"'#' INPUT_ALB
sed -i 's#@raw_data_aluvd_filename@#'"$raw_data_aluvd"'#' INPUT_ALB
sed -i 's#@raw_data_t_clim_filename@#'"$raw_data_tclim_fine"'#' INPUT_TCLIM
sed -i 's#@raw_data_lu_filename@#'"$raw_data_glc2000"'#' INPUT_LU
sed -i 's#@raw_data_glcc_filename@#'"$raw_data_glcc"'#' INPUT_LU
sed -i 's#@aster_prefix@#'"$aster_prefix"'#g' INPUT_ORO
sed -i 's#@raw_data_flake_filename@#'"$raw_data_flake"'#' INPUT_FLAKE
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define paths and variables independent from host

# directories
currentdir=$(pwd)

# Names of executables

# fortran executables
binary_alb=extpar_alb_to_buffer.exe
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_tclim=extpar_cru_to_buffer.exe
binary_ndvi=extpar_ndvi_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_sgsl=extpar_sgsl_to_buffer.exe
binary_ahf=extpar_ahf_to_buffer.exe
binary_isa=extpar_isa_to_buffer.exe
binary_consistency_check=extpar_consistency_check.exe

# link raw data files to local workdir
ln -s -f ${data_dir}/*.nc .
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define test-specific paths and variables 

type_of_test=`echo $currentdir | rev | cut -d"/" -f2 | rev`
name_of_test=`echo $currentdir | rev | cut -d"/" -f1 | rev`
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# launch extpar executables

echo ">>>> Data will be processed and produced in `pwd` <<<<"

# dwd
if [[ $type_of_test == dwd ]]; then
    # run tclim twice
    run_sequential ${binary_tclim}

    # modify namelist INPUT_TCLIM
    cat > INPUT_TCLIM << EOF_tclim
&t_clim_raw_data
  raw_data_t_clim_path='${data_dir}',
  raw_data_t_clim_filename='${raw_data_tclim_fine}',
  it_cl_type = 1
/  

&t_clim_io_extpar
  t_clim_buffer_file='crutemp_climF_extpar_BUFFER.nc',
  t_clim_output_file='crutemp_climF_extpar_BUFFER.nc'
/  
EOF_tclim

fi

run_parallel ${binary_alb}
run_parallel ${binary_aot}
run_parallel ${binary_tclim}
run_parallel ${binary_lu}
run_parallel ${binary_topo} 
run_parallel ${binary_ndvi} 
run_parallel ${binary_soil} 
run_parallel ${binary_flake}

if [ -f INPUT_AHF ] ; then
    run_parallel ${binary_ahf}
fi
if [ -f INPUT_ISA ] ; then
    run_parallel ${binary_isa}
fi
if [ -f INPUT_SGSL ] ; then
    run_parallel ${binary_sgsl}
fi

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

if [ -f INPUT_AHF ] ; then
    check_exit_status ${binary_ahf} error_count
fi
if [ -f INPUT_ISA ] ; then
    check_exit_status ${binary_isa} error_count
fi
if [ -f INPUT_SGSL ] ; then
    check_exit_status ${binary_sgsl} error_count
fi

# if execution of some Extpar executables failed exit script
if [[ $error_count > 0 ]]; then

    echo "*****************************************"
    echo ""
    echo "Some Extpar executables did not terminate correctly!"
    echo "See ${logfile} for more information"
    echo ""
    echo "*****************************************"
    exit 1 

fi

# dwd
if [[ $type_of_test == dwd ]]; then

    # modify namelist TCLIM_FINAL for consistency check
    cat > INPUT_TCLIM_FINAL << EOF_tclim
&t_clim_raw_data
  raw_data_t_clim_path='${data_dir}',
  raw_data_t_clim_filename='${raw_data_tclim_fine}',
  it_cl_type = 1
/  

&t_clim_io_extpar
  t_clim_buffer_file='crutemp_climF_extpar_BUFFER.nc',
  t_clim_output_file='crutemp_climC_extpar_BUFFER.nc'
/  
EOF_tclim

fi

run_sequential ${binary_consistency_check}

#--------------------------------------------------------------------------------
# clean-up
rm exit_status_*
rm time_*

echo ">>>> External parameters for COSMO model generated <<<<"
