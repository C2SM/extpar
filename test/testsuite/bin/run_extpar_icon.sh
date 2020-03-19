#!/bin/ksh

# import functions to launch Extpar executables
. ./runcontrol_functions.sh

ulimit -s unlimited
ulimit -c 0

# get hostname
hostname="`echo $HOSTNAME`"
logfile="extpar_runscript.log"

rm ${logfile}

#--------------------------------------------------------------------------------
# define host-dependent paths and variables

# mistral
if [[ $hostname == m* ]]; then

    export OMP_NUM_THREADS=8

    # directories
    datadir=/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/extpar_input.2016/
    dir_during_test=./

    # variables for albedo script
    raw_data_alb='month_alb_new.nc'
    raw_data_alnid='month_alnid_new.nc'
    raw_data_aluvd='month_aluvd_new.nc'
    buffer_alb='month_alb_BUFFER.nc'
    output_alb='month_alb_ICON.nc'

    # variables for ndvi script
    raw_data_ndvi='NDVI_1998_2003.nc'
    buffer_ndvi='ndvi_BUFFER.nc'
    output_ndvi='ndvi_ICON.nc'

    # variables for tclim script
    raw_data_tclim_coarse='absolute_hadcrut3.nc'
    raw_data_tclim_fine='CRU_T2M_SURF_clim.nc'
    buffer_tclim='crutemp_clim_extpar_BUFFER.nc'
    output_tclim='crutemp_clim_extpar_ICON.nc'


# unkown host
else

    # exit script in case of unknown host
    echo ERROR: Unkown host: $hostname >> ${logfile}
    exit 1
fi
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define paths and variables independent from host

# directories
currentdir=$(pwd)
rootdir=${currentdir}/../../../../..

# Names of executables

# python and cdo executables
binary_alb=extpar_alb_to_buffer.sh
binary_ndvi=extpar_ndvi_to_buffer.sh
binary_tclim=extpar_cru_to_buffer.sh

# fortran executables
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_emiss=extpar_emiss_to_buffer.exe
binary_consistency_check=extpar_consistency_check.exe

# link raw data files to local workdir
ln -s -f ${datadir}/*.nc .
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define test-specific paths and variables 

type_of_test=`echo $currentdir | rev | cut -d"/" -f2 | rev`
name_of_test=`echo $currentdir | rev | cut -d"/" -f1 | rev`
icon_grid_dir=$rootdir/test/testsuite/data/$type_of_test/$name_of_test/
icon_grid_file=icon_grid*

# mpim
if [[ $type_of_test == mpim ]]; then
    ln -sf ${icon_grid_dir}/ei_sst_an1986-2015_0013_R02B04_G_BUFFER.nc .
    ln -sf ${icon_grid_dir}/ei_t2m_an1986-2015_0013_R02B04_G_BUFFER.nc .

# unknowm test
else

    # exit script in case of unknown host
    echo ERROR: Unkown test: $type_of_test >> ${logfile}
    exit 1
fi

ln -sf ${icon_grid_dir}/${icon_grid_file} .
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# launch extpar executables

echo ">>>> Data will be processed and produced in `pwd` <<<<"

# 1) topography needs to be processed first - result is input fro the
#    CRU data processing

run_sequential ${binary_topo}

#________________________________________________________________________________
# 2) drive the cdo repacement scripts of the failing extpar routines
# because of algorithmic problems for high res output with respect to
# low res source data

run_sequential "${binary_alb} -r ${raw_data_alb} -u ${raw_data_aluvd} -i ${raw_data_alnid} -g ${icon_grid_file} -b ${buffer_alb} -p ${dir_during_test}"

run_sequential "${binary_ndvi} -r ${raw_data_ndvi} -g ${icon_grid_file} -b ${buffer_ndvi} -p ${dir_during_test}"

run_sequential "${binary_tclim} -c ${raw_data_tclim_coarse} -f ${raw_data_tclim_fine} -g ${icon_grid_file} -b ${buffer_tclim} -p ${dir_during_test}"

#________________________________________________________________________________
# 3) handle all the remaining files

run_sequential ${binary_aot}

run_sequential ${binary_lu}

run_sequential ${binary_soil}

run_sequential ${binary_flake}

run_sequential ${binary_emiss}

run_sequential ${binary_consistency_check}
#________________________________________________________________________________

echo ">>>> External parameters for ICON model generated <<<<"
