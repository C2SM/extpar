#!/bin/ksh
set -x
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
    currentdir=$(pwd)
    rootdir=${currentdir}/../../../../..
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

# Names of executables


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
# define test-specific paths and variables 

type_of_test=`echo $currentdir | rev | cut -d"/" -f2 | rev`
name_of_test=`echo $currentdir | rev | cut -d"/" -f1 | rev`
icon_grid_dir=$rootdir/test/testsuite/data/$type_of_test/$name_of_test/
icon_grid_file=icon_grid*

# mpim
if [[ $type_of_test == mpim ]]; then

    # python and cdo executables
    binary_alb=extpar_alb_to_buffer.sh
    binary_ndvi=extpar_ndvi_to_buffer.sh
    binary_tclim=extpar_cru_to_buffer.sh

    ln -sf ${icon_grid_dir}/ei_sst_an1986-2015_0013_R02B04_G_BUFFER.nc .
    ln -sf ${icon_grid_dir}/ei_t2m_an1986-2015_0013_R02B04_G_BUFFER.nc .

# ecmwf
elif [[ $type_of_test == ecmwf ]]; then

    # fortran executables
    binary_alb=extpar_alb_to_buffer.exe
    binary_ndvi=extpar_ndvi_to_buffer.exe
    binary_tclim=extpar_cru_to_buffer.exe

    # run TCLIM with COARSE and FINE
    cp INPUT_TCLIM_COARSE INPUT_TCLIM

    ln -sf ${icon_grid_dir}/ei_an1986-2015_0099_R19B10_BUFFER.nc
    ln -sf ${icon_grid_dir}/ei_2t_an1986-2015_0099_R19B10_BUFFER.nc

# unknowm test
else

    # exit script in case of unknown host
    echo ERROR: Unkown test: $type_of_test >> ${logfile}
    exit 1
fi

ln -sf ${icon_grid_dir}/${icon_grid_file} .

#--------------------------------------------------------------------------------
# launch extpar executables

echo ">>>> Data will be processed and produced in `pwd` <<<<"

# 1) topography needs to be processed first - result is input fro the
#    CRU data processing

run_sequential ${binary_topo}

#________________________________________________________________________________

# mpim
if [[ $type_of_test == mpim ]]; then

    # 2) drive the cdo repacement scripts of the failing extpar routines
    # because of algorithmic problems for high res output with respect to
    # low res source data
    run_sequential "${binary_alb} -r ${raw_data_alb} -u ${raw_data_aluvd} -i ${raw_data_alnid} -g ${icon_grid_file} -b ${buffer_alb} -p ${dir_during_test}"

    run_sequential "${binary_ndvi} -r ${raw_data_ndvi} -g ${icon_grid_file} -b ${buffer_ndvi} -p ${dir_during_test}"

    run_sequential "${binary_tclim} -c ${raw_data_tclim_coarse} -f ${raw_data_tclim_fine} -g ${icon_grid_file} -b ${buffer_tclim} -p ${dir_during_test}"

# ecmwf
elif [[ $type_of_test == ecmwf ]]; then

    # 2) use conventional fortran executables
    run_sequential ${binary_alb}

    run_sequential ${binary_tclim}

    cp INPUT_TCLIM_FINE INPUT_TCLIM
    run_sequential ${binary_tclim}

    run_sequential ${binary_ndvi}
    
# other tests
else
    # 2) use conventional fortran executables
    run_sequential ${binary_alb}

    run_sequential ${binary_ndvi}

    run_sequential ${binary_tclim}

fi

#________________________________________________________________________________
# 3) handle all the remaining files

run_sequential ${binary_aot}

run_sequential ${binary_lu}

run_sequential ${binary_soil}

run_sequential ${binary_flake}

if [ -f INPUT_EMISS ] ; then
    run_sequential ${binary_emiss}
fi

# ecmwf
if [[ $type_of_test == ecmwf ]]; then

    cp INPUT_TCLIM_FINAL INPUT_TCLIM
    run_sequential ${binary_consistency_check}

# other tests
else

    run_sequential ${binary_consistency_check}
fi
#________________________________________________________________________________

echo ">>>> External parameters for ICON model generated <<<<"
