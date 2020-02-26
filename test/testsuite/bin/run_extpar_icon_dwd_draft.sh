#!/bin/ksh -l

# exit on error:
set -e

module load unsupported
module unload intel/14.0.0.080
module load intel/15.0.3.187
module unload python
module load netcdf
module load cdo
module load nco/4.6.2

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
    (( runtime = end - start ))
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

scriptpath=$0
scriptname=${scriptpath##*/}
logfile=${scriptname%.*}_$(date +%Y%m%d%H%M%S).log

#________________________________________________________________________________
# Names of executables

# Scripts
binary_alb=extpar_alb_to_buffer.exe
binary_ndvi=extpar_ndvi_to_buffer.exe
binary_emiss=extpar_emiss_to_buffer.exe
binary_tclim=extpar_cru_to_buffer.exe
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_consistency_check=extpar_consistency_check.exe

#________________________________________________________________________________

export OMP_NUM_THREADS=8

# Assume ASTER orography as default
cp INPUT_ORO_ASTER INPUT_ORO

cp INPUT_TCLIM_COARSE INPUT_TCLIM
run_command ${binary_tclim}
cp INPUT_TCLIM_FINE INPUT_TCLIM

run_command ${binary_tclim}

run_command ${binary_emiss}

run_command ${binary_topo}

run_command ${binary_alb}
run_command ${binary_ndvi}
run_command ${binary_lu}
run_command ${binary_aot}
run_command ${binary_soil}
run_command ${binary_flake}

run_command ${binary_consistency_check}
