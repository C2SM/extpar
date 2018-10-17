#! /bin/bash

# This is a script for compilation of Extpar by Jenkins slaves

case "$(hostname)" in
    # CSCS machines
    daint*)
        module swap PrgEnv-cray PrgEnv-pgi
        module load cray-netcdf
        make clean
        make
        ;;
    # DKRZ machines    
    mlogin*)
        case "$compiler" in
            gcc)
                export MACH=mistral.gcc
                module unload gcc
                module load gcc/6.2.0
                ;;
            nag)
                export MACH=mistral.nag
                module unload nag
                module load nag/6.2
                ;;
            intel)
                export MACH=mistral.intel
                module unload gcc
                module load gcc/6.2.0
                module unload intel
                module load intel/18.0.2
                ;;
        esac
        make clean
        make -j 4
        ;;
esac 


