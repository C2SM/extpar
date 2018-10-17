#! /bin/bash

# This is a script for compilation of Extpar by Jenkins slaves

case "$(hostname)" in
    daint*)
        module swap PrgEnv-cray PrgEnv-pgi
        module load cray-netcdf
        make clean
        make
        ;;
    mlogin*)
        case "$compiler" in
            gcc)
                export MACH=mistral.gcc
                module purge gcc
                module load gcc/6.2.0
                ;;
            nag)
                export MACH=mistral.nag
                module purge nag
                module load nag/6.2
                ;;
            intel)
                export MACH=mistral.intel
                module purge intel
                module load intel/18.0.2
                ;;
        esac
        make clean
        make -j 4
        ;;
esac 


