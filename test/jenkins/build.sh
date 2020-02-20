#! /bin/bash

# This is a script for compilation of Extpar by Jenkins slaves

case "$(hostname)" in
    # CSCS machines
    daint*)
        module swap PrgEnv-cray PrgEnv-pgi
        module load cray-netcdf
        module list
        make clean
        make &> compile.log
        ;;
    kesch*)
        export MODULEPATH=$MODULEPATH:/oprusers/owm/modules/RH7.5/modulefiles
        module load PE/17.06
        module load gcc
        module load netcdf/4.4.1.1-gmvolf-17.02
        module load cdo
        module list
        make clean
        make &> compile.log

        ;;
    # DKRZ machines    
    mlogin*)
        if [[ -r /sw/rhel6-x64/etc/profile.mistral ]]
        then
           source /sw/rhel6-x64/etc/profile.mistral
        fi
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
        module list
        make clean
        make &> compile.log

        ;;
esac 


