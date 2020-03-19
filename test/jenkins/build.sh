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
        echo compile extpar...
        make &> compile.log
        echo          ...done
        echo See compile.log for more information!

        ;;
    tsa*)
        source /oprusers/osm/.opr_setup_dir
        export MODULEPATH=$MODULEPATH\:$OPR_SETUP_DIR/modules/modulefiles
        module load PrgEnv-gnu/19.2
        module load netcdf-fortran/4.4.4-gnu-8.3.0-with-system-zlib
        make clean
        echo compile extpar...
        make &> compile.log
        echo          ...done
        echo See compile.log for more information!

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
                module load intel/18.0.4
                ;;
        esac
        module list
        make clean
        echo compile extpar...
        make &> compile.log
        echo          ...done
        echo See compile.log for more information!

        ;;
esac 


