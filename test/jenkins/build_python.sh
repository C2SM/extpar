#! /bin/bash

# This is a script for compilation of Extpar by Jenkins slaves

# Define run_command function
function run_command {
        "$@"
        local status=$?
        if [ $status -ne 0 ]; then
           echo "error with $1" >&2
           exit 1
        fi
        return $status
}


##############################################################
# Begin script

case "$(hostname)" in
    # CSCS machines
    daint*)
        run_command git submodule init
        run_command git submodule update
        run_command ./configure.daint.gcc
        run_command source modules.env
        run_command make clean
        echo compile extpar...
        run_command make &> compile.log
        echo          ...done
        echo See compile.log for more information!
        echo Building python package
        for f in $(find bin -type l); do cp --remove-destination $(readlink $f) $f; done
        run_command python setup.py bdist_wheel
        ;;

    tsa*)
        run_command git submodule init
        run_command git submodule update
        run_command ./configure.tsa.gcc
        run_command source modules.env
        run_command make clean
        echo compile extpar...
        run_command make &> compile.log
        echo          ...done
        echo See compile.log for more information!
        echo Building python package
        for f in $(find bin -type l); do cp --remove-destination $(readlink $f) $f; done
        run_command python setup.py bdist_wheel
        ;;

    # DKRZ machines
    *levante*)
        if [[ -r /sw/etc/profile.levante ]]
        then
           source /sw/etc/profile.levante
        fi
        run_command git submodule init
        run_command git submodule update
        run_command ./configure.levante.$compiler
        run_command source modules.env
        run_command make clean
        echo compile extpar...
        run_command make &> compile.log
        echo          ...done
        echo See compile.log for more information!
        echo Building python package
        for f in $(find bin -type l); do cp --remove-destination $(readlink $f) $f; done
        run_command python setup.py bdist_wheel
        ;;
esac
