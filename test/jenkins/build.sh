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
    balfrin*)
        run_command ./configure.balfrin.gcc
        run_command source modules.env
        echo compile extpar...
        run_command make &> compile.log
        echo          ...done
        echo See compile.log for more information!
        ;;

    # DKRZ machines    
    *levante*)
        if [[ -r /sw/etc/profile.levante ]]
        then
           source /sw/etc/profile.levante
        fi
        run_command ./configure.levante.gcc
        run_command source modules.env
        run_command make clean
        echo compile extpar...
        run_command make &> compile.log
        echo          ...done
        echo See compile.log for more information!

esac 

if [[ $compiler == 'python-package' ]]; then
    run_command source modules.env
    run_command python -m venv .venv
    run_command source .venv/bin/activate
    run_command python setup.py sdist
    run_command pip install dist/extpar-*.tar.gz
    run_command python -m extpar.WrapExtpar -h
    deactivate
fi
