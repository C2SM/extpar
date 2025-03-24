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
    # DKRZ machines    
    *levante*)
        if [[ -r /sw/etc/profile.levante ]]; then
            source /sw/etc/profile.levante
        fi
        run_command ./configure.levante.gcc
        run_command source modules.env
        run_command make clean
        echo compile extpar...
        run_command make &> compile.log
        echo          ...done
        echo See compile.log for more information!
        ;;

    # container build run at CO2
    *co2* | *iacdipl-7*)
        run_command podman build -t extpar-base:latest -f Dockerfile.base .
        run_command podman build -t extpar:$ghprbPullId -f Dockerfile.extpar .
        ;;
esac
