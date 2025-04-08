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

case "$(hostname)" in
    *co2* | *iacdipl-7*)
        run_command podman build -t extpar-base:latest -f Dockerfile.base .
        run_command podman build -t extpar:$ghprbPullId -f Dockerfile.extpar .
        ;;
esac
