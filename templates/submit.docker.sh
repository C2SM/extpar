#!/bin/bash
export PYTHONPATH=@PYTHONPATH@:$PYTHONPATH

logfile=extpar.log
rm $logfile

# source /etc/profile if this has not yet been done
if ! command -v module &> /dev/null; then source /etc/profile; fi

if [ -e modules.env ]; then source modules.env; fi
source runcontrol_functions.sh

for exe in @EXTPAR_EXECUTABLES@
do
    run_sequential $exe
done
