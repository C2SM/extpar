#!/bin/bash
export PYTHONPATH=@PYTHONPATH@:$PYTHONPATH

logfile=extpar.log
rm $logfile

source runcontrol_functions.sh

for exe in @EXTPAR_EXECUTABLES@
do
    run_sequential $exe
done
