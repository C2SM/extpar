#!/bin/bash
set -e

VERBOSITY=${1:-"2"}

./src/testsuite.py -a --exe=run_extpar_cosmo.sh -v ${VERBOSITY}  --testlist=testlist_cosmo.xml --mpicmd='sleep 1 &&'
./src/testsuite.py -a --exe=run_extpar_icon.sh -v ${VERBOSITY}  --testlist=testlist_icon.xml --mpicmd='sleep 1 &&'
./src/testsuite.py -a --exe=run_extpar_icon.sh -v ${VERBOSITY}  --testlist=testlist_landuse.xml --mpicmd='sleep 1 &&'
