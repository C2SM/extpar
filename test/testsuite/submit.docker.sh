#!/bin/bash
set -e

./src/testsuite.py -a --exe=run_extpar_cosmo.sh -v 2  --testlist=testlist_cosmo.xml --mpicmd='sleep 1 &&'  
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2  --testlist=testlist_icon.xml --mpicmd='sleep 1 &&'  
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2  --testlist=testlist_landuse.xml --mpicmd='sleep 1 &&'  
