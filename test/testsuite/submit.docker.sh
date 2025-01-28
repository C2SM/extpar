#!/bin/bash
set -e

./src/testsuite.py --exe=run_extpar_cosmo.sh -v 2 -o testsuite.out --testlist=testlist_cosmo.xml --mpicmd='sleep 1 &&' 
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_icon.xml --mpicmd='sleep 1 &&'  
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_landuse.xml --mpicmd='sleep 1 &&'  
./src/testsuite.py -a --exe=run_extpar_icon.sh -v 2 -o testsuite.out --testlist=testlist_art.xml --mpicmd='sleep 1 &&'  
