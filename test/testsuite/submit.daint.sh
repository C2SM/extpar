#!/bin/bash
module load daint-gpu
module load CDO
./src/testsuite.py --exe=run_extpar_daint.sh -v 1 -o testsuite.out 
