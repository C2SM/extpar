#!/bin/bash

# This script runs the Extpar testsuite on Daint
#


cd test/testsuite

cd data
./get_data.sh
cd ..

./submit.daint.sh
