#!/bin/bash

# This is a script for compilation of Extpar by Jenkins on Daint
#

module swap PrgEnv-cray PrgEnv-pgi
module load cray-netcdf

make clean

make


