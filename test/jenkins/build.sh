#!/bin/bash

# This is a script for compilation of Extpar by Jenkins on Daint
#

case "$(hostname)" in
daint*)
   module swap PrgEnv-cray PrgEnv-pgi
   module load cray-netcdf
   make clean
   make
   ;;
mlogin*)
   module load gcc/6.2.0
   MACH=mistral.gcc make clean
   MACH=mistral.gcc make -j 4
   ;;
esac 


