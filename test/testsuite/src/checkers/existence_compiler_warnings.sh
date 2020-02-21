#!/bin/bash
# COSMO TECHNICAL TESTSUITE
#
# This script checks whether the compilation of Extpar caused any
# warnings for the NAG and the GCC compiler
#
# Author       Jonas Jucker 
# Maintainer   katherine.osterried@env.ethz.ch

# root dir of Extpar repo
compiledir="../../"

# define compiler warnings to be checked
warnings=("Wunused-dummy-argument" "Wconversion" "Wunused-variable" "Line longer than 132 characters" "Unused dummy variable")

if [ ! -f "$compiledir/compile.log" ] ; then
  echo "No compilation log-file found"  1>&1
  exit 20 # FAIL
fi

for warning in "${warnings[@]}"; do
    echo $warning
    grep $warning $compiledir/compile.log
    if [ $? -ne 1 ] ; then
       echo "Compiler warnings found for $warning" 1>&1
       exit 20 # FAIL
    fi
done

# goodbye
exit 0 # MATCH

