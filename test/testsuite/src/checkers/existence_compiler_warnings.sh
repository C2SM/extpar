#!/bin/bash
# COSMO TECHNICAL TESTSUITE
#
# This script checks whether the Extpar code output is different
# from the reference using the cdo diffv function.  
#
# Author       Katie Osterried
# Maintainer   katherine.osterried@env.ethz.ch

# check environment variables
RUNDIR=${TS_RUNDIR}
VERBOSE=${TS_VERBOSE}
REFOUTDIR=${TS_REFOUTDIR}
compiledir="../../"

# define compiler warnings to be checked
compiler_warnings=("Wunused-dummy-argument" \
                   "Wconversion" \
                   "Wunused-variable" \
                   "Unused dummy variable" \
                   "Line longer than 132 characters")


if [ -z "${VERBOSE}" ] ; then
  echo "Environment variable TS_VERBOSE is not set" 1>&1
  exit 20 # FAIL
fi
if [ -z "${RUNDIR}" ] ; then
  echo "Environment variable TS_RUNDIR is not set" 1>&1
  exit 20 # FAIL
fi
if [ -z "${REFOUTDIR}" ] ; then
  echo "Environment variable TS_REFOUTDIR is not set" 1>&1
  exit 20 # FAIL
fi
if [ ! -d "${RUNDIR}" ] ; then
  echo "Directory TS_RUNDIR=${RUNDIR} does not exist" 1>&1
  exit 20 # FAIL
fi

if [ ! -f "$compiledir/compile.log"* ] ; then
  echo "No compilation log-file found"  1>&1
  exit 20 # FAIL
fi

for warning in ${compiler_warnings[@]}; do
    grep $warning $compiledir/compile.log
    if [ $? -ne 1 ] ; then
       echo "Compiler warnings found for $warning" 1>&1
       exit 20 # FAIL
    fi
done

# goodbye
exit 0 # MATCH

