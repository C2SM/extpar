#!/bin/bash

# COSMO TECHNICAL TESTSUITE
#
# This script checks whether the run produced correctly a netCDF
# file at hour zero

# Author       Burkhardt Rockel
# Maintainer   burkhardt.rockel@hzg.de

# check environment variables
RUNDIR=${TS_RUNDIR}
VERBOSE=${TS_VERBOSE}
if [ -z "${VERBOSE}" ] ; then
  echo "Environment variable TS_VERBOSE is not set" 1>&1
  exit 20 # FAIL
fi
if [ -z "${RUNDIR}" ] ; then
  echo "Environment variable TS_RUNDIR is not set" 1>&1
  exit 20 # FAIL
fi
if [ ! -d "${RUNDIR}" ] ; then
  echo "Directory TS_RUNDIR=${RUNDIR} does not exist" 1>&1
  exit 20 # FAIL
fi

FILELIST=$(ls -1 ${RUNDIR}/output/l[bf]f*00.nc 2>/dev/null)
if [ $? -ne 0 ] ; then
  echo "No netCDF output file found"  1>&1
  exit 20 # FAIL
fi
for FILE in ${FILELIST}
do
  if [ ! -s "$FILE" ]; then
    if [ "$VERBOSE" -gt 0 ]; then
      echo "File $FILE is zero size"  1>&1
    fi
    exit 20 # FAIL
  fi
done

# goodbye
exit 0 # MATCH

