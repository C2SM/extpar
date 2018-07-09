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

if [ ! -f "$REFOUTDIR/external_parameter"* ] ; then
  echo "No netCDF reference file found"  1>&1
  exit 20 # FAIL
fi

FILE=$(ls -1 ${RUNDIR}/external_parameter*.nc 2>/dev/null)
if [ $? -ne 0 ] ; then
  echo "No netCDF output file found"  1>&1
  exit 20 # FAIL
fi

cdo diffv $REFOUTDIR/external_parameter*.nc $FILE 2>&1 | grep differ
if [ $? -ne 1 ] ; then
   echo "Differences from reference file found" 1>&1
   exit 20 # FAIL
fi

cdo diffv $REFOUTDIR/external_parameter*.nc $FILE 2>&1 | grep Abort
if [ $? -ne 1 ] ; then
   echo "CDO abort during file comparison" 1>&1
   exit 20 # FAIL
fi
# goodbye
exit 0 # MATCH

