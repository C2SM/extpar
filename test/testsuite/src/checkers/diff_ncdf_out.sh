#!/bin/bash

# COSMO TECHNICAL TESTSUITE
#
# This script checks whether the Extpar code output is different
# from the reference using the cdo diffv function.  

# Fields allowed to be not bit-identical can be declared in fields_with_roundoff
#
# The following fields are not bit-identical for Tsa and Daint
# in comparison to Mistral
# because of different versions of Python, CDO and GCC:
#
#   -NDVI
#   -NDVI_MAX
#   -NDVI_MRAT
#   -ALB_DIF12
#   -ALNID12
#
# ############ Remarks ##############

#  Fields with round-offs are allowed to have 
#  round-off in all tests on all machines!

#  Test can go wrong if CDO outputs another output
#  than defined in pattern_to_ignore

#
# Author       Katie Osterried, Jonas Jucker
# Maintainer   katherine.osterried@env.ethz.ch

log='diffv.log'

# define fields that can have round-offs smaller than 0.001
# names that occur in other names (NDVI -> NDVI_MRAT) must
# have a leading and trailing whitespace
fields_with_roundoff=(" NDVI " "NDVI_MRAT" "NDVI_MAX" "ALB_DIF12" "ALNID12" )

# CDO output that is not interesting for this checker
pattern_to_ignore="diffn\|differ\|nhori\|Warning"

# clean-up old stuff
[[ ! -f $log ]] ||  rm ${log}

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


# compare fields and write output to logfile
cdo --sortname diffv $REFOUTDIR/external_parameter*.nc $FILE &> $log 

# check for Abort
grep Abort $log > /dev/null
if [ $? -ne 1 ] ; then
   echo "CDO abort during file comparison" 1>&1
   exit 20 # FAIL
fi

# check if fields differ
grep 'differ'  $log > /dev/null

# fields do not differ
if [ $? -ne 0 ] ; then
    exit 0 # MATCH

else # fields differ

    # check if fields differ more than 0.001
    grep '0 of'  $log > /dev/null

    # fields differ more than 0.001
    if [ $? -ne 0 ] ; then
        echo "Differences from reference file found" 1>&1
        [[ $VERBOSE -gt 1 ]] && grep -v "$pattern_to_ignore" $log 1>&1
        exit 20 # FAIL

    else # fields differ less than 0.001

        # remove all fields allowed to be not bit-identical
        counter=0
        for field in "${fields_with_roundoff[@]}"; do
            if [[ $counter == "0" ]]; then

                # remove all lines with uninteresting content for this checker
                # remove one field from log
                grep -v "$pattern_to_ignore" $log | grep -v "$field"  > field_step_${counter}
            else
                # remove another field from log
                grep -v "$field" field_step_${last_iteration} > field_step_${counter}
            fi
            last_iteration=$counter
            (( counter++ ))
        done

        # count lines of log -> should be only 1
        # can go wrong, if CDO outputs unexpected output
        # that is not removed earlier in the code
        should_be_one=`wc -l field_step_${last_iteration} | awk '{print $1}'`

        # field differs that is not allowed
        if [[ $should_be_one -ne 1 ]]; then
            echo "Differences from reference file found" 1>&1
            head field_step_${last_iteration} 1>&1
            exit 20 #FAIL

        else # field not bit-identical, but OK

            # count number of fields not bit-identical
            nr_of_fields=`grep -v "$pattern_to_ignore\|Date" $log | wc -l | awk '{print $1}'`

            # print info to logfile of testsuite
            echo "$nr_of_fields fields not bit-identical" 1>&1
            [[ $VERBOSE -gt 1 ]] && grep -v "$pattern_to_ignore" $log 1>&1
            rm field_step_*
            exit 10 # OK
        fi

    fi # end of fields differ more than 0.001

fi # end of fields differ
