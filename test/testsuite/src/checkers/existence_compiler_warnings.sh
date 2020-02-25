#!/bin/bash
# COSMO TECHNICAL TESTSUITE
#
# This script checks whether the compilation of Extpar caused any
# warnings for the NAG and the GCC compiler
#
# Author       Jonas Jucker 
# Maintainer   katherine.osterried@env.ethz.ch
set -x
# root dir of Extpar repo
compiledir="../../"

# warning identifier for different compilers
warning_gcc="Warning:"
warning_intel=
warning_nag=
# define compiler warnings to be checked

# define compiler warnings to be checked
to_ignore_gcc=("Wmaybe-uninitialized")
to_ignore_intel=("Wmaybe-uninitialized")
to_ignore_nag=("Wmaybe-uninitialized")

if [ ! -f "$compiledir/compile.log" ] ; then
  echo "No compilation log-file found"  1>&1
  exit 20 # FAIL
fi

# GCC 
counter=0
for ignore in "${to_ignore_gcc[@]}"; do
    if [[ $counter == "0" ]]; then
        # look for warnings, but ignore 1 of those defined above
        grep  "$warning_gcc" $compiledir/compile.log | grep -v "$ignore"  > ignore_step_${counter}
    else

        # look for warnings, but ignore 1 of those defined above
        grep "$warning_gcc" ignore_step_$last_iteration | grep -v "$ignore"> ignore_step_${counter}
    fi
    last_iteration=$counter
    (( counter++ ))
done

# look for warnings, if still some present -> exit
grep "$warning_gcc" ignore_step_$last_iteration > /dev/null
if [ $? -ne 1 ] ; then
   echo "Compiler warnings found for GCC" 1>&1
   exit 20 # FAIL
fi

# NAG 
counter=0
for ignore in "${to_ignore_nag[@]}"; do
    if [[ $counter == "0" ]]; then
        # look for warnings, but ignore 1 of those defined above
        grep  "$warning_nag" $compiledir/compile.log | grep -v "$ignore"  > ignore_step_${counter}
    else

        # look for warnings, but ignore 1 of those defined above
        grep "$warning_nag" ignore_step_$last_iteration | grep -v "$ignore"> ignore_step_${counter}
    fi
    last_iteration=$counter
    (( counter++ ))
done

# look for warnings, if still some present -> exit
grep "$warning_nag" ignore_step_$last_iteration > /dev/null
if [ $? -ne 1 ] ; then
   echo "Compiler warnings found for NAG" 1>&1
   exit 20 # FAIL
fi

# INTEL 
counter=0
for ignore in "${to_ignore_intel[@]}"; do
    if [[ $counter == "0" ]]; then
        # look for warnings, but ignore 1 of those defined above
        grep  "$warning_intel" $compiledir/compile.log | grep -v "$ignore"  > ignore_step_${counter}
    else

        # look for warnings, but ignore 1 of those defined above
        grep "$warning_intel" ignore_step_$last_iteration | grep -v "$ignore"> ignore_step_${counter}
    fi
    last_iteration=$counter
    (( counter++ ))
done

# look for warnings, if still some present -> exit
grep "$warning_intel" ignore_step_$last_iteration > /dev/null
if [ $? -ne 1 ] ; then
   echo "Compiler warnings found for INTEL" 1>&1
   exit 20 # FAIL
fi

#  cleanup
rm ignore_step_*

# goodbye
exit 0 # MATCH

