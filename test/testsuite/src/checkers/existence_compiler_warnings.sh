#!/bin/bash
# COSMO TECHNICAL TESTSUITE
#
# This script checks whether the compilation of Extpar caused any
# warnings for the NAG, INTEL and the GCC compiler
#
# Author       Jonas Jucker 
# Maintainer   katherine.osterried@env.ethz.ch

# root dir of Extpar repository
compiledir="../../"

# warning identifier for different compilers
warning_gcc="Warning:"
warning_nag="Warning\|Extension\|Questionable"
warning_intel="warning #6843"

# define compiler warnings to be ignored
to_ignore_gcc=("Wmaybe-uninitialized")
to_ignore_nag=("set but never referenced" "OpenMP" "STAT")
to_ignore_intel=("explicit")

if [ ! -f "$compiledir/compile.log" ] ; then
  echo "No compilation log-file found"  1>&1
  exit 20 # FAIL
fi

# to find out compiler the ORDER of the grep searches must ALWAYS be
# GCC, INTEL and last NAG because NAG compilation 
# also requires to load the GCC compiler

# $compiler is "GCC" for GCC and GNU
grep -i "GCC\|GNU" $compiledir/modules.env > /dev/null
if [ $? -ne 1 ] ; then
 compiler="GCC" 
fi

# INTEL
grep -i "INTEL" $compiledir/modules.env > /dev/null
if [ $? -ne 1 ] ; then
 compiler="INTEL" 
fi

# if NAG is found in modules.env, replace $compiler with NAG
grep -i "NAG" $compiledir/modules.env> /dev/null
if [ $? -ne 1 ] ; then
 compiler="NAG" 
fi

# GCC 
if [[ $compiler == "GCC" ]]; then
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
elif [[ $compiler == "NAG" ]]; then
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
elif [[ $compiler == "INTEL" ]]; then
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
    
# no compiler found in logfile
else
    echo "Could not determine compiler from compile.log" 1>&1
    exit 20 # FAIL
fi

#  cleanup
rm ignore_step_*

# goodbye
exit 0 # MATCH

