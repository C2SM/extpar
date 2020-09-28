#!/bin/bash

# COSMO TECHNICAL TESTSUITE
#
# This script checks whether the Extpar Python-code violates
# the Pep8 coding standards
#
# Author       Jonas Jucker
# Maintainer   katherine.osterried@env.ethz.ch

path_to_lib="../../python/lib"
path_to_buffer_scripts="../../python"

# check if necessary scripts exist
if [ ! -f "$path_to_buffer_scripts/pep8_checker.sh" ] ; then
  echo "Script pep8_checker.sh not found" 1>&1
  exit 20 # FAIL
fi

if [ ! -f "$path_to_lib/pycodestyle.py" ] ; then
  echo "Script pycodestyle not found" 1>&1
  exit 20 # FAIL
fi

# count number of py-files
py_files_in_lib=`find $path_to_lib/*.py | wc -l`
py_buffer_files=`find $path_to_buffer_scripts/*.py | wc -l`

if [ $py_files_in_lib -lt 1 ] ; then
  echo "No Python files found in directory $path_to_lib" 1>&1
  exit 20 # FAIL
fi

if [ $py_buffer_files -lt 1 ] ; then
  echo "No Python files found in directory $path_to_buffer_scripts" 1>&1
  exit 20 # FAIL
fi

# run Pep8 test
pep8_violations=`${path_to_buffer_scripts}/pep8_checker.sh ${path_to_lib}`

# no violations found
if [ ! -z "${pep8_violations}" ]; then
    echo "Violations of Pep8-standard found!" 1>&1
    exit 20 # FAIL
fi

# goodbye
exit 0 # MATCH
