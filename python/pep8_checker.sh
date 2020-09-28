#!/bin/ksh

# Scripts to check if Python files allign with the Pep8 Coding standard
# All files in the "python" and "python/lib" are checked

# The scripts has two use-cases:

#   1. ./check_PEP8_standard.sh to see if your Python code in
#      the directories "python" and "python/lib" 
#      alligns with the Pep8 standard

#   2. ./check_PEP8_standard.sh $PYTHONPATH to check the Pep8 standard 
#      from a checker during execution of testsuite.py 

# Not all coding rules are required du to simple practical reasons
# or (E722) because no other way was possible in the code.

# The following picky coding rules are ignored:

############# ERROR CODES ####################

# E265: block comment should start with '#'
# E201: whitespace after '('
# E221: multiple spaces before operator
# E231: missing whitespace after ',' ';' or ':'
# E241: multiples spaces after ','
# E722: do not use bare except, specify exception instead

############# WARNINGS #######################

# W291: trailing whitespace
# W504: line break after binary operator

# check arguments

# assume that argument is the path to python/lib

# use case 2:
if [ ${#} -gt 0 ]; then
    path_to_lib=$1
    path_to_buffer_scripts=$path_to_lib/..

# use case 1:
else
    path_to_lib="lib/"
    path_to_buffer_scripts="./"
fi

# python directory
$path_to_lib/pycodestyle.py --ignore=E265,W291,E201,E221,E231,E241,E722,W504 $path_to_buffer_scripts/*.py

# python/lib directory
$path_to_lib/pycodestyle.py --exclude=namelist* --ignore=E265,W291,E201,E221,E231,E241,E722,W504 $path_to_lib/*.py
