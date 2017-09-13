#!/bin/bash

# This script runs the Extpar testsuite on Daint
#

# Error function
exitError()
{
        echo "ERROR $1: $2" 1>&2
        echo "ERROR     LOCATION=$0" 1>&2
        exit $1
}


# Code body

cd test/testsuite

# Get the input data
cd data
./get_data.sh
cd ..

# Copy the executables
cp ../../bin/*.exe bin

test -f submit.daint.sh || exitError 1260 "submit script submit.daint.sh does not exist" 

echo "Running submit script"
./submit.daint.sh
echo "Finished with submit script"


# echo output to stdout
test -f testsuite.out || exitError 1261 "output file testsuite.out does not exist"
echo "=== testsuite.out BEGIN ==="
cat testsuite.out | /bin/sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"
echo "=== testsuite.out END ==="

# check result of testsuite
grep RESULT testsuite.out | egrep 'FAIL|CRASH' > /dev/null
if [ $? -eq 0 ] ; then
  exitError 1271 "testsuite did not complete successfully"
fi
