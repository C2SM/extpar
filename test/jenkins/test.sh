#!/bin/bash

# This script runs the Extpar testsuite on Daint
#


cd test/testsuite

cd data
./get_data.sh
cd ..

./submit.daint.sh

# echo output to stdout
test -f testsuite.out || exitError 1261 ${LINENO} "output file testsuite.out does not exist"
echo "=== testsuite.out BEGIN ==="
cat testsuite.out | /bin/sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"
echo "=== testsuite.out END ==="

# check result of testsuite
grep RESULT testsuite.out | egrep 'FAIL|CRASH' > /dev/null
if [ $? -eq 0 ] ; then
  exitError 1271 ${LINENO} "testsuite did not complete successfully"
fi
