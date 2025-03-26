cd /workspace/test/testsuite

# Get the input data
cd data
./get_data.sh
cd ..

# Copy the executables
cp ../../bin/* bin

./submit.docker.sh 2>&1 | tee -a testsuite.out

grep RESULT testsuite.out | egrep 'FAIL|CRASH' > /dev/null
if [ $? -eq 0 ] ; then
   echo "testsuite did not complete successfully"
   exit 1

fi
