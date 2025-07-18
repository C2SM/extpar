cd /workspace/test/testsuite

# Get the input data
cd data
./get_data.sh
cd ..

# Copy the executables
cp ../../bin/* bin

./submit.docker.sh ${1:-"2"} 2>&1 | tee -a testsuite.out

cp testsuite.out /net/co2/c2sm-services/extpar/test/.
for dir in work/*/; do
   subdir=$(basename $dir)
   mkdir -p /net/co2/c2sm-services/extpar/test/$subdir
   cp $dir/*.log /net/co2/c2sm-services/extpar/test/$subdir/.
   cp $dir/INPUT_* /net/co2/c2sm-services/extpar/test/$subdir/.
   cp $dir/TEST_RES /net/co2/c2sm-services/extpar/test/$subdir/.
done

grep RESULT testsuite.out | egrep 'FAIL|CRASH' > /dev/null
if [ $? -eq 0 ] ; then
   echo "testsuite did not complete successfully"
   exit 1
fi
