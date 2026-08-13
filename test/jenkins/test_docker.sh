MODE=${1:-'NORMAL'}

cd /workspace/test/testsuite

# Get the input data
cd data
./get_data.sh
cd ..

# Copy the executables
cp ../../bin/* bin

./submit.docker.sh 2>&1 | tee -a testsuite.out

if [[ ${MODE} == "DEBUG" ]]; then
    cp testsuite.out /artifacts/.
    for dir in work/*/*/; do
        subdir=${dir#"work/"}
        mkdir -p /artifacts/$subdir
        cp $dir/*.log /artifacts/$subdir/.
        cp $dir/INPUT_* /artifacts/$subdir/.
        cp $dir/TEST_RES /artifacts/$subdir/.
    done
fi

grep RESULT testsuite.out | egrep 'FAIL|CRASH' > /dev/null
if [ $? -eq 0 ] ; then
    echo "testsuite did not complete successfully"
    exit 1
fi
