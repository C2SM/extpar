cd /workspace/test/testsuite

# Get the input data
cd data
./get_data.sh
cd ..

# Copy the executables
cp ../../bin/* bin

./submit.docker.sh
