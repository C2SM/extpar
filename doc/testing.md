# Testing
Extpar is tested with an adapted version of the [COSMO technical testsuite](https://github.com/C2SM-RCM/testsuite).

## Run Tests
First step is to compile the code following the instructions in [compile_run](compile_run.md).

#### Daint

```
cp bin/* test/testsuite/bin/.
cd test/testsuite
./data/get_data.sh
./bin/extract_inputfiles_from_namelist.py
sbatch --wait ./submit.daint.transfer.sh
sbatch submit.daint.sh
```

#### Tsa
```
cp bin/* test/testsuite/bin/.
cd test/testsuite
./data/get_data.sh
./bin/extract_inputfiles_from_namelist.py
sbatch --wait ./submit.tsa.transfer.sh
sbatch submit.tsa.sh
```

#### Mistral

```
cp bin/* test/testsuite/bin/.
cd test/testsuite
./data/get_data.sh
sbatch submit.mistral.sh # or submit.mistral.intel.sh
```

The results of the testsuite can be found in file _testsuite.out_

An example output could look as follows:  
<img width="400" alt="Testsuite_log" src="https://user-images.githubusercontent.com/39263956/132489686-1987e5f7-3d26-4994-8af2-ae459a908b44.png">

## Testlists
There are many different testlist, each containing a set of tests for different setups, compiler or models:

#### GCC
* [COSMO](../test/testsuite/testlist_cosmo.xml)
* [ICON](../test/testsuite/testlist_icon.xml)
* [Landuse](../test/testsuite/testlist_landuse.xml)

#### Intel
* [COSMO](../test/testsuite/testlist_cosmo_intel.xml)
* [ICON](../test/testsuite/testlist_icon_intel.xml)
 
 
 ## Tolerances
 It is possible to define an optional tolerance threshold for each test and each field.
 To allow deviations for the test _icon_d2_ for example, just dit the the [tolerances-file](../test/testsuite/data/dwd/icon_d2/tolerances).
 The syntax is as follows:
 ```
 PARAMETER, abs_diff
 NDVI, 9.0e-08
 W_SNOW, 5.0e-05
 ```
 
 ## Add a new test

  1. Modify the testlist.xml file to add the new test.  Alternatively, you could also add a new 
     testlist XML file (with a new name).  The testlist which is run can be chosen from the testsuite
     command line.  

  2. Make a folder in the data folder for the new test containing the INPUT_* files and the namelist.py for the Python-CDO modules.  

  3. Send any binary reference files to upload to the ftp site to the source code administrator. 
