# Testing

EXTPAR is tested with an adapted version of the [COSMO technical testsuite :material-open-in-new:](https://github.com/C2SM-RCM/testsuite){:target="_blank"}.

## Run Tests

First step is to compile the code following the instructions in [Compile and Run](compile_run.md).

### Docker

```bash
docker run extpar bash -c "/workspace/test/jenkins/test_docker.sh"
```

### Levante

```bash
cp bin/* test/testsuite/bin/.
cd test/testsuite
./data/get_data.sh
sbatch submit.levante.sh
```

The results of the testsuite can be found in file _testsuite.out_

An example output could look as follows:  
<img width="400" alt="Testsuite_log" src="https://user-images.githubusercontent.com/39263956/132489686-1987e5f7-3d26-4994-8af2-ae459a908b44.png">

## Testlists

There are many different testlist, each containing a set of tests for different setups, compiler or models:

### GCC

* [COSMO :material-open-in-new:](https://github.com/C2SM/extpar/tree/master/test/testsuite/testlist_cosmo.xml){:target="_blank"}
* [ICON :material-open-in-new:](https://github.com/C2SM/extpar/tree/master/test/testsuite/testlist_icon.xml){:target="_blank"}
* [Landuse :material-open-in-new:](https://github.com/C2SM/extpar/tree/master/test/testsuite/testlist_landuse.xml){:target="_blank"}

## Tolerances

It is possible to define an optional tolerance threshold for each test and each field.
To allow deviations for the test _icon_d2_ for example, just dit the the [tolerances file :material-open-in-new:](https://github.com/C2SM/extpar/tree/master/test/testsuite/data/dwd/icon_d2/tolerances){:target="_blank"}.

The syntax is as follows:

```
PARAMETER, abs_diff
NDVI, 9.0e-08
W_SNOW, 5.0e-05
```
 
## Add a New Test

 1. Modify the testlist.xml file to add the new test.  Alternatively, you could also add a new 
    testlist XML file (with a new name).  The testlist which is run can be chosen from the testsuite
    command line.  

 2. Make a folder in the data folder for the new test containing the `INPUT_*` files and the namelist.py for the Python-CDO modules.  

 3. Send any binary reference files to upload to the ftp site to the source code administrator. 
