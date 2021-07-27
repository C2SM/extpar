# EXTPAR

# General Information
EXTPAR (External Parameters for Numerical Weather Prediction and Climate Application) is an official software of the [COSMO Consortium](www.cosmo-model.org).  It is used to prepare the external parameter data files that are used as input for the COSMO model, and additionally now the ICON model.

The code is written in Fortran 90 and in Python. The Python scripts use CDO for the most compute-intensive parts.  Currently the code is tested regularly using the gcc, NAG, and Intel Fortran compilers.  The code is also accelerated in some places with OpenMP parallelization.  

The code once compiled generates 8 Fortran executables and 4 Python scripts, which can be run simultaneously except for the final extpar_consistency_check.exe, which is used to tie together all the external parameter results into one output file.  

In order to run Extpar, input data files for the external parameter variables are needed.  These input files are currently available on CSCS machines, DWD machines, and on mistral at DKRZ.  The input data files are also stored in a git-LFS data repository found at: https://gitlab.dkrz.de/extpar-data/extpar-input-data.  Instructions to download or update the input data files can be found in this repository.  To gain access to the git-LFS input data repository, contact the Extpar source code administrator (Jonas Jucker) at jonas.jucker@c2sm.ethz.ch.

A full documentation of code can be found in the [manual](doc/user_and_implementation_manual.pdf)

# Getting started
For compilation instructions see: [README.compile_run](doc/README.compile_run.md)

Once you have a compiled version, there are several sample run scripts available in the run_scripts folder of this repository.  You should be able to adapt one of these for your needs.  

# Testing
The extpar code comes with a technical testsuite to ensure the accuracy of the results.  More information about the testsuite can be found [here](test/testsuite/README.md)

# Information for developers
The coding rules and development workflow for Extpar can be found [here](doc/development.md)

# Support 
In the case of issues or questions, please contact the current source code administrator (Jonas Jucker) at jonas.jucker@c2sm.ethz.ch.  


