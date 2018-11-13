# README for MPI-M DYAMOND preprocessing

*Luis Kornblueh, 2018-11-13, MPI-M*

This set of scripts is for creating the input data for a global ICON simulation with NWP physics on mistral@DKRZ.
Processing input data consists of several steps:
1. interpolation and averaging of the high frequency SST and seaice data from IFS operational forecasts,
2. interpolation and averaging of the SST, seaice, and snow data from ERA interim for providing climatologies,
3.  run extpar for generating all other boundary condition data and finish with a consitency check, and
4.  create the initial data from an IFS analysis
First step is to check that all scripts are using the same grid. Next one has to take care to set the work and output directories. 
The sequence to call the single steps is:
```
dkrz-prepare-sst-seaice.sh
dkrz-monmean-sst-seaice-snow.sh
dkrz-monmean-t2m.sh
dkrz-extpar.sh
dkrz-ifs2icon.sh
```
The two available python scripts are for preparing the results of step 2 for use in the consistency check.
In case of problems contact *luis.kornblueh@mpimet.mpg.de*
