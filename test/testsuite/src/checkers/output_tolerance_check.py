#!/usr/bin/env python

"""
COSMO TECHNICAL TESTSUITE

This script checks whether the output of two runs is identical by comparing the
YUCHKDAT files.
"""

# built-in modules
import os, sys

# private modules
sys.path.append(os.path.join(os.path.dirname(__file__), "../tools"))
from ts_utilities import read_environ, dir_path, str_to_bool
from ts_fortran_nl import get_param
from sys import maxsize
import ts_yuchdat
import ts_thresholds

# information
__author__     = "Xavier Lapillonne, David Leutwyler, Santiago Moreno"
__email__      = "cosmo-wg6@cosmo.org"
__maintainer__ = "xavier.lapillonne@meteoswiss.ch"

# some global definitions
yufile   = 'YUCHKDAT'    # name of special testsuite output
yuswitch = 'lcheck'      # namelist switch controlling YUPRDBG output
lnout    = 'ngribout'    # namelist entry which specifies the nr of output lists

def check():

    # get name of myself
    myname = os.path.basename(__file__)
    header = myname+': '

    # get environment variables
    env = read_environ()
    verbose = int(env['VERBOSE'])
    rundir = dir_path(env['RUNDIR'])
    refoutdir = dir_path(env['REFOUTDIR'])
    namelistdir = dir_path(env['NAMELISTDIR'])
    switch = env['DT_FILE']
    tolerance = env['TOLERANCE']
    forcematch = int(env['FORCEMATCH']) == 1
    tune_thresholds = str_to_bool(env['TUNE_THRESHOLDS'])
    
    # defines the 2 file that belongs logically to the checker
    yufile1 = rundir + yufile
    yufile2 = refoutdir + yufile

    if not switch: # no need to extract timestep and to check output lists
        nlfile  = env['NL_TS_SWITCH']
    else:
        nlfile  = 'INPUT_IO'    # namelist file containing yuswitch
        nlfile2 = 'INPUT_ORG'   # namelist file containing dt

        # extract timestep
        try:
            dt = float(get_param(rundir+nlfile2, 'dt'))
        except:
            if verbose:
                print(header+'failed to extract dt from '+rundir+nlfile2)
            return 20 # FAIL

        # check for output lists
        try:
            nout_list = get_param(rundir+nlfile, lnout)
        except:
            if verbose:
                print(header+'no output lists in '+rundir+nlfile)
            return 20 # FAIL
    
        # check if special testsuite output was activated in every output list
        if get_param(rundir+nlfile, yuswitch, occurrence=nout_list) in ['.FALSE.', '.false.']:
            if verbose:
                print(yuswitch+' is set to .false. in '+rundir+nlfile+' for this simulation')
            return 20 # FAIL

    #check if tolerance file exists in namelistdir
    tolerance_path=namelistdir + tolerance
    if os.path.exists(tolerance_path):
        tolerance_file=tolerance_path   #in namelist dir
        ltol_file=True
    else:
        tolerance_file=''
        ltol_file=False

    # Create Thresholds object
    thres = """
 minval = 1e-12
  steps =         1        20         60        100
      * =   1.00e-7   1.00e-4   1.00e-06   1.00e-02
"""
    threshold = ts_thresholds.Thresholds(thres)
    threshold_var = '*'
    # define tolerances for comparing YUCHKDAT files
    # Use tolerance file if exists
    if ltol_file:
        if verbose==2:
            print(header + 'Using tolerance values from file')
        elif verbose>2:
            print(header + 'Using tolerance values from file '+tolerance_file)
        try:
            threshold = ts_thresholds.Thresholds(tolerance_file)
            if "CHKDAT" in threshold.variables:
                threshold_var = "CHKDAT"
        except:
            if verbose:
                print(header+'Error while reading '+tolerance_file)
                print(header+'Cannot read one of the following parameter: tol_times,tol_out,minval')
            return 20 # FAIL
    
    # Identical thresholds

    thres = """
 minval = 0.0
  steps =   0
      * =   0
"""
    threshold_identical = ts_thresholds.Thresholds(thres)

    # Override the tolerances in case FORCEMATCH is enabled
    if forcematch:
        threshold = threshold_identical

    try:
        # check for bit identical results
        if verbose>1:
            print(header + 'Checking first if results are bit identical')
        err_count_identical = ts_yuchdat.compare(yufile1, yufile2, threshold_identical, threshold_var, v_level=-1)
        
        if verbose>1:
            if err_count_identical==0:
                print(header + 'Results are bit identical')
            else:
                print(header + 'Results are not bit identical')
        
        # check if results are within tolerance values
        if verbose>1:
            print(header + 'Checking if results are within tolerance')
        
        error_count = ts_yuchdat.compare(yufile1, yufile2, threshold, threshold_var, tune_thresholds, v_level=0)
        
        if (tune_thresholds):
            threshold.to_file(tolerance_path)
        if verbose>1:
            if error_count==0:
                print(header + 'Results are within thresholds')
            else:
                print(header + 'Results are not within thresholds')

    except Exception as e:
        if verbose:
            print(e)
        return 20 # FAIL

    if err_count_identical == 0:
        return 0 # MATCH
    elif error_count == 0:
        return 10 # OK
    else:
        return 20 # FAIL

if __name__ == "__main__":
    sys.exit(check())



