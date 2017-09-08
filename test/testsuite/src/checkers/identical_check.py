#!/usr/bin/env python2

"""
COSMO TECHNICAL TESTSUITE

This script checks whether two runs are bit-identical by comparing the
YUPRTEST files.
"""

# built-in modules
import os, re, sys

# private modules
sys.path.append(os.path.join(os.path.dirname(__file__), "../tools")) # this is the generic folder for subroutines
from ts_utilities import read_environ, dir_path, str_to_bool
from ts_fortran_nl import get_param
from sys import maxsize
import comp_yuprtest

# information
__author__     = "Nicolo Lardelli, Oliver Fuhrer, Santiago Moreno"
__email__      = "cosmo-wg6@cosmo.org"
__maintainer__ = "xavier.lapillonne@meteoswiss.ch"

# some global definitions
yufile = 'YUPRTEST'   # name of special testsuite output

def check():

    # get name of myself
    myname = os.path.basename(__file__)
    header = myname+': '

    # get environment variables
    env = read_environ()
    verbose = int(env['VERBOSE'])
    rundir = dir_path(env['RUNDIR'])
    refoutdir = dir_path(env['REFOUTDIR'])
    switch = env['DT_FILE']
    tune_thresholds = env['TUNE_THRESHOLDS']
    if(str_to_bool(tune_thresholds)):
        print ("WARNING: Skipped as identical checker as it isn't usable while tuning thresholds")
        return 15 #SKIP

    # defines the 2 file that belongs logically to the checker
    yufile1 = rundir + yufile
    yufile2 = refoutdir + yufile

    if switch: # extract timestep
        nlfile = switch  # namelist file containing dt
        try:
            dt = float(get_param(rundir+nlfile, 'dt'))
        except:
            if verbose:
                print header+'failed to extract dt from '+rundir+nlfile
            return 20 # FAIL
    
    # define tolerances for comparing YUPRTEST files
    nts = [maxsize]
    tols = [0.0]
    tas = [0.0]

    try:
        # check for bit identical results
        if verbose>1:
            print header + 'Checking if results are bit identical'
        if verbose>2:
            print header + 'comp_yuprtest()'
            print header + '  file1 = '+yufile1
            print header + '  file2 = '+yufile2
            print header + '  minval = '+str(-1)
            print header + '  nts = [%s]' % ','.join(map(str, nts))
            print header + '  tols = [%s]' % ','.join(map(str, tols))
            print header + '  tas = [%s]' % ','.join(map(str, tas))
        error_count = comp_yuprtest.cmp_(yufile1, yufile2, 0, -1, nts, tols, tas)
        if verbose>1:
            if error_count==0:
                print header + 'Results are bit identical'
            else:
                print header + 'Results are not bit identical'

    except RuntimeError as e:
        if verbose:
            print e
        return 20 # FAIL

    if error_count:
        return 20 # FAIL
    else:
        return 0 # MATCH

if __name__ == "__main__":
    sys.exit(check())

