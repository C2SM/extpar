#!/usr/bin/env python

"""
COSMO TECHNICAL TESTSUITE

This script checks whether two runs are bit-identical by comparing the
YUPRTEST files.
"""

# built-in modules
import os, sys

# private modules
sys.path.append(os.path.join(os.path.dirname(__file__), "../tools"))
from ts_utilities import read_environ, dir_path, str_to_bool
from ts_fortran_nl import get_param
from ts_yuprtest import *
from sys import maxsize
import traceback

# information
__author__      = "Santiago Moreno"
__email__       = "cosmo-wg6@cosmo.org"
__maintainer__  = "xavier.lapillonne@meteoswiss.ch"

# some global definitions
yufile = 'YUPRTEST'     # name of special testsuite output
yuswitch = 'ltestsuite' # namelist switch controlling YUPRTEST output

def check():
    # get name of myself
    myname = os.path.basename(__file__)
    header = myname+': '

    # get environment variables
    env = read_environ()
    verbosity = int(env['VERBOSE'])
    rundir = dir_path(env['RUNDIR'])
    refoutdir = dir_path(env['REFOUTDIR'])
    namelistdir = dir_path(env['NAMELISTDIR'])
    tolerance = env['TOLERANCE']
    switch = env['NL_TS_SWITCH']
    forcematch = int(env['FORCEMATCH']) == 1
    tune_thresholds = str_to_bool(env['TUNE_THRESHOLDS'])
    tune_times = int(env['TUNING_ITERATIONS'])
    reset_thresholds = str_to_bool(env['RESET_THRESHOLDS'])

    #check if namelist file with switch exists in namelistdir
    switch_path = namelistdir + switch
    if not os.path.exists(switch_path):
        if verbosity>0:
            print header + "unable to find namelist file with switch in " + switch_path
        return 20 # FAIL

    # defines the 1 file that belongs logically to the checker
    yufile1 = rundir + yufile
    yufile2 = refoutdir + yufile
    # check if special testsuite output was activated
    if get_param(rundir+switch, yuswitch) in ['.FALSE.', '.false.']:
        if verbosity:
            print yuswitch +' is set to .false. in '+ rundir + switch +' for this simulation'
        return 20 # FAIL

    #check if tolerance file exists in namelistdir or type dir
    tolerance_path = namelistdir + tolerance

    # The forcematch option gives precedence over the tolerance file
    if forcematch:
        thresh = """
        minval =    0.0
        steps =   {maxsteps}
            * =     0.0""".format(maxsteps=maxsize)
    elif os.path.exists(tolerance_path):
        thresh=tolerance_path   #in namelist dir
    else:
        if verbosity>0:
            print header + "unable to find tolerance file at " + tolerance
        return 20 # FAIL

    try:
        c = Compare(yufile1, yufile2, thresh)
        if reset_thresholds:
            c.reset_thresholds()
        if tune_thresholds:
            c.update_thresholds()
            c.write_threshold_to_file(tolerance_path)
        result = c.compare_data()

        if (verbosity > 1) or (result >= 2):
            print(header)
            c.print_results()
        elif verbosity > 1:
            print(header)

    except Exception as e :
        traceback.print_exc(file=sys.stdout)
        print(header + str(e))
        return 30 # CRASH
    if (result == 0):
        if verbosity>1:
            print header + "Results are within the thresholds and bit identical"
        return 0 # MATCH
    if (result == 1):
        if verbosity>1:
            print header + "Results are within thresholds, but are not bit identical"
        return 10 # OK
    if (result == 2):
        if verbosity>1:
            print header + "Some or all Results are not within thresholds"
        return 20 # FAIL


if __name__ == "__main__":
    sys.exit(check())
