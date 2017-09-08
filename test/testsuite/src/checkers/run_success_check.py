#!/usr/bin/env python

"""
COSMO TECHNICAL TESTSUITE

This script checks whether the standard output of the COSMO simulation
contains the given patterns which are indicative of a fully successful
simulation

New patterns can be defined directly in the script
the tuple patternlist contains all the patterns.
"""
# built-in modules
import os, sys
# private modules
sys.path.append(os.path.join(os.path.dirname(__file__), "../tools"))  # this is the generic folder for subroutines
from ts_utilities import read_environ, dir_path
from filechecker import *

# information
__author__      = "Santiago Moreno, Nicolo Lardelli, Oliver Fuhrer"
__email__       = "cosmo-wg6@cosmo.org"
__maintainer__  = "xavier.lapillonne@meteoswiss.ch"

def run_checker():
    # get environment variables
    env = read_environ()
    verbose = int(env['VERBOSE'])
    rundir = env['RUNDIR']
    log_output = env['LOGFILE']
    # construct stdout filename
    working_dir = dir_path(rundir).replace("./", "", 1) 
    logfile = os.path.join(working_dir, log_output)

    cosmo_patterns = [
    #   Class/Type                  Name                    RegularExpression
        WarningPattern(             "CFL pattern",          "CFL"                       ),
        OccurrenceCrashPattern(     "Cleanup pattern",      "(.*)^(.*)CLEAN(\s*)UP(.*)" )
    ]
    cosmo_filechecker = FileChecker()
    cosmo_filechecker.add_pattern_list(cosmo_patterns)
    return cosmo_filechecker.check(logfile, verbose)

if __name__ == "__main__":
    sys.exit(run_checker())
