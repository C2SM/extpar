#!/usr/bin/env python

'''
Simple tolerance checker for EXTPAR

Parse the output of cdo diffv to see, whether differences of fields
are still within tolerances.

User-defined tolerances can be passed for each field with the
tolerance file.
It requires the following structure:

    all_ABS_DIFF: 1.0e-14

    PARAMETER, abs_diff
    T2M,1.0e-7
    SO2,1.0e-12

all_ABS_DIFF defines the maximal absolute difference allowed for ALL fields.
Parameters requiring higher tolerance thresholds can be defined below.
all_ABS_DIFF nor single parameters need to be defined, default behaviour is not allowing
differences from the reference.
It is possible to only define all_ABS_DIFF or only define PARAMETERS or both.

Current limitations:
    - In case a tolerance file is present for a test, the test cannot detect,
      whether the results are bit-identical or within the tolerances.
      In that case the result will always be "OK".
      The checker diff_ncdf_out.sh should be used in that case.


# Author       Jonas Jucker
# Maintainer   jonas.jucker@env.ethz.ch

'''

import os
import glob
import sys
import subprocess

def filter_string_list(string_list, filter_strings):
    '''
    Keep all entries of a list of strings, that do not
    contain an entry of filter_strings
    '''
    filtered_string = []
    for line in string_list:
        if not any(filter_string in line for filter_string in filter_strings):
            filtered_string.append(line)

    return filtered_string

def read_tolerances_from_text_file(tolerance_file):
    '''
    Infer the all_ABS_DIFF parameter and the tolerances per
    parameter by reading the tolerance-file
    '''
    tolerance_per_param = {}
    all_abs_diff = 0.0
    with open(tolerance_file) as f:
        for line in f:
            if 'PARAMETER' not in line:
                line = line.strip('\n')
                if 'all_ABS_DIFF' in line:
                    
                    all_abs_diff = float(line.split()[-1])
                else:
                    split = line.split(',')
                    try:
                        param = split[0]
                        abs_diff = float(split[1])
                    except IndexError:
                        param = 'parsing failed'
                        abs_diff = 0.0

                    absval_dict = {'abs_diff': abs_diff}
                    tolerance_per_param[param] = absval_dict

    return tolerance_per_param, all_abs_diff

def get_cdo_version():
    '''
    get CDO version from `cdo -V`                                                                                                                                                                                                 
    '''

    cdo_cmd = ['cdo', '-V']
    try:

        process = subprocess.run(cdo_cmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=True,
            universal_newlines=True)

        output = process.stdout + process.stderr

    except FileNotFoundError:
        print('No CDO available')
        sys.exit(20)

    cdo_version = output.split()[4]

    return cdo_version

def versiontuple(v):
    return tuple(map(int, (v.split("."))))


'''
Start main script
'''

# important environment variables
try:
    verbose = int(os.environ['TS_VERBOSE'])
    rundir = os.environ['TS_RUNDIR']
    refoutdir = os.environ['TS_REFOUTDIR']
except KeyError:
    print('Testsuite environment variables not set')
    sys.exit(20)

cdo_version = get_cdo_version()    

# get all netCDF present in refoutdir
ref_files = glob.glob('{}/external_parameter*.nc'.format(refoutdir))

if not ref_files:
    print('netCDF reference files found')
    sys.exit(20)
else:
    ref_file = ref_files[0]

# read tolerance file if present
tolerance_file = '{}/tolerances'.format(refoutdir)
if os.path.isfile(tolerance_file):
    tolerance_dict, cdo_abs_diff = read_tolerances_from_text_file(tolerance_file)

# assume zero tolerances for all fields
else:
    cdo_abs_diff = 0.0
    tolerance_dict = {}

# compare references with test output
file_to_test = glob.glob('{}/external_parameter*.nc'.format(rundir))[0]

if not os.path.isfile(file_to_test):
  print('No netCDF output found')
  sys.exit(20)
else:
    # compare fields using CDO
    if versiontuple(cdo_version) > versiontuple('1.9.5'):                                                                                                                                                                     
        diffv_cmd = 'diffv,abslim={}'.format(cdo_abs_diff) 
    else:
        diffv_cmd = 'diffv,{}'.format(cdo_abs_diff)
    cdo_cmd = ['cdo', '--sortname', diffv_cmd, ref_file, file_to_test]

try:
    process = subprocess.run(cdo_cmd, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=True,
        universal_newlines=True)

    output = process.stdout + process.stderr

except FileNotFoundError:
    print('No CDO available')
    sys.exit(20)

except subprocess.CalledProcessError as e:
    output = e.stdout + e.stderr
    crash_indicators = ['Abort', 'failed']                                                                                                                                                                                    
    if any(x in output for x in crash_indicators): 
        print('CDO abort during file comparison')
        sys.exit(20)

test_fail = False
test_ok = False

# there are differences greater that cdo_abs_diff
if 'differ' in output:

    bad_words = ['Warning','diffn','differ','Parameter','nhori']
    clean_output = filter_string_list(output.split('\n'), bad_words)

    # false alarm -> remaining lines not related to field differences
    if not clean_output:
        if verbose > 0:
            print('Differences OK')
        test_ok = True

    # some fields differ
    else:
        failed_parameters = []

        # check line by line of remaining CDO output
        for line in clean_output:
            if line:
                split = line.split()
                parameter = split[-1]
                abs_diff = split[11]
                diff_count = split[7]
                total_count = split[5]
                info_string = '{} {} {}/{}'.format(parameter, abs_diff, diff_count, total_count)

                # parameter has non-default tolerance threshold
                try:
                    parameter_tolerances = tolerance_dict[parameter]

                    if parameter_tolerances['abs_diff'] < float(abs_diff):
                        failed_parameters.append(info_string)

                # parameter has default tolerance threshold
                except KeyError:
                        failed_parameters.append(info_string)

        # all parameters have differences within thresholds
        if not failed_parameters:
            if verbose > 0:
                print('Differences OK')
            test_ok = True

        # info prints
        else:
            counter = 0
            if verbose > 1:
                print('  PARAM ABSDIFF DIFF/TOTAL')
                for entry in failed_parameters:
                    if counter < 50:
                        print('  ' + entry)
                    counter += 1

            if verbose > 0:
                print('Differences from reference found')
            test_fail = True

# no differences greater then value defined in tolerances-file
else:
    if verbose > 0:
        if (cdo_abs_diff <= 1.0e-30) and not bool(tolerance_dict):
            print('Bit-identical')
        else:
            print('Differences OK')

if test_fail:
    sys.exit(20)

elif (test_ok or cdo_abs_diff >= 1.0e-30 or bool(tolerance_dict)):
    sys.exit(10)

else:
    sys.exit(0)
