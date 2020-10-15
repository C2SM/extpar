import sys
import os
import logging
import subprocess

'''
Module environment provides functions that interact with
the system Extpar is running on, it contains:

-get_cdo_version: get CDO version from environment variable

-check_environment_for_extpar: check if all required modules are loaded

-get_omp_num_threads: get value of environment variable for OpenMP
'''


def get_cdo_version(extpar_programme, host):

#    cdo_version = None
    try:
        cdo_version_text = subprocess.run([ "cdo", "--version" ], capture_output=True, encoding="utf-8")
    except KeyError:
        logging.error(f'CDO is not available on host {host}.')
        raise

    stderr_lines = cdo_version_text.stderr
    for line in stderr_lines.split("\n"):
        if "Climate Data Operators version" in line:
            (_dum1, _dum2, _dum3, _dum4, cdo_version, _dum5) = line.split()

    return cdo_version


def check_environment_for_extpar(extpar_programme):
    '''
    get hostname, python version, pythonpath and cdo version

    put all together into an info print for the logfile
    if CDO is not loaded, exit
    '''

    hostname = os.uname()[1]
    python_version = sys.version
    pythonpath = os.environ['PYTHONPATH']
    cdo_version = get_cdo_version(extpar_programme,hostname)

    logging.info('')
    logging.info('============= listen to environment ============')
    logging.info('')
    logging.info(f'Hostname -> {hostname}')
    logging.info('')
    logging.info(f'Python version -> {python_version}')
    logging.info('')
    logging.info(f'PYTHONPATH -> {pythonpath}')
    logging.info('')
    logging.info(f'CDO version -> {cdo_version}')
    logging.info('')


def get_omp_num_threads():
    '''
    get environment variables for OMP,
    if not set, assume 1 as default
    '''

    try:
        omp = os.environ['OMP_NUM_THREADS']
    except KeyError:
        omp = 1
        logging.warning('OMP_NUM_THREADS not set -> '
                        f'use OMP_NUM_THREADS = {omp} instead')
    return omp
