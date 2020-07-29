import sys
import os
import logging

'''
Module environment provides functions that interact with
the system Extpar is running on, it contains:

-get_cdo_version: get CDO version from environment variable

-check_environment_for_extpar: check if all required modules are loaded

-get_omp_num_threads: get value of environment variable for OpenMP
'''


def get_cdo_version(extpar_programme, host):
    '''
    get CDO version from environment variable EBVERSIONCDO

    If EBVERSIONCDO is not set, exit and print a message.
    For Mistral (CDO loaded per default), print a warning

    **Caution: For hosts containing an 'm' in the name, this
               can go wrong -> Extpar will then stop without
               error message at the first launch_shell()
               involving CDO.
    '''

    # get version of CDO
    try:
        import subprocess
        cdoinfo = subprocess.run(['cdo', '-V'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cdo_version=cdoinfo.stderr.split()[4].decode("utf-8")
    except KeyError:
        hostname = os.uname()[1]

        # host is not mistral
        if 'm' not in hostname: 
            logging.error(f'CDO is not loaded on host {host}. Please '
                          f'load CDO before you run {extpar_programme}')
            sys.exit(1)

        # host is mistral
        else:
            logging.warning('Assume that host is mistral@dkrz.de '
                            '-> CDO is already loaded as default')

            cdo_version = 'mistral default'

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
