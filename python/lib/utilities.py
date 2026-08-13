import logging
import os
import subprocess

try:
    from extpar.lib.fortran_namelist import read_variable
except ImportError:  # package not installed -> use PYTHONPATH
    from fortran_namelist import read_variable
'''
Module utilities provides a bunch of helpful functions for Extpar,
it contains:

-launch_shell: wrapper to launch a programm in the shell

-remove: remove a file from the system

-clean_path: make a clean path from a dir and file

-check_gridtype : check whether itype_grid from namelist is correct

-check_itype_cru : check whether itype_cru from namelist is correct

-check_eratype: check whether iera_type from namelist is correct

-check_albtype: check whether ialb_type from namelist is correct

-check_ahftype: check whether iahf_type from namelist is correct

-check_isatype: check whether isa_type from namelist is correct

-check_emisstype: check whether emiss_type from namelist is correct

-determine_albedo_varnames: assign correct varnames for different ialb_type
'''


def launch_shell(bin, *args):
    '''
    wrapper to launch an external programme on the system

    bin is the executable to run
    *args are the arguments for bin, need to be convertable to string
    stdout/stderr of bin is written to logfile
    '''

    #convert *args to string
    arg_list = []
    arg_list.insert(0, str(bin))
    for arg in args:
        if arg:  # Prevents empty strings from being written into list
            arg_list.append(str(arg))

    args_for_logger = ' '.join(arg_list)

    logging.info(f'Launch shell command: {args_for_logger}')
    logging.info('')

    try:
        process = subprocess.run(arg_list,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 check=True,
                                 universal_newlines=True)
        output = process.stdout + process.stderr
    except FileNotFoundError:
        logging.warning(f'Problems with shell command: {args_for_logger} \n'
                        '-> it appears your shell does not know this command')

        logging.error('Shell command failed', exc_info=True)
        raise

    except subprocess.CalledProcessError as e:
        output = e.stdout + e.stderr
        logging.warning(f'Problems with shell command: {args_for_logger} \n'
                        '-> the output returned to the shell is:')
        logging.warning(f'{output}')

        logging.error('Shell command failed', exc_info=True)
        raise

    logging.info('Output:')
    logging.info(f'{output}')

    return output


def remove(file):
    '''
    remove file from system if exists
    '''

    if os.path.exists(file):
        os.remove(file)
    else:
        logging.warning(f'Could not remove {file}')

    return


def clean_path(dir, file):
    '''
    returns a clean path from a dir and file

    used to check if all files from namelist
    exist, if not exit programme
    '''

    clean_path = os.path.join(dir, file)
    try:
        f = open(clean_path)

    except FileNotFoundError:
        logging.error('File not found', exc_info=True)
        raise

    f.close()

    return clean_path


def check_eratype(era_type):
    '''
    check era_type for correctness and return value,
    if not exit programme
    '''

    if (era_type > 2 or era_type < 1):
        logging.error(f'iera_type {era_type} does not exist')
        raise ValueError(f'iera_type {era_type} does not exist')

    if (era_type == 1):
        logging.info('process ERA5 data')

    if (era_type == 2):
        logging.info('process ERA-I data')

    return era_type


def check_albtype(alb_type):
    '''
    check alb_type for correctnes and return value,
    if not exit programme
    '''

    if (alb_type > 3 or alb_type < 1):
        logging.error(f'ialb_type {alb_type} does not exist.')
        raise ValueError(f'ialb_type {alb_type} does not exist.')

    if (alb_type == 1):
        logging.info('process albedo data  for VIS, NIR and UV spectra')

    if (alb_type == 2):
        logging.info('process albedo data  for dry and saturated soil')

    if (alb_type == 3):
        logging.info('process albedo data  for VIS only')

    return alb_type


def check_ahftype(ahf_type):
    '''
    check ahf_type for correctnes and return value,
    if not exit programme
    '''

    if (ahf_type > 2 or ahf_type < 1):
        logging.error(f'iahf_type {ahf_type} does not exist.')
        raise ValueError(f'iahf_type {ahf_type} does not exist.')

    if (ahf_type == 1):
        logging.info('process ahf data with spatial resolution of 2.5 min')

    if (ahf_type == 2):
        logging.info('process ahf data with spatial resolution of 30 sec')

    return ahf_type


def check_isatype(isa_type):
    '''
    check isa_type for correctnes and return value,
    if not exit programme
    '''

    if (isa_type > 2 or isa_type < 1):
        logging.error(f'isa_type {isa_type} does not exist.')
        raise ValueError(f'isa_type {isa_type} does not exist.')

    if (isa_type == 1):
        logging.info('process isa data with spatial resolution of 30sec')

    if (isa_type == 2):
        logging.info('process isa data with spatial resolution of 10sec')

    return isa_type


def check_aottype(aot_type):
    '''
    check aot_type for correctnes and return value,
    if not exit programme
    '''

    if (aot_type > 2 or aot_type < 1):
        logging.error(f'aot_type {aot_type} does not exist.')
        raise ValueError(f'aot_type {aot_type} does not exist.')

    if (aot_type == 1):
        logging.info('process aot data with spatial resolution of 30sec')

    if (aot_type == 2):
        logging.info('process aot data with spatial resolution of 10sec')

    return aot_type


def check_emisstype(emiss_type):
    '''
    check emiss_type for correctness and return value,
    if not exit programme
    '''
    if (emiss_type < 1 or emiss_type > 2):
        logging.error(f'iemiss_type {emiss_type} does not exist. '
                      'Use 1 (full-range) or 2 (longwave) instead!')

    if (emiss_type == 1):
        logging.info('process full-range emissivity data')

    if (emiss_type == 2):
        logging.info('process long-wave emissivity data only')

    return emiss_type


def check_cdnctype(cdnc_type):
    '''
    check cdnc_type for correctness and return value,
    if not exit programme
    '''

    if (cdnc_type == 1):
        logging.info('process cdnc data from sampling method Q06')
    elif (cdnc_type == 2):
        logging.info('process cdnc data from sampling method G18')
    elif (cdnc_type == 3):
        logging.info('process cdnc data from sampling method BR17')
    else:
        logging.error(f'cdnc_type {cdnc_type} does not exist.')
        raise ValueError(f'cdnc_type {cdnc_type} does not exist.')

    return cdnc_type


def check_gridtype(input_grid_org):
    '''
    check gridtype read from input_grid_org

    read gridtype and domain_def_namelist from Fortran namelist
    if gridtype is valid, return gridtype and name of domain_def_namelist
    exit if gridtype is wrong or any namelist does not exist
    '''

    grid_org = clean_path('', input_grid_org)

    grid_type = read_variable(grid_org, 'igrid_type', int)

    def_domain_namelist = read_variable(grid_org, 'domain_def_namelist', str)

    grid_fortran_namelist = clean_path('', def_domain_namelist)

    if (grid_type < 1 or grid_type > 2):
        logging.error(f'grid_type {grid_type} does not exist. '
                      f'Use 1 (Icon) or 2 (Cosmo) instead!')
        raise ValueError(f'grid_type {grid_type} does not exist. '
                         f'Use 1 (Icon) or 2 (Cosmo) instead!')

    return grid_type, grid_fortran_namelist


def check_itype_cru(itype_cru):
    '''
    check itype_cru for correctness and return value,
    if not exit programme
    '''

    if (itype_cru > 2 or itype_cru < 1):
        logging.error(f'itype_cru {itype_cru} does not exist. '
                      f'Use 1 (fine) or 2 (coarse and fine) instead!')

        raise ValueError(f'itype_cru {itype_cru} does not exist. '
                         f'Use 1 (fine) or 2 (coarse and fine) instead!')

    if (itype_cru == 1):
        logging.info('Process fine resolution for land')

    if (itype_cru == 2):
        logging.info('Process fine resolution for land, '
                     'coarse resolution for sea')

    return itype_cru


def check_tcorr_lapse_rate(tcorr_lapse_rate):
    '''
    check tcorr_lapse_rate for correctness and return value,
    if not exit programme
    '''

    if (tcorr_lapse_rate > 0.0098 or tcorr_lapse_rate < 0.005):
        logging.error(f'Invalid value entered for temperature lapse rate. '
                      f'Use a value between 0.005 and 0.0098 K/m instead!')

        raise ValueError(f'Invalid value entered for temperature lapse rate. '
                         f'Use a value between 0.005 and 0.0098 K/m instead!')

    logging.info(f'Temperate lapse rate is set to {tcorr_lapse_rate} K/m')

    return tcorr_lapse_rate


def check_tcorr_offset(tcorr_offset):
    '''
    check tcorr_offset for correctness and return value,
    if not exit programme
    '''

    if (tcorr_offset > 10.0 or tcorr_offset < -10.0):
        logging.error(f'Invalid value entered for temperature offset. '
                      f'Use a value between -10.0 and +10.0 K instead!')

        raise ValueError(f'Invalid value entered for temperature offset. '
                         f'Use a value between -10.0 and +10.0 K instead!')

    logging.info(f'Chosen temperature offset: {tcorr_offset} K')

    return tcorr_offset


def determine_emiss_varnames(iemiss_type):
    '''
    assign the correct variable names for different iemiss_type
    '''

    if (iemiss_type == 1):
        var = 'bbemis_full'

    if (iemiss_type == 2):
        var = 'bbemis_longwave'

    return var


def determine_albedo_varnames(ialb_type):
    '''
    assign the correct variable names for different ialb_type
    '''

    if (ialb_type == 1):
        var_1 = 'al'
        var_2 = 'alnid'
        var_3 = 'aluvd'

    elif (ialb_type == 2):
        var_1 = 'ALB_SAT'
        var_2 = 'ALB_DRY'
        var_3 = ''

    elif (ialb_type == 3):
        var_1 = 'al'
        var_2 = ''
        var_3 = ''

    return var_1, var_2, var_3
