import logging
import sys

try:
    import extpar.lib.utilities as utils
except ImportError:  # package not installed -> use PYTHONPATH
    import utilities as utils

'''
Module providing function and classes needed for writing/reading
Fortran-namelists with the python version of Extpar,
it contains:

    -read_variable: read a variable from given Fortran namelist

    -write_fortran_namelist: write a Fortran namelist using corr. class

    -classes defining the groups and its variables of a Fortran namelist
        -InputTclim
        -InputAlb
        -InputEmiss
        -InputNdvi
        -InputEra
        -InputAhf
        -InputIsa
'''


def read_variable(namelist, variable, type_to_convert):
    '''
    read variable from existing Fortran namelist

    namelist is read line by line, the first occurence of variable
    is taken as the value and returned with type int, float or string
    '''

    namelist = utils.clean_path('', namelist)

    with open(namelist, 'r') as f:

        # read line by line
        for line in f:
            line = line.rstrip().lstrip()

            # line is commented
            if line.startswith("!"):
                logging.warning('Ignore commented line: '
                                f'{line}')

            # valid entry in namelist
            else:

                if variable in line:
                    split = line.split('=')

                    # return last element of split
                    raw_variable = split[-1].strip()

                    characters_to_strip = ["'", ",", '"']
                    for character in characters_to_strip:
                        raw_variable = raw_variable.strip(character)

                    # convert string to type_to_convert

                    # string
                    if (type_to_convert == str):
                        try:
                            converted_variable = str(raw_variable)
                        except ValueError:
                            logging.error('Could not convert string '
                                          f'{raw_variable} to type '
                                          f'{type_to_convert}')
                            sys.exit(1)

                    # integer
                    elif (type_to_convert == int):
                        try:
                            converted_variable = int(raw_variable)
                        except ValueError:
                            logging.error('Could not convert string '
                                          f'{raw_variable} to type '
                                          f'{type_to_convert}')
                            sys.exit(1)

                    # float
                    elif (type_to_convert == float):
                        try:
                            converted_variable = float(raw_variable)
                        except ValueError:
                            logging.error('Could not convert string '
                                          f'{raw_variable} to type '
                                          f'{type_to_convert}')
                            sys.exit(1)

                    #unsupported
                    else:
                        logging.error(f'Unsupported type {type_to_convert} '
                                      f'to read from Fortran namelist')
                        sys.exit(1)

                    return converted_variable

    # variable not found in namelist
    logging.error(f'Could not find {variable} in {namelist}')
    sys.exit(1)


def write_fortran_namelist(name, namelist, nl_class):
    '''
    write Fortran namelist with name for extpar_consistency_check

    the groups of the namelist and its variables are defined in nl_class,
    the values of the variables are stored in namelist
    '''

    with open(name, 'w') as f:
        for group in nl_class.variables.keys():
            f.write(f'{group}\n')
            for var in nl_class.variables[group]:

                if isinstance(namelist[var], str):
                    if (str(namelist[var]) and str(namelist[var]).strip()):
                        f.write(f"  {var} = '{namelist[var]}',\n")
                    else:
                        f.write(f"  {var} = '""'\n")
                else:
                    f.write(f'  {var} = {namelist[var]},\n')

            f.write('/\n')


class InputTclim:
    '''
    define structure of former namelist "INPUT_TCLIM"
    '''

    def __init__(self):

        self.variables = {'&t_clim_raw_data':{'it_cl_type'}}

        self.variables.update({'&t_clim_io_extpar':{'t_clim_buffer_file'}})


class InputAlb:
    '''
    define structure of former namelist "INPUT_ALB"
    '''

    def __init__(self):

        self.variables = {'&alb_raw_data':{'raw_data_alb_path',
                                           'raw_data_alb_filename',
                                           'ialb_type'}}

        self.variables.update({'&alnid_raw_data':{'raw_data_alb_path',
                                                  'raw_data_alnid_filename',
                                                  }})

        self.variables.update({'&aluvd_raw_data':{'raw_data_alb_path',
                                                  'raw_data_aluvd_filename',
                                                  }})

        self.variables.update({'&alb_io_extpar':{'alb_buffer_file'}})


class InputNdvi:
    '''
    define structure of former namelist "INPUT_NDVI"
    '''

    def __init__(self):

        self.variables = {'&ndvi_raw_data':{'raw_data_ndvi_path',
                                            'raw_data_ndvi_filename'}}

        self.variables.update({'&ndvi_io_extpar':{'ndvi_buffer_file'}})


class InputEmiss:
    '''
    define structure of former namelist "INPUT_EMISS"
    '''

    def __init__(self):

        self.variables = {'&emiss_raw_data':{'raw_data_emiss_path',
                                             'raw_data_emiss_filename'}}

        self.variables.update({'&emiss_io_extpar': {'emiss_buffer_file'}})


class InputEra:
    '''
    define structure of namelist "INPUT_ERA"
    '''

    def __init__(self):

        self.variables = {'&era_raw_data':{'iera_type'}}

        # required due to strange bug on Piz Daint
        dict = {'&era_io_extpar':{'era_buffer_file'}}

        self.variables.update(dict)


class InputAhf:
    '''
    define structure of namelist "INPUT_AHF"
    '''

    def __init__(self):

        self.variables = {'&ahf_raw_data': {'iahf_type'}}

        self.variables.update({'&ahf_io_extpar': {'ahf_buffer_file'}})


class InputIsa:
    '''
    define structure of namelist "INPUT_ISA"
    '''

    def __init__(self):

        self.variables = {'&isa_raw_data': {'isa_type'}}

        self.variables.update({'&isa_io_extpar': {'isa_buffer_file'}})
