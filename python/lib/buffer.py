import logging
import sys
import netCDF4 as nc
import numpy as np
'''
Module providing functions to write the buffer file,
it contains:

    -init_netcdf: create buffer file and init dimensions for albedo

    -open_netcdf: open netcdf, if exception occurs exit extpar

    -close_netcdf: close netcdf, if exception occurs exit extpar

    -write_field_to_buffer: interface to call the correct writing-function

    -write_3d_field: write 3d field to buffer file

    -write_4d_field: write 4d field to buffer file
'''


def init_netcdf(buffer_name, len_je, len_ie):
    '''
    create netcdf and init base dimensions ie, je ,ke, mlev
    '''
    # create extpar BUFFER
    buffer = open_netcdf(buffer_name)

    buffer.createDimension('ie', len_ie)
    buffer.createDimension('je', len_je)
    buffer.createDimension('ke', 1)
    buffer.createDimension('mlev', 1)

    extpar_mlev = buffer.createVariable('mlev', np.float32, ('mlev', ))
    extpar_mlev[:] = np.array([1], dtype=np.float32)

    return buffer


def add_dimension_month(buffer):
    '''
    add 12 months as dimension to netCDF
    '''

    buffer.createDimension('time', None)

    extpar_time = buffer.createVariable('time', np.float32, ('time', ))
    extpar_time[:] = np.array([
        11110101, 11110201, 11110301, 11110401, 11110501, 11110601, 11110701,
        11110801, 11110901, 11111001, 11111101, 11111201
    ],
                              dtype=np.float32)

    return buffer


def open_netcdf(buffer_name):
    '''
    open netcdf with buffer_name
    '''
    try:
        buffer = nc.Dataset(buffer_name, "w", format='NETCDF4')
    except:
        logging.error(f'Could not open {buffer_name}', exc_info=True)
        raise

    logging.info(f'create netCDF {buffer_name}')

    return buffer


def close_netcdf(buffer):
    '''
    close netcdf-file buffer
    '''
    try:
        buffer.close()
    except:
        logging.error('Could not close netCDF', exc_info=True)
        raise


def write_field_to_buffer(buffer, field, field_meta):
    '''
    determine number of dimensions of field and call resp. writing-function
    '''

    # determine number of dimensions
    dim_nr = len(field_meta.dim)

    # call the correct writing-function

    # 3d-field
    if (dim_nr == 3):
        write_3d_field(buffer, field, field_meta)

    # 4d-field
    elif (dim_nr == 4):
        write_4d_field(buffer, field, field_meta)

    # unsupported
    else:
        logging.error(f'field {field_meta.name} has {dim_nr} dimensions,'
                      ' currently supported are 3 or 4 dimensions')
        raise


def write_3d_field(buffer, field_3d, meta):
    '''
    write 3d field.data to buffer

    buffer is the netCDF file to write field_3d
    field_3d needs to have the same shape as netcdf_var
    the metadata for each variable is stored in meta
    meta is defined in module metadata
    '''
    netcdf_var = buffer.createVariable(meta.name, meta.type,
                                       (meta.dim[0], meta.dim[1], meta.dim[2]))

    netcdf_var.standard_name = meta.long
    netcdf_var.long_name = meta.long
    netcdf_var.units = meta.units

    try:
        netcdf_var[:, :, :] = field_3d.data
    except ValueError:
        logging.error('Error during netCDF IO', exc_info=True)
        raise

    logging.info(f'3D-field {meta.name} written')


def write_4d_field(buffer, field_4d, meta):
    '''
    write 4d field.data to buffer

    buffer is the netCDF file to write field_4d
    field_4d needs to have the same shape as netcdf_var
    the metadata for each variable is stored in meta
    meta is defined in module metadata
    '''
    netcdf_var = buffer.createVariable(
        meta.name, meta.type,
        (meta.dim[0], meta.dim[1], meta.dim[2], meta.dim[3]))

    netcdf_var.standard_name = meta.long
    netcdf_var.long_name = meta.long
    netcdf_var.units = meta.units

    try:
        netcdf_var[:, :, :, :] = field_4d.data
    except ValueError:
        logging.error('Error during netCDF IO', exc_info=True)
        raise

    logging.info(f'4D-field {meta.name} written')
