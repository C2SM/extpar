#!/usr/bin/env python3
import logging
import os
import sys
import subprocess
import netCDF4 as nc
import numpy as np

# extpar modules from lib
try:
    from extpar.lib import (
        utilities as utils,
        grid_def,
        buffer,
        metadata,
        fortran_namelist,
        environment as env,
    )
except ImportError:
    import utilities as utils
    import grid_def
    import buffer
    import metadata
    import fortran_namelist
    import environment as env
from namelist import input_edgar as iedgar

# initialize logger
logging.basicConfig(filename='extpar_edgar_to_buffer.log',
                    level=logging.INFO,
                    format='%(message)s',
                    filemode='w')

logging.info('============= start extpar_edgar_to_buffer ======')
logging.info('')

# print a summary of the environment
env.check_environment_for_extpar(__file__)

# check HDF5
lock = env.check_hdf5_threadsafe()

# get number of OpenMP threads for CDO
omp = env.get_omp_num_threads()

# unique names for files written to system to allow parallel execution
grid = 'grid_description_edgar'  # name for grid description file
reduced_grid = 'reduced_icon_grid_edgar.nc'  # name for reduced icon grid
weights = 'weights_edgar'  # name for weights of spatial interpolation

# names for output of CDO
edgar_bc_cdo = 'edgar_bc_ycon.nc'
edgar_oc_cdo = 'edgar_oc_ycon.nc'
edgar_so2_cdo = 'edgar_so2_ycon.nc'

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= delete files from old runs =======')
logging.info('')

utils.remove(grid)
utils.remove(reduced_grid)
utils.remove(weights)
utils.remove(edgar_bc_cdo)
utils.remove(edgar_oc_cdo)
utils.remove(edgar_so2_cdo)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= init variables from namelist =====')
logging.info('')

igrid_type, grid_namelist = utils.check_gridtype('INPUT_grid_org')

if (igrid_type == 1):
    path_to_grid = \
        fortran_namelist.read_variable(grid_namelist,
                                       'icon_grid_dir',
                                       str)

    icon_grid = \
        fortran_namelist.read_variable(grid_namelist,
                                       'icon_grid_nc_file',
                                       str)

    icon_grid = utils.clean_path(path_to_grid, icon_grid)

    tg = grid_def.IconGrid(icon_grid)

    grid = tg.reduce_grid(reduced_grid)

elif (igrid_type == 2):
    raise exception("EDGAR emission data only works with ICON")

raw_data_edgar_bc = utils.clean_path(iedgar['raw_data_edgar_path'],
                                     iedgar['raw_data_edgar_filename_bc'])
raw_data_edgar_oc = utils.clean_path(iedgar['raw_data_edgar_path'],
                                     iedgar['raw_data_edgar_filename_oc'])
raw_data_edgar_so2 = utils.clean_path(iedgar['raw_data_edgar_path'],
                                      iedgar['raw_data_edgar_filename_so2'])

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= initialize metadata ==============')
logging.info('')

lat_meta = metadata.Lat()
lon_meta = metadata.Lon()

edgarbc_meta = metadata.EdgarBC()
edgaroc_meta = metadata.EdgarOC()
edgarso2_meta = metadata.EdgarSO2()

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= write FORTRAN namelist ===========')
logging.info('')

input_edgar = fortran_namelist.InputEdgar()
fortran_namelist.write_fortran_namelist('INPUT_edgar', iedgar, input_edgar)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= CDO: remap to target grid ========')
logging.info('')

# calculate weights
utils.launch_shell('cdo', lock, '-f', 'nc4', '-P', omp, f'genycon,{grid}',
                   tg.cdo_sellonlat(), raw_data_edgar_bc, weights)

# regrid
utils.launch_shell('cdo', lock, '-f',
                   'nc4', '-P', omp, f'-remap,{grid},{weights}',
                   tg.cdo_sellonlat(), raw_data_edgar_bc, edgar_bc_cdo)
utils.launch_shell('cdo', lock, '-f',
                   'nc4', '-P', omp, f'-remap,{grid},{weights}',
                   tg.cdo_sellonlat(), raw_data_edgar_oc, edgar_oc_cdo)
utils.launch_shell('cdo', lock, '-f',
                   'nc4', '-P', omp, f'-remap,{grid},{weights}',
                   tg.cdo_sellonlat(), raw_data_edgar_so2, edgar_so2_cdo)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= reshape CDO output ===============')
logging.info('')

edgar_bc_nc = nc.Dataset(edgar_bc_cdo, "r")
edgar_oc_nc = nc.Dataset(edgar_oc_cdo, "r")
edgar_so2_nc = nc.Dataset(edgar_so2_cdo, "r")

# infer coordinates/dimensions form CDO file
ie_tot = len(edgar_bc_nc.dimensions['cell'])
je_tot = 1
ke_tot = 1
lon = np.rad2deg(
    np.reshape(edgar_bc_nc.variables['clon'][:], (ke_tot, je_tot, ie_tot)))
lat = np.rad2deg(
    np.reshape(edgar_bc_nc.variables['clat'][:], (ke_tot, je_tot, ie_tot)))

edgar_bc = np.reshape(edgar_bc_nc.variables['emi_bc'][:],
                      (ke_tot, je_tot, ie_tot))
edgar_oc = np.reshape(edgar_oc_nc.variables['emi_oc'][:],
                      (ke_tot, je_tot, ie_tot))
edgar_so2 = np.reshape(edgar_so2_nc.variables['emi_so2'][:],
                       (ke_tot, je_tot, ie_tot))

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= write to buffer file =============')
logging.info('')

# init buffer file
buffer_file = buffer.init_netcdf(iedgar['edgar_buffer_file'], je_tot, ie_tot)

# write lat/lon
buffer.write_field_to_buffer(buffer_file, lon, lon_meta)
buffer.write_field_to_buffer(buffer_file, lat, lat_meta)

# write edgar fields
buffer.write_field_to_buffer(buffer_file, edgar_bc, edgarbc_meta)
buffer.write_field_to_buffer(buffer_file, edgar_oc, edgaroc_meta)
buffer.write_field_to_buffer(buffer_file, edgar_so2, edgarso2_meta)

buffer.close_netcdf(buffer_file)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= clean up =========================')
logging.info('')

utils.remove(weights)
utils.remove(edgar_bc_cdo)
utils.remove(edgar_oc_cdo)
utils.remove(edgar_so2_cdo)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= extpar_edgar_to_buffer done =======')
logging.info('')
