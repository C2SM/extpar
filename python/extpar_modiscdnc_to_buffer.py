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
from namelist import input_modis_cdnc as icdnc

# initialize logger
logging.basicConfig(filename='extpar_modis_cdnc_to_buffer.log',
                    level=logging.INFO,
                    format='%(message)s',
                    filemode='w')

logging.info('============= start extpar_modis_cdnc_to_buffer ======')
logging.info('')

# print a summary of the environment
env.check_environment_for_extpar(__file__)

# check HDF5
lock = env.check_hdf5_threadsafe()

# get number of OpenMP threads for CDO
omp = env.get_omp_num_threads()

# unique names for files written to system to allow parallel execution
grid = 'grid_description_modis_cdnc'  # name for grid description file
reduced_grid = 'reduced_icon_grid_modis_cdnc.nc'  # name for reduced icon grid
weights = 'weights_modis_cdnc'  # name for weights of spatial interpolation

# names for output of CDO
modis_cdnc_cdo = 'modis_cdnc_ycon.nc'

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= delete files from old runs =======')
logging.info('')

utils.remove(grid)
utils.remove(reduced_grid)
utils.remove(weights)
utils.remove(modis_cdnc_cdo)

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
    raise exception("MODIS cdnc data only works with ICON")

raw_data_modis_cdnc = utils.clean_path(icdnc['raw_data_modis_cdnc_path'],
                                       icdnc['raw_data_modis_cdnc_filename'])

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= initialize metadata ==============')
logging.info('')

lat_meta = metadata.Lat()
lon_meta = metadata.Lon()

modiscdnc_meta = metadata.ModisCdnc()

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= write FORTRAN namelist ===========')
logging.info('')

input_modis_cdnc = fortran_namelist.InputModisCdnc()
fortran_namelist.write_fortran_namelist('INPUT_modis_cdnc', icdnc,
                                        input_modis_cdnc)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= CDO: remap to target grid ========')
logging.info('')

# calculate weights
utils.launch_shell('cdo', lock, '-f', 'nc4', '-P', omp, f'genycon,{grid}',
                   tg.cdo_sellonlat(), raw_data_modis_cdnc, weights)

# regrid
utils.launch_shell('cdo', lock, '-f',
                   'nc4', '-P', omp, f'-remap,{grid},{weights}',
                   tg.cdo_sellonlat(), raw_data_modis_cdnc, modis_cdnc_cdo)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= reshape CDO output ===============')
logging.info('')

modis_cdnc_nc = nc.Dataset(modis_cdnc_cdo, "r")

# infer coordinates/dimensions form CDO file
ie_tot = len(modis_cdnc_nc.dimensions['cell'])
je_tot = 1
ke_tot = 1
lon = np.rad2deg(
    np.reshape(modis_cdnc_nc.variables['clon'][:], (ke_tot, je_tot, ie_tot)))
lat = np.rad2deg(
    np.reshape(modis_cdnc_nc.variables['clat'][:], (ke_tot, je_tot, ie_tot)))

modis_cdnc = np.reshape(modis_cdnc_nc.variables['modis_cdnc'][:],
                        (ke_tot, je_tot, ie_tot))

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= write to buffer file =============')
logging.info('')

# init buffer file
buffer_file = buffer.init_netcdf(icdnc['modis_cdnc_buffer_file'], je_tot,
                                 ie_tot)

# write lat/lon
buffer.write_field_to_buffer(buffer_file, lon, lon_meta)
buffer.write_field_to_buffer(buffer_file, lat, lat_meta)

# write modis cdnc field
buffer.write_field_to_buffer(buffer_file, modis_cdnc, modis_cdnc_meta)

buffer.close_netcdf(buffer_file)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= clean up =========================')
logging.info('')

utils.remove(weights)
utils.remove(modis_cdnc_cdo)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= extpar_modis_cdnc_to_buffer done =======')
logging.info('')
