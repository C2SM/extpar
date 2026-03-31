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
from namelist import input_aot as iaot

# initialize logger
logging.basicConfig(filename='extpar_aot_to_buffer.log',
                    level=logging.INFO,
                    format='%(message)s',
                    filemode='w')

logging.info('============= start extpar_aot_to_buffer =======')
logging.info('')

# print a summary of the environment
env.check_environment_for_extpar(__file__)

# check HDF5
lock = env.check_hdf5_threadsafe()

# get number of OpenMP threads for CDO
omp = env.get_omp_num_threads()

# unique names for files written to system to allow parallel execution
grid = 'grid_description_aot'  # name for grid description file
reduced_grid = 'reduced_icon_grid_aot.nc'  # name for reduced icon grid
weights = 'weights_aot.nc'  # name for weights of spatial interpolation

# names for output of CDO
aot_cdo = 'aot_ycon.nc'
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= delete files from old runs =======')
logging.info('')

utils.remove(grid)
utils.remove(reduced_grid)
utils.remove(weights)
utils.remove(aot_cdo)

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
    tg = grid_def.CosmoGrid(grid_namelist)
    tg.create_grid_description(grid)

aot_type = utils.check_aottype(iaot['iaot_type'])

raw_data_aot = utils.clean_path(iaot['raw_data_aot_path'],
                                iaot['raw_data_aot_filename'])

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= initialize metadata ==============')
logging.info('')

lat_meta = metadata.Lat()
lon_meta = metadata.Lon()

if (aot_type == 1):
    aot_meta = metadata.AotTegen()
elif (aot_type == 2):
    aot_meta = metadata.AotAeroCom()

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= write FORTRAN namelist ===========')
logging.info('')

input_aot = fortran_namelist.InputAot()
fortran_namelist.write_fortran_namelist('INPUT_AOT', iaot, input_aot)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= CDO: remap to target grid ========')
logging.info('')

# Compute longitude and latitude resolution of the raw dataset
raw_data_aot_nc = nc.Dataset(raw_data_aot, "r")
lon = raw_data_aot_nc.variables['lon']
lat = raw_data_aot_nc.variables['lat']
dlon = lon[1] - lon[0]
dlat = lat[1] - lat[0]

# calculate weights
utils.launch_shell('cdo', '-f', 'nc4', lock, '-P', omp, f'genbil,{grid}',
                   tg.cdo_sellonlat(dlon, dlat), raw_data_aot, weights)

# regrid aot
utils.launch_shell('cdo', '-f', 'nc4', lock, '-P', omp,
                   f'settaxis,1111-01-01,0,1mo', f'-remap,{grid},{weights}',
                   tg.cdo_sellonlat(dlon, dlat), raw_data_aot, aot_cdo)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= reshape CDO output ===============')
logging.info('')

aot_nc = nc.Dataset(aot_cdo, "r")

if (igrid_type == 1):

    # infer coordinates/dimensions from CDO file
    ie_tot = len(aot_nc.dimensions['cell'])
    je_tot = 1
    ke_tot = 1
    lon = np.rad2deg(
        np.reshape(aot_nc.variables['clon'][:], (ke_tot, je_tot, ie_tot)))
    lat = np.rad2deg(
        np.reshape(aot_nc.variables['clat'][:], (ke_tot, je_tot, ie_tot)))

else:

    # infer coordinates/dimensions from tg
    lat, lon = tg.latlon_cosmo_to_latlon_regular()
    ie_tot = tg.ie_tot
    je_tot = tg.je_tot
    ke_tot = tg.ke_tot

aot = np.empty((12, 5, ke_tot, je_tot, ie_tot), dtype=aot_meta.type)
aerosol_names = ['black_carbon', 'dust', 'organic', 'sulfate', 'sea_salt']

for i in range(5):
    aerosol_name = aerosol_names[i]
    aot[:, i, :, :, :] = np.reshape(aot_nc.variables[aerosol_name][:, :],
                                    (12, ke_tot, je_tot, ie_tot))

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= write to buffer file =============')
logging.info('')

# init buffer file
buffer_file = buffer.init_netcdf(iaot['aot_buffer_file'], je_tot, ie_tot)

# add 12 months as additional dimension
buffer_file = buffer.add_dimension_month(buffer_file)

# add 5 aerosol types as additional dimension
buffer_file = buffer.add_dimension_aerosols(buffer_file)

# write lat/lon
buffer.write_field_to_buffer(buffer_file, lon, lon_meta)
buffer.write_field_to_buffer(buffer_file, lat, lat_meta)

# write aot fields
buffer.write_field_to_buffer(buffer_file, aot, aot_meta)

buffer.close_netcdf(buffer_file)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= clean up =========================')
logging.info('')

utils.remove(weights)
utils.remove(aot_cdo)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= extpar_aot_to_buffer done ========')
logging.info('')
