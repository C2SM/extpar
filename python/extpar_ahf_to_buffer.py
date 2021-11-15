#!/usr/bin/env python3
import logging
import os
import sys
import subprocess
import netCDF4 as nc
import numpy as np

# extpar modules from lib
import utilities as utils
import grid_def
import buffer
import metadata
import fortran_namelist
import environment as env
from namelist import input_ahf as iahf

# initialize logger
logging.basicConfig(filename='extpar_ahf_to_buffer.log',
                    level=logging.INFO,
                    format='%(message)s',
                    filemode='w')


logging.info('============= start extpar_ahf_to_buffer =======')
logging.info('')

# print a summary of the environment
env.check_environment_for_extpar(__file__)

# get number of OpenMP threads for CDO
omp = env.get_omp_num_threads()

# unique names for files written to system to allow parallel execution
grid = 'grid_description_ahf'  # name for grid description file
reduced_grid = 'reduced_icon_grid_ahf.nc'  # name for reduced icon grid
weights = 'weights_ahf.nc'        # name for weights of spatial interpolation

# names for output of CDO
ahf_cdo = 'ahf_ycon.nc'
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= delete files from old runs =======')
logging.info('')

utils.remove(grid)
utils.remove(reduced_grid)
utils.remove(weights)
utils.remove(ahf_cdo)

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

    icon_grid = utils.clean_path(path_to_grid,icon_grid)

    tg = grid_def.IconGrid(icon_grid)

    grid = tg.reduce_grid(reduced_grid)

elif(igrid_type == 2):
    tg = grid_def.CosmoGrid(grid_namelist)
    tg.create_grid_description(grid)

iahf_type = utils.check_ahftype(iahf['iahf_type'])

raw_data_ahf  = utils.clean_path(iahf['raw_data_ahf_path'],
                                 iahf['raw_data_ahf_filename'])

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= initialize metadata ==============')
logging.info('')

lat_meta   = metadata.Lat()
lon_meta   = metadata.Lon()

if (iahf_type == 1):
    ahf_meta  = metadata.Ahf_2min()

elif (iahf_type == 2):
    ahf_meta  = metadata.Ahf_30sec()

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= write FORTRAN namelist ===========')
logging.info( '')

input_ahf = fortran_namelist.InputAhf()
fortran_namelist.write_fortran_namelist('INPUT_AHF', iahf, input_ahf)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= CDO: remap to target grid ========')
logging.info('')

lock = env.check_hdf5_threadsafe()

# calculate weights
utils.launch_shell('cdo', '-f', 'nc4', lock, '-P', omp, f'genbil,{grid}',
                   tg.cdo_sellonlat(),
                   raw_data_ahf, weights)

# regrid AHF
utils.launch_shell('cdo', '-f', 'nc4', lock, '-P', omp, 
                   f'settaxis,1111-01-01,0,1mo',
                   f'-remap,{grid},{weights}',
                   tg.cdo_sellonlat(),
                   raw_data_ahf, ahf_cdo)


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= reshape CDO output ===============')
logging.info('')

ahf_nc = nc.Dataset(ahf_cdo, "r")

if (igrid_type == 1):

    # infer coordinates/dimensions form CDO file
    ie_tot = len(ahf_nc.dimensions['cell'])
    je_tot = 1
    ke_tot = 1
    lon    = np.rad2deg(np.reshape(ahf_nc.variables['clon'][:], 
                                   (ke_tot, je_tot, ie_tot)))
    lat    = np.rad2deg(np.reshape(ahf_nc.variables['clat'][:],
                                   (ke_tot, je_tot, ie_tot)))

else:

    # infer coordinates/dimensions from tg
    lat, lon = tg.latlon_cosmo_to_latlon_regular()
    ie_tot   = tg.ie_tot
    je_tot   = tg.je_tot
    ke_tot   = tg.ke_tot

ahf  = np.reshape(ahf_nc.variables['AHF'][:], 
                  (1, ke_tot, je_tot, ie_tot))

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= write to buffer file =============')
logging.info( '')

# init buffer file
buffer_file = buffer.init_netcdf(iahf['ahf_buffer_file'], je_tot, ie_tot)


# write lat/lon
buffer.write_field_to_buffer(buffer_file, lon, lon_meta)
buffer.write_field_to_buffer(buffer_file, lat, lat_meta)

# write AHF fields
buffer.write_field_to_buffer(buffer_file, ahf, ahf_meta)

buffer.close_netcdf(buffer_file)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= clean up =========================')
logging.info('')

utils.remove(weights)
utils.remove(ahf_cdo)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= extpar_ahf_to_buffer done ========')
logging.info( '')
