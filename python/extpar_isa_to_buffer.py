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
from namelist import input_isa as iisa

# initialize logger
logging.basicConfig(filename='extpar_isa_to_buffer.log',
                    level=logging.INFO,
                    format='%(message)s',
                    filemode='w')


logging.info('============= start extpar_isa_to_buffer =======')
logging.info('')

# print a summary of the environment
env.check_environment_for_extpar(__file__)

# get number of OpenMP threads for CDO
omp = env.get_omp_num_threads()

# unique names for files written to system to allow parallel execution
grid = 'grid_description_isa'  # name for grid description file
reduced_grid = 'reduced_icon_grid_isa.nc'  # name for reduced icon grid
weights = 'weights_isa.nc'        # name for weights of spatial interpolation

# names for output of CDO
isa_cdo = 'isa_ycon.nc'
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= delete files from old runs =======')
logging.info('')

utils.remove(grid)
utils.remove(reduced_grid)
utils.remove(weights)
utils.remove(isa_cdo)

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

    grid = utils.reduce_icon_grid(icon_grid, reduced_grid)

elif(igrid_type == 2):
    tg = grid_def.CosmoGrid(grid_namelist)
    tg.create_grid_description(grid)


isa_type = utils.check_isatype(iisa['isa_type'])

raw_data_isa  = utils.clean_path(iisa['raw_data_isa_path'],
                                 iisa['raw_data_isa_filename'])

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= initialize metadata ==============')
logging.info('')

lat_meta   = metadata.Lat()
lon_meta   = metadata.Lon()

if (isa_type == 1):
    isa_meta  = metadata.Isa_30sec()
elif (isa_type == 1):
    isa_meta  = metadata.Isa_10sec()

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= write FORTRAN namelist ===========')
logging.info( '')

input_isa = fortran_namelist.InputIsa()
fortran_namelist.write_fortran_namelist('INPUT_ISA', iisa, input_isa)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= CDO: remap to target grid ========')
logging.info('')

# calculate weights
utils.launch_shell('cdo', '-f', 'nc4', '-L', '-P', omp,  f'genycon,{grid}',
                   raw_data_isa, weights)

# regrid ISA
utils.launch_shell('cdo', '-f', 'nc4', '-L', '-P', omp, 
                   f'settaxis,1111-01-01,0,1mo',
                   f'-remap,{grid},{weights}', raw_data_isa, isa_cdo)


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= reshape CDO output ===============')
logging.info('')

isa_nc = nc.Dataset(isa_cdo, "r")

if (igrid_type == 1):

    # infer coordinates/dimensions form CDO file
    ie_tot = len(isa_nc.dimensions['cell'])
    je_tot = 1
    ke_tot = 1
    lon    = np.rad2deg(np.reshape(isa_nc.variables['clon'][:], 
                                   (ke_tot, je_tot, ie_tot)))
    lat    = np.rad2deg(np.reshape(isa_nc.variables['clat'][:],
                                   (ke_tot, je_tot, ie_tot)))

else:

    # infer coordinates/dimensions from tg
    lat, lon = tg.latlon_cosmo_to_latlon_regular()
    ie_tot   = tg.ie_tot
    je_tot   = tg.je_tot
    ke_tot   = tg.ke_tot

isa  = np.reshape(isa_nc.variables['ISA'][:], 
                  (1, ke_tot, je_tot, ie_tot))

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= write to buffer file =============')
logging.info( '')

# init buffer file
buffer_file = buffer.init_netcdf(iisa['isa_buffer_file'], je_tot, ie_tot)


# write lat/lon
buffer.write_field_to_buffer(buffer_file, lon, lon_meta)
buffer.write_field_to_buffer(buffer_file, lat, lat_meta)

# write ISA fields
buffer.write_field_to_buffer(buffer_file, isa, isa_meta)

buffer.close_netcdf(buffer_file)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= clean up =========================')
logging.info('')

utils.remove(weights)
utils.remove(isa_cdo)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= extpar_isa_to_buffer done ========')
logging.info( '')
