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
from namelist import input_era as iera

# initialize logger
logging.basicConfig(filename='extpar_era_to_buffer.log',
                    level=logging.INFO,
                    format='%(message)s',
                    filemode='w')


logging.info('============= start extpar_era_to_buffer =======')
logging.info('')

# print a summary of the environment
env.check_environment_for_extpar(__file__)

# get number of OpenMP threads for CDO
omp = env.get_omp_num_threads()

# unique names for files written to system to allow parallel execution
grid = 'grid_description_era'  # name for grid description file
reduced_grid = 'reduced_icon_grid_era.nc'  # name for reduced icon grid
weights = 'weights_era.nc'        # name for weights of spatial interpolation

# names for output of CDO
sst_cdo = 'sst_ycon.nc'
t2m_cdo = 't2m_ycon.nc'
oro_cdo = 'oro_ycon.nc'
sd_cdo = 'sd_ycon.nc'

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= delete files from old runs =======')
logging.info('')

utils.remove(grid)
utils.remove(reduced_grid)
utils.remove(weights)
utils.remove(sst_cdo)
utils.remove(t2m_cdo)
utils.remove(oro_cdo)
utils.remove(sd_cdo)

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


raw_data_sst  = utils.clean_path(iera['raw_data_era_path'],
                                  iera['raw_data_era_SST'])

raw_data_t2m  = utils.clean_path(iera['raw_data_era_path'],
                                  iera['raw_data_era_T2M'])

raw_data_oro  = utils.clean_path(iera['raw_data_era_path'],
                                  iera['raw_data_era_ORO'])

raw_data_sd   = utils.clean_path(iera['raw_data_era_path'],
                                  iera['raw_data_era_SD'])
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= initialize metadata ==============')
logging.info('')

lat_meta   = metadata.Lat()
lon_meta   = metadata.Lon()

sst_meta  = metadata.SstEra()
t2m_meta  = metadata.T2mEra()
oro_meta  = metadata.OroEra()
sd_meta   = metadata.SdEra()

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= write FORTRAN namelist ===========')
logging.info( '')

input_era = fortran_namelist.InputEra()
fortran_namelist.write_fortran_namelist('INPUT_ERA', iera, input_era)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= CDO: remap to target grid ========')
logging.info('')

# calculate weights
utils.launch_shell('cdo', '-f', 'nc4', '-L', '-P', omp, f'genycon,{grid}',
                   raw_data_sst, weights)

# regrid SST
utils.launch_shell('cdo', '-f', 'nc4', '-L', '-P', omp, 
                   f'settaxis,1111-01-01,0,1mo',
                   f'-remap,{grid},{weights}', raw_data_sst, sst_cdo)

# regrid T2M
utils.launch_shell('cdo', '-f', 'nc4', '-L', '-P', omp, 
                   f'settaxis,1111-01-01,0,1mo',
                   f'-remap,{grid},{weights}', raw_data_t2m, t2m_cdo)

# regrid ORO
utils.launch_shell('cdo', '-f', 'nc4', '-L', '-P', omp, 
                   f'settaxis,1111-01-01,0,1mo',
                   f'-remap,{grid},{weights}', raw_data_oro, oro_cdo)

# regrid SD
utils.launch_shell('cdo', '-f', 'nc4', '-L', '-P', omp, 
                   f'settaxis,1111-01-01,0,1mo',
                   f'-remap,{grid},{weights}', raw_data_sd, sd_cdo)


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= reshape CDO output ===============')
logging.info('')

sst_nc = nc.Dataset(sst_cdo, "r")
t2m_nc = nc.Dataset(t2m_cdo, "r")
oro_nc = nc.Dataset(oro_cdo, "r")
sd_nc = nc.Dataset(sd_cdo, "r")

if (igrid_type == 1):

    # infer coordinates/dimensions form CDO file
    ie_tot = len(sst_nc.dimensions['cell'])
    je_tot = 1
    ke_tot = 1
    lon    = np.rad2deg(np.reshape(sst_nc.variables['clon'][:], 
                                   (ke_tot, je_tot, ie_tot)))
    lat    = np.rad2deg(np.reshape(sst_nc.variables['clat'][:],
                                   (ke_tot, je_tot, ie_tot)))

else:

    # infer coordinates/dimensions from tg
    lat, lon = tg.latlon_cosmo_to_latlon_regular()
    ie_tot   = tg.ie_tot
    je_tot   = tg.je_tot
    ke_tot   = tg.ke_tot

sst  = np.reshape(sst_nc.variables['sst'][:,:], 
                   (12, ke_tot, je_tot, ie_tot))

t2m  = np.reshape(t2m_nc.variables['2t'][:,:], 
                   (12, ke_tot, je_tot, ie_tot))

oro  = np.reshape(oro_nc.variables['HSURF'][:,:], 
                   (ke_tot, je_tot, ie_tot))

sd  = np.reshape(sd_nc.variables['sd'][:,:], 
                   (12, ke_tot, je_tot, ie_tot))

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= compute W_SNOW ===================')
logging.info( '')

sd = sd * 1000.

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= write to buffer file =============')
logging.info( '')

# init buffer file
buffer_file = buffer.init_netcdf(iera['era_buffer_file'], je_tot, ie_tot)

# add 12 months as additional dimension
buffer_file = buffer.add_dimension_month(buffer_file)

# write lat/lon
buffer.write_field_to_buffer(buffer_file, lon, lon_meta)
buffer.write_field_to_buffer(buffer_file, lat, lat_meta)

# write ERA fields
buffer.write_field_to_buffer(buffer_file, sst, sst_meta)
buffer.write_field_to_buffer(buffer_file, t2m, t2m_meta)
buffer.write_field_to_buffer(buffer_file, oro, oro_meta)
buffer.write_field_to_buffer(buffer_file, sd, sd_meta)

buffer.close_netcdf(buffer_file)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= clean up =========================')
logging.info('')

utils.remove(weights)
utils.remove(sst_cdo)
utils.remove(t2m_cdo)
utils.remove(oro_cdo)
utils.remove(sd_cdo)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= extpar_era_to_buffer done ========')
logging.info( '')
