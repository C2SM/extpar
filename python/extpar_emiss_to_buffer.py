#!/usr/bin/env python3
import logging
import os
import sys
import subprocess
import netCDF4 as nc
import numpy as np

# extpar modules from lib
import grid_def
import buffer
import metadata
import fortran_namelist
import utilities as utils
import environment as env
from namelist import input_emiss as iemiss

# initialize logger
logging.basicConfig(filename='extpar_emiss_to_buffer.log',
                    level=logging.INFO,
                    format='%(message)s',
                    filemode='w')


logging.info('============= start extpar_emiss_to_buffer ======')
logging.info('')

# print a summary of the environment
env.check_environment_for_extpar(__file__)

# get number of OpenMP threads for CDO
omp = env.get_omp_num_threads()

# unique names for files written to system to allow parallel execution
grid = 'grid_description_emiss'  # name for grid description file
reduced_grid = 'reduced_icon_grid_emiss.nc'  # name for reduced icon grid
weights = 'weights_emiss'        # name for weights of spatial interpolation

# names for output of CDO
emiss_cdo_1 = 'emiss_1.nc'
emiss_cdo_2 = 'emiss_2.nc'
emiss_cdo_3 = 'emiss_3.nc'
emiss_cdo_4 = 'emiss_4.nc'
emiss_cdo_5 = 'emiss_5.nc'

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= delete files from old runs =======')
logging.info('')

utils.remove(grid)
utils.remove(reduced_grid)
utils.remove(weights)
utils.remove(emiss_cdo_1)
utils.remove(emiss_cdo_2)
utils.remove(emiss_cdo_3)
utils.remove(emiss_cdo_4)
utils.remove(emiss_cdo_5)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= init variables from namelist =====')
logging.info('')

iemiss_type = utils.check_emisstype(iemiss['iemiss_type'])

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

raw_data_emiss  = utils.clean_path(iemiss['raw_data_emiss_path'],
                                   iemiss['raw_data_emiss_filename'])

# determine varname for postprocessing of CDO output
var = utils.determine_emiss_varnames(iemiss_type)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= initialize metadata ==============')
logging.info('')

lat_meta   = metadata.Lat()
lon_meta   = metadata.Lon()

emiss_meta  = metadata.EmissMean()
max_meta   = metadata.EmissMax()
mrat_meta  = metadata.EmissMrat()

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= write FORTRAN namelist ===========')
logging.info( '')

input_emiss = fortran_namelist.InputEmiss()
fortran_namelist.write_fortran_namelist('INPUT_EMISS', iemiss, input_emiss)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= CDO: remap to target grid ========')
logging.info('')

utils.launch_shell('cp', raw_data_emiss, emiss_cdo_1)

# calculate weights
utils.launch_shell('cdo', '-f', 'nc4', '-P', omp, '--silent',f'genycon,{grid}',
                   emiss_cdo_1, weights)

# Check of artificial low values
# (useful range is between approx. 0.6 and 1. for earth surface)
utils.launch_shell('cdo',  '-f', 'nc4', '-P', omp,
                   f'-expr,{var}=({var}<0.5)'
                   f'?-999:{var};',
                   emiss_cdo_1, emiss_cdo_2)

# Ensure artificial low values are set to missing
utils.launch_shell('cdo', '-f', 'nc4', '-P', omp,
                   'setmissval,-999', emiss_cdo_2, emiss_cdo_3)

# Set missing values to nearest neighbors -> useful values for high-res grids
utils.launch_shell('cdo', '-f', 'nc4', '-P', omp,
                   'setmisstonn', emiss_cdo_3, emiss_cdo_4)

# regrid 1 
# -L option prevents crash for non-threads save compilations of HDF5-library
utils.launch_shell('cdo', '-f', 'nc4', '-P', omp, '-L', 
                   f'settaxis,1111-01-01,0,1mo',
                   f'-remap,{grid},{weights}', emiss_cdo_4, emiss_cdo_5)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= reshape CDO output ===============')
logging.info('')

emiss_nc = nc.Dataset(emiss_cdo_5, "r")

if (igrid_type == 1):

    # infer coordinates/dimensions form CDO file
    ie_tot = len(emiss_nc.dimensions['cell'])
    je_tot = 1
    ke_tot = 1
    lon    = np.rad2deg(np.reshape(emiss_nc.variables['clon'][:], 
                                   (ke_tot, je_tot, ie_tot)))
    lat    = np.rad2deg(np.reshape(emiss_nc.variables['clat'][:],
                                   (ke_tot, je_tot, ie_tot)))

else:

    # infer coordinates/dimensions from tg
    lat, lon = tg.latlon_cosmo_to_latlon_regular()
    ie_tot   = tg.ie_tot
    je_tot   = tg.je_tot
    ke_tot   = tg.ke_tot

emiss  = np.reshape(emiss_nc.variables[var][:,:], 
                    (12, ke_tot, je_tot, ie_tot))

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= compute EMISS_MAX and EMISS_MRAT ==')
logging.info('')

# calculate maxval over 12 month
emiss_max = np.amax(np.reshape(emiss_nc.variables[var][:,:], 
                               (12, ke_tot, je_tot, ie_tot)), axis=0)

# calculate ratio of emiss/emiss_max per month and set 'missing value' to -1
emiss_mrat = np.empty((12, ke_tot, je_tot, ie_tot), dtype=mrat_meta.type)

for t in np.arange(12):
    emiss_mrat[t,:,:,:] = np.divide(emiss[t,:,:,:], emiss_max[:,:,:],
                                    where=emiss_max[:,:,:] != 0.0)
    emiss_mrat[t,:,:,:] = np.where(emiss_max[:,:,:] <= 0.0, -1.0, 
                                   emiss_mrat[t,:,:,:])

# debug -> print these statistics setting level=logging.DEBUG
logging.debug('Diagnostics:')

logging.debug('   EMISS max')
logging.debug('   min: {0:8.5f} mean: {1:8.5f} max: {2:8.5f}'
              .format(np.min(emiss_max),
                      np.mean(emiss_max),
                      np.max(emiss_max)))

logging.debug('   EMISS')
for t in np.arange(12):
    logging.debug('   min: {0:8.5f} mean: {1:8.5f} max: {2:8.5f}'
                  .format(np.min(emiss[t,:,:,:]),
                          np.mean(emiss[t,:,:,:]),
                          np.max(emiss[t,:,:,:])))

logging.debug('   EMISS mrat')
for t in np.arange(12):
    logging.debug('   min: {0:8.5f} mean: {1:8.5f} max: {2:8.5f}'
                  .format(np.min(emiss_mrat[t,:,:,:]),
                          np.mean(emiss_mrat[t,:,:,:]),
                          np.max(emiss_mrat[t,:,:,:])))

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= write to buffer file =============')
logging.info( '')

# init buffer file
buffer_file = buffer.init_netcdf(iemiss['emiss_buffer_file'], je_tot, ie_tot)

# add 12 months as additional dimension
buffer_file = buffer.add_dimension_month(buffer_file)

# write lat/lon
buffer.write_field_to_buffer(buffer_file, lon, lon_meta)
buffer.write_field_to_buffer(buffer_file, lat, lat_meta)

# write EMISS fields
buffer.write_field_to_buffer(buffer_file, emiss_max, max_meta)
buffer.write_field_to_buffer(buffer_file, emiss, emiss_meta)
buffer.write_field_to_buffer(buffer_file, emiss_mrat, mrat_meta)

buffer.close_netcdf(buffer_file)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= clean up =========================')
logging.info('')

utils.remove(weights)
utils.remove(emiss_cdo_1)
utils.remove(emiss_cdo_2)
utils.remove(emiss_cdo_3)
utils.remove(emiss_cdo_4)
utils.remove(emiss_cdo_5)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info( '')
logging.info( '============= extpar_emiss_to_buffer done =======')
logging.info( '')
