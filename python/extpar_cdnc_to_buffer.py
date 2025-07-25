#!/usr/bin/env python3
import logging
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
from namelist import input_cdnc as icdnc

# initialize logger
logging.basicConfig(filename='extpar_cdnc_to_buffer.log',
                    level=logging.INFO,
                    format='%(message)s',
                    filemode='w')

logging.info('============= start extpar_cdnc_to_buffer ======')
logging.info('')

# print a summary of the environment
env.check_environment_for_extpar(__file__)

# check HDF5
lock = env.check_hdf5_threadsafe()

# get number of OpenMP threads for CDO
omp = env.get_omp_num_threads()

# unique names for files written to system to allow parallel execution
grid = 'grid_description_cdnc'  # name for grid description file
reduced_grid = 'reduced_icon_grid_cdnc.nc'  # name for reduced icon grid
weights = 'weights_cdnc'  # name for weights of spatial interpolation

# names for output of CDO
cdnc_cdo = 'cdnc_ycon.nc'
cdnc_tmp = 'cdnc_tmp.nc'

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= delete files from old runs =======')
logging.info('')

utils.remove(grid)
utils.remove(reduced_grid)
utils.remove(weights)
utils.remove(cdnc_cdo)
utils.remove(cdnc_tmp)

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
    raise exception("cdnc data only works with ICON")

cdnc_type = utils.check_cdnctype(icdnc['icdnc_type'])

raw_data_cdnc = utils.clean_path(icdnc['raw_data_cdnc_path'],
                                 icdnc['raw_data_cdnc_filename'])

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= initialize metadata ==============')
logging.info('')

lat_meta = metadata.Lat()
lon_meta = metadata.Lon()

if (cdnc_type == 1):
    cdnc_meta = metadata.CdncQ06()
elif (cdnc_type == 2):
    cdnc_meta = metadata.CdncBR17Adjusted()

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= write FORTRAN namelist ===========')
logging.info('')

input_cdnc = fortran_namelist.InputCdnc()
fortran_namelist.write_fortran_namelist('INPUT_CDNC', icdnc, input_cdnc)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= CDO: remap to target grid ========')
logging.info('')

# calculate weights
utils.launch_shell('cdo', lock, '-f', 'nc4', '-P', omp, f'genycon,{grid}',
                   tg.cdo_sellonlat(), raw_data_cdnc, weights)

# regrid
utils.launch_shell('cdo', lock, '-f',
                   'nc4', '-P', omp, f'-remap,{grid},{weights}',
                   tg.cdo_sellonlat(), raw_data_cdnc, cdnc_cdo)

# set missing values to a minimal value
utils.launch_shell('cdo', '-f', 'nc4', '-P', omp, lock, f'setmisstoc,-1',
                   cdnc_cdo, cdnc_tmp)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= reshape CDO output ===============')
logging.info('')

cdnc_nc = nc.Dataset(cdnc_tmp, "r")

# infer coordinates/dimensions form CDO file
ie_tot = len(cdnc_nc.dimensions['cell'])
je_tot = 1
ke_tot = 1
lon = np.rad2deg(
    np.reshape(cdnc_nc.variables['clon'][:], (ke_tot, je_tot, ie_tot)))
lat = np.rad2deg(
    np.reshape(cdnc_nc.variables['clat'][:], (ke_tot, je_tot, ie_tot)))

cdnc = np.reshape(cdnc_nc.variables['cdnc'][:, :],
                  (12, ke_tot, je_tot, ie_tot))

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= write to buffer file =============')
logging.info('')

# init buffer file
buffer_file = buffer.init_netcdf(icdnc['cdnc_buffer_file'], je_tot, ie_tot)

# add 12 months as additional dimension
buffer_file = buffer.add_dimension_month(buffer_file)

# write lat/lon
buffer.write_field_to_buffer(buffer_file, lon, lon_meta)
buffer.write_field_to_buffer(buffer_file, lat, lat_meta)

# write cdnc field
buffer.write_field_to_buffer(buffer_file, cdnc, cdnc_meta)

buffer.close_netcdf(buffer_file)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= clean up =========================')
logging.info('')

utils.remove(grid)
utils.remove(reduced_grid)
utils.remove(weights)
utils.remove(cdnc_cdo)
utils.remove(cdnc_tmp)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= extpar_cdnc_to_buffer done =======')
logging.info('')
