#!/usr/bin/env python3
import logging
import netCDF4 as nc
import numpy as np

from joblib import Parallel, delayed, dump, load
from tqdm import tqdm
from sklearn.neighbors import BallTree

# extpar modules from lib
try:
    from extpar.lib import (
        grid_def,
        buffer,
        metadata,
        fortran_namelist,
        utilities as utils,
        environment as env,
    )
    from extpar.lib.namelist import input_art as iart
except ImportError:
    import grid_def
    import buffer
    import metadata
    import fortran_namelist
    import utilities as utils
    import environment as env
    from namelist import input_art as iart


def get_nearest_neighbors(lon_array,
                          lat_array,
                          lon_point,
                          lat_point,
                          metric='haversine'):
    """
    Find the indices of n closest cells in a grid, relative to a given latitude/longitude.    
    The function builds a BallTree from the provided 1D or 2D array of latitude and longitude and queries it to find n nearest neghbors to a given point.
    Borrowed from iconarray utilities    
    """
    lon_array, lat_array, lon_point, lat_point = [
        np.deg2rad(arr)
        for arr in [lon_array, lat_array, lon_point, lat_point]
    ]
    lon_lat_array = np.column_stack((lon_array.flatten(), lat_array.flatten()))
    points = np.column_stack((lon_point, lat_point))
    indices = BallTree(lon_lat_array, metric=metric,
                       leaf_size=3).query(points, k=1)[1].squeeze()
    return indices


def get_neighbor_index(index, lons, lats, hlon, hlat, idxs):
    xx = np.ones((hlon.size))
    idxs[index, :] = get_nearest_neighbors(lons, lats, hlon, hlat[index] * xx)


def get_memory_map(mat, filename_mmap):
    dump(mat, filename_mmap)
    mat = load(filename_mmap, mmap_mode='r+')
    return mat


def generate_memory_map(raw_lus, soiltype_memmap_filename,
                        nearest_neighbor_index_memmap_filename):
    idxs = -1 * np.ones(raw_lus.shape, int)
    idxs = get_memory_map(idxs, nearest_neighbor_index_memmap_filename)
    lus = get_memory_map(raw_lus, soiltype_memmap_filename)
    return lus, idxs


def calculate_soil_fraction(tg, lus, idxs, ncpu=2):
    """
    lus: LU classes from HWSD data
    idxs: indices corrsponding to icon grid for each grid in HWSD data
    tg: ICON grid
    """
    soil_types = np.arange(1, 14)
    fracs = np.zeros((tg.lons.size, soil_types.size))
    grid_ids, grid_counts = np.unique(idxs, return_counts=True)

    def get_fraction_per_soil_type(lu):
        grid_class, grid_count = np.unique(np.where(lus == lu, idxs, -1),
                                           return_counts=True)
        for grid_id in np.arange(tg.lons.size):
            frac = np.array(grid_count[grid_class == grid_id] /
                            grid_counts[grid_ids == grid_id])
            if len(frac) != 0:
                fracs[grid_id, lu - 1] = frac

    Parallel(n_jobs=ncpu,
             max_nbytes='100M',
             mmap_mode='w+',
             backend='threading')(delayed(get_fraction_per_soil_type)(lu)
                                  for lu in tqdm(soil_types))
    return fracs


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# initialize logger
logging.basicConfig(
    filename="extpar_art_to_buffer.log",
    level=logging.INFO,
    format="%(message)s",
    filemode="w",
)

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
logging.info("============= start extpar_art_to_buffer =======")
logging.info("")

# Use all available CPUs
omp = int(env.get_omp_num_threads())

# unique names for files written to system to allow parallel execution
soiltype_memmap_filename = 'memmap_soiltype'
nearest_neighbor_index_memmap_filename = 'memmap_index'

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info(
    '============= delete files from old runs =========================')
logging.info('')

utils.remove(soiltype_memmap_filename)
utils.remove(nearest_neighbor_index_memmap_filename)

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
logging.info("")
logging.info("============= init variables from namelist =====")
logging.info("")

igrid_type, grid_namelist = utils.check_gridtype('INPUT_grid_org')
raw_data_art = utils.clean_path(iart['raw_data_art_path'],
                                iart['raw_data_art_filename'])

if igrid_type == 1:
    path_to_grid = fortran_namelist.read_variable(grid_namelist,
                                                  "icon_grid_dir", str)
    icon_grid = fortran_namelist.read_variable(grid_namelist,
                                               "icon_grid_nc_file", str)
    icon_grid = utils.clean_path(path_to_grid, icon_grid)
    tg = grid_def.IconGrid(icon_grid)
else:
    logging.error('COSMO grid not supported')
    raise

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= write FORTRAN namelist ===========')
logging.info('')

input_art = fortran_namelist.InputArt()
fortran_namelist.write_fortran_namelist("INPUT_ART", iart, input_art)

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
logging.info('')
logging.info('============= Reading Raw HWSD Land use data ===============')
logging.info('')

hwsd_nc = nc.Dataset(raw_data_art, "r")
raw_lon = hwsd_nc.variables['lon'][:]
raw_lat = hwsd_nc.variables['lat'][:]
raw_lus = hwsd_nc.variables['LU'][:]

lons = np.array(tg.lons)
lats = np.array(tg.lats)

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
logging.info("")
logging.info(
    "============= Mapping raw pixel data to memory map for efficient use ========"
)
logging.info("")

soil_types, neighbor_ids = generate_memory_map(
    raw_lus, soiltype_memmap_filename, nearest_neighbor_index_memmap_filename)

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
logging.info("")
logging.info(
    "============= Find cells in target grid nearest to grid boxes in original HWSD dataset ========"
)
logging.info("")

nrows = np.arange(raw_lat.size)
Parallel(n_jobs=omp, max_nbytes='100M', mmap_mode='w+')(
    delayed(get_neighbor_index)(i, lons, lats, raw_lon, raw_lat, neighbor_ids)
    for i in tqdm(nrows))

# --------------------------------------------------------------------------
logging.info("")
logging.info("============= Calculate LU Fraction for target grid ========")
logging.info("")

fracs = calculate_soil_fraction(tg, soil_types, neighbor_ids, ncpu=2)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= initialize metadata ==============')
logging.info('')

# infer coordinates/dimensions from CDO file
ie_tot = len(tg.lons)
je_tot = 1
ke_tot = 1
lon = np.rad2deg(np.reshape(lons, (ke_tot, je_tot, ie_tot)))
lat = np.rad2deg(np.reshape(lats, (ke_tot, je_tot, ie_tot)))

lat_meta = metadata.Lat()
lon_meta = metadata.Lon()
hcla_meta = metadata.ART_hcla()
silc_meta = metadata.ART_silc()
lcla_meta = metadata.ART_lcla()
sicl_meta = metadata.ART_sicl()
cloa_meta = metadata.ART_cloa()
silt_meta = metadata.ART_silt()
silo_meta = metadata.ART_silo()
scla_meta = metadata.ART_scla()
loam_meta = metadata.ART_loam()
sclo_meta = metadata.ART_sclo()
sloa_meta = metadata.ART_sloa()
lsan_meta = metadata.ART_lsan()
sand_meta = metadata.ART_sand()
udef_meta = metadata.ART_udef()

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= write to buffer file =============')
logging.info('')

buffer_file = buffer.init_netcdf(iart['art_buffer_file'], je_tot, ie_tot)
buffer.write_field_to_buffer(buffer_file, lon, lon_meta)
buffer.write_field_to_buffer(buffer_file, lat, lat_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 0], hcla_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 1], silc_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 2], lcla_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 3], sicl_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 4], cloa_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 5], silt_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 6], silo_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 7], scla_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 8], loam_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 9], sclo_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 10], sloa_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 11], lsan_meta)
buffer.write_field_to_buffer(buffer_file, fracs[:, 12], sand_meta)
buffer.write_field_to_buffer(buffer_file, 1 - fracs.sum(axis=1), udef_meta)
buffer.close_netcdf(buffer_file)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= clean up =========================')
logging.info('')

utils.remove(soiltype_memmap_filename)
utils.remove(nearest_neighbor_index_memmap_filename)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
logging.info('')
logging.info('============= extpar_art_to_buffer done =======')
logging.info('')
