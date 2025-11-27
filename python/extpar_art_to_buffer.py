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


def get_neighbor_index(index, hlon, hlat, idxs, ones, balltree):
    points = np.column_stack((hlon, hlat[index] * ones))
    idxs[index, :] = balltree.query(points, k=1)[1].squeeze()


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


def calculate_soil_fraction(target_grid,
                            soil_types_raw,
                            nearest_target_cell_to_raw_cells,
                            ncpu=13):
    """
    target_grid: target ICON grid
    soil_types_raw: landuse class for each cell from the HWSD dataset (LU variable)
    nearest_target_cell_to_raw_cell: indices of the cell from the target ICON grid which is nearest to each cell of the raw grid (from HWSD dataset)
    """
    ncells_target = target_grid.lons.size
    nsoil_types = 13
    nthreads = min(nsoil_types, ncpu)

    soil_ids = np.arange(1, nsoil_types + 1)
    soil_fractions_target = np.zeros((ncells_target, nsoil_types))

    n_nearest_raw_cells = np.bincount(nearest_target_cell_to_raw_cells.ravel(),
                                      minlength=ncells_target)

    def get_fraction_per_soil_type(soil_id):
        n_nearest_raw_cells_with_soil_type = np.bincount(
            nearest_target_cell_to_raw_cells[soil_types_raw == soil_id],
            minlength=ncells_target)

        np.divide(n_nearest_raw_cells_with_soil_type,
                  n_nearest_raw_cells,
                  out=soil_fractions_target[:, soil_id - 1],
                  where=n_nearest_raw_cells != 0)

    Parallel(n_jobs=nthreads,
             max_nbytes='100M',
             mmap_mode='w+',
             backend='threading')(delayed(get_fraction_per_soil_type)(soil_id)
                                  for soil_id in tqdm(soil_ids))

    return soil_fractions_target


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

vlons = np.array(tg.vlons)
vlats = np.array(tg.vlats)

lon_min = max(np.min(vlons), -180.0)
lon_max = min(np.max(vlons), 180.0)
lat_min = max(np.min(vlats), -90.0)
lat_max = min(np.max(vlats), 90.0)

lon_mask = (raw_lon >= lon_min) & (raw_lon <= lon_max)
lat_mask = (raw_lat >= lat_min) & (raw_lat <= lat_max)

raw_lon = raw_lon[lon_mask]
raw_lat = raw_lat[lat_mask]
raw_lus = raw_lus[np.ix_(lat_mask, lon_mask)]

lons = np.deg2rad(lons)
lats = np.deg2rad(lats)
raw_lon = np.deg2rad(raw_lon)
raw_lat = np.deg2rad(raw_lat)

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

lon_lat_array = np.column_stack((lons.ravel(), lats.ravel()))
balltree = BallTree(lon_lat_array, metric="haversine", leaf_size=3)

ones = np.ones((raw_lon.size))

nrows = np.arange(raw_lat.size)
Parallel(n_jobs=omp, max_nbytes='100M', mmap_mode='w+')(delayed(
    get_neighbor_index)(i, raw_lon, raw_lat, neighbor_ids, ones, balltree)
                                                        for i in tqdm(nrows))

# --------------------------------------------------------------------------
logging.info("")
logging.info("============= Calculate LU Fraction for target grid ========")
logging.info("")

fracs = calculate_soil_fraction(tg, soil_types, neighbor_ids, ncpu=omp)

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
