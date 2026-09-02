#!/usr/bin/env python3
import logging
import sys
import numpy as np
import xarray as xr

import horayzon_extpar as hray

# extpar modules from lib
try:
    from extpar.lib import (
        utilities as utils,
        buffer,
        metadata,
        fortran_namelist,
        radtopo,
    )
except ImportError:
    import utilities as utils
    import buffer
    import metadata
    import fortran_namelist
    import radtopo
from namelist import input_radtopo as iradtopo
from namelist import input_oro as ioro

if (not iradtopo["lradtopo"]) or (iradtopo["radtopo_type"] == 1):
    sys.exit()

# initialize logger
logging.basicConfig(filename='extpar_radtopo_to_buffer.log',
                    level=logging.INFO,
                    format='%(message)s',
                    filemode='w')

logging.info('============= start extpar_radtopo_to_buffer ======')
logging.info('')

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

elif (igrid_type != 1):
    error_message = "RADTOPO (HORAYZON) only works with ICON"
    logging.error(error_message)
    raise ValueError(error_message)

#--------------------------------------------------------------------------
# Compute radiation-topography parameters on grid-scale
#--------------------------------------------------------------------------
if iradtopo["radtopo_type"] == 2:

    logging.info('')
    logging.info('==== compute grid-scale radtopo parameters =====')
    logging.info('')

    # Settings
    num_azim_agg = iradtopo["nhori"]
    refine_factor = 10
    num_azim = refine_factor * num_azim_agg
    azim_offset = -(360.0 / (num_azim_agg * 2.0)) + (360 / (num_azim * 2.0))
    # offset of first azimuth position from 0.0 [deg]
    dist_search = float(iradtopo["radius"]) # horizon search distance [m]
    ray_origin_elev = 0.2  # [m]
    horizon_acc = 0.25  # horizon accuracy [deg]
    elev_angle_min = -10.0  # threshold for sampling in negative elevation 
    # angle direction (relevant for 'void regions' at edge of mesh) [deg]
    size_horizon_max = 1.0  # [GB]

    # Load ICON grid
    with xr.open_dataset(icon_grid) as ds:
        clon = ds["clon"].values  # (num_cell) [rad]
        clat = ds["clat"].values  # (num_cell) [rad]
        vlon = ds["vlon"].values  # (num_vertex) [rad]
        vlat = ds["vlat"].values  # (num_vertex) [rad]
        cells_of_vertex = ds["cells_of_vertex"].values - 1  # (6, num_vertex)

    # Load ICON topography (elevation of cell circumcenters)
    with xr.open_dataset("topography_buffer.nc") as ds:
        elevation = ds["HSURF"].values.squeeze()  # (num_cell, float32) [m]
        logging.info(f"DEM source: {ds["HSURF"].data_set}")

    # Generate Embree triangle mesh from ICON grid
    tri_vert, tri_face = radtopo.build_tri_mesh_circ_vert(
        clon,
        clat,
        elevation.astype(np.float64),
        vlon,
        vlat,
        cells_of_vertex,
        neigh_min=6,
    )  # currently only works correctly with 6!
    logging.info(f"Number of ICON triangles: {clon.size}")
    logging.info(f"Number of Embree triangles: {tri_face.shape[0]}")
    if (tri_vert.nbytes / 1e9) > 16.0:
        raise ValueError("Triangle vertices array is larger than 16 GB")

    # Coordinate transformation (lon/lat to ENU)
    radtopo.lonlat2ecef(tri_vert)
    coord_origin = tri_vert.mean(axis=0)
    radius = np.sqrt((coord_origin ** 2).sum())
    lon_origin = np.arctan2(coord_origin[1], coord_origin[0])  # y / x
    lat_origin = np.arcsin(coord_origin[2] / radius)  # z / radius
    # works correctly for ICON domains containing the North/South Pole and/or
    # crossing the +/- 180 deg meridian
    logging.info(f"Origin of ENU system: "
        f"lon = {np.rad2deg(lon_origin):.3f} deg, "
        f"lat = {np.rad2deg(lat_origin):.3f} deg")
    radtopo.ecef2enu(tri_vert, lon_origin=lon_origin, lat_origin=lat_origin)

    # Compute the earth centre and North Pole in ENU coordinates
    earth_centre = np.array(
        [0.0, 0.0, 0.0], dtype=np.float64).reshape((1, 3))
    radtopo.ecef2enu(earth_centre, lon_origin, lat_origin)
    radius_earth = 6_371_229.0
    north_pole = np.array(
        [0.0, 0.0, radius_earth], dtype=np.float64).reshape((1, 3))
    radtopo.ecef2enu(north_pole, lon_origin, lat_origin)

    # Type casting to float32 for Embree
    tri_vert = tri_vert.astype(np.float32)
    earth_centre = earth_centre[0, :].astype(np.float32)
    north_pole = north_pole[0, :].astype(np.float32)

    # Build Embree BVH from terrain triangle mesh
    terrain = hray.Terrain()
    terrain.initialise(
        tri_vert,
        tri_face,
        earth_centre,
        north_pole,
    )

    # Compute horizon in blocks to limit memory allocation
    size_horizon = (clon.size * num_azim * 4) / 1e9  # [GB]
    num_iter = int(np.ceil(size_horizon / size_horizon_max))
    logging.info(f"Number of iterations: {num_iter}")
    slice_locs = np.linspace(0, clon.size, num_iter + 1, dtype=np.uint32)
    horizon_agg = np.empty((clon.size, num_azim_agg), dtype=np.float32)
    sky_view_factor = np.empty(clon.size, dtype=np.float32)
    for idx_iter in range(num_iter):

        # Compute horizon
        slice_loc = slice_locs[idx_iter:(idx_iter + 2)]
        if (slice_loc[0] < 0) or (slice_loc[1] > tri_vert.shape[0]):
            raise ValueError(
                "Indices in 'slice_loc' must be in the range [0, num_vert]"
            )
        # -> move this check later into 'horizon_comp.cpp' and 'horizon.pyx'
        horizon = terrain.horizon_vertex(
            num_azim,
            azim_offset,
            slice_loc,
            dist_search,
            ray_origin_elev,
            horizon_acc,
            elev_angle_min,
        )
        if ((horizon.min() < (elev_angle_min - 2.0 * horizon_acc))
            or (horizon.max() >= 90.0)):
            raise ValueError("Horizon value(s) out of bounds")
        horizon = horizon.clip(min=0.0)

        logging.info(f"Size of horizon array: {horizon.nbytes / 1e9:.2f} GB")

        # Aggregate horizon ('refine_factor' neighbouring values)
        horizon_agg[slice_loc[0]:slice_loc[1], :] = horizon.reshape(
            horizon.shape[0], num_azim_agg, refine_factor).mean(axis=2)

        # Compute the geometric sky view factor
        svf = radtopo.geometric_svf(
            np.deg2rad(horizon),
            scaling=iradtopo["itype_scaling"],
        )
        sky_view_factor[slice_loc[0]:slice_loc[1]] = svf

    del terrain
    del tri_vert, tri_face, earth_centre, north_pole

    # Save terrain horizon and sky view factor to buffer
    logging.info('')
    logging.info('============= write to buffer file =============')
    logging.info('')
    topo_data_set = {1: "GLOBE", 2: "ASTER", 3: "MERIT", 4: "COPERNICUS"}
    ie_tot = clon.size
    lon = np.rad2deg(clon)
    lat = np.rad2deg(clat)
    je_tot = 1
    buffer_file = buffer.init_netcdf(
        iradtopo['radtopo_buffer_file'], je_tot, ie_tot
    )
    buffer_file = buffer.add_dimension_azimuth(buffer_file)
    buffer.write_field_to_buffer(buffer_file, lon, metadata.Lon())
    buffer.write_field_to_buffer(buffer_file, lat, metadata.Lat())
    horizon_meta = metadata.Horizon()
    buffer.write_field_to_buffer(
        buffer_file, horizon_agg.transpose(), horizon_meta
    )
    buffer_file["HORIZON"].data_set = topo_data_set[ioro["itopo_type"]]
    svf_meta = metadata.SVF(scaling=iradtopo["itype_scaling"])
    buffer.write_field_to_buffer(buffer_file, sky_view_factor, svf_meta)
    buffer_file["SKYVIEW"].data_set = topo_data_set[ioro["itopo_type"]]
    buffer.close_netcdf(buffer_file)

#--------------------------------------------------------------------------
# Compute radiation-topography parameters on subgrid-scale
#--------------------------------------------------------------------------
else:

    logging.info('')
    logging.info('=== compute subgrid-scale radtopo parameters ===')
    logging.info('')

    # Get relevant raw topography files -> when ASTER is selected,
    # always use Copernicus DEM !!!

    logging.info('Not yet implemented')

logging.info('')
logging.info('============= extpar_radtopo_to_buffer done =======')
logging.info('')
