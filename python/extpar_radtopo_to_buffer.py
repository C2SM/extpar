#!/usr/bin/env python3
import logging
from time import perf_counter
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

if (not iradtopo["lradtopo"]) or (iradtopo.get("radtopo_type", 1) == 1):
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

topo_data_set = {1: "GLOBE", 2: "ASTER", 3: "MERIT", 4: "COPERNICUS"}

# Constants
radius_earth = 6_371_229.0  # ICON/COSMO earth radius [m]

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
    orography_buffer_file = ioro.get("orography_buffer_file",
                                     "topography_buffer.nc")
    with xr.open_dataset(orography_buffer_file) as ds:
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

    logging.info('')
    logging.info('============= write to buffer file =============')
    logging.info('')

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

    raw_data_path = ioro.get("raw_data_path", "")

    # Check input topography
    itopo_type = ioro["itopo_type"]
    if itopo_type == 1:
        sys.exit("GLOBE not supported for subgrid-scale radtopo")
    elif itopo_type == 2:
        logging.info(
            'ASTER not supported for subgrid-scale radtopo due to artefacts.\n'
            'Copernicus DEM is used instead.'
        )
        itopo_type = 4

    # Input DEM settings
    if itopo_type == 3: # MERIT
        dem_spacing = 3.0 / 3600.0  # [deg]
        num_tile_lon, num_tile_lat = 12, 6
    else: # COPERNICUS
        dem_spacing = 1.0 / 3600.0  # [deg]
        num_tile_lon, num_tile_lat = 18, 18

    # Resolution of DEM
    earth_circ = 40_075_000.0  # equatorial circumference of Earth [m]
    dem_spacing_m = (earth_circ / 360.0) * dem_spacing  # ~max. spacing [m]
    logging.info(f"DEM resolution: {(dem_spacing_m):.1f} m")
    cell_area_dem = dem_spacing_m ** 2  # [m2]

    # Load ICON grid with specific resolution
    with xr.open_dataset(icon_grid) as ds:
        clon = ds["clon"].values  # (num_cell) [rad]
        clat = ds["clat"].values  # (num_cell) [rad]
        vlon = ds["vlon"].values
        vlat = ds["vlat"].values
        vertex_of_cell = ds["vertex_of_cell"].values - 1  # (3, num_cell)
        # ordered counter-clockwise
        edge_of_cell = ds["edge_of_cell"].values - 1  # (3, num_cell)
        # ordered counter-clockwise
        edge_vertices = ds["edge_vertices"].values - 1  # (2, num_edge)
        grid_level = int(ds.attrs["grid_level"])
        grid_root = int(ds.attrs["grid_root"])
    res_icon = (5050.0 / (grid_root * 2 ** grid_level)) * 1000.0  # [m]
    logging.info(f"ICON resolution: {res_icon:.1f} m")
    cell_area_icon = res_icon ** 2  # [m2]

    # Cartesian coordinates (unit sphere)
    vertices = np.column_stack((
        np.cos(vlat) * np.cos(vlon),
        np.cos(vlat) * np.sin(vlon),
        np.sin(vlat)
    ))

    # Compute refinement level
    logging.info(" Computed refinement level ".center(60, "-"))
    # n ** 2 = cell_area_icon / cell_area_dem -> solve for n
    n_sel = np.sqrt(cell_area_icon / cell_area_dem)
    n_sel = round(n_sel)  # closest to DEM resolution
    logging.info(f"Division steps (n): {n_sel}")
    res_icon_fine = np.sqrt(cell_area_icon / (n_sel ** 2))
    logging.info(f"Refined mesh resolution: {res_icon_fine:.1f} m")
    num_tri_fine = vertex_of_cell.shape[1] * (n_sel ** 2)
    logging.info(
        f"Number of resulting triangles: {num_tri_fine:,}".replace(",", "'")
    )
    logging.info(
        f"Size of 'faces_child' array: {(12 * num_tri_fine / 1e9):.2f} GB"
    )
    logging.info("-" * 60)

    # Refine ICON triangle mesh
    vertices_child, faces_child = radtopo.refine_tri_mesh(
        vertices,
        vertex_of_cell,
        edge_of_cell,
        edge_vertices,
        n_sel,
    )

    # Check size of vertices array
    size_float32 = vertices_child.nbytes / (2.0 * 1e9) # [GB]
    logging.info(
        f"Size of 'vertices_child' array: {size_float32:.2f} GB (32-bit float)"
    )
    if size_float32 > 16.0:
        raise ValueError("Triangle vertices array is larger than 16 GB")

    # Values relevant for child-parent triangle relation
    num_cell_parent = vertex_of_cell.shape[1]
    num_cell_child_per_parent = n_sel ** 2

    # Compute spherical coordinates (longitude/latitude) of child vertices
    t_beg = perf_counter()
    np.arctan2(vertices_child[:, 1], vertices_child[:, 0],
            out=vertices_child[:, 0])
    np.arcsin(vertices_child[:, 2], out=vertices_child[:, 1])
    vertices_child[:, 2] = np.nan  # defined later...
    t_end = perf_counter()
    logging.info(f"Coordinate transformation: {t_end - t_beg:.1f} s")

    # Interpolate elevation data from raw DEM on refined mesh vertices
    logging.info(f"Input DEM: {topo_data_set[itopo_type]}")
    idx_lin_tile = radtopo.assign_points_to_tiles(
        vertices_child[:, 0],
        vertices_child[:, 1],
        num_tile_lon,
        num_tile_lat,
    )
    elevation_interp = np.empty(idx_lin_tile.size, dtype=np.float32)
    num_tile = (num_tile_lon, num_tile_lat)
    for idx_lin in np.unique(idx_lin_tile):

        # Get relevant points for tile
        mask = (idx_lin_tile == idx_lin)
        lon_vert_tile = vertices_child[mask, 0]
        lat_vert_tile = vertices_child[mask, 1]

        # Bounding box for DEM domain
        lon_min = np.rad2deg(lon_vert_tile.min()) - dem_spacing
        lon_max = np.rad2deg(lon_vert_tile.max()) + dem_spacing
        lat_min = np.rad2deg(lat_vert_tile.min()) - dem_spacing
        lat_max = np.rad2deg(lat_vert_tile.max()) + dem_spacing

        # Longitude (i) and latitude (j) tile index
        i = idx_lin % num_tile_lon
        j = idx_lin // num_tile_lon

        # Get a list of required DEM tiles
        i_left = (i + 1) % num_tile_lon
        dem_tiles = [
            radtopo.get_tile_name(i, j, *num_tile),
            radtopo.get_tile_name(i_left, j, *num_tile),
        ]
        if itopo_type == 4:  # COPERNICUS
            if j + 1 < 18:
                j_below = j + 1
                dem_tiles.append(
                    radtopo.get_tile_name(i, j_below, *num_tile)
                )
                dem_tiles.append(
                    radtopo.get_tile_name(i_left, j_below, *num_tile)
                )
            lon_slice = slice(0, 3_600 * 20 + 1)
            lat_slice = slice(0, 3_600 * 10 + 1)
            dem_tiles = [f"COPERNICUS_{tile}.nc" for tile in dem_tiles]
            var_elevation = "elevation"
        else: # MERIT
            if j - 1 >= 0:
                j_above = j - 1
                dem_tiles.append(
                    radtopo.get_tile_name(i, j_above, *num_tile)
                )
                dem_tiles.append(
                    radtopo.get_tile_name(i_left, j_above, *num_tile)
                )
            lon_slice = slice(0, 1_200 * 30 + 1)
            lat_slice = slice(1_200 * 30 - 1, 2 * 1_200 * 30)
            dem_tiles = [f"MERIT_{tile}.nc" if tile[:3] != "S60"
                         else f"REMA_BKG_{tile}.nc" for tile in dem_tiles]
            var_elevation = "Elevation"
        dem_tiles = [
            utils.clean_path(raw_data_path, tile) for tile in dem_tiles
        ]
        logging.info("\n".join(dem_tiles))

        # Load DEM data
        with xr.open_mfdataset(
            [tile for tile in dem_tiles], mask_and_scale=False
        ) as ds:
            ds = ds.isel(lon=lon_slice, lat=lat_slice)
            ds = ds.sel(
                lon=slice(lon_min, lon_max),
                lat=slice(lat_max, lat_min)
            )
            lon_dem = np.deg2rad(ds["lon"].values)  # [rad]
            lat_dem = np.deg2rad(ds["lat"].values)  # [rad]
            elevation_dem = ds[var_elevation].values  # [m]
            if itopo_type == 3:  # MERIT
                elevation_dem[elevation_dem == -32767] = 0
            elevation_dem = elevation_dem.astype(np.float32)
        logging.info(f"Elevation range: {elevation_dem.min():.1f} "
                     f"– {elevation_dem.max():.1f} m")

        # Interpolate elevation bilinearly from DEM to triangle mesh vertices
        if ((lon_vert_tile.min() < lon_dem.min())
            or (lon_vert_tile.max() > lon_dem.max())
            or (lat_vert_tile.min() < lat_dem.min())
            or (lat_vert_tile.max() > lat_dem.max())):
            raise ValueError("Interpolation point(s) outside of source grid")
            # -> this check has to be removed for the North/South Pole because
            #    there, extrapolation is required (-> clamping to grid)
        atol_grid = np.deg2rad(1.0e-8)  # ca. 1 mm (max) for degree
        elevation_interp[mask] = radtopo.interp_bilinear(
            elevation_dem,
            lon_dem,
            lat_dem,
            lon_vert_tile,
            lat_vert_tile,
            atol_grid,
        )
    vertices_child[:, 2] = elevation_interp
    del lon_dem, lat_dem, elevation_dem
    del lon_vert_tile, lat_vert_tile, elevation_interp

    # Coordinate transformation (lon/lat to ENU)
    radtopo.lonlat2ecef(vertices_child)
    coord_origin = vertices_child.mean(axis=0)
    radius = np.sqrt((coord_origin ** 2).sum())
    lon_origin = np.arctan2(coord_origin[1], coord_origin[0])  # y / x
    lat_origin = np.arcsin(coord_origin[2] / radius)  # z / radius
    # works correctly for ICON domains containing the North/South Pole and/or
    # crossing the +/- 180 deg meridian
    logging.info(f"Origin of ENU system: "
                 f"lon = {np.rad2deg(lon_origin):.3f} deg, "
                 f"lat = {np.rad2deg(lat_origin):.3f} deg")
    radtopo.ecef2enu(
        vertices_child, lon_origin=lon_origin, lat_origin=lat_origin
    )

    # Compute the earth centre and North Pole in ENU coordinates
    earth_centre = np.array(
        [0.0, 0.0, 0.0], dtype=np.float64).reshape((1, 3))
    radtopo.ecef2enu(
        earth_centre, lon_origin=lon_origin, lat_origin=lat_origin
    )
    north_pole = np.array(
        [0.0, 0.0, radius_earth], dtype=np.float64).reshape((1, 3))
    radtopo.ecef2enu(north_pole, lon_origin=lon_origin, lat_origin=lat_origin)

    # Type casting to float32 for Embree
    vertices_child = vertices_child.astype(np.float32)
    earth_centre = earth_centre[0, :].astype(np.float32)
    north_pole = north_pole[0, :].astype(np.float32)

    # Build Embree BVH from terrain triangle mesh
    terrain = hray.Terrain()
    terrain.initialise(
        vertices_child,
        faces_child,
        earth_centre,
        north_pole,
        )

    # Settings
    num_azim = iradtopo["nhori"]
    dist_search = float(iradtopo["radius"])  # horizon search distance [m]
    ray_origin_elev = 0.2  # 0.1, 0.2 [m]
    horizon_acc = 0.25  # horizon accuracy [deg]
    elev_angle_min = -10.0  # -85.0
    size_horizon_max = 2.0  # [GB] (0.5, 1.0, 2.0)
    num_elev = 181  # 91, 181
    sw_dir_cor_max = 25.0  # maximum for individual values
    sw_dir_cor_agg_max = 10.0  # maximum for aggregated values
    num_nodes = 7
    eta_sel = 2.0  # set hard-coded value in ICON code accordingly

    # Compute block size and number of iterations
    size_horizon_ppt = (num_cell_child_per_parent * num_azim * 4) / 1e9
    # horizon array size per parent triangle [GB]
    num_ptpi = int(size_horizon_max / size_horizon_ppt)
    # number of parent triangles per iteration
    if num_ptpi == 0:
        raise ValueError(
            "'size_horizon_max' too small for a single parent triangle"
        )
    num_iter = int(np.ceil(num_cell_parent / num_ptpi))
    logging.info(f"Number of iterations: {num_iter}")
    slice_locs_child = np.linspace(
        0, num_cell_parent, num_iter + 1, dtype=np.uint32
    ) * num_cell_child_per_parent

    # Allocate output arrays
    svf_agg_all = np.empty(num_cell_parent, dtype=np.float32)
    sw_dir_cor_sparse_all = np.empty((num_cell_parent, num_azim, num_nodes),
                                     dtype=np.float32)
    terrain_normal_all = np.empty((num_cell_parent, 3), dtype=np.float32)
    elev_shadow_all = np.empty((num_cell_parent, num_azim, 3),
                               dtype=np.float32)

    # Loop
    for idx_iter in range(num_iter):

        logging.info(f" Iteration {idx_iter + 1}/{num_iter} ".center(60, "-"))

        # Compute horizon
        slice_loc_child = slice_locs_child[idx_iter:(idx_iter + 2)]
        if ((slice_loc_child[0] < 0)
            or (slice_loc_child[1] > faces_child.shape[0])):
            raise ValueError(
                "Indices in 'slice_loc' must be in the range [0, num_face]"
            )
        # -> move this check later into 'horizon_comp.cpp' and 'horizon.pyx'
        horizon = terrain.horizon_centroid(
            num_azim,
            slice_loc_child,
            dist_search,
            ray_origin_elev,
            horizon_acc,
            elev_angle_min,
        )
        if ((horizon.min() < (elev_angle_min - 2.0 * horizon_acc))
            or (horizon.max() >= 90.0)):
            raise ValueError("Horizon value(s) out of bounds")
        horizon = horizon.clip(min=0.0)

        # Compute the geometric sky view factor (spatially aggregated)
        svf = radtopo.geometric_svf(np.deg2rad(horizon), scaling=2)
        slice_loc_parent = (
            slice_loc_child / num_cell_child_per_parent
        ).astype(int)
        svf_agg = np.nanmean(
            svf.reshape(-1, num_cell_child_per_parent), axis=1
        )
        svf_agg_all[slice_loc_parent[0]:slice_loc_parent[1]] = svf_agg

        elev = np.linspace(0.0, 90.0, num_elev)  # [deg]
        faces_child_sel = faces_child[slice_loc_child[0]:slice_loc_child[1], :]

        radtopo.process_block(
            num_cell_child_per_parent,
            num_nodes,
            sw_dir_cor_max,
            sw_dir_cor_agg_max,
            eta_sel,
            slice_loc_parent,
            elev,
            vertices_child,
            faces_child_sel,
            earth_centre,
            north_pole,
            horizon,
            terrain_normal_all,
            elev_shadow_all,
            sw_dir_cor_sparse_all,
        )

    del terrain
    del vertices_child, faces_child, earth_centre, north_pole

    logging.info('')
    logging.info('============= write to buffer file =============')
    logging.info('')

    ie_tot = clon.size
    lon = np.rad2deg(clon)
    lat = np.rad2deg(clat)
    je_tot = 1
    buffer_file = buffer.init_netcdf(
        iradtopo['radtopo_buffer_file'], je_tot, ie_tot
    )
    buffer_file = buffer.add_dimension_azimuth(buffer_file)
    buffer_file = buffer.add_dimension_vector_component(buffer_file)
    buffer_file = buffer.add_dimension_element(buffer_file)
    buffer.write_field_to_buffer(buffer_file, lon, metadata.Lon())
    buffer.write_field_to_buffer(buffer_file, lat, metadata.Lat())
    shp = (num_cell_parent, num_azim * 3)
    horizon = elev_shadow_all.reshape(shp).transpose()
    horizon_meta = metadata.Horizon()
    buffer.write_field_to_buffer(buffer_file, horizon, horizon_meta)
    buffer_file["HORIZON"].data_set = topo_data_set[itopo_type]
    svf_meta = metadata.SVF(scaling=iradtopo["itype_scaling"])
    buffer.write_field_to_buffer(buffer_file, svf_agg_all, svf_meta)
    buffer_file["SKYVIEW"].data_set = topo_data_set[itopo_type]
    shp = (num_cell_parent, num_azim * (num_nodes - 2))
    swdir_cor = sw_dir_cor_sparse_all[:, :, 1:-1].reshape(shp).transpose()
    swdir_cor_meta = metadata.SWDIR_COR()
    buffer.write_field_to_buffer(buffer_file, swdir_cor, swdir_cor_meta)
    buffer_file["SWDIR_COR"].data_set = topo_data_set[itopo_type]
    terrain_normal_meta = metadata.Terrain_normal()
    buffer.write_field_to_buffer(buffer_file, terrain_normal_all.transpose(),
                                 terrain_normal_meta)
    buffer_file["TERRAIN_NORMAL"].data_set = topo_data_set[itopo_type]
    buffer.close_netcdf(buffer_file)

logging.info('')
logging.info('============= extpar_radtopo_to_buffer done =======')
logging.info('')
