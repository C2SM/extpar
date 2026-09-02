import functools
import time

import numpy as np
from numba import njit, prange
from numba import float32, float64, int32, int64, uint32, void
from numba import types # type: ignore
'''
Module utilities with auxiliary functions to compute grid- and subgrid-scale
radiation-topography correction parameters.
'''

# -----------------------------------------------------------------------------
# Shared (grid- and subgrid-scale)
# -----------------------------------------------------------------------------

def measure_time(func):
    """
    Decorator to measure the execution time of a function.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        time_start = time.perf_counter()
        result = func(*args, **kwargs)
        time_end = time.perf_counter()
        print(f"{func.__name__}: {time_end - time_start:.1f} s")
        return result
    return wrapper


@njit(void(
    float64[:, :],
), parallel=True, cache=True)
def lonlat2ecef(
    coord,
):
    """
    Transform spherical longitude/latitude coordinates to earth-centered,
    earth-fixed (ECEF) coordinates in-place.

    Parameters
    ----------
    coord : ndarray of float64 (num_points, 3)
        Input: longitude [rad], latitude [rad], elevation above surface [m]
        Output: ECEF x, y and z-coordinates [m]

    Notes
    -----
    - This and the below function could be improved by accepting an xarray
      DataArray instead of a numpy.ndarray, which allows to specify the
      coordinate system.
    """
    radius_earth = 6_371_229.0  # ICON/COSMO earth radius [m]
    for i in prange(coord.shape[0]):
        lon = coord[i, 0]
        lat = coord[i, 1]
        radius = radius_earth + coord[i, 2]
        coord[i, 0] = radius * np.cos(lat) * np.cos(lon)
        coord[i, 1] = radius * np.cos(lat) * np.sin(lon)
        coord[i, 2] = radius * np.sin(lat)


@njit(void(
    float64[:, :],
    float64,
    float64,
), parallel=True, cache=True)
def ecef2enu(
    coord,
    lon_origin,
    lat_origin,
):
    """
    Transform earth-centered, earth-fixed (ECEF) coordinates to local
    east-north-up (ENU) coordinates in-place. The origin of ENU coordinate
    system is bound to the surface of the sphere.

    Parameters
    ----------
    coord : ndarray of float64 (num_points, 3)
        Input: ECEF x, y and z-coordinates [m]
        Output: ENU x, y and z-coordinates [m]
    lon_origin : float64
        Longitude of the ENU coordinate system's origin [rad]
    lat_origin : float64
        Latitude of the ENU coordinate system's origin [rad]
    """
    radius_earth = 6_371_229.0  # ICON/COSMO earth radius [m]
    x_origin = radius_earth * np.cos(lat_origin) * np.cos(lon_origin)
    y_origin = radius_earth * np.cos(lat_origin) * np.sin(lon_origin)
    z_origin = radius_earth * np.sin(lat_origin)
    sin_lon_origin = np.sin(lon_origin)
    cos_lon_origin = np.cos(lon_origin)
    sin_lat_origin = np.sin(lat_origin)
    cos_lat_origin = np.cos(lat_origin)
    for i in prange(coord.shape[0]):
        dx = coord[i, 0] - x_origin
        dy = coord[i, 1] - y_origin
        dz = coord[i, 2] - z_origin
        coord[i, 0] = -sin_lon_origin * dx \
            + cos_lon_origin * dy
        coord[i, 1] = -sin_lat_origin * cos_lon_origin * dx \
            - sin_lat_origin * sin_lon_origin * dy \
            + cos_lat_origin * dz
        coord[i, 2] = cos_lat_origin * cos_lon_origin * dx \
            + cos_lat_origin * sin_lon_origin * dy \
            + sin_lat_origin * dz


@njit(float32[:](
    float32[:, :],
    int64,
), parallel=True, cache=True)
def geometric_svf(
    horizon,
    scaling,
):
    """
    Compute the geometric sky view factor (SVF) from the terrain horizon.

    Parameters
    ----------
    horizon : ndarray of float32 (num_cell, num_azim)
        Terrain horizon angles [rad]
    scaling : int64
        Scaling exponent; (scaling + 1) is applied to sin(horizon)

    Returns
    -------
    svf : ndarray of float32 (num_cell)
        Geometric sky view factor for each grid cell [-]
    """
    if scaling < 0 or scaling > 2:
        raise ValueError("Scaling must be 0, 1 or 2")
    num_cell, num_azim = horizon.shape
    svf = np.empty(num_cell, dtype=horizon.dtype)
    for i in prange(num_cell):
        s = 0.0
        for j in range(num_azim):
            x = np.sin(horizon[i, j])
            s += 1.0 - x ** (scaling + 1)
        svf[i] = s / num_azim
    return svf


# -----------------------------------------------------------------------------
# Grid-scale
# -----------------------------------------------------------------------------

@njit(types.Tuple((
    float64[:, :],
    uint32[:, :],
))(
    float64[:],
    float64[:],
    float64[:],
    float64[:],
    float64[:],
    int32[:, :],
    int64,
), cache=True)
def build_tri_mesh_circ_vert(
    lon_circ,
    lat_circ,
    elevation_circ,
    lon_vert,
    lat_vert,
    cells_of_vertex,
    neigh_min=6,
):
    """
    Build a triangle mesh from the ICON grid cell circumcenters and vertices.

    Parameters
    ----------
    lon_circ : ndarray of float64 (num_cell)
        Longitude of ICON grid cell circumcenters [rad]
    lat_circ : ndarray of float64 (num_cell)
        Latitude of ICON grid cell circumcenters [rad]
    elevation_circ : ndarray of float64 (num_cell)
        Elevation of ICON grid cell circumcenters [m]
    lon_vert : ndarray of float64 (num_vert)
        Longitude of ICON grid cell vertices [rad]
    lat_vert : ndarray of float64 (num_vert)
        Latitude of ICON grid cell vertices [rad]
    cells_of_vertex : ndarray of int32 (6, num_vert)
        Indices of the (up to six) ICON grid cells adjoining each grid
        vertex, padded with -2 for vertices with fewer than six adjoining
        cells
    neigh_min : int64
        Minimum number of neighboring ICON cells to include the polygon for
        triangulation (default: 4)

    Returns
    -------
    tri_vert : ndarray of float64 (num_cell, 3)
        Longitude, latitude and elevation of the Embree triangle mesh vertices
        [rad, rad, m]
    tri_face : ndarray of uint32 (num_face, 3)
        Indices of the Embree triangle mesh faces

    Notes
    -----
    - limitation of algorithm: works only correctly for 'neigh_min = 6' and as
      long as the ICON base icosahedron vertices, which only have 5 ICON cell
      neighbours, are not in the domain. For smaller 'neigh_min', the ICON
      grid cells around a centre vertex do not form a closed cycle (except at
      the vertices of the base icosahedron) and correct ordering would be
      necessary
   - the above limitation could be resolved by ordering the ICON circumcenters
      based on walking through the cells utilising 'neighbor_cell_index'
    """
    if (neigh_min < 2) or (neigh_min > 6):
        raise ValueError("Argument 'neigh_min' must be in the range [2, 6]")
    num_cell = lon_circ.size  # number of ICON grid cells
    num_vert = lat_vert.size  # number of ICON grid cell vertices
    tri_vert = np.empty((num_cell + num_vert, 3), dtype=np.float64)
    # Embree buffer (allocate maximal possible size)
    for i in range(num_cell):
        tri_vert[i, 0] = lon_circ[i]
        tri_vert[i, 1] = lat_circ[i]
        tri_vert[i, 2] = elevation_circ[i]
    num_tri_vert = num_cell
    idx_add = num_cell
    angles = np.empty(6, dtype=np.float64)
    tri_face = np.empty((num_vert * 6, 3), dtype=np.uint32)  # Embree buffer
    # (allocate maximal possible size; 6 triangles per ICON grid vertex)
    idx_face = 0
    for idx_vert in range(num_vert):
        elevation_mean = 0.0
        num_angle = 0
        for j in range(6):
            idx_cell = cells_of_vertex[j, idx_vert]
            if idx_cell != -2:
                angle = np.arctan2(lat_circ[idx_cell] - lat_vert[idx_vert],
                                   lon_circ[idx_cell] - lon_vert[idx_vert])
                # anti-clockwise angle from positive x-axis
                if angle < 0.0:
                    angle += 2.0 * np.pi
                angles[num_angle] = angle
                num_angle += 1
                elevation_mean += elevation_circ[idx_cell]
        if num_angle >= neigh_min:
            tri_vert[num_tri_vert, 0] = lon_vert[idx_vert]
            tri_vert[num_tri_vert, 1] = lat_vert[idx_vert]
            tri_vert[num_tri_vert, 2] = elevation_mean / float(num_angle)
            num_tri_vert += 1
            idx_sort = np.argsort(angles[:num_angle])
            # num_iter = num_angle
            num_iter = 6 if (num_angle == 6) else num_angle - 1
            for j in range(num_iter):
                tri_face[idx_face, 0] \
                    = cells_of_vertex[idx_sort[j], idx_vert]
                tri_face[idx_face, 1] \
                    = cells_of_vertex[idx_sort[(j + 1) % num_angle], idx_vert]
                tri_face[idx_face, 2] \
                    = idx_add
                idx_face += 1
            idx_add += 1
    tri_vert = tri_vert[:num_tri_vert, :]
    tri_face = tri_face[:idx_face, :]
    return tri_vert, tri_face


# -----------------------------------------------------------------------------
# Subgrid-scale
# -----------------------------------------------------------------------------

@measure_time
@njit(types.Tuple((
    float64[:, :],
    uint32[:, :],
))(
    float64[:, :],
    int32[:, :],
    int32[:, :],
    int32[:, :],
    int32,
), cache=True)
def refine_tri_mesh(
    vertices,
    vertex_of_cell,
    edge_of_cell,
    edge_vertices,
    n,
):
    """
    Refine a triangle mesh on a sphere surface by subdividing each parent
    triangle n ** 2 child triangles.

    Parameters
    ----------
    vertices : ndarray of float64 (num_vert, 3)
        Cartesian coordinates (x, y, z) of vertices [-] (unit sphere)
    vertex_of_cell : ndarray of int32 (3, num_cells)
        Indices of triangle faces
    edge_of_cell : ndarray of int32 (3, num_cells)
        Indices of triangle edges
    edge_vertices : ndarray of int32 (2, num_edges)
        Indices of vertices of edges
    n : int32
        Triangle division number

    Returns
    -------
    vertices_child : ndarray of float64 (num_vertex_child, 3)
        Cartesian coordinates (x, y, z) of refined vertices [-] (unit sphere)
    faces_child : ndarray of uint32 (num_faces_child, 3)
        Indices of refined triangle faces. The winding order (orientation) is
        the same as in the input (vertex_of_cell)
    """

    # Compute number of new vertices
    num_vertex_in = vertices.shape[0]
    num_vertex_edge = edge_vertices.shape[1] * (n - 1)
    num_vertex_interior_pgc = 0  # number of interior vertices per grid cell
    for i in range(n - 1):
        num_vertex_interior_pgc += i
    num_vertex_interior = vertex_of_cell.shape[1] * num_vertex_interior_pgc
    num_vertex_child = num_vertex_in + num_vertex_edge + num_vertex_interior

    # Mapping of triangle vertex indices
    num_vert_per_tri = 3 + 3 * (n - 1) + num_vertex_interior_pgc
    # print(f"Number of vertices per triangle: {num_vert_per_tri}")
    mapping = np.empty(num_vert_per_tri, dtype=np.uint32)
    idx = 0
    mapping[idx] = 0
    idx += 1
    mapping[idx:(idx + (n - 1))] \
        = np.arange(3 + 2 * (n - 1), 3 + 2 * (n - 1) + (n - 1),
                    dtype=np.uint32)
    idx += (n - 1)
    mapping[idx] = 2
    idx += 1
    start = n * 3
    for i in range(0, n - 1):
        mapping[idx] = 3 + i
        idx += 1
        mapping[idx:(idx + n - 2 - i)] \
            = np.arange(start, start + n - 2 - i, dtype=np.uint32)
        start = start + n - 2 - i
        idx += (n - 2 - i)
        mapping[idx] = 3 + (n - 1) + i
        idx += 1
    mapping[idx] = 1
    idx += 1
    if not np.all(np.diff(np.unique(mapping)) == 1) \
        or (num_vert_per_tri != idx):
        raise ValueError("Error while computing the remapping array")

    # Compute faces
    index_2d = np.empty((n + 1, n + 1), dtype=np.uint32)
    index_2d.fill(-999)
    idx_vertex = 0
    for i in range(n + 1):
        for j in range(n + 1 - i):
            index_2d[i, j] = idx_vertex
            idx_vertex += 1
    faces = np.empty((n ** 2, 3), dtype=np.uint32)
    idx_face = 0
    for i in range(n):
        for j in range(n - i):
            idx_v0 = index_2d[i, j]
            idx_v1 = index_2d[i + 1, j]
            idx_v2 = index_2d[i, j + 1]
            faces[idx_face, :] = (idx_v0, idx_v1, idx_v2)
            idx_face += 1
            if i + j + 1 < n:
                idx_v3 = index_2d[i + 1, j + 1]
                faces[idx_face, :] = (idx_v1, idx_v3, idx_v2)
                idx_face += 1

    # Allocate arrays for refined mesh
    vertices_child = np.empty((num_vertex_child, 3), dtype=np.float64)
    faces_child = np.empty((n ** 2 * vertex_of_cell.shape[1], 3),
                           dtype=np.uint32)

    # Add vertices from base mesh
    vertices_child[:num_vertex_in, :] = vertices

    # Add vertices located on the edge of bash mesh (shared)
    t = np.linspace(0.0, 1.0, num=(n + 1))[1:-1]
    idx_vertex = num_vertex_in
    for i in range(edge_vertices.shape[1]):  # loop through all edges
        vertex_0 = vertices[edge_vertices[0, i], :]
        vertex_1 = vertices[edge_vertices[1, i], :]
        for j in range(n - 1):
            vertices_child[idx_vertex, :] \
                = vertex_0 + t[j] * (vertex_1 - vertex_0)
            idx_vertex += 1

    # Add vertices located in the interior of base mesh triangles
    for idx_cell in range(vertex_of_cell.shape[1]): # loop through all cells
        vertex_0 = vertices[vertex_of_cell[0, idx_cell]]
        vertex_1 = vertices[vertex_of_cell[1, idx_cell]]
        vertex_2 = vertices[vertex_of_cell[2, idx_cell]]
        for i in range(1, n):
            for j in range(1, n - i):
                k = n - i - j
                vertices_child[idx_vertex, :] \
                    = (k * vertex_0 + i * vertex_1 + j * vertex_2) / n
                idx_vertex += 1
    for i in range(vertices_child.shape[0]):
        vertices_child[i, :] /= np.linalg.norm(vertices_child[i, :])
    # unit vectors

    # Connect vertices into child triangle
    indices = np.empty(num_vert_per_tri, dtype=np.uint32)
    idx_face = 0
    for idx_cell in range(vertex_of_cell.shape[1]):
    # for ind_cell in range(500):
        indices[:3] = vertex_of_cell[:, idx_cell]
        # counter-clockwise ordered
        idx_vertex = 0
        for idx in edge_of_cell[:, idx_cell]:
            indices_edge = np.arange(num_vertex_in + idx * (n - 1),
                                    num_vertex_in + (idx + 1) * (n - 1))
            # ordering (clockwise vs. counter-clockwise) not consistent
            if idx_vertex == 0: # order: 0 -> 1
                if vertex_of_cell[idx_vertex, idx_cell] \
                    != edge_vertices[0, edge_of_cell[idx_vertex, idx_cell]]:
                    indices_edge = indices_edge[::-1]
            else: # order: 2 -> 1, 0 -> 2
                if vertex_of_cell[idx_vertex, idx_cell] \
                    == edge_vertices[0, edge_of_cell[idx_vertex, idx_cell]]:
                    indices_edge = indices_edge[::-1]
            slice_v = slice(3 + idx_vertex * (n - 1),
                            3 + (idx_vertex + 1) * (n - 1))
            indices[slice_v] = indices_edge
            idx_vertex += 1
        indices[slice_v.stop:] \
            = np.arange(num_vertex_in + num_vertex_edge
                        + idx_cell * num_vertex_interior_pgc,
                        num_vertex_in + num_vertex_edge
                        + (idx_cell + 1) * num_vertex_interior_pgc)
        for i in range(n ** 2):
            faces_child[idx_face, :] = indices[mapping[faces[i, :]]]
            idx_face += 1

    return vertices_child, faces_child


@measure_time
def assign_points_to_tiles(
    lon,
    lat,
    num_tile_lon,
    num_tile_lat,
):
    """
    Assign points to tiles covering the entire surface of the Earth. Tile
    indices increase eastward from -180 deg longitude and southward from +90
    deg latitude.

    Parameters
    ----------
    lon : ndarray of float64 (num_points)
        Longitude coordinates of points [rad]
    lat : ndarray of float64 (num_points)
        Latitude coordinates of points [rad]
    num_tile_lon : int
        Number of tiles along the longitudinal direction
    num_tile_lat : int
        Number of tiles along the latitudinal direction

    Returns
    -------
    idx_lin_tile : ndarray of int16 (num_points)
        Linear index of tile associated with point
    """

    # Check validity of input arguments
    if lon.size != lat.size:
        raise ValueError("'lon' and 'lat' must have equal size")
    if (lon.min() < -np.pi) or (lon.max() > +np.pi):
        raise ValueError("Longitude value(s) out of range [-pi, +pi]")
    if (lat.min() < -(np.pi / 2.0)) or (lat.max() > +(np.pi / 2.0)):
        raise ValueError("Latitude value(s) out of range [-pi/2, +pi/2]")
    if (num_tile_lon <= 0) or (360 % num_tile_lon != 0):
        raise ValueError("Invalid number of longitudinal tiles")
    if (num_tile_lat <= 0) or (180 % num_tile_lat != 0):
        raise ValueError("Invalid number of latitudinal tiles")

    tile_extent_lon = np.deg2rad(360.0 / num_tile_lon)
    tile_extent_lat = np.deg2rad(180.0 / num_tile_lat)
    idx_tile_lon = np.floor((lon + np.pi) / tile_extent_lon).astype(np.int16)
    idx_tile_lat = np.floor(
        ((np.pi / 2.0) - lat) / tile_extent_lat
    ).astype(np.int16)
    idx_tile_lon = np.clip(idx_tile_lon, 0, num_tile_lon - 1)
    idx_tile_lat = np.clip(idx_tile_lat, 0, num_tile_lat - 1)

    return idx_tile_lat * num_tile_lon + idx_tile_lon


def get_tile_name(
        idx_tile_lon,
        idx_tile_lat,
        num_tile_lon,
        num_tile_lat,
):
    """
    Return coordinates part of tile name based on longitudinal and latitudinal
    tile indices. Tile indices increase eastward from -180 deg longitude and
    southward from +90 deg latitude.

    Parameters
    ----------
    idx_tile_lon : int
        Longitudinal tile index
    idx_tile_lat : int
        Latitudinal tile index
    num_tile_lon : int
        Number of tiles along the longitudinal direction
    num_tile_lat : int
        Number of tiles along the latitudinal direction

    Returns
    -------
    tile_name_coord : str
        Coordinates part of tile name
    """

    # Check validity of input arguments
    if (num_tile_lon <= 0) or (360 % num_tile_lon != 0):
        raise ValueError("Invalid number of longitudinal tiles")
    if (num_tile_lat <= 0) or (180 % num_tile_lat != 0):
        raise ValueError("Invalid number of latitudinal tiles")
    if (idx_tile_lon < 0) or (idx_tile_lon >= num_tile_lon):
        raise ValueError("Longitudinal tile index out of bounds")
    if (idx_tile_lat < 0) or (idx_tile_lat >= num_tile_lat):
        raise ValueError("Latitudinal tile index out of bounds")

    tile_extent_lon = 360 // num_tile_lon
    lon_west = -180 + idx_tile_lon * tile_extent_lon
    lon_east = lon_west + tile_extent_lon
    letter_west = "E" if lon_west >= 0 else "W"
    letter_east = "E" if lon_east >= 0 else "W"

    tile_extent_lat = 180 // num_tile_lat
    lat_north = 90 - idx_tile_lat * tile_extent_lat
    lat_south = lat_north - tile_extent_lat
    letter_north = "N" if lat_north >= 0 else "S"
    letter_south = "N" if lat_south >= 0 else "S"

    tile_name_coord = (
        f"{letter_north}{abs(lat_north):02d}-"
        f"{letter_south}{abs(lat_south):02d}_"
        f"{letter_west}{abs(lon_west):03d}-"
        f"{letter_east}{abs(lon_east):03d}"
    )

    return tile_name_coord


@measure_time
@njit(float32[:](
    float32[:, :],
    float64[::1],
    float64[::1],
    float64[::1],
    float64[::1],
    float64,
), parallel=True, cache=True)
def interp_bilinear(
    data,
    x_axis,
    y_axis,
    x_interp,
    y_interp,
    atol_grid,
):
    """
    Bilinear interpolation on a regular grid. Interpolation points outside the
    grid domain are clamped to the nearest grid boundary.

    Parameters
    ----------
    data : ndarray of float32 (len_y, len_x)
        Input data
    x_axis : ndarray of float64 (len_x)
        Regularly spaced x-coordinates of the input data
    y_axis : ndarray of float64 (len_y)
        Regularly spaced y-coordinates of the input data
    x_interp : ndarray of float64 (num_points)
        x-coordinates for bilinear interpolation
    y_interp : ndarray of float64 (num_points)
        y-coordinates for bilinear interpolation
    atol_grid : float64
        Absolute tolerance used to check whether the grid spacing along the
        x- and y-axes is regular

    Returns
    -------
    data_interp : ndarray of float32 (num_points)
        Bilinearly interpolated data
    """

    # Check validity of input arguments
    len_y, len_x = data.shape
    if (len_x != x_axis.size) or (len_y != y_axis.size):
        raise ValueError("Shape of 'data' does not match coordinate axes")
    if x_interp.size != y_interp.size:
        raise ValueError("'x_interp' and 'y_interp' must have equal size")
    if (len_x < 2) or (len_y < 2):
        raise ValueError("Coordinate axes must contain at least two points")
    if atol_grid < 0.0:
        raise ValueError("'atol_grid' must be non-negative")
    delta_x = (x_axis[-1] - x_axis[0]) / (len_x - 1)
    delta_y = (y_axis[-1] - y_axis[0]) / (len_y - 1)
    if np.max(np.abs(np.diff(x_axis) - delta_x)) > atol_grid:
        raise ValueError("x-axis must be regularly spaced")
    if np.max(np.abs(np.diff(y_axis) - delta_y)) > atol_grid:
        raise ValueError("y-axis must be regularly spaced")

    # Loop through interpolation points
    data_interp = np.empty(x_interp.size, dtype=np.float32)

    for k in prange(x_interp.size):

        x = (x_interp[k] - x_axis[0]) / delta_x
        y = (y_interp[k] - y_axis[0]) / delta_y
        x = min(max(x, 0.0), len_x - 1.0)  # clamp x
        y = min(max(y, 0.0), len_y - 1.0)  # clamp y
        i_0 = min(int(x), len_x - 2)
        j_0 = min(int(y), len_y - 2)
        i_1 = i_0 + 1
        j_1 = j_0 + 1
        weight_x = x - i_0
        weight_y = y - j_0
        data_interp[k] = (
            (1.0 - weight_x) * (1.0 - weight_y) * data[j_0, i_0]
            + weight_x * (1.0 - weight_y) * data[j_0, i_1]
            + (1.0 - weight_x) * weight_y * data[j_1, i_0]
            + weight_x * weight_y * data[j_1, i_1]
        )

    return data_interp


@njit(types.Tuple((
    float64[:, :],
    uint32[:, :],
    float64[:],
))(
    float32[:, :],
    uint32[:, :],
    float32[:],
    float32[:],
    float32[:, :],
    int32,
    int32,
    int32,
    float64,
), cache=True)
def compute_sw_dir_cor_agg(
    tri_vert,
    tri_face,
    earth_centre,
    north_pole,
    horizon,
    idx_tri_start,
    idx_tri_stop,
    num_elev,
    sw_dir_cor_max,
):
    """
    Compute direct shortwave correction factors and fractional illumination
    for a hemisphere, aggregated over multiple triangles. Additionally, the
    averaged terrain normal is returned.

    Parameters
    ----------
    tri_vert : ndarray of float32 (num_vert, 3)
        Cartesian coordinates of triangle mesh vertices [m]
    tri_face : ndarray of uint32 (num_face, 3)
        Indices of the triangle mesh faces
    earth_centre : ndarray of float32 (3)
        Earth centre in cartesian coordinates [m]
    north_pole : ndarray of float32 (3)
        North Pole in cartesian coordinates [m]
    horizon : ndarray of float32 (num_face, num_azim)
        Terrain horizon of individual triangles [deg]
    idx_tri_start : int32 (num_points, 3)
        Start indices of considered triangles
    idx_tri_stop : int32(num_points, 3)
        Stop indices of considered triangles
    num_elev : int32
        Number of elevation angles for correction factor computation
    sw_dir_cor_max : float64
        Maximum value for individual correction factor [-]

    Returns
    -------
    sw_dir_cor : ndarray of float64 (num_azim, num_elev)
        Averaged direct shortwave correction factor [-]
    illuminated : ndarray of uint32 (num_azim, num_elev)
        Illumination fraction [-]
    terrain_normal : ndarray of float64 (3)
        Averaged terrain normal (not normalised) [-]

    Notes
    -----
    - The below code guarantees the following conditions:
      - np.all(sw_dir_cor[:, 0] == 0.0)
      - np.all(illuminated[:, 0] == 0)
      - Range of 'sw_dir_cor'-values: [0.0, sw_dir_cor_max]
    - Strictly speaking, the spatial averaging of 'sw_dir_cor', 'illuminated'
      and 'terrain_normal' is not correct because these quantities are all
      computed in slightly different reference systems (local ENU coordinates).
      However, this effect is negligible for parent grid cells with surfaces
      areas in the square kilometre range and below.
    """

    # Check validity of input arguments
    if (idx_tri_start < 0) or (idx_tri_stop > tri_face.shape[0]):
        raise ValueError("'idx_tri_start' and/or 'idx_tri_stop' ouf of bounds")
    if tri_face.shape[0] != horizon.shape[0]:
        raise ValueError("'tri_face' and 'horizon' have incompatible shapes")

    num_azim = horizon.shape[1]
    area_factor = 1.0  # currently constant
    sw_dir_cor = np.zeros((num_azim, num_elev), dtype=np.float64)
    illuminated = np.zeros((num_azim, num_elev), dtype=np.uint32)
    terrain_normal = np.zeros(3, dtype=np.float64)

    # Evaluated trigonometric functions
    azim = np.deg2rad(np.arange(0.0, 360.0, 360.0 / num_azim))  # [rad]
    azim_sin = np.sin(azim)
    azim_cos = np.cos(azim)
    elev = np.deg2rad(np.linspace(0.0, 90.0, num_elev))  # [rad]
    elev_sin = np.sin(elev)
    elev_cos = np.cos(elev)

    for idx_tri in range(idx_tri_start, idx_tri_stop):

        # Rotation matrix to transform from global to local ENU coordinates
        vert_0, vert_1, vert_2 = tri_vert[tri_face[idx_tri, :]]
        # winding order in 'tri_face': counter-clockwise -> relevant for
        # below 'tri_normal' calculation
        centroid = (vert_0 + vert_1 + vert_2) / np.float32(3.0)
        sphere_normal = centroid - earth_centre
        sphere_normal /= np.linalg.norm(sphere_normal)  # [-]
        v_n = north_pole - centroid
        dot_prod = np.dot(v_n, sphere_normal)
        north_dir = v_n - sphere_normal * dot_prod
        north_dir /= np.linalg.norm(north_dir)  # [-]
        east_dir = np.cross(north_dir, sphere_normal)  # [-]

        # Triangle face normal (global ENU coordinates)
        tri_normal = np.cross(vert_2 - vert_1, vert_0 - vert_1)
        tri_normal /= np.linalg.norm(tri_normal)  # [-]

        # Triangle face normal (local ENU coordinates)
        tri_normal_x = np.dot(east_dir, tri_normal)
        tri_normal_y = np.dot(north_dir, tri_normal)
        tri_normal_z = np.dot(sphere_normal, tri_normal)
        terrain_normal[0] += tri_normal_x
        terrain_normal[1] += tri_normal_y
        terrain_normal[2] += tri_normal_z

        # Loop through sun azimuth and elevation angles
        for idx_azim in range(num_azim):
            for idx_elev in range(1, num_elev):
                # for idx_elev = 0, set sw_dir_cor to 0.0 because
                # dot_prod_hs = dot(sun_dir, up_local) = 0.0
                dir_x = elev_cos[idx_elev] * azim_sin[idx_azim]
                dir_y = elev_cos[idx_elev] * azim_cos[idx_azim]
                dir_z = elev_sin[idx_elev]
                dot_prod_ts = (
                    dir_x * tri_normal_x +
                    dir_y * tri_normal_y +
                    dir_z * tri_normal_z
                )
                if dot_prod_ts <= 0.0: # sw_dir_cor += 0.0, illuminated += 0
                    continue
                dot_prod_hs = dir_z
                is_illum = (
                    elev[idx_elev] > np.deg2rad(horizon[idx_tri, idx_azim])
                ) # True (1): illuminated, False (0): shadow
                sw_dir_cor[idx_azim, idx_elev] += min(
                    (1.0 / dot_prod_hs) * area_factor
                    * np.float64(is_illum) * dot_prod_ts,
                    sw_dir_cor_max
                )
                illuminated[idx_azim, idx_elev] += np.uint32(is_illum)

    # Average shortwave correction factors and terrain normal
    num_tri = idx_tri_stop - idx_tri_start
    sw_dir_cor /= float(num_tri)
    terrain_normal /= float(num_tri)

    return sw_dir_cor, illuminated, terrain_normal


@njit(float64[:, :](
    uint32[:, :],
    float64[:],
    int32,
), cache=True)
def get_elev_shadow(
    illuminated,
    elev,
    num_cell,
):
    """
    Compute relevant elevation angles for fractional shadow.

    Parameters
    ----------
    illuminated : ndarray of uint32 (num_azim, num_elev)
        Number of illuminated cells, increasing monotonically with elevation
        angle
    elev : ndarray of float64 (num_elev)
        Elevation angle [deg]
    num_cell : int32
        Total number of cells

    Returns
    -------
    elev_shadow : ndarray of float64 (num_azim, 3)
        0: last elevation angle with total shadow [deg]
        1: first elevation angle with full illumination [deg]
        2: first elevation angle with the majority of cells illuminated [deg]
    """
    num_azim, num_elev = illuminated.shape
    elev_shadow = np.empty((num_azim, 3), dtype=np.float64)
    num_cell_half = int(num_cell / 2.0)
    for idx_azim in range(num_azim):

        idx_0 = -1
        idx_1 = -1
        idx_2 = -1

        for idx_elev in range(num_elev):
            if illuminated[idx_azim, idx_elev] == 0:
                idx_0 = idx_elev
            if idx_1 == -1 and illuminated[idx_azim, idx_elev] == num_cell:
                idx_1 = idx_elev
            if idx_2 == -1 and illuminated[idx_azim, idx_elev] > num_cell_half:
                idx_2 = idx_elev

        elev_shadow[idx_azim, 0] = elev[idx_0]
        elev_shadow[idx_azim, 1] = elev[idx_1]
        elev_shadow[idx_azim, 2] = elev[idx_2]

    return elev_shadow


@njit(float64[:](
    float64,
    float64,
    int32,
    float64,
), cache=True)
def spacing_exp(
    x_start,
    x_stop,
    num_nodes,
    eta,
):
    """
    Compute spacing between x_start and x_stop with increasing spacing. The
    output array starts/ends exactly with x_start/x_stop.

    Parameters
    ----------
    x_start : float64
        Start of the spacing
    x_stop : float64
        End of the spacing.
    num_nodes : int32
        Number of points in the spacing
    eta : float64
        Exponent for the spacing. Must be >= 1.0.

    Returns
    -------
    x_spacing : ndarray of float64 (num_nodes)
        Array with spacing
    """
    x_spacing = np.empty(num_nodes, dtype=np.float64)
    x_spacing[0] = x_start
    for i in range(1, num_nodes - 1):
        x_spacing[i] = x_start + (x_stop - x_start) \
            * (float(i) / float(num_nodes - 1)) ** eta
    x_spacing[num_nodes - 1] = x_stop
    return x_spacing


@njit(float64[:, :](
    float64[:, :],
    float64[:],
    float64[:, :],
    int32,
    float64,
), cache=True)
def compress_sw_dir_cor(
    sw_dir_cor,
    elev,
    elev_shadow,
    num_nodes,
    eta,
):
    """
    Compress correction information by interpolating the values to an
    exponentially increasing elevation angle spacing.

    Parameters
    ----------
    sw_dir_cor : ndarray of float64 (num_azim, num_elev)
        Direct shortwave correction factor [-]
    elev : ndarray of float64 (num_elev)
        Elevation angles [deg]
    elev_shadow : ndarray of float64 (num_azim, 2)
        0: last elevation angle with total shadow [deg]
        1: first elevation angle with full illumination [deg]
    num_nodes : int32
        Number of nodes along the elevation angle
    eta : float64
        Exponent for the spacing. Must be >= 1.0.

    Returns
    -------
    sw_dir_cor_sparse : ndarray of float64 (num_azim, num_nodes)
        Compressed direct shortwave correction factor [-]
    """
    sw_dir_cor_sparse = np.empty(
        (sw_dir_cor.shape[0], num_nodes), dtype=np.float64
    )
    for idx_azim in range(sw_dir_cor.shape[0]):
        elev_start = elev_shadow[idx_azim, 0]
        elev_stop = elev_shadow[idx_azim, 1]
        elev_sparse = spacing_exp(elev_start, elev_stop, num_nodes, eta)
        sw_dir_cor_sparse[idx_azim, :] = np.interp(
            x=elev_sparse, xp=elev, fp=sw_dir_cor[idx_azim, :]
        )
    return sw_dir_cor_sparse


@measure_time
@njit(void(
    int32,             # num_cell_child_per_parent
    int32,             # num_nodes
    float64,           # sw_dir_cor_max
    float64,           # sw_dir_cor_agg_max
    float64,           # eta_sel
    int64[:],          # slice_loc_parent
    float64[:],        # elev
    float32[:, :],     # vertices_child
    uint32[:, :],      # faces_child_sel
    float32[:],        # earth_centre
    float32[:],        # north_pole
    float32[:, :],     # horizon
    float32[:, :],     # terrain_normal_all
    float32[:, :, :],  # elev_shadow_all
    float32[:, :, :],  # sw_dir_cor_sparse_all
), parallel=True, cache=True)
def process_block(
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
):
    """
    Parallel block processing of parent grid cells. First, spatially aggregated
    values over a child cell are computed and relevant elevation angles for
    fractional shadow are derived. Subsequently, the correction information is
    compressed.

    Parameters
    ----------
    num_cell_child_per_parent : int32
        Number of child cells per parent cell
    num_nodes : int32
        Number of nodes along the elevation angle
    sw_dir_cor_max : float64
        Maximum correction for individual values [-]
    sw_dir_cor_agg_max : float64
        Maximum correction for spatially aggregated values [-]
    eta_sel : float64
        Exponent for the spacing. Must be >= 1.0
    slice_loc_parent : ndarray of int64 (2)
        Slice indices for parent grid
    elev : ndarray of float64 (num_elev)
        Elevation angles [deg]
    vertices_child : ndarray of float32 (num_cell, 3)
        Vertices of child triangle mesh [m]
    faces_child_sel : ndarray of uint32 (num_face, 3)
        Faces of child triangle mesh (same sub-selection as for horizon)
    earth_centre : ndarray of float32 (3)
        Earth centre in cartesian coordinates [m]
    north_pole : ndarray of float32 (3)
        North Pole in cartesian coordinates [m]
    horizon : ndarray of float32 (num_face, num_azim)
        Terrain horizon of child triangles [deg]
    terrain_normal_all : ndarray of float32 (num_parent, 3) (output)
        Terrain normal [-]
    elev_shadow_all : ndarray of float32 (num_parent, num_azim, 3) (output)
        Relevant shadow angles (no, full and half illumination) [deg]
    sw_dir_cor_sparse_all : ndarray of float32
                            (num_parent, num_azim, num_nodes) (output)
        Compressed direct shortwave correction factor [-]
    """
    num_elev = elev.size
    num_parent = slice_loc_parent[1] - slice_loc_parent[0]
    invalid_illuminated = np.zeros(num_parent, dtype=np.bool_)

    # Loop through parent cells (ICON grid cells)
    for idx_parent_rel in prange(num_parent):

        idx_loc_start = num_cell_child_per_parent * idx_parent_rel
        idx_loc_stop = num_cell_child_per_parent * (idx_parent_rel + 1)

        # Aggregate values over child cells
        sw_dir_cor, illuminated, terrain_normal = compute_sw_dir_cor_agg(
            vertices_child,
            faces_child_sel,
            earth_centre,
            north_pole,
            horizon,
            idx_loc_start,
            idx_loc_stop,
            num_elev,
            sw_dir_cor_max,
        )

        idx_parent = slice_loc_parent[0] + idx_parent_rel
        terrain_normal_all[idx_parent, :] = terrain_normal

        # Get relevant elevation angles for fractional shadow
        for i in range(illuminated.shape[0]):
            if illuminated[i, -1] != num_cell_child_per_parent:
                invalid_illuminated[idx_parent_rel] = True
            for j in range(1, illuminated.shape[1]):
                if illuminated[i, j] < illuminated[i, j - 1]:
                    invalid_illuminated[idx_parent_rel] = True
        elev_shadow = get_elev_shadow(
            illuminated,
            elev,
            num_cell_child_per_parent,
        )
        # 0: last elevation angle with total shadow
        # 1: first elevation angle with full illumination
        # 2: first elevation angle with the majority of cells illuminated

        elev_shadow_all[idx_parent, :, :] = elev_shadow

        # Set upper limit for correction values
        sw_dir_cor = np.minimum(sw_dir_cor, sw_dir_cor_agg_max)

        # Compress correction values
        sw_dir_cor_sparse = compress_sw_dir_cor(
            sw_dir_cor,  # C-contiguous
            elev,  # C-contiguous
            elev_shadow[:, :2],  # not C-contiguous!
            num_nodes,
            eta_sel
        )
        sw_dir_cor_sparse_all[idx_parent, :, :] = sw_dir_cor_sparse

    if np.any(invalid_illuminated):
        raise ValueError("Array 'illuminated' is invalid")
