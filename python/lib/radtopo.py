import functools
import time
import logging

import numpy as np
from numba import njit, prange
from numba import boolean, float32, float64, int32, int64, uint32, void
from numba import types  # type: ignore
'''
Module utilities with auxiliary functions to compute grid- and subgrid-scale
radiation-topography correction parameters.
'''

RADIUS_EARTH = 6_371_229.0  # ICON/COSMO earth radius [m]

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
        logging.info(f"{func.__name__}: {time_end - time_start:.2f} s")
        return result

    return wrapper


@njit(void(float64[:, :]), parallel=True, cache=True)
def lonlat2ecef(coord):
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
    if coord.shape[1] != 3:
        raise ValueError("'coord' must have shape (num_points, 3)")

    for i in prange(coord.shape[0]):
        lon = coord[i, 0]
        lat = coord[i, 1]
        radius = RADIUS_EARTH + coord[i, 2]
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
    if coord.shape[1] != 3:
        raise ValueError("'coord' must have shape (num_points, 3)")

    sin_lon_origin = np.sin(lon_origin)
    cos_lon_origin = np.cos(lon_origin)
    sin_lat_origin = np.sin(lat_origin)
    cos_lat_origin = np.cos(lat_origin)

    x_origin = RADIUS_EARTH * cos_lat_origin * cos_lon_origin
    y_origin = RADIUS_EARTH * cos_lat_origin * sin_lon_origin
    z_origin = RADIUS_EARTH * sin_lat_origin

    for i in prange(coord.shape[0]):
        dx = coord[i, 0] - x_origin
        dy = coord[i, 1] - y_origin
        dz = coord[i, 2] - z_origin
        # component of offset perpendicular to the earth's axis (in meridian
        # plane of origin)
        t = cos_lon_origin * dx + sin_lon_origin * dy
        coord[i, 0] = -sin_lon_origin * dx + cos_lon_origin * dy
        coord[i, 1] = -sin_lat_origin * t + cos_lat_origin * dz
        coord[i, 2] = cos_lat_origin * t + sin_lat_origin * dz


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
            s += 1.0 - x**(scaling + 1)
        svf[i] = s / num_azim
    return svf


# -----------------------------------------------------------------------------
# Grid-scale
# -----------------------------------------------------------------------------


@measure_time
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
    int32[:, :],
),
      cache=True)
def build_tri_mesh_circ_vert(
    lon_circ,
    lat_circ,
    elevation_circ,
    lon_vert,
    lat_vert,
    cells_of_vertex,
    neighbor_cell_index,
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
        Indices of ICON grid cells (up to six) adjoining each grid vertex,
        padded with -1 for vertices with fewer than six adjoining cells
    neighbor_cell_index : ndarray of int32 (3, num_cell)
        Indices of ICON grid cells (up to three) adjoining each grid cell,
        padded with -1 for cells with fewer than three adjoining cells

    Returns
    -------
    tri_vert : ndarray of float64 (num_vert_out, 3)
        Longitude, latitude and elevation of the Embree triangle mesh vertices
        [rad, rad, m]
    tri_face : ndarray of uint32 (num_face, 3)
        Indices of the Embree triangle mesh faces. The winding order of the
        vertices is counter-clockwise
    """

    num_cell = lon_circ.size  # number of ICON grid cells
    num_vert = lon_vert.size  # number of ICON grid cell vertices

    # Check validity of input arguments
    if (lat_circ.size != num_cell) or (elevation_circ.size != num_cell):
        raise ValueError("Inconsistent sizes of cell circumcenter arrays")
    if lat_vert.size != num_vert:
        raise ValueError("Inconsistent sizes of vertex arrays")
    if ((cells_of_vertex.shape[0] != 6)
            or (cells_of_vertex.shape[1] != num_vert)):
        raise ValueError("'cells_of_vertex' must have shape (6, num_vert)")
    if ((neighbor_cell_index.shape[0] != 3)
            or (neighbor_cell_index.shape[1] != num_cell)):
        raise ValueError("'neighbor_cell_index' must have shape (3, num_cell)")
    if (cells_of_vertex.min() < -1) or (cells_of_vertex.max() >= num_cell):
        raise ValueError("'cells_of_vertex' contains invalid indices")
    if ((neighbor_cell_index.min() < -1)
            or (neighbor_cell_index.max() >= num_cell)):
        raise ValueError("'neighbor_cell_index' contains invalid indices")
    for i in range(num_cell):
        if ((neighbor_cell_index[0, i] == -1)
                and (neighbor_cell_index[1, i] == -1)
                and (neighbor_cell_index[2, i] == -1)):
            raise ValueError("ICON grid cell(s) without any neighbour cell")

    # Compute cartesian coordinates of ICON cell circumcenters and vertices
    c_pts = np.empty((num_cell, 3), dtype=np.float64)
    for i in range(num_cell):
        c_pts[i, 0] = np.cos(lat_circ[i]) * np.cos(lon_circ[i])
        c_pts[i, 1] = np.cos(lat_circ[i]) * np.sin(lon_circ[i])
        c_pts[i, 2] = np.sin(lat_circ[i])
    v_pts = np.empty((num_vert, 3), dtype=np.float64)
    for i in range(num_vert):
        v_pts[i, 0] = np.cos(lat_vert[i]) * np.cos(lon_vert[i])
        v_pts[i, 1] = np.cos(lat_vert[i]) * np.sin(lon_vert[i])
        v_pts[i, 2] = np.sin(lat_vert[i])

    # Allocate output arrays
    tri_vert = np.empty((num_cell + num_vert, 3), dtype=np.float64)
    for i in range(num_cell):
        tri_vert[i, 0] = lon_circ[i]
        tri_vert[i, 1] = lat_circ[i]
        tri_vert[i, 2] = elevation_circ[i]
    tri_face = np.empty((num_vert * 6, 3), dtype=np.uint32)  # Embree buffer
    # (allocate maximal possible size; 6 triangles per ICON grid vertex)

    # Loop through ICON grid vertices
    indices_cell = np.empty(6, dtype=np.int32)
    idx_face = 0
    for idx_vert in range(num_vert):

        num_cell_adj = 0
        elevation_mean = 0.0
        for i in range(6):
            idx_cell = cells_of_vertex[i, idx_vert]
            if idx_cell != -1:
                indices_cell[num_cell_adj] = idx_cell
                num_cell_adj += 1
                elevation_mean += elevation_circ[idx_cell]
        if num_cell_adj == 0:
            raise ValueError("Grid vertex has no adjacent cells")

        idx_tri_vert = idx_vert + num_cell
        tri_vert[idx_tri_vert, 0] = lon_vert[idx_vert]
        tri_vert[idx_tri_vert, 1] = lat_vert[idx_vert]
        tri_vert[idx_tri_vert, 2] = elevation_mean / float(num_cell_adj)

        pt_0 = v_pts[idx_vert, :]
        face_idx_start = idx_face

        # Loop through cells adjacent to vertex
        for idx_cell in indices_cell[:num_cell_adj]:
            pt_1 = c_pts[idx_cell, :]
            indices_cell_neigh = neighbor_cell_index[:, idx_cell]

            # Loop through neighbour cells
            for idx_cell_neigh in indices_cell_neigh:
                if idx_cell_neigh in indices_cell[:num_cell_adj]:
                    pt_2 = c_pts[idx_cell_neigh, :]

                    # d = np.dot(np.cross(pt_1 - pt_0, pt_2 - pt_0), pt_0)
                    a_x = pt_1[0] - pt_0[0]
                    a_y = pt_1[1] - pt_0[1]
                    a_z = pt_1[2] - pt_0[2]
                    b_x = pt_2[0] - pt_0[0]
                    b_y = pt_2[1] - pt_0[1]
                    b_z = pt_2[2] - pt_0[2]
                    d = ((a_y * b_z - a_z * b_y) * pt_0[0] +
                         (a_z * b_x - a_x * b_z) * pt_0[1] +
                         (a_x * b_y - a_y * b_x) * pt_0[2])
                    if d > 0.0:  # counter-clockwise
                        idx_1, idx_2 = idx_cell, idx_cell_neigh
                    elif d < 0.0:  # clockwise
                        idx_1, idx_2 = idx_cell_neigh, idx_cell
                    else:
                        raise ValueError("Degenerate triangle")
                    # Add triangle if not already present for this vertex
                    if idx_1 not in tri_face[face_idx_start:idx_face, 1]:
                        tri_face[idx_face, 0] = idx_tri_vert
                        tri_face[idx_face, 1] = idx_1
                        tri_face[idx_face, 2] = idx_2
                        idx_face += 1

    tri_face = tri_face[:idx_face, :]

    # Check that all ICON grid cell circumcenters are part of the mesh
    circ_used = np.zeros(num_cell, dtype=np.bool_)
    for i in range(idx_face):
        circ_used[tri_face[i, 1]] = True
        circ_used[tri_face[i, 2]] = True
    if not np.all(circ_used):
        raise ValueError("ICON grid cell circumcenter(s) not part of the "
                         "output triangle mesh")

    return tri_vert, tri_face


# -----------------------------------------------------------------------------
# Subgrid-scale
# -----------------------------------------------------------------------------


@njit(int32(int32[:], int32), cache=True)
def _local_corner(idx_corners, idx_vertex):
    """Return local corner number (0, 1, 2) of a global vertex index."""
    for i in range(3):
        if idx_corners[i] == idx_vertex:
            return i
    raise ValueError("Edge vertex is not a vertex of the cell")


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
),
      cache=True)
def refine_tri_mesh(
    vertices,
    vertex_of_cell,
    edge_of_cell,
    edge_vertices,
    n,
):
    """
    Refine a triangle mesh on a sphere surface by subdividing each parent
    triangle into n ** 2 child triangles.

    Parameters
    ----------
    vertices : ndarray of float64 (num_vert, 3)
        Cartesian coordinates (x, y, z) of vertices on unit sphere [-]
    vertex_of_cell : ndarray of int32 (3, num_cell)
        Indices of triangle vertices
    edge_of_cell : ndarray of int32 (3, num_cell)
        Indices of triangle edges
    edge_vertices : ndarray of int32 (2, num_edge)
        Indices of vertices of edges
    n : int32
        Triangle division number

    Returns
    -------
    vertices_child : ndarray of float64 (num_vert_child, 3)
        Cartesian coordinates (x, y, z) of refined vertices on unit sphere [-]
    faces_child : ndarray of uint32 (num_cell * n ** 2, 3)
        Indices of refined triangle faces. The winding order is the same as in
        the input (vertex_of_cell)

    Notes
    -----
    - child vertices are stored in the order: (1) parent vertices, (2) vertices
      on parent edges (n - 1 per edge, shared by adjacent cells), (3) vertices
      in the interior of parent cells ((n - 1) * (n - 2) / 2 per cell).
    - Within a parent cell, a child vertex is addressed by the lattice position
      (i, j) with i + j <= n, which corresponds to the point
      ((n - i - j) * v_0 + i * v_1 + j * v_2) / n (before projection onto the
      unit sphere).
    - The child triangles of parent cell 'i' are stored contiguously in rows
      [i * n ** 2, (i + 1) * n ** 2).
    """

    # Number of parent vertices, edges and cells
    num_vert_p = vertices.shape[0]
    num_edge_p = edge_vertices.shape[1]
    num_cell_p = vertex_of_cell.shape[1]

    # Check validity of input arguments
    if n < 1:
        raise ValueError("Triangle division number 'n' must be at least 1")
    if vertices.shape[1] != 3:
        raise ValueError("'vertices' must have shape (num_vert, 3)")
    if vertex_of_cell.shape[0] != 3:
        raise ValueError("'vertex_of_cell' must have shape (3, num_cell)")
    if (edge_of_cell.shape[0] != 3) or (edge_of_cell.shape[1] != num_cell_p):
        raise ValueError("'edge_of_cell' must have shape (3, num_cell)")
    if edge_vertices.shape[0] != 2:
        raise ValueError("'edge_vertices' must have shape (2, num_edge)")
    if (vertex_of_cell.min() < 0) or (vertex_of_cell.max() >= num_vert_p):
        raise ValueError("'vertex_of_cell' contains invalid indices")
    if (edge_of_cell.min() < 0) or (edge_of_cell.max() >= num_edge_p):
        raise ValueError("'edge_of_cell' contains invalid indices")
    if (edge_vertices.min() < 0) or (edge_vertices.max() >= num_vert_p):
        raise ValueError("'edge_vertices' contains invalid indices")

    # Compute number of child vertices and cells
    num_vert_c_edge = num_edge_p * (n - 1)
    num_vert_c_interior = (n - 1) * (n - 2) // 2 * num_cell_p
    num_vert_c = num_vert_p + num_vert_c_edge + num_vert_c_interior
    num_cell_c = num_cell_p * n**2

    # Allocate output arrays
    vertices_child = np.empty((num_vert_c, 3), dtype=np.float64)
    faces_child = np.empty((num_cell_c, 3), dtype=np.uint32)

    # Parent vertices (keep original indices for connectivity)
    vertices_child[:num_vert_p, :] = vertices

    # Vertices on parent edges (shared between adjacent cells)
    # Interpolate n - 1 points along each straight edge between its endpoints.
    # Compute them once per global edge so adjacent cells share these vertices.
    # Project these points onto the sphere after constructing the whole mesh.
    idx_vert = num_vert_p
    for idx_edge in range(num_edge_p):
        v_0 = vertices[edge_vertices[0, idx_edge], :]
        v_1 = vertices[edge_vertices[1, idx_edge], :]
        for s in range(1, n):
            vertices_child[idx_vert, :] = ((n - s) * v_0 + s * v_1) / n
            idx_vert += 1

    corner_i = np.array([0, n, 0], dtype=np.int32)
    corner_j = np.array([0, 0, n], dtype=np.int32)
    # Use a triangular lattice with v_0 at (0, 0), v_1 at (n, 0),
    # and v_2 at (0, n). Before projection, each step in i is
    # (v_1 - v_0) / n; each step in j is (v_2 - v_0) / n.
    lattice = np.empty((n + 1, n + 1), dtype=np.uint32)  # child vertex indices
    # Map lattice positions (i, j) to global indices in vertices_child.
    # Refill the corners, edges and interior for each parent in this array.
    # Use only positions with i + j <= n; the rest lie outside the triangle.
    idx_face = 0
    for idx_cell in range(num_cell_p):
        idx_corners = vertex_of_cell[:, idx_cell]
        v_0 = vertices[idx_corners[0], :]
        v_1 = vertices[idx_corners[1], :]
        v_2 = vertices[idx_corners[2], :]

        # Corners
        for c in range(3):
            lattice[corner_i[c], corner_j[c]] = idx_corners[c]

        # Edges (direction of edge determined from 'edge_vertices')
        # Locate each edge's endpoints among this cell's corners, then walk
        # between them in lattice coordinates. Follow the global endpoint order
        # to match the precomputed vertex indices, regardless of cell winding.
        for idx_edge in edge_of_cell[:, idx_cell]:
            c_0 = _local_corner(idx_corners, edge_vertices[0, idx_edge])
            c_1 = _local_corner(idx_corners, edge_vertices[1, idx_edge])
            # One subdivision step along the edge in lattice coordinates.
            di = (corner_i[c_1] - corner_i[c_0]) // n
            dj = (corner_j[c_1] - corner_j[c_0]) // n
            # Each edge owns n - 1 entries after the parent vertex block.
            idx_vert_edge = num_vert_p + idx_edge * (n - 1)
            for s in range(1, n):
                lattice[corner_i[c_0] + s * di, corner_j[c_0] + s * dj] \
                    = idx_vert_edge + s - 1

        # Interior (vertices not shared -> computed here)
        # Compute interior coordinates as a weighted average of the corners.
        # Append new indices: these vertices belong only to this parent.
        for i in range(1, n - 1):
            for j in range(1, n - i):
                vertices_child[idx_vert, :] \
                    = ((n - i - j) * v_0 + i * v_1 + j * v_2) / n
                lattice[i, j] = idx_vert
                idx_vert += 1

        # Child faces ('upward' and 'downward' triangles)
        for i in range(n):
            for j in range(n - i):
                # Four neighboring lattice points can form two triangles.
                # Use three for the first; both triangles share the edge
                # between (i + 1, j) and (i, j + 1), with parent winding.
                faces_child[idx_face, 0] = lattice[i, j]
                faces_child[idx_face, 1] = lattice[i + 1, j]
                faces_child[idx_face, 2] = lattice[i, j + 1]
                idx_face += 1
                # Add the second only if (i + 1, j + 1) is inside.
                # On the boundary i + j = n - 1, it lies outside.
                if i + j + 1 < n:
                    faces_child[idx_face, 0] = lattice[i + 1, j]
                    faces_child[idx_face, 1] = lattice[i + 1, j + 1]
                    faces_child[idx_face, 2] = lattice[i, j + 1]
                    idx_face += 1

    # Project vertices onto unit sphere
    # Interpolation placed new points on flat edges and parent triangle planes.
    # Normalize their vectors to project them radially onto the unit sphere;
    # also normalize the copied parent vertices.
    for i in range(num_vert_c):
        vertices_child[i, :] /= np.linalg.norm(vertices_child[i, :])

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
    idx_lin_tile : ndarray of int32 (num_points)
        Linear index of tile associated with point
    """

    # Check validity of input arguments
    if lon.size != lat.size:
        raise ValueError("'lon' and 'lat' must have equal size")
    if (lon.min() < -np.pi) or (lon.max() > +np.pi):
        raise ValueError("Longitude value(s) out of range [-pi, +pi]")
    if (lat.min() < -(np.pi / 2.0)) or (lat.max() > +(np.pi / 2.0)):
        raise ValueError("Latitude value(s) out of range [-pi/2, +pi/2]")
    if not isinstance(num_tile_lon, (int, np.integer)):
        raise TypeError("Number of longitudinal tiles must be an integer")
    if not isinstance(num_tile_lat, (int, np.integer)):
        raise TypeError("Number of latitudinal tiles must be an integer")
    if (num_tile_lon <= 0) or (360 % num_tile_lon != 0):
        raise ValueError("Invalid number of longitudinal tiles")
    if (num_tile_lat <= 0) or (180 % num_tile_lat != 0):
        raise ValueError("Invalid number of latitudinal tiles")

    tile_extent_lon = 2.0 * np.pi / num_tile_lon
    tile_extent_lat = np.pi / num_tile_lat
    idx_tile_lon = np.floor((lon + np.pi) / tile_extent_lon).astype(np.int16)
    idx_tile_lat = np.floor(
        ((np.pi / 2.0) - lat) / tile_extent_lat).astype(np.int16)
    idx_tile_lon = np.clip(idx_tile_lon, 0, num_tile_lon - 1)
    idx_tile_lat = np.clip(idx_tile_lat, 0, num_tile_lat - 1)

    return idx_tile_lat.astype(np.int32) * num_tile_lon + idx_tile_lon


def get_tile_name(
    idx_tile_lon,
    idx_tile_lat,
    num_tile_lon,
    num_tile_lat,
):
    """
    Return coordinates part of tile name based on longitudinal and latitudinal
    tile indices. Tile indices increase eastward from -180 deg longitude and
    southward from +90 deg latitude. Wrap index around longitude.

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
    lon_add : float
        Additive term for tile longitude to wrap [deg]
    """

    # Check validity of input arguments
    args = (idx_tile_lon, idx_tile_lat, num_tile_lon, num_tile_lat)
    if not all(isinstance(arg, (int, np.integer)) for arg in args):
        raise TypeError("All input arguments must be integers")
    if (num_tile_lon <= 0) or (360 % num_tile_lon != 0):
        raise ValueError("Invalid number of longitudinal tiles")
    if (num_tile_lat <= 0) or (180 % num_tile_lat != 0):
        raise ValueError("Invalid number of latitudinal tiles")
    if (idx_tile_lat < 0) or (idx_tile_lat >= num_tile_lat):
        raise ValueError("Latitudinal tile index out of bounds")

    tile_extent_lon = 360 // num_tile_lon
    lon_add = idx_tile_lon // num_tile_lon * 360.0
    idx_tile_lon %= num_tile_lon  # wrap index around longitude
    lon_west = -180 + idx_tile_lon * tile_extent_lon
    lon_east = lon_west + tile_extent_lon
    letter_west = "E" if lon_west >= 0 else "W"
    letter_east = "E" if lon_east >= 0 else "W"

    tile_extent_lat = 180 // num_tile_lat
    lat_north = 90 - idx_tile_lat * tile_extent_lat
    lat_south = lat_north - tile_extent_lat
    letter_north = "N" if lat_north >= 0 else "S"
    letter_south = "N" if lat_south >= 0 else "S"

    tile_name_coord = (f"{letter_north}{abs(lat_north):02d}-"
                       f"{letter_south}{abs(lat_south):02d}_"
                       f"{letter_west}{abs(lon_west):03d}-"
                       f"{letter_east}{abs(lon_east):03d}")

    return tile_name_coord, lon_add


@measure_time
@njit(float32[:](
    float32[:, :],
    float64[::1],
    float64[::1],
    float64[::1],
    float64[::1],
    float64,
),
      parallel=True,
      cache=True)
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
    if not np.isfinite(atol_grid) or atol_grid < 0.0:
        raise ValueError("'atol_grid' must be finite and non-negative")
    if not np.all(np.isfinite(x_axis)) or not np.all(np.isfinite(y_axis)):
        raise ValueError("Coordinate axes must contain only finite values")

    dx = np.diff(x_axis)
    dy = np.diff(y_axis)
    if not (np.all(dx > 0.0) or np.all(dx < 0.0)):
        raise ValueError("x-axis must be strictly ascending or descending")
    if not (np.all(dy > 0.0) or np.all(dy < 0.0)):
        raise ValueError("y-axis must be strictly ascending or descending")

    delta_x = (x_axis[-1] - x_axis[0]) / (len_x - 1)
    delta_y = (y_axis[-1] - y_axis[0]) / (len_y - 1)
    if not np.isfinite(delta_x) or delta_x == 0.0:
        raise ValueError("x-axis spacing must be finite and nonzero")
    if not np.isfinite(delta_y) or delta_y == 0.0:
        raise ValueError("y-axis spacing must be finite and nonzero")
    if np.max(np.abs(dx - delta_x)) > atol_grid:
        raise ValueError("x-axis must be regularly spaced")
    if np.max(np.abs(dy - delta_y)) > atol_grid:
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
        data_interp[k] = ((1.0 - weight_x) *
                          (1.0 - weight_y) * data[j_0, i_0] + weight_x *
                          (1.0 - weight_y) * data[j_0, i_1] +
                          (1.0 - weight_x) * weight_y * data[j_1, i_0] +
                          weight_x * weight_y * data[j_1, i_1])

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
    float64[:],
),
      cache=True)
def _compute_sw_dir_cor_agg(
    tri_vert,
    tri_face,
    earth_centre,
    north_pole,
    horizon,
    elev,
):
    """
    Compute direct shortwave correction factors and fractional illumination
    for a hemisphere, aggregated over all passed triangles. Additionally, the
    averaged terrain normal is returned.

    Parameters
    ----------
    tri_vert : ndarray of float32 (num_vert, 3)
        Cartesian coordinates of triangle mesh vertices [m]
    tri_face : ndarray of uint32 (num_face, 3)
        Indices of the considered triangle mesh faces (counter-clockwise
        winding order)
    earth_centre : ndarray of float32 (3)
        Earth centre in cartesian coordinates [m]
    north_pole : ndarray of float32 (3)
        North Pole in cartesian coordinates [m]
    horizon : ndarray of float32 (num_face, num_azim)
        Terrain horizon of the considered triangles [deg]
    elev : ndarray of float64 (num_elev)
        Elevation angles, strictly ascending from 0.0 to 90.0 [deg]

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
      - Range of 'sw_dir_cor'-values:
        [0.0, 1.0 / np.sin(np.deg2rad(elev[1]))]
    - Strictly speaking, the spatial averaging of 'sw_dir_cor', 'illuminated'
      and 'terrain_normal' is not correct because these quantities are all
      computed in slightly different reference systems (local ENU coordinates).
      However, this effect is negligible for parent grid cells with surfaces
      areas in the square kilometre range and below.
    """

    num_azim = horizon.shape[1]
    num_elev = elev.size
    area_factor = 1.0  # currently constant
    sw_dir_cor = np.zeros((num_azim, num_elev), dtype=np.float64)
    illuminated = np.zeros((num_azim, num_elev), dtype=np.uint32)
    terrain_normal = np.zeros(3, dtype=np.float64)

    # Evaluated trigonometric functions
    azim = np.deg2rad(np.arange(0.0, 360.0, 360.0 / num_azim))  # [rad]
    azim_sin = np.sin(azim)
    azim_cos = np.cos(azim)
    elev_rad = np.deg2rad(elev)  # [rad]
    elev_sin = np.sin(elev_rad)
    elev_cos = np.cos(elev_rad)

    num_tri = tri_face.shape[0]
    for idx_tri in range(num_tri):

        # Rotation matrix to transform from global to local ENU coordinates
        vert_0, vert_1, vert_2 = tri_vert[tri_face[idx_tri, :]]
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
        # points 'upward' for counter-clockwise oriented triangles
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
            horizon_rad = np.deg2rad(np.float64(horizon[idx_tri, idx_azim]))
            # for idx_elev = 0, set sw_dir_cor to 0.0 because
            # dot_prod_hs = dir_z = dot(sun_dir, up_local) = 0.0
            for idx_elev in range(1, num_elev):
                if elev_rad[idx_elev] <= horizon_rad:
                    # sw_dir_cor += 0.0, illuminated += 0
                    continue
                dir_x = elev_cos[idx_elev] * azim_sin[idx_azim]
                dir_y = elev_cos[idx_elev] * azim_cos[idx_azim]
                dir_z = elev_sin[idx_elev]
                dot_prod_ts = (dir_x * tri_normal_x + dir_y * tri_normal_y +
                               dir_z * tri_normal_z)
                if dot_prod_ts <= 0.0:  # sw_dir_cor += 0.0, illuminated += 0
                    continue
                sw_dir_cor[idx_azim, idx_elev] += \
                    ((1.0 / dir_z) * area_factor * dot_prod_ts)
                # range of '(1.0 / dir_z)':
                # [1.0, 1.0 / np.sin(elev_rad[1])]
                # range of 'dot_prod_ts': [>0.0, 1.0]
                # -> range of above additive term is restricted to:
                # [0.0, 1.0 / np.sin(elev_rad[1])]
                # (ca. 114.6 for elev[1] = 0.5 deg)
                illuminated[idx_azim, idx_elev] += np.uint32(1)

    # Average shortwave correction factors and terrain normal
    sw_dir_cor /= float(num_tri)
    terrain_normal /= float(num_tri)

    return sw_dir_cor, illuminated, terrain_normal


@njit(types.Tuple((
    float64[:, :],
    boolean,
))(
    uint32[:, :],
    float64[:],
    int32,
),
      cache=True)
def _get_elev_shadow(
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
        angle from 0 (horizon) to 'num_cell' (zenith)
    elev : ndarray of float64 (num_elev)
        Elevation angles, strictly ascending from 0.0 to 90.0 [deg]
    num_cell : int32
        Total number of cells

    Returns
    -------
    elev_shadow : ndarray of float64 (num_azim, 3)
        0: last elevation angle with total shadow [deg]
        1: first elevation angle with full illumination [deg]
        2: first elevation angle with the majority of cells illuminated [deg]
    valid : bool
        False if 'illuminated' violates the above conditions (in this case,
        'elev_shadow' is meaningless)
    """
    num_azim, num_elev = illuminated.shape
    elev_shadow = np.empty((num_azim, 3), dtype=np.float64)
    num_cell_half = num_cell // 2
    valid = True
    for idx_azim in range(num_azim):

        # No illumination at the horizon, full illumination at zenith
        if (illuminated[idx_azim, 0] != 0) or \
                (illuminated[idx_azim, -1] != num_cell):
            valid = False

        idx_0 = -1
        idx_1 = -1
        idx_2 = -1

        for idx_elev in range(num_elev):
            # Number of illuminated cells must increase monotonically
            if (idx_elev > 0) and (illuminated[idx_azim, idx_elev]
                                   < illuminated[idx_azim, idx_elev - 1]):
                valid = False
            if illuminated[idx_azim, idx_elev] == 0:
                idx_0 = idx_elev
            if idx_1 == -1 and illuminated[idx_azim, idx_elev] == num_cell:
                idx_1 = idx_elev
            if idx_2 == -1 and illuminated[idx_azim, idx_elev] > num_cell_half:
                idx_2 = idx_elev

        elev_shadow[idx_azim, 0] = elev[idx_0]
        elev_shadow[idx_azim, 1] = elev[idx_1]
        elev_shadow[idx_azim, 2] = elev[idx_2]

    return elev_shadow, valid


@njit(float64[:](
    float64,
    float64,
    int32,
    float64,
), cache=True)
def _spacing_power(
    x_start,
    x_stop,
    num_nodes,
    eta,
):
    """
    Generate power-law-spaced nodes between x_start and x_stop.
    For an increasing interval, eta = 1 gives uniform spacing and
    eta > 1 gives increasing interval sizes (up to floating-point
    rounding). Both endpoints are assigned exactly.

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
),
      cache=True)
def _compress_sw_dir_cor(
    sw_dir_cor,
    elev,
    elev_shadow,
    num_nodes,
    eta,
):
    """
    Compress correction information by interpolating the values to a
    power-law-spaced elevation angle grid.

    Parameters
    ----------
    sw_dir_cor : ndarray of float64 (num_azim, num_elev)
        Direct shortwave correction factor [-]
    elev : ndarray of float64 (num_elev)
        Elevation angles, strictly ascending from 0.0 to 90.0 [deg]
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
    sw_dir_cor_sparse = np.empty((sw_dir_cor.shape[0], num_nodes),
                                 dtype=np.float64)
    for idx_azim in range(sw_dir_cor.shape[0]):
        elev_start = elev_shadow[idx_azim, 0]
        elev_stop = elev_shadow[idx_azim, 1]
        elev_sparse = _spacing_power(elev_start, elev_stop, num_nodes, eta)
        sw_dir_cor_sparse[idx_azim, :] = np.interp(x=elev_sparse,
                                                   xp=elev,
                                                   fp=sw_dir_cor[idx_azim, :])
    return sw_dir_cor_sparse


@measure_time
@njit(
    void(
        int32,  # num_cell_child_per_parent
        int32,  # num_nodes
        float64,  # sw_dir_cor_agg_max
        float64,  # eta_sel
        int64[:],  # slice_loc_parent
        float64[:],  # elev
        float32[:, :],  # vertices_child
        uint32[:, :],  # faces_child_sel
        float32[:],  # earth_centre
        float32[:],  # north_pole
        float32[:, :],  # horizon
        float32[:, :],  # terrain_normal_all
        float32[:, :, :],  # elev_shadow_all
        float32[:, :, :],  # sw_dir_cor_sparse_all
    ),
    parallel=True,
    cache=True)
def process_block(
    num_cell_child_per_parent,
    num_nodes,
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
    sw_dir_cor_agg_max : float64
        Maximum correction for spatially aggregated values [-]
    eta_sel : float64
        Exponent for the spacing. Must be >= 1.0
    slice_loc_parent : ndarray of int64 (2)
        Slice indices for parent grid
    elev : ndarray of float64 (num_elev)
        Elevation angles, strictly ascending from 0.0 to 90.0 [deg]
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
    num_parent = slice_loc_parent[1] - slice_loc_parent[0]
    num_azim = horizon.shape[1]

    # Check validity of input arguments (no errors can be raised in the
    # parallel region below)
    if num_cell_child_per_parent < 1:
        raise ValueError("'num_cell_child_per_parent' must be positive")
    if num_nodes < 2:
        raise ValueError("'num_nodes' must be at least 2")
    if (elev.size < 2) or (elev[0] != 0.0) or (elev[-1] != 90.0) or \
            np.any(elev[1:] <= elev[:-1]):
        raise ValueError("'elev' must be strictly ascending from 0.0 to "
                         "90.0 deg")
    num_child = num_parent * num_cell_child_per_parent
    if (faces_child_sel.shape[0] != num_child) or \
            (faces_child_sel.shape[1] != 3) or \
            (horizon.shape[0] != num_child):
        raise ValueError("'faces_child_sel' and/or 'horizon' have "
                         "incompatible shapes")
    if (slice_loc_parent[0] < 0) or \
            (slice_loc_parent[1] > terrain_normal_all.shape[0]) or \
            (slice_loc_parent[1] > elev_shadow_all.shape[0]) or \
            (slice_loc_parent[1] > sw_dir_cor_sparse_all.shape[0]):
        raise ValueError("'slice_loc_parent' out of bounds")
    if (terrain_normal_all.shape[1] != 3) or \
            (elev_shadow_all.shape[1] != num_azim) or \
            (elev_shadow_all.shape[2] != 3) or \
            (sw_dir_cor_sparse_all.shape[1] != num_azim) or \
            (sw_dir_cor_sparse_all.shape[2] != np.int64(num_nodes)):
        raise ValueError("Output arrays have incompatible shapes")

    invalid_illuminated = np.zeros(num_parent, dtype=np.bool_)

    # Loop through parent cells (ICON grid cells)
    for idx_parent_rel in prange(num_parent):

        idx_loc_start = num_cell_child_per_parent * idx_parent_rel
        idx_loc_stop = num_cell_child_per_parent * (idx_parent_rel + 1)

        # Aggregate values over child cells
        sw_dir_cor, illuminated, terrain_normal = _compute_sw_dir_cor_agg(
            vertices_child,
            faces_child_sel[idx_loc_start:idx_loc_stop, :],
            earth_centre,
            north_pole,
            horizon[idx_loc_start:idx_loc_stop, :],
            elev,
        )

        idx_parent = slice_loc_parent[0] + idx_parent_rel
        terrain_normal_all[idx_parent, :] = terrain_normal

        # Get relevant elevation angles for fractional shadow
        elev_shadow, valid = _get_elev_shadow(
            illuminated,
            elev,
            num_cell_child_per_parent,
        )
        invalid_illuminated[idx_parent_rel] = not valid
        # 0: last elevation angle with total shadow
        # 1: first elevation angle with full illumination
        # 2: first elevation angle with the majority of cells illuminated

        elev_shadow_all[idx_parent, :, :] = elev_shadow

        # Set upper limit for correction values
        sw_dir_cor = np.minimum(sw_dir_cor, sw_dir_cor_agg_max)

        # Compress correction values
        sw_dir_cor_sparse = _compress_sw_dir_cor(
            sw_dir_cor,  # C-contiguous
            elev,  # C-contiguous
            elev_shadow[:, :2],  # not C-contiguous
            num_nodes,
            eta_sel)
        sw_dir_cor_sparse_all[idx_parent, :, :] = sw_dir_cor_sparse

    if np.any(invalid_illuminated):
        raise ValueError("Array 'illuminated' is invalid")
