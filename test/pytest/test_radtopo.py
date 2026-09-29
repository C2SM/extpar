import os
import sys
import itertools

import numpy as np
import pytest

try:
    from extpar.lib import radtopo
except ImportError:
    sys.path.append(
        os.path.join(os.path.dirname(__file__), '..', '..', 'python', 'lib'))
    import radtopo

RADIUS_EARTH = 6_371_229.0  # ICON/COSMO earth radius [m]

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def lonlat2cart(lon, lat):
    """Longitude/latitude [rad] to cartesian coordinates on the unit sphere."""
    return np.column_stack(
        (np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)))


def cart2lonlat(pts):
    """Cartesian coordinates to longitude/latitude [rad]."""
    pts = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    return np.arctan2(pts[:, 1], pts[:, 0]), np.arcsin(pts[:, 2])


def is_ccw(pts, faces):
    """Check if triangles are counter-clockwise when viewed from outside."""
    p0, p1, p2 = pts[faces[:, 0]], pts[faces[:, 1]], pts[faces[:, 2]]
    return np.einsum("ij,ij->i", np.cross(p1 - p0, p2 - p0), p0) > 0.0


def is_closed_manifold(faces):
    """Check if every directed edge occurs exactly once, together with its
    reverse (-> closed, consistently oriented surface)."""
    edges = {}
    for face in faces:
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2],
                                                              face[0])):
            edges[(int(a), int(b))] = edges.get((int(a), int(b)), 0) + 1
    return all(count == 1 and edges.get((b, a)) == 1
               for (a, b), count in edges.items())


@pytest.fixture
def icosahedron():
    """
    ICON-like description of an icosahedron (closed triangle mesh). Cells
    are ordered counter-clockwise and edge k of a cell joins its vertices k
    and k + 1 (like in the ICON grids). The direction of the edges in
    'edge_vertices' is randomised.
    """
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    vertices = []
    for s1, s2 in itertools.product((-1.0, 1.0), repeat=2):
        vertices += [(0.0, s1, s2 * phi), (s1, s2 * phi, 0.0),
                     (s2 * phi, 0.0, s1)]
    vertices = np.array(vertices)
    vertices /= np.linalg.norm(vertices, axis=1, keepdims=True)

    # Faces: triples of mutually adjacent vertices (shortest distance)
    dist = np.linalg.norm(vertices[:, None, :] - vertices[None, :, :], axis=2)
    dist_edge = dist[dist > 0.0].min()
    adjacent = np.isclose(dist, dist_edge)
    faces = []
    for i, j, k in itertools.combinations(range(12), 3):
        if adjacent[i, j] and adjacent[j, k] and adjacent[i, k]:
            face = [i, j, k]
            if not is_ccw(vertices, np.array([face]))[0]:
                face = [i, k, j]
            faces.append(face)
    vertex_of_cell = np.array(faces, dtype=np.int32).T  # (3, 20)

    # Edges
    rng = np.random.default_rng(42)
    edge_index = {}
    edge_vertices = []
    edge_of_cell = np.empty_like(vertex_of_cell)
    for idx_cell in range(vertex_of_cell.shape[1]):
        for k in range(3):
            a = vertex_of_cell[k, idx_cell]
            b = vertex_of_cell[(k + 1) % 3, idx_cell]
            key = (min(a, b), max(a, b))
            if key not in edge_index:
                edge_index[key] = len(edge_vertices)
                edge_vertices.append(key if rng.random() < 0.5 else key[::-1])
            edge_of_cell[k, idx_cell] = edge_index[key]
    edge_vertices = np.array(edge_vertices, dtype=np.int32).T  # (2, 30)

    # Adjacency (padded with -2, like the ICON indices after subtracting 1)
    num_vert, num_cell = vertices.shape[0], vertex_of_cell.shape[1]
    cells_of_vertex = np.full((6, num_vert), -2, dtype=np.int32)
    for idx_vert in range(num_vert):
        cells = np.where((vertex_of_cell == idx_vert).any(axis=0))[0]
        slots = rng.permutation(6)[:cells.size]  # padding at random slot
        cells_of_vertex[slots, idx_vert] = cells
    neighbor_cell_index = np.full((3, num_cell), -2, dtype=np.int32)
    for idx_cell in range(num_cell):
        for k in range(3):
            cells = np.where(
                (edge_of_cell == edge_of_cell[k, idx_cell]).any(axis=0))[0]
            neighbor_cell_index[k, idx_cell] = cells[cells != idx_cell][0]

    return {
        "vertices": vertices,
        "vertex_of_cell": vertex_of_cell,
        "edge_of_cell": edge_of_cell,
        "edge_vertices": edge_vertices,
        "cells_of_vertex": cells_of_vertex,
        "neighbor_cell_index": neighbor_cell_index,
    }


# -----------------------------------------------------------------------------
# lonlat2ecef
# -----------------------------------------------------------------------------


def test_lonlat2ecef_known_points():
    coord = np.array([
        [0.0, 0.0, 0.0],
        [np.pi / 2.0, 0.0, 0.0],
        [np.pi, 0.0, 100.0],
        [0.3, np.pi / 2.0, 0.0],
        [-1.2, -np.pi / 2.0, -50.0],
    ])
    radtopo.lonlat2ecef(coord)
    expected = np.array([
        [RADIUS_EARTH, 0.0, 0.0],
        [0.0, RADIUS_EARTH, 0.0],
        [-(RADIUS_EARTH + 100.0), 0.0, 0.0],
        [0.0, 0.0, RADIUS_EARTH],
        [0.0, 0.0, -(RADIUS_EARTH - 50.0)],
    ])
    np.testing.assert_allclose(coord, expected, rtol=0.0, atol=1e-6)


def test_lonlat2ecef_random_points():
    rng = np.random.default_rng(0)
    num = 1000
    lon = rng.uniform(-np.pi, np.pi, num)
    lat = rng.uniform(-np.pi / 2.0, np.pi / 2.0, num)
    elevation = rng.uniform(-400.0, 8_800.0, num)
    coord = np.column_stack((lon, lat, elevation))
    radtopo.lonlat2ecef(coord)
    # Radius
    np.testing.assert_allclose(np.linalg.norm(coord, axis=1),
                               RADIUS_EARTH + elevation,
                               rtol=1e-12)
    # Direction
    lon_back, lat_back = cart2lonlat(coord)
    np.testing.assert_allclose(lon_back, lon, atol=1e-10)
    np.testing.assert_allclose(lat_back, lat, atol=1e-10)


# -----------------------------------------------------------------------------
# ecef2enu
# -----------------------------------------------------------------------------


@pytest.mark.parametrize("lon_origin, lat_origin", [(0.0, 0.0), (0.15, 0.81),
                                                    (-2.5, -0.6),
                                                    (np.pi, 1.2)])
def test_ecef2enu_local_directions(lon_origin, lat_origin):
    d = 1e-5  # [rad] (ca. 64 m)
    h = 250.0  # [m]
    coord = np.array([
        [lon_origin, lat_origin, 0.0],  # origin
        [lon_origin, lat_origin, h],  # above origin
        [lon_origin + d, lat_origin, 0.0],  # east of origin
        [lon_origin, lat_origin + d, 0.0],  # north of origin
    ])
    radtopo.lonlat2ecef(coord)
    radtopo.ecef2enu(coord, lon_origin, lat_origin)

    np.testing.assert_allclose(coord[0], [0.0, 0.0, 0.0], atol=1e-6)
    np.testing.assert_allclose(coord[1], [0.0, 0.0, h], atol=1e-6)
    dist_east = RADIUS_EARTH * np.cos(lat_origin) * d
    dist_north = RADIUS_EARTH * d
    np.testing.assert_allclose(coord[2, :2], [dist_east, 0.0], atol=1e-3)
    np.testing.assert_allclose(coord[3, :2], [0.0, dist_north], atol=1e-3)
    # curvature of the earth -> points below tangent plane
    assert coord[2, 2] <= 0.0
    assert coord[3, 2] < 0.0


def test_ecef2enu_preserves_distances():
    rng = np.random.default_rng(1)
    num = 500
    coord = np.column_stack(
        (rng.uniform(-np.pi, np.pi,
                     num), rng.uniform(-np.pi / 2.0, np.pi / 2.0,
                                       num), rng.uniform(0.0, 5_000.0, num)))
    radtopo.lonlat2ecef(coord)
    coord_ecef = coord.copy()
    radtopo.ecef2enu(coord, 0.7, -0.4)
    # Rotation + translation -> pairwise distances are preserved
    idx = rng.integers(0, num, size=(200, 2))
    dist_ecef = np.linalg.norm(coord_ecef[idx[:, 0]] - coord_ecef[idx[:, 1]],
                               axis=1)
    dist_enu = np.linalg.norm(coord[idx[:, 0]] - coord[idx[:, 1]], axis=1)
    np.testing.assert_allclose(dist_enu, dist_ecef, rtol=1e-9, atol=1e-6)


def test_ecef2enu_earth_centre():
    lon_origin, lat_origin = -0.9, 0.5
    coord = np.zeros((1, 3))  # earth centre
    radtopo.ecef2enu(coord, lon_origin, lat_origin)
    np.testing.assert_allclose(coord[0], [0.0, 0.0, -RADIUS_EARTH], atol=1e-6)


# -----------------------------------------------------------------------------
# geometric_svf
# -----------------------------------------------------------------------------


@pytest.mark.parametrize("scaling", [0, 1, 2])
def test_geometric_svf_flat_and_closed(scaling):
    num_azim = 24
    horizon = np.zeros((2, num_azim), dtype=np.float32)
    horizon[1, :] = np.pi / 2.0
    svf = radtopo.geometric_svf(horizon, scaling)
    assert svf.dtype == np.float32
    np.testing.assert_allclose(svf, [1.0, 0.0], atol=1e-6)


@pytest.mark.parametrize("scaling", [0, 1, 2])
def test_geometric_svf_random_horizon(scaling):
    rng = np.random.default_rng(2)
    horizon = rng.uniform(0.0, np.deg2rad(60.0), (50, 36)).astype(np.float32)
    svf = radtopo.geometric_svf(horizon, scaling)
    expected = np.mean(1.0 - np.sin(horizon.astype(np.float64))**(scaling + 1),
                       axis=1)
    np.testing.assert_allclose(svf, expected, rtol=1e-5)
    assert np.all((svf >= 0.0) & (svf <= 1.0))


def test_geometric_svf_scaling_order():
    # Higher scaling exponent -> larger sky view factor (0 < sin(h) < 1)
    horizon = np.full((1, 16), np.deg2rad(30.0), dtype=np.float32)
    svf = [radtopo.geometric_svf(horizon, s)[0] for s in (0, 1, 2)]
    np.testing.assert_allclose(svf, [0.5, 0.75, 0.875], rtol=1e-6)


@pytest.mark.parametrize("scaling", [-1, 3])
def test_geometric_svf_invalid_scaling(scaling):
    horizon = np.zeros((2, 4), dtype=np.float32)
    with pytest.raises(ValueError):
        radtopo.geometric_svf(horizon, scaling)


# -----------------------------------------------------------------------------
# build_tri_mesh_circ_vert
# -----------------------------------------------------------------------------


def test_build_tri_mesh_circ_vert(icosahedron):
    ico = icosahedron
    vertices = ico["vertices"]
    vertex_of_cell = ico["vertex_of_cell"]
    num_cell = vertex_of_cell.shape[1]
    num_vert = vertices.shape[0]

    # Circumcenters (equilateral triangles -> normalised centroid)
    circ = vertices[vertex_of_cell].mean(axis=0)
    lon_circ, lat_circ = cart2lonlat(circ)
    lon_vert, lat_vert = cart2lonlat(vertices)
    rng = np.random.default_rng(3)
    elevation_circ = rng.uniform(0.0, 3_000.0, num_cell)

    tri_vert, tri_face = radtopo.build_tri_mesh_circ_vert(
        lon_circ, lat_circ, elevation_circ, lon_vert, lat_vert,
        ico["cells_of_vertex"], ico["neighbor_cell_index"])

    # Vertices: circumcenters followed by ICON vertices
    assert tri_vert.shape == (num_cell + num_vert, 3)
    np.testing.assert_array_equal(tri_vert[:num_cell, 0], lon_circ)
    np.testing.assert_array_equal(tri_vert[:num_cell, 1], lat_circ)
    np.testing.assert_array_equal(tri_vert[:num_cell, 2], elevation_circ)
    np.testing.assert_array_equal(tri_vert[num_cell:, 0], lon_vert)
    np.testing.assert_array_equal(tri_vert[num_cell:, 1], lat_vert)
    for idx_vert in range(num_vert):
        cells = ico["cells_of_vertex"][:, idx_vert]
        cells = cells[cells != -2]
        np.testing.assert_allclose(tri_vert[num_cell + idx_vert, 2],
                                   elevation_circ[cells].mean())

    # Faces: one triangle per pair of neighbouring cells around each vertex
    assert tri_face.dtype == np.uint32
    assert tri_face.shape == (3 * num_cell, 3)
    for idx_vert in range(num_vert):
        faces_vert = tri_face[tri_face[:, 0] == num_cell + idx_vert]
        assert faces_vert.shape[0] == 5
        cells = set(ico["cells_of_vertex"][:, idx_vert]) - {-2}
        assert set(faces_vert[:, 1:].ravel()) == cells
    assert np.all(tri_face[:, 1:] < num_cell)

    # Orientation and closedness
    pts = lonlat2cart(tri_vert[:, 0], tri_vert[:, 1])
    assert np.all(is_ccw(pts, tri_face))
    assert is_closed_manifold(tri_face)


# -----------------------------------------------------------------------------
# refine_tri_mesh
# -----------------------------------------------------------------------------


@pytest.mark.parametrize("n", [1, 2, 3, 6])
def test_refine_tri_mesh(icosahedron, n):
    ico = icosahedron
    vertices = ico["vertices"]
    vertex_of_cell = ico["vertex_of_cell"]
    num_vert, num_cell = vertices.shape[0], vertex_of_cell.shape[1]
    num_edge = ico["edge_vertices"].shape[1]

    vertices_child, faces_child = radtopo.refine_tri_mesh(
        vertices, vertex_of_cell, ico["edge_of_cell"], ico["edge_vertices"], n)

    # Sizes
    num_vert_expected = (num_vert + num_edge * (n - 1) + num_cell * (n - 1) *
                         (n - 2) // 2)
    assert vertices_child.shape == (num_vert_expected, 3)
    assert faces_child.shape == (num_cell * n**2, 3)
    assert faces_child.dtype == np.uint32
    assert faces_child.max() == num_vert_expected - 1

    # Vertices: parent vertices unchanged, all on unit sphere, no duplicates
    np.testing.assert_allclose(vertices_child[:num_vert], vertices)
    np.testing.assert_allclose(np.linalg.norm(vertices_child, axis=1), 1.0)
    dist = np.linalg.norm(vertices_child[:, None, :] -
                          vertices_child[None, :, :],
                          axis=2)
    np.fill_diagonal(dist, np.inf)
    assert dist.min() > 1e-6
    assert np.unique(faces_child).size == num_vert_expected

    # Faces: orientation preserved, closed surface, Euler characteristic
    assert np.all(is_ccw(vertices_child, faces_child))
    assert is_closed_manifold(faces_child)
    num_edge_child = faces_child.shape[0] * 3 // 2
    assert (vertices_child.shape[0] - num_edge_child +
            faces_child.shape[0]) == 2

    # Child triangles are located within their parent triangle
    for idx_cell in range(num_cell):
        p0, p1, p2 = vertices[vertex_of_cell[:, idx_cell]]
        normals = np.array(
            [np.cross(p0, p1),
             np.cross(p1, p2),
             np.cross(p2, p0)])  # inward normals of great circle planes
        faces = faces_child[idx_cell * n**2:(idx_cell + 1) * n**2]
        pts = vertices_child[np.unique(faces)]
        assert pts.shape[0] == (n + 1) * (n + 2) // 2
        assert np.all(pts @ normals.T > -1e-12)

    # All child triangles have a similar size (projection of the linearly
    # subdivided parent triangle onto the sphere causes some variation)
    p0, p1, p2 = (vertices_child[faces_child[:, i]] for i in range(3))
    area = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
    assert area.max() / area.min() < 2.0


def test_refine_tri_mesh_edges_shared(icosahedron):
    # Refined vertices on a parent edge are shared by both adjacent cells
    ico = icosahedron
    n = 4
    num_vert = ico["vertices"].shape[0]
    vertices_child, faces_child = radtopo.refine_tri_mesh(
        ico["vertices"], ico["vertex_of_cell"], ico["edge_of_cell"],
        ico["edge_vertices"], n)
    for idx_edge in range(ico["edge_vertices"].shape[1]):
        idx_edge_vert = np.arange(num_vert + idx_edge * (n - 1),
                                  num_vert + (idx_edge + 1) * (n - 1))
        cells = np.where((ico["edge_of_cell"] == idx_edge).any(axis=0))[0]
        assert cells.size == 2
        for idx_cell in cells:
            faces = faces_child[idx_cell * n**2:(idx_cell + 1) * n**2]
            assert np.isin(idx_edge_vert, faces).all()
        # Edge vertices lie on great circle between the edge's end points
        v0, v1 = ico["vertices"][ico["edge_vertices"][:, idx_edge]]
        np.testing.assert_allclose(
            vertices_child[idx_edge_vert] @ np.cross(v0, v1), 0.0, atol=1e-12)


# -----------------------------------------------------------------------------
# assign_points_to_tiles
# -----------------------------------------------------------------------------


def test_assign_points_to_tiles_known_points():
    lon = np.deg2rad([-180.0, -170.0, 179.9, 180.0, 0.0, 15.0, -160.0])
    lat = np.deg2rad([90.0, 85.0, -89.9, -90.0, 0.0, 47.0, 50.0])
    idx = radtopo.assign_points_to_tiles(lon, lat, 18, 18)
    idx_lon = [0, 0, 17, 17, 9, 9, 1]
    idx_lat = [0, 0, 17, 17, 9, 4, 4]
    np.testing.assert_array_equal(idx, np.array(idx_lat) * 18 + idx_lon)


@pytest.mark.parametrize("num_tile_lon, num_tile_lat", [(18, 18), (12, 6),
                                                        (1, 1), (180, 90)])
def test_assign_points_to_tiles_consistent_with_tile_name(
        num_tile_lon, num_tile_lat):
    rng = np.random.default_rng(4)
    lon = rng.uniform(-np.pi, np.pi, 500)
    lat = rng.uniform(-np.pi / 2.0, np.pi / 2.0, 500)
    idx = radtopo.assign_points_to_tiles(lon, lat, num_tile_lon, num_tile_lat)
    assert idx.min() >= 0
    assert idx.max() < num_tile_lon * num_tile_lat
    extent_lon = 360.0 / num_tile_lon
    extent_lat = 180.0 / num_tile_lat
    for k in range(lon.size):
        i, j = idx[k] % num_tile_lon, idx[k] // num_tile_lon
        lon_west = -180.0 + i * extent_lon
        lat_north = 90.0 - j * extent_lat
        assert lon_west <= np.rad2deg(lon[k]) <= lon_west + extent_lon
        assert lat_north - extent_lat <= np.rad2deg(lat[k]) <= lat_north


@pytest.mark.parametrize("lon, lat, num_tile_lon, num_tile_lat", [
    ([0.0, 1.0], [0.0], 18, 18),
    ([3.2], [0.0], 18, 18),
    ([0.0], [-1.6], 18, 18),
    ([0.0], [0.0], 7, 18),
    ([0.0], [0.0], 18, 0),
])
def test_assign_points_to_tiles_invalid(lon, lat, num_tile_lon, num_tile_lat):
    with pytest.raises(ValueError):
        radtopo.assign_points_to_tiles(np.array(lon), np.array(lat),
                                       num_tile_lon, num_tile_lat)


# -----------------------------------------------------------------------------
# get_tile_name
# -----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "idx_tile_lon, idx_tile_lat, num_tile_lon, num_tile_lat, name, lon_add", [
        (0, 0, 18, 18, "N90-N80_W180-W160", 0.0),
        (17, 17, 18, 18, "S80-S90_E160-E180", 0.0),
        (9, 9, 18, 18, "N00-S10_E000-E020", 0.0),
        (8, 4, 18, 18, "N50-N40_W020-E000", 0.0),
        (0, 1, 12, 6, "N60-N30_W180-W150", 0.0),
        (11, 5, 12, 6, "S60-S90_E150-E180", 0.0),
        (18, 3, 18, 18, "N60-N50_W180-W160", 360.0),
        (-1, 3, 18, 18, "N60-N50_E160-E180", -360.0),
        (12, 2, 12, 6, "N30-N00_W180-W150", 360.0),
        (25, 2, 12, 6, "N30-N00_W150-W120", 720.0),
    ])
def test_get_tile_name(idx_tile_lon, idx_tile_lat, num_tile_lon, num_tile_lat,
                       name, lon_add):
    result = radtopo.get_tile_name(idx_tile_lon, idx_tile_lat, num_tile_lon,
                                   num_tile_lat)
    assert result[0] == name
    assert result[1] == lon_add


def test_get_tile_name_all_tiles_unique():
    names = {
        radtopo.get_tile_name(i, j, 18, 18)[0]
        for i in range(18)
        for j in range(18)
    }
    assert len(names) == 18 * 18


@pytest.mark.parametrize(
    "idx_tile_lon, idx_tile_lat, num_tile_lon, "
    "num_tile_lat", [
        (0, -1, 18, 18),
        (0, 18, 18, 18),
        (0, 0, 7, 18),
        (0, 0, 18, 7),
        (0, 0, 0, 18),
    ])
def test_get_tile_name_invalid(idx_tile_lon, idx_tile_lat, num_tile_lon,
                               num_tile_lat):
    with pytest.raises(ValueError):
        radtopo.get_tile_name(idx_tile_lon, idx_tile_lat, num_tile_lon,
                              num_tile_lat)


# -----------------------------------------------------------------------------
# interp_bilinear
# -----------------------------------------------------------------------------


def bilinear_function(x, y):
    return 3.0 + 2.0 * x - 1.5 * y + 0.5 * x * y


@pytest.fixture
def regular_grid():
    x_axis = np.linspace(-2.0, 3.0, 26)
    y_axis = np.linspace(4.0, 1.0, 16)  # descending (like DEM latitudes)
    xx, yy = np.meshgrid(x_axis, y_axis)
    data = bilinear_function(xx, yy).astype(np.float32)
    return data, x_axis, y_axis


def test_interp_bilinear_grid_points(regular_grid):
    data, x_axis, y_axis = regular_grid
    xx, yy = np.meshgrid(x_axis, y_axis)
    result = radtopo.interp_bilinear(data, x_axis, y_axis, xx.ravel(),
                                     yy.ravel(), 1e-10)
    assert result.dtype == np.float32
    np.testing.assert_allclose(result, data.ravel(), rtol=1e-6)


def test_interp_bilinear_exact_for_bilinear_function(regular_grid):
    data, x_axis, y_axis = regular_grid
    rng = np.random.default_rng(5)
    x = rng.uniform(x_axis[0], x_axis[-1], 1000)
    y = rng.uniform(y_axis[-1], y_axis[0], 1000)
    result = radtopo.interp_bilinear(data, x_axis, y_axis, x, y, 1e-10)
    np.testing.assert_allclose(result, bilinear_function(x, y), rtol=1e-5)


def test_interp_bilinear_clamping(regular_grid):
    data, x_axis, y_axis = regular_grid
    x = np.array([-10.0, 10.0, 0.5, 0.5, -10.0, 10.0])
    y = np.array([2.0, 2.0, 10.0, -10.0, 10.0, -10.0])
    x_clamped = np.clip(x, x_axis.min(), x_axis.max())
    y_clamped = np.clip(y, y_axis.min(), y_axis.max())
    result = radtopo.interp_bilinear(data, x_axis, y_axis, x, y, 1e-10)
    np.testing.assert_allclose(result,
                               bilinear_function(x_clamped, y_clamped),
                               rtol=1e-5)


def test_interp_bilinear_minimal_grid():
    data = np.array([[0.0, 1.0], [2.0, 3.0]], dtype=np.float32)
    x_axis = np.array([0.0, 1.0])
    y_axis = np.array([0.0, 1.0])
    x = np.array([0.5, 0.0, 1.0, 0.25])
    y = np.array([0.5, 1.0, 1.0, 0.75])
    result = radtopo.interp_bilinear(data, x_axis, y_axis, x, y, 1e-10)
    np.testing.assert_allclose(result, [1.5, 2.0, 3.0, 1.75], rtol=1e-6)


def test_interp_bilinear_invalid(regular_grid):
    data, x_axis, y_axis = regular_grid
    x = np.array([0.0])
    y = np.array([2.0])
    # Shape mismatch
    with pytest.raises(ValueError):
        radtopo.interp_bilinear(data, x_axis[:-1].copy(), y_axis, x, y, 1e-10)
    # Unequal size of interpolation coordinates
    with pytest.raises(ValueError):
        radtopo.interp_bilinear(data, x_axis, y_axis, np.array([0.0, 1.0]), y,
                                1e-10)
    # Less than two points along an axis
    with pytest.raises(ValueError):
        radtopo.interp_bilinear(data[:1, :].copy(), x_axis, y_axis[:1].copy(),
                                x, y, 1e-10)
    # Negative tolerance
    with pytest.raises(ValueError):
        radtopo.interp_bilinear(data, x_axis, y_axis, x, y, -1.0)
    # Irregular spacing
    x_irregular = x_axis.copy()
    x_irregular[5] += 0.01
    with pytest.raises(ValueError):
        radtopo.interp_bilinear(data, x_irregular, y_axis, x, y, 1e-10)
