// C++ program to compute topographic horizon and sky view factor

#define _USE_MATH_DEFINES
#include <cstdio>
#include <iostream>
#include <sstream>  // Fortran interface
#include <cstring>  // Fortran interface
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <embree4/rtcore.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

// Namespace
#if defined(RTC_NAMESPACE_USE)
    RTC_NAMESPACE_USE
#endif

//-----------------------------------------------------------------------------
// Definition of geometries
//-----------------------------------------------------------------------------

// Point in 3D space
struct geom_point{
    double x, y, z;
    // also used for geographic coordinates: lon (x), lat (y), elevation (z)
};

// Vector in 3D space
struct geom_vector{
    double x, y, z;
};

// Vertex (for Embree)
struct Vertex{
    float x, y, z;
};

// Triangle specified by vertex indices (for Embree)
struct Triangle{
    int v0, v1, v2;
};
// Indices should be 32-bit unsigned integers according to the Embree
// documentation. However, until 2'147'483'647, the binary representation
// between signed/unsigned integers is identical.

//-----------------------------------------------------------------------------
// Functions (not dependent on Embree)
//-----------------------------------------------------------------------------

/**
 * @brief Converts degree to radian.
 * @param ang Input angle [deg].
 * @return Output angle [rad].
 */
inline double deg2rad(double ang) {
	return ((ang / 180.0) * M_PI);
}

/**
 * @brief Converts radian to degree.
 * @param ang Input angle [rad].
 * @return Output angle [deg].
 */
inline double rad2deg(double ang) {
	return ((ang / M_PI) * 180.0);
}

/**
 * @brief Computes the dot product between two vectors.
 * @param a Vector a.
 * @param b Vector b.
 * @return Resulting dot product.
 */
inline double dot_product(geom_vector a, geom_vector b) {
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

/**
 * @brief Computes cross dot product between two vectors.
 * @param a Vector a.
 * @param b Vector b.
 * @return Resulting cross product.
 */
inline geom_vector cross_product(geom_vector a, geom_vector b) {
    geom_vector c = {a.y * b.z - a.z * b.y,
                     a.z * b.x - a.x * b.z,
                     a.x * b.y - a.y * b.x};
    return c;
}

/**
 * @brief Computes the unit vector (normalised vector) of a vector in-place.
 * @param a Vector a.
 */
void unit_vector(geom_vector& a) {
    double vector_mag = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    a.x /= vector_mag;
    a.y /= vector_mag;
    a.z /= vector_mag;
}

/**
 * @brief Rotates vector v around unit vector k with a given angle.
 *
 * This function rotates vector v around a unit vector k with a given angle
 * according to the Rodrigues' rotation formula. For performance reasons,
 * trigonometric function have to be pre-computed.
 *
 * @param v Vector that should be rotated.
 * @param k Unit vector specifying the rotation axis.
 * @param ang_rot_sin Sine of the rotation angle.
 * @param ang_rot_cos Cosine of the rotation angle.
 * @return Rotated vector.
 */
inline geom_vector vector_rotation(geom_vector v, geom_vector k,
    double ang_rot_sin, double ang_rot_cos) {
    geom_vector v_rot;
    double term = dot_product(k, v) * (1.0 - ang_rot_cos);
    v_rot.x = v.x * ang_rot_cos + (k.y * v.z - k.z * v.y) * ang_rot_sin
        + k.x * term;
    v_rot.y = v.y * ang_rot_cos + (k.z * v.x - k.x * v.z) * ang_rot_sin
        + k.y * term;
    v_rot.z = v.z * ang_rot_cos + (k.x * v.y - k.y * v.x) * ang_rot_sin
        + k.z * term;
    return v_rot;
}

/**
 * @brief Returns indices that would sort an array in ascending order.
 * @param values Input values.
 * @return Indices that would sort the array.
 */
std::vector<int> sort_index(std::vector<double>& values){
	std::vector<int> index(values.size());
     for (size_t i = 0 ; i < index.size() ; i++) {
        index[i] = i;
    }
    std::sort(index.begin(), index.end(), [&](const int& a, const int& b){
        return (values[a] < values[b]);
    });
	return index;
}

/**
 * @brief Transforms geographic to ECEF coordinates in-place.
 *
 * This function transforms geographic longitude/latitude to earth-centered,
 * earth-fixed (ECEF) coordinates. A spherical Earth is assumed.
 *
 * @param points Points (lon, lat elevation) in geographic coordinates
 *               [rad, rad, m].
 * @param rad_earth Radius of Earth [m].
 */
void lonlat2ecef(std::vector<geom_point>& points, double rad_earth){
    for (size_t i = 0; i < points.size(); i++){
        double sin_lon = sin(points[i].x);
        double cos_lon = cos(points[i].x);
        double sin_lat = sin(points[i].y);
        double cos_lat = cos(points[i].y);
        double elevation = points[i].z;
        points[i].x = (rad_earth + elevation) * cos_lat * cos_lon;
        points[i].y = (rad_earth + elevation) * cos_lat * sin_lon;
        points[i].z = (rad_earth + elevation) * sin_lat;
    }
}

/**
 * @brief Transforms points from ECEF to ENU coordinates in-place.
 * @param points Points (x, y, z) in ECEF coordinates [m].
 * @param lon_orig Longitude of ENU coordinate system origin [rad].
 * @param lat_orig Latitude of ENU coordinate system origin [rad].
 * @param rad_earth Radius of Earth [m].
 */
void ecef2enu_point(std::vector<geom_point>& points, double lon_orig,
    double lat_orig, double rad_earth){
    double sin_lon = sin(lon_orig);
    double cos_lon = cos(lon_orig);
    double sin_lat = sin(lat_orig);
    double cos_lat = cos(lat_orig);
    double x_ecef_orig = rad_earth * cos(lat_orig) * cos(lon_orig);
    double y_ecef_orig = rad_earth * cos(lat_orig) * sin(lon_orig);
    double z_ecef_orig = rad_earth * sin(lat_orig);
    double x_enu, y_enu, z_enu;
    for (size_t i = 0; i < points.size(); i++){
        x_enu = - sin_lon * (points[i].x - x_ecef_orig)
            + cos_lon * (points[i].y - y_ecef_orig);
        y_enu = - sin_lat * cos_lon * (points[i].x - x_ecef_orig)
            - sin_lat * sin_lon * (points[i].y - y_ecef_orig)
            + cos_lat * (points[i].z - z_ecef_orig);
        z_enu = + cos_lat * cos_lon * (points[i].x - x_ecef_orig)
            + cos_lat * sin_lon * (points[i].y - y_ecef_orig)
            + sin_lat * (points[i].z - z_ecef_orig);
        points[i].x = x_enu;
        points[i].y = y_enu;
        points[i].z = z_enu;
    }
}

/**
 * @brief Builds the triangle mesh from the ICON grid.
 *
 * This function builds the triangle mesh from the ICON grid cell circumcenters
 * (and vertices). Two options are available:
 * - 0: Build triangle mesh solely from ICON grid cell circumcenters
 *      (non-unique triangulation of hexa- and pentagons; relatively long
 *      triangle edges can cause artefacts in horizon computation)
 * - 1: Build triangle mesh from ICON grid cell circumcenters and vertices
 *      (elevation at vertices is computed as mean from adjacent cell
 *      circumcenters; triangulation is unique and artefacts are reduced)
 *
 * @param clon Longitude of ICON grid cell circumcenters [rad].
 * @param clat Latitudes of ICON grid cell circumcenters [rad].
 * @param hsurf Elevation of ICON grid cell circumcenters [m].
 * @param vlon Longitude of ICON grid cell vertices [rad].
 * @param vlat Latitudes of ICON grid cell vertices [rad].
 * @param cells_of_vertex Indices of ICON cells adjacent to ICON vertices.
 * @param num_cell Number of ICON grid cells.
 * @param num_vertex Number of ICON grid vertices.
 * @param grid_type Grid type option for building mesh.
 * @param vertices Vertices of build triangle mesh.
 * @param vertex_of_triangle Indices of triangles' vertices.
 * @return Number of triangles in build mesh.
 */
int build_triangle_mesh(double* clon, double* clat, double* hsurf,
    double* vlon, double* vlat, int* cells_of_vertex,
    int num_cell, int num_vertex, int grid_type,
    std::vector<geom_point>& vertices, std::vector<int>& vertex_of_triangle){
    int ind_cell;
    int ind_cov;
    for (int i = 0; i < num_cell; i++){
        vertices[i].x = clon[i];
        vertices[i].y = clat[i];
        vertices[i].z = hsurf[i];
    }
    if (grid_type == 0) {
        std::cout << "Build triangle mesh solely from ICON grid cell"
            << " circumcenters\n (non-unique triangulation)" << std::endl;
        int ind_1, ind_2;
        for (int ind_vertex = 0; ind_vertex < num_vertex; ind_vertex++){
            std::vector<double> angles;
            angles.reserve(6);
            for (int j = 0; j < 6; j++){
                ind_cell = cells_of_vertex[num_vertex * j + ind_vertex];
                if (ind_cell != -2) {
                    double angle = atan2(clon[ind_cell] - vlon[ind_vertex],
                                         clat[ind_cell] - vlat[ind_vertex]);
                    // clockwise angle from positive latitude-axis (y-axis)
                    if (angle < 0.0) {
                        angle += 2.0 * M_PI;
                    }
                    angles.push_back(angle);
                }
            }
            if (angles.size() >= 3){
                // at least 3 vertices are needed to create one or multiple
                // triangles(s) from the polygon
                std::vector<int> ind_sort = sort_index(angles);
                ind_1 = 1;
                ind_2 = 2;
                for (size_t j = 0; j < (angles.size() - 2); j++){
                    ind_cov = num_vertex * ind_sort[0] + ind_vertex;
                    vertex_of_triangle.push_back(cells_of_vertex[ind_cov]);
                    ind_cov = num_vertex * ind_sort[ind_1] + ind_vertex;
                    vertex_of_triangle.push_back(cells_of_vertex[ind_cov]);
                    ind_1 ++;
                    ind_cov = num_vertex * ind_sort[ind_2] + ind_vertex;
                    vertex_of_triangle.push_back(cells_of_vertex[ind_cov]);
                    ind_2 ++;
                    // add indices of triangle's vertices in clockwise order
                }
            }
        }
    }  else {
        std::cout << "Build triangle mesh from ICON grid cell circumcenters"
            << " and vertices\n (unique triangulation)" << std::endl;
        int ind_add = num_cell;
        int ind[7] = {0, 1, 2, 3, 4, 5, 0};
        for (int ind_vertex = 0; ind_vertex < num_vertex; ind_vertex++){
            std::vector<double> angles;
            angles.reserve(6);
            double hsurf_mean = 0.0;
            for (int j = 0; j < 6; j++){
                ind_cell = cells_of_vertex[num_vertex * j + ind_vertex];
                if (ind_cell != -2) {
                    double angle = atan2(clon[ind_cell] - vlon[ind_vertex],
                                         clat[ind_cell] - vlat[ind_vertex]);
                    // clockwise angle from positive latitude-axis (y-axis)
                    if (angle < 0.0) {
                        angle += 2.0 * M_PI;
                    }
                    angles.push_back(angle);
                    hsurf_mean += hsurf[ind_cell];
                }
            }
            if (angles.size() == 6){
                vertices.push_back({vlon[ind_vertex], vlat[ind_vertex],
                                    hsurf_mean / 6.0});
                std::vector<int> ind_sort = sort_index(angles);
                for (int j = 0; j < 6; j++){
                    ind_cov = num_vertex * ind_sort[ind[j]] + ind_vertex;
                    vertex_of_triangle.push_back(cells_of_vertex[ind_cov]);
                    ind_cov = num_vertex * ind_sort[ind[j + 1]] + ind_vertex;
                    vertex_of_triangle.push_back(cells_of_vertex[ind_cov]);
                    vertex_of_triangle.push_back(ind_add);
                }
                ind_add += 1;
            }
        }
    }
    int num_triangle = vertex_of_triangle.size() / 3;
    std::cout << "Number of triangles in mesh: " << num_triangle << std::endl;
    return num_triangle;
}

/**
 * @brief Computes the sky view factor for a horizontally aligned plane.
 *
 * This function computes the sky view factor (SVF) for a horizontally aligned
 * plane. Three methods are available:
 * - Visible sky fraction / pure geometric sky view factor
 * - Sky view factor / geometric scaled with sin(horizon)
 * - Sky view factor additionally scaled with sin(horizon) /
 *   geometric scaled with sin(horizon)**2
 *
 * @param horizon_cell Horizon array [rad].
 * @param horizon_cell_len Length of the horizon array.
 * @return Sky view factor [-].
 */
double (*function_pointer)(double* horizon_cell, int horizon_cell_len);
double pure_geometric_svf(double* horizon_cell, int horizon_cell_len){
    double svf = 0.0;
    for(int i = 0; i < horizon_cell_len; i++){
        svf += (1.0 - sin(horizon_cell[i]));
    }
    svf /= (double)horizon_cell_len;
    return svf;
}
double geometric_svf_scaled_1(double* horizon_cell, int horizon_cell_len){
    double svf = 0.0;
    for(int i = 0; i < horizon_cell_len; i++){
        svf += (1.0 - (sin(horizon_cell[i]) * sin(horizon_cell[i])));
    }
    svf /= (double)horizon_cell_len;
    return svf;
}
double geometric_svf_scaled_2(double* horizon_cell, int horizon_cell_len){
    double svf = 0.0;
    for(int i = 0; i < horizon_cell_len; i++){
        svf += (1.0 - (sin(horizon_cell[i]) * sin(horizon_cell[i])
            * sin(horizon_cell[i])));
    }
    svf /= (double)horizon_cell_len;
    return svf;
}

//-----------------------------------------------------------------------------
// Functions (Embree related)
//-----------------------------------------------------------------------------

/**
 * @brief Error function for device initialiser.
 * @param userPtr
 * @param error
 * @param str
 */
void errorFunction(void* userPtr, enum RTCError error, const char* str) {
    printf("error %d: %s\n", error, str);
}

/**
 * @brief Initialises device and registers error handler
 * @return Device instance.
 */
RTCDevice initializeDevice() {
    RTCDevice device = rtcNewDevice(NULL);
    if (!device) {
        printf("error %d: cannot create device\n", rtcGetDeviceError(NULL));
    }
    rtcSetDeviceErrorFunction(device, errorFunction, NULL);
    return device;
}

/**
 * @brief Initialises the Embree scene.
 * @param device Initialised device.
 * @param vertex_of_triangle Indices of the triangle vertices.
 * @param num_triangle Number of triangles.
 * @param vertices Vertices of the triangles [m].
 * @return Embree scene.
 */
RTCScene initializeScene(RTCDevice device, int* vertex_of_triangle,
    int num_triangle, std::vector<geom_point>& vertices){

    RTCScene scene = rtcNewScene(device);
    rtcSetSceneFlags(scene, RTC_SCENE_FLAG_ROBUST);
    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

    // Vertices
    Vertex* vertices_embree = (Vertex*) rtcSetNewGeometryBuffer(geom,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(Vertex),
        vertices.size());
    for (size_t i = 0; i < vertices.size(); i++) {
        vertices_embree[i].x = (float)vertices[i].x;
        vertices_embree[i].y = (float)vertices[i].y;
        vertices_embree[i].z = (float)vertices[i].z;
    }

    // Cell (triangle) indices to vertices
    rtcSetSharedGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0,
        RTC_FORMAT_UINT3, vertex_of_triangle, 0, 3*sizeof(int), num_triangle);

    auto start = std::chrono::high_resolution_clock::now();

    // Commit geometry and scene
    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    std::cout << std::setprecision(2) << std::fixed;
    std::cout << "Building bounding volume hierarchy (BVH): " << time.count()
        << " s" << std::endl;

    return scene;

}

/**
 * @brief Ray casting with occlusion testing (hit / no hit).
 * @param scene Embree scene.
 * @param ox x-coordinate of the ray origin [m].
 * @param oy y-coordinate of the ray origin [m].
 * @param oz z-coordinate of the ray origin [m].
 * @param dx x-component of the ray direction [m].
 * @param dy y-component of the ray direction [m].
 * @param dz z-component of the ray direction [m].
 * @param dist_search Search distance for potential collision [m].
 * @return Collision status (true: hit, false: no hit).
 */
bool castRay_occluded1(RTCScene scene, float ox, float oy, float oz, float dx,
    float dy, float dz, float dist_search){
    struct RTCRay ray;
    ray.org_x = ox;
    ray.org_y = oy;
    ray.org_z = oz;
    ray.dir_x = dx;
    ray.dir_y = dy;
    ray.dir_z = dz;
    ray.tnear = 0.0;
    ray.tfar = dist_search;
    ray.mask = 1;
    rtcOccluded1(scene, &ray); // intersect ray with scene
    return (ray.tfar < 0.0);
}

/**
 * @brief Computes the terrain horizon for a specific point.
 *
 * This function computes the terrain horizon for a specific point on the
 * triangle mesh. It iteratively samples a certain azimuth direction with rays
 * until the horizon is found. For all but the first azimuth direction, the
 * elevation angle for the search is initialised with a value equal to the
 * horizon from the previous azimuth direction +/- the horizon accuracy value.
 *
 * @param ray_org_x x-coordinate of the ray origin [m].
 * @param ray_org_y y-coordinate of the ray origin [m].
 * @param ray_org_z z-coordinate of the ray origin [m].
 * @param hori_acc Horizon accuracy [rad].
 * @param dist_search Search distance for potential collision [m].
 * @param elev_ang_thresh Threshold angle for sampling in negative elevation
 *                        angle direction [rad].
 * @param scene Embree scene.
 * @param num_rays Number of rays casted.
 * @param horizon_cell Horizon array [rad].
 * @param horizon_cell_len Length of the horizon array.
 * @param azim_shift Azimuth shift for the first azimuth sector [rad].
 * @param sphere_normal Sphere normal at the point location [m].
 * @param north_direction North direction at the point location [m].
 * @param azim_sin Sine of the azimuth angle spacing.
 * @param azim_cos Cosine of the azimuth angle spacing.
 * @param elev_sin_2ha Sine of the double elevation angle spacing.
 * @param elev_cos_2ha Cosine of the double elevation angle spacing.
 */
void terrain_horizon(float ray_org_x, float ray_org_y, float ray_org_z,
    double hori_acc, float dist_search, double elev_ang_thresh,
    RTCScene scene, size_t &num_rays,
    double* horizon_cell, int horizon_cell_len,
    double azim_shift,
    geom_vector sphere_normal, geom_vector north_direction,
    double azim_sin, double azim_cos,
    double elev_sin_2ha, double elev_cos_2ha){

    // Initial ray direction
    geom_vector ray_dir;
    ray_dir.x = north_direction.x;
    ray_dir.y = north_direction.y;
    ray_dir.z = north_direction.z;

    // Shift azimuth angle in case of 'refine_factor' > 1 so that first
    // azimuth sector is centred around 0.0 deg (North)
    ray_dir = vector_rotation(ray_dir, sphere_normal, sin(-azim_shift),
        cos(-azim_shift));

    // Sample along azimuth
    double elev_ang = 0.0;
    for (int i = 0; i < horizon_cell_len; i++){

        // Rotation axis
        geom_vector rot_axis = cross_product(ray_dir, sphere_normal);
        unit_vector(rot_axis);
        // not necessarily a unit vector because vectors are mostly not
        // perpendicular

        // Find terrain horizon by iterative ray sampling
        bool hit = castRay_occluded1(scene, ray_org_x, ray_org_y,
            ray_org_z, (float)ray_dir.x, (float)ray_dir.y, (float)ray_dir.z,
            dist_search);
        num_rays += 1;
        if (hit) { // terrain hit -> increase elevation angle
            while (hit){
                elev_ang += (2.0 * hori_acc);
                ray_dir = vector_rotation(ray_dir, rot_axis, elev_sin_2ha,
                    elev_cos_2ha);
                hit = castRay_occluded1(scene, ray_org_x, ray_org_y,
                ray_org_z, (float)ray_dir.x, (float)ray_dir.y,
                (float)ray_dir.z, dist_search);
                num_rays += 1;
            }
            horizon_cell[i] = elev_ang - hori_acc;
        } else { // terrain not hit -> decrease elevation angle
            while ((!hit) && (elev_ang > elev_ang_thresh)){
                elev_ang -= (2.0 * hori_acc);
                ray_dir = vector_rotation(ray_dir, rot_axis, -elev_sin_2ha,
                    elev_cos_2ha); // sin(-x) == -sin(x), cos(x) == cos(-x)
                hit = castRay_occluded1(scene, ray_org_x, ray_org_y,
                ray_org_z, (float)ray_dir.x, (float)ray_dir.y,
                (float)ray_dir.z, dist_search);
                num_rays += 1;
            }
            horizon_cell[i] = elev_ang + hori_acc;
        }

        // Azimuthal rotation of ray direction (clockwise; first to east)
        ray_dir = vector_rotation(ray_dir, sphere_normal, -azim_sin,
            azim_cos);  // sin(-x) == -sin(x), cos(x) == cos(-x)

    }

}

//-----------------------------------------------------------------------------
// Main function
//-----------------------------------------------------------------------------

extern "C" { // Fortran interface
void horizon_svf_comp(double* clon, double* clat, double* hsurf,
    double* vlon, double* vlat,
    int* cells_of_vertex,
    double* horizon, double* skyview,
    int num_cell, int num_vertex, int azim_num,
    int grid_type, double dist_search_dp,
    double ray_org_elev, int refine_factor,
    // int svf_type){ // Fortran interface
    int svf_type, char* buffer, int* buffer_len){ // Fortran interface

    // Redirect std::cout (Fortran interface)
    std::stringstream string_stream;
    auto* old_buf = std::cout.rdbuf(string_stream.rdbuf());

    // Fixed settings
    double hori_acc = deg2rad(0.25); // horizon accuracy [rad]
    double elev_ang_thresh = deg2rad(-85.0);
    // threshold for sampling in negative elevation angle direction [rad]
    // - relevant for 'void sampling directions' at edge of mesh
    // - necessary requirement: (elev_ang_thresh - (2.0 * hori_acc)) > -90.0

    // Constants
    double rad_earth = 6371229.0;  // ICON/COSMO earth radius [m]

    // Type casting
    float dist_search = (float)dist_search_dp;

    std::cout << "------------------------------------------------------------"
        << "-------------------" << std::endl;
    std::cout << "Horizon and SVF computation with Intel Embree (v1.1)"
        << std::endl;
    std::cout << "------------------------------------------------------------"
        << "-------------------" << std::endl;

    // Build triangle mesh from ICON grid
    auto start_mesh = std::chrono::high_resolution_clock::now();
    std::vector<geom_point> vertices(num_cell);
    std::vector<int> vertex_of_triangle;
    int num_triangle = build_triangle_mesh(clon, clat, hsurf,
        vlon, vlat, cells_of_vertex,
        num_cell, num_vertex, grid_type,
        vertices, vertex_of_triangle);
    auto end_mesh = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_mesh = end_mesh - start_mesh;
    std::cout << std::setprecision(2) << std::fixed;
    std::cout << "Triangle mesh building: " << time_mesh.count() << " s"
        << std::endl;

    // In-place transformation from geographic to ECEF coordinates
    lonlat2ecef(vertices, rad_earth);

    // Earth center and North Pole in ECEF coordinates
    std::vector<geom_point> earth_centre(1);
    earth_centre[0].x = 0.0;
    earth_centre[0].y = 0.0;
    earth_centre[0].z = 0.0;
    std::vector<geom_point> north_pole(1);
    north_pole[0].x = 0.0;
    north_pole[0].y = 0.0;
    north_pole[0].z = rad_earth;

    // Origin of ENU coordinate system
    double x_orig = 0.0;
    double y_orig = 0.0;
    double z_orig = 0.0;
    for (int i = 0; i < num_cell; i++){
        x_orig += vertices[i].x;
        y_orig += vertices[i].y;
        z_orig += vertices[i].z;
    }
    double radius = sqrt(x_orig * x_orig + y_orig * y_orig + z_orig * z_orig);
    double lon_orig = atan2(y_orig, x_orig);
    double lat_orig = asin(z_orig / radius);
    // works correctly for ICON domains containing the North/South Pole and/or
    // crossing the +/- 180 deg meridian

    // In-place transformation from ECEF to ENU coordinates
    std::cout << std::setprecision(4) << std::fixed;
    std::cout << "Origin of ENU coordinate system: " << rad2deg(lat_orig)
        << " deg lat, "  << rad2deg(lon_orig) << " deg lon" << std::endl;
    ecef2enu_point(vertices, lon_orig, lat_orig, rad_earth);
    ecef2enu_point(earth_centre, lon_orig, lat_orig, rad_earth);
    ecef2enu_point(north_pole, lon_orig, lat_orig, rad_earth);

    // Build bounding volume hierarchy (BVH)
    RTCDevice device = initializeDevice();
    RTCScene scene = initializeScene(device, vertex_of_triangle.data(),
        num_triangle, vertices);

    // Evaluated trigonometric functions for rotation along azimuth/elevation
    // angle
    int horizon_cell_len = azim_num * refine_factor;
    double azim_sin = sin(deg2rad(360.0) / (double)horizon_cell_len);
    double azim_cos = cos(deg2rad(360.0) / (double)horizon_cell_len);
    double elev_sin_2ha = sin(2.0 * hori_acc);
    double elev_cos_2ha = cos(2.0 * hori_acc);
    // Note: sin(-x) == -sin(x), cos(x) == cos(-x)

    // Compute shift for azimuth angle so that first azimuth sector is
    // centred around 0.0 deg (North) in case of 'refine_factor' > 1
    double azim_shift;
    if (refine_factor == 1) {
        azim_shift = 0.0;
    } else {
        azim_shift = -(deg2rad(360.0) / (2.0 * azim_num))
            + (deg2rad(360.0) / (2.0 * (double)horizon_cell_len));
    }

    // Select algorithm for sky view factor computation
    std::cout << "Sky View Factor computation algorithm: ";
    if (svf_type == 0) {
        std::cout << "pure geometric SVF" << std::endl;
        function_pointer = pure_geometric_svf;
    } else if (svf_type == 1) {
        std::cout << "geometric scaled with sin(horizon)" << std::endl;
        function_pointer = geometric_svf_scaled_1;
    } else if (svf_type == 2) {
        std::cout << "geometric scaled with sin(horizon)**2" << std::endl;
        function_pointer = geometric_svf_scaled_2;
    }

    auto start_ray = std::chrono::high_resolution_clock::now();
    size_t num_rays = 0;

    num_rays += tbb::parallel_reduce(
    tbb::blocked_range<size_t>(0, num_cell), 0.0,
    [&](tbb::blocked_range<size_t> r, size_t num_rays) {  // parallel

    // for (size_t i = 0; i < (size_t)num_cell; i++){ // serial
    for (size_t i=r.begin(); i<r.end(); ++i) {  // parallel

        // Compute sphere normal
        geom_vector sphere_normal = {
            (vertices[i].x - earth_centre[0].x),
            (vertices[i].y - earth_centre[0].y),
            (vertices[i].z - earth_centre[0].z)
        };
        unit_vector(sphere_normal);

        // Compute north direction (orthogonal to sphere normal)
        geom_vector north_direction;
        geom_vector v_n;
        double dot_prod;
        v_n.x = north_pole[0].x - vertices[i].x;
        v_n.y = north_pole[0].y - vertices[i].y;
        v_n.z = north_pole[0].z - vertices[i].z;
        dot_prod = dot_product(v_n, sphere_normal);
        north_direction.x = v_n.x - dot_prod * sphere_normal.x;
        north_direction.y = v_n.y - dot_prod * sphere_normal.y;
        north_direction.z = v_n.z - dot_prod * sphere_normal.z;
        unit_vector(north_direction);

        // Elevate origin for ray tracing by 'safety margin'
        float ray_org_x = (float)(vertices[i].x
            + sphere_normal.x * ray_org_elev);
        float ray_org_y = (float)(vertices[i].y
            + sphere_normal.y * ray_org_elev);
        float ray_org_z = (float)(vertices[i].z
            + sphere_normal.z * ray_org_elev);
        // The origin of the ray is slightly elevated to avoid potential ray-
        // terrain collisions near the origin due to numerical imprecisions.

        double* horizon_cell = new double[horizon_cell_len];  // [rad]

        // Compute terrain horizon
        terrain_horizon(ray_org_x, ray_org_y, ray_org_z,
            hori_acc, dist_search, elev_ang_thresh,
            scene, num_rays,
            horizon_cell, horizon_cell_len,
            azim_shift,
            sphere_normal, north_direction,
            azim_sin, azim_cos,
            elev_sin_2ha, elev_cos_2ha);

        // Clip lower limit of terrain horizon values to 0.0
        for(int j = 0; j < horizon_cell_len; j++){
            if (horizon_cell[j] < 0.0){
                horizon_cell[j] = 0.0;
            }
        }

        // Compute mean horizon for sector and save in 'horizon' buffer
        for(int j = 0; j < azim_num; j++){
            double horizon_mean = 0.0;
            for(int k = 0; k < refine_factor; k++){
                horizon_mean += horizon_cell[(j * refine_factor) + k];
            }
            horizon[(j * num_cell) + i] = (rad2deg(horizon_mean)
                / (double)refine_factor);
        }

        // Compute sky view factor and save in 'skyview' buffer
        skyview[i] = function_pointer(horizon_cell, horizon_cell_len);

        delete[] horizon_cell;

    }

    return num_rays;  // parallel
    }, std::plus<size_t>());  // parallel

    auto end_ray = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_ray = end_ray - start_ray;
    std::cout << std::setprecision(2) << std::fixed;
    std::cout << "Ray tracing: " << time_ray.count() << " s" << std::endl;

    // Print number of rays needed for location and azimuth direction
    std::cout << "Number of rays shot: " << num_rays << std::endl;
    double ratio = (double)num_rays / (double)(num_cell * azim_num);
    std::cout << std::setprecision(2) << std::fixed;
    std::cout << "Average number of rays per cell and azimuth sector: "
        << ratio << std::endl;

    std::cout << "------------------------------------------------------------"
        << "-------------------" << std::endl;

    // Restore original std::cout and copy output to buffer (Fortran interface)
    std::cout.rdbuf(old_buf);
    std::string output = string_stream.str();
    int copy_len = std::min(*buffer_len, (int)output.size());
    std::memcpy(buffer, output.c_str(), copy_len);
    *buffer_len = copy_len;

}
} // Fortran interface
