#include "geometry/mesher.hpp"
#include "core/constants.hpp"
#include <gmsh.h>
#include <cmath>
#include <iostream>
#include <filesystem>
#include <stdexcept>

namespace paultrap::geometry {

namespace {

struct GmshSessionGuard {
    GmshSessionGuard() {
        if (!gmsh::isInitialized()) {
            gmsh::initialize();
            gmsh::option::setNumber("General.Terminal", 0); // quiet by default
        }
    }
    ~GmshSessionGuard() {
        if (gmsh::isInitialized()) {
            gmsh::finalize();
        }
    }
};

} // anonymous namespace

Mesher::Mesher(const core::TrapConfig& config)
    : geom_(config.geometry), mesh_cfg_(config.mesh) {}

MeshStats Mesher::generate(bool surface_only) {
    stats_ = generate(geom_, mesh_cfg_, surface_only);
    return stats_;
}

MeshStats Mesher::generate(const core::GeometryConfig& geom, 
                          const core::MeshConfig& mesh_cfg,
                          bool surface_only) {
    if (geom.r0 <= 0.0 || geom.z0 <= 0.0) {
        throw std::runtime_error("Mesher: Geometric dimensions r0 and z0 must be positive.");
    }
    if (mesh_cfg.cl_scale <= 0.0) {
        throw std::runtime_error("Mesher: Mesh characteristic length cl_scale must be positive.");
    }

    // Ensure output directory exists
    std::filesystem::path out_path(mesh_cfg.mesh_file);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    GmshSessionGuard guard;

    gmsh::model::add("paul_trap_3d");

    const double r0 = geom.r0;
    const double z0 = geom.z0;
    const double lc = mesh_cfg.cl_scale;

    // Set mesh options
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", mesh_cfg.cl_min);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", mesh_cfg.cl_max);
    gmsh::option::setNumber("Mesh.ElementOrder", mesh_cfg.mesh_order);
    gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);

    // Center points on bottom (-z0) and top (+z0) planes
    int p_bot_center = gmsh::model::geo::addPoint(0.0, 0.0, -z0, lc);
    int p_top_center = gmsh::model::geo::addPoint(0.0, 0.0,  z0, lc);

    // Azimuthal angles for the 4 radial quadrants centered at +x, +y, -x, -y
    // Angles for boundaries between quadrants: -pi/4, pi/4, 3pi/4, 5pi/4
    const std::array<double, 4> boundary_angles = {
        -core::PI / 4.0,
         core::PI / 4.0,
        3.0 * core::PI / 4.0,
        5.0 * core::PI / 4.0
    };

    std::array<int, 4> p_bot;
    std::array<int, 4> p_top;
    for (int i = 0; i < 4; ++i) {
        double x = r0 * std::cos(boundary_angles[i]);
        double y = r0 * std::sin(boundary_angles[i]);
        p_bot[i] = gmsh::model::geo::addPoint(x, y, -z0, lc);
        p_top[i] = gmsh::model::geo::addPoint(x, y,  z0, lc);
    }

    // Bottom and top circle arcs, plus vertical corner lines
    std::array<int, 4> arc_bot;
    std::array<int, 4> arc_top;
    std::array<int, 4> line_vert;
    for (int i = 0; i < 4; ++i) {
        int next = (i + 1) % 4;
        arc_bot[i] = gmsh::model::geo::addCircleArc(p_bot[i], p_bot_center, p_bot[next]);
        arc_top[i] = gmsh::model::geo::addCircleArc(p_top[i], p_top_center, p_top[next]);
        line_vert[i] = gmsh::model::geo::addLine(p_bot[i], p_top[i]);
    }

    // 4 radial lateral cylindrical surfaces
    std::array<int, 4> surf_radial;
    for (int i = 0; i < 4; ++i) {
        int next = (i + 1) % 4;
        // Curve loop orientation: bottom arc, vertical forward, top arc backward, vertical backward
        int loop = gmsh::model::geo::addCurveLoop({arc_bot[i], line_vert[next], -arc_top[i], -line_vert[i]});
        surf_radial[i] = gmsh::model::geo::addSurfaceFilling({loop});
    }

    // Bottom and top endcap planar disk surfaces
    int loop_bot = gmsh::model::geo::addCurveLoop({arc_bot[0], arc_bot[1], arc_bot[2], arc_bot[3]});
    int surf_bot = gmsh::model::geo::addPlaneSurface({loop_bot});

    int loop_top = gmsh::model::geo::addCurveLoop({arc_top[0], arc_top[1], arc_top[2], arc_top[3]});
    int surf_top = gmsh::model::geo::addPlaneSurface({loop_top});

    int vol = 0;
    if (!surface_only) {
        // 3D Volume enclosure
        int surf_loop = gmsh::model::geo::addSurfaceLoop({
            surf_bot, surf_top,
            surf_radial[0], surf_radial[1], surf_radial[2], surf_radial[3]
        });
        vol = gmsh::model::geo::addVolume({surf_loop});
    }

    // Synchronize geometry to Gmsh CAD model
    gmsh::model::geo::synchronize();

    // Assign Physical Groups matching core boundary constants
    // Radial 1: +x quadrant (i=0 spans [-pi/4, pi/4])
    gmsh::model::addPhysicalGroup(2, {surf_radial[0]}, core::BOUNDARY_ID_RADIAL_1, "Radial1");
    // Radial 2: +y quadrant (i=1 spans [pi/4, 3pi/4])
    gmsh::model::addPhysicalGroup(2, {surf_radial[1]}, core::BOUNDARY_ID_RADIAL_2, "Radial2");
    // Radial 3: -x quadrant (i=2 spans [3pi/4, 5pi/4])
    gmsh::model::addPhysicalGroup(2, {surf_radial[2]}, core::BOUNDARY_ID_RADIAL_3, "Radial3");
    // Radial 4: -y quadrant (i=3 spans [5pi/4, 7pi/4])
    gmsh::model::addPhysicalGroup(2, {surf_radial[3]}, core::BOUNDARY_ID_RADIAL_4, "Radial4");
    // Endcap 1: +z top planar disc
    gmsh::model::addPhysicalGroup(2, {surf_top}, core::BOUNDARY_ID_ENDCAP_1, "Endcap1");
    // Endcap 2: -z bottom planar disc
    gmsh::model::addPhysicalGroup(2, {surf_bot}, core::BOUNDARY_ID_ENDCAP_2, "Endcap2");

    if (!surface_only) {
        // 3D domain volume
        gmsh::model::addPhysicalGroup(3, {vol}, core::DOMAIN_VOLUME_ID, "Domain");
        // Generate 3D tetrahedral mesh
        gmsh::model::mesh::generate(3);
    } else {
        // Generate 2D surface triangle mesh only
        gmsh::model::mesh::generate(2);
    }

    // Save mesh file
    gmsh::write(mesh_cfg.mesh_file);

    // Gather mesh statistics
    std::vector<std::size_t> node_tags;
    std::vector<double> coords;
    std::vector<double> param_coords;
    gmsh::model::mesh::getNodes(node_tags, coords, param_coords);

    size_t num_elements = 0;
    if (!surface_only) {
        std::vector<int> elem_types;
        std::vector<std::vector<std::size_t>> elem_tags;
        std::vector<std::vector<std::size_t>> elem_node_tags;
        gmsh::model::mesh::getElements(elem_types, elem_tags, elem_node_tags, 3, vol);

        for (const auto& tags : elem_tags) {
            num_elements += tags.size();
        }
    } else {
        std::vector<int> elem_types;
        std::vector<std::vector<std::size_t>> elem_tags;
        std::vector<std::vector<std::size_t>> elem_node_tags;
        gmsh::model::mesh::getElements(elem_types, elem_tags, elem_node_tags, 2);

        for (const auto& tags : elem_tags) {
            num_elements += tags.size();
        }
    }

    MeshStats stats;
    stats.num_vertices = node_tags.size();
    stats.num_elements = num_elements;
    stats.num_surfaces = 6;
    stats.mesh_path = mesh_cfg.mesh_file;
    stats.is_surface_only = surface_only;

    if (surface_only) {
        std::cout << "[Mesher] Generated 2D surface mesh: " << stats.num_vertices << " vertices, "
                  << stats.num_elements << " triangles -> " << stats.mesh_path << std::endl;
    } else {
        std::cout << "[Mesher] Generated 3D volume mesh: " << stats.num_vertices << " vertices, "
                  << stats.num_elements << " tetrahedra -> " << stats.mesh_path << std::endl;
    }

    return stats;
}

} // namespace paultrap::geometry
