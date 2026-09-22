#include "geometry/mesher.hpp"
#include "core/config_parser.hpp"
#include "mfem.hpp"
#include <iostream>
#include <cassert>
#include <filesystem>

using namespace paultrap;

void test_mesh_generation() {
    std::cout << "[TEST] Running Mesher with config.toml..." << std::endl;
    auto cfg = core::ConfigParser::load("config.toml");
    
    // Override with slightly coarser scale for quick testing
    cfg.mesh.cl_scale = 0.001;
    cfg.mesh.mesh_file = "data/meshes/test_trap.msh";

    geometry::Mesher mesher(cfg);
    auto stats = mesher.generate();

    assert(stats.num_vertices > 50);
    assert(stats.num_elements > 100);
    assert(stats.num_surfaces == 6);
    assert(std::filesystem::exists(stats.mesh_path));
    assert(std::filesystem::file_size(stats.mesh_path) > 1000);

    std::cout << "[PASS] Mesher successfully created mesh with "
              << stats.num_vertices << " vertices and "
              << stats.num_elements << " elements." << std::endl;

    // Verify MFEM can load it directly
    std::cout << "[TEST] Verifying MFEM compatibility of generated mesh..." << std::endl;
    mfem::Mesh mfem_mesh(stats.mesh_path.c_str(), 1, 1);
    assert(mfem_mesh.Dimension() == 3);
    assert(mfem_mesh.bdr_attributes.Size() == 6);

    for (int i = 0; i < 6; ++i) {
        assert(mfem_mesh.bdr_attributes[i] == (i + 1));
    }

    std::cout << "[PASS] MFEM loaded mesh with all 6 boundary attributes verified." << std::endl;
}

void test_mesher_validation() {
    std::cout << "[TEST] Testing Mesher invalid parameter checks..." << std::endl;
    core::GeometryConfig bad_geom;
    bad_geom.r0 = -0.01;
    core::MeshConfig mesh_cfg;

    try {
        geometry::Mesher::generate(bad_geom, mesh_cfg);
        assert(false && "Should have thrown for negative r0");
    } catch (const std::runtime_error& e) {
        std::cout << "  Caught expected exception: " << e.what() << std::endl;
    }

    std::cout << "[PASS] Mesher parameter validation verified." << std::endl;
}

void test_surface_mesh_generation() {
    std::cout << "[TEST] Running Mesher with surface_only = true..." << std::endl;
    auto cfg = core::ConfigParser::load("config.toml");
    
    cfg.mesh.cl_scale = 0.001;
    cfg.mesh.mesh_file = "data/meshes/test_surface_trap.msh";

    geometry::Mesher mesher(cfg);
    auto stats = mesher.generate(true);

    assert(stats.is_surface_only);
    assert(stats.num_vertices > 20);
    assert(stats.num_elements > 50);
    assert(stats.num_surfaces == 6);
    assert(std::filesystem::exists(stats.mesh_path));
    assert(std::filesystem::file_size(stats.mesh_path) > 500);

    std::cout << "[PASS] Mesher successfully created 2D surface mesh with "
              << stats.num_vertices << " vertices and "
              << stats.num_elements << " triangles." << std::endl;
}

int main() {
    try {
        test_mesh_generation();
        test_surface_mesh_generation();
        test_mesher_validation();
        std::cout << "=== ALL MESHER TESTS PASSED ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
