#include "core/config_parser.hpp"
#include "core/constants.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace paultrap;

void test_load_valid_config() {
    std::cout << "[TEST] Loading config.toml..." << std::endl;
    auto cfg = core::ConfigParser::load("config.toml");

    assert(std::abs(cfg.geometry.r0 - 0.005) < 1e-9);
    assert(std::abs(cfg.geometry.z0 - 0.005) < 1e-9);
    assert(cfg.mesh.mesh_file == "data/meshes/trap.msh");
    assert(cfg.fem.order >= 1);
    assert(cfg.interpolation.grid_nx == 21);
    assert(cfg.interpolation.grid_ny == 21);
    assert(cfg.interpolation.grid_nz == 21);

    auto v = cfg.voltages.to_vector();
    assert(v.size() == 6);
    assert(std::abs(v[0] - 100.0) < 1e-9);
    assert(std::abs(v[1] - (-100.0)) < 1e-9);
    assert(std::abs(v[4] - 50.0) < 1e-9);

    assert(cfg.solver.method == "bem");
    assert(cfg.bem.element_type == "triangle");
    assert(cfg.bem.collocation_method == "centroid");
    assert(cfg.bem.quadrature_order == 4);
    assert(cfg.output.vtk_bem_file == "output/vtk/bem_solution.vtu");

    std::cout << "[PASS] Valid config loaded successfully." << std::endl;
}

void test_invalid_config_handling() {
    std::cout << "[TEST] Testing invalid config handling..." << std::endl;
    try {
        core::ConfigParser::load("non_existent_file.toml");
        assert(false && "Should have thrown for missing file");
    } catch (const std::runtime_error& e) {
        std::cout << "  Caught expected exception: " << e.what() << std::endl;
    }

    core::TrapConfig bad_cfg;
    bad_cfg.geometry.r0 = -1.0;
    try {
        core::ConfigParser::validate(bad_cfg);
        assert(false && "Should have thrown for negative r0");
    } catch (const std::runtime_error& e) {
        std::cout << "  Caught expected exception: " << e.what() << std::endl;
    }

    core::TrapConfig bad_solver_cfg;
    bad_solver_cfg.solver.method = "invalid_solver";
    try {
        core::ConfigParser::validate(bad_solver_cfg);
        assert(false && "Should have thrown for invalid solver method");
    } catch (const std::runtime_error& e) {
        std::cout << "  Caught expected exception: " << e.what() << std::endl;
    }

    std::cout << "[PASS] Error handling verified." << std::endl;
}

int main() {
    try {
        test_load_valid_config();
        test_invalid_config_handling();
        std::cout << "=== ALL CONFIG TESTS PASSED ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with error: " << e.what() << std::endl;
        return 1;
    }
}
