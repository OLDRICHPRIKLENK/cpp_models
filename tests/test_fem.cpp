#include "physics/fem_solver.hpp"
#include "geometry/mesher.hpp"
#include "core/config_parser.hpp"
#include <mfem.hpp>
#include <iostream>
#include <cassert>
#include <cmath>
#include <filesystem>

using namespace paultrap;

void test_fem_laplace_solve() {
    std::cout << "[TEST] Setting up mesh and running FEM Laplace solver..." << std::endl;
    auto cfg = core::ConfigParser::load("config.toml");
    
    // Use test mesh configuration for fast test execution
    cfg.mesh.mesh_file = "data/meshes/test_fem_trap.msh";
    cfg.mesh.cl_scale = 0.001;

    // Generate test mesh
    geometry::Mesher mesher(cfg);
    mesher.generate();

    physics::FemSolver solver(cfg);
    solver.solve_basis_fields();

    assert(solver.is_solved());

    // Check each basis field satisfies maximum principle: 0.0 <= phi_k <= 1.0
    for (int k = 0; k < 6; ++k) {
        const auto& phi_k = solver.get_basis_field(k);
        assert(phi_k.Min() >= -1.0e-4);
        assert(phi_k.Max() <= 1.0 + 1.0e-4);
        assert(phi_k.Max() > 0.99); // boundary reached 1.0V
    }

    std::cout << "[PASS] All 6 basis fields satisfy maximum principle (0 <= phi <= 1)." << std::endl;

    // Test Partition of Unity: When all 6 electrodes are at 1V, potential must be identically 1.0V
    std::cout << "[TEST] Testing partition of unity (all electrodes at 1.0V)..." << std::endl;
    Eigen::Matrix<double, 6, 1> unity_voltages;
    unity_voltages.setOnes();
    solver.compute_combined_field(unity_voltages);

    const auto& u_unity = solver.get_combined_field();
    const double u_min = u_unity.Min();
    const double u_max = u_unity.Max();
    std::cout << "  Unity field min: " << u_min << ", max: " << u_max << std::endl;
    assert(std::abs(u_min - 1.0) < 1.0e-3);
    assert(std::abs(u_max - 1.0) < 1.0e-3);

    std::cout << "[PASS] Partition of unity verified: sum(phi_k) == 1.0 everywhere." << std::endl;

    // Test extraction of nodal data
    std::cout << "[TEST] Extracting nodal data to dense Eigen structures..." << std::endl;
    auto nodal = solver.extract_nodal_data();
    const size_t n_nodes = nodal.coordinates.cols();
    assert(n_nodes > 50);
    assert(nodal.potential.size() == static_cast<Eigen::Index>(n_nodes));
    assert(nodal.basis_potentials.cols() == static_cast<Eigen::Index>(n_nodes));
    assert(nodal.basis_potentials.rows() == 6);

    std::cout << "[PASS] Extracted " << n_nodes << " nodal data points successfully." << std::endl;

    // Test VTK export
    std::cout << "[TEST] Testing VTK export..." << std::endl;
    const std::string vtk_test_path = "output/vtk/test_fem_solution.vtu";
    solver.export_vtk(vtk_test_path);
    assert(std::filesystem::exists(vtk_test_path));
    assert(std::filesystem::file_size(vtk_test_path) > 1000);

    std::cout << "[PASS] VTK file successfully generated." << std::endl;

    // Clean up test files
    std::filesystem::remove(cfg.mesh.mesh_file);
    std::filesystem::remove(vtk_test_path);
}

int main() {
    try {
        test_fem_laplace_solve();
        std::cout << "=== ALL FEM SOLVER TESTS PASSED ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
