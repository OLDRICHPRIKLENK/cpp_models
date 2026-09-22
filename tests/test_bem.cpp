#include "physics/bem_solver.hpp"
#include "geometry/mesher.hpp"
#include "core/config_parser.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <filesystem>

using namespace paultrap;

void test_bem_solve_and_partition_of_unity() {
    std::cout << "[TEST] Setting up surface mesh and running BEM Laplace solver..." << std::endl;
    auto cfg = core::ConfigParser::load("config.toml");

    // Use test mesh configuration for fast test execution
    cfg.mesh.mesh_file = "data/meshes/test_bem_trap.msh";
    cfg.mesh.cl_scale = 0.001; // Coarser mesh for rapid test execution

    // Generate 2D surface mesh
    geometry::Mesher mesher(cfg);
    auto stats = mesher.generate(true); // surface_only = true
    assert(stats.is_surface_only);
    assert(stats.num_elements > 0);

    physics::BemSolver solver(cfg);
    solver.solve_basis_fields();

    assert(solver.is_solved());
    assert(solver.num_elements() == stats.num_elements);

    // Test Partition of Unity: When all 6 electrodes are at 1.0V, potential must be identically 1.0V
    std::cout << "[TEST] Testing partition of unity (all electrodes at 1.0V)..." << std::endl;
    Eigen::Matrix<double, 6, 1> unity_voltages;
    unity_voltages.setOnes();
    solver.compute_combined_field(unity_voltages);

    // Sample interior points in central region [-1 mm, +1 mm]
    const int n_pts = 27;
    Eigen::Matrix3Xd test_pts(3, n_pts);
    int idx = 0;
    for (double x : {-0.001, 0.0, 0.001}) {
        for (double y : {-0.001, 0.0, 0.001}) {
            for (double z : {-0.001, 0.0, 0.001}) {
                test_pts.col(idx++) = Eigen::Vector3d(x, y, z);
            }
        }
    }

    auto nodal_unity = solver.extract_nodal_data(test_pts);
    assert(nodal_unity.potential.size() == n_pts);
    assert(nodal_unity.basis_potentials.cols() == n_pts);

    double max_unity_dev = 0.0;
    for (int i = 0; i < n_pts; ++i) {
        double dev = std::abs(nodal_unity.potential(i) - 1.0);
        if (dev > max_unity_dev) max_unity_dev = dev;
        assert(dev < 1.0e-3 && "Partition of unity violated");
    }
    std::cout << "  Max partition-of-unity deviation: " << max_unity_dev << " V" << std::endl;
    std::cout << "[PASS] Partition of unity verified: sum(phi_k) == 1.0 everywhere." << std::endl;

    // Check maximum principle: 0.0 <= phi_k <= 1.0 for each basis field
    std::cout << "[TEST] Checking maximum principle for all 6 basis fields..." << std::endl;
    for (int k = 0; k < 6; ++k) {
        for (int i = 0; i < n_pts; ++i) {
            double val = nodal_unity.basis_potentials(k, i);
            assert(val >= -1.0e-4 && "Potential cannot be negative");
            assert(val <= 1.0 + 1.0e-4 && "Potential cannot exceed 1.0V");
        }
    }
    std::cout << "[PASS] All 6 basis fields satisfy maximum principle (0 <= phi_k <= 1)." << std::endl;

    // Test symmetry: Radial 1 (+x) and Radial 3 (-x) at (x, 0, 0) vs (-x, 0, 0)
    std::cout << "[TEST] Testing geometric symmetry..." << std::endl;
    Eigen::Matrix3Xd symm_pts(3, 2);
    symm_pts.col(0) = Eigen::Vector3d( 0.001, 0.0, 0.0);
    symm_pts.col(1) = Eigen::Vector3d(-0.001, 0.0, 0.0);
    auto symm_data = solver.extract_nodal_data(symm_pts);

    double phi1_pos = symm_data.basis_potentials(0, 0); // Radial 1 at (+x)
    double phi3_neg = symm_data.basis_potentials(2, 1); // Radial 3 at (-x)
    assert(std::abs(phi1_pos - phi3_neg) < 1.0e-4 && "Radial 1 and Radial 3 symmetry violated");

    double phi5_pos = symm_data.basis_potentials(4, 0); // Top endcap at (+x)
    double phi6_pos = symm_data.basis_potentials(5, 0); // Bottom endcap at (+x)
    assert(std::abs(phi5_pos - phi6_pos) < 1.0e-4 && "Top and bottom endcap symmetry violated");
    std::cout << "[PASS] Trap geometric symmetry verified." << std::endl;

    // Test VTK Export
    std::cout << "[TEST] Testing BEM surface VTK export..." << std::endl;
    std::string test_vtk = "output/vtk/test_bem_surface.vtu";
    solver.export_vtk(test_vtk);
    assert(std::filesystem::exists(test_vtk));
    assert(std::filesystem::file_size(test_vtk) > 500);
    std::cout << "[PASS] VTK surface file exported successfully." << std::endl;
}

int main() {
    try {
        test_bem_solve_and_partition_of_unity();
        std::cout << "=== ALL BEM TESTS PASSED ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
