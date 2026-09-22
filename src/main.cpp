#include "core/config_parser.hpp"
#include "geometry/mesher.hpp"
#include "physics/analytical.hpp"
#include "physics/fem_solver.hpp"
#include "interpolation/regular_grid.hpp"
#include "interpolation/rbf_interpolator.hpp"
#include "validation/metrics.hpp"

#include <chrono>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    try {
        std::cout << "=======================================================\n";
        std::cout << "   3D PAUL TRAP FIELD SOLVER (NATIVE C++20 / MFEM / VTK)\n";
        std::cout << "=======================================================\n";

        const auto start_total = std::chrono::high_resolution_clock::now();

        // 1. Configuration parsing
        std::string config_path = "config.toml";
        if (argc > 1) {
            config_path = argv[1];
        }
        std::cout << "[Step 1/7] Loading configuration from: " << config_path << std::endl;
        const auto config = paultrap::core::ConfigParser::load(config_path);
        std::cout << "  Loaded parameters: r0 = " << config.geometry.r0 * 1000.0 << " mm, z0 = " 
                  << config.geometry.z0 * 1000.0 << " mm" << std::endl;
        std::cout << "  Applied voltages: [" 
                  << config.voltages.v1 << ", " << config.voltages.v2 << ", "
                  << config.voltages.v3 << ", " << config.voltages.v4 << ", "
                  << config.voltages.v5 << ", " << config.voltages.v6 << "] V" << std::endl;
        std::cout << "  Polynomial fit degree: " << config.interpolation.polynomial_degree << std::endl;

        // 2. 3D Mesh Generation (Gmsh)
        std::cout << "\n[Step 2/7] Generating 3D Paul trap geometry and mesh..." << std::endl;
        paultrap::geometry::Mesher mesher(config);
        const auto mesh_stats = mesher.generate();
        std::cout << "  Generated mesh: " << mesh_stats.num_vertices << " vertices, "
                  << mesh_stats.num_elements << " tetrahedra" << std::endl;

        // 3. FEM Laplace Equation Solve (MFEM 6 Basis Fields)
        std::cout << "\n[Step 3/7] Solving 6-basis Laplace boundary value problems with MFEM..." << std::endl;
        const auto t_fem_start = std::chrono::high_resolution_clock::now();
        paultrap::physics::FemSolver fem_solver(config);
        fem_solver.solve_basis_fields();
        fem_solver.compute_combined_field(config.voltages);
        const auto t_fem_end = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> fem_elapsed = t_fem_end - t_fem_start;
        std::cout << "  FEM solution completed in " << fem_elapsed.count() << " seconds." << std::endl;

        // Export FEM Solution to VTK (SRP Dataset 1)
        fem_solver.export_vtk(config.output.vtk_fem_file);

        // 4. Extract Nodal Data for Interpolation
        std::cout << "\n[Step 4/7] Extracting dense unstructured nodal data..." << std::endl;
        const auto nodal_data = fem_solver.extract_nodal_data();
        std::cout << "  Extracted " << nodal_data.coordinates.cols() << " spatial nodal coordinates." << std::endl;

        // 5. Polynomial & RBF Spatial Interpolation
        std::cout << "\n[Step 5/7] Fitting spatial interpolators..." << std::endl;

        // 5a. Polynomial Interpolator (Quadrupole expansion)
        paultrap::interpolation::PolynomialInterpolator poly(config.interpolation.polynomial_degree);
        poly.fit(nodal_data.coordinates, nodal_data.potential);

        // Print equation to stdout and save report file
        poly.print_equation();
        poly.save_equation_report(config.output.report_polynomial_file);

        // Create Regular 3D Cartesian Grid
        paultrap::interpolation::RegularGrid reg_grid(config.interpolation);
        std::cout << "  Interpolating onto regular 3D grid (" 
                  << reg_grid.nx() << " x " << reg_grid.ny() << " x " << reg_grid.nz() 
                  << " = " << reg_grid.total_points() << " points)..." << std::endl;
        poly.interpolate_grid(reg_grid);

        // Export Interpolated Grid to VTK (SRP Dataset 2)
        reg_grid.export_vtk(config.output.vtk_interpolated_file);

        // 5b. RBF Interpolator (ALGLIB + nanoflann)
        std::cout << "  Fitting ALGLIB RBF interpolator..." << std::endl;
        paultrap::interpolation::RbfInterpolator rbf(config.interpolation);
        rbf.fit(nodal_data.coordinates, nodal_data.potential);

        // 6. Analytical Ground Truth Baseline
        std::cout << "\n[Step 6/7] Evaluating analytical 3D quadrupole baseline..." << std::endl;
        paultrap::physics::AnalyticalTrapField analytical(config);
        std::cout << "  Analytical multipole modes: V_rq = " << analytical.radial_quadrupole_voltage() 
                  << " V, V_aq = " << analytical.axial_quadrupole_voltage() 
                  << " V, V_offset = " << analytical.offset_voltage() << " V" << std::endl;

        // Export Analytical Ideal Field to VTK (SRP Dataset 3)
        paultrap::validation::FieldValidator::export_analytical_vtk(
            reg_grid, analytical, config.output.vtk_analytical_file
        );

        // 7. Validation & Rigorous Error Metrics
        std::cout << "\n[Step 7/7] Computing validation metrics and residual error map..." << std::endl;
        const auto metrics = paultrap::validation::FieldValidator::compute(reg_grid, analytical);

        // Print and save metrics report
        metrics.print_summary();
        metrics.save_report(config.output.report_metrics_file);

        // Export Dedicated Residual Error Map to VTK (SRP Dataset 4)
        paultrap::validation::FieldValidator::export_residual_vtk(
            reg_grid, metrics, config.output.vtk_residual_file
        );

        const auto end_total = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> total_elapsed = end_total - start_total;

        std::cout << "=======================================================\n";
        std::cout << "               SIMULATION COMPLETED IN " << total_elapsed.count() << " s\n";
        std::cout << "=======================================================\n";
        std::cout << "Generated Output Files (Strict SRP):\n";
        std::cout << "  1. FEM Mesh Solution:    " << config.output.vtk_fem_file << "\n";
        std::cout << "  2. Interpolated Grid:    " << config.output.vtk_interpolated_file << "\n";
        std::cout << "  3. Analytical Baseline:  " << config.output.vtk_analytical_file << "\n";
        std::cout << "  4. Residual Error Map:   " << config.output.vtk_residual_file << "\n";
        std::cout << "  5. Polynomial Equation:  " << config.output.report_polynomial_file << "\n";
        std::cout << "  6. Validation Metrics:   " << config.output.report_metrics_file << "\n";
        std::cout << "=======================================================\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n[FATAL ERROR] " << e.what() << std::endl;
        return 1;
    }
}
