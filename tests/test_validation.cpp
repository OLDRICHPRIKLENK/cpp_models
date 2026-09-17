#include "validation/metrics.hpp"
#include "interpolation/regular_grid.hpp"
#include "physics/analytical.hpp"
#include "core/config_parser.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <filesystem>

using namespace paultrap;

void test_metrics_calculation() {
    std::cout << "[TEST] Testing FieldValidator error norms calculation..." << std::endl;

    const int n = 100;
    Eigen::VectorXd ideal = Eigen::VectorXd::LinSpaced(n, 10.0, 20.0);
    Eigen::VectorXd comp = ideal;

    auto m_zero = validation::FieldValidator::compute(comp, ideal);
    assert(m_zero.relative_l2_error == 0.0);
    assert(m_zero.absolute_l2_error == 0.0);
    assert(m_zero.max_infinity_error == 0.0);

    // Add constant offset of 0.5V
    comp = ideal.array() + 0.5;
    auto m_offset = validation::FieldValidator::compute(comp, ideal);
    assert(std::abs(m_offset.absolute_l2_error - 0.5) < 1e-12);
    assert(std::abs(m_offset.max_infinity_error - 0.5) < 1e-12);
    assert(std::abs(m_offset.mean_absolute_error - 0.5) < 1e-12);

    std::string summary = m_offset.summary_string();
    assert(!summary.empty());

    std::cout << "[PASS] Metric calculations and error norms verified." << std::endl;
}

void test_validation_vtk_exports() {
    std::cout << "[TEST] Testing standalone VTK exports for residual and analytical fields..." << std::endl;
    auto cfg = core::ConfigParser::load("config.toml");
    physics::AnalyticalTrapField analytical(cfg);

    interpolation::RegularGrid grid(-0.001, 0.001, 5,
                                   -0.001, 0.001, 5,
                                   -0.001, 0.001, 5);

    // Fill grid potential with analytical + small artificial noise
    for (size_t i = 0; i < grid.total_points(); ++i) {
        auto pt = grid.points().col(i);
        grid.potential()(i) = analytical.potential(pt) + 0.01;
    }

    auto metrics = validation::FieldValidator::compute(grid, analytical);

    const std::string res_path = "output/vtk/test_residual_error.vtu";
    const std::string ana_path = "output/vtk/test_analytical_field.vtu";

    validation::FieldValidator::export_residual_vtk(grid, metrics, res_path);
    validation::FieldValidator::export_analytical_vtk(grid, analytical, ana_path);

    assert(std::filesystem::exists(res_path));
    assert(std::filesystem::file_size(res_path) > 500);

    assert(std::filesystem::exists(ana_path));
    assert(std::filesystem::file_size(ana_path) > 500);

    std::filesystem::remove(res_path);
    std::filesystem::remove(ana_path);

    std::cout << "[PASS] Standalone VTK export verified." << std::endl;
}

int main() {
    try {
        test_metrics_calculation();
        test_validation_vtk_exports();
        std::cout << "=== ALL VALIDATION TESTS PASSED ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
