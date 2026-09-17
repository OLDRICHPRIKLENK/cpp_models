#include "physics/analytical.hpp"
#include "core/config_parser.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace paultrap;

void test_laplacian_zero() {
    std::cout << "[TEST] Verifying Laplace condition div(grad(V)) == 0..." << std::endl;
    auto cfg = core::ConfigParser::load("config.toml");
    physics::AnalyticalTrapField field(cfg);

    const double h = 1.0e-4;
    const double h2 = h * h;

    // Test across several interior spatial points
    const std::vector<Eigen::Vector3d> test_points = {
        {0.0, 0.0, 0.0},
        {0.001, 0.001, 0.001},
        {-0.0015, 0.0012, -0.0008},
        {0.0005, -0.001, 0.0015}
    };

    for (const auto& pt : test_points) {
        double v = field.potential(pt.x(), pt.y(), pt.z());
        double v_xp = field.potential(pt.x() + h, pt.y(), pt.z());
        double v_xm = field.potential(pt.x() - h, pt.y(), pt.z());
        double v_yp = field.potential(pt.x(), pt.y() + h, pt.z());
        double v_ym = field.potential(pt.x(), pt.y() - h, pt.z());
        double v_zp = field.potential(pt.x(), pt.y(), pt.z() + h);
        double v_zm = field.potential(pt.x(), pt.y(), pt.z() - h);

        double laplacian = (v_xp + v_xm - 2.0 * v) / h2
                         + (v_yp + v_ym - 2.0 * v) / h2
                         + (v_zp + v_zm - 2.0 * v) / h2;

        assert(std::abs(laplacian) < 1.0e-5);
    }

    std::cout << "[PASS] Laplace equation satisfied everywhere (Laplacian == 0)." << std::endl;
}

void test_electric_field_gradient() {
    std::cout << "[TEST] Verifying E == -grad(V)..." << std::endl;
    auto cfg = core::ConfigParser::load("config.toml");
    physics::AnalyticalTrapField field(cfg);

    const double h = 1.0e-6;
    Eigen::Vector3d pt(0.0008, -0.0007, 0.0011);

    Eigen::Vector3d e_analytical = field.electric_field(pt);

    double dV_dx = (field.potential(pt.x() + h, pt.y(), pt.z()) - field.potential(pt.x() - h, pt.y(), pt.z())) / (2.0 * h);
    double dV_dy = (field.potential(pt.x(), pt.y() + h, pt.z()) - field.potential(pt.x(), pt.y() - h, pt.z())) / (2.0 * h);
    double dV_dz = (field.potential(pt.x(), pt.y(), pt.z() + h) - field.potential(pt.x(), pt.y(), pt.z() - h)) / (2.0 * h);

    Eigen::Vector3d e_fd(-dV_dx, -dV_dy, -dV_dz);

    assert((e_analytical - e_fd).norm() < 1.0e-5);
    std::cout << "[PASS] Electric field matches -grad(V) finite difference gradient." << std::endl;
}

void test_batch_evaluation() {
    std::cout << "[TEST] Verifying batch evaluation consistency..." << std::endl;
    auto cfg = core::ConfigParser::load("config.toml");
    physics::AnalyticalTrapField field(cfg);

    Eigen::Matrix3Xd pts(3, 4);
    pts << 0.0, 0.001, -0.001,  0.0005,
           0.0, 0.001,  0.000, -0.0010,
           0.0, 0.000,  0.001,  0.0015;

    auto v_batch = field.evaluate_potential(pts);
    auto e_batch = field.evaluate_field(pts);

    assert(v_batch.size() == 4);
    assert(e_batch.cols() == 4);

    for (int i = 0; i < 4; ++i) {
        double v_single = field.potential(pts.col(i));
        Eigen::Vector3d e_single = field.electric_field(pts.col(i));
        assert(std::abs(v_batch(i) - v_single) < 1.0e-12);
        assert((e_batch.col(i) - e_single).norm() < 1.0e-12);
    }

    std::cout << "[PASS] Batch evaluations match single-point evaluations exactly." << std::endl;
}

int main() {
    try {
        test_laplacian_zero();
        test_electric_field_gradient();
        test_batch_evaluation();
        std::cout << "=== ALL ANALYTICAL TESTS PASSED ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
