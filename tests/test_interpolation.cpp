#include "interpolation/regular_grid.hpp"
#include "interpolation/rbf_interpolator.hpp"
#include "core/config_parser.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <filesystem>

using namespace paultrap;

void test_regular_grid() {
    std::cout << "[TEST] Testing RegularGrid..." << std::endl;
    interpolation::RegularGrid grid(-0.002, 0.002, 5,
                                   -0.002, 0.002, 5,
                                   -0.002, 0.002, 5);

    assert(grid.total_points() == 125);
    assert(grid.nx() == 5 && grid.ny() == 5 && grid.nz() == 5);
    assert(std::abs(grid.dx() - 0.001) < 1e-12);

    auto p0 = grid.point(0, 0, 0);
    assert(std::abs(p0.x() - (-0.002)) < 1e-12);
    assert(std::abs(p0.y() - (-0.002)) < 1e-12);
    assert(std::abs(p0.z() - (-0.002)) < 1e-12);

    auto p_mid = grid.point(2, 2, 2);
    assert(std::abs(p_mid.x()) < 1e-12);
    assert(std::abs(p_mid.y()) < 1e-12);
    assert(std::abs(p_mid.z()) < 1e-12);

    std::cout << "[PASS] RegularGrid creation and coordinates verified." << std::endl;
}

void test_polynomial_interpolator() {
    std::cout << "[TEST] Testing PolynomialInterpolator fitting known quadratic field..." << std::endl;

    // Test quadratic potential: V(x,y,z) = 15.0 + 200.0*x - 300.0*y + 4000.0*x^2 - 4000.0*y^2 + 8000.0*z^2
    const int n_pts = 64;
    Eigen::Matrix3Xd pts(3, n_pts);
    Eigen::VectorXd vals(n_pts);

    int idx = 0;
    for (int ix = -2; ix <= 1; ++ix) {
        for (int iy = -2; iy <= 1; ++iy) {
            for (int iz = -2; iz <= 1; ++iz) {
                double x = ix * 0.0005;
                double y = iy * 0.0005;
                double z = iz * 0.0005;
                pts(0, idx) = x;
                pts(1, idx) = y;
                pts(2, idx) = z;
                vals(idx) = 15.0 + 200.0 * x - 300.0 * y + 4000.0 * x * x - 4000.0 * y * y + 8000.0 * z * z;
                ++idx;
            }
        }
    }

    interpolation::PolynomialInterpolator poly(2);
    poly.fit(pts, vals);

    assert(poly.is_fitted());

    // Check evaluate on arbitrary point
    double tx = 0.0003, ty = -0.0004, tz = 0.0007;
    double exact = 15.0 + 200.0 * tx - 300.0 * ty + 4000.0 * tx * tx - 4000.0 * ty * ty + 8000.0 * tz * tz;
    double fitted = poly.evaluate(tx, ty, tz);
    assert(std::abs(fitted - exact) < 1e-8);

    // Check electric field: Ex = -(200 + 8000*x), Ey = -(-300 - 8000*y), Ez = -16000*z
    Eigen::Vector3d e_field = poly.electric_field(tx, ty, tz);
    assert(std::abs(e_field.x() - (-(200.0 + 8000.0 * tx))) < 1e-6);
    assert(std::abs(e_field.y() - (-(-300.0 - 8000.0 * ty))) < 1e-6);
    assert(std::abs(e_field.z() - (-16000.0 * tz)) < 1e-6);

    std::string eq_str = poly.get_equation_string();
    assert(!eq_str.empty());

    std::cout << "[PASS] PolynomialInterpolator correctly fits and computes electric fields." << std::endl;
}

void test_polynomial_interpolator_degree4() {
    std::cout << "[TEST] Testing PolynomialInterpolator fitting degree 4 polynomial..." << std::endl;

    // Test degree 4 potential with mixed and quartic terms
    const int n_side = 5;
    const int n_pts = n_side * n_side * n_side; // 125 points > 35 terms
    Eigen::Matrix3Xd pts(3, n_pts);
    Eigen::VectorXd vals(n_pts);

    int idx = 0;
    for (int ix = -2; ix <= 2; ++ix) {
        for (int iy = -2; iy <= 2; ++iy) {
            for (int iz = -2; iz <= 2; ++iz) {
                double x = ix * 0.001;
                double y = iy * 0.001;
                double z = iz * 0.001;
                pts(0, idx) = x;
                pts(1, idx) = y;
                pts(2, idx) = z;
                vals(idx) = 5.0 + 10.0 * x - 20.0 * y + 30.0 * z
                          + 100.0 * x * x - 200.0 * y * y + 300.0 * z * z
                          + 1e6 * (x * x * x * x) - 2e6 * (y * y * y * y) + 3e6 * (x * x * y * y);
                ++idx;
            }
        }
    }

    interpolation::PolynomialInterpolator poly(4);
    assert(poly.degree() == 4);
    poly.fit(pts, vals);
    assert(poly.is_fitted());
    assert(poly.coefficients().size() == 35);

    // Test evaluation at an arbitrary test point
    double tx = 0.0006, ty = -0.0008, tz = 0.0004;
    double exact = 5.0 + 10.0 * tx - 20.0 * ty + 30.0 * tz
                 + 100.0 * tx * tx - 200.0 * ty * ty + 300.0 * tz * tz
                 + 1e6 * std::pow(tx, 4) - 2e6 * std::pow(ty, 4) + 3e6 * (tx * tx * ty * ty);
    double fitted = poly.evaluate(tx, ty, tz);
    assert(std::abs(fitted - exact) < 1e-6);

    // Electric field: Ex = -dV/dx, Ey = -dV/dy, Ez = -dV/dz
    double dV_dx = 10.0 + 200.0 * tx + 4e6 * std::pow(tx, 3) + 6e6 * tx * (ty * ty);
    double dV_dy = -20.0 - 400.0 * ty - 8e6 * std::pow(ty, 3) + 6e6 * (tx * tx) * ty;
    double dV_dz = 30.0 + 600.0 * tz;
    Eigen::Vector3d e_field = poly.electric_field(tx, ty, tz);
    assert(std::abs(e_field.x() - (-dV_dx)) < 1e-5);
    assert(std::abs(e_field.y() - (-dV_dy)) < 1e-5);
    assert(std::abs(e_field.z() - (-dV_dz)) < 1e-5);

    // Test fallback check on invalid degree > 4
    interpolation::PolynomialInterpolator poly_invalid(5);
    assert(poly_invalid.degree() == 2); // Falls back to default 2

    std::cout << "[PASS] PolynomialInterpolator degree 4 fitting, derivatives, and degree fallback verified." << std::endl;
}

void test_rbf_interpolator() {
    std::cout << "[TEST] Testing RbfInterpolator with ALGLIB and nanoflann..." << std::endl;

    const int n_pts = 64;
    Eigen::Matrix3Xd pts(3, n_pts);
    Eigen::VectorXd vals(n_pts);

    int idx = 0;
    for (int ix = -2; ix <= 1; ++ix) {
        for (int iy = -2; iy <= 1; ++iy) {
            for (int iz = -2; iz <= 1; ++iz) {
                double x = ix * 0.0005;
                double y = iy * 0.0005;
                double z = iz * 0.0005;
                pts(0, idx) = x;
                pts(1, idx) = y;
                pts(2, idx) = z;
                vals(idx) = 100.0 * (x * x - y * y);
                ++idx;
            }
        }
    }

    interpolation::RbfInterpolator rbf("multiquadric", 1e-5);
    rbf.fit(pts, vals);

    assert(rbf.is_fitted());

    // Test KD-tree
    auto neighbors = rbf.find_k_nearest(Eigen::Vector3d(0.0, 0.0, 0.0), 4);
    assert(neighbors.size() == 4);

    // Test evaluation
    double tx = 0.00025, ty = 0.00025, tz = 0.0;
    double val_rbf = rbf.evaluate(tx, ty, tz);
    double val_exact = 100.0 * (tx * tx - ty * ty);
    assert(std::abs(val_rbf - val_exact) < 1e-4);

    std::cout << "[PASS] RbfInterpolator successfully evaluated and tested." << std::endl;
}

int main() {
    try {
        test_regular_grid();
        test_polynomial_interpolator();
        test_polynomial_interpolator_degree4();
        test_rbf_interpolator();
        std::cout << "=== ALL INTERPOLATION TESTS PASSED ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
