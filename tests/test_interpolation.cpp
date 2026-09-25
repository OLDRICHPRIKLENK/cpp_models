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
    std::cout << "[TEST] Testing PolynomialInterpolator fitting known quadratic harmonic field..." << std::endl;

    // Test quadratic harmonic potential (9 degrees of freedom):
    // V(x,y,z) = 15.0 + 200.0*x - 300.0*y + 50.0*z 
    //          + 120.0*(x*y) - 80.0*(y*z) + 60.0*(z*x)
    //          + 4000.0*(x^2 - y^2) + 2500.0*(2*z^2 - x^2 - y^2)
    // Note that d^2V = 4000*(2 - 2) + 2500*(-2 - 2 + 4) = 0 identically everywhere.
    const int n_pts = 64;
    Eigen::Matrix3Xd pts(3, n_pts);
    Eigen::VectorXd vals(n_pts);

    auto eval_exact_pot = [](double x, double y, double z) {
        return 15.0 + 200.0 * x - 300.0 * y + 50.0 * z
               + 120.0 * (x * y) - 80.0 * (y * z) + 60.0 * (z * x)
               + 4000.0 * (x * x - y * y)
               + 2500.0 * (2.0 * z * z - x * x - y * y);
    };

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
                vals(idx) = eval_exact_pot(x, y, z);
                ++idx;
            }
        }
    }

    interpolation::PolynomialInterpolator poly(2);
    poly.fit(pts, vals);

    assert(poly.is_fitted());
    // Verification criterion: degree 2 must have exactly 9 coefficients (reduced from 10)
    assert(poly.coefficients().size() == 9);

    // Verify recovery of the harmonic expansion coefficients
    assert(std::abs(poly.coefficients()(0) - 15.0) < 1e-6);
    assert(std::abs(poly.coefficients()(1) - 200.0) < 1e-6);
    assert(std::abs(poly.coefficients()(2) - (-300.0)) < 1e-6);
    assert(std::abs(poly.coefficients()(3) - 50.0) < 1e-6);
    assert(std::abs(poly.coefficients()(4) - 120.0) < 1e-6);
    assert(std::abs(poly.coefficients()(5) - (-80.0)) < 1e-6);
    assert(std::abs(poly.coefficients()(6) - 60.0) < 1e-6);
    assert(std::abs(poly.coefficients()(7) - 4000.0) < 1e-6);
    assert(std::abs(poly.coefficients()(8) - 2500.0) < 1e-6);

    // Check evaluate on arbitrary point inside volume
    double tx = 0.0003, ty = -0.0004, tz = 0.0007;
    double exact = eval_exact_pot(tx, ty, tz);
    double fitted = poly.evaluate(tx, ty, tz);
    assert(std::abs(fitted - exact) < 1e-10);

    // Check exact gradient and electric field:
    // dV/dx = 200 + 120*y + 60*z + 8000*x - 5000*x = 200 + 120*y + 60*z + 3000*x
    // dV/dy = -300 + 120*x - 80*z - 8000*y - 5000*y = -300 + 120*x - 80*z - 13000*y
    // dV/dz = 50 - 80*y + 60*x + 10000*z
    double dV_dx = 200.0 + 120.0 * ty + 60.0 * tz + 3000.0 * tx;
    double dV_dy = -300.0 + 120.0 * tx - 80.0 * tz - 13000.0 * ty;
    double dV_dz = 50.0 - 80.0 * ty + 60.0 * tx + 10000.0 * tz;

    Eigen::Vector3d grad = poly.gradient(tx, ty, tz);
    assert(std::abs(grad.x() - dV_dx) < 1e-7);
    assert(std::abs(grad.y() - dV_dy) < 1e-7);
    assert(std::abs(grad.z() - dV_dz) < 1e-7);

    Eigen::Vector3d e_field = poly.electric_field(tx, ty, tz);
    assert(std::abs(e_field.x() - (-dV_dx)) < 1e-7);
    assert(std::abs(e_field.y() - (-dV_dy)) < 1e-7);
    assert(std::abs(e_field.z() - (-dV_dz)) < 1e-7);

    // Verify analytical Laplacian is identically 0
    assert(std::abs(poly.laplacian(tx, ty, tz)) < 1e-14);

    std::string eq_str = poly.get_equation_string();
    assert(!eq_str.empty());

    // Test saving equation report
    const std::string rep_path = "output/reports/test_poly_report.txt";
    poly.save_equation_report(rep_path);
    assert(std::filesystem::exists(rep_path));
    std::filesystem::remove(rep_path);

    std::cout << "[PASS] PolynomialInterpolator degree 2 (9 harmonic terms) correctly fits, satisfies Laplace, and computes electric fields." << std::endl;
}

void test_polynomial_interpolator_degree4() {
    std::cout << "[TEST] Testing PolynomialInterpolator fitting degree 4 harmonic polynomial..." << std::endl;

    // Test degree 4 harmonic potential:
    // V(x,y,z) = 5.0 + 10.0*x - 20.0*y + 30.0*z
    //          + 100.0*(x^2 - y^2) + 200.0*(2*z^2 - x^2 - y^2)
    //          + 50.0*(3*x^2*y - y^3) + 70.0*(x^3 - 3*x*y^2)
    //          + 1e5*(x^4 - 6*x^2*y^2 + y^4)
    //          + 2e5*(8*z^4 - 24*z^2*(x^2 + y^2) + 3*(x^2 + y^2)^2)
    const int n_side = 5;
    const int n_pts = n_side * n_side * n_side; // 125 points > 25 terms
    Eigen::Matrix3Xd pts(3, n_pts);
    Eigen::VectorXd vals(n_pts);

    auto eval_exact_d4 = [](double x, double y, double z) {
        double x2 = x * x, y2 = y * y, z2 = z * z;
        double r2 = x2 + y2;
        return 5.0 + 10.0 * x - 20.0 * y + 30.0 * z
               + 100.0 * (x2 - y2) + 200.0 * (2.0 * z2 - r2)
               + 50.0 * (3.0 * x2 * y - y * y2) + 70.0 * (x * x2 - 3.0 * x * y2)
               + 1e5 * (x2 * x2 - 6.0 * x2 * y2 + y2 * y2)
               + 2e5 * (8.0 * z2 * z2 - 24.0 * z2 * r2 + 3.0 * r2 * r2);
    };

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
                vals(idx) = eval_exact_d4(x, y, z);
                ++idx;
            }
        }
    }

    interpolation::PolynomialInterpolator poly(4);
    assert(poly.degree() == 4);
    poly.fit(pts, vals);
    assert(poly.is_fitted());
    // Verification criterion: degree 4 must have exactly 25 coefficients (reduced from 35)
    assert(poly.coefficients().size() == 25);

    // Test evaluation at an arbitrary test point
    double tx = 0.0006, ty = -0.0008, tz = 0.0004;
    double exact = eval_exact_d4(tx, ty, tz);
    double fitted = poly.evaluate(tx, ty, tz);
    assert(std::abs(fitted - exact) < 1e-8);

    // Numerical gradient verification via central finite differences (step h = 1e-6)
    constexpr double h = 1e-6;
    double dV_dx_fd = (eval_exact_d4(tx + h, ty, tz) - eval_exact_d4(tx - h, ty, tz)) / (2.0 * h);
    double dV_dy_fd = (eval_exact_d4(tx, ty + h, tz) - eval_exact_d4(tx, ty - h, tz)) / (2.0 * h);
    double dV_dz_fd = (eval_exact_d4(tx, ty, tz + h) - eval_exact_d4(tx, ty, tz - h)) / (2.0 * h);

    Eigen::Vector3d e_field = poly.electric_field(tx, ty, tz);
    assert(std::abs(e_field.x() - (-dV_dx_fd)) < 1e-4);
    assert(std::abs(e_field.y() - (-dV_dy_fd)) < 1e-4);
    assert(std::abs(e_field.z() - (-dV_dz_fd)) < 1e-4);

    // Verify analytical Laplacian is identically zero (machine precision)
    assert(std::abs(poly.laplacian(tx, ty, tz)) < 1e-12);

    // Verify numerical Laplacian (finite difference of E-field) is zero within FD truncation
    constexpr double h_fd = 1e-6;
    Eigen::Vector3d ep_x = poly.electric_field(tx + h_fd, ty, tz);
    Eigen::Vector3d em_x = poly.electric_field(tx - h_fd, ty, tz);
    Eigen::Vector3d ep_y = poly.electric_field(tx, ty + h_fd, tz);
    Eigen::Vector3d em_y = poly.electric_field(tx, ty - h_fd, tz);
    Eigen::Vector3d ep_z = poly.electric_field(tx, ty, tz + h_fd);
    Eigen::Vector3d em_z = poly.electric_field(tx, ty, tz - h_fd);
    double div_E = ((ep_x.x() - em_x.x()) + (ep_y.y() - em_y.y()) + (ep_z.z() - em_z.z())) / (2.0 * h_fd);
    assert(std::abs(div_E) < 1e-3);

    // Test fallback check on invalid degree > 4 or < 1
    interpolation::PolynomialInterpolator poly_invalid(5);
    assert(poly_invalid.degree() == 2); // Falls back to default 2
    assert(poly_invalid.coefficients().size() == 0); // Not fitted yet

    interpolation::PolynomialInterpolator poly_invalid_zero(0);
    assert(poly_invalid_zero.degree() == 2); // Falls back to default 2

    std::cout << "[PASS] PolynomialInterpolator degree 4 fitting (25 harmonic terms), exact Laplace vanishing, and degree fallback verified." << std::endl;
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

void test_solid_harmonic_properties() {
    std::cout << "[TEST] Testing solid harmonic term counts and space charge vanishing across all degrees..." << std::endl;

    // Check term counts for degrees 1, 2, 3, 4: (degree + 1)^2
    // Degree 1: 4
    // Degree 2: 9 (decreased from 10)
    // Degree 3: 16 (decreased from 20)
    // Degree 4: 25 (decreased from 35)
    const int n_pts = 64;
    Eigen::Matrix3Xd pts(3, n_pts);
    Eigen::VectorXd vals(n_pts);
    for (int i = 0; i < n_pts; ++i) {
        double x = (i % 4 - 1.5) * 0.001;
        double y = ((i / 4) % 4 - 1.5) * 0.001;
        double z = (i / 16 - 1.5) * 0.001;
        pts(0, i) = x;
        pts(1, i) = y;
        pts(2, i) = z;
        // Linear combination of harmonic terms
        vals(i) = 10.0 + 2.0 * x - 3.0 * y + 4.0 * z + 50.0 * (x * x - y * y);
    }

    for (int deg = 1; deg <= 4; ++deg) {
        interpolation::PolynomialInterpolator poly(deg);
        poly.fit(pts, vals);
        const int expected_terms = (deg + 1) * (deg + 1);
        assert(poly.coefficients().size() == expected_terms);

        // Verify that analytical Laplacian is strictly zero at various test locations
        const double test_locs[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0008, -0.0005, 0.0003},
            {-0.0007, 0.0009, -0.0006}
        };

        for (const auto& loc : test_locs) {
            double lap = poly.laplacian(loc[0], loc[1], loc[2]);
            assert(std::abs(lap) < 1e-14);
            // Space charge rho = -epsilon_0 * laplacian = 0
            double rho = -8.8541878128e-12 * lap;
            assert(std::abs(rho) < 1e-25);
        }
    }

    std::cout << "[PASS] Solid harmonic term counts and space charge vanishing verified across degrees 1 to 4." << std::endl;
}

int main() {
    try {
        test_regular_grid();
        test_polynomial_interpolator();
        test_polynomial_interpolator_degree4();
        test_solid_harmonic_properties();
        test_rbf_interpolator();
        std::cout << "=== ALL INTERPOLATION TESTS PASSED ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
