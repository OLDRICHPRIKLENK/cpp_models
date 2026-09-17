#include "physics/analytical.hpp"
#include <cmath>
#include <stdexcept>

namespace paultrap::physics {

AnalyticalTrapField::AnalyticalTrapField(const core::TrapConfig& config)
    : geom_(config.geometry), voltages_(config.voltages) {
    compute_multipole_coefficients();
}

AnalyticalTrapField::AnalyticalTrapField(const core::GeometryConfig& geom, 
                                         const core::VoltagesConfig& voltages)
    : geom_(geom), voltages_(voltages) {
    compute_multipole_coefficients();
}

void AnalyticalTrapField::compute_multipole_coefficients() {
    if (geom_.r0 <= 0.0 || geom_.z0 <= 0.0) {
        throw std::runtime_error("AnalyticalTrapField: Dimensions r0 and z0 must be positive.");
    }

    const double v1 = voltages_.v1;
    const double v2 = voltages_.v2;
    const double v3 = voltages_.v3;
    const double v4 = voltages_.v4;
    const double v5 = voltages_.v5;
    const double v6 = voltages_.v6;

    // Radial quadrupole amplitude: ((+x and -x) - (+y and -y)) / 4
    v_rad_quad_ = ((v1 + v3) - (v2 + v4)) / 4.0;

    // Axial quadrupole amplitude: endcap avg - radial ring avg
    const double v_endcap_avg = (v5 + v6) / 2.0;
    const double v_radial_avg = (v1 + v2 + v3 + v4) / 4.0;
    v_axial_quad_ = v_endcap_avg - v_radial_avg;

    // Monopole offset
    v_offset_ = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;

    // Dipole steering terms
    v_dipole_x_ = (v1 - v3) / 2.0;
    v_dipole_y_ = (v2 - v4) / 2.0;
    v_dipole_z_ = (v5 - v6) / 2.0;

    // Axial geometry factor
    geom_factor_axial_ = 1.0 / (2.0 * geom_.z0 * geom_.z0 + geom_.r0 * geom_.r0);
}

double AnalyticalTrapField::potential(double x, double y, double z) const noexcept {
    const double r0 = geom_.r0;
    const double z0 = geom_.z0;

    // Dipole contribution
    const double v_dip = (v_dipole_x_ * x / r0) + 
                         (v_dipole_y_ * y / r0) + 
                         (v_dipole_z_ * z / z0);

    // Radial quadrupole contribution: (x^2 - y^2) / r0^2
    const double v_rq = v_rad_quad_ * (x * x - y * y) / (r0 * r0);

    // Axial quadrupole contribution: (2z^2 - x^2 - y^2) / (2*z0^2 + r0^2)
    const double v_aq = v_axial_quad_ * geom_factor_axial_ * (2.0 * z * z - x * x - y * y);

    return v_offset_ + v_dip + v_rq + v_aq;
}

double AnalyticalTrapField::potential(const Eigen::Vector3d& pos) const noexcept {
    return potential(pos.x(), pos.y(), pos.z());
}

Eigen::Vector3d AnalyticalTrapField::electric_field(double x, double y, double z) const noexcept {
    const double r0 = geom_.r0;
    const double z0 = geom_.z0;
    const double r0_sq = r0 * r0;

    // Ex = - dV/dx
    const double dV_dx = (v_dipole_x_ / r0) 
                       + (2.0 * v_rad_quad_ / r0_sq) * x 
                       - (2.0 * v_axial_quad_ * geom_factor_axial_) * x;

    // Ey = - dV/dy
    const double dV_dy = (v_dipole_y_ / r0) 
                       - (2.0 * v_rad_quad_ / r0_sq) * y 
                       - (2.0 * v_axial_quad_ * geom_factor_axial_) * y;

    // Ez = - dV/dz
    const double dV_dz = (v_dipole_z_ / z0) 
                       + (4.0 * v_axial_quad_ * geom_factor_axial_) * z;

    return {-dV_dx, -dV_dy, -dV_dz};
}

Eigen::Vector3d AnalyticalTrapField::electric_field(const Eigen::Vector3d& pos) const noexcept {
    return electric_field(pos.x(), pos.y(), pos.z());
}

Eigen::VectorXd AnalyticalTrapField::evaluate_potential(const Eigen::Matrix3Xd& points) const {
    const Eigen::Index n = points.cols();
    Eigen::VectorXd potentials(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        potentials(i) = potential(points(0, i), points(1, i), points(2, i));
    }
    return potentials;
}

Eigen::Matrix3Xd AnalyticalTrapField::evaluate_field(const Eigen::Matrix3Xd& points) const {
    const Eigen::Index n = points.cols();
    Eigen::Matrix3Xd fields(3, n);
    for (Eigen::Index i = 0; i < n; ++i) {
        fields.col(i) = electric_field(points(0, i), points(1, i), points(2, i));
    }
    return fields;
}

} // namespace paultrap::physics
