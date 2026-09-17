#pragma once

#include "core/config_parser.hpp"
#include <Eigen/Dense>

namespace paultrap::physics {

class AnalyticalTrapField {
public:
    explicit AnalyticalTrapField(const core::TrapConfig& config);
    AnalyticalTrapField(const core::GeometryConfig& geom, const core::VoltagesConfig& voltages);

    // Evaluate ideal potential V(x,y,z) at a single spatial point
    [[nodiscard]] double potential(double x, double y, double z) const noexcept;
    [[nodiscard]] double potential(const Eigen::Vector3d& pos) const noexcept;

    // Evaluate ideal electric field E = -grad(V) at a single point
    [[nodiscard]] Eigen::Vector3d electric_field(double x, double y, double z) const noexcept;
    [[nodiscard]] Eigen::Vector3d electric_field(const Eigen::Vector3d& pos) const noexcept;

    // Batch evaluations for spatial points (3 x N matrix of positions)
    [[nodiscard]] Eigen::VectorXd evaluate_potential(const Eigen::Matrix3Xd& points) const;
    [[nodiscard]] Eigen::Matrix3Xd evaluate_field(const Eigen::Matrix3Xd& points) const;

    // Accessors for analytical mode coefficients
    [[nodiscard]] double radial_quadrupole_voltage() const noexcept { return v_rad_quad_; }
    [[nodiscard]] double axial_quadrupole_voltage() const noexcept { return v_axial_quad_; }
    [[nodiscard]] double offset_voltage() const noexcept { return v_offset_; }

private:
    void compute_multipole_coefficients();

    core::GeometryConfig geom_;
    core::VoltagesConfig voltages_;

    double v_rad_quad_{0.0};
    double v_axial_quad_{0.0};
    double v_offset_{0.0};
    double v_dipole_x_{0.0};
    double v_dipole_y_{0.0};
    double v_dipole_z_{0.0};
    double geom_factor_axial_{0.0}; // 1 / (2 * z0^2 + r0^2)
};

} // namespace paultrap::physics
