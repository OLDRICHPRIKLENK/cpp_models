#pragma once

#include "interpolation/regular_grid.hpp"
#include "physics/analytical.hpp"
#include <Eigen/Dense>
#include <string>
#include <iosfwd>

namespace paultrap::validation {

struct ValidationMetrics {
    double relative_l2_error{0.0};
    double absolute_l2_error{0.0};
    double max_infinity_error{0.0};
    double mean_absolute_error{0.0};
    size_t num_samples{0};

    Eigen::VectorXd residuals; // Pointwise signed difference: V_computed - V_ideal
    Eigen::VectorXd absolute_residuals; // Pointwise absolute difference: |V_computed - V_ideal|

    [[nodiscard]] std::string summary_string() const;
    void print_summary(std::ostream& os) const;
    void print_summary() const;
    void save_report(const std::string& filepath) const;
};

class FieldValidator {
public:
    // Compute metrics comparing computed field against ideal baseline vectors
    [[nodiscard]] static ValidationMetrics compute(
        const Eigen::VectorXd& computed_potentials,
        const Eigen::VectorXd& ideal_potentials
    );

    // Compute metrics on a RegularGrid against AnalyticalTrapField
    [[nodiscard]] static ValidationMetrics compute(
        const interpolation::RegularGrid& grid,
        const physics::AnalyticalTrapField& analytical
    );

    // Export standalone residual error map to .vtu file (Strict SRP)
    static void export_residual_vtk(
        const interpolation::RegularGrid& grid,
        const ValidationMetrics& metrics,
        const std::string& filepath
    );

    // Export standalone analytical ideal field to .vtu file (Strict SRP)
    static void export_analytical_vtk(
        const interpolation::RegularGrid& grid,
        const physics::AnalyticalTrapField& analytical,
        const std::string& filepath
    );
};

} // namespace paultrap::validation
