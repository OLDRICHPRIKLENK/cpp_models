#pragma once

#include "core/config_parser.hpp"
#include "physics/fem_solver.hpp" // For NodalFieldData
#include <Eigen/Dense>
#include <memory>
#include <string>

namespace paultrap::physics {

class BemSolver {
public:
    explicit BemSolver(const core::TrapConfig& config);
    ~BemSolver();

    BemSolver(const BemSolver&) = delete;
    BemSolver& operator=(const BemSolver&) = delete;
    BemSolver(BemSolver&&) noexcept;
    BemSolver& operator=(BemSolver&&) noexcept;

    // Solves the 6 independent basis Boundary Integral problems
    void solve_basis_fields();

    // Computes total surface potential and charge density for given voltages
    void compute_combined_field(const Eigen::Matrix<double, 6, 1>& voltages);
    void compute_combined_field(const core::VoltagesConfig& voltages);

    // Direct evaluation of BEM solution at arbitrary 3D spatial points (e.g. RegularGrid points)
    [[nodiscard]] NodalFieldData extract_nodal_data(const Eigen::Matrix3Xd& eval_points) const;

    // Native VTK .vtu file exporter for 2D surface triangle mesh
    void export_vtk(const std::string& filepath) const;

    // Accessors
    [[nodiscard]] bool is_solved() const noexcept;
    [[nodiscard]] size_t num_elements() const noexcept;
    [[nodiscard]] const Eigen::VectorXd& get_surface_charges() const;
    [[nodiscard]] const Eigen::VectorXd& get_surface_voltages() const;
    [[nodiscard]] const Eigen::Matrix<double, 6, Eigen::Dynamic>& get_basis_charges() const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

} // namespace paultrap::physics
