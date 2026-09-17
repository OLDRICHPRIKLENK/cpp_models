#pragma once

#include "core/config_parser.hpp"
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

// Forward declarations for MFEM types
namespace mfem {
class Mesh;
class GridFunction;
} // namespace mfem

namespace paultrap::physics {

struct NodalFieldData {
    Eigen::Matrix3Xd coordinates;      // 3 x N nodal coordinates (x, y, z)
    Eigen::VectorXd potential;          // N combined potential values
    Eigen::Matrix<double, 6, Eigen::Dynamic> basis_potentials; // 6 x N individual basis fields
};

class FemSolver {
public:
    explicit FemSolver(const core::TrapConfig& config);
    ~FemSolver();

    FemSolver(const FemSolver&) = delete;
    FemSolver& operator=(const FemSolver&) = delete;
    FemSolver(FemSolver&&) noexcept;
    FemSolver& operator=(FemSolver&&) noexcept;

    // Solves the 6 independent basis Dirichlet boundary value problems
    void solve_basis_fields();

    // Evaluates linear combination of basis fields for given voltages
    void compute_combined_field(const Eigen::Matrix<double, 6, 1>& voltages);
    void compute_combined_field(const core::VoltagesConfig& voltages);

    // Extracts nodal coordinates and potential values into dense Eigen structures
    [[nodiscard]] NodalFieldData extract_nodal_data() const;

    // Native VTK .vtu file exporter
    void export_vtk(const std::string& filepath) const;

    // Accessors
    [[nodiscard]] const mfem::Mesh& get_mesh() const;
    [[nodiscard]] const mfem::GridFunction& get_basis_field(int index) const; // 0..5
    [[nodiscard]] const mfem::GridFunction& get_combined_field() const;
    [[nodiscard]] bool is_solved() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

} // namespace paultrap::physics
