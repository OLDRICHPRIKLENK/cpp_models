#pragma once

#include "core/config_parser.hpp"
#include "interpolation/regular_grid.hpp"
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

namespace paultrap::interpolation {

class RbfInterpolator {
public:
    explicit RbfInterpolator(const core::InterpolationConfig& config);
    RbfInterpolator(const std::string& basis = "multiquadric", double smoothing = 1.0e-5);
    ~RbfInterpolator();

    RbfInterpolator(const RbfInterpolator&) = delete;
    RbfInterpolator& operator=(const RbfInterpolator&) = delete;
    RbfInterpolator(RbfInterpolator&&) noexcept;
    RbfInterpolator& operator=(RbfInterpolator&&) noexcept;

    // Fits RBF model to unstructured 3D FEM nodal data
    void fit(const Eigen::Matrix3Xd& points, const Eigen::VectorXd& values);

    // Evaluate potential at single spatial coordinate
    [[nodiscard]] double evaluate(double x, double y, double z) const noexcept;
    [[nodiscard]] double evaluate(const Eigen::Vector3d& pos) const noexcept;

    // Evaluate Electric field E = -grad(V)
    [[nodiscard]] Eigen::Vector3d electric_field(double x, double y, double z) const noexcept;
    [[nodiscard]] Eigen::Vector3d electric_field(const Eigen::Vector3d& pos) const noexcept;

    // Evaluates potential and E-field onto entire RegularGrid
    void interpolate_grid(RegularGrid& grid) const;

    // Fast KD-tree nearest neighbor query via nanoflann
    [[nodiscard]] std::vector<size_t> find_k_nearest(const Eigen::Vector3d& query, size_t k) const;

    [[nodiscard]] bool is_fitted() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

} // namespace paultrap::interpolation
