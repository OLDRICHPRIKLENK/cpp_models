#include "interpolation/rbf_interpolator.hpp"

#include <interpolation.h>
#include <nanoflann.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace paultrap::interpolation {

namespace {

// Nanoflann adaptor for Eigen::Matrix3Xd point cloud
struct EigenMatrixPointCloud {
    const Eigen::Matrix3Xd& points;

    inline size_t kdtree_get_point_count() const noexcept {
        return static_cast<size_t>(points.cols());
    }

    inline double kdtree_get_pt(const size_t idx, const size_t dim) const noexcept {
        return points(static_cast<Eigen::Index>(dim), static_cast<Eigen::Index>(idx));
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const noexcept {
        return false;
    }
};

using KDTree_t = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, EigenMatrixPointCloud>,
    EigenMatrixPointCloud,
    3,
    size_t
>;

} // anonymous namespace

struct RbfInterpolator::Impl {
    std::string basis{"multiquadric"};
    double smoothing{1.0e-5};
    alglib::rbfmodel model;
    bool is_fitted{false};

    Eigen::Matrix3Xd stored_points;
    std::unique_ptr<EigenMatrixPointCloud> cloud_adaptor;
    std::unique_ptr<KDTree_t> kdtree;
};

RbfInterpolator::RbfInterpolator(const core::InterpolationConfig& config)
    : RbfInterpolator(config.rbf_basis, config.rbf_smoothing) {}

RbfInterpolator::RbfInterpolator(const std::string& basis, double smoothing)
    : pimpl_(std::make_unique<Impl>()) {
    pimpl_->basis = basis;
    pimpl_->smoothing = smoothing;
}

RbfInterpolator::~RbfInterpolator() = default;

RbfInterpolator::RbfInterpolator(RbfInterpolator&&) noexcept = default;
RbfInterpolator& RbfInterpolator::operator=(RbfInterpolator&&) noexcept = default;

void RbfInterpolator::fit(const Eigen::Matrix3Xd& points, const Eigen::VectorXd& values) {
    const Eigen::Index n = points.cols();
    if (n < 4) {
        throw std::runtime_error("RbfInterpolator: At least 4 points required for 3D RBF interpolation.");
    }
    if (values.size() != n) {
        throw std::runtime_error("RbfInterpolator: Number of points and values must match.");
    }

    // 1. Build nanoflann KD-Tree for spatial queries
    pimpl_->stored_points = points;
    pimpl_->cloud_adaptor = std::make_unique<EigenMatrixPointCloud>(EigenMatrixPointCloud{pimpl_->stored_points});
    pimpl_->kdtree = std::make_unique<KDTree_t>(
        3, *pimpl_->cloud_adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10)
    );
    pimpl_->kdtree->buildIndex();

    // 2. Estimate average inter-point distance via nearest neighbors
    double mean_dist = 0.0;
    const size_t sample_count = std::min(static_cast<size_t>(n), static_cast<size_t>(100));
    for (size_t i = 0; i < sample_count; ++i) {
        double query_pt[3] = {points(0, i), points(1, i), points(2, i)};
        size_t ret_index[2];
        double out_dist_sqr[2];
        (void)pimpl_->kdtree->knnSearch(query_pt, 2, ret_index, out_dist_sqr);
        mean_dist += std::sqrt(out_dist_sqr[1]);
    }
    mean_dist /= static_cast<double>(sample_count);
    const double r_base = std::max(mean_dist * 2.5, 1.0e-4);

    // 3. Configure ALGLIB RBF model
    alglib::rbfcreate(3, 1, pimpl_->model);

    alglib::real_2d_array xy;
    xy.setlength(n, 4);
    for (Eigen::Index i = 0; i < n; ++i) {
        xy[i][0] = points(0, i);
        xy[i][1] = points(1, i);
        xy[i][2] = points(2, i);
        xy[i][3] = values(i);
    }

    alglib::rbfsetpoints(pimpl_->model, xy);

    // Multilayer RBF model
    alglib::rbfsetalgomultilayer(pimpl_->model, r_base, 3, pimpl_->smoothing);

    alglib::rbfreport rep;
    alglib::rbfbuildmodel(pimpl_->model, rep);

    if (rep.terminationtype <= 0) {
        throw std::runtime_error("RbfInterpolator: ALGLIB RBF fitting failed with code " + 
                                 std::to_string(rep.terminationtype));
    }

    pimpl_->is_fitted = true;
    std::cout << "[RbfInterpolator] Successfully fitted RBF model with " << n 
              << " points (r_base = " << r_base << ", smoothing = " << pimpl_->smoothing 
              << ", iterations = " << rep.iterationscount << ")" << std::endl;
}

double RbfInterpolator::evaluate(double x, double y, double z) const noexcept {
    if (!pimpl_->is_fitted) return 0.0;
    try {
        return alglib::rbfcalc3(pimpl_->model, x, y, z);
    } catch (...) {
        return 0.0;
    }
}

double RbfInterpolator::evaluate(const Eigen::Vector3d& pos) const noexcept {
    return evaluate(pos.x(), pos.y(), pos.z());
}

Eigen::Vector3d RbfInterpolator::electric_field(double x, double y, double z) const noexcept {
    if (!pimpl_->is_fitted) return Eigen::Vector3d::Zero();

    // Central finite difference gradient E = -grad(V)
    constexpr double h = 1.0e-6; // 1 micrometer step
    constexpr double inv_2h = 1.0 / (2.0 * h);

    const double dV_dx = (evaluate(x + h, y, z) - evaluate(x - h, y, z)) * inv_2h;
    const double dV_dy = (evaluate(x, y + h, z) - evaluate(x, y - h, z)) * inv_2h;
    const double dV_dz = (evaluate(x, y, z + h) - evaluate(x, y, z - h)) * inv_2h;

    return {-dV_dx, -dV_dy, -dV_dz};
}

Eigen::Vector3d RbfInterpolator::electric_field(const Eigen::Vector3d& pos) const noexcept {
    return electric_field(pos.x(), pos.y(), pos.z());
}

void RbfInterpolator::interpolate_grid(RegularGrid& grid) const {
    if (!pimpl_->is_fitted) {
        throw std::runtime_error("RbfInterpolator: Cannot interpolate grid before fitting.");
    }

    const size_t n = grid.total_points();
    const auto& pts = grid.points();

    for (size_t i = 0; i < n; ++i) {
        const double x = pts(0, i);
        const double y = pts(1, i);
        const double z = pts(2, i);

        grid.potential()(i) = evaluate(x, y, z);
        grid.electric_field().col(i) = electric_field(x, y, z);
    }
}

std::vector<size_t> RbfInterpolator::find_k_nearest(const Eigen::Vector3d& query, size_t k) const {
    if (!pimpl_->kdtree) {
        throw std::runtime_error("RbfInterpolator: KDTree not built yet.");
    }

    std::vector<size_t> indices(k);
    std::vector<double> dist_sqr(k);
    double q_pt[3] = {query.x(), query.y(), query.z()};

    (void)pimpl_->kdtree->knnSearch(q_pt, k, indices.data(), dist_sqr.data());
    return indices;
}

bool RbfInterpolator::is_fitted() const noexcept {
    return pimpl_ && pimpl_->is_fitted;
}

} // namespace paultrap::interpolation
