#pragma once

#include "core/config_parser.hpp"
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iosfwd>

namespace paultrap::interpolation {

class RegularGrid {
public:
    explicit RegularGrid(const core::InterpolationConfig& config);
    RegularGrid(double x_min, double x_max, int nx,
                double y_min, double y_max, int ny,
                double z_min, double z_max, int nz);

    [[nodiscard]] double x_min() const noexcept { return x_min_; }
    [[nodiscard]] double x_max() const noexcept { return x_max_; }
    [[nodiscard]] double y_min() const noexcept { return y_min_; }
    [[nodiscard]] double y_max() const noexcept { return y_max_; }
    [[nodiscard]] double z_min() const noexcept { return z_min_; }
    [[nodiscard]] double z_max() const noexcept { return z_max_; }

    [[nodiscard]] int nx() const noexcept { return nx_; }
    [[nodiscard]] int ny() const noexcept { return ny_; }
    [[nodiscard]] int nz() const noexcept { return nz_; }
    [[nodiscard]] size_t total_points() const noexcept { return total_points_; }

    [[nodiscard]] double dx() const noexcept { return dx_; }
    [[nodiscard]] double dy() const noexcept { return dy_; }
    [[nodiscard]] double dz() const noexcept { return dz_; }

    // Coordinates of all grid points: 3 x N matrix
    [[nodiscard]] const Eigen::Matrix3Xd& points() const noexcept { return points_; }
    [[nodiscard]] Eigen::Vector3d point(int ix, int iy, int iz) const;
    [[nodiscard]] size_t flat_index(int ix, int iy, int iz) const noexcept;

    // Field storage
    [[nodiscard]] const Eigen::VectorXd& potential() const noexcept { return potential_; }
    [[nodiscard]] Eigen::VectorXd& potential() noexcept { return potential_; }

    [[nodiscard]] const Eigen::Matrix3Xd& electric_field() const noexcept { return electric_field_; }
    [[nodiscard]] Eigen::Matrix3Xd& electric_field() noexcept { return electric_field_; }

    // Direct export of grid to VTK .vtu file (UnstructuredGrid with VOXEL cells)
    void export_vtk(const std::string& filepath) const;

private:
    void generate_coordinates();

    double x_min_, x_max_;
    double y_min_, y_max_;
    double z_min_, z_max_;
    int nx_, ny_, nz_;
    size_t total_points_{0};
    double dx_{0.0}, dy_{0.0}, dz_{0.0};

    Eigen::Matrix3Xd points_;
    Eigen::VectorXd potential_;
    Eigen::Matrix3Xd electric_field_;
};

class PolynomialInterpolator {
public:
    explicit PolynomialInterpolator(int degree = 2);

    // Fit polynomial to unstructured 3D points and corresponding potential values
    void fit(const Eigen::Matrix3Xd& points, const Eigen::VectorXd& values);

    // Evaluate potential at (x, y, z)
    [[nodiscard]] double evaluate(double x, double y, double z) const noexcept;
    [[nodiscard]] double evaluate(const Eigen::Vector3d& pos) const noexcept;

    // Evaluate gradient grad(V)
    [[nodiscard]] Eigen::Vector3d gradient(double x, double y, double z) const noexcept;
    [[nodiscard]] Eigen::Vector3d gradient(const Eigen::Vector3d& pos) const noexcept;

    // Evaluate Electric field E = -grad(V)
    [[nodiscard]] Eigen::Vector3d electric_field(double x, double y, double z) const noexcept;
    [[nodiscard]] Eigen::Vector3d electric_field(const Eigen::Vector3d& pos) const noexcept;

    // Evaluate on entire regular grid
    void interpolate_grid(RegularGrid& grid) const;

    // Formatting and equation report logging
    [[nodiscard]] std::string get_equation_string() const;
    void print_equation(std::ostream& os) const;
    void print_equation() const;
    void save_equation_report(const std::string& filepath) const;

    [[nodiscard]] const Eigen::VectorXd& coefficients() const noexcept { return coefficients_; }
    [[nodiscard]] int degree() const noexcept { return degree_; }
    [[nodiscard]] bool is_fitted() const noexcept { return is_fitted_; }

private:
    int degree_{2};
    bool is_fitted_{false};
    Eigen::VectorXd coefficients_; // [c0, c_x, c_y, c_z, c_xx, c_yy, c_zz, c_xy, c_yz, c_zx]
};

} // namespace paultrap::interpolation
