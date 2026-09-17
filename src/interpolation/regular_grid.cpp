#include "interpolation/regular_grid.hpp"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkVoxel.h>
#include <vtkCellArray.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace paultrap::interpolation {

// ---------------------------------------------------------------------------
// RegularGrid Implementation
// ---------------------------------------------------------------------------

RegularGrid::RegularGrid(const core::InterpolationConfig& config)
    : RegularGrid(config.grid_x_min, config.grid_x_max, config.grid_nx,
                  config.grid_y_min, config.grid_y_max, config.grid_ny,
                  config.grid_z_min, config.grid_z_max, config.grid_nz) {}

RegularGrid::RegularGrid(double x_min, double x_max, int nx,
                         double y_min, double y_max, int ny,
                         double z_min, double z_max, int nz)
    : x_min_(x_min), x_max_(x_max),
      y_min_(y_min), y_max_(y_max),
      z_min_(z_min), z_max_(z_max),
      nx_(nx), ny_(ny), nz_(nz) {
    if (nx < 2 || ny < 2 || nz < 2) {
        throw std::runtime_error("RegularGrid: Resolution in each dimension must be at least 2.");
    }
    if (x_min >= x_max || y_min >= y_max || z_min >= z_max) {
        throw std::runtime_error("RegularGrid: Min bounds must be strictly less than max bounds.");
    }

    total_points_ = static_cast<size_t>(nx_) * ny_ * nz_;
    dx_ = (x_max_ - x_min_) / (nx_ - 1);
    dy_ = (y_max_ - y_min_) / (ny_ - 1);
    dz_ = (z_max_ - z_min_) / (nz_ - 1);

    generate_coordinates();

    potential_.resize(total_points_);
    potential_.setZero();

    electric_field_.resize(3, total_points_);
    electric_field_.setZero();
}

void RegularGrid::generate_coordinates() {
    points_.resize(3, total_points_);

    for (int iz = 0; iz < nz_; ++iz) {
        const double z = z_min_ + iz * dz_;
        for (int iy = 0; iy < ny_; ++iy) {
            const double y = y_min_ + iy * dy_;
            for (int ix = 0; ix < nx_; ++ix) {
                const double x = x_min_ + ix * dx_;
                const size_t idx = flat_index(ix, iy, iz);
                points_(0, idx) = x;
                points_(1, idx) = y;
                points_(2, idx) = z;
            }
        }
    }
}

size_t RegularGrid::flat_index(int ix, int iy, int iz) const noexcept {
    return static_cast<size_t>(ix) + static_cast<size_t>(nx_) * (iy + ny_ * iz);
}

Eigen::Vector3d RegularGrid::point(int ix, int iy, int iz) const {
    const size_t idx = flat_index(ix, iy, iz);
    return points_.col(idx);
}

void RegularGrid::export_vtk(const std::string& filepath) const {
    std::filesystem::path out_path(filepath);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(total_points_);

    for (size_t i = 0; i < total_points_; ++i) {
        points->SetPoint(i, points_(0, i), points_(1, i), points_(2, i));
    }

    // Connect voxels
    auto cells = vtkSmartPointer<vtkCellArray>::New();

    for (int iz = 0; iz < nz_ - 1; ++iz) {
        for (int iy = 0; iy < ny_ - 1; ++iy) {
            for (int ix = 0; ix < nx_ - 1; ++ix) {
                auto voxel = vtkSmartPointer<vtkVoxel>::New();
                // VTK_VOXEL vertex ordering
                voxel->GetPointIds()->SetId(0, flat_index(ix,     iy,     iz));
                voxel->GetPointIds()->SetId(1, flat_index(ix + 1, iy,     iz));
                voxel->GetPointIds()->SetId(2, flat_index(ix,     iy + 1, iz));
                voxel->GetPointIds()->SetId(3, flat_index(ix + 1, iy + 1, iz));
                voxel->GetPointIds()->SetId(4, flat_index(ix,     iy,     iz + 1));
                voxel->GetPointIds()->SetId(5, flat_index(ix + 1, iy,     iz + 1));
                voxel->GetPointIds()->SetId(6, flat_index(ix,     iy + 1, iz + 1));
                voxel->GetPointIds()->SetId(7, flat_index(ix + 1, iy + 1, iz + 1));
                cells->InsertNextCell(voxel);
            }
        }
    }

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    grid->SetCells(VTK_VOXEL, cells);

    // Potential scalar field
    auto pot_array = vtkSmartPointer<vtkDoubleArray>::New();
    pot_array->SetName("Potential");
    pot_array->SetNumberOfValues(total_points_);
    for (size_t i = 0; i < total_points_; ++i) {
        pot_array->SetValue(i, potential_(i));
    }
    grid->GetPointData()->AddArray(pot_array);
    grid->GetPointData()->SetActiveScalars("Potential");

    // Electric field vector array
    auto ef_array = vtkSmartPointer<vtkDoubleArray>::New();
    ef_array->SetName("ElectricField");
    ef_array->SetNumberOfComponents(3);
    ef_array->SetNumberOfTuples(total_points_);
    for (size_t i = 0; i < total_points_; ++i) {
        double e_vec[3] = {electric_field_(0, i), electric_field_(1, i), electric_field_(2, i)};
        ef_array->SetTuple(i, e_vec);
    }
    grid->GetPointData()->AddArray(ef_array);
    grid->GetPointData()->SetActiveVectors("ElectricField");

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filepath.c_str());
    writer->SetInputData(grid);
    if (writer->Write() != 1) {
        throw std::runtime_error("RegularGrid: Failed to write VTK file: " + filepath);
    }

    std::cout << "[RegularGrid] Exported grid (" << nx_ << "x" << ny_ << "x" << nz_ 
              << ") to " << filepath << std::endl;
}

// ---------------------------------------------------------------------------
// PolynomialInterpolator Implementation
// ---------------------------------------------------------------------------

PolynomialInterpolator::PolynomialInterpolator(int degree) : degree_(degree) {
    if (degree_ != 2) {
        // Degree 2 is primary for 3D quadrupole; higher degrees can be added if needed
        degree_ = 2;
    }
}

void PolynomialInterpolator::fit(const Eigen::Matrix3Xd& points, const Eigen::VectorXd& values) {
    const Eigen::Index n = points.cols();
    if (n < 10) {
        throw std::runtime_error("PolynomialInterpolator: Need at least 10 data points for 3D quadratic fit.");
    }
    if (values.size() != n) {
        throw std::runtime_error("PolynomialInterpolator: Point count and value count mismatch.");
    }

    // Degree 2 model: V = c0 + c1*x + c2*y + c3*z + c4*x^2 + c5*y^2 + c6*z^2 + c7*x*y + c8*y*z + c9*z*x
    // Design matrix A (N x 10)
    Eigen::MatrixXd A(n, 10);
    for (Eigen::Index i = 0; i < n; ++i) {
        const double x = points(0, i);
        const double y = points(1, i);
        const double z = points(2, i);

        A(i, 0) = 1.0;
        A(i, 1) = x;
        A(i, 2) = y;
        A(i, 3) = z;
        A(i, 4) = x * x;
        A(i, 5) = y * y;
        A(i, 6) = z * z;
        A(i, 7) = x * y;
        A(i, 8) = y * z;
        A(i, 9) = z * x;
    }

    // Solve least-squares problem via Column-Pivoting Householder QR
    coefficients_ = A.colPivHouseholderQr().solve(values);
    is_fitted_ = true;
}

double PolynomialInterpolator::evaluate(double x, double y, double z) const noexcept {
    if (!is_fitted_) return 0.0;
    const auto& c = coefficients_;
    return c(0) 
         + c(1) * x + c(2) * y + c(3) * z 
         + c(4) * x * x + c(5) * y * y + c(6) * z * z 
         + c(7) * x * y + c(8) * y * z + c(9) * z * x;
}

double PolynomialInterpolator::evaluate(const Eigen::Vector3d& pos) const noexcept {
    return evaluate(pos.x(), pos.y(), pos.z());
}

Eigen::Vector3d PolynomialInterpolator::gradient(double x, double y, double z) const noexcept {
    if (!is_fitted_) return Eigen::Vector3d::Zero();
    const auto& c = coefficients_;
    // dV/dx = c1 + 2*c4*x + c7*y + c9*z
    const double dV_dx = c(1) + 2.0 * c(4) * x + c(7) * y + c(9) * z;
    // dV/dy = c2 + 2*c5*y + c7*x + c8*z
    const double dV_dy = c(2) + 2.0 * c(5) * y + c(7) * x + c(8) * z;
    // dV/dz = c3 + 2*c6*z + c8*y + c9*x
    const double dV_dz = c(3) + 2.0 * c(6) * z + c(8) * y + c(9) * x;

    return {dV_dx, dV_dy, dV_dz};
}

Eigen::Vector3d PolynomialInterpolator::gradient(const Eigen::Vector3d& pos) const noexcept {
    return gradient(pos.x(), pos.y(), pos.z());
}

Eigen::Vector3d PolynomialInterpolator::electric_field(double x, double y, double z) const noexcept {
    return -gradient(x, y, z);
}

Eigen::Vector3d PolynomialInterpolator::electric_field(const Eigen::Vector3d& pos) const noexcept {
    return -gradient(pos);
}

void PolynomialInterpolator::interpolate_grid(RegularGrid& grid) const {
    if (!is_fitted_) {
        throw std::runtime_error("PolynomialInterpolator: Cannot interpolate grid before fitting.");
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

std::string PolynomialInterpolator::get_equation_string() const {
    if (!is_fitted_) {
        return "Polynomial: Not fitted yet.";
    }

    const auto& c = coefficients_;
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << "V(x, y, z) = " << c(0) << "\n"
       << "  + (" << c(1) << ") * x\n"
       << "  + (" << c(2) << ") * y\n"
       << "  + (" << c(3) << ") * z\n"
       << "  + (" << c(4) << ") * x^2\n"
       << "  + (" << c(5) << ") * y^2\n"
       << "  + (" << c(6) << ") * z^2\n"
       << "  + (" << c(7) << ") * x*y\n"
       << "  + (" << c(8) << ") * y*z\n"
       << "  + (" << c(9) << ") * z*x";
    return ss.str();
}

void PolynomialInterpolator::print_equation(std::ostream& os) const {
    os << "\n=======================================================\n";
    os << "           FITTED 3D POLYNOMIAL FIELD EQUATION        \n";
    os << "=======================================================\n";
    os << get_equation_string() << "\n";
    os << "=======================================================\n" << std::endl;
}

void PolynomialInterpolator::print_equation() const {
    print_equation(std::cout);
}

void PolynomialInterpolator::save_equation_report(const std::string& filepath) const {
    std::filesystem::path out_path(filepath);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    std::ofstream out(filepath);
    if (!out.is_open()) {
        throw std::runtime_error("PolynomialInterpolator: Could not open report file: " + filepath);
    }

    out << "=======================================================\n";
    out << "       3D PAUL TRAP FIELD SOLVER: POLYNOMIAL FIT       \n";
    out << "=======================================================\n\n";
    out << get_equation_string() << "\n\n";

    if (is_fitted_) {
        const auto& c = coefficients_;
        out << "--- Analytical Multipole Mode Equivalents ---\n";
        out << "Radial Quadrupole (c_xx - c_yy)/2: " << (c(4) - c(5)) / 2.0 << " V/m^2\n";
        out << "Axial Quadrupole c_zz:             " << c(6) << " V/m^2\n";
        out << "Laplacian residual (2*c_xx + 2*c_yy + 2*c_zz): "
            << 2.0 * (c(4) + c(5) + c(6)) << " V/m^2\n";
    }

    std::cout << "[PolynomialInterpolator] Saved equation report to " << filepath << std::endl;
}

} // namespace paultrap::interpolation
