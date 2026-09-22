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

namespace {

struct PolynomialTerm {
    int px{0};
    int py{0};
    int pz{0};
    std::string name;
};

const std::vector<PolynomialTerm>& get_terms_for_degree(int degree) {
    // Degree 1: 4 terms
    static const std::vector<PolynomialTerm> terms1 = {
        {0, 0, 0, "1"},
        {1, 0, 0, "x"},
        {0, 1, 0, "y"},
        {0, 0, 1, "z"}
    };

    // Degree 2: 10 terms (quadrupole expansion)
    static const std::vector<PolynomialTerm> terms2 = {
        {0, 0, 0, "1"},
        {1, 0, 0, "x"},
        {0, 1, 0, "y"},
        {0, 0, 1, "z"},
        {2, 0, 0, "x^2"},
        {0, 2, 0, "y^2"},
        {0, 0, 2, "z^2"},
        {1, 1, 0, "x*y"},
        {0, 1, 1, "y*z"},
        {1, 0, 1, "z*x"}
    };

    // Degree 3: 20 terms
    static const std::vector<PolynomialTerm> terms3 = {
        {0, 0, 0, "1"},
        {1, 0, 0, "x"},
        {0, 1, 0, "y"},
        {0, 0, 1, "z"},
        {2, 0, 0, "x^2"},
        {0, 2, 0, "y^2"},
        {0, 0, 2, "z^2"},
        {1, 1, 0, "x*y"},
        {0, 1, 1, "y*z"},
        {1, 0, 1, "z*x"},
        {3, 0, 0, "x^3"},
        {0, 3, 0, "y^3"},
        {0, 0, 3, "z^3"},
        {2, 1, 0, "x^2*y"},
        {1, 2, 0, "x*y^2"},
        {2, 0, 1, "x^2*z"},
        {1, 0, 2, "x*z^2"},
        {0, 2, 1, "y^2*z"},
        {0, 1, 2, "y*z^2"},
        {1, 1, 1, "x*y*z"}
    };

    // Degree 4: 35 terms (octupole & dodecapole terms)
    static const std::vector<PolynomialTerm> terms4 = {
        {0, 0, 0, "1"},
        {1, 0, 0, "x"},
        {0, 1, 0, "y"},
        {0, 0, 1, "z"},
        {2, 0, 0, "x^2"},
        {0, 2, 0, "y^2"},
        {0, 0, 2, "z^2"},
        {1, 1, 0, "x*y"},
        {0, 1, 1, "y*z"},
        {1, 0, 1, "z*x"},
        {3, 0, 0, "x^3"},
        {0, 3, 0, "y^3"},
        {0, 0, 3, "z^3"},
        {2, 1, 0, "x^2*y"},
        {1, 2, 0, "x*y^2"},
        {2, 0, 1, "x^2*z"},
        {1, 0, 2, "x*z^2"},
        {0, 2, 1, "y^2*z"},
        {0, 1, 2, "y*z^2"},
        {1, 1, 1, "x*y*z"},
        {4, 0, 0, "x^4"},
        {0, 4, 0, "y^4"},
        {0, 0, 4, "z^4"},
        {3, 1, 0, "x^3*y"},
        {1, 3, 0, "x*y^3"},
        {3, 0, 1, "x^3*z"},
        {1, 0, 3, "x*z^3"},
        {0, 3, 1, "y^3*z"},
        {0, 1, 3, "y*z^3"},
        {2, 2, 0, "x^2*y^2"},
        {0, 2, 2, "y^2*z^2"},
        {2, 0, 2, "z^2*x^2"},
        {2, 1, 1, "x^2*y*z"},
        {1, 2, 1, "x*y^2*z"},
        {1, 1, 2, "x*y*z^2"}
    };

    switch (degree) {
        case 1: return terms1;
        case 3: return terms3;
        case 4: return terms4;
        case 2:
        default:
            return terms2;
    }
}

} // anonymous namespace

PolynomialInterpolator::PolynomialInterpolator(int degree) : degree_(degree) {
    if (degree_ < 1 || degree_ > 4) {
        // Degree up to 4 supported; fallback to standard degree 2 quadrupole
        degree_ = 2;
    }
}

void PolynomialInterpolator::fit(const Eigen::Matrix3Xd& points, const Eigen::VectorXd& values) {
    const Eigen::Index n = points.cols();
    const auto& terms = get_terms_for_degree(degree_);
    const Eigen::Index num_terms = static_cast<Eigen::Index>(terms.size());

    if (n < num_terms) {
        throw std::runtime_error("PolynomialInterpolator: Need at least " + 
                                 std::to_string(num_terms) + " data points for degree " + 
                                 std::to_string(degree_) + " fit.");
    }
    if (values.size() != n) {
        throw std::runtime_error("PolynomialInterpolator: Point count and value count mismatch.");
    }

    Eigen::MatrixXd A(n, num_terms);
    for (Eigen::Index i = 0; i < n; ++i) {
        const double x = points(0, i);
        const double y = points(1, i);
        const double z = points(2, i);

        const double x_pow[5] = {1.0, x, x * x, x * x * x, x * x * x * x};
        const double y_pow[5] = {1.0, y, y * y, y * y * y, y * y * y * y};
        const double z_pow[5] = {1.0, z, z * z, z * z * z, z * z * z * z};

        for (Eigen::Index j = 0; j < num_terms; ++j) {
            const auto& t = terms[j];
            A(i, j) = x_pow[t.px] * y_pow[t.py] * z_pow[t.pz];
        }
    }

    // Solve least-squares problem via Column-Pivoting Householder QR
    coefficients_ = A.colPivHouseholderQr().solve(values);
    is_fitted_ = true;
}

double PolynomialInterpolator::evaluate(double x, double y, double z) const noexcept {
    if (!is_fitted_) return 0.0;
    const auto& terms = get_terms_for_degree(degree_);
    const double x_pow[5] = {1.0, x, x * x, x * x * x, x * x * x * x};
    const double y_pow[5] = {1.0, y, y * y, y * y * y, y * y * y * y};
    const double z_pow[5] = {1.0, z, z * z, z * z * z, z * z * z * z};

    double val = 0.0;
    const Eigen::Index n = std::min(static_cast<Eigen::Index>(terms.size()), coefficients_.size());
    for (Eigen::Index j = 0; j < n; ++j) {
        const auto& t = terms[j];
        val += coefficients_(j) * x_pow[t.px] * y_pow[t.py] * z_pow[t.pz];
    }
    return val;
}

double PolynomialInterpolator::evaluate(const Eigen::Vector3d& pos) const noexcept {
    return evaluate(pos.x(), pos.y(), pos.z());
}

Eigen::Vector3d PolynomialInterpolator::gradient(double x, double y, double z) const noexcept {
    if (!is_fitted_) return Eigen::Vector3d::Zero();
    const auto& terms = get_terms_for_degree(degree_);
    const double x_pow[5] = {1.0, x, x * x, x * x * x, x * x * x * x};
    const double y_pow[5] = {1.0, y, y * y, y * y * y, y * y * y * y};
    const double z_pow[5] = {1.0, z, z * z, z * z * z, z * z * z * z};

    double dV_dx = 0.0;
    double dV_dy = 0.0;
    double dV_dz = 0.0;

    const Eigen::Index n = std::min(static_cast<Eigen::Index>(terms.size()), coefficients_.size());
    for (Eigen::Index j = 0; j < n; ++j) {
        const auto& t = terms[j];
        const double c = coefficients_(j);

        if (t.px > 0) {
            dV_dx += c * t.px * x_pow[t.px - 1] * y_pow[t.py] * z_pow[t.pz];
        }
        if (t.py > 0) {
            dV_dy += c * t.py * x_pow[t.px] * y_pow[t.py - 1] * z_pow[t.pz];
        }
        if (t.pz > 0) {
            dV_dz += c * t.pz * x_pow[t.px] * y_pow[t.py] * z_pow[t.pz - 1];
        }
    }

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

    const auto& terms = get_terms_for_degree(degree_);
    const Eigen::Index n = std::min(static_cast<Eigen::Index>(terms.size()), coefficients_.size());
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << "V(x, y, z) = " << coefficients_(0);
    for (Eigen::Index j = 1; j < n; ++j) {
        ss << "\n  + (" << coefficients_(j) << ") * " << terms[j].name;
    }
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
