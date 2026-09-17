#include "validation/metrics.hpp"

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
#include <cmath>

namespace paultrap::validation {

std::string ValidationMetrics::summary_string() const {
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << "-------------------------------------------------------\n";
    ss << "            PAUL TRAP FIELD VALIDATION METRICS        \n";
    ss << "-------------------------------------------------------\n";
    ss << "Sample Points:         " << num_samples << "\n";
    ss << "Relative L2 Error:     " << relative_l2_error << "\n";
    ss << "Absolute L2 (RMS) Err: " << absolute_l2_error << " V\n";
    ss << "Max (L_inf) Error:     " << max_infinity_error << " V\n";
    ss << "Mean Absolute Error:   " << mean_absolute_error << " V\n";
    ss << "-------------------------------------------------------";
    return ss.str();
}

void ValidationMetrics::print_summary(std::ostream& os) const {
    os << "\n" << summary_string() << "\n" << std::endl;
}

void ValidationMetrics::print_summary() const {
    print_summary(std::cout);
}

void ValidationMetrics::save_report(const std::string& filepath) const {
    std::filesystem::path out_path(filepath);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    std::ofstream out(filepath);
    if (!out.is_open()) {
        throw std::runtime_error("ValidationMetrics: Could not open report file: " + filepath);
    }

    out << summary_string() << "\n\n";
    out << "Detailed Residual Statistics:\n";
    if (residuals.size() > 0) {
        out << "Min Signed Residual: " << residuals.minCoeff() << " V\n";
        out << "Max Signed Residual: " << residuals.maxCoeff() << " V\n";
        out << "Residual Variance:   " 
            << (residuals.array() - residuals.mean()).square().mean() << " V^2\n";
    }

    std::cout << "[Validation] Saved metrics report to " << filepath << std::endl;
}

ValidationMetrics FieldValidator::compute(
    const Eigen::VectorXd& computed,
    const Eigen::VectorXd& ideal) {
    if (computed.size() == 0 || computed.size() != ideal.size()) {
        throw std::runtime_error("FieldValidator: Mismatch or zero size in potential vectors.");
    }

    const Eigen::Index n = computed.size();
    ValidationMetrics m;
    m.num_samples = static_cast<size_t>(n);
    m.residuals = computed - ideal;
    m.absolute_residuals = m.residuals.cwiseAbs();

    const double diff_l2 = m.residuals.norm();
    const double ideal_l2 = ideal.norm();

    m.relative_l2_error = (ideal_l2 > 1.0e-14) ? (diff_l2 / ideal_l2) : diff_l2;
    m.absolute_l2_error = std::sqrt(m.residuals.squaredNorm() / static_cast<double>(n));
    m.max_infinity_error = m.absolute_residuals.maxCoeff();
    m.mean_absolute_error = m.absolute_residuals.mean();

    return m;
}

ValidationMetrics FieldValidator::compute(
    const interpolation::RegularGrid& grid,
    const physics::AnalyticalTrapField& analytical) {
    const size_t n = grid.total_points();
    Eigen::VectorXd ideal(n);
    const auto& pts = grid.points();

    for (size_t i = 0; i < n; ++i) {
        ideal(i) = analytical.potential(pts(0, i), pts(1, i), pts(2, i));
    }

    return compute(grid.potential(), ideal);
}

void FieldValidator::export_residual_vtk(
    const interpolation::RegularGrid& grid,
    const ValidationMetrics& metrics,
    const std::string& filepath) {
    std::filesystem::path out_path(filepath);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    const size_t n = grid.total_points();
    if (metrics.residuals.size() != static_cast<Eigen::Index>(n)) {
        throw std::runtime_error("FieldValidator: Metrics residuals size does not match grid points.");
    }

    const auto& pts = grid.points();
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(n);
    for (size_t i = 0; i < n; ++i) {
        points->SetPoint(i, pts(0, i), pts(1, i), pts(2, i));
    }

    const int nx = grid.nx();
    const int ny = grid.ny();
    const int nz = grid.nz();

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    for (int iz = 0; iz < nz - 1; ++iz) {
        for (int iy = 0; iy < ny - 1; ++iy) {
            for (int ix = 0; ix < nx - 1; ++ix) {
                auto voxel = vtkSmartPointer<vtkVoxel>::New();
                voxel->GetPointIds()->SetId(0, grid.flat_index(ix,     iy,     iz));
                voxel->GetPointIds()->SetId(1, grid.flat_index(ix + 1, iy,     iz));
                voxel->GetPointIds()->SetId(2, grid.flat_index(ix,     iy + 1, iz));
                voxel->GetPointIds()->SetId(3, grid.flat_index(ix + 1, iy + 1, iz));
                voxel->GetPointIds()->SetId(4, grid.flat_index(ix,     iy,     iz + 1));
                voxel->GetPointIds()->SetId(5, grid.flat_index(ix + 1, iy,     iz + 1));
                voxel->GetPointIds()->SetId(6, grid.flat_index(ix,     iy + 1, iz + 1));
                voxel->GetPointIds()->SetId(7, grid.flat_index(ix + 1, iy + 1, iz + 1));
                cells->InsertNextCell(voxel);
            }
        }
    }

    auto vtk_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtk_grid->SetPoints(points);
    vtk_grid->SetCells(VTK_VOXEL, cells);

    // Absolute residual error array
    auto abs_res_array = vtkSmartPointer<vtkDoubleArray>::New();
    abs_res_array->SetName("AbsoluteResidual");
    abs_res_array->SetNumberOfValues(n);
    for (size_t i = 0; i < n; ++i) {
        abs_res_array->SetValue(i, metrics.absolute_residuals(i));
    }
    vtk_grid->GetPointData()->AddArray(abs_res_array);
    vtk_grid->GetPointData()->SetActiveScalars("AbsoluteResidual");

    // Signed residual array
    auto res_array = vtkSmartPointer<vtkDoubleArray>::New();
    res_array->SetName("SignedResidual");
    res_array->SetNumberOfValues(n);
    for (size_t i = 0; i < n; ++i) {
        res_array->SetValue(i, metrics.residuals(i));
    }
    vtk_grid->GetPointData()->AddArray(res_array);

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filepath.c_str());
    writer->SetInputData(vtk_grid);
    if (writer->Write() != 1) {
        throw std::runtime_error("FieldValidator: Failed to write residual VTK file: " + filepath);
    }

    std::cout << "[Validation] Exported residual map to " << filepath << std::endl;
}

void FieldValidator::export_analytical_vtk(
    const interpolation::RegularGrid& grid,
    const physics::AnalyticalTrapField& analytical,
    const std::string& filepath) {
    std::filesystem::path out_path(filepath);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    const size_t n = grid.total_points();
    const auto& pts = grid.points();

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(n);
    for (size_t i = 0; i < n; ++i) {
        points->SetPoint(i, pts(0, i), pts(1, i), pts(2, i));
    }

    const int nx = grid.nx();
    const int ny = grid.ny();
    const int nz = grid.nz();

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    for (int iz = 0; iz < nz - 1; ++iz) {
        for (int iy = 0; iy < ny - 1; ++iy) {
            for (int ix = 0; ix < nx - 1; ++ix) {
                auto voxel = vtkSmartPointer<vtkVoxel>::New();
                voxel->GetPointIds()->SetId(0, grid.flat_index(ix,     iy,     iz));
                voxel->GetPointIds()->SetId(1, grid.flat_index(ix + 1, iy,     iz));
                voxel->GetPointIds()->SetId(2, grid.flat_index(ix,     iy + 1, iz));
                voxel->GetPointIds()->SetId(3, grid.flat_index(ix + 1, iy + 1, iz));
                voxel->GetPointIds()->SetId(4, grid.flat_index(ix,     iy,     iz + 1));
                voxel->GetPointIds()->SetId(5, grid.flat_index(ix + 1, iy,     iz + 1));
                voxel->GetPointIds()->SetId(6, grid.flat_index(ix,     iy + 1, iz + 1));
                voxel->GetPointIds()->SetId(7, grid.flat_index(ix + 1, iy + 1, iz + 1));
                cells->InsertNextCell(voxel);
            }
        }
    }

    auto vtk_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtk_grid->SetPoints(points);
    vtk_grid->SetCells(VTK_VOXEL, cells);

    // Analytical potential array
    auto pot_array = vtkSmartPointer<vtkDoubleArray>::New();
    pot_array->SetName("IdealPotential");
    pot_array->SetNumberOfValues(n);

    // Analytical electric field array
    auto ef_array = vtkSmartPointer<vtkDoubleArray>::New();
    ef_array->SetName("IdealElectricField");
    ef_array->SetNumberOfComponents(3);
    ef_array->SetNumberOfTuples(n);

    for (size_t i = 0; i < n; ++i) {
        const double x = pts(0, i);
        const double y = pts(1, i);
        const double z = pts(2, i);

        pot_array->SetValue(i, analytical.potential(x, y, z));
        Eigen::Vector3d ef = analytical.electric_field(x, y, z);
        double e_vec[3] = {ef.x(), ef.y(), ef.z()};
        ef_array->SetTuple(i, e_vec);
    }

    vtk_grid->GetPointData()->AddArray(pot_array);
    vtk_grid->GetPointData()->SetActiveScalars("IdealPotential");

    vtk_grid->GetPointData()->AddArray(ef_array);
    vtk_grid->GetPointData()->SetActiveVectors("IdealElectricField");

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filepath.c_str());
    writer->SetInputData(vtk_grid);
    if (writer->Write() != 1) {
        throw std::runtime_error("FieldValidator: Failed to write analytical VTK file: " + filepath);
    }

    std::cout << "[Validation] Exported ideal analytical field to " << filepath << std::endl;
}

} // namespace paultrap::validation
