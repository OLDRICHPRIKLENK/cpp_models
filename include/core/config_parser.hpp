#pragma once

#include <string>
#include <array>
#include <vector>
#include <stdexcept>
#include <filesystem>
#include <Eigen/Dense>
#include <toml++/toml.hpp>

namespace paultrap::core {

struct GeometryConfig {
    double r0{0.005};            // Cylinder radius [m]
    double z0{0.005};            // Half-length along z [m]
    double electrode_gap{0.0002}; // Gap between segmented electrodes [m]
};

struct MeshConfig {
    std::string mesh_file{"data/meshes/trap.msh"};
    double cl_scale{0.0006};     // Characteristic mesh element size [m]
    double cl_min{0.0003};
    double cl_max{0.0010};
    int mesh_order{1};
};

struct FemConfig {
    int order{1};
    int refine_steps{0};
    std::string solver_type{"cg"};
    int max_iter{3000};
    double rel_tol{1.0e-12};
    double abs_tol{1.0e-16};
};

struct VoltagesConfig {
    double v1{100.0};  // Radial electrode 1 (+x)
    double v2{-100.0}; // Radial electrode 2 (+y)
    double v3{100.0};  // Radial electrode 3 (-x)
    double v4{-100.0}; // Radial electrode 4 (-y)
    double v5{50.0};   // Endcap 1 (+z)
    double v6{50.0};   // Endcap 2 (-z)

    [[nodiscard]] Eigen::Matrix<double, 6, 1> to_vector() const noexcept {
        Eigen::Matrix<double, 6, 1> vec;
        vec << v1, v2, v3, v4, v5, v6;
        return vec;
    }

    [[nodiscard]] std::array<double, 6> to_array() const noexcept {
        return {v1, v2, v3, v4, v5, v6};
    }
};

struct InterpolationConfig {
    int polynomial_degree{2};
    std::string rbf_basis{"multiquadric"};
    double rbf_smoothing{1.0e-5};

    double grid_x_min{-0.002};
    double grid_x_max{0.002};
    double grid_y_min{-0.002};
    double grid_y_max{0.002};
    double grid_z_min{-0.002};
    double grid_z_max{0.002};

    int grid_nx{21};
    int grid_ny{21};
    int grid_nz{21};

    [[nodiscard]] size_t total_points() const noexcept {
        return static_cast<size_t>(grid_nx) * grid_ny * grid_nz;
    }
};

struct OutputConfig {
    std::string vtk_fem_file{"output/vtk/fem_solution.vtu"};
    std::string vtk_interpolated_file{"output/vtk/interpolated_grid.vtu"};
    std::string vtk_analytical_file{"output/vtk/analytical_field.vtu"};
    std::string vtk_residual_file{"output/vtk/residual_error.vtu"};
    std::string report_polynomial_file{"output/reports/polynomial_equation.txt"};
    std::string report_metrics_file{"output/reports/validation_metrics.txt"};
};

struct TrapConfig {
    GeometryConfig geometry;
    MeshConfig mesh;
    FemConfig fem;
    VoltagesConfig voltages;
    InterpolationConfig interpolation;
    OutputConfig output;
};

class ConfigParser {
public:
    static TrapConfig load(const std::string& filepath) {
        if (!std::filesystem::exists(filepath)) {
            throw std::runtime_error("ConfigParser: Configuration file not found at " + filepath);
        }

        toml::table tbl;
        try {
            tbl = toml::parse_file(filepath);
        } catch (const toml::parse_error& err) {
            throw std::runtime_error("ConfigParser: Failed to parse TOML config: " + std::string(err.description()));
        }

        TrapConfig cfg;

        // [geometry]
        if (auto geom = tbl["geometry"].as_table()) {
            cfg.geometry.r0 = geom->get("r0")->value_or(cfg.geometry.r0);
            cfg.geometry.z0 = geom->get("z0")->value_or(cfg.geometry.z0);
            cfg.geometry.electrode_gap = geom->get("electrode_gap")->value_or(cfg.geometry.electrode_gap);
        } else {
            throw std::runtime_error("ConfigParser: Missing [geometry] section in " + filepath);
        }

        // [mesh]
        if (auto mesh = tbl["mesh"].as_table()) {
            cfg.mesh.mesh_file = mesh->get("mesh_file")->value_or(cfg.mesh.mesh_file);
            cfg.mesh.cl_scale = mesh->get("cl_scale")->value_or(cfg.mesh.cl_scale);
            cfg.mesh.cl_min = mesh->get("cl_min")->value_or(cfg.mesh.cl_min);
            cfg.mesh.cl_max = mesh->get("cl_max")->value_or(cfg.mesh.cl_max);
            cfg.mesh.mesh_order = mesh->get("mesh_order")->value_or(cfg.mesh.mesh_order);
        } else {
            throw std::runtime_error("ConfigParser: Missing [mesh] section in " + filepath);
        }

        // [fem]
        if (auto fem = tbl["fem"].as_table()) {
            cfg.fem.order = fem->get("order")->value_or(cfg.fem.order);
            cfg.fem.refine_steps = fem->get("refine_steps")->value_or(cfg.fem.refine_steps);
            cfg.fem.solver_type = fem->get("solver_type")->value_or(cfg.fem.solver_type);
            cfg.fem.max_iter = fem->get("max_iter")->value_or(cfg.fem.max_iter);
            cfg.fem.rel_tol = fem->get("rel_tol")->value_or(cfg.fem.rel_tol);
            cfg.fem.abs_tol = fem->get("abs_tol")->value_or(cfg.fem.abs_tol);
        } else {
            throw std::runtime_error("ConfigParser: Missing [fem] section in " + filepath);
        }

        // [voltages]
        if (auto volt = tbl["voltages"].as_table()) {
            cfg.voltages.v1 = volt->get("v1")->value_or(cfg.voltages.v1);
            cfg.voltages.v2 = volt->get("v2")->value_or(cfg.voltages.v2);
            cfg.voltages.v3 = volt->get("v3")->value_or(cfg.voltages.v3);
            cfg.voltages.v4 = volt->get("v4")->value_or(cfg.voltages.v4);
            cfg.voltages.v5 = volt->get("v5")->value_or(cfg.voltages.v5);
            cfg.voltages.v6 = volt->get("v6")->value_or(cfg.voltages.v6);
        } else {
            throw std::runtime_error("ConfigParser: Missing [voltages] section in " + filepath);
        }

        // [interpolation]
        if (auto interp = tbl["interpolation"].as_table()) {
            cfg.interpolation.polynomial_degree = interp->get("polynomial_degree")->value_or(cfg.interpolation.polynomial_degree);
            cfg.interpolation.rbf_basis = interp->get("rbf_basis")->value_or(cfg.interpolation.rbf_basis);
            cfg.interpolation.rbf_smoothing = interp->get("rbf_smoothing")->value_or(cfg.interpolation.rbf_smoothing);
            cfg.interpolation.grid_x_min = interp->get("grid_x_min")->value_or(cfg.interpolation.grid_x_min);
            cfg.interpolation.grid_x_max = interp->get("grid_x_max")->value_or(cfg.interpolation.grid_x_max);
            cfg.interpolation.grid_y_min = interp->get("grid_y_min")->value_or(cfg.interpolation.grid_y_min);
            cfg.interpolation.grid_y_max = interp->get("grid_y_max")->value_or(cfg.interpolation.grid_y_max);
            cfg.interpolation.grid_z_min = interp->get("grid_z_min")->value_or(cfg.interpolation.grid_z_min);
            cfg.interpolation.grid_z_max = interp->get("grid_z_max")->value_or(cfg.interpolation.grid_z_max);
            cfg.interpolation.grid_nx = interp->get("grid_nx")->value_or(cfg.interpolation.grid_nx);
            cfg.interpolation.grid_ny = interp->get("grid_ny")->value_or(cfg.interpolation.grid_ny);
            cfg.interpolation.grid_nz = interp->get("grid_nz")->value_or(cfg.interpolation.grid_nz);
        } else {
            throw std::runtime_error("ConfigParser: Missing [interpolation] section in " + filepath);
        }

        // [output]
        if (auto out = tbl["output"].as_table()) {
            cfg.output.vtk_fem_file = out->get("vtk_fem_file")->value_or(cfg.output.vtk_fem_file);
            cfg.output.vtk_interpolated_file = out->get("vtk_interpolated_file")->value_or(cfg.output.vtk_interpolated_file);
            cfg.output.vtk_analytical_file = out->get("vtk_analytical_file")->value_or(cfg.output.vtk_analytical_file);
            cfg.output.vtk_residual_file = out->get("vtk_residual_file")->value_or(cfg.output.vtk_residual_file);
            cfg.output.report_polynomial_file = out->get("report_polynomial_file")->value_or(cfg.output.report_polynomial_file);
            cfg.output.report_metrics_file = out->get("report_metrics_file")->value_or(cfg.output.report_metrics_file);
        } else {
            throw std::runtime_error("ConfigParser: Missing [output] section in " + filepath);
        }

        validate(cfg);
        return cfg;
    }

    static void validate(const TrapConfig& cfg) {
        if (cfg.geometry.r0 <= 0.0) throw std::runtime_error("Invalid geometry: r0 must be positive.");
        if (cfg.geometry.z0 <= 0.0) throw std::runtime_error("Invalid geometry: z0 must be positive.");
        if (cfg.mesh.cl_scale <= 0.0) throw std::runtime_error("Invalid mesh: cl_scale must be positive.");
        if (cfg.fem.order < 1) throw std::runtime_error("Invalid fem: order must be >= 1.");
        if (cfg.interpolation.grid_nx < 2 || cfg.interpolation.grid_ny < 2 || cfg.interpolation.grid_nz < 2) {
            throw std::runtime_error("Invalid interpolation: grid resolutions must be >= 2.");
        }
        if (cfg.interpolation.grid_x_min >= cfg.interpolation.grid_x_max ||
            cfg.interpolation.grid_y_min >= cfg.interpolation.grid_y_max ||
            cfg.interpolation.grid_z_min >= cfg.interpolation.grid_z_max) {
            throw std::runtime_error("Invalid interpolation: grid min bounds must be strictly less than max bounds.");
        }
    }
};

} // namespace paultrap::core
