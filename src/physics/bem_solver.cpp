#include "physics/bem_solver.hpp"
#include "core/constants.hpp"

#include <gmsh.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkTriangle.h>
#include <cmath>
#include <array>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <filesystem>
#include <stdexcept>

namespace paultrap::physics {

namespace {

struct GmshSessionGuard {
    GmshSessionGuard() {
        if (!gmsh::isInitialized()) {
            gmsh::initialize();
            gmsh::option::setNumber("General.Terminal", 0);
        }
    }
    ~GmshSessionGuard() {
        if (gmsh::isInitialized()) {
            gmsh::finalize();
        }
    }
};

struct SurfaceTriangle {
    std::array<size_t, 3> node_indices; // 0-based indices into mesh vertex list
    int physical_id{0};                  // 1..6 (Radial1..4, Endcap1..2)
    Eigen::Vector3d centroid;
    Eigen::Vector3d normal;              // Outward-pointing unit normal
    double area{0.0};
};

struct SurfaceMesh {
    Eigen::Matrix3Xd vertices;            // 3 x N_vertices
    std::vector<SurfaceTriangle> elements; // N_elements
};

// 4-point symmetric Gaussian quadrature rule on triangle (Dunavant / Hammer rule, degree 3)
struct TriangleQuadRule {
    static constexpr int NUM_POINTS = 4;
    std::array<double, 4> weights;
    std::array<std::array<double, 3>, 4> bary_coords;

    TriangleQuadRule() {
        weights[0] = -9.0 / 16.0;
        bary_coords[0] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};

        const double w_perm = 25.0 / 48.0;
        weights[1] = w_perm;
        bary_coords[1] = {3.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0};

        weights[2] = w_perm;
        bary_coords[2] = {1.0 / 5.0, 3.0 / 5.0, 1.0 / 5.0};

        weights[3] = w_perm;
        bary_coords[3] = {1.0 / 5.0, 1.0 / 5.0, 3.0 / 5.0};
    }
};

const TriangleQuadRule QUAD_RULE;

// Analytical closed-form evaluation of weakly singular integral:
// \int_{\Gamma_i} \frac{1}{4\pi \|\mathbf{x}_i - \mathbf{y}\|} d\Gamma(\mathbf{y})
// where \mathbf{x}_i is the centroid of flat triangle \Gamma_i.
double compute_singular_G_ii(const Eigen::Vector3d& v0, 
                             const Eigen::Vector3d& v1, 
                             const Eigen::Vector3d& v2, 
                             const Eigen::Vector3d& centroid) {
    const std::array<Eigen::Vector3d, 3> verts = {v0, v1, v2};
    double total_integral = 0.0;

    for (int k = 0; k < 3; ++k) {
        const Eigen::Vector3d a = verts[k] - centroid;
        const Eigen::Vector3d b = verts[(k + 1) % 3] - centroid;
        const Eigen::Vector3d edge = b - a;
        const double L = edge.norm();
        if (L < 1.0e-14) continue;

        const Eigen::Vector3d t = edge / L; // Unit tangent along edge
        const double sA = a.dot(t);
        const double sB = b.dot(t); // sB = sA + L

        const Eigen::Vector3d p_foot = a - sA * t;
        const double h = p_foot.norm(); // Perpendicular distance to edge line
        if (h < 1.0e-14) continue;

        const double rA = a.norm();
        const double rB = b.norm();

        const double num = rB + sB;
        const double den = rA + sA;
        if (den > 1.0e-14 && num > 1.0e-14) {
            total_integral += h * std::log(num / den);
        }
    }

    return total_integral / (4.0 * core::PI);
}

SurfaceMesh load_surface_mesh_from_gmsh(const std::string& mesh_file) {
    if (!std::filesystem::exists(mesh_file)) {
        throw std::runtime_error("BemSolver: Mesh file does not exist: " + mesh_file);
    }

    GmshSessionGuard guard;
    gmsh::clear();
    gmsh::open(mesh_file);

    std::vector<std::size_t> node_tags;
    std::vector<double> coords;
    std::vector<double> param_coords;
    gmsh::model::mesh::getNodes(node_tags, coords, param_coords);

    if (node_tags.empty()) {
        throw std::runtime_error("BemSolver: No nodes found in mesh file: " + mesh_file);
    }

    std::unordered_map<std::size_t, size_t> tag_to_idx;
    tag_to_idx.reserve(node_tags.size());

    SurfaceMesh mesh;
    mesh.vertices.resize(3, node_tags.size());
    for (size_t i = 0; i < node_tags.size(); ++i) {
        tag_to_idx[node_tags[i]] = i;
        mesh.vertices(0, i) = coords[3 * i + 0];
        mesh.vertices(1, i) = coords[3 * i + 1];
        mesh.vertices(2, i) = coords[3 * i + 2];
    }

    const std::array<int, 6> boundary_ids = {
        core::BOUNDARY_ID_RADIAL_1, core::BOUNDARY_ID_RADIAL_2,
        core::BOUNDARY_ID_RADIAL_3, core::BOUNDARY_ID_RADIAL_4,
        core::BOUNDARY_ID_ENDCAP_1, core::BOUNDARY_ID_ENDCAP_2
    };

    for (int bdr_id : boundary_ids) {
        std::vector<int> entity_tags;
        gmsh::model::getEntitiesForPhysicalGroup(2, bdr_id, entity_tags);

        for (int ent : entity_tags) {
            std::vector<int> elem_types;
            std::vector<std::vector<std::size_t>> elem_tags;
            std::vector<std::vector<std::size_t>> elem_node_tags;
            gmsh::model::mesh::getElements(elem_types, elem_tags, elem_node_tags, 2, ent);

            for (size_t t = 0; t < elem_types.size(); ++t) {
                if (elem_types[t] == 2) { // 3-node triangle
                    const auto& nodes = elem_node_tags[t];
                    const size_t num_tris = elem_tags[t].size();
                    for (size_t j = 0; j < num_tris; ++j) {
                        size_t idx0 = tag_to_idx.at(nodes[3 * j + 0]);
                        size_t idx1 = tag_to_idx.at(nodes[3 * j + 1]);
                        size_t idx2 = tag_to_idx.at(nodes[3 * j + 2]);

                        const Eigen::Vector3d v0 = mesh.vertices.col(idx0);
                        const Eigen::Vector3d v1 = mesh.vertices.col(idx1);
                        const Eigen::Vector3d v2 = mesh.vertices.col(idx2);

                        Eigen::Vector3d cross = (v1 - v0).cross(v2 - v0);
                        double area = 0.5 * cross.norm();
                        if (area < 1.0e-15) continue;

                        Eigen::Vector3d normal = cross.normalized();
                        Eigen::Vector3d centroid = (v0 + v1 + v2) / 3.0;

                        // Ensure outward-pointing normal based on electrode ID
                        if (bdr_id >= core::BOUNDARY_ID_RADIAL_1 && bdr_id <= core::BOUNDARY_ID_RADIAL_4) {
                            // Cylindrical radial surface: normal must point away from z-axis
                            if (normal.dot(Eigen::Vector3d(centroid.x(), centroid.y(), 0.0)) < 0.0) {
                                normal = -normal;
                                std::swap(idx1, idx2);
                            }
                        } else if (bdr_id == core::BOUNDARY_ID_ENDCAP_1) {
                            // Top endcap (+z): normal must point +z
                            if (normal.z() < 0.0) {
                                normal = -normal;
                                std::swap(idx1, idx2);
                            }
                        } else if (bdr_id == core::BOUNDARY_ID_ENDCAP_2) {
                            // Bottom endcap (-z): normal must point -z
                            if (normal.z() > 0.0) {
                                normal = -normal;
                                std::swap(idx1, idx2);
                            }
                        }

                        SurfaceTriangle tri;
                        tri.node_indices = {idx0, idx1, idx2};
                        tri.physical_id = bdr_id;
                        tri.centroid = centroid;
                        tri.normal = normal;
                        tri.area = area;

                        mesh.elements.push_back(tri);
                    }
                }
            }
        }
    }

    if (mesh.elements.empty()) {
        throw std::runtime_error("BemSolver: No 2D triangle boundary elements found in mesh: " + mesh_file);
    }

    return mesh;
}

} // anonymous namespace

struct BemSolver::Impl {
    core::TrapConfig config;
    SurfaceMesh mesh;
    bool solved{false};

    // Dense BEM matrices
    Eigen::MatrixXd G;
    Eigen::MatrixXd H;

    // Basis solutions: Q is N x 6 matrix of surface charge densities for basis problems 1..6
    Eigen::MatrixXd basis_charges;  // N x 6
    Eigen::MatrixXd basis_voltages; // N x 6

    // Combined solution
    Eigen::VectorXd surface_charges;  // N
    Eigen::VectorXd surface_voltages; // N

    explicit Impl(const core::TrapConfig& cfg) : config(cfg) {}
};

BemSolver::BemSolver(const core::TrapConfig& config)
    : pimpl_(std::make_unique<Impl>(config)) {}

BemSolver::~BemSolver() = default;

BemSolver::BemSolver(BemSolver&&) noexcept = default;
BemSolver& BemSolver::operator=(BemSolver&&) noexcept = default;

void BemSolver::solve_basis_fields() {
    const std::string& mesh_path = pimpl_->config.mesh.mesh_file;
    pimpl_->mesh = load_surface_mesh_from_gmsh(mesh_path);

    const size_t N = pimpl_->mesh.elements.size();
    std::cout << "[BemSolver] Loaded surface mesh with " << pimpl_->mesh.vertices.cols() 
              << " vertices and " << N << " boundary triangles." << std::endl;

    pimpl_->G.resize(N, N);
    pimpl_->H.resize(N, N);

    // Pre-calculate quadrature point coordinates on all triangles
    // Shape: N x 4 x 3
    std::vector<std::array<Eigen::Vector3d, TriangleQuadRule::NUM_POINTS>> quad_points(N);
    for (size_t j = 0; j < N; ++j) {
        const auto& tri = pimpl_->mesh.elements[j];
        const Eigen::Vector3d v0 = pimpl_->mesh.vertices.col(tri.node_indices[0]);
        const Eigen::Vector3d v1 = pimpl_->mesh.vertices.col(tri.node_indices[1]);
        const Eigen::Vector3d v2 = pimpl_->mesh.vertices.col(tri.node_indices[2]);

        for (int p = 0; p < TriangleQuadRule::NUM_POINTS; ++p) {
            const auto& b = QUAD_RULE.bary_coords[p];
            quad_points[j][p] = b[0] * v0 + b[1] * v1 + b[2] * v2;
        }
    }

    // Assemble dense matrices G and H
    std::cout << "[BemSolver] Assembling " << N << " x " << N << " collocation matrices H and G..." << std::endl;

    for (size_t i = 0; i < N; ++i) {
        const Eigen::Vector3d& x_i = pimpl_->mesh.elements[i].centroid;

        for (size_t j = 0; j < N; ++j) {
            if (i == j) {
                // Diagonal self-element
                const auto& tri = pimpl_->mesh.elements[i];
                const Eigen::Vector3d v0 = pimpl_->mesh.vertices.col(tri.node_indices[0]);
                const Eigen::Vector3d v1 = pimpl_->mesh.vertices.col(tri.node_indices[1]);
                const Eigen::Vector3d v2 = pimpl_->mesh.vertices.col(tri.node_indices[2]);

                pimpl_->G(i, i) = compute_singular_G_ii(v0, v1, v2, x_i);
                pimpl_->H(i, i) = 0.5; // Theoretical flat-element jump term
            } else {
                // Off-diagonal smooth element
                const auto& tri_j = pimpl_->mesh.elements[j];
                const double Aj = tri_j.area;
                const Eigen::Vector3d& nj = tri_j.normal;

                double g_sum = 0.0;
                double h_sum = 0.0;

                for (int p = 0; p < TriangleQuadRule::NUM_POINTS; ++p) {
                    const Eigen::Vector3d diff = x_i - quad_points[j][p];
                    const double r2 = diff.squaredNorm();
                    const double r = std::sqrt(r2);
                    const double inv_r = 1.0 / r;
                    const double inv_r3 = inv_r / r2;
                    const double w = QUAD_RULE.weights[p];

                    g_sum += w * inv_r;
                    h_sum += w * (diff.dot(nj) * inv_r3);
                }

                pimpl_->G(i, j) = (Aj / (4.0 * core::PI)) * g_sum;
                pimpl_->H(i, j) = (Aj / (4.0 * core::PI)) * h_sum;
            }
        }

        // Rigid-body condition regularization: H(i, i) = -sum_{j != i} H(i, j)
        // Guarantees exact equipotential conservation H * 1 = 0
        double off_diag_sum = 0.0;
        for (size_t j = 0; j < N; ++j) {
            if (i != j) {
                off_diag_sum += pimpl_->H(i, j);
            }
        }
        pimpl_->H(i, i) = -off_diag_sum;
    }

    // Set up Dirichlet basis potentials (N x 6)
    // Electrode k has 1.0 V on boundary k+1 and 0.0 V elsewhere
    pimpl_->basis_voltages = Eigen::MatrixXd::Zero(N, core::NUM_ELECTRODES);
    for (size_t j = 0; j < N; ++j) {
        int phys_id = pimpl_->mesh.elements[j].physical_id;
        int k = phys_id - 1;
        if (k >= 0 && k < core::NUM_ELECTRODES) {
            pimpl_->basis_voltages(j, k) = 1.0;
        }
    }

    // Right-hand side matrix B = H * U (N x 6)
    Eigen::MatrixXd B = pimpl_->H * pimpl_->basis_voltages;

    // Solve G * Q = B for unknown surface normal charges Q (N x 6)
    std::cout << "[BemSolver] Performing dense LU factorization of matrix G..." << std::endl;
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(pimpl_->G);
    if (lu.info() != Eigen::Success) {
        throw std::runtime_error("BemSolver: LU decomposition of BEM matrix G failed.");
    }

    pimpl_->basis_charges = lu.solve(B);
    pimpl_->solved = true;
    std::cout << "[BemSolver] Successfully solved 6 basis field surface distributions." << std::endl;
}

void BemSolver::compute_combined_field(const Eigen::Matrix<double, 6, 1>& voltages) {
    if (!pimpl_->solved) {
        throw std::runtime_error("BemSolver: solve_basis_fields() must be called before compute_combined_field().");
    }

    pimpl_->surface_voltages = pimpl_->basis_voltages * voltages;
    pimpl_->surface_charges = pimpl_->basis_charges * voltages;
}

void BemSolver::compute_combined_field(const core::VoltagesConfig& voltages) {
    compute_combined_field(voltages.to_vector());
}

NodalFieldData BemSolver::extract_nodal_data(const Eigen::Matrix3Xd& eval_points) const {
    if (!pimpl_->solved) {
        throw std::runtime_error("BemSolver: solve_basis_fields() must be called before extract_nodal_data().");
    }

    const size_t N = pimpl_->mesh.elements.size();
    const size_t M = eval_points.cols();

    NodalFieldData data;
    data.coordinates = eval_points;
    data.potential.resize(M);
    data.basis_potentials.resize(core::NUM_ELECTRODES, M);

    // Precompute triangle quad points
    std::vector<std::array<Eigen::Vector3d, TriangleQuadRule::NUM_POINTS>> quad_points(N);
    for (size_t j = 0; j < N; ++j) {
        const auto& tri = pimpl_->mesh.elements[j];
        const Eigen::Vector3d v0 = pimpl_->mesh.vertices.col(tri.node_indices[0]);
        const Eigen::Vector3d v1 = pimpl_->mesh.vertices.col(tri.node_indices[1]);
        const Eigen::Vector3d v2 = pimpl_->mesh.vertices.col(tri.node_indices[2]);

        for (int p = 0; p < TriangleQuadRule::NUM_POINTS; ++p) {
            const auto& b = QUAD_RULE.bary_coords[p];
            quad_points[j][p] = b[0] * v0 + b[1] * v1 + b[2] * v2;
        }
    }

    // Direct boundary representation evaluation at all internal points
    // phi(x) = sum_j [ q_j * G_j(x) - u_j * H_j(x) ]
    for (size_t m = 0; m < M; ++m) {
        const Eigen::Vector3d x = eval_points.col(m);

        Eigen::Matrix<double, core::NUM_ELECTRODES, 1> basis_phi = Eigen::Matrix<double, core::NUM_ELECTRODES, 1>::Zero();
        double combined_phi = 0.0;

        for (size_t j = 0; j < N; ++j) {
            const auto& tri_j = pimpl_->mesh.elements[j];
            const double Aj = tri_j.area;
            const Eigen::Vector3d& nj = tri_j.normal;

            double g_sum = 0.0;
            double h_sum = 0.0;

            for (int p = 0; p < TriangleQuadRule::NUM_POINTS; ++p) {
                const Eigen::Vector3d diff = x - quad_points[j][p];
                const double r2 = diff.squaredNorm();
                const double r = std::sqrt(r2);
                const double inv_r = 1.0 / r;
                const double inv_r3 = inv_r / r2;
                const double w = QUAD_RULE.weights[p];

                g_sum += w * inv_r;
                h_sum += w * (diff.dot(nj) * inv_r3);
            }

            const double gj = (Aj / (4.0 * core::PI)) * g_sum;
            const double hj = (Aj / (4.0 * core::PI)) * h_sum;

            // Combined potential
            combined_phi += pimpl_->surface_charges(j) * gj - pimpl_->surface_voltages(j) * hj;

            // 6 Basis potentials
            for (int k = 0; k < core::NUM_ELECTRODES; ++k) {
                basis_phi(k) += pimpl_->basis_charges(j, k) * gj - pimpl_->basis_voltages(j, k) * hj;
            }
        }

        data.potential(m) = combined_phi;
        for (int k = 0; k < core::NUM_ELECTRODES; ++k) {
            data.basis_potentials(k, m) = basis_phi(k);
        }
    }

    return data;
}

void BemSolver::export_vtk(const std::string& filepath) const {
    if (!pimpl_->solved) {
        throw std::runtime_error("BemSolver: solve_basis_fields() must be called before export_vtk().");
    }

    // Ensure output directory exists
    std::filesystem::path out_path(filepath);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    auto points = vtkSmartPointer<vtkPoints>::New();
    const size_t num_verts = pimpl_->mesh.vertices.cols();
    for (size_t i = 0; i < num_verts; ++i) {
        points->InsertNextPoint(pimpl_->mesh.vertices(0, i),
                                pimpl_->mesh.vertices(1, i),
                                pimpl_->mesh.vertices(2, i));
    }

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);

    auto pot_array = vtkSmartPointer<vtkDoubleArray>::New();
    pot_array->SetName("Potential");

    auto q_array = vtkSmartPointer<vtkDoubleArray>::New();
    q_array->SetName("SurfaceChargeDensity");

    auto bdr_array = vtkSmartPointer<vtkIntArray>::New();
    bdr_array->SetName("BoundaryID");

    std::array<vtkSmartPointer<vtkDoubleArray>, 6> basis_q_arrays;
    for (int k = 0; k < 6; ++k) {
        basis_q_arrays[k] = vtkSmartPointer<vtkDoubleArray>::New();
        basis_q_arrays[k]->SetName(("BasisCharge_" + std::to_string(k + 1)).c_str());
    }

    const size_t num_elements = pimpl_->mesh.elements.size();
    for (size_t j = 0; j < num_elements; ++j) {
        const auto& tri = pimpl_->mesh.elements[j];
        auto cell = vtkSmartPointer<vtkTriangle>::New();
        cell->GetPointIds()->SetId(0, tri.node_indices[0]);
        cell->GetPointIds()->SetId(1, tri.node_indices[1]);
        cell->GetPointIds()->SetId(2, tri.node_indices[2]);
        grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds());

        pot_array->InsertNextValue(pimpl_->surface_voltages(j));
        q_array->InsertNextValue(pimpl_->surface_charges(j));
        bdr_array->InsertNextValue(tri.physical_id);
        for (int k = 0; k < 6; ++k) {
            basis_q_arrays[k]->InsertNextValue(pimpl_->basis_charges(j, k));
        }
    }

    grid->GetCellData()->AddArray(bdr_array);
    grid->GetCellData()->AddArray(pot_array);
    grid->GetCellData()->AddArray(q_array);
    grid->GetCellData()->SetActiveScalars("Potential");
    for (int k = 0; k < 6; ++k) {
        grid->GetCellData()->AddArray(basis_q_arrays[k]);
    }

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filepath.c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToBinary();

    if (writer->Write() != 1) {
        throw std::runtime_error("BemSolver: Failed to write VTK surface solution to " + filepath);
    }

    std::cout << "[BemSolver] Exported surface BEM solution to: " << filepath << std::endl;
}

bool BemSolver::is_solved() const noexcept {
    return pimpl_ && pimpl_->solved;
}

size_t BemSolver::num_elements() const noexcept {
    return pimpl_ ? pimpl_->mesh.elements.size() : 0;
}

const Eigen::VectorXd& BemSolver::get_surface_charges() const {
    return pimpl_->surface_charges;
}

const Eigen::VectorXd& BemSolver::get_surface_voltages() const {
    return pimpl_->surface_voltages;
}

const Eigen::Matrix<double, 6, Eigen::Dynamic>& BemSolver::get_basis_charges() const {
    // basis_charges is stored as N x 6, so transposing to 6 x N
    static Eigen::Matrix<double, 6, Eigen::Dynamic> transposed;
    transposed = pimpl_->basis_charges.transpose();
    return transposed;
}

} // namespace paultrap::physics
