#include "physics/fem_solver.hpp"
#include "core/constants.hpp"

#include <mfem.hpp>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkCellArray.h>

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace paultrap::physics {

struct FemSolver::Impl {
    core::TrapConfig config;
    std::unique_ptr<mfem::Mesh> mesh;
    std::unique_ptr<mfem::H1_FECollection> fec;
    std::unique_ptr<mfem::FiniteElementSpace> fespace;
    std::vector<std::unique_ptr<mfem::GridFunction>> basis_fields; // 6 fields
    std::unique_ptr<mfem::GridFunction> combined_field;
    bool solved{false};

    explicit Impl(const core::TrapConfig& cfg) : config(cfg) {
        basis_fields.resize(core::NUM_ELECTRODES);
    }
};

FemSolver::FemSolver(const core::TrapConfig& config)
    : pimpl_(std::make_unique<Impl>(config)) {}

FemSolver::~FemSolver() = default;

FemSolver::FemSolver(FemSolver&&) noexcept = default;
FemSolver& FemSolver::operator=(FemSolver&&) noexcept = default;

void FemSolver::solve_basis_fields() {
    const std::string& mesh_path = pimpl_->config.mesh.mesh_file;
    if (!std::filesystem::exists(mesh_path)) {
        throw std::runtime_error("FemSolver: Mesh file does not exist: " + mesh_path);
    }

    // 1. Load mesh from Gmsh MSH 2.2
    pimpl_->mesh = std::make_unique<mfem::Mesh>(mesh_path.c_str(), 1, 1);
    if (pimpl_->mesh->Dimension() != 3) {
        throw std::runtime_error("FemSolver: Expected 3D mesh, found dimension " + 
                                 std::to_string(pimpl_->mesh->Dimension()));
    }

    // Uniform mesh refinements if configured
    for (int lev = 0; lev < pimpl_->config.fem.refine_steps; ++lev) {
        pimpl_->mesh->UniformRefinement();
    }

    // 2. Build H1 finite element space
    pimpl_->fec = std::make_unique<mfem::H1_FECollection>(pimpl_->config.fem.order, 3);
    pimpl_->fespace = std::make_unique<mfem::FiniteElementSpace>(pimpl_->mesh.get(), pimpl_->fec.get());

    std::cout << "[FemSolver] Mesh loaded: " << pimpl_->mesh->GetNV() << " vertices, "
              << pimpl_->mesh->GetNE() << " elements, "
              << pimpl_->fespace->GetNDofs() << " H1 DOFs." << std::endl;

    // 3. Assemble Laplace diffusion operator (div grad u = 0)
    auto a = std::make_unique<mfem::BilinearForm>(pimpl_->fespace.get());
    mfem::ConstantCoefficient diffusion_coef(1.0);
    a->AddDomainIntegrator(new mfem::DiffusionIntegrator(diffusion_coef));
    a->Assemble();

    // 4. Mark all boundary attributes as essential Dirichlet boundaries
    const int max_bdr = pimpl_->mesh->bdr_attributes.Max();
    mfem::Array<int> ess_bdr(max_bdr);
    ess_bdr = 1;
    mfem::Array<int> ess_tdof_list;
    pimpl_->fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // 5. Solve for each of the 6 basis fields
    for (int k = 0; k < core::NUM_ELECTRODES; ++k) {
        pimpl_->basis_fields[k] = std::make_unique<mfem::GridFunction>(pimpl_->fespace.get());
        *pimpl_->basis_fields[k] = 0.0;

        // Apply 1V to electrode k+1, 0V to all others
        mfem::Vector bdr_values(max_bdr);
        bdr_values = 0.0;
        if (k < max_bdr) {
            bdr_values(k) = 1.0;
        }
        mfem::PWConstCoefficient bdr_coeff(bdr_values);
        pimpl_->basis_fields[k]->ProjectBdrCoefficient(bdr_coeff, ess_bdr);

        // Zero RHS for Laplace equation
        mfem::LinearForm b(pimpl_->fespace.get());
        b = 0.0;

        mfem::OperatorPtr A;
        mfem::Vector B, X;
        a->FormLinearSystem(ess_tdof_list, *pimpl_->basis_fields[k], b, A, X, B);

        // Solve linear system with Conjugate Gradient
        mfem::CGSolver solver;
        solver.SetRelTol(pimpl_->config.fem.rel_tol);
        solver.SetAbsTol(pimpl_->config.fem.abs_tol);
        solver.SetMaxIter(pimpl_->config.fem.max_iter);
        solver.SetPrintLevel(0);
        solver.SetOperator(*A);
        solver.Mult(B, X);

        a->RecoverFEMSolution(X, b, *pimpl_->basis_fields[k]);
    }

    pimpl_->solved = true;
    std::cout << "[FemSolver] Successfully solved all " << core::NUM_ELECTRODES 
              << " basis Laplace fields." << std::endl;

    // Compute initial combined field using configuration voltages
    compute_combined_field(pimpl_->config.voltages.to_vector());
}

void FemSolver::compute_combined_field(const Eigen::Matrix<double, 6, 1>& voltages) {
    if (!pimpl_->solved) {
        throw std::runtime_error("FemSolver: Cannot compute combined field before solving basis fields.");
    }

    if (!pimpl_->combined_field) {
        pimpl_->combined_field = std::make_unique<mfem::GridFunction>(pimpl_->fespace.get());
    }
    *pimpl_->combined_field = 0.0;

    for (int k = 0; k < core::NUM_ELECTRODES; ++k) {
        pimpl_->combined_field->Add(voltages(k), *pimpl_->basis_fields[k]);
    }
}

void FemSolver::compute_combined_field(const core::VoltagesConfig& voltages) {
    compute_combined_field(voltages.to_vector());
}

NodalFieldData FemSolver::extract_nodal_data() const {
    if (!pimpl_->solved || !pimpl_->combined_field) {
        throw std::runtime_error("FemSolver: extract_nodal_data called before solving.");
    }

    const int num_vertices = pimpl_->mesh->GetNV();
    NodalFieldData data;
    data.coordinates.resize(3, num_vertices);
    data.potential.resize(num_vertices);
    data.basis_potentials.resize(core::NUM_ELECTRODES, num_vertices);

    for (int v = 0; v < num_vertices; ++v) {
        const double* coord = pimpl_->mesh->GetVertex(v);
        data.coordinates(0, v) = coord[0];
        data.coordinates(1, v) = coord[1];
        data.coordinates(2, v) = coord[2];

        mfem::Array<int> dofs;
        pimpl_->fespace->GetVertexDofs(v, dofs);
        const int dof = dofs[0];
        data.potential(v) = (*pimpl_->combined_field)(dof);

        for (int k = 0; k < core::NUM_ELECTRODES; ++k) {
            data.basis_potentials(k, v) = (*pimpl_->basis_fields[k])(dof);
        }
    }

    return data;
}

void FemSolver::export_vtk(const std::string& filepath) const {
    if (!pimpl_->solved || !pimpl_->combined_field) {
        throw std::runtime_error("FemSolver: export_vtk called before solving.");
    }

    std::filesystem::path out_path(filepath);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    const int num_vertices = pimpl_->mesh->GetNV();
    const int num_elements = pimpl_->mesh->GetNE();

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(num_vertices);

    for (int v = 0; v < num_vertices; ++v) {
        const double* coord = pimpl_->mesh->GetVertex(v);
        points->SetPoint(v, coord[0], coord[1], coord[2]);
    }

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    for (int e = 0; e < num_elements; ++e) {
        mfem::Array<int> v;
        pimpl_->mesh->GetElementVertices(e, v);
        if (v.Size() == 4) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            for (int i = 0; i < 4; ++i) {
                tetra->GetPointIds()->SetId(i, v[i]);
            }
            cells->InsertNextCell(tetra);
        }
    }

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    grid->SetCells(VTK_TETRA, cells);

    // Potential scalar array
    auto pot_array = vtkSmartPointer<vtkDoubleArray>::New();
    pot_array->SetName("Potential");
    pot_array->SetNumberOfValues(num_vertices);

    for (int v = 0; v < num_vertices; ++v) {
        mfem::Array<int> dofs;
        pimpl_->fespace->GetVertexDofs(v, dofs);
        const int dof = dofs[0];
        pot_array->SetValue(v, (*pimpl_->combined_field)(dof));
    }

    grid->GetPointData()->AddArray(pot_array);
    grid->GetPointData()->SetActiveScalars("Potential");

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filepath.c_str());
    writer->SetInputData(grid);
    if (writer->Write() != 1) {
        throw std::runtime_error("FemSolver: Failed to write VTK file: " + filepath);
    }

    std::cout << "[FemSolver] Exported FEM solution to " << filepath << std::endl;
}

const mfem::Mesh& FemSolver::get_mesh() const {
    if (!pimpl_->mesh) throw std::runtime_error("FemSolver: Mesh not loaded.");
    return *pimpl_->mesh;
}

const mfem::GridFunction& FemSolver::get_basis_field(int index) const {
    if (index < 0 || index >= core::NUM_ELECTRODES) {
        throw std::out_of_range("FemSolver: Basis index out of range [0, 5].");
    }
    if (!pimpl_->basis_fields[index]) {
        throw std::runtime_error("FemSolver: Basis fields not solved yet.");
    }
    return *pimpl_->basis_fields[index];
}

const mfem::GridFunction& FemSolver::get_combined_field() const {
    if (!pimpl_->combined_field) {
        throw std::runtime_error("FemSolver: Combined field not computed yet.");
    }
    return *pimpl_->combined_field;
}

bool FemSolver::is_solved() const noexcept {
    return pimpl_ && pimpl_->solved;
}

} // namespace paultrap::physics
