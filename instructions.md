Markdown
# C++ 3D Paul Trap Field Solver: Architecture & Implementation Guide

## Objective
Act as an expert computational physicist and C++ software architect. Your task is to build a highly optimized, native C++ codebase for precise 3D quadrupole field calculations. 

The architecture must adhere strictly to the Single Responsibility Principle, use atomic composable methods, and be optimized for raw speed using modern C++20 features and dense matrix operations.

## Core Physics & Mathematical Approach
*   **Geometry:** 3D Paul trap with 4 radial electrodes and 2 endcap electrodes. Cylindrical surfaces defined by `r0` and `z0`.
*   **Basis Fields:** The solver must compute 6 independent basis potentials by applying 1V to a single electrode while grounding the other 5. The final spatial potential is the dot product of these 6 basis fields and a 6-element time-varying voltage control vector.
*   **Analytical Ground Truth:** Include a module with the ideal analytical 3D quadrupole potential equations for V(x,y,z) to serve as a baseline.

## Required C++ Library Stack
Do not attempt to write from-scratch implementations for standard numerical tasks. Rely on the following industry-standard C++ libraries:
*   **Linear Algebra & Vectors:** `Eigen` (Dense operations, fast matrix math).
*   **Configuration:** `tomlplusplus` (Parsing runtime variables).
*   **Mesh Generation:** `Gmsh` (Native C++ API via `gmsh.h`).
*   **FEM Solver:** `MFEM` (Scalable C++ finite element library to solve Laplace's equation).
*   **Interpolation:** `ALGLIB` (for RBF) and `nanoflann` (KD-tree for fast regular-grid spatial lookups).
*   **Visualization:** Native `VTK` output (writing `.vtu` files for ParaView). Do not use Matplotlib wrappers.

## Directory Structure
Generate the files according to the following strict CMake-based directory structure:

```text
paul_trap_sim/
├── CMakeLists.txt              
├── vcpkg.json                  
├── config.toml                 
├── include/                    
│   ├── core/
│   │   ├── config_parser.hpp   
│   │   └── constants.hpp       
│   ├── geometry/
│   │   └── mesher.hpp          
│   ├── physics/
│   │   ├── analytical.hpp      
│   │   └── fem_solver.hpp      
│   ├── interpolation/
│   │   ├── regular_grid.hpp    
│   │   └── rbf_interpolator.hpp
│   └── validation/
│       └── metrics.hpp         
├── src/                        
│   ├── main.cpp                
│   ├── geometry/
│   │   └── mesher.cpp
│   ├── physics/
│   │   ├── analytical.cpp
│   │   └── fem_solver.cpp
│   ├── interpolation/
│   │   ├── regular_grid.cpp
│   │   └── rbf_interpolator.cpp
│   └── validation/
│       └── metrics.cpp
├── data/                       
│   ├── meshes/                 
│   └── fields/                 
└── output/                     
    ├── vtk/                    
    └── reports/ 
```


## Module Requirements
config.toml / config_parser: A centralized configuration setup containing all geometric parameters, mesh resolutions, interpolation bounds, and physical constants. Do not hardcode or duplicate configs in other files.

Mesher (geometry/): Interfaces with Gmsh to generate the 3D domain, assign physical groups to the 4 radial and 2 endcap electrodes, and export to data/meshes/.
FEM Solver (physics/): Uses MFEM to load the mesh and solve Laplace's equation 6 times (once for each basis field). Must expose methods to extract the solved nodal data.

Interpolator (interpolation/): Maps unstructured FEM data onto a regular 3D Cartesian grid. Implement both a fast simple polynomial interpolator and an RBF interpolator. Both must support spatial gradient evaluation (Electric field mapping).

Equation Logging: The module must print the resulting polynomial equation to the standard terminal output and save it as a text file in output/reports/.

Visualization (vtk): Generate .vtu files directly.
Strict SRP for Outputs: The Single Responsibility Principle applies to VTK files as well. Each .vtu file must represent exactly one object or dataset. Do not conflate multiple meshes or fields into a single combined file.

Generate separate .vtu files for: (a) The FEM mesh data, (b) The interpolated data grids, (c) The analytical ideal field, and (d) A dedicated residual error map (see Validation below).

Validation (validation/): A rigorous verification module.
Error Metrics: Must calculate the relative L2 error, absolute L2 error, and maximum (infinity norm) error comparing the simple polynomial interpolation results against the analytical ideal trap geometry equation across the central trapping region.

Residual Map: Generate a distinct dataset representing the spatial residual error map (interpolated vs. ideal) to be exported as its own standalone .vtu file.

## C++ Best Practices & Code Quality Requirements
To ensure the codebase is DRY, composable, and easily extensible, you must strictly adhere to the following C++ paradigms:

Always perform unittests. After any implementation, first test out the module functionality before moving on to the next task.

RAII & Memory Management: Never use raw new or delete. Use std::unique_ptr, std::shared_ptr, or standard library containers (std::vector) to manage lifetimes.

Const-Correctness & References: Avoid unnecessary memory copies. Pass heavy objects (like Eigen vectors, strings, and MFEM GridFunctions) by constant reference (e.g., const Eigen::VectorXd& vector). Mark class methods that do not modify state as const.

Compile-Time Optimization: Utilize constexpr for mathematical constants and fixed geometry bounds to offload work to the compiler.

Namespaces: Wrap all modules in a root namespace (e.g., namespace paultrap { ... }) and use sub-namespaces (e.g., paultrap::physics) to prevent global scope pollution.

Decoupled Interfaces: Keep header files (.hpp) clean. Put implementation details in the .cpp files. Avoid #include bloat in headers by using forward declarations where possible.

Error Handling: Do not fail silently. Throw std::runtime_error with descriptive messages if file I/O fails (e.g., missing mesh file) or if matrix assembly fails.