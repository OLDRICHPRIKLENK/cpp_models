cmake_minimum_required(VERSION 3.20)
project(paul_trap_sim LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add Homebrew prefix path for macOS
if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew")
endif()

# Compiler warnings and optimization
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-Wall -Wextra -O3)
endif()

# Suppress deprecated warnings from third-party headers
add_compile_options(-Wno-deprecated-declarations -Wno-deprecated-literal-operator)

# -----------------------------------------------------------------------------
# Find Required Packages
# -----------------------------------------------------------------------------

# Eigen3
find_package(Eigen3 REQUIRED)

# tomlplusplus
find_package(tomlplusplus REQUIRED)

# nanoflann
find_package(nanoflann REQUIRED)

# MFEM
find_package(mfem REQUIRED)

# VTK
find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel IOXML)

# Gmsh
find_package(gmsh QUIET)
if (NOT TARGET gmsh AND NOT TARGET gmsh::gmsh)
    find_path(GMSH_INCLUDE_DIR NAMES gmsh.h HINTS /opt/homebrew/include /usr/local/include)
    find_library(GMSH_LIBRARY NAMES gmsh HINTS /opt/homebrew/lib /usr/local/lib)
    if (GMSH_INCLUDE_DIR AND GMSH_LIBRARY)
        add_library(gmsh_imported SHARED IMPORTED)
        set_target_properties(gmsh_imported PROPERTIES
            IMPORTED_LOCATION "${GMSH_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMSH_INCLUDE_DIR}"
        )
        add_library(gmsh::gmsh ALIAS gmsh_imported)
    else()
        message(FATAL_ERROR "Gmsh headers or library not found!")
    endif()
endif()

# ALGLIB (via FetchContent)
include(FetchContent)
FetchContent_Declare(
    alglib
    GIT_REPOSITORY https://github.com/hjkuijf/ALGLIB.git
    GIT_TAG master
)
FetchContent_MakeAvailable(alglib)
if(NOT TARGET alglib)
    file(GLOB ALGLIB_SOURCES "${alglib_SOURCE_DIR}/src/*.cpp")
    add_library(alglib STATIC ${ALGLIB_SOURCES})
    target_include_directories(alglib PUBLIC "${alglib_SOURCE_DIR}/src")
    target_compile_options(alglib PRIVATE -Wno-deprecated-declarations)
endif()

# -----------------------------------------------------------------------------
# Core Library: paul_trap_lib
# -----------------------------------------------------------------------------

set(LIB_SOURCES
    src/geometry/mesher.cpp
    src/physics/analytical.cpp
    src/physics/fem_solver.cpp
    src/interpolation/regular_grid.cpp
    src/interpolation/rbf_interpolator.cpp
    src/validation/metrics.cpp
)

add_library(paul_trap_lib STATIC ${LIB_SOURCES})

target_include_directories(paul_trap_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(paul_trap_lib PUBLIC
    Eigen3::Eigen
    tomlplusplus::tomlplusplus
    nanoflann::nanoflann
    mfem
    gmsh::gmsh
    alglib
    ${VTK_LIBRARIES}
)

# -----------------------------------------------------------------------------
# Main Application Executable
# -----------------------------------------------------------------------------

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
    add_executable(paul_trap_sim src/main.cpp)
    target_link_libraries(paul_trap_sim PRIVATE paul_trap_lib)
endif()

# -----------------------------------------------------------------------------
# Testing Suite
# -----------------------------------------------------------------------------

enable_testing()

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/config.toml" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

function(add_paultrap_test TEST_NAME TEST_SOURCE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    target_link_libraries(${TEST_NAME} PRIVATE paul_trap_lib)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    set_tests_properties(${TEST_NAME} PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
endfunction()

add_paultrap_test(test_config tests/test_config.cpp)
add_paultrap_test(test_mesher tests/test_mesher.cpp)
add_paultrap_test(test_analytical tests/test_analytical.cpp)
add_paultrap_test(test_fem tests/test_fem.cpp)
add_paultrap_test(test_interpolation tests/test_interpolation.cpp)
add_paultrap_test(test_validation tests/test_validation.cpp)
