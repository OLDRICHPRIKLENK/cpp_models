#pragma once

#include <numbers>

namespace paultrap::core {

// Mathematical constants
inline constexpr double PI = std::numbers::pi_v<double>;
inline constexpr double TWO_PI = 2.0 * PI;
inline constexpr double HALF_PI = PI / 2.0;

// Fundamental physical constants (SI units)
inline constexpr double EPSILON_0 = 8.8541878128e-12;       // Vacuum permittivity [F/m]
inline constexpr double ELEMENTARY_CHARGE = 1.602176634e-19; // Elementary charge [C]
inline constexpr double ATOMIC_MASS_UNIT = 1.66053906660e-27; // 1 amu / 1 Da [kg]
inline constexpr double SPEED_OF_LIGHT = 299792458.0;       // Speed of light [m/s]
inline constexpr double BOLTZMANN_CONSTANT = 1.380649e-23;  // Boltzmann constant [J/K]

// Trap system constants
inline constexpr int NUM_ELECTRODES = 6;
inline constexpr int NUM_RADIAL_ELECTRODES = 4;
inline constexpr int NUM_ENDCAP_ELECTRODES = 2;

// Physical Group Boundary IDs for Gmsh and MFEM
inline constexpr int BOUNDARY_ID_RADIAL_1 = 1;
inline constexpr int BOUNDARY_ID_RADIAL_2 = 2;
inline constexpr int BOUNDARY_ID_RADIAL_3 = 3;
inline constexpr int BOUNDARY_ID_RADIAL_4 = 4;
inline constexpr int BOUNDARY_ID_ENDCAP_1 = 5;
inline constexpr int BOUNDARY_ID_ENDCAP_2 = 6;
inline constexpr int DOMAIN_VOLUME_ID     = 1;

} // namespace paultrap::core
