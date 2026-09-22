#pragma once

#include "core/config_parser.hpp"
#include <string>
#include <vector>

namespace paultrap::geometry {

struct MeshStats {
    size_t num_vertices{0};
    size_t num_elements{0};
    size_t num_surfaces{0};
    std::string mesh_path;
    bool is_surface_only{false};
};

class Mesher {
public:
    explicit Mesher(const core::TrapConfig& config);
    ~Mesher() = default;

    // Generates 3D cylindrical Paul trap mesh (or 2D surface-only if surface_only=true)
    MeshStats generate(bool surface_only = false);

    // Direct generation with custom parameters
    static MeshStats generate(const core::GeometryConfig& geom, 
                              const core::MeshConfig& mesh_cfg,
                              bool surface_only = false);

    [[nodiscard]] const MeshStats& get_stats() const noexcept { return stats_; }

private:
    core::GeometryConfig geom_;
    core::MeshConfig mesh_cfg_;
    MeshStats stats_;
};

} // namespace paultrap::geometry
