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
};

class Mesher {
public:
    explicit Mesher(const core::TrapConfig& config);
    ~Mesher() = default;

    // Generates 3D cylindrical Paul trap mesh with 4 radial and 2 endcap electrodes
    MeshStats generate();

    // Direct generation with custom parameters
    static MeshStats generate(const core::GeometryConfig& geom, 
                              const core::MeshConfig& mesh_cfg);

    [[nodiscard]] const MeshStats& get_stats() const noexcept { return stats_; }

private:
    core::GeometryConfig geom_;
    core::MeshConfig mesh_cfg_;
    MeshStats stats_;
};

} // namespace paultrap::geometry
