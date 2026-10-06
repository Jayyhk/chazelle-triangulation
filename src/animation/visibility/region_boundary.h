#pragma once

#include "up_phase.h"

#include <array>
#include <memory>
#include <vector>

namespace chazelle::animation {

struct OriginalBoundaryLocation {
    std::size_t edge = NONE;
    Side side = LEFT;
    std::size_t arc = NONE;
};

struct RegionBoundaryPiece {
    Polygon curve;
    std::size_t first_original_vertex = NONE;
    bool reversed = false;
    Side original_side = LEFT;
    std::size_t original_arc = NONE;
    std::vector<OriginalBoundaryLocation> edge_locations;
    std::vector<SymbolicY> vertex_levels;
    std::shared_ptr<const Submap> canonical_submap{};
};

struct RegionBoundary {
    explicit RegionBoundary(std::vector<RegionBoundaryPiece> boundary_pieces,
                            std::size_t original_offset = 0);

    std::vector<RegionBoundaryPiece> pieces;
    Polygon curve;
    std::size_t original_offset;

    OriginalBoundaryLocation original_location(std::size_t edge, Side side) const;
    SymbolicY original_level(std::size_t tag, const Polygon& original) const;
};

struct CanonicalBoundaryChain {
    Polygon curve;
    Submap submap;
    std::unique_ptr<RayShootingStructure> rays;
};

class RegionBoundaryGeometry {
public:
    RegionBoundaryGeometry(const UpPhase& up_phase, const Polygon& curve);

    RegionBoundary boundary(const Polygon& curve, const Submap& submap, std::size_t region);

    UpPhase::PortionResult canonical_piece(const RegionBoundaryPiece& piece) const;

    UpPhase::PortionResult canonical_original_piece(std::size_t first, std::size_t last, Side side,
                                                    bool reversed) const;

    CanonicalBoundaryChain& canonical_chain(std::size_t grade, std::size_t index, Side side,
                                            bool reversed) const;

private:
    const UpPhase* up_phase_;
    std::size_t first_original_vertex_;
    std::size_t original_edge_count_;
    std::array<std::unique_ptr<Polygon>, 2> sides_;
    std::size_t next_tag_;
    mutable std::array<std::vector<std::vector<std::unique_ptr<CanonicalBoundaryChain>>>, 4>
        chains_;

    Point arc_endpoint(const Polygon& original, std::size_t global_edge, Side side,
                       const SymbolicY& level, bool first);
};

}
