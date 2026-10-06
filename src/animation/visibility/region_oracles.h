#pragma once

#include "region_boundary.h"

#include <vector>

namespace chazelle::animation {

class RegionBoundaryOracles final : public RayShootingOracle, public ArcCuttingOracle {
public:
    RegionBoundaryOracles(const RegionBoundaryGeometry& geometry, const RegionBoundary& boundary,
                          const Polygon& input_curve, std::size_t input_offset, std::size_t grade);

    RayHit shoot(Point origin, Side direction, std::size_t arc, const Subarc& target,
                 SourceOffset offset = SOURCE_OFFSET_NONE) const override;

    std::vector<ArcPiece> cut(std::size_t arc, const Subarc& target) const override;

    static std::size_t piece_count_bound(std::size_t grade) {
        return 192 * (2 + 6 * grade);
    }

private:
    struct Piece {
        Subarc subarc;
        CanonicalBoundaryChain* chain = nullptr;
        std::size_t granularity = 0;
    };
    struct SingleEdge {
        std::size_t edge;
        std::unique_ptr<CanonicalBoundaryChain> canonical;
    };

    const RegionBoundaryGeometry* geometry_;
    const RegionBoundary* boundary_;
    const Polygon* input_curve_;
    std::size_t input_offset_;
    std::size_t grade_;
    mutable std::vector<SingleEdge> edges_;

    std::vector<Piece> partition(const Subarc& target) const;
    CanonicalBoundaryChain& single_edge(std::size_t edge) const;
};

UpPhase::PortionResult canonical_region_boundary(const RegionBoundaryGeometry& geometry,
                                                 const RegionBoundary& boundary, std::size_t grade);

}
