#pragma once

#include "../polygon/perturbation.h"
#include "../polygon/polygon.h"
#include "../submap/boundary_geometry.h"
#include "../submap/submap.h"
#include "oracle.h"

#include <cstddef>
#include <vector>

namespace chazelle::animation {

struct FusionVertex {
    SymbolicY y;
    std::size_t edge;
    Side side;
    std::size_t chord_idx;
    bool is_left_endpoint;
    bool is_companion;
};

RayHit local_shoot(const Point& p, Side direction, std::size_t region, const Submap& submap,
                   const Polygon& curve, const RayShootingOracle& oracle, bool require_hit = true,
                   const SourceOffset& source_x_offset = SOURCE_OFFSET_NONE, bool record = true);

struct FusionState {
    std::vector<FusionVertex> sequence;
    std::size_t current_stop = 0;

    bool junction_at_end = true;

    Point p{0.0, 0.0, NONE};
    std::size_t p_edge = NONE;
    Side p_side = LEFT;

    SymbolicY p_y{};

    std::size_t current_region = NONE;

    std::vector<std::size_t> arc_starts;

    struct DiscoveredChord {
        SymbolicY y;
        std::size_t left_edge;
        Side left_side;
        std::size_t right_edge;
        Side right_side;
        bool left_on_first_curve = true;
        bool right_on_first_curve = true;
    };
    std::vector<DiscoveredChord> chords;

    std::vector<bool> invalidated_first_chords;
    std::vector<bool> invalidated_second_chords;
};

std::size_t fusion_startup(FusionState& state, const Submap& first_submap,
                           const Polygon& first_curve, const Submap& second_submap,
                           const Polygon& second_curve, const RayShootingOracle& oracle1,
                           const RayShootingOracle& oracle2);

void fuse_submaps(FusionState& state, const Submap& first_submap, const Polygon& first_curve,
                  const Submap& second_submap, const Polygon& second_curve,
                  const RayShootingOracle& oracle1, const RayShootingOracle& oracle2);

void build_fusion_sequence(FusionState& state, const Submap& submap, const Polygon& curve);

void rebuild_submap(Submap& submap, const Polygon& curve, const Submap& first_submap,
                    const Polygon& first_curve, const Submap& second_submap,
                    const Polygon& second_curve, const FusionState& state1,
                    const FusionState& state2);
}
