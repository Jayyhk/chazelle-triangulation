#pragma once

#include "../polygon/perturbation.h"
#include "../polygon/polygon.h"
#include "../submap/submap.h"
#include "oracle.h"

#include <cstddef>
#include <vector>

namespace chazelle::animation {

struct ArcSource {
    bool on_first_curve = true;
    std::size_t input_arc = NONE;
};

std::vector<ArcSource> identify_arc_sources(const Submap& submap, const Polygon& curve,
                                            const Submap& first_submap, const Polygon& first_curve,
                                            const Submap& second_submap,
                                            const Polygon& second_curve);

struct CycleArc {
    std::size_t arc = NONE;

    bool is_zero_length = false;
};

struct FusedRegionCycle {
    static constexpr std::size_t MAX_ARCS = 20;
    CycleArc arcs[MAX_ARCS];
    std::size_t count = 0;
};

FusedRegionCycle fused_region_cycle(const Submap& submap, const Polygon& curve, std::size_t region,
                                    const std::vector<std::size_t>& arcs_of_region);

struct FusedShootContext {
    const Submap* submap = nullptr;
    const Polygon* curve = nullptr;
    const Polygon* first_curve = nullptr;
    const Polygon* second_curve = nullptr;
    const RayShootingOracle* first_ray_shooter = nullptr;
    const RayShootingOracle* second_ray_shooter = nullptr;
    const std::vector<ArcSource>* arc_sources = nullptr;
};

RayHit local_shoot_fused(Point p, const SymbolicY& p_y, Side direction,
                         const FusedRegionCycle& cycle, const FusedShootContext& ctx,
                         bool require_hit = true, std::size_t source_edge_c = NONE);

struct VisiblePoint {
    bool found = false;

    std::size_t p_table_arc = NONE;
    std::size_t p_edge = NONE;
    Side p_side = LEFT;
    Exact p_x = 0.0;
    SymbolicY y{};

    std::size_t q_table_arc = NONE;
    std::size_t q_edge = NONE;
    Side q_side = LEFT;
    Exact q_x = 0.0;
};

struct ConformalityOracles {
    const Submap* first_submap = nullptr;
    const Submap* second_submap = nullptr;
    const Polygon* first_curve = nullptr;
    const Polygon* second_curve = nullptr;
    const RayShootingOracle* first_ray_shooter = nullptr;
    const RayShootingOracle* second_ray_shooter = nullptr;
    const ArcCuttingOracle* first_arc_cutter = nullptr;
    const ArcCuttingOracle* second_arc_cutter = nullptr;

    std::size_t first_piece_count_bound = 0;
    std::size_t second_piece_count_bound = 0;
    std::size_t first_piece_granularity_bound = 0;
    std::size_t second_piece_granularity_bound = 0;
};

VisiblePoint find_visible_point(const Submap& submap, const Polygon& curve, std::size_t region,
                                std::size_t A1, std::size_t A2, const FusedRegionCycle& cycle,
                                const std::vector<ArcSource>& arc_sources,
                                const ConformalityOracles& oracles);

void restore_conformality(Submap& submap, const Polygon& curve, const ConformalityOracles& oracles);

void restore_conformality(Submap& submap, const Polygon& curve,
                          const RayShootingOracle& ray_shooter, const ArcCuttingOracle& arc_cutter,
                          std::size_t piece_count_bound, std::size_t piece_granularity_bound);

}
