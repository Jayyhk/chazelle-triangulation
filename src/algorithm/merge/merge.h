#pragma once

#include "../polygon/polygon.h"
#include "../submap/submap.h"
#include "oracle.h"

#include <cstddef>
#include <utility>

namespace chazelle {

struct MergeInput {
    const Polygon* first_curve = nullptr;
    const Polygon* second_curve = nullptr;

    const Submap* first_submap = nullptr;
    const Submap* second_submap = nullptr;

    std::size_t first_granularity = 0;
    std::size_t second_granularity = 0;
    std::size_t granularity = 0;

    const RayShootingOracle* first_ray_shooter = nullptr;
    const RayShootingOracle* second_ray_shooter = nullptr;
    const ArcCuttingOracle* first_arc_cutter = nullptr;
    const ArcCuttingOracle* second_arc_cutter = nullptr;

    std::size_t first_piece_count_bound = 0;
    std::size_t second_piece_count_bound = 0;
    std::size_t first_piece_granularity_bound = 0;
    std::size_t second_piece_granularity_bound = 0;
};

struct MergeResult {
    Polygon curve;
    Submap submap;
    explicit MergeResult(Polygon merged_curve) : curve(std::move(merged_curve)) {}
};

void assert_merge_preconditions(const MergeInput& input);

MergeResult merge(const MergeInput& input);

}
