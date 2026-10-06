#pragma once

#include "../polygon/polygon.h"
#include "../submap/submap.h"

#include <cstddef>

namespace chazelle::animation {

struct RegionCompletionWork {
    std::size_t vertex_occurrences = 0;
    std::size_t height_comparisons = 0;
    std::size_t ray_edge_tests = 0;
    std::size_t endpoint_visits = 0;
};

struct RegionCompletion {
    Submap visibility_map;
    RegionCompletionWork work;
};

inline constexpr std::size_t REGION_COMPLETION_GRANULARITY = 4;

RegionCompletion complete_bounded_regions(const Polygon& curve, const Submap& submap,
                                          std::size_t granularity);

}
