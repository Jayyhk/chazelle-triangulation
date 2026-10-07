#pragma once

#include "../merge/oracle.h"
#include "../polygon/polygon.h"
#include "../submap/submap.h"

#include <cstddef>

namespace chazelle::animation {

RayHit naive_first_contact(const Polygon& curve, const Point& p, const SymbolicY& sy, Side dir,
                           std::size_t source_edge = NONE, std::size_t trace_query = NONE);

Submap build_full_visibility_map(const Polygon& curve);

std::size_t canonical_granularity(std::size_t num_edges);

Submap build_canonical_submap_naive(const Polygon& curve);

}
