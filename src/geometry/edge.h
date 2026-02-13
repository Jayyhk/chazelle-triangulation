#pragma once

#include <cstddef>

namespace chazelle {

/// A directed edge of the polygon, stored in the read-only input table.
/// Arc-structures reference edges by index into this table.
struct Edge {
    std::size_t index;     ///< Position in the input table (edge i connects vertex i to vertex i+1).
    std::size_t start_idx; ///< Index of the start vertex in the polygon.
    std::size_t end_idx;   ///< Index of the end vertex in the polygon.
};

} // namespace chazelle
