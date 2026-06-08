#pragma once

/// [C91 §2.1]: Directed edge of a simple polygonal curve.

#include <cstddef>

namespace chazelle {

/// A directed edge from vertex start_idx to vertex end_idx.
struct Edge {
    std::size_t start_idx; ///< Index of the start vertex.
    std::size_t end_idx;   ///< Index of the end vertex.
};

} // namespace chazelle
