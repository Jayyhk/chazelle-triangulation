#pragma once

/// Intrusive doubly-linked polygon vertex list for Fournier-Montuno.
///
/// Algorithm 2 rewires prev/next pointers to split the polygon into
/// sub-polygons; Algorithm 3 removes vertices from this list.

#include "geometry/point.h"

#include <cstddef>
#include <limits>
#include <vector>

namespace chazelle {

struct VertexNode {
    std::size_t index = 0;   ///< Vertex index in the polygon.
    double x = 0.0;
    double y = 0.0;

    std::size_t prev = 0;    ///< Index into the VertexNode array.
    std::size_t next = 0;    ///< Index into the VertexNode array.

    /// Per FM Algorithm 1 output: "some vertices (of type 2) may
    /// point to two trapezoids."  We store up to 2 trapezoid indices.
    std::size_t trapezoid_idx  = std::numeric_limits<std::size_t>::max();
    std::size_t trapezoid_idx2 = std::numeric_limits<std::size_t>::max();
    bool done = false;       ///< Visited flag for Algorithm 2.
};

/// Build a circular doubly-linked list of VertexNodes from polygon vertices.
/// Returns the array; node i corresponds to vertex i.
/// prev/next form a circular chain: 0 → 1 → … → n−1 → 0.
inline std::vector<VertexNode>
build_vertex_list(const std::vector<Point>& vertices) {
    const std::size_t n = vertices.size();
    std::vector<VertexNode> nodes(n);
    for (std::size_t i = 0; i < n; ++i) {
        nodes[i].index = vertices[i].index;
        nodes[i].x     = vertices[i].x;
        nodes[i].y     = vertices[i].y;
        nodes[i].prev  = (i == 0) ? n - 1 : i - 1;
        nodes[i].next  = (i + 1 == n) ? 0 : i + 1;
    }
    return nodes;
}

} // namespace chazelle
