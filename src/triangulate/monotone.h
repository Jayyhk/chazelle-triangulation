#pragma once

/// Algorithm 3: triangulate a unimonotone polygon.
///
/// A unimonotone polygon has one edge connecting the y-max and y-min
/// vertices, and one monotone chain between them.  The algorithm walks
/// the chain, removing convex vertices and emitting triangles.
///
/// Fournier & Montuno 1984, Algorithm 3.

#include "vertex_node.h"

#include <array>
#include <cstddef>
#include <vector>

namespace chazelle {

/// A triangle: three vertex indices.
using Triangle = std::array<std::size_t, 3>;

/// Triangulate a unimonotone sub-polygon delimited by nodes
/// `first` and `last` in the VertexNode list.
///
/// `first` and `last` are indices into `nodes`.  The sub-polygon
/// is the circular chain first → next → … → last → first.
///
/// Emits triangles into `out`.
void triangulate_monotone(std::vector<VertexNode>& nodes,
                          std::size_t first,
                          std::size_t last,
                          std::vector<Triangle>& out);

} // namespace chazelle
