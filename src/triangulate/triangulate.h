#pragma once

/// Fournier-Montuno Algorithm 2 + 3 wrapper:
///   trapezoidized polygon → triangulation.
///
/// Algorithm 2 classifies trapezoids (Class A / Class B), inserts
/// diagonals for Class B, and recurses.  When no more Class B
/// trapezoids remain, the sub-polygon is unimonotone and Algorithm 3
/// handles it.

#include "trapezoid.h"
#include "vertex_node.h"
#include "monotone.h"

#include <cstddef>
#include <vector>

namespace chazelle {

/// Full triangulation from a trapezoidized polygon.
///
/// @param nodes   Doubly-linked vertex list (prev/next form a circular chain).
/// @param traps   Trapezoid records; each vertex's `trapezoid_idx` field
///                points into this array (or FM_NONE if no trapezoid).
/// @param first   Starting node index for the (sub-)polygon.
/// @param last    Ending node index (polygon wraps: … → last → first).
/// @return        List of triangles (vertex index triples).
std::vector<Triangle>
fm_triangulate(std::vector<VertexNode>& nodes,
               std::vector<Trapezoid>& traps,
               std::size_t first,
               std::size_t last);

} // namespace chazelle
