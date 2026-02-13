#pragma once

/// Bridging step: convert V(P) (a normal-form Submap with all chords)
/// into the Fournier-Montuno trapezoid + vertex-list format.
///
/// Each region of V(P) corresponds to one trapezoid.  The trapezoid's
/// top/bottom vertices are determined by the bounding chords' y-
/// coordinates, and its left/right edges come from the polygon edges
/// forming the region's arcs.
///
/// This is a linear-time traversal of the submap tree.

#include "trapezoid.h"
#include "vertex_node.h"
#include "geometry/polygon.h"
#include "visibility/submap.h"

#include <cstddef>
#include <vector>

namespace chazelle {

struct ConvertResult {
    std::vector<VertexNode>  nodes;
    std::vector<Trapezoid>   traps;
    std::size_t              first = 0;
    std::size_t              last  = 0;
};

/// Convert a fully-refined visibility map (all chords retained) into
/// the Fournier-Montuno format.
///
/// @param vp       The complete visibility map (Submap with all chords).
/// @param polygon  The input polygon.
/// @return         ConvertResult with vertex-node list and trapezoids.
ConvertResult convert_submap_to_fm(const Submap& vp,
                                   const Polygon& polygon);

} // namespace chazelle
