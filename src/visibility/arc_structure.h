#pragma once

/// ArcStructure — represents an arc of ∂C within a single region of a submap.
///
/// An arc is a connected piece of ∂C (the double boundary of the input curve)
/// between two consecutive chord endpoints, belonging to a single region.
///
/// Per Chazelle §2.3:
/// - Single-edge arc (t = 1): one input-table edge pointer + one side flag.
/// - Multi-edge arc (t > 1): two input-table edge pointers (first and last
///   edges of P containing the arc's first and last sub-edges, in clockwise
///   order) + two side flags.
/// - Null-length arcs occur at local extrema.
///
/// The endpoints of the arc are NOT stored — chords handle that.

#include "common.h"
#include "geometry/point.h"   // for Side enum

#include <cstddef>

namespace chazelle {

struct ArcStructure {
    /// Index of the first (or only) edge of P containing this arc.
    std::size_t first_edge = NONE;
    /// Side of the double boundary for the first edge.
    Side first_side = Side::LEFT;

    /// Index of the last edge of P containing this arc (= first_edge for t = 1).
    std::size_t last_edge = NONE;
    /// Side of the double boundary for the last edge.
    Side last_side = Side::LEFT;

    /// Back-pointer to the tree node (region) this arc belongs to.
    std::size_t region_node = NONE;

    /// Number of non-null-length polygon edges composing this arc.
    /// Used for weight computation.
    std::size_t edge_count = 0;

    /// Is this a null-length (zero-length) arc?
    bool is_null() const { return first_edge == last_edge && edge_count == 0; }

    /// Is this a single-edge arc?
    bool is_single_edge() const { return first_edge == last_edge && edge_count <= 1; }
};

} // namespace chazelle
