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

#include <cmath>
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

    /// Virtual-edge y-coordinate (§4.2).
    ///
    /// When finite (not NaN), this arc represents a "tilted" exit-chord
    /// edge rather than a real polygon edge.  The chord is a nearly-
    /// horizontal segment at height virtual_y connecting a point on
    /// polygon edge first_edge to a point on polygon edge last_edge.
    /// During ray-shooting, intersection is computed with this tilted
    /// segment instead of with the polygon edges.
    ///
    /// Per Chazelle §4.2: "tilting the edges a₁b₁ and a₂b₂
    /// symbolically to keep the merging algorithm from complaining
    /// later."
    double virtual_y = NAN;

    /// Is this a virtual (exit-chord-turned-edge) arc?
    bool is_virtual() const { return !std::isnan(virtual_y); }

    /// Is this a null-length (zero-length) arc?
    bool is_null() const { return first_edge == last_edge && edge_count == 0; }

    /// Is this a single-edge arc?
    bool is_single_edge() const { return first_edge == last_edge && edge_count <= 1; }
};

} // namespace chazelle
