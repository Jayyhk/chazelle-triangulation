#pragma once

/// Chord — a horizontal segment connecting two mutually visible points on ∂C.
///
/// Each chord is also a tree edge in the submap tree. The two regions it
/// separates are the tree nodes at its endpoints.
///
/// Per Chazelle §2.3, a chord record stores pointers to the 2, 3, or 4
/// adjacent arc-structures:
///   - 2 arcs: standard case (one arc on each side of the chord).
///   - 3 arcs: one chord endpoint is at a local y-extremum, producing
///     companion/duplicate vertices → two arcs on one side.
///   - 4 arcs: both chord endpoints are at local extrema.

#include "common.h"

#include <cstddef>
#include <array>

namespace chazelle {

struct Chord {
    /// The two tree-node indices (regions) that this chord separates.
    /// These are indices into Submap::nodes_.
    std::size_t region[2] = {NONE, NONE};

    /// Adjacent arc-structure indices (into Submap::arc_sequence_).
    /// Up to 4; unused slots are NONE.
    std::array<std::size_t, 4> adj_arcs = {NONE, NONE, NONE, NONE};

    /// Number of valid entries in adj_arcs (2, 3, or 4).
    std::size_t num_adj_arcs = 0;

    /// Whether this chord is null-length (occurs at local extrema).
    bool is_null_length = false;

    /// Y-coordinate of the chord (for horizontal chords, this is exact).
    /// For implementation convenience; strictly, the chord's position is
    /// implicit from the arc-structures.
    double y = 0.0;

    /// The two endpoint vertex indices on ∂C (into the polygon's vertex array).
    /// left_vertex is the left endpoint, right_vertex is the right endpoint.
    std::size_t left_vertex  = NONE;
    std::size_t right_vertex = NONE;
};

} // namespace chazelle
