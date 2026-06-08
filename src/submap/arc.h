#pragma once

/// [C91 §2.2] (tex 96): region-arc — "pieces of C and not of C"
/// alternating with exit chords along a region's boundary.
///
/// [C91 §2.4] (tex 133): Arc-structure — the in-memory representation.
/// "Each arc is represented by a separate arc-structure." Points
/// directly into the input table (Polygon's edge array).
///
/// §2.4(ii): "each arc-structure has a pointer to the node of the tree
/// whose corresponding region is incident upon the arc."

#include "../polygon/point.h"
#include "../polygon/perturbation.h"
#include "../common.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>

namespace chazelle {

struct Arc {
    /// Index of the first edge (e₁) in the input table.
    std::size_t first_edge = NONE;
    /// Which side of the double boundary for e₁.
    Side first_side = LEFT;

    /// Index of the last edge (eₜ) in the input table.
    /// For a single-edge arc (t=1): last_edge == first_edge.
    std::size_t last_edge = NONE;
    /// Which side of the double boundary for eₜ.
    Side last_side = LEFT;

    /// [C91 §2.4(ii)]: "a pointer to the node of the tree whose
    /// corresponding region is incident upon the arc."
    std::size_t region_node = NONE;

    /// [C91 §2.2]: Cached nonnull edge count for weight computation.
    std::size_t edge_count = 0;

    /// [C91 §2.4]: Symbolic y-coordinate of the arc's starting
    /// position on ∂C.  Used for disambiguation in double_identify
    /// when multiple arcs share the same edge ("sorting the endpoints
    /// falling within the same edges by considering y-coordinates").
    double key_y = 0.0;
    std::size_t key_y_tag = NONE;

    SymbolicY key_symbolic_y() const noexcept {
        return {key_y, key_y_tag};
    }

    /// [C91 §3]: The underlying range of C edges for this arc (ᾱ).
    ///
    /// "We introduce the notation ᾱ to refer to the portion of C to
    /// which an arc α of ∂C corresponds.  Recall that an arc may
    /// double-back around an endpoint of C, so ᾱ may not always be
    /// as 'long' as α."
    ///
    /// Non-wrapped (first_side == last_side):
    ///   ᾱ = [min(first_edge, last_edge), max(first_edge, last_edge)]
    ///
    /// Wrapped (first_side != last_side, double-backing):
    ///   The arc wraps around a C endpoint.  ᾱ is the union of both
    ///   legs, which meet at the turnaround edge.  The range on C
    ///   covers both legs and is strictly shorter than the ∂C span.
    ///
    /// @param c_start  Index of C's start vertex (for LEFT→RIGHT wrap)
    /// @param c_end    Index of C's end vertex (for RIGHT→LEFT wrap)
    /// @return (lo_edge, hi_edge) inclusive range on C.
    std::pair<std::size_t, std::size_t> underlying_edge_range(
            std::size_t c_start, std::size_t c_end) const noexcept {
        if (first_side == last_side) {
            // Non-wrapped: straightforward min/max.
            return {std::min(first_edge, last_edge),
                    std::max(first_edge, last_edge)};
        }
        // Wrapped (double-backing): both legs map to the same C edges.
        // LEFT→RIGHT wraps at c_end vertex: turnaround = c_end - 1.
        // RIGHT→LEFT wraps at c_start vertex: turnaround = c_start.
        if (first_side == LEFT) {
            // LEFT→RIGHT: first_edge ascending, last_edge descending,
            // meet at c_end - 1.
            // [C91 §2.1 tex 72]: wrap occurs at a C endpoint; both legs
            // must lie within the C-side bounded by the turnaround edge.
            assert(c_end >= 1 &&
                   "wrap requires c_end >= 1 (C has at least 2 vertices)");
            assert(first_edge <= c_end - 1 && last_edge <= c_end - 1 &&
                   "LEFT→RIGHT wrap: both legs must end at or before "
                   "the turnaround edge (c_end - 1)");
            std::size_t turnaround = c_end - 1;
            return {std::min(first_edge, last_edge), turnaround};
        } else {
            // RIGHT→LEFT: first_edge descending, last_edge ascending,
            // meet at c_start.
            // [C91 §2.1 tex 72]: wrap occurs at a C endpoint; both legs
            // must lie within the C-side bounded by the turnaround edge.
            assert(first_edge >= c_start && last_edge >= c_start &&
                   "RIGHT→LEFT wrap: both legs must start at or after "
                   "the turnaround edge (c_start)");
            std::size_t turnaround = c_start;
            return {turnaround, std::max(first_edge, last_edge)};
        }
    }

    /// Tombstone flag for O(1) removal (§3.3).  Dead arcs remain in
    /// the arc-sequence table to keep indices stable; they are stripped
    /// by compact().
    bool dead = false;
};

} // namespace chazelle
