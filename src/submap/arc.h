#pragma once

// [C91 §2.2 tex 96]: region-arc — "pieces of C and not of C" alternating
// with exit chords along a region's boundary.
// [C91 §2.4 tex 133]: Arc-structure points directly into Polygon's edge
// array.  [C91 §2.4(ii)]: pointer to the tree node whose region the arc bounds.

#include "../polygon/point.h"
#include "../polygon/perturbation.h"
#include "../common.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>

namespace chazelle {

struct Arc {
    std::size_t first_edge = NONE;  // First edge (e₁) in input table.
    Side first_side = LEFT;
    std::size_t last_edge = NONE;   // Last edge (eₜ); == first_edge for single-edge arcs.
    Side last_side = LEFT;

    // [C91 §2.4(ii)]: Tree node of the region this arc bounds.
    std::size_t region_node = NONE;

    // [C91 §2.2]: Cached nonnull edge count for weight computation.
    std::size_t edge_count = 0;

    // [C91 §3]: Underlying range of C edges (ᾱ) — "the portion of C to which
    // an arc α of ∂C corresponds.  An arc may double-back around an
    // endpoint of C, so ᾱ may not always be as 'long' as α."
    //
    // Non-wrapped (first_side == last_side): ᾱ = [min,max] of edges.
    // Wrapped: both legs map to the same C edges, meeting at the
    // turnaround edge ([C91 §2.1 tex 72]: wraps occur at C endpoints).
    //
    // @param c_start  C's start vertex index (RIGHT→LEFT wrap turnaround).
    // @param c_end    C's end vertex index   (LEFT→RIGHT wrap turnaround at c_end - 1).
    std::pair<std::size_t, std::size_t> underlying_edge_range(
            std::size_t c_start, std::size_t c_end) const noexcept {
        if (first_side == last_side) {
            return {std::min(first_edge, last_edge),
                    std::max(first_edge, last_edge)};
        }
        if (first_side == LEFT) {
            // LEFT→RIGHT: legs meet at c_end - 1.
            assert(c_end >= 1 &&
                   "wrap requires c_end >= 1 (C has ≥ 2 vertices)");
            assert(first_edge <= c_end - 1 && last_edge <= c_end - 1 &&
                   "LEFT→RIGHT wrap: both legs ≤ turnaround edge (c_end - 1)");
            return {std::min(first_edge, last_edge), c_end - 1};
        } else {
            // RIGHT→LEFT: legs meet at c_start.
            assert(first_edge >= c_start && last_edge >= c_start &&
                   "RIGHT→LEFT wrap: both legs ≥ turnaround edge (c_start)");
            return {c_start, std::max(first_edge, last_edge)};
        }
    }

    // [C91 §3.3]: Tombstone for O(1) removal; stripped by compact().
    bool dead = false;
};

} // namespace chazelle
