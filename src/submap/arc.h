#pragma once

// [C91 §2.2 tex 96]: region-arc — "pieces of C and not of C" alternating
// with exit chords along a region's boundary.
// [C91 §2.4 tex 133]: Arc-structure points directly into Polygon's edge
// array.  [C91 §2.4(ii)]: pointer to the tree node whose region the arc bounds.
//
// [C91 §2.4 tex 142]: "caution must be used since an arc might wrap
// around both sides of C, something we call double-backing."  A wrap-
// spanning arc is stored as ONE arc-structure — never split at C's
// endpoints — so the arc-structure count of a region equals its paper
// arc count (= node degree, ≤ 4 for conformal submaps, [C91 §2.3 tex
// 114]).  The (first, last) pointer pair encodes every class:
//
//   no wrap        first_side == last_side; LEFT ascends (first ≤ last),
//                  RIGHT descends (first ≥ last)  [C91 §2.4 tex 133].
//   end wrap       first_side == LEFT, last_side == RIGHT — double-backs
//                  around C's END vertex.
//   start wrap     first_side == RIGHT, last_side == LEFT — double-backs
//                  around C's START vertex.
//   double wrap    first_side == last_side with the traversal order
//                  VIOLATED (LEFT with first > last, RIGHT with
//                  first < last) — wraps around BOTH endpoints of C,
//                  covering the whole opposite side.  (Same-edge
//                  ambiguity cannot arise: a double-wrap arc's two
//                  endpoints are chord endpoints, every chord joins two
//                  distinct (edge, side) positions, and all other chord
//                  endpoints would lie interior to the arc — so its two
//                  bounding endpoints sit on distinct (edge, side)
//                  positions.)
//   closed         the single arc of a chordless submap (the entire
//                  ∂C).  By convention it is cut at the canonical
//                  traversal origin — C's start turnaround ([C91
//                  §2.4(iii) tex 138]) — and stored as an end-wrap arc
//                  with first_edge == last_edge == C's first edge.

#include "../polygon/point.h"
#include "../polygon/perturbation.h"
#include "../polygon/polygon.h"
#include "../common.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>

namespace chazelle {

// One maximal single-side piece of an arc's ∂C coverage, at edge
// granularity ([C91 §2.4 tex 142]: reconstruction "between the locations
// indicated by the two pointers", detecting double-backing at C's
// endpoint edges).
struct ArcLeg {
    Side side = LEFT;
    std::size_t lo = NONE;   // inclusive edge range, ascending
    std::size_t hi = NONE;
};

struct Arc {
    std::size_t first_edge = NONE;  // First edge (e₁) in input table.
    Side first_side = LEFT;
    std::size_t last_edge = NONE;   // Last edge (eₜ); == first_edge for single-edge arcs.
    Side last_side = LEFT;

    // [C91 §2.4(ii)]: Tree node of the region this arc bounds.
    std::size_t region_node = NONE;

    // [C91 §2.2]: Cached nonnull edge count for weight computation.
    std::size_t edge_count = 0;

    // [C91 §2.4 tex 142]: does the (first, last) encoding double-back
    // around an endpoint of C?  (Closed arcs report true: they wrap
    // both endpoints.)
    bool wraps() const noexcept {
        if (first_side != last_side) return true;
        // Double wrap: traversal order violated ([C91 §2.4 tex 133]:
        // LEFT ascends, RIGHT descends).
        if (first_side == LEFT)  return first_edge > last_edge;
        return first_edge < last_edge;
    }

    // Double-backs around C's END vertex ([C91 §2.4 tex 142]).
    bool wraps_end() const noexcept {
        if (first_side == LEFT && last_side == RIGHT) return true;   // end wrap / closed
        return first_side == last_side && wraps();                    // double wrap
    }

    // Double-backs around C's START vertex.  Caveat: the closed arc of
    // a chordless submap shares the end-wrap encoding with first_edge
    // == last_edge == C's first edge, and only IT passes the start
    // turnaround — callers distinguish by context (a chordless submap
    // has exactly one arc).
    bool wraps_start() const noexcept {
        if (first_side == RIGHT && last_side == LEFT) return true;   // start wrap
        return first_side == last_side && wraps();                    // double wrap
    }

    // [C91 §2.4 tex 142]: decompose the arc into its 1–3 single-side
    // legs.  @param c_start  C's start vertex index (= C's first edge).
    //        @param c_end    C's end vertex index (last edge = c_end−1).
    // Legs are emitted in clockwise ∂C traversal order.
    std::size_t legs(std::size_t c_start, std::size_t c_end,
                     ArcLeg out[3]) const noexcept {
        assert(c_end >= 1 && "[C91 §2.1]: C has ≥ 2 vertices");
        const std::size_t last_c_edge = c_end - 1;
        if (first_side == last_side) {
            if (!wraps()) {
                out[0] = ArcLeg{first_side,
                                std::min(first_edge, last_edge),
                                std::max(first_edge, last_edge)};
                return 1;
            }
            // Double wrap ([C91 §2.4 tex 142]: around both endpoints).
            if (first_side == LEFT) {
                assert(first_edge > last_edge);
                out[0] = ArcLeg{LEFT, first_edge, last_c_edge};
                out[1] = ArcLeg{RIGHT, c_start, last_c_edge};
                out[2] = ArcLeg{LEFT, c_start, last_edge};
            } else {
                assert(first_edge < last_edge);
                out[0] = ArcLeg{RIGHT, c_start, first_edge};
                out[1] = ArcLeg{LEFT, c_start, last_c_edge};
                out[2] = ArcLeg{RIGHT, last_edge, last_c_edge};
            }
            return 3;
        }
        if (first_side == LEFT) {
            // End wrap: LEFT leg ascends to C's end, RIGHT leg returns.
            assert(first_edge <= last_c_edge && last_edge <= last_c_edge);
            out[0] = ArcLeg{LEFT, first_edge, last_c_edge};
            out[1] = ArcLeg{RIGHT, last_edge, last_c_edge};
            return 2;
        }
        // Start wrap: RIGHT leg descends to C's start, LEFT leg leaves it.
        assert(first_edge >= c_start && last_edge >= c_start);
        out[0] = ArcLeg{RIGHT, c_start, first_edge};
        out[1] = ArcLeg{LEFT, c_start, last_edge};
        return 2;
    }

    // Does the arc's coverage include ∂C position (edge, side), at edge
    // granularity?  ([C91 §3.1 tex 220]: "arc-structures encode on which
    // side(s) of the double boundary the arcs lie.")
    bool covers(std::size_t edge, Side side,
                std::size_t c_start, std::size_t c_end) const noexcept {
        ArcLeg lg[3];
        std::size_t n = legs(c_start, c_end, lg);
        for (std::size_t i = 0; i < n; ++i)
            if (lg[i].side == side && lg[i].lo <= edge && edge <= lg[i].hi)
                return true;
        return false;
    }

    // [C91 §3]: Underlying range of C edges (ᾱ) — "the portion of C to
    // which an arc α of ∂C corresponds.  An arc may double-back around
    // an endpoint of C, so ᾱ may not always be as 'long' as α"
    // ([C91 §3.0 tex 166]).
    //
    // Non-wrapped: ᾱ = [min,max] of edges.  Wrapped: the legs meet at
    // the turnaround ([C91 §2.4 tex 142]); a double-wrap arc covers a
    // full ∂C side, so ᾱ is all of C.
    std::pair<std::size_t, std::size_t> underlying_edge_range(
            std::size_t c_start, std::size_t c_end) const noexcept {
        assert(c_end >= 1 &&
               "wrap ranges require c_end >= 1 (C has ≥ 2 vertices)");
        if (first_side == last_side) {
            if (!wraps())
                return {std::min(first_edge, last_edge),
                        std::max(first_edge, last_edge)};
            // Double wrap: covers all of one side ⟹ ᾱ = all of C.
            return {c_start, c_end - 1};
        }
        if (first_side == LEFT) {
            // End wrap: legs meet at c_end − 1.
            assert(first_edge <= c_end - 1 && last_edge <= c_end - 1 &&
                   "[C91 §2.4 tex 142]: end wrap legs ≤ C's last edge");
            return {std::min(first_edge, last_edge), c_end - 1};
        }
        // Start wrap: legs meet at c_start.
        assert(first_edge >= c_start && last_edge >= c_start &&
               "[C91 §2.4 tex 142]: start wrap legs ≥ C's first edge");
        return {c_start, std::max(first_edge, last_edge)};
    }

    // [C91 §3.3]: Tombstone for O(1) removal; stripped by compact().
    bool dead = false;
};

// [C91 §2.2 tex 106]: "the weight of a region [is] ... the maximum
// number of nonnull length edges in any of its arcs."  An arc is a
// piece of ∂C — the DOUBLE boundary — so a double-backing arc counts
// the ∂C edges of EVERY leg: the two sides of a doubled-over C edge
// are distinct ∂C edges.  (Lemma 2.3's proof, tex 129, divides the
// ≤ 4-fold per-C-vertex multiplicity back out — a step that
// presupposes this double-boundary count; tex 108 likewise treats
// weight as a per-side notion.)
//
// Partial end pieces of nonzero length count as one edge each ([C91
// §2.1 tex 70]: the subdivision's pieces are its edges); end pieces of
// ZERO length are excluded ("NONNULL length edges") — these arise when
// an arc endpoint sits exactly on a turnaround corner, so the
// (first, last) encoding carries a leg with no geometric extent
// ([C91 §2.4 tex 142]).  The corner tests need the arc's endpoint ys,
// which the encoding does not record ([C91 §2.4 tex 133]: "chords take
// care of that") — callers pass them in.  Zero-length arcs ([C91 §2.2
// tex 96]) carry edge_count = 0 by explicit construction and must not
// be routed through this helper.
inline std::size_t arc_boundary_edge_count(
        const Arc& a, const Polygon& C,
        std::size_t c_start, std::size_t c_end,
        SymbolicY start_y, SymbolicY end_y) {
    ArcLeg lg[3];
    const std::size_t n = a.legs(c_start, c_end, lg);
    std::size_t total = 0;
    for (std::size_t i = 0; i < n; ++i)
        total += C.count_nonnull_edges(lg[i].lo, lg[i].hi);
    if (n == 1) return total;

    // Zero-coverage end edges ([C91 §2.2 tex 106]: an arc's weight
    // counts an edge only where the arc covers nonzero length of its
    // face).  The FIRST edge contributes nothing when the arc starts
    // exactly at the vertex where that face's clockwise traversal
    // EXITS it; the LAST edge contributes nothing when the arc ends
    // exactly at the vertex where the traversal ENTERS it.  The
    // classic case is an endpoint at one of C's turnaround corners
    // ([C91 §2.4 tex 142]); after a null-tail extension
    // ([C91 §4 tex 316], up_phase.cpp) the identical geometry appears
    // at the OLD corner, by then an ordinary interior vertex — the
    // rule depends only on the endpoint's position within its edge,
    // not on the vertex being a corner.
    {
        const auto& fe = C.edge(a.first_edge);
        const SymbolicY f_exit = symbolic_y_of(C.vertex(
            a.first_side == LEFT ? fe.end_idx : fe.start_idx));
        if (symbolic_y_equal(start_y, f_exit))
            total -= C.count_nonnull_edges(a.first_edge, a.first_edge);
        const auto& le = C.edge(a.last_edge);
        const SymbolicY l_entry = symbolic_y_of(C.vertex(
            a.last_side == LEFT ? le.start_idx : le.end_idx));
        if (symbolic_y_equal(end_y, l_entry))
            total -= C.count_nonnull_edges(a.last_edge, a.last_edge);
    }
    return total;
}

} // namespace chazelle
