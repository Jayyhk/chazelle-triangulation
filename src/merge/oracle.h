#pragma once

// [C91 §3]: Oracle primitives for submap merging — ray-shooter (§3.0(i))
// and arc-cutter (§3.0(ii)).  Abstract interfaces; implementations in
// §3.4 and the up-phase (§4).

#include "../polygon/point.h"
#include "../polygon/polygon.h"
#include "../submap/submap.h"

#include <cassert>
#include <cstddef>
#include <vector>

namespace chazelle {

// [C91 §3.0(i)(ii)]: A contiguous piece of a region arc — endpoints
// specified by (edge of P, ∂P side).
//
// Invariant: (first_edge, first_side) precedes (last_edge, last_side) in
// clockwise ∂C order.  Holds for both oracle contexts: §3.0(i) subarcs
// inherit α's clockwise traversal (§2.4 tex 133/138); §3.0(ii) cut() pieces
// require it explicitly (tex 170 condition (1)).
struct Subarc {
    std::size_t first_edge;
    Side first_side;
    std::size_t last_edge;
    Side last_side;
};

// [C91 §3.0(i)] (tex 169): Ray-hit report — point + edge of P containing it.
struct RayHit {
    bool hit = false;
    double x = 0.0;
    double y = 0.0;                 // = p.y (horizontal ray).
    std::size_t edge = 0;
    Side side = LEFT;
    // Index of the struck region arc.  Not part of the paper's spec
    // (tex 169 only mandates point + edge); set by local_shoot() after
    // the oracle call.  Oracle implementations MUST NOT set this.
    std::size_t hit_arc_idx = NONE;
};

// [C91 §3.0(ii)] (tex 170): One piece αⱼ returned by cut().
// Conditions (1)–(3) summarized; full quotes in `assert_cut_postconditions`.
struct ArcPiece {
    Subarc subarc;                  // (1): location on ∂C, clockwise.
    const Submap* submap = nullptr; // (3): h(γᵢ)-granular conformal V(ᾱⱼ); null for boundary pieces.
    const Polygon* curve = nullptr; // (3): underlying ᾱⱼ; null for boundary pieces.
    bool is_boundary_piece = false; // true ⟹ single-edge ᾱ₁/ᾱ₂ at α''s endpoints.
};

// [C91 §3.0(i)] (tex 166–169): Ray-shooter — "given any point p along with
// a horizontal direction and any subarc α' of α, reports the single point
// of α' (if any) that a ray from p would hit."  O(f(γᵢ)) per query.
// Implementation: §3.4.
struct RayShootingOracle {
    virtual ~RayShootingOracle() = default;
    // @param arc_idx  region arc α's index in Sᵢ's arc-sequence
    //                 (tex 166: "specified by a pointer to its arc-structure").
    virtual RayHit shoot(Point p, Side direction,
                         std::size_t arc_idx,
                         const Subarc& target) const = 0;
};

// [C91 §3.0(ii)] (tex 170): Arc-cutter — subdivides α' into at most g(γᵢ)
// pieces in O(g(γᵢ)) time, satisfying conditions (1)–(3).
// Implementation: §4 (up-phase).
struct ArcCuttingOracle {
    virtual ~ArcCuttingOracle() = default;
    // Callers MUST validate via `assert_cut_postconditions` before consuming.
    // @param arc_idx  region arc α's index in Sᵢ's arc-sequence (tex 166).
    virtual std::vector<ArcPiece> cut(std::size_t arc_idx,
                                      const Subarc& target) const = 0;
};

// [C91 §3.0(ii)] (tex 170): Enforces cut()'s post-conditions.
//
// Bound (verbatim): "in O(g(γᵢ)) time, subdivides the subarc α' into at
//                    most g(γᵢ) subarcs α₁, α₂, ..."
// (1) (verbatim):   "the pair should be ordered to reflect a clockwise
//                    traversal along ∂P and two flags should be included
//                    to indicate on which sides of ∂P these endpoints fall"
// (2) (verbatim):   "the relative interior of no ᾱⱼ should contain a point
//                    of ∂Cᵢ that corresponds to an endpoint of Cᵢ, that is,
//                    each subarc must lie entirely on one side of C
//                    (no double-backing)"
// (3) (verbatim):   "except for ᾱ₁ and ᾱ₂ (in the case where these are
//                    single-edge pieces attached to the points of C
//                    corresponding to the endpoints of α'), all the ᾱⱼ are
//                    vertex-to-vertex subchains of Cᵢ ... and, for each of
//                    them, an h(γᵢ)-granular conformal submap of V(ᾱⱼ) is
//                    available in normal form."
//
// Plus the implicit subdivision identity α' = α₁ ∪ ... ∪ αₖ:
//   first piece's first endpoint = α'.first;
//   last piece's last endpoint   = α'.last.
//
// Callers in §3.2 / §3.3 must invoke this on every `cut` result to ensure
// the oracle's contract holds before consuming the pieces.
inline void assert_cut_postconditions(
        const Subarc& target,
        const ArcPiece* pieces, std::size_t count,
        std::size_t max_pieces,
        [[maybe_unused]] std::size_t h_gamma) {
    assert(count >= 1 && "§3.0(ii) tex 170: cut() must produce ≥1 piece");
    assert(count <= max_pieces &&
           "§3.0(ii) tex 170: cut() must produce ≤ g(γᵢ) pieces");

    // Subdivision identity: α' = α₁ ∪ ... ∪ αₖ.
    assert(pieces[0].subarc.first_edge == target.first_edge &&
           pieces[0].subarc.first_side == target.first_side &&
           "§3.0(ii) tex 170: first piece must start at α'.first");
    assert(pieces[count - 1].subarc.last_edge == target.last_edge &&
           pieces[count - 1].subarc.last_side == target.last_side &&
           "§3.0(ii) tex 170: last piece must end at α'.last");

    for (std::size_t j = 0; j < count; ++j) {
        const ArcPiece& p = pieces[j];

        // (2): no double-backing.
        assert(p.subarc.first_side == p.subarc.last_side &&
               "§3.0(ii)(2) tex 170: each piece must lie on one side of C");

        // (1): clockwise — LEFT ascends, RIGHT descends (§2.4 tex 138).
        if (p.subarc.first_side == LEFT) {
            assert(p.subarc.first_edge <= p.subarc.last_edge &&
                   "§3.0(ii)(1) tex 170: LEFT piece must ascend in edge");
        } else {
            assert(p.subarc.first_edge >= p.subarc.last_edge &&
                   "§3.0(ii)(1) tex 170: RIGHT piece must descend in edge");
        }

        // (3): boundary pieces only at j=0 or j=count-1 (the at-most-two
        // endpoints of α').
        const bool at_endpoint = (j == 0 || j + 1 == count);
        if (!at_endpoint) {
            assert(!p.is_boundary_piece &&
                   "§3.0(ii)(3) tex 170: only first/last pieces may be boundary");
        }

        if (!p.is_boundary_piece) {
            // (3): h(γᵢ)-granular conformal submap of V(ᾱⱼ) in normal form.
            assert(p.submap != nullptr &&
                   "§3.0(ii)(3) tex 170: non-boundary piece requires a submap");
            assert(p.curve != nullptr &&
                   "§3.0(ii)(3) tex 170: non-boundary piece requires its curve");
            // §2.4(iv) tex 139: conformal ⟹ tree decomposition available.
            assert(!p.submap->tree_decomposition().empty() &&
                   "§2.4(iv) tex 139: normal-form conformal submap "
                   "needs its tree decomposition");
#ifdef CHAZELLE_EXPENSIVE_ASSERTS
            // Full normal-form / conformal / granular validation — O(m), gated
            // to keep cut() within its O(g(γᵢ)) paper budget.
            p.submap->check_invariants(*p.curve);
            assert(p.submap->is_conformal() &&
                   "§3.0(ii)(3) tex 170: non-boundary piece must be conformal");
            assert(p.submap->is_granular(h_gamma, *p.curve) &&
                   "§3.0(ii)(3) tex 170: non-boundary piece must be h(γᵢ)-granular");
#endif
        } else {
            // (3): boundary pieces are single-edge at α''s endpoints.
            assert(p.subarc.first_edge == p.subarc.last_edge &&
                   "§3.0(ii)(3) tex 170: boundary piece must be single-edge");
        }
    }
}

} // namespace chazelle
