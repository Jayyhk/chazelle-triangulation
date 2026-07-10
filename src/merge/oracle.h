#pragma once

// [C91 §3]: Oracle primitives for submap merging — ray-shooter ([C91 §3.0(i)])
// and arc-cutter ([C91 §3.0(ii)]).  Abstract interfaces; implementations in
// [C91 §3.4] and the up-phase ([C91 §4]).

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
// clockwise ∂C order.  Holds for both oracle contexts: [C91 §3.0(i)] subarcs
// inherit α's clockwise traversal ([C91 §2.4 tex 133/138]); [C91 §3.0(ii)] cut() pieces
// require it explicitly (tex 170 condition (1)).
struct Subarc {
    std::size_t first_edge;
    Side first_side;
    std::size_t last_edge;
    Side last_side;
    // [C91 §3.0(i) tex 169]: α' is "specified by its two endpoints
    // (along with two pointers to the input table to indicate the names
    // of the edges of P that contain these two endpoints as well as two
    // flags indicating which side of ∂P is to be understood)" — the
    // endpoints are exact ∂C positions; the edge pointers + side flags
    // above only locate them.  Every subarc passed to
    // RayShootingOracle::shoot MUST carry real endpoint ys
    // (tag != SOS_NONE).  ArcPiece.subarc endpoints may keep default
    // ys — pieces are never shot today (their ys arrive with §4.1's
    // production cutter).
    SymbolicY first_y{};
    SymbolicY last_y{};
};

// [C91 §2.4 tex 138]: A Subarc respects clockwise ∂C ordering — LEFT
// side ascends in edge, RIGHT side descends.  Wrap subarcs inherit α's
// double-backing per [C91 §2.4 tex 142]: first_side ≠ last_side wraps
// one endpoint of C; a same-side pair VIOLATING the traversal order
// wraps both (the subarc covers the whole opposite side), so no
// same-side edge-index check is possible — the encoding is exactly
// Arc's (arc.h).
inline void assert_subarc_clockwise(const Subarc&) {}

// [C91 §2.4 tex 142]: decompose a Subarc into its 1–3 single-side legs —
// the Subarc (first, last) encoding is Arc's, so reuse Arc::legs.
inline std::size_t subarc_legs(const Subarc& s, std::size_t c_start,
                               std::size_t c_end, ArcLeg out[3]) {
    Arc a;
    a.first_edge = s.first_edge;
    a.first_side = s.first_side;
    a.last_edge = s.last_edge;
    a.last_side = s.last_side;
    a.region_node = 0;   // unused by legs()
    return a.legs(c_start, c_end, out);
}

// [C91 §3.0(i) tex 169]: position-exact membership — is the ∂C point
// (edge, side, y) on α′?  Leg coverage as subarc_covers_position, plus
// within-edge clipping at α′'s two endpoint faces: subarc_legs returns
// legs in traversal order, so leg 0 begins at (first_edge, first_side)
// and the final leg ends at (last_edge, last_side); within an edge the
// ∂C traversal follows the edge on LEFT and reverses it on RIGHT
// ([C91 §2.4(iii) tex 138]), so the within-face order of two symbolic
// ys is the y order iff the face traversal ascends.
inline bool subarc_contains_point(const Subarc& s, const Polygon& C,
                                  std::size_t edge, Side side, SymbolicY y,
                                  std::size_t c_start, std::size_t c_end) {
    ArcLeg legs[3];
    std::size_t n = subarc_legs(s, c_start, c_end, legs);
    auto face_trav_ascends = [&](std::size_t e_, Side s_) {
        const auto& ed = C.edge(e_);
        bool asc = symbolic_y_less(symbolic_y_of(C.vertex(ed.start_idx)),
                                   symbolic_y_of(C.vertex(ed.end_idx)));
        return (s_ == LEFT) ? asc : !asc;
    };
    for (std::size_t g = 0; g < n; ++g) {
        if (legs[g].side != side || edge < legs[g].lo || edge > legs[g].hi)
            continue;
        if (g == 0 && edge == s.first_edge && side == s.first_side) {
            assert(s.first_y.tag != SOS_NONE &&
                   "[C91 §3.0(i) tex 169]: shot subarcs carry exact endpoints");
            bool ok = face_trav_ascends(edge, side)
                          ? symbolic_y_geq(y, s.first_y)
                          : symbolic_y_leq(y, s.first_y);
            if (!ok) continue;
        }
        if (g + 1 == n && edge == s.last_edge && side == s.last_side) {
            assert(s.last_y.tag != SOS_NONE &&
                   "[C91 §3.0(i) tex 169]: shot subarcs carry exact endpoints");
            bool ok = face_trav_ascends(edge, side)
                          ? symbolic_y_leq(y, s.last_y)
                          : symbolic_y_geq(y, s.last_y);
            if (!ok) continue;
        }
        return true;
    }
    return false;
}

// [C91 §3.0(i) tex 169]: Ray-hit report — point + edge of P containing it.
struct RayHit {
    bool hit = false;
    double x = 0.0;
    double y = 0.0;                 // = p.y (horizontal ray).
    std::size_t edge = 0;
    Side side = LEFT;
    // [C91 §2.1 tex 70]: "If it were to go to infinity in the Cartesian
    // plane, then it would actually wrap around in the spherical plane
    // until it hits C again."  True iff the reported hit was reached
    // through the point at infinity — its x lies BEHIND the source in
    // the travel direction, and every direct (non-wrapped) hit is
    // nearer than every wrapped one.
    bool wrapped = false;
    // Index of the struck region arc.  Not part of the paper's spec
    // (tex 169 only mandates point + edge); set by local_shoot() after
    // the oracle call.  Oracle implementations MUST NOT set this.
    std::size_t hit_arc_idx = NONE;
};

// [C91 §3.0(ii) tex 170]: One piece αⱼ returned by cut().
// Conditions (1)–(3) summarized; full quotes in `assert_cut_postconditions`.
struct ArcPiece {
    Subarc subarc;                  // (1): location on ∂C, clockwise.
    const Submap* submap = nullptr; // (3): conformal V(ᾱⱼ); null for boundary pieces.
    const Polygon* curve = nullptr; // (3): underlying ᾱⱼ; null for boundary pieces.
    bool is_boundary_piece = false; // true ⟹ single-edge ᾱ₁/ᾱ₂ at α''s endpoints.
    // (3): the γⱼ for which `submap` is γⱼ-granular, with γⱼ ≤ h(γᵢ).
    // tex 170 says "an h(γᵢ)-granular conformal submap of V(ᾱⱼ) is
    // available", but granularity is not monotone in γ (criterion (ii)
    // of [C91 §2.3 tex 121] demands MORE contraction at larger γ), and
    // the paper's own oracle hands out pieces of differing
    // granularities: [C91 §4.1 tex 352] "we have conformal submaps of
    // granularity AT MOST 2^⌈β⌈βλ⌉⌉" — the up-phase pieces are chains
    // in grades μ ≤ ⌈βλ⌉ whose canonical submaps are 2^⌈βμ⌉-granular
    // ([C91 §4.1 tex 338]).  h(γᵢ) is the uniform bound (tex 353:
    // "h(γ) ≤ 2^⌈β⌈βλ⌉⌉") entering the §3.2 complexity ("there are only
    // O(h(γ₂)) vertices in the region", tex 252), which needs region
    // weights ≤ h — implied by γⱼ-granularity with γⱼ ≤ h.
    std::size_t granularity = 0;
};

// [C91 §3.0(i) tex 166–169]: Ray-shooter — "given any point p along with
// a horizontal direction and any subarc α' of α, reports the single point
// of α' (if any) that a ray from p would hit."  O(f(γᵢ)) per query.
// Implementation: [C91 §3.4].
struct RayShootingOracle {
    virtual ~RayShootingOracle() = default;
    // @param arc_idx  region arc α's index in Sᵢ's arc-sequence
    //                 (tex 166: "specified by a pointer to its arc-structure").
    // @param source_x_offset  the SOURCE wall's perturbed x-offset
    //                 (polygon.h::perturbed_x_offset — how far the SoS
    //                 perturbation really moves the point off its raw
    //                 position, computed once by the original caller in
    //                 the source's own frame; the offset is plain
    //                 geometry, identical in every subchain frame), or
    //                 SOURCE_OFFSET_NONE.  [C91 §2 tex 47]: a contact
    //                 at raw distance 0 shares the source's raw point;
    //                 whether it is strictly FORWARD (direct) or met
    //                 after a full wrap is decided by comparing the two
    //                 perturbed x-offsets
    //                 (polygon.h::perturbed_hit_forward).  Without the
    //                 source's offset the historic convention applies
    //                 (d = 0 ⟹ wrapped), exact when p is a wall of its
    //                 own duplicate pair.
    virtual RayHit shoot(Point p, Side direction,
                         std::size_t arc_idx,
                         const Subarc& target,
                         double source_x_offset = SOURCE_OFFSET_NONE) const = 0;
};

// [C91 §3.0(ii) tex 170]: Arc-cutter — subdivides α' into at most g(γᵢ)
// pieces in O(g(γᵢ)) time, satisfying conditions (1)–(3).
// Implementation: [C91 §4] (up-phase).
struct ArcCuttingOracle {
    virtual ~ArcCuttingOracle() = default;
    // Callers MUST validate via `assert_cut_postconditions` before consuming.
    // @param arc_idx  region arc α's index in Sᵢ's arc-sequence (tex 166).
    virtual std::vector<ArcPiece> cut(std::size_t arc_idx,
                                      const Subarc& target) const = 0;
};

// [C91 §3.0(ii) tex 170]: Enforces cut()'s post-conditions.
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
// Callers in [C91 §3.2 / §3.3] must invoke this on every `cut` result to
// verify the §3.0(ii) tex 170 post-conditions before consuming the pieces.
inline void assert_cut_postconditions(
        [[maybe_unused]] const Polygon& Ci,
        [[maybe_unused]] const Subarc& target,
        const ArcPiece* pieces, std::size_t count,
        [[maybe_unused]] std::size_t max_pieces,
        [[maybe_unused]] std::size_t h_gamma) {
    assert(count >= 1 && "[C91 §3.0(ii) tex 170]: cut() must produce ≥1 piece");
    assert(count <= max_pieces &&
           "[C91 §3.0(ii) tex 170]: cut() must produce ≤ g(γᵢ) pieces");

    // Subdivision identity: α' = α₁ ∪ ... ∪ αₖ.
    assert(pieces[0].subarc.first_edge == target.first_edge &&
           pieces[0].subarc.first_side == target.first_side &&
           "[C91 §3.0(ii) tex 170]: first piece must start at α'.first");
    assert(pieces[count - 1].subarc.last_edge == target.last_edge &&
           pieces[count - 1].subarc.last_side == target.last_side &&
           "[C91 §3.0(ii) tex 170]: last piece must end at α'.last");

    for (std::size_t j = 0; j < count; ++j) {
        const ArcPiece& p = pieces[j];

        // (2): no double-backing.
        assert(p.subarc.first_side == p.subarc.last_side &&
               "[C91 §3.0(ii)(2) tex 170]: each piece must lie on one side of C");

        // (1): clockwise — LEFT ascends, RIGHT descends ([C91 §2.4 tex 138]).
        if (p.subarc.first_side == LEFT) {
            assert(p.subarc.first_edge <= p.subarc.last_edge &&
                   "[C91 §3.0(ii)(1) tex 170]: LEFT piece must ascend in edge");
        } else {
            assert(p.subarc.first_edge >= p.subarc.last_edge &&
                   "[C91 §3.0(ii)(1) tex 170]: RIGHT piece must descend in edge");
        }

        // (3): boundary pieces only at j=0 or j=count-1 (the at-most-two
        // endpoints of α').
        const bool at_endpoint = (j == 0 || j + 1 == count);
        if (!at_endpoint) {
            assert(!p.is_boundary_piece &&
                   "[C91 §3.0(ii)(3) tex 170]: only first/last pieces may be boundary");
        }

        if (!p.is_boundary_piece) {
            // (3): h(γᵢ)-granular conformal submap of V(ᾱⱼ) in normal form.
            assert(p.submap != nullptr &&
                   "[C91 §3.0(ii)(3) tex 170]: non-boundary piece requires a submap");
            assert(p.curve != nullptr &&
                   "[C91 §3.0(ii)(3) tex 170]: non-boundary piece requires its curve");
            // (3): "vertex-to-vertex subchains ... do not stop in the middle of
            // an edge."  ᾱⱼ = p.curve must cover exactly the piece's polygon
            // edge range.  (2) above guarantees first_side == last_side, so the
            // piece is single-side and its edge count is |last - first| + 1.
            {
                [[maybe_unused]] std::size_t lo =
                    std::min(p.subarc.first_edge, p.subarc.last_edge);
                [[maybe_unused]] std::size_t hi =
                    std::max(p.subarc.first_edge, p.subarc.last_edge);
                assert(p.curve->num_edges() == hi - lo + 1 &&
                       "[C91 §3.0(ii)(3) tex 170]: non-boundary piece is "
                       "vertex-to-vertex; p.curve must cover exactly the "
                       "piece's polygon edge range");
                // (3): "vertex-to-vertex subchains OF Cᵢ" — identity, not
                // just length: ᾱⱼ must be THE subchain of Cᵢ spanning the
                // piece's edge range ([C91 §2.4 tex 133]: SoS indices are
                // consecutive along Cᵢ, so matching both ends pins it).
                assert(p.curve->vertex(0).index == Ci.vertex(lo).index &&
                       p.curve->vertex(p.curve->num_vertices() - 1).index ==
                           Ci.vertex(hi + 1).index &&
                       "[C91 §3.0(ii)(3) tex 170]: ᾱⱼ must be the "
                       "vertex-to-vertex subchain of Cᵢ at the piece's "
                       "edge range");
            }
            // [C91 §2.4(iv) tex 139]: conformal ⟹ tree decomposition available.
            assert(!p.submap->tree_decomposition().empty() &&
                   "[C91 §2.4(iv) tex 139]: normal-form conformal submap "
                   "needs its tree decomposition");
            // (3) + [C91 §4.1 tex 352]: each piece submap is γⱼ-granular
            // for a declared γⱼ ≤ h(γᵢ) (see ArcPiece::granularity).
            assert(p.granularity >= 1 && p.granularity <= h_gamma &&
                   "[C91 §4.1 tex 352]: piece granularity is at most h(γᵢ)");
#ifdef CHAZELLE_EXPENSIVE_ASSERTS
            // Gated: O(m) per piece × g(γᵢ) pieces would blow cut()'s
            // O(g(γᵢ)) budget at [C91 §3.0(ii) tex 170].
            p.submap->check_invariants(*p.curve);
            assert(p.submap->is_conformal() &&
                   "[C91 §3.0(ii)(3) tex 170]: non-boundary piece must be conformal");
            assert(p.submap->is_granular(p.granularity, *p.curve) &&
                   "[C91 §3.0(ii)(3) tex 170]: non-boundary piece must be "
                   "γⱼ-granular for its declared γⱼ");
            assert(p.submap->is_semigranular(h_gamma) &&
                   "[C91 §3.2 tex 252]: O(h(γᵢ)) vertices per piece region "
                   "— weights bounded by h(γᵢ)");
#endif
        } else {
            // (3): boundary pieces are single-edge at α''s endpoints.
            assert(p.subarc.first_edge == p.subarc.last_edge &&
                   "[C91 §3.0(ii)(3) tex 170]: boundary piece must be single-edge");
        }
    }
}

} // namespace chazelle
