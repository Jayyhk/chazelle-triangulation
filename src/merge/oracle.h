#pragma once

/// [C91 §3]: Oracle primitives for submap merging.
///
/// "To facilitate the exposition we assume that we have at our
/// disposal two primitives."
///
/// These are abstract interfaces — implementations are provided
/// by §3.4 and the up-phase (§4).

#include "../polygon/point.h"
#include "../polygon/polygon.h"
#include "../submap/submap.h"

#include <cassert>
#include <cstddef>
#include <vector>

namespace chazelle {

/// [C91 §3.0(i)(ii)]: A subarc specification — a contiguous piece of
/// a region arc α of Sᵢ, "specified by its two endpoints (along with
/// two pointers to the input table to indicate the names of the edges
/// of P that contain these two endpoints as well as two flags
/// indicating which side of ∂P is to be understood)."
///
/// Clockwise ordering invariant: `first_edge`/`first_side` precedes
/// `last_edge`/`last_side` in clockwise ∂C order.  This holds for both
/// oracle contexts the struct is shared by, but for distinct paper reasons:
///   - §3.0(i) ray-shooter subarcs of α inherit α's traversal direction,
///     which is itself clockwise by construction (§2.4 tex 138 canonical
///     traversal of ∂C; tex 133 "Let e₁, ..., eₜ be the edges of an arc
///     in clockwise order along the double boundary").
///   - §3.0(ii) cut() pieces are required to be clockwise explicitly by
///     condition (1) (tex 170: "the pair should be ordered to reflect a
///     clockwise traversal along ∂P").
struct Subarc {
    std::size_t first_edge;   ///< Edge of P containing the first endpoint (clockwise-first).
    Side first_side;          ///< Which side of ∂P for the first endpoint.
    std::size_t last_edge;    ///< Edge of P containing the last endpoint (clockwise-last).
    Side last_side;           ///< Which side of ∂P for the last endpoint.
};

/// [C91 §3.0(i)]: Result of a ray-shooting query.
///
/// "reports the single point of α' (if any) that a ray of light shot
/// from p in the given direction would hit in the absence of any
/// obstacle except α'.  In addition to the point hit, the report
/// should also include the name of the edge of P that contains it."
struct RayHit {
    bool hit = false;               ///< True if the ray hit the subarc.
    double x = 0.0;                 ///< x-coordinate of the hit point.
    double y = 0.0;                 ///< y-coordinate of the hit point (= p.y).
    std::size_t edge = 0;           ///< "the name of the edge of P that contains it."
    Side side = LEFT;               ///< Which side of ∂P the hit is on.
    /// Index of the struck region arc.
    /// NOT part of the paper's oracle specification (tex 169 specifies only
    /// hit point + edge name).  Set by local_shoot() after the oracle call;
    /// oracle implementations MUST NOT set this field.
    std::size_t hit_arc_idx = NONE;
};

/// [C91 §3.0(ii)]: Result of an arc-cutting query — one piece αⱼ.
///
/// (1) "each αⱼ is specified by its two endpoints and a pair of
///     pointers into the input table to indicate which edges of P
///     contain these endpoints; the pair should be ordered to reflect
///     a clockwise traversal along ∂P and two flags should be
///     included to indicate on which sides of ∂P these endpoints fall"
///
/// (2) "the relative interior of no ᾱⱼ should contain a point of ∂Cᵢ
///     that corresponds to an endpoint of Cᵢ, that is, each subarc
///     must lie entirely on one side of C (no double-backing)"
///
/// (3) "except for ᾱ₁ and ᾱ₂ (in the case where these are single-
///     edge pieces attached to the points of C corresponding to the
///     endpoints of α'), all the ᾱⱼ are vertex-to-vertex subchains
///     of Cᵢ [...] and, for each of them, an h(γᵢ)-granular conformal
///     submap of V(ᾱⱼ) is available in normal form."
struct ArcPiece {
    /// (1): The piece's location on ∂C, endpoints ordered clockwise.
    Subarc subarc;

    /// (3): h(γᵢ)-granular conformal submap of V(ᾱⱼ) in normal form.
    /// nullptr for the first/last single-edge boundary pieces.
    const Submap* submap = nullptr;

    /// (3): The underlying curve ᾱⱼ (vertex-to-vertex subchain of Cᵢ).
    /// nullptr for the first/last single-edge boundary pieces.
    const Polygon* curve = nullptr;

    /// True if this is a single-edge boundary piece (ᾱ₁ or ᾱ₂ from
    /// condition (3)) that need not have a precomputed submap.
    bool is_boundary_piece = false;
};

/// [C91 §3.0(i)]: Ray-shooting oracle.
///
/// "For each region arc α of Sᵢ (i=1,2) specified by a pointer to
/// its arc-structure: (i) There exists a ray-shooter which, given
/// any point p along with a horizontal direction (left or right) and
/// any subarc α' of α [...], reports the single point of α' (if any)
/// that a ray of light shot from p in the given direction would hit
/// in the absence of any obstacle except α'."
///
/// "The report should take O(f(γᵢ)) time, where f is a nondecreasing
/// function."
///
/// Implementation: §3.4.
struct RayShootingOracle {
    virtual ~RayShootingOracle() = default;

    /// Shoot a horizontal ray from p in the given direction toward
    /// subarc α' of arc α.  Returns the single hit point (if any)
    /// plus the edge of P containing it.
    ///
    /// @param p          The ray origin.
    /// @param direction  LEFT (shoot left) or RIGHT (shoot right).
    /// @param arc_idx    Index of region arc α in Sᵢ's arc-sequence
    ///                   table — "specified by a pointer to its
    ///                   arc-structure" (§3.0, tex 166).
    /// @param target     The subarc α' to shoot toward.
    virtual RayHit shoot(Point p, Side direction,
                         std::size_t arc_idx,
                         const Subarc& target) const = 0;
};

/// [C91 §3.0(ii)]: Arc-cutting oracle.
///
/// "For each region arc α of Sᵢ specified by a pointer to its
/// arc-structure: (ii) There exists an arc-cutter which, in
/// O(g(γᵢ)) time, subdivides the subarc α' into at most g(γᵢ)
/// subarcs α₁, α₂, ..., such that
/// (1) [each specified by endpoints + edges + sides, clockwise order],
/// (2) [no double-backing — each piece on one side of C], and
/// (3) [except boundary pieces, all are vertex-to-vertex subchains
///     with h(γᵢ)-granular conformal submaps in normal form]."
///
/// Implementation: §4 (up-phase).
struct ArcCuttingOracle {
    virtual ~ArcCuttingOracle() = default;

    /// Subdivide subarc α' of arc α into at most g(γᵢ) pieces
    /// satisfying conditions (1)–(3).
    ///
    /// Callers MUST validate the result via `assert_cut_postconditions`
    /// (defined below) before consuming the pieces — that helper enforces
    /// the four §3.0(ii) post-conditions the paper guarantees.
    ///
    /// @param arc_idx Index of region arc α in Sᵢ's arc-sequence
    ///                table — "specified by a pointer to its
    ///                arc-structure" (§3.0, tex 166).
    /// @param target  The subarc α' to cut.
    /// @return Pieces α₁, α₂, ..., ordered along ∂C.
    ///         result.size() must be ≤ g(γᵢ) (tex 175).
    virtual std::vector<ArcPiece> cut(std::size_t arc_idx,
                                      const Subarc& target) const = 0;
};

/// [C91 §3.0(ii)] (tex 170): Asserts the paper-mandated post-conditions on
/// the output of `ArcCuttingOracle::cut`.
///
/// The four conditions from tex 170 (verbatim):
///   - bound: "in O(g(γᵢ)) time, subdivides the subarc α' into at most
///     g(γᵢ) subarcs α₁, α₂, ..."
///   - (1):   "the pair should be ordered to reflect a clockwise traversal
///     along ∂P and two flags should be included to indicate on which sides
///     of ∂P these endpoints fall"
///   - (2):   "the relative interior of no ᾱⱼ should contain a point of
///     ∂Cᵢ that corresponds to an endpoint of Cᵢ, that is, each subarc
///     must lie entirely on one side of C (no double-backing)"
///   - (3):   "except for ᾱ₁ and ᾱ₂ (in the case where these are single-
///     edge pieces attached to the points of C corresponding to the
///     endpoints of α'), all the ᾱⱼ are vertex-to-vertex subchains of Cᵢ
///     ... and, for each of them, an h(γᵢ)-granular conformal submap of
///     V(ᾱⱼ) is available in normal form."
///
/// Plus the implicit subdivision identity α' = α₁ ∪ ... ∪ αₖ:
///   - first piece's first endpoint = α'.first
///   - last  piece's last  endpoint = α'.last
///
/// Callers in §3.2 / §3.3 must invoke this on every `cut` result to ensure
/// the oracle's contract holds before consuming the pieces.
///
/// @param target       The α' that was passed to `cut`.
/// @param pieces       The pieces returned by `cut` (data() + size()).
/// @param count        Number of pieces.
/// @param max_pieces   Upper bound g(γᵢ) on the piece count.
/// @param h_gamma      h(γᵢ) — granularity bound for non-boundary submaps.
inline void assert_cut_postconditions(
        const Subarc& target,
        const ArcPiece* pieces, std::size_t count,
        std::size_t max_pieces,
        [[maybe_unused]] std::size_t h_gamma) {
    // Bound (tex 170): "subdivides ... into at most g(γᵢ) subarcs"
    // implies count ≥ 1 (a subarc cannot be subdivided into zero pieces)
    // and count ≤ g(γᵢ).
    assert(count >= 1 &&
           "§3.0(ii) tex 170: cut() must produce at least one piece");
    assert(count <= max_pieces &&
           "§3.0(ii) tex 170: cut() must produce at most g(γᵢ) pieces");

    // Subdivision identity α' = α₁ ∪ ... ∪ αₖ.
    assert(pieces[0].subarc.first_edge == target.first_edge &&
           pieces[0].subarc.first_side == target.first_side &&
           "§3.0(ii) tex 170: first piece must start at α'.first "
           "(subdivision identity α' = α₁ ∪ ... ∪ αₖ)");
    assert(pieces[count - 1].subarc.last_edge == target.last_edge &&
           pieces[count - 1].subarc.last_side == target.last_side &&
           "§3.0(ii) tex 170: last piece must end at α'.last "
           "(subdivision identity α' = α₁ ∪ ... ∪ αₖ)");

    for (std::size_t j = 0; j < count; ++j) {
        const ArcPiece& p = pieces[j];

        // (2) (tex 170): "each subarc must lie entirely on one side of C
        // (no double-backing)".
        assert(p.subarc.first_side == p.subarc.last_side &&
               "§3.0(ii)(2) tex 170: each piece must lie entirely on one "
               "side of C (no double-backing)");

        // (1) (tex 170): "the pair should be ordered to reflect a clockwise
        // traversal along ∂P".  LEFT side traverses edges ascending; RIGHT
        // side traverses edges descending (§2.4 tex 138 canonical traversal).
        if (p.subarc.first_side == LEFT) {
            assert(p.subarc.first_edge <= p.subarc.last_edge &&
                   "§3.0(ii)(1) tex 170: LEFT-side piece must have "
                   "first_edge ≤ last_edge (clockwise = ascending on LEFT)");
        } else {
            assert(p.subarc.first_edge >= p.subarc.last_edge &&
                   "§3.0(ii)(1) tex 170: RIGHT-side piece must have "
                   "first_edge ≥ last_edge (clockwise = descending on RIGHT)");
        }

        // (3) (tex 170): "except for ᾱ₁ and ᾱ₂ ... all the ᾱⱼ are
        // vertex-to-vertex subchains of Cᵢ ..."  Boundary pieces are
        // permitted only at the sequence endpoints (j=0 or j=count-1)
        // since they correspond to the (at most two) endpoints of α'.
        const bool at_endpoint = (j == 0 || j + 1 == count);
        if (!at_endpoint) {
            assert(!p.is_boundary_piece &&
                   "§3.0(ii)(3) tex 170: only the first and last pieces "
                   "may be boundary pieces (ᾱ₁ and ᾱ₂ in the paper's "
                   "labeling correspond to α''s two endpoints)");
        }

        if (!p.is_boundary_piece) {
            // (3) (tex 170): "an h(γᵢ)-granular conformal submap of V(ᾱⱼ)
            // is available in normal form."  "Available" implies non-null
            // submap + non-null underlying curve.
            assert(p.submap != nullptr &&
                   "§3.0(ii)(3) tex 170: non-boundary piece requires an "
                   "h(γᵢ)-granular conformal submap (non-null)");
            assert(p.curve != nullptr &&
                   "§3.0(ii)(3) tex 170: non-boundary piece requires its "
                   "underlying curve V(ᾱⱼ) (non-null)");
            // §2.4(iv) tex 139: "If the submap is conformal, then its tree
            // decomposition should be available."  Part of "normal form" —
            // O(1) check, kept unguarded.
            assert(!p.submap->tree_decomposition().empty() &&
                   "§2.4(iv) tex 139: non-boundary piece's conformal submap "
                   "must have its tree decomposition available (normal form)");
#ifdef CHAZELLE_EXPENSIVE_ASSERTS
            // "in normal form" — full structural validation (§2.4(i)–(iv));
            // O(m) so gated to avoid inflating cut() beyond paper's
            // O(g(γᵢ)) per call (§3.0(ii) tex 170).
            p.submap->check_invariants(*p.curve);
            // "conformal submap of V(ᾱⱼ)" — O(num_nodes).
            assert(p.submap->is_conformal() &&
                   "§3.0(ii)(3) tex 170: non-boundary piece's submap must "
                   "be conformal");
            // "h(γᵢ)-granular ... submap of V(ᾱⱼ)" — O(num_chords +
            // simulated_contraction_weight per chord).
            assert(p.submap->is_granular(h_gamma, *p.curve) &&
                   "§3.0(ii)(3) tex 170: non-boundary piece's submap must "
                   "be h(γᵢ)-granular");
#endif
        } else {
            // (3) (tex 170): boundary pieces are "single-edge pieces
            // attached to the points of C corresponding to the endpoints
            // of α'".
            assert(p.subarc.first_edge == p.subarc.last_edge &&
                   "§3.0(ii)(3) tex 170: boundary piece must be a single-"
                   "edge piece (attached to an endpoint of α')");
        }
    }
}

} // namespace chazelle
