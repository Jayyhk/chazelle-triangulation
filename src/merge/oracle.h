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

#include <cstddef>
#include <vector>

namespace chazelle {

/// [C91 §3.0(i)(ii)]: A subarc specification — a contiguous piece of
/// a region arc α of Sᵢ, "specified by its two endpoints (along with
/// two pointers to the input table to indicate the names of the edges
/// of P that contain these two endpoints as well as two flags
/// indicating which side of ∂P is to be understood)."
///
/// §3.0(ii) condition (1): "the pair should be ordered to reflect a
/// clockwise traversal along ∂P."  That is, first_edge/first_side
/// comes before last_edge/last_side in the clockwise ∂C order.
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
    bool hit = false;          ///< True if the ray hit the subarc.
    double x = 0.0;            ///< x-coordinate of the hit point.
    double y = 0.0;            ///< y-coordinate of the hit point (= p.y).
    std::size_t edge = 0;      ///< "the name of the edge of P that contains it."
    Side side = LEFT;          ///< Which side of ∂P the hit is on.
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
    /// @param arc_idx Index of region arc α in Sᵢ's arc-sequence
    ///                table — "specified by a pointer to its
    ///                arc-structure" (§3.0, tex 166).
    /// @param target  The subarc α' to cut.
    /// @return Pieces α₁, α₂, ..., ordered along ∂C.
    virtual std::vector<ArcPiece> cut(std::size_t arc_idx,
                                       const Subarc& target) const = 0;
};

} // namespace chazelle
