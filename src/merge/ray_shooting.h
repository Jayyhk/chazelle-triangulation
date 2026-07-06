#pragma once

// [C91 §3.4 tex 282–312]: Implementing the Oracles — the ray-shooting
// structure of Lemma 3.6.
//
// "Let C be a connected vertex-to-vertex piece of the input polygonal
// curve P and let m be its number of vertices.  Let S be a normal-form
// γ-granular conformal submap of V(C).  Then it is possible to
// preprocess S in O(m(log m)/γ + 1) time so that ray-shooting within S
// can be done in time O(γ^{1/3} m^{2/3})."  [C91 Lemma 3.6 tex 310–312]
//
// Pieces ([C91 §3.4]):
//   S*  — S as a planar subdivision with the double boundary collapsed;
//         faces ↔ nonempty regions (tex 286).
//   G   — the dual graph of S*: nodes = faces, edges = shared S*-edges;
//         planar, ≤ 3μ−6 edges (tex 295–297).
//   D*  — recursive Lipton–Tarjan separators until pieces ≤ μ^{2/3},
//         plus the D_i partition (tex 297–304; planar_separator.h).
//   The vertical line right of all vertices, split by the chords that
//   cross it into segments identified with regions (tex 306).
//   Query — shoot within each D* region naively, identify the last
//   region R before the first hit by double identification, then scan
//   R's D_i; if no D* region is hit, locate the D_i through the
//   vertical line (tex 306–308).
//
// The arc-cutting oracle is NOT implemented here: "The arc-cutter is
// implemented by using the divide-and-conquer structure of the up-phase
// of the visibility algorithm ... we postpone the discussion of its
// implementation" ([C91 §3.4 tex 284]; owned by §4).

#include "../polygon/point.h"
#include "../polygon/perturbation.h"
#include "../polygon/polygon.h"
#include "../submap/submap.h"
#include "oracle.h"
#include "planar_separator.h"

#include <cstddef>
#include <vector>

namespace chazelle {

// [C91 §3.4 Lemma 3.6]: preprocessing of one normal-form γ-granular
// conformal submap S of V(C) for ray-shooting toward ∂C.
class RayShootingStructure {
public:
    // Preconditions (asserted): S in normal form ([C91 §3 tex 166]),
    // conformal, γ-semigranular (weights ≤ γ; [C91 §2.3]).
    // Time: O(μ log m) for μ > 1 ([C91 §3.4 tex 297]), trivial for μ = 1.
    RayShootingStructure(const Submap& S, const Polygon& C,
                         std::size_t gamma);

    // [C91 §3.4 tex 306–308]: shoot the horizontal ray from p toward ∂C
    // and report the first contact (edge of C + struck ∂C side + point),
    // in the wrap metric of [C91 §2.1 tex 70] (direct contacts precede
    // through-infinity ones; RayHit.wrapped marks the latter).  No
    // contact at all (RayHit.hit = false) occurs exactly when p's
    // perturbed y lies outside C's y-range.  hit_arc_idx is NOT set
    // (oracle.h contract).
    // Time: O(m) if μ = 1, else O(γ μ^{2/3}) (tex 297/308).
    RayHit shoot_toward_boundary(Point p, Side direction) const;

    // ── Introspection (tests) ─────────────────────────────────────
    std::size_t mu() const noexcept { return mu_; }
    // region → face id, NONE for empty regions ([C91 §3.4 tex 286]:
    // "empty regions have no associated faces").
    std::size_t face_of_region(std::size_t r) const {
        return face_of_region_[r];
    }
    std::size_t region_of_face(std::size_t f) const {
        return region_of_face_[f];
    }
    // Deduplicated dual-graph edges as face-id pairs.
    const std::vector<std::pair<std::size_t, std::size_t>>& dual_edges()
        const noexcept { return dual_edges_; }
    const SeparatorDecomposition& decomposition() const noexcept {
        return decomp_;
    }
    // [C91 §3.4 tex 306]: the vertical line's crossings, bottom to top;
    // the crossing chords are exactly the through-infinity chords of S.
    struct LineCrossing {
        std::size_t chord = NONE;
        SymbolicY y{};
        std::size_t region_below = NONE;
        std::size_t region_above = NONE;
    };
    const std::vector<LineCrossing>& vertical_line() const noexcept {
        return line_;
    }
    // Region owning the whole vertical line when no chord crosses it.
    std::size_t region_at_infinity() const noexcept {
        return region_infinity_;
    }

private:
    const Submap* S_;
    const Polygon* C_;
    std::size_t gamma_;

    std::size_t mu_ = 0;
    std::vector<std::size_t> face_of_region_;
    std::vector<std::size_t> region_of_face_;
    std::vector<std::vector<std::size_t>> arcs_of_region_;
    std::vector<std::pair<std::size_t, std::size_t>> dual_edges_;
    SeparatorDecomposition decomp_;
    std::vector<std::size_t> dstar_faces_;
    // [C91 §3.4 tex 308]: member faces of each D_i, built once at
    // preprocessing so a query scans |D_i| ≤ μ^{2/3} faces, not μ.
    std::vector<std::vector<std::size_t>> subset_faces_;
    std::vector<LineCrossing> line_;
    std::size_t region_infinity_ = NONE;

    // One side of ∂C as sorted positive-length boundary intervals (the
    // same lists that drive the dual-graph merge, [C91 §3.4 tex 295]);
    // retained for the query's double identification of the first hit
    // ([C91 §3.4 tex 306–308]: "To identify R can be done by double
    // identification ... might require a binary search").
    struct BoundaryInterval {
        std::size_t lo_edge = NONE, hi_edge = NONE;
        SymbolicY lo_y{}, hi_y{};
        std::size_t region = NONE;
    };
    std::vector<BoundaryInterval> left_ivals_, right_ivals_;

    // Regions whose boundary covers ∂C position (edge, side, y):
    // one region at interval interiors, both flanking regions at
    // interval junctions ([C91 §2 tex 47] point coincidences).
    void regions_at_boundary(std::size_t edge, Side side, SymbolicY y,
                             std::vector<std::size_t>& out) const;

    void build_faces();
    void build_dual_graph_and_decomposition();
    void build_vertical_line();
};

// [C91 §3.0(i) tex 169]-shaped oracle over the [C91 Lemma 3.6]
// structure: shoot(p, dir, α, α') reports the first contact of the ray
// with ∂C filtered to the target subarc α'.
//
// Semantics note.  [C91 §3.0(i) tex 169] specifies the report "in the
// absence of any obstacle except α'"; the paper realizes that contract
// in the up-phase by decomposing α' into pieces and shooting in each
// piece's own Lemma 3.6 structure, whose curve IS the piece
// ([C91 §4.1 tex 341]: "to shoot a ray toward α', we shoot toward each
// of the O(λ) subarcs of its decomposition and determine the closest
// hit").  At the §3.4 stage the available structure is the one over C
// itself, so the obstacle is C ⊇ ᾱ'.  The two reports agree whenever
// the tex-169 hit is visible with respect to C, and every §3.1/§3.2
// consumer discards non-C-visible candidates anyway:
//  · local shooting ([C91 §3.1 tex 181]) takes the nearest hit over the
//    arcs of the region containing p, which IS the first C-contact;
//  · the case (i) test ([C91 §3.1 tex 220]) infers a_j ∈ R from the
//    local orientation at the first hit — under obstacle-C semantics a
//    report on an arc of R exists iff a_j lies in R (a horizontal ray
//    can enter a region's interior only through its arcs, never
//    through its horizontal chords);
//  · the case (ii) test ([C91 §3.1 tex 222]) confirms every candidate
//    with the back-shot mutual-visibility test, which rejects exactly
//    the pairs blocked by C.
class SubmapRayShooter final : public RayShootingOracle {
public:
    SubmapRayShooter(const Submap& S, const Polygon& C, std::size_t gamma)
        : S_(&S), C_(&C), impl_(S, C, gamma) {}

    RayHit shoot(Point p, Side direction, std::size_t arc_idx,
                 const Subarc& target) const override;

    const RayShootingStructure& structure() const noexcept { return impl_; }

private:
    const Submap* S_;
    const Polygon* C_;
    RayShootingStructure impl_;
};

} // namespace chazelle
