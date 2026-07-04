#pragma once

// [C91 §3.2 tex 236–272]: Restoring Conformality.
//
// After fusion, "there is no reason to believe that the fusion S should
// be conformal" (tex 238), but no region has more than a bounded number
// of arcs: the arcs of any region partition into at most TWO runs (one
// per ∂Cᵢ, by Lemma 2.2), and a run "cannot have more than four exit
// chords in its midst, not counting the new chords incident upon a₀ or
// a_{m+1}" (tex 238).  This module reduces every region to at most four
// arcs by adding visibility chords (Lemma 3.2 finds them; Lemma 3.3
// proves they exist; Lemma 3.4 bounds the total cost).

#include "../polygon/polygon.h"
#include "../polygon/perturbation.h"
#include "../submap/submap.h"
#include "oracle.h"

#include <cstddef>
#include <vector>

namespace chazelle {

// ── [C91 §3.2 tex 244]: Arc provenance ──────────────────────────
//
// "A₁ and A₂ are arcs or subarcs of either S₁ or S₂ but cannot overlap
// both ∂C₁ and ∂C₂.  The reason is that all chord endpoints in S₁ and
// S₂ are still chord endpoints in S ... and that we added chords
// incident upon the vertices of ∂C resulting from the duplication of
// the vertex C₁ ∩ C₂."  Every arc of the fused S therefore lies inside
// exactly one arc of S₁ or S₂; the paper's procedures ("the arc-cutter
// associated with the arc of S₁ or S₂ containing A₁", tex 248, and
// local shooting within S, tex 244) require that containing arc.

struct ArcProvenance {
    bool on_c1 = true;              // arc lies on ∂C₁ (else ∂C₂)
    std::size_t arc_in_si = NONE;   // containing arc's index in Sᵢ
};

// Compute provenance for every live arc of S.  One double_identify per
// arc — O(M log(n₁+n₂)) total, within the [C91 §3.1 tex 226] rebuild
// budget (computed once per merge, before any §3.2 mutation).
std::vector<ArcProvenance> compute_arc_provenance(
    const Submap& S, const Polygon& C,
    const Submap& S1, const Polygon& C1,
    const Submap& S2, const Polygon& C2);

// ── [C91 §3.2 tex 238]: Fused region cycles ─────────────────────
//
// The paper's arcs are maximal ∂C pieces between chord endpoints; the
// fused table additionally splits them at C's two endpoint wraps
// (rebuild_submap emits single-side arcs only).  A LogicalArc re-glues
// wrap-adjacent table pieces so that the cycle counts PAPER arcs —
// Lemma 3.3's "k > 4" and the conformality target "four arcs or less"
// (tex 264) are stated for paper arcs.

struct LogicalArc {
    static constexpr std::size_t MAX_PIECES = 3;  // both wraps + between
    std::size_t pieces[MAX_PIECES] = {NONE, NONE, NONE};
    std::size_t piece_count = 0;
    // [C91 §2.2 tex 96] "some arcs may be of zero length": a zero-length
    // arc is a single ∂C point.  Its visibility is settled — it is (or
    // duplicates) a chord endpoint, whose horizontal visibility is
    // unique ([C91 §2.1 tex 70/72]) — so it is never a Lemma 3.2
    // candidate, though it still counts toward k.
    bool is_zero_length = false;
};

struct FusedRegionCycle {
    // [C91 §3.2 tex 238]: ≤ 2 runs × (4 midst exit chords + ≤ 4 new
    // junction chords + 1 junction null chord + 1) = 20 paper arcs;
    // generous headroom for the ≤ 2 extra wrap pieces at table level.
    static constexpr std::size_t MAX_ARCS = 24;
    LogicalArc arcs[MAX_ARCS];
    std::size_t count = 0;
};

// Build the region's boundary cycle from its (unordered) live arc list.
// [C91 §2.2 tex 96] Lemma 2.2: sorting by clockwise ∂C position yields
// the region's boundary order.  Wrap-adjacent table pieces (no chord
// endpoint between them) are merged into one logical arc.
FusedRegionCycle fused_region_cycle(const Submap& S, const Polygon& C,
                                    std::size_t region,
                                    const std::vector<std::size_t>& arcs_of_region);

// ── [C91 §3.2 tex 244/246]: Local shooting within fused S ───────
//
// "Because of the bounded number of arcs per region, it is still
// possible to do local shooting within any region of S ... this takes
// O(f(γ₂)) time" (tex 244).  Each region arc is shot via its provenance
// (Sᵢ arc + the fused arc's span as the target subarc), and each hit is
// validated against the arc's true extent and ∂C side — the paper's
// "local checking" for the companion ambiguity (tex 246: "local
// shooting reports edges of P and does not tell us if the point hit is
// on the desired arc or is the companion").
//
// On return, RayHit.edge is in C's frame and RayHit.hit_arc_idx is the
// struck arc's index in S's table.

struct FusedShootContext {
    const Submap* S = nullptr;
    const Polygon* C = nullptr;
    const Polygon* C1 = nullptr;
    const Polygon* C2 = nullptr;
    const RayShootingOracle* ray1 = nullptr;
    const RayShootingOracle* ray2 = nullptr;
    const std::vector<ArcProvenance>* provenance = nullptr;
};

RayHit local_shoot_fused(Point p, SymbolicY p_y, Side direction,
                         const FusedRegionCycle& cycle,
                         const FusedShootContext& ctx,
                         bool require_hit = true);

// ── [C91 §3.2 Lemma 3.2]: Finding a visible point ───────────────

struct VisiblePoint {
    bool found = false;
    // The point of A₁ ([C91 §3.2 Lemma 3.2]: "a point of A₁ (not
    // necessarily a vertex of ∂C) that sees A₂").
    std::size_t p_table_arc = NONE;  // table arc of S containing p
    std::size_t p_edge = NONE;       // C-frame edge
    Side p_side = LEFT;
    double p_x = 0.0;
    SymbolicY y{};                   // the chord's symbolic y
    // The point of A₂ that p sees (the local-shoot hit).
    std::size_t q_table_arc = NONE;
    std::size_t q_edge = NONE;       // C-frame edge
    Side q_side = LEFT;
    double q_x = 0.0;
};

struct ConformalityOracles {
    const Submap* S1 = nullptr;
    const Submap* S2 = nullptr;
    const Polygon* C1 = nullptr;
    const Polygon* C2 = nullptr;
    const RayShootingOracle* ray1 = nullptr;
    const RayShootingOracle* ray2 = nullptr;
    const ArcCuttingOracle* cut1 = nullptr;
    const ArcCuttingOracle* cut2 = nullptr;
    // [C91 §3.0(ii) tex 170]: cut() bounds — ≤ g(γᵢ) pieces, each with an
    // h(γᵢ)-granular conformal submap; enforced by
    // assert_cut_postconditions at every §3.2 cut() call.
    std::size_t g_gamma1 = 0;
    std::size_t g_gamma2 = 0;
    std::size_t h_gamma1 = 0;
    std::size_t h_gamma2 = 0;
};

// [C91 §3.2 Lemma 3.2]: if A₁ has a vertex of C that sees A₂ (with
// respect to C), find a point of A₁ that sees A₂ in
// O(f(γ₂)g(γ₂)(h(γ₂) + log γ₂)) time.  Sound also when the premise
// fails: success is only reported after an actual local shoot lands on
// A₂, so a not-found result simply means this ordered pair contributes
// no chord.
VisiblePoint find_visible_point(const Submap& S, const Polygon& C,
                                std::size_t region,
                                const LogicalArc& A1, const LogicalArc& A2,
                                const FusedRegionCycle& cycle,
                                const std::vector<ArcProvenance>& provenance,
                                const ConformalityOracles& oracles);

// ── [C91 §3.2 Lemma 3.4]: Driver ────────────────────────────────
//
// "For any region with more than four arcs, let us apply Lemma 3.2 to
// every pair of nonconsecutive arcs until we find a chord which we can
// add to S.  We iterate on this process until no region has more than
// four arcs."  (tex 264)  Lemma 3.3 guarantees each iteration finds a
// chord; each addition splits a k-arc region into two regions of k₁+1
// and k₂+1 arcs with k₁+k₂ = k and k₁,k₂ ≥ 2 (nonconsecutive split), so
// Σ max(0, arcs−4) drops by ≥ 2 per chord and the total number of
// additions is O(n₁/γ₁ + n₂/γ₂ + 1) — Lemma 3.4's bound
// O((n₁/γ₁ + n₂/γ₂ + 1)·f(γ₂)g(γ₂)(h(γ₂) + log γ₂)).
//
// Postcondition: every region of S has at most four arcs and S is
// conformal ([C91 §2.3 tex 114]).  S is NOT left in normal form — the
// appended arc halves break canonical table order; [C91 §3.3 tex 276]
// re-normalizes.
void restore_conformality(Submap& S, const Polygon& C,
                          const ConformalityOracles& oracles);

} // namespace chazelle
