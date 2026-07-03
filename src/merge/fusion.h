#pragma once

// [C91 §3.1]: Fusion of two submaps — "determining the points of ∂C
// that are seen by the endpoints of the exit chords of S₁ and by the
// companion vertices resulting from the duplication of C₁ ∩ C₂."  The
// reverse direction (fuse S₂ into S₁) is the same algorithm with
// arguments swapped per tex 179 "by symmetry."

#include "../polygon/polygon.h"
#include "../polygon/perturbation.h"
#include "../submap/submap.h"
#include "oracle.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

namespace chazelle {

// [C91 §3.1]: A stop in the fusion sequence a₀, a₁, ..., aₘ, a_{m+1}.
//   a₀, a_{m+1}    — junction companion vertices (start/end of cw tour).
//   a₁, ..., aₘ    — S₁ exit chord endpoints in clockwise ∂C₁ order.
struct FusionVertex {
    SymbolicY y;
    std::size_t edge;
    Side side;
    std::size_t chord_idx;     // Chord index in S₁ (NONE for companions).
    bool is_left_endpoint;     // True if this is the chord's LEFT endpoint.
    bool is_companion;         // True for a₀ or a_{m+1}.
};

// [C91 §3.1 tex 181]: Fixed-capacity result for collect_region_arcs —
// conformality bounds the count: "at most four arcs need to be checked."
struct RegionArcs {
    static constexpr std::size_t MAX = 4;
    std::array<std::size_t, MAX> arcs = {};
    std::size_t count = 0;
    void push(std::size_t arc_idx) {
        assert(count < MAX);
        arcs[count++] = arc_idx;
    }
    const std::size_t* begin() const { return arcs.data(); }
    const std::size_t* end()   const { return arcs.data() + count; }
};

// [C91 §3.1]: All arc indices of a region.  O(1) for conformal submaps
// (degree ≤ 4 ⟹ at most ~8 arc candidates via chord adjacency).
RegionArcs collect_region_arcs(const Submap& S, std::size_t region);

// [C91 §3.1]: Determine which point of ∂Cᵢ point p sees in Sᵢ, in
// O(f(γᵢ)).  At an exit-chord endpoint this is O(1); else shoot into
// each of p's region's arcs (closed under visibility by Lemma 2.1) and
// pick the nearest hit.  p need not lie on ∂Cᵢ as long as `region` is
// the region of Sᵢ containing p.
//
// `require_hit = true` (default): asserts the ray hits (Lemma 2.1 — p
// is inside the region by precondition).  `require_hit = false`: returns
// `RayHit{hit=false}` when no arc is struck — needed by the `[C91 §3.1
// tex 220]` case (i) test, which uses no-hit as evidence that a_j ∉ R.
RayHit local_shoot(Point p, Side direction,
                    std::size_t region,
                    const Submap& S, const Polygon& C,
                    const RayShootingOracle& oracle,
                    bool require_hit = true);

// [C91 §3.1]: Traversal state.  "We let a variable p run through ∂C₁
// in clockwise order, stopping at a₀, ..., a_{m+1} ... determining what
// p sees while keeping track of the current region of S₂."
//
// [C91 §3.1 tex 179] "by symmetry" requires a second pass fusing S₂
// into S₁.  NOTE (limitation): the current implementation supports only
// walks whose junction (C₁ ∩ C₂) is the walked curve's LAST vertex —
// the S₁→S₂ direction.  The symmetric pass walks ∂C₂, whose junction is
// C₂'s FIRST vertex, and needs the mirrored tour (a₀ = LEFT companion
// at vertex 0, LEFT half first); build_fusion_sequence / fusion_startup
// / fuse_submaps do not implement that orientation yet.  This is one of
// the blockers recorded at merge.cpp::fuse.
struct FusionState {
    std::vector<FusionVertex> sequence;
    std::size_t current_stop = 0;

    Point p{0.0, 0.0, NONE};
    std::size_t p_edge = NONE;
    Side p_side = LEFT;
    // [C91 §2 tex 47]: SoS tag identifying p.y's perturbed source —
    // the polygon vertex or chord whose y-event placed p.  `Point.index`
    // is NONE for mid-edge p, so we carry the SoS tag separately to
    // keep downstream symbolic comparisons (cw_position, p ≠ a_l)
    // tie-break-correct.
    SymbolicY p_y{};

    // The region of S₂ that currently contains p.
    std::size_t s2_region = NONE;

    // Maps cw_pos(arc_index) → start block in `sequence`; O(1) lookup
    // for context propagation during the traversal.
    std::vector<std::size_t> arc_starts;

    // Chords discovered during the traversal — feed into the fused submap.
    // [C91 §3.1 tex 224]: discovered chords connect ∂C₁ and ∂C₂ — endpoints
    // live in DIFFERENT curve frames.  Flags are WALKER-frame: `true` means
    // the slot's `edge` indexes into the curve this pass WALKED (the first
    // curve argument of fuse_submaps), `false` the target curve.
    // rebuild_submap translates to the merged-C frame per pass.  (For the
    // case (ii) startup chord both endpoints can sit on the walked curve;
    // both flags would be true.)
    struct DiscoveredChord {
        SymbolicY y;
        std::size_t left_edge;
        Side left_side;
        std::size_t right_edge;
        Side right_side;
        bool left_on_walker  = true;
        bool right_on_walker = true;
    };
    std::vector<DiscoveredChord> chords;

    // [C91 §3.1 tex 224]: chord indices invalidated by a case (i)/(ii)
    // discovery in this pass — dropped by rebuild_submap's visibility
    // filter.  `invalidated_self` indexes the walked submap (S₁) and
    // `invalidated_other` indexes the target (S₂).
    std::vector<bool> invalidated_self;
    std::vector<bool> invalidated_other;
};

// [C91 §3.1]: The unique shoot direction at a ∂C point ("because of
// the double boundary the shooting direction is always uniquely defined").
Side shooting_direction(std::size_t edge, Side side,
                         const Polygon& C);

// [C91 §3.1 Start-Up]: Initialize p and the current S₂ region.
//
// Compute c₀ = the point of ∂C that a₀ sees, then:
//   case (i)  (c₀ ∈ ∂C₂): p = a₀, current = S₂ region crossed by a₀c₀.
//   case (ii) (c₀ ∈ ∂C₁): skip to c₀, p = c₀, current = S₂ region of a₀.
//
// Returns the index into state.sequence at which the main loop starts
// (0 for case (i), the index of c₀ for case (ii)).
std::size_t fusion_startup(FusionState& state,
                            const Submap& S1, const Polygon& C1,
                            const Submap& S2, const Polygon& C2,
                            const RayShootingOracle& oracle1,
                            const RayShootingOracle& oracle2);

// [C91 §3.1]: Fuse S₁ into S₂ — walk ∂C₁ clockwise, shoot into S₂ at
// each stop to discover the chords seen.  [C91 §3.1 Lemma 3.1]: each
// pass runs in O((n₁/γ₁ + n₂/γ₂ + 1)(f(γ₂) + log(n₁+n₂))).
//
// LIMITATION: only the S₁→S₂ direction (junction = walked curve's LAST
// vertex) is implemented; see the FusionState note above.  The
// [C91 §3.1 tex 179] symmetric pass (S₂→S₁, junction = walked curve's
// FIRST vertex) is pending — tracked at merge.cpp::fuse.
void fuse_submaps(FusionState& state,
                const Submap& S1, const Polygon& C1,
                const Submap& S2, const Polygon& C2,
                const RayShootingOracle& oracle1,
                const RayShootingOracle& oracle2);

// [C91 §3.1]: Build the fusion vertex sequence for S₁→S₂.
//
// a₀ and a_{m+1} are the two junction companions; a₁, ..., aₘ are S₁'s
// canonical vertex enumeration in clockwise ∂C₁ order.  The cw tour from
// the junction goes RIGHT side first (c_end→c_start, descending edges),
// then LEFT side (c_start→c_end, ascending edges).
void build_fusion_sequence(FusionState& state, const Submap& S,
                           const Polygon& C);

// [C91 §3.1 tex 226]: Set up the fused submap S in normal form from the
// chord inventory.  Sorts chord endpoints along ∂C, builds the arc-sequence
// table + region tree + chord-arc pointers.  Skips tree decomposition
// (§3.2 restores conformality first).
//
// Time: O((n₁/γ₁ + n₂/γ₂ + 1) log(n₁+n₂)) per [tex 226]: dominated by the
// endpoint sort; everything else is linear in the chord count.
//
// `state1.chords` and `state2.chords` carry the newly discovered chords
// (with per-endpoint `left_on_c1` / `right_on_c1` flags for curve-frame
// translation).  Old chords of S₁ and S₂ — filtered for V(C) visibility —
// and null-length chords are gathered from S₁ and S₂ here.
void rebuild_submap(Submap& out_S,
                     const Polygon& C,
                     const Submap& S1, const Polygon& C1,
                     const Submap& S2, const Polygon& C2,
                     const FusionState& state1,
                     const FusionState& state2);
} // namespace chazelle
