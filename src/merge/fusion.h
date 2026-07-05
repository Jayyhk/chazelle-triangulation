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

// [C91 §3.1 tex 181]: Fixed-capacity result for collect_region_arcs.
// Conformality bounds a region to ≤ 4 arcs ([C91 §2.3 tex 114]: degree
// ≤ 4; [C91 §2.2 tex 96]: the boundary alternates arcs and exit
// chords, so arc-structure count == degree) — a wrap-spanning arc is
// ONE arc-structure ([C91 §2.4 tex 142]), never split.
struct RegionArcs {
    static constexpr std::size_t MAX = 4;
    std::array<std::size_t, MAX> arcs = {};
    std::size_t count = 0;
    void push(std::size_t arc_idx) {
        assert(count < MAX &&
               "[C91 §2.3 tex 114]: conformal region has ≤ 4 arcs");
        arcs[count++] = arc_idx;
    }
    const std::size_t* begin() const { return arcs.data(); }
    const std::size_t* end()   const { return arcs.data() + count; }
};

// [C91 §3.1]: All arc indices of a region.  O(1) for conformal submaps
// (degree ≤ 4 ⟹ ≤ 16 slot candidates via chord adjacency).
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
// [C91 §3.1 tex 179] "by symmetry": the second pass fuses S₂ into S₁ by
// running the same algorithm with the walked curve's junction at its
// FIRST vertex instead of its LAST.  `junction_at_end` selects the tour:
//   true  (pass 1, walk C₁): a₀ = RIGHT companion at the last vertex;
//         cw tour = RIGHT half (descending edges), wrap at the start
//         vertex, LEFT half (ascending), ends at a_{m+1} = LEFT companion.
//   false (pass 2, walk C₂): a₀ = LEFT companion at vertex 0; cw tour =
//         LEFT half (ascending), wrap at the end vertex, RIGHT half
//         (descending), ends at a_{m+1} = RIGHT companion.
// Set it before build_fusion_sequence / fuse_submaps.
struct FusionState {
    std::vector<FusionVertex> sequence;
    std::size_t current_stop = 0;

    // [C91 §3.1 tex 179]: where C₁ ∩ C₂ sits in the WALKED curve's frame.
    bool junction_at_end = true;

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

// [C91 §2.1 tex 70/72]: does this non-null chord run through infinity?
// A chord wraps exactly when the chord direction at its endpoints
// points away from the other endpoint — the left (smaller-x) slot
// shoots LEFT — and a non-null chord whose endpoints share their x (an
// outside duplicate pair at a y-extremum) always wraps: the direct
// zero-length segment is realized by the null-length INSIDE pair, so
// the outside pair's visibility runs the other way around the sphere.
bool chord_runs_through_infinity(const Polygon& C, const Chord& c);

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
// [C91 §3.1 tex 179] "by symmetry": the reverse pass (fuse S₂ into S₁)
// is this same function with the arguments swapped AND
// state.junction_at_end = false (the junction is the walked curve's
// FIRST vertex there):
//   pass 1: st1.junction_at_end = true;  fuse_submaps(st1, S₁,C₁, S₂,C₂, o1,o2)
//   pass 2: st2.junction_at_end = false; fuse_submaps(st2, S₂,C₂, S₁,C₁, o2,o1)
void fuse_submaps(FusionState& state,
                const Submap& S1, const Polygon& C1,
                const Submap& S2, const Polygon& C2,
                const RayShootingOracle& oracle1,
                const RayShootingOracle& oracle2);

// [C91 §3.1]: Build the fusion vertex sequence for the walked submap.
//
// a₀ and a_{m+1} are the two junction companions; a₁, ..., aₘ are the
// walked submap's canonical vertex enumeration in clockwise ∂C order,
// starting and ending at the junction.  `state.junction_at_end` selects
// the tour orientation (see FusionState).
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
