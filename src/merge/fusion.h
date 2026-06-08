#pragma once

// [C91 §3.1]: Fuse S₁ into S₂ — "determining the points of ∂C that are
// seen by the endpoints of the exit chords of S₁ and by the companion
// vertices resulting from the duplication of C₁ ∩ C₂."  By symmetry,
// covers the reverse direction too.

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

// [C91 §3.1] (tex 181): Fixed-capacity result for collect_region_arcs —
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
RayHit local_shoot(Point p, Side direction,
                    std::size_t region,
                    const Submap& S, const Polygon& C,
                    const RayShootingOracle& oracle);

// [C91 §3.1]: Traversal state.  "We let a variable p run through ∂C₁
// in clockwise order, stopping at a₀, ..., a_{m+1} ... determining what
// p sees while keeping track of the current region of S₂."
struct FusionState {
    std::vector<FusionVertex> sequence;
    std::size_t current_stop = 0;

    Point p{0.0, 0.0, NONE};
    std::size_t p_edge = NONE;
    Side p_side = LEFT;

    // The region of S₂ that currently contains p.
    std::size_t s2_region = NONE;

    // Maps cw_pos(arc_index) → start block in `sequence`; O(1) lookup
    // for context propagation during the traversal.
    std::vector<std::size_t> arc_starts;

    // Chords discovered during the traversal — feed into the fused submap.
    struct DiscoveredChord {
        SymbolicY y;
        std::size_t left_edge;
        Side left_side;
        std::size_t right_edge;
        Side right_side;
    };
    std::vector<DiscoveredChord> chords;
};

// [C91 §3.1]: The unique shoot direction at a ∂C point ("because of
// the double boundary the shooting direction is always uniquely defined").
Side shooting_direction(std::size_t edge, Side side,
                         const Polygon& C);

// [C91 §3.1 Start-Up]: Initialize p and the current S₂ region.
//
// Compute c₀ = the point of ∂C that a₀ sees, then:
//   Case 1 (c₀ ∈ ∂C₂): p = a₀, current = S₂ region crossed by a₀c₀.
//   Case 2 (c₀ ∈ ∂C₁): skip to c₀, p = c₀, current = S₂ region of a₀.
//
// Returns the index into state.sequence at which the main loop starts
// (0 for Case 1, the index of c₀ for Case 2).
std::size_t fusion_startup(FusionState& state,
                            const Submap& S1, const Polygon& C1,
                            const Submap& S2, const Polygon& C2,
                            const RayShootingOracle& oracle1,
                            const RayShootingOracle& oracle2);

// [C91 §3.1]: Main loop — traverse ∂C₁ clockwise, shooting into S₂ at
// each stop to discover new chords.
// TODO: (§3.1) implement main loop body.
void fuse_s1_into_s2(FusionState& state,
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

} // namespace chazelle
