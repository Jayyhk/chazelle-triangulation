#pragma once

/// [C91 §3.1]: Fusion of two submaps.
///
/// "By symmetry, we may limit our discussion to the problem of
/// fusing S₁ into S₂, that is, determining the points of ∂C that
/// are seen by the endpoints of the exit chords of S₁ and by the
/// companion vertices resulting from the duplication of C₁ ∩ C₂."

#include "../polygon/polygon.h"
#include "../polygon/perturbation.h"
#include "../submap/submap.h"
#include "oracle.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

namespace chazelle {

/// [C91 §3.1]: A vertex in the fusion sequence.
///
/// The sequence is a₀, a₁, ..., aₘ, a_{m+1} where:
///   - a₀ = companion vertex at junction (start of clockwise tour)
///   - a₁, ..., aₘ = exit chord endpoints of S₁ in clockwise ∂C₁ order
///   - a_{m+1} = companion vertex at junction (end of clockwise tour)
struct FusionVertex {
    SymbolicY y;               ///< The vertex's symbolic y-coordinate.
    std::size_t edge;          ///< Edge of P containing this vertex.
    Side side;                 ///< Which side of ∂C₁.
    std::size_t chord_idx;     ///< Chord index in S₁ (NONE for companions).
    bool is_left_endpoint;     ///< True if this is the chord's LEFT endpoint.
    bool is_companion;         ///< True for a₀ or a_{m+1}.
};

/// [C91 §3.1]: Fixed-capacity result for collect_region_arcs.
///
/// §3.1 (tex 181): "the conformality of Sᵢ, which ensures that at most
/// four arcs need to be checked."  A degree-d region has exactly one arc
/// per chord-gap around its boundary, so degree ≤ 4 ⟹ at most 4 arcs.
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

/// [C91 §3.1]: Collect all arc indices belonging to a region.
///
/// Uses chord→arc adjacency + start_arc/end_arc.  O(1) for conformal
/// submaps (degree ≤ 4 → at most ~8 arc candidates).
///
/// @param S         The submap.
/// @param region    The region index.
RegionArcs collect_region_arcs(const Submap& S, std::size_t region);

/// [C91 §3.1]: Local shooting.
///
/// "Given any point p of ∂Cᵢ and the arc to which it belongs, we can
/// determine which point of ∂Cᵢ it sees (with respect to Cᵢ) in
/// O(f(γᵢ)) time."
///
/// "If p is an endpoint of an exit chord we can easily do that (even
/// in constant time).  If not, then p belongs to a unique region of
/// Sᵢ [...] and the point of ∂Cᵢ that it sees lies on one of the
/// region's arcs.  This is because regions are closed under visibility,
/// which is a corollary of Lemma 2.1.  Using the appropriate
/// ray-shooters, we can find that point by checking each arc in turn
/// and finding the nearest hit."
///
/// "Note that local shooting is still possible even if p does not lie
/// on ∂Cᵢ; it can lie anywhere in the spherical plane as long as a
/// horizontal direction (left or right) has been specified and we know
/// in which region of Sᵢ the point lies."
///
/// @param p          The point to shoot from.
/// @param direction  LEFT (shoot left) or RIGHT (shoot right).
/// @param region     The region of Sᵢ containing p.
/// @param S          The submap Sᵢ.
/// @param C          The curve Cᵢ.
/// @param oracle     The ray-shooting oracle for Sᵢ.
/// @return The hit point on ∂Cᵢ (always hits — regions are bounded).
RayHit local_shoot(Point p, Side direction,
                    std::size_t region,
                    const Submap& S, const Polygon& C,
                    const RayShootingOracle& oracle);

/// [C91 §3.1]: State maintained during the fusion traversal.
///
/// "To fuse S₁ into S₂ we let a variable p run through ∂C₁ in
/// clockwise order, stopping at a₀, ..., a_{m+1} as well as at some
/// other places to be specified.  We determine what p sees along the
/// way, while keeping track of the current region of S₂ in which p
/// lies."
struct FusionState {
    /// The fusion vertex sequence a₀, ..., a_{m+1}.
    std::vector<FusionVertex> sequence;

    /// Index into `sequence` of the current stop.
    std::size_t current_stop = 0;

    /// [C91 §3.1]: "a variable p" — current position on ∂C₁.
    Point p{0.0, 0.0, NONE};
    std::size_t p_edge = NONE;
    Side p_side = LEFT;

    /// Current region of S₂ that p lies in.
    std::size_t s2_region = NONE;

    /// [C91 §3.1]: O(1) array mapping to support context propagation.
    /// Maps a `cw_pos(arc_index)` to its start block in `sequence`.
    std::vector<std::size_t> arc_starts;

    /// Discovered chords (points of ∂C seen by fusion vertices).
    /// Populated during the traversal; used to build the fused submap.
    struct DiscoveredChord {
        SymbolicY y;           ///< y-coordinate of the chord.
        std::size_t left_edge; ///< Edge on LEFT side of ∂C.
        Side left_side;
        std::size_t right_edge; ///< Edge on RIGHT side of ∂C.
        Side right_side;
    };
    std::vector<DiscoveredChord> chords;
};

/// [C91 §3.1]: Determine the shooting direction for a point on ∂C.
///
/// "Because of the double boundary the shooting direction is always
/// uniquely defined."  A point on the LEFT side of an edge shoots
/// in the direction that crosses C; similarly for RIGHT.
Side shooting_direction(std::size_t edge, Side side,
                         const Polygon& C);

/// [C91 §3.1 Start-Up]: Initialize p and the current S₂ region.
///
/// Finds c₀ = the point of ∂C that a₀ sees, then:
///   Case 1: c₀ ∈ ∂C₂ → p = a₀, current = S₂ region crossed by a₀c₀.
///   Case 2: c₀ ∈ ∂C₁ → skip to c₀, p = c₀, current = S₂ region
///           containing a₀.
///
/// @return The index into state.sequence where the main loop should
///         start (0 for case 1, the index of c₀ for case 2).
std::size_t fusion_startup(FusionState& state,
                            const Submap& S1, const Polygon& C1,
                            const Submap& S2, const Polygon& C2,
                            const RayShootingOracle& oracle1,
                            const RayShootingOracle& oracle2);

/// [C91 §3.1]: Fuse S₁ into S₂ — traverse ∂C₁ clockwise, shooting
/// into S₂ at each fusion vertex to discover new chords.
///
/// TODO: (§3.1) implement main loop body.
void fuse_s1_into_s2(FusionState& state,
                      const Submap& S1, const Polygon& C1,
                      const Submap& S2, const Polygon& C2,
                      const RayShootingOracle& oracle1,
                      const RayShootingOracle& oracle2);

/// [C91 §3.1]: Build the fusion vertex sequence for fusing S₁ into S₂.
///
/// "Let a_{m+1} and a₀ be the companion vertices, as they appear next
/// to each other clockwise around ∂C₁, resulting from the duplication
/// of the vertex C₁ ∩ C₂ in ∂C₁.  Let a₁, a₂, ..., aₘ be the
/// canonical vertex enumeration of S₁."
///
/// "A clockwise tour around ∂C₁ that begins at a₀ thus ends at a_{m+1}."
///
/// The clockwise ∂C₁ order starting from the junction (c_end of C₁) is:
///   RIGHT side (c_end → c_start, descending edges), then
///   LEFT side (c_start → c_end, ascending edges).
///
/// @param state   The FusionState to populate with sequence and arc mappings.
/// @param S       The submap S₁.
/// @param C       The curve C₁.
void build_fusion_sequence(FusionState& state, const Submap& S,
                           const Polygon& C);

} // namespace chazelle
