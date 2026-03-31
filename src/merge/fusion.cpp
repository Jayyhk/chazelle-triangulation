/// src/merge/fusion.cpp

#include "fusion.h"

#include <algorithm>

namespace chazelle {

// ── collect_region_arcs ─────────────────────────────────────────

void collect_region_arcs(const Submap& S, std::size_t region,
                          std::vector<std::size_t>& out) {
    assert(region < S.num_nodes() && !S.node(region).dead);

    // Collect via chord→arc adjacency (same traversal as region_weight).
    auto check_adj = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            std::size_t ai = adj.arcs[k];
            assert(ai < S.num_arcs() && !S.arc(ai).dead);
            if (S.arc(ai).region_node == region) {
                // Avoid duplicates: check if already added.
                bool dup = false;
                for (std::size_t x : out)
                    if (x == ai) { dup = true; break; }
                if (!dup) out.push_back(ai);
            }
        }
    };

    for (std::size_t ci : S.node(region).incident_chords) {
        assert(ci < S.num_chords());
        if (S.chord(ci).dead) continue;
        check_adj(S.chord(ci).left_adj);
        check_adj(S.chord(ci).right_adj);
    }

    // Also check start_arc / end_arc.
    if (S.start_arc != NONE && S.start_arc < S.num_arcs() &&
        !S.arc(S.start_arc).dead &&
        S.arc(S.start_arc).region_node == region) {
        bool dup = false;
        for (std::size_t x : out)
            if (x == S.start_arc) { dup = true; break; }
        if (!dup) out.push_back(S.start_arc);
    }
    if (S.end_arc != NONE && S.end_arc < S.num_arcs() &&
        !S.arc(S.end_arc).dead &&
        S.arc(S.end_arc).region_node == region) {
        bool dup = false;
        for (std::size_t x : out)
            if (x == S.end_arc) { dup = true; break; }
        if (!dup) out.push_back(S.end_arc);
    }
}

// ── local_shoot ─────────────────────────────────────────────────

RayHit local_shoot(Point p, Side direction,
                    std::size_t region,
                    const Submap& S, const Polygon& /*C*/,
                    const RayShootingOracle& oracle) {
    assert(S.is_conformal() &&
           "§3.1: local shooting requires conformal submap "
           "(at most 4 arcs per region)");

    // [C91 §3.1]: "Using the appropriate ray-shooters, we can find
    // that point by checking each arc in turn and finding the nearest
    // hit.  The claim on the time follows from the conformality of Sᵢ,
    // which ensures that at most four arcs need to be checked."
    std::vector<std::size_t> arcs;
    collect_region_arcs(S, region, arcs);

    RayHit best;
    best.hit = false;

    for (std::size_t ai : arcs) {
        const auto& a = S.arc(ai);

        // Build a Subarc covering the full arc.
        Subarc sub;
        sub.first_edge = a.first_edge;
        sub.first_side = a.first_side;
        sub.last_edge = a.last_edge;
        sub.last_side = a.last_side;

        RayHit hit = oracle.shoot(p, direction, sub);
        if (!hit.hit) continue;

        if (!best.hit) {
            best = hit;
        } else {
            // Nearest in the shooting direction.
            // LEFT: largest x < p.x → closest to p.
            // RIGHT: smallest x > p.x → closest to p.
            double best_dist = (direction == LEFT)
                ? (p.x - best.x) : (best.x - p.x);
            double hit_dist = (direction == LEFT)
                ? (p.x - hit.x) : (hit.x - p.x);
            if (hit_dist < best_dist)
                best = hit;
        }
    }

    return best;
}

// ── shooting_direction ──────────────────────────────────────────

Side shooting_direction(std::size_t edge, Side side,
                         const Polygon& C) {
    // [C91 §3.1]: "because of the double boundary the shooting
    // direction is always uniquely defined."
    //
    // A point on the LEFT side of edge e sees across C to the RIGHT.
    // The horizontal direction depends on the edge's slope:
    //   Edge going up (start.y < end.y): LEFT side is to the left
    //     of the upward direction → shoot RIGHT.
    //   Edge going down: LEFT side is to the right → shoot LEFT.
    // For RIGHT side, the direction is reversed.
    assert(edge < C.num_edges());
    const auto& e = C.edge(edge);
    SymbolicY start_y = symbolic_y_of(C.vertex(e.start_idx));
    SymbolicY end_y = symbolic_y_of(C.vertex(e.end_idx));
    bool edge_ascending = symbolic_y_less(start_y, end_y);

    if (side == LEFT)
        return edge_ascending ? RIGHT : LEFT;
    else
        return edge_ascending ? LEFT : RIGHT;
}

// ── fusion_startup ──────────────────────────────────────────────

std::size_t fusion_startup(FusionState& state,
                            const Submap& S1, const Polygon& C1,
                            const Submap& S2, const Polygon& C2,
                            const RayShootingOracle& oracle1,
                            const RayShootingOracle& oracle2) {
    assert(!state.sequence.empty());
    const auto& a0 = state.sequence[0];
    assert(a0.is_companion && a0.side == RIGHT);

    // [C91 §3.1]: "Using local shooting, we find the point of ∂C₁
    // that a₀ sees with respect to C₁."
    //
    // a₀ is at c_end of C₁, RIGHT side.  The region in S₁ is the
    // one containing the end_arc (or its RIGHT-side neighbor).
    std::size_t a0_region_s1 = S1.arc(S1.end_arc).region_node;
    Side a0_dir = shooting_direction(a0.edge, a0.side, C1);
    Point a0_point{0.0, a0.y.y, a0.y.tag}; // x doesn't matter for horizontal
    // Compute x from the edge geometry at the junction.
    {
        const auto& e = C1.edge(a0.edge);
        double y0 = C1.vertex(e.start_idx).y;
        double y1 = C1.vertex(e.end_idx).y;
        double x0 = C1.vertex(e.start_idx).x;
        double x1 = C1.vertex(e.end_idx).x;
        if (y0 != y1) {
            double t = (a0.y.y - y0) / (y1 - y0);
            a0_point.x = x0 + t * (x1 - x0);
        } else {
            a0_point.x = (x0 + x1) / 2.0;
        }
    }

    RayHit hit_c1 = local_shoot(a0_point, a0_dir, a0_region_s1,
                                 S1, C1, oracle1);

    // [C91 §3.1]: "using the information about the endpoints of C₂
    // encoded in the normal-form representation of S₂ (namely,
    // pointers to incident arcs), we can find, in constant time, in
    // which region of S₂ the point a₀ lies."
    //
    // a₀ corresponds to C₂'s start vertex → S₂.start_arc's region.
    assert(S2.start_arc != NONE);
    std::size_t a0_region_s2 = S2.arc(S2.start_arc).region_node;

    // [C91 §3.1]: "This allows us to do local shooting and find the
    // point of ∂C₂ that a₀ sees with respect to C₂."
    RayHit hit_c2 = local_shoot(a0_point, a0_dir, a0_region_s2,
                                 S2, C2, oracle2);

    // [C91 §3.1]: "These two pieces of information combine to give us
    // the unique point c₀ of ∂C that a₀ sees with respect to C."
    //
    // c₀ is the nearer of hit_c1 and hit_c2 in the shooting direction.
    bool c0_on_c2;
    RayHit c0;
    if (!hit_c1.hit && !hit_c2.hit) {
        // No hit — shouldn't happen in a valid visibility map.
        assert(false && "§3.1: a₀ must see some point on ∂C");
    } else if (!hit_c1.hit) {
        c0 = hit_c2; c0_on_c2 = true;
    } else if (!hit_c2.hit) {
        c0 = hit_c1; c0_on_c2 = false;
    } else {
        double d1 = (a0_dir == LEFT)
            ? (a0_point.x - hit_c1.x) : (hit_c1.x - a0_point.x);
        double d2 = (a0_dir == LEFT)
            ? (a0_point.x - hit_c2.x) : (hit_c2.x - a0_point.x);
        if (d2 <= d1) {
            c0 = hit_c2; c0_on_c2 = true;
        } else {
            c0 = hit_c1; c0_on_c2 = false;
        }
    }

    if (c0_on_c2) {
        // [C91 §3.1 Case 1]: "If c₀ belongs to ∂C₂, then we set
        // p = a₀ and we call the region of S₂ crossed by a₀c₀
        // current: the start-up phase is over."
        state.s2_region = a0_region_s2;
        state.chords.push_back({a0.y, a0.edge, a0.side,
                                c0.edge, c0.side});
        return 0;
    } else {
        // [C91 §3.1 Case 2]: "If c₀ belongs to ∂C₁... We can
        // therefore skip all the way to c₀.  Now, however, c₀ sees
        // a point of ∂C₂, namely a₀, so we set p = c₀ and call the
        // region of S₂ containing a₀ current."
        state.s2_region = a0_region_s2;
        state.chords.push_back({a0.y, c0.edge, c0.side,
                                a0.edge, a0.side});

        // Find the fusion sequence index closest to c₀ to skip to.
        // All vertices between a₀ and c₀ (clockwise) can be skipped
        // because "the points of ∂C that its exit chord endpoints see
        // all belong to ∂C₁ and thus are available directly from S₁."
        //
        // TODO: (§3.1) precise skip index computation.
        // For now, return 0 (conservative — don't skip).
        return 0;
    }
}

// ── fuse_s1_into_s2 ─────────────────────────────────────────────

void fuse_s1_into_s2(FusionState& state,
                      const Submap& S1, const Polygon& C1,
                      const Submap& S2, const Polygon& C2,
                      const RayShootingOracle& oracle1,
                      const RayShootingOracle& oracle2) {
    // [C91 §3.1]: Build the fusion vertex sequence.
    state.sequence = build_fusion_sequence(S1, C1);
    state.current_stop = 0;

    // [C91 §3.1 Start-Up]: Initialize p and current S₂ region.
    std::size_t start_idx = fusion_startup(state, S1, C1, S2, C2,
                                            oracle1, oracle2);

    // [C91 §3.1]: "We let a variable p run through ∂C₁ in clockwise
    // order, stopping at a₀, ..., a_{m+1}."
    for (std::size_t i = start_idx; i < state.sequence.size(); ++i) {
        state.current_stop = i;

        // TODO: (§3.1) At each stop:
        //   1. Determine what p sees (local_shoot into S₂).
        //   2. Record discovered chord.
        //   3. Update s2_region for the next stop.
    }
}

// ── build_fusion_sequence ───────────────────────────────────────

std::vector<FusionVertex> build_fusion_sequence(const Submap& S,
                                                 const Polygon& C) {
    std::size_t n = C.num_vertices();
    assert(n >= 2);

    // [C91 §3.1]: The junction vertex is C₁ ∩ C₂ = last vertex of C₁.
    std::size_t junction_v = n - 1;
    std::size_t junction_edge = junction_v - 1; // last edge of C₁

    // [C91 §3.1]: "Let a_{m+1} and a₀ be the companion vertices,
    // as they appear next to each other clockwise around ∂C₁,
    // resulting from the duplication of the vertex C₁ ∩ C₂."
    //
    // Clockwise ∂C₁: LEFT ascending → junction → RIGHT descending.
    // a_{m+1} = LEFT companion at junction (end of LEFT traversal).
    // a₀ = RIGHT companion at junction (start of RIGHT traversal).
    //
    // Both have the junction vertex's y-coordinate and SoS tag.
    SymbolicY junction_y = symbolic_y_of(C.vertex(junction_v));

    FusionVertex a_0;
    a_0.y = junction_y;
    a_0.edge = junction_edge;
    a_0.side = RIGHT;  // RIGHT companion: start of clockwise tour
    a_0.chord_idx = NONE;
    a_0.is_left_endpoint = false;
    a_0.is_companion = true;

    FusionVertex a_m1;
    a_m1.y = junction_y;
    a_m1.edge = junction_edge;
    a_m1.side = LEFT;   // LEFT companion: end of clockwise tour
    a_m1.chord_idx = NONE;
    a_m1.is_left_endpoint = false;
    a_m1.is_companion = true;

    // [C91 §3.1]: "Let a₁, a₂, ..., aₘ be the canonical vertex
    // enumeration of S₁.  Recall that this enumerates the exit chord
    // endpoints in S₁ as we encounter them going clockwise around ∂C₁."
    //
    // Each chord has two endpoints (LEFT and RIGHT).  Collect them all.
    std::vector<FusionVertex> seq;
    seq.reserve(2 * S.num_live_chords() + 2);

    for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
        if (S.chord(ci).dead) continue;
        const auto& c = S.chord(ci);

        FusionVertex left_ep;
        left_ep.y = c.symbolic_y();
        left_ep.edge = c.left_edge;
        left_ep.side = c.left_side;
        left_ep.chord_idx = ci;
        left_ep.is_left_endpoint = true;
        left_ep.is_companion = false;
        seq.push_back(left_ep);

        FusionVertex right_ep;
        right_ep.y = c.symbolic_y();
        right_ep.edge = c.right_edge;
        right_ep.side = c.right_side;
        right_ep.chord_idx = ci;
        right_ep.is_left_endpoint = false;
        right_ep.is_companion = false;
        seq.push_back(right_ep);
    }

    // Sort in clockwise ∂C₁ order starting from a₀.
    //
    // Clockwise ∂C₁ from a₀:
    //   RIGHT side: junction_edge descending to edge 0
    //   LEFT side: edge 0 ascending to junction_edge
    //
    // Assign a traversal position for sorting:
    //   RIGHT vertices: pos = (junction_edge - edge)  [0 at junction, increasing toward start]
    //   LEFT vertices: pos = (junction_edge + 1) + edge  [continues after RIGHT]
    //
    // Within the same edge, sort by y along the ∂C traversal direction.
    std::size_t right_half_len = junction_edge + 1; // edges 0..junction_edge

    auto traversal_pos = [&](const FusionVertex& v) -> std::size_t {
        if (v.side == RIGHT) {
            return (junction_edge - v.edge);
        } else {
            return right_half_len + v.edge;
        }
    };

    std::sort(seq.begin(), seq.end(),
        [&](const FusionVertex& a, const FusionVertex& b) {
            std::size_t pa = traversal_pos(a);
            std::size_t pb = traversal_pos(b);
            if (pa != pb) return pa < pb;
            // Same edge: sort by y along ∂C traversal direction.
            // RIGHT side descends in edge index but the y-order within
            // an edge depends on the edge's geometric direction.
            // LEFT side: ∂C traversal follows edge direction (ascending).
            // RIGHT side: ∂C traversal is reversed.
            // For vertices on the same edge, use symbolic y:
            // on LEFT, earlier in traversal = lower y on ascending edges
            // (but depends on edge direction — just use y for now,
            // ties broken by SoS).
            return symbolic_y_less(a.y, b.y);
        });

    // [C91 §3.1]: "we assume that [a₀ and a_{m+1}] are not [in the
    // sequence] and therefore add them to the sequence."
    // Prepend a₀, append a_{m+1}.
    std::vector<FusionVertex> result;
    result.reserve(seq.size() + 2);
    result.push_back(a_0);
    result.insert(result.end(), seq.begin(), seq.end());
    result.push_back(a_m1);

    return result;
}

} // namespace chazelle
