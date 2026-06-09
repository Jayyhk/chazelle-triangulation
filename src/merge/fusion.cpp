// src/merge/fusion.cpp

#include "fusion.h"

namespace chazelle {

// ── collect_region_arcs ─────────────────────────────────────────

RegionArcs collect_region_arcs(const Submap& S, std::size_t region) {
    assert(region < S.num_nodes() && !S.node(region).dead);

    RegionArcs out;
    // [C91 §3.1 tex 181]: conformality bounds the iteration to O(1).
    assert(S.node(region).degree() <= 4 &&
           "[C91 §2.3]: conformal regions MUST have degree ≤ 4");
    auto check_adj = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            std::size_t ai = adj.arcs[k];
            assert(ai < S.num_arcs() && !S.arc(ai).dead);
            if (S.arc(ai).region_node == region) {
                // Avoid duplicates: O(count) scan, count ≤ MAX = O(1).
                bool dup = false;
                for (std::size_t i = 0; i < out.count; ++i)
                    if (out.arcs[i] == ai) { dup = true; break; }
                if (!dup) out.push(ai);
            }
        }
    };

    for (std::size_t ci : S.node(region).incident_chords) {
        assert(ci < S.num_chords());
        assert(!S.chord(ci).dead &&
               "[C91 §2.4]: normal-form (compacted) submap must have no "
               "dead chords in incident_chords");
        check_adj(S.chord(ci).left_adj);
        check_adj(S.chord(ci).right_adj);
    }

    // [C91 §2.4(iii) tex 138]: normal-form ⟹ start_arc/end_arc are set.
    auto check_endpoint_arc = [&](std::size_t arc_idx) {
        assert(arc_idx != NONE && arc_idx < S.num_arcs() &&
               "[C91 §2.4(iii)]: endpoint arc must be set for normal-form submaps");
        assert(!S.arc(arc_idx).dead &&
               "[C91 §2.4(iii)]: endpoint arc must be live");
        if (S.arc(arc_idx).region_node != region) return;
        bool dup = false;
        for (std::size_t i = 0; i < out.count; ++i)
            if (out.arcs[i] == arc_idx) { dup = true; break; }
        if (!dup) out.push(arc_idx);
    };
    check_endpoint_arc(S.start_arc);
    check_endpoint_arc(S.end_arc);

    // [C91 §3.1 tex 181]: degree ≤ 4 ⟹ at most 4 arcs (one per chord-gap).
    assert(out.count <= 4 &&
           "[C91 §3.1 tex 181]: conformal region has at most 4 arcs (degree ≤ 4)");

    return out;
}

// ── local_shoot ─────────────────────────────────────────────────

RayHit local_shoot(Point p, Side direction,
                    std::size_t region,
                    const Submap& S, const Polygon& C,
                    const RayShootingOracle& oracle) {
    // [C91 §3.1 tex 181]: check each region arc, take the nearest hit.
    // Conformality bounds the loop to ≤ 4 arcs.
    RegionArcs arcs = collect_region_arcs(S, region);

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

        RayHit hit = oracle.shoot(p, direction, ai, sub);
        if (!hit.hit) continue;

        // [C91 §3.0(i) tex 169]: oracle reports forward hits only.  A backward
        // hit (negative signed distance) would beat every legitimate hit
        // and silently corrupt the nearest-of selection.
        double hit_signed_dist = (direction == LEFT)
            ? (p.x - hit.x) : (hit.x - p.x);
        assert(hit_signed_dist >= 0.0 &&
               "[C91 §3.0(i) tex 169]: oracle must report a forward hit");

        hit.hit_arc_idx = ai;  // propagate O(1) context

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
            
            if (hit_dist < best_dist) {
                best = hit;
            } else if (hit_dist == best_dist) {
                // [C91 §2.1 tex 72]: Double-boundary disambiguation.  ∂C's two
                // sides share geometric location but are topologically
                // distinct; the ray hits the "first face" depending on its
                // travel direction and the edge's orientation.
                assert(hit.edge < C.num_edges());
                const auto& e = C.edge(hit.edge);
                bool edge_ascending = symbolic_y_less(
                    symbolic_y_of(C.vertex(e.start_idx)),
                    symbolic_y_of(C.vertex(e.end_idx)));

                Side minus_x_face = edge_ascending ? LEFT : RIGHT;
                Side plus_x_face  = edge_ascending ? RIGHT : LEFT;

                Side struck_first_face;
                if (best_dist == 0.0) {
                    // d = 0: Crossing the immediate zero-width canal boundary.
                    // Ray traveling RIGHT (crosses from -X to +X) hits +X face next.
                    // Ray traveling LEFT (crosses from +X to -X) hits -X face next.
                    struck_first_face = (direction == RIGHT) ? plus_x_face : minus_x_face;
                } else {
                    // d > 0: Hitting a distant obstacle.
                    // Ray traveling RIGHT (from -X) hits the -X face first.
                    // Ray traveling LEFT (from +X) hits the +X face first.
                    struck_first_face = (direction == RIGHT) ? minus_x_face : plus_x_face;
                }

                if (hit.side == struck_first_face && best.side != struck_first_face) {
                    best = hit;
                }
            }
        }
    }

    // [C91 §3.1 tex 181]: regions are closed under visibility (Lemma 2.1
    // corollary) ⟹ shooting from inside a region always hits.
    assert(best.hit &&
           "[C91 §3.1 tex 181]: local shoot inside a region must hit (Lemma 2.1)");
    return best;
}

// ── shooting_direction ──────────────────────────────────────────

Side shooting_direction(std::size_t edge, Side side,
                         const Polygon& C) {
    // [C91 §3.1]: "because of the double boundary the shooting direction is
    // always uniquely defined."  Edge going up (start.y < end.y): LEFT
    // side is on the west, shoots east (RIGHT).  Edge going down: LEFT
    // shoots west (LEFT).  RIGHT side is the reverse.
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
    // [C91 §3.1 tex 179]: sequence starts at a₀ (RIGHT companion of junction)
    // and ends at a_{m+1} (LEFT companion) — the cw tour of ∂C₁.
    assert(state.sequence.size() >= 2 &&
           "[C91 §3.1]: fusion sequence needs at least a₀ and a_{m+1}");
    const auto& a0 = state.sequence[0];
    assert(a0.is_companion && a0.side == RIGHT &&
           "[C91 §3.1 tex 179]: sequence[0] = a₀ (RIGHT companion at junction)");
    const auto& a_m1 = state.sequence.back();
    assert(a_m1.is_companion && a_m1.side == LEFT &&
           "[C91 §3.1 tex 179]: sequence.back() = a_{m+1} (LEFT companion)");

    // [C91 §3.1]: "Using local shooting, we find the point of ∂C₁ that a₀
    // sees with respect to C₁."  a₀ is at c_end of C₁, RIGHT side.
    assert(S1.end_arc != NONE &&
           "[C91 §2.4(iii)]: S₁ must have end_arc set (normal form)");
    std::size_t a0_region_s1 = S1.arc(S1.end_arc).region_node;
    Side a0_dir = shooting_direction(a0.edge, a0.side, C1);

    // a₀ is the junction vertex; use its coords directly.
    std::size_t junction_v = C1.num_vertices() - 1;
    Point a0_point = C1.vertex(junction_v);

    RayHit hit_c1 = local_shoot(a0_point, a0_dir, a0_region_s1,
                                 S1, C1, oracle1);

    // [C91 §3.1]: a₀ "most often" coincides with a ∂C₂ point but not at
    // y-extrema due to duplication.  SoS ([C91 §2 tex 47]) gives every vertex
    // a unique symbolic y, so visibility is well-defined either way.
    //
    // [C91 §3.1]: "the normal-form representation of S₂ ... lets us find,
    // in constant time, in which region of S₂ the point a₀ lies."
    assert(S2.start_arc != NONE &&
           "[C91 §2.4(iii) tex 138]: S₂ must have start_arc set (normal form)");
    std::size_t a0_region_s2 = S2.arc(S2.start_arc).region_node;

    RayHit hit_c2 = local_shoot(a0_point, a0_dir, a0_region_s2,
                                 S2, C2, oracle2);

    // [C91 §3.1 tex 185]: c₀ = the closer of hit_c1 and hit_c2.  Both are
    // guaranteed to hit (Lemma 2.1).
    assert(hit_c1.hit && "[C91 §3.1]: a₀ must see ∂C₁ (Lemma 2.1)");
    assert(hit_c2.hit && "[C91 §3.1]: a₀ must see ∂C₂ (Lemma 2.1)");
    bool c0_on_c2;
    RayHit c0;
    {
        double d1 = (a0_dir == LEFT)
            ? (a0_point.x - hit_c1.x) : (hit_c1.x - a0_point.x);
        double d2 = (a0_dir == LEFT)
            ? (a0_point.x - hit_c2.x) : (hit_c2.x - a0_point.x);
        // [C91 §3.1 tex 191]: "for simplicity we still go on saying that
        // c₀ sees a point of ∂C₂ with respect to C."  The paper allows
        // edge cases (a₀ at a y-extremum where c₀ technically cannot see
        // ∂C₂) and explicitly mandates defaulting to "c₀ on ∂C₂".  The
        // `d2 <= d1` (rather than `d2 < d1`) selects Case 1 in ties,
        // implementing the paper's tex 191 default.
        if (d2 <= d1) {
            c0 = hit_c2; c0_on_c2 = true;
        } else {
            c0 = hit_c1; c0_on_c2 = false;
        }
    }

    // [C91 §3.1]: "Another minor problem is that a₀c₀ might lie on
    // an exit chord of S₂ and thus there might be more than one
    // candidate for the status of current region.  We break ties by
    // electing the region that we locally enter as we leave p in a
    // clockwise traversal of ∂C₁."
    auto resolve_s2_region = [&](std::size_t initial_region, bool leaving_downward) -> std::size_t {
        // [C91 §3.1 tex 191]: If a₀ lies on an S₂ chord, that chord
        // must be incident on initial_region.  [C91 §2.3 tex 114]: conformal
        // submaps have node-degree ≤ 4, so the incident_chords scan is
        // O(1).  (The paper does not separately bound the tie-break;
        // startup as a whole is O(f(γ₂)) per tex 220, so this O(1)
        // sub-operation fits well within budget.)
        SymbolicY a0y = a0.y;
        for (std::size_t ci : S2.node(initial_region).incident_chords) {
            assert(!S2.chord(ci).dead &&
                   "[C91 §2.4]: normal-form (compacted) submap must have no "
                   "dead chords in incident_chords");
            const auto& ch = S2.chord(ci);
            if (!symbolic_y_equal(ch.symbolic_y(), a0y))
                continue;

            // [C91 §2.1 tex 72]: junction has TWO ∂C₂ companions, each
            // incident on one chord → up to two distinct S₂ chords share
            // a₀'s symbolic y.  Pick the SPECIFIC chord that a₀c₀ lies on:
            // one endpoint at a₀'s position on a₀'s ∂C side.  In C₂'s
            // frame, junction is C₂.vertex(0), so its edge is C₂.edge(0).
            // (a0.edge is in C₁'s frame — that's why we hardcode 0 here.)
            const bool incident_on_a0 =
                (ch.left_edge  == 0 && ch.left_side  == a0.side) ||
                (ch.right_edge == 0 && ch.right_side == a0.side);
            if (!incident_on_a0) continue;

            // [C91 §3.1 tex 191]: a₀c₀ lies on this S₂ chord.  Elect "the
            // region we locally enter as we leave p clockwise on ∂C₁."
            // `leaving_downward` encodes the y-direction of motion from p
            // (set differently for Case 1 leaving a₀ vs Case 2 leaving c₀).
            //
            // For each adj arc at the chord endpoint, look at the polygon
            // vertex adjacent to the chord along the arc's traversal; its
            // SoS y vs chord.symbolic_y tells us which side the arc's
            // region is on.  Prefer a mid-edge endpoint (its two adj arcs
            // straddle the chord per [C91 §2.2 tex 96]).  If both endpoints
            // are polygon vertices, fall through to the single-arc-per-side
            // case below.
            const Chord::AdjArcs* adj = nullptr;
            std::size_t ch_edge_at_endpoint = NONE;
            Side ch_side_at_endpoint = LEFT;
            if (ch.left_adj.count == 2) {
                adj = &ch.left_adj;
                ch_edge_at_endpoint = ch.left_edge;
                ch_side_at_endpoint = ch.left_side;
            } else if (ch.right_adj.count == 2) {
                adj = &ch.right_adj;
                ch_edge_at_endpoint = ch.right_edge;
                ch_side_at_endpoint = ch.right_side;
            }

            // Returns true iff the arc's region is strictly ABOVE the
            // chord under SoS.  Looks at the polygon vertex adjacent to
            // the chord along the arc's traversal direction.
            auto arc_region_above_chord = [&](const Arc& a) -> bool {
                // Structural: does the arc START at the chord (first_edge
                // matches the chord's slot) or END there (last_edge matches)?
                bool first_matches =
                    (a.first_edge == ch_edge_at_endpoint &&
                     a.first_side == ch_side_at_endpoint);
                bool last_matches =
                    (a.last_edge == ch_edge_at_endpoint &&
                     a.last_side == ch_side_at_endpoint);
                bool arc_starts_at_chord;
                if (first_matches && !last_matches) {
                    arc_starts_at_chord = true;
                } else if (!first_matches && last_matches) {
                    arc_starts_at_chord = false;
                } else {
                    // Single-edge arc on chord's edge — disambiguate via
                    // key_y_tag (matches chord.y_tag ⟺ arc starts here).
                    assert(first_matches && last_matches &&
                           "adj arc must touch the chord endpoint via "
                           "first_edge or last_edge");
                    arc_starts_at_chord = (a.key_y_tag == ch.y_tag);
                }

                // Vertex adjacent to the chord along the arc's traversal
                // (i.e. the next vertex AWAY from the chord).  LEFT
                // ascends, RIGHT descends.
                std::size_t adj_v_idx;
                if (arc_starts_at_chord) {
                    adj_v_idx = (a.first_side == LEFT)
                              ? C2.edge(a.first_edge).end_idx
                              : C2.edge(a.first_edge).start_idx;
                } else {
                    adj_v_idx = (a.last_side == LEFT)
                              ? C2.edge(a.last_edge).start_idx
                              : C2.edge(a.last_edge).end_idx;
                }

                // SoS ⟹ every vertex has a unique perturbed y, so this
                // comparison is definitive.
                return symbolic_y_greater(
                    symbolic_y_of(C2.vertex(adj_v_idx)), ch.symbolic_y());
            };

            if (adj) {
                // Mid-edge endpoint: the two adj arcs are on OPPOSITE
                // sides of the chord (one above, one below per
                // [C91 §2.2 tex 96]).  Classify each via the helper.
                bool arc0_above = arc_region_above_chord(S2.arc(adj->arcs[0]));
                bool arc1_above = arc_region_above_chord(S2.arc(adj->arcs[1]));
                assert(arc0_above != arc1_above &&
                       "[C91 §2.2 tex 96]: the two adj arcs at a mid-edge "
                       "chord endpoint must lie on opposite sides of the chord");
                std::size_t r0 = S2.arc(adj->arcs[0]).region_node;
                std::size_t r1 = S2.arc(adj->arcs[1]).region_node;
                std::size_t above_r = arc0_above ? r0 : r1;
                std::size_t below_r = arc0_above ? r1 : r0;
                return leaving_downward ? below_r : above_r;
            }

            // Vertex-to-vertex case: both endpoints are polygon vertices
            // (e.g., endpoint companion chord at the junction per
            // [C91 §2.1 tex 72]).  Each endpoint has one adj arc; classify it
            // using the same helper with the appropriate endpoint slot.
            ch_edge_at_endpoint = ch.left_edge;
            ch_side_at_endpoint = ch.left_side;
            bool left_arc_above =
                arc_region_above_chord(S2.arc(ch.left_adj.arcs[0]));
            std::size_t r_arc = S2.arc(ch.left_adj.arcs[0]).region_node;
            std::size_t r_other = (ch.region[0] == r_arc)
                                ? ch.region[1] : ch.region[0];
            std::size_t above_r = left_arc_above ? r_arc : r_other;
            std::size_t below_r = left_arc_above ? r_other : r_arc;
            return leaving_downward ? below_r : above_r;
        }
        return initial_region;
    };

    // Record a₀c₀ — same in both cases (fields depend only on a₀, c₀).
    // DiscoveredChord slot convention (chord.h, fusion.h): LEFT slot holds
    // the LEFT-side endpoint; RIGHT slot holds the RIGHT-side endpoint.
    assert(c0.side == LEFT && a0.side == RIGHT &&
           "DiscoveredChord slot convention: c0 on LEFT, a0 on RIGHT");
    state.chords.push_back({a0.y, c0.edge, c0.side, a0.edge, a0.side});

    if (c0_on_c2) {
        // [C91 §3.1 Case 1]: "set p = a₀, call the region of S₂ crossed by
        // a₀c₀ current."  Main loop starts at k=1 (A_k = arc a₀→a₁).
        state.p = a0_point;
        state.p_edge = a0.edge;
        state.p_side = a0.side;

        // CW step from a₀ goes from junction_v to junction_v-1; its
        // y-direction is the topological "leaving p" vector.
        SymbolicY y_here = symbolic_y_of(C1.vertex(junction_v));
        SymbolicY y_next = symbolic_y_of(C1.vertex(junction_v - 1));
        bool a0_leaving_downward = symbolic_y_less(y_next, y_here);

        state.s2_region = resolve_s2_region(a0_region_s2, a0_leaving_downward);
        return 1;
    } else {
        // [C91 §3.1 Case 2]: "skip all the way to c₀ ... set p = c₀, call
        // the region of S₂ containing a₀ current."  c₀ is on an edge
        // interior, not a named vertex (index = NONE).
        state.p = Point{c0.x, c0.y, NONE};
        state.p_edge = c0.edge;
        state.p_side = c0.side;

        // Topological slope leaving c₀.
        assert(c0.edge < C1.num_edges());
        const auto& e = C1.edge(c0.edge);
        bool edge_ascending = symbolic_y_less(
            symbolic_y_of(C1.vertex(e.start_idx)),
            symbolic_y_of(C1.vertex(e.end_idx)));
        bool c0_leaving_downward =
            (c0.side == LEFT) ? !edge_ascending : edge_ascending;

        // Always invoke the tie-break: a₀ may be a reflex vertex generated
        // by an interior cut, allowing an S₂ chord at the exact same y.
        state.s2_region = resolve_s2_region(a0_region_s2, c0_leaving_downward);

        // Skip to c₀: find first fusion vertex at or past c₀'s ∂C₁ position.
        // O(1) via arc_starts (local_shoot's hit arc → sequence block) plus
        // a constant within-bucket scan (degree ≤ 4 ⟹ O(1) endpoints/bucket).
        std::size_t junction_edge = C1.num_edges() - 1;
        std::size_t right_half_len = junction_edge + 1;
        auto trav_pos = [&](std::size_t edge, Side side) -> std::size_t {
            // CW traversal: RIGHT half descends in edge, LEFT half ascends.
            return (side == RIGHT) ? junction_edge - edge
                                   : right_half_len + edge;
        };

        // c₀ shares its raw y with a₀ (horizontal visibility chord), but
        // it's mid-edge so it has no vertex tag of its own.  We borrow a₀'s
        // tag (= junction's tag).  Downstream comparisons (vertex_past_c0)
        // only pair c₀ against non-companion chord endpoints whose tags
        // differ from the junction's, so the SoS tie-break is never
        // ambiguous and raw-y equality is unproblematic.
        SymbolicY c0_y = a0.y;
        std::size_t c0_pos = trav_pos(c0.edge, c0.side);

        // Within an edge, ∂C₁ traversal direction:
        //   LEFT side  follows the edge (start→end).
        //   RIGHT side reverses it (end→start).
        auto vertex_past_c0 = [&](const FusionVertex& v) -> bool {
            std::size_t v_pos = trav_pos(v.edge, v.side);
            if (v_pos != c0_pos) return v_pos > c0_pos;
            assert(v.edge < C1.num_edges());
            const auto& e = C1.edge(v.edge);
            bool edge_ascending = symbolic_y_less(
                symbolic_y_of(C1.vertex(e.start_idx)),
                symbolic_y_of(C1.vertex(e.end_idx)));
            bool traversal_ascending =
                (v.side == LEFT) ? edge_ascending : !edge_ascending;
            return traversal_ascending ? symbolic_y_geq(v.y, c0_y)
                                       : symbolic_y_leq(v.y, c0_y);
        };

        // O(1) skip: c0's hit arc indexes directly into the sequence.
        assert(hit_c1.hit_arc_idx != NONE && "[C91 §3.1]: c0 must carry hit context");
        std::size_t c0_arc_idx = hit_c1.hit_arc_idx;
        
        std::size_t cw_key = 0;
        {
            std::size_t lrb = S1.left_right_boundary();
            std::size_t num_right = S1.num_arcs() - lrb;
            cw_key = (c0_arc_idx >= lrb) ? (c0_arc_idx - lrb) : (num_right + c0_arc_idx);
        }
        
        std::size_t lo = state.arc_starts[cw_key];

        // The exact match within the arc's block (degree <= 4 → O(1) iteration).
        while (lo < state.sequence.size() && !vertex_past_c0(state.sequence[lo])) {
            lo++;
        }

        // [C91 §3.1 tex 179, 188]: c₀ strictly precedes a_{m+1} on the cw
        // tour, and a_{m+1} sits at the maximal trav_pos → the while loop
        // always terminates within bounds.
        assert(lo < state.sequence.size() &&
               "[C91 §3.1 tex 179, 188]: skip-to-c₀ must land within sequence");

        // [C91 §3.1 Case 2 / Lemma 2.1]: skipped vertices see ∂C₁ — paper-
        // guaranteed (their chord ends on ∂C₁ before c₀ in cw order).
        for (std::size_t skipped_i = 1; skipped_i < lo; ++skipped_i) {
            if (!state.sequence[skipped_i].is_companion) {
                assert(state.sequence[skipped_i].chord_idx != NONE &&
                       "[C91 §3.1 Case 2 (Lemma 2.1)]: skipped vertex must be a "
                       "chord endpoint of S₁ (sees ∂C₁ by construction)");
            }
        }

        return lo;
    }
}

// ── fuse_s1_into_s2 ─────────────────────────────────────────────

void fuse_s1_into_s2(FusionState& state,
                      const Submap& S1, const Polygon& C1,
                      const Submap& S2, const Polygon& C2,
                      const RayShootingOracle& oracle1,
                      const RayShootingOracle& oracle2) {
    build_fusion_sequence(state, S1, C1);

    // [C91 §3.1 tex 197]: k = index of arc A_k of S₁ containing p, tie-break
    // p ≠ a_k.  A_j runs from a_{j-1} to a_j cw on ∂C₁ (A_1 = a_0→a_1,
    // A_{m+1} = a_m→a_{m+1}).
    std::size_t k = fusion_startup(state, S1, C1, S2, C2, oracle1, oracle2);

    while (true) {
        // [C91 §3.1 tex 197]: p == a_{m+1} ⟹ terminate (no A_k defined).
        if (k >= state.sequence.size()) return;

        // [C91 §3.1 tex 195]: loop invariants (A) and (B).
        assert(!state.chords.empty() &&
               "[C91 §3.1 invariant (A)]: chord list non-empty");
        assert(state.s2_region != NONE &&
               state.s2_region < S2.num_nodes() &&
               !S2.node(state.s2_region).dead &&
               "[C91 §3.1 invariant (B)]: current S₂ region must be valid");
        assert(local_shoot(state.p,
                   shooting_direction(state.p_edge, state.p_side, C1),
                   state.s2_region, S2, C2, oracle2).hit &&
               "[C91 §3.1 invariant (B)]: p must see ∂C₂");

        // [C91 §3.1 tex 197]: R = current region; snapshot fixed for the
        // inner walk; updates take effect next outer iteration.
        [[maybe_unused]] const std::size_t R = state.s2_region;

        // [C91 §3.1 tex 200(i)]: a_j ∈ R AND a_j sees ∂C₂.  O(f(γ₂)) test
        // at [C91 §3.1 tex 220]: shoot a_j against each arc of R; no hit
        // ⟹ a_j ∉ R; else nearest hit's side+orientation confirms in-R
        // and gives s ∈ ∂C₂; compare to t = a_j's ∂C₁ hit (O(1) when a_j
        // is an S₁ chord endpoint, else local_shoot in S₁).
        // TODO ([C91 §3.1] — actions paragraph): wire predicate + action.
        auto case_i_fires = [&](std::size_t /*j*/) -> bool {
            return false;
        };

        // [C91 §3.1 tex 202(ii)]: (i) fails, but R has an exit chord whose
        // endpoint sees a point on A_j strictly after p.  O(f(γ₁)) test
        // at [C91 §3.1 tex 222]: per exit-chord endpoint a, shoot via
        // oracle1 toward A_j — disqualify on miss / off-chord / before-p /
        // wrong orientation; else verify by local_shoot from hit s back
        // toward a in S₁ (passes through a ⟺ s and a see each other).
        // TODO ([C91 §3.1] — actions paragraph): wire predicate + action.
        auto case_ii_fires = [&](std::size_t /*j*/) -> bool {
            return false;
        };

        for (std::size_t j = k; ; ++j) {
            // [C91 §3.1 tex 203(iii) + tex 206 action]: j == m+2 ⟹ stop;
            // remaining a_k..a_{m+1} see ∂C₁ (informational — no chords
            // recorded, only S₁'s existing chords carry over in rebuild).
            if (j == state.sequence.size()) return;

            state.current_stop = j;

            if (case_i_fires(j)) {
                // TODO ([C91 §3.1 tex 206 case (i) action]):
                //   - Find q ∈ ∂C₂ seen by a_j (local_shoot from a_j in R).
                //   - Record chord a_j → q in state.chords.
                //   - p ← a_j, R unchanged.
                //   - k ← j+1 (per A_k tie-break p ≠ a_k).
                // Intervening a_i (k ≤ i < j) see ∂C₁ — informational.
                break;
            }
            if (case_ii_fires(j)) {
                // TODO ([C91 §3.1 tex 206 case (ii) action]):
                //   - Among (ii) candidates, pick endpoint q' whose seen
                //     point p' is the LAST one cw on ∂C₁ starting from p.
                //   - Record chord q' → p'.
                //   - p ← p'; s2_region ← region entered when locally
                //     crossing the exit chord at p' along ∂C₁.
                //   - k ← j (p' lies in A_j's interior; tie-break not
                //     triggered since p' isn't a fusion vertex).
                // Intervening a_i (k ≤ i < j) see ∂C₁ — informational.
                break;
            }
        }
    }
}

// ── build_fusion_sequence ───────────────────────────────────────

void build_fusion_sequence(FusionState& state, const Submap& S,
                           const Polygon& C) {
    std::size_t n = C.num_vertices();
    assert(n >= 2);

    // Junction = C₁ ∩ C₂ = last vertex of C₁.
    std::size_t junction_v = n - 1;
    std::size_t junction_edge = junction_v - 1;

    // [C91 §3.1]: a₀ = RIGHT companion at junction (start of cw tour);
    // a_{m+1} = LEFT companion (end).  Both inherit junction's symbolic y.
    SymbolicY junction_y = symbolic_y_of(C.vertex(junction_v));

    FusionVertex a_0;
    a_0.y = junction_y;
    a_0.edge = junction_edge;
    a_0.side = RIGHT;
    a_0.chord_idx = NONE;
    a_0.is_left_endpoint = false;
    a_0.is_companion = true;

    FusionVertex a_m1;
    a_m1.y = junction_y;
    a_m1.edge = junction_edge;
    a_m1.side = LEFT;
    a_m1.chord_idx = NONE;
    a_m1.is_left_endpoint = false;
    a_m1.is_companion = true;

    // [C91 §3.1]: a₁, ..., aₘ = canonical enumeration of S₁'s exit chord
    // endpoints in cw ∂C₁ order.  [C91 §2.4(iii) tex 142]: arc-sequence is in
    // ∂C order (LEFT ascending, then RIGHT descending).  cw ∂C₁ from a₀
    // is RIGHT then LEFT; we counting-sort endpoints by cw-arc position
    // in O(m).  Each endpoint is associated with one arc via adj_arcs:
    // vertex endpoint ⟹ the single adj arc; mid-edge endpoint ⟹ the
    // "starting" adj arc whose key_y == chord.y.

    std::size_t num_arcs = S.num_arcs();
    std::size_t lrb = S.left_right_boundary();
    std::size_t num_right = num_arcs - lrb;

    // Arc-table index → cw position.  RIGHT arcs [lrb..end) come first,
    // then LEFT arcs [0..lrb).
    auto cw_pos = [&](std::size_t arc_idx) -> std::size_t {
        return (arc_idx >= lrb) ? arc_idx - lrb
                                : num_right + arc_idx;
    };

    struct KeyedVertex {
        FusionVertex v;
        std::size_t key;        // cw arc position
    };
    std::vector<KeyedVertex> endpoints;
    endpoints.reserve(2 * S.num_live_chords());

    auto make_vertex = [](const Chord& c, std::size_t ci,
                           bool is_left) -> FusionVertex {
        FusionVertex fv;
        fv.y = c.symbolic_y();
        fv.edge = is_left ? c.left_edge : c.right_edge;
        fv.side = is_left ? c.left_side : c.right_side;
        fv.chord_idx = ci;
        fv.is_left_endpoint = is_left;
        fv.is_companion = false;
        return fv;
    };

    // Mid-edge endpoint: of the two adj arcs, exactly one starts at the
    // chord (its key_y matches the chord's y).
    auto starting_arc = [&](const Chord::AdjArcs& adj,
                             SymbolicY chord_y) -> std::size_t {
        assert(adj.count == 2);
        if (symbolic_y_equal(S.arc(adj.arcs[0]).key_symbolic_y(), chord_y))
            return adj.arcs[0];
        assert(symbolic_y_equal(S.arc(adj.arcs[1]).key_symbolic_y(), chord_y) &&
               "[C91 §2.4]: one adj arc must start at the chord endpoint");
        return adj.arcs[1];
    };

    for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
        assert(!S.chord(ci).dead && "[C91 §2.4]: normal-form submap has no dead chords");
        const auto& c = S.chord(ci);
        // [C91 §3.1 tex 179]: only EXIT chords enter the enumeration.
        // Null-length chords ([C91 §2.2 tex 96]) are carried over to the fused
        // submap directly (tex 224); including them would blow Lemma 3.1's
        // O(n/γ + 1) bound (tex 209–210).
        if (c.is_null_length) continue;
        SymbolicY cy = c.symbolic_y();

        // LEFT endpoint → associated arc.
        std::size_t left_arc = (c.left_adj.count == 2)
            ? starting_arc(c.left_adj, cy)
            : c.left_adj.arcs[0];
        endpoints.push_back({make_vertex(c, ci, true), cw_pos(left_arc)});

        // RIGHT endpoint → associated arc.
        std::size_t right_arc = (c.right_adj.count == 2)
            ? starting_arc(c.right_adj, cy)
            : c.right_adj.arcs[0];
        endpoints.push_back({make_vertex(c, ci, false), cw_pos(right_arc)});
    }

    // Counting sort by clockwise arc key.  O(m) time, O(m) space.
    // [C91 §2.4 tex 142]: "canonical vertex enumerations in optimal
    // time" — O(m) where m = number of chord endpoints.
    std::vector<KeyedVertex> sorted(endpoints.size());
    std::vector<std::size_t> bucket_starts;
    if (!endpoints.empty()) {
        std::vector<std::size_t> counts(num_arcs + 1, 0);
        for (const auto& ep : endpoints)
            ++counts[ep.key + 1];
        for (std::size_t i = 1; i <= num_arcs; ++i)
            counts[i] += counts[i - 1];
            
        bucket_starts.assign(counts.begin(), counts.begin() + static_cast<std::ptrdiff_t>(num_arcs));
        
        for (const auto& ep : endpoints)
            sorted[counts[ep.key]++] = ep;
    } else {
        bucket_starts.assign(num_arcs, 0);
    }

    // Within each bucket (same arc): O(1) elements (conformal degree ≤ 4).
    // Insertion-sort by y along the arc's traversal direction; O(m) total.
    std::size_t right_half_len = junction_edge + 1;
    auto trav_pos = [&](std::size_t edge, Side side) -> std::size_t {
        if (side == RIGHT)
            return junction_edge - edge;
        else
            return right_half_len + edge;
    };
    auto vertex_before = [&](const FusionVertex& u, const FusionVertex& v) -> bool {
        std::size_t u_pos = trav_pos(u.edge, u.side);
        std::size_t v_pos = trav_pos(v.edge, v.side);
        if (u_pos != v_pos) return u_pos < v_pos;
        assert(u.edge < C.num_edges());
        const auto& e = C.edge(u.edge);
        bool edge_ascending = symbolic_y_less(
            symbolic_y_of(C.vertex(e.start_idx)),
            symbolic_y_of(C.vertex(e.end_idx)));
        bool trav_asc = (u.side == LEFT) ? edge_ascending : !edge_ascending;
        return trav_asc ? symbolic_y_less(u.y, v.y) : symbolic_y_greater(u.y, v.y);
    };

    {
        std::size_t i = 0;
        while (i < sorted.size()) {
            std::size_t j = i + 1;
            while (j < sorted.size() && sorted[j].key == sorted[i].key)
                ++j;
            // sorted[i..j) share the same arc key.
            if (j - i > 1) {
                // Insertion sort — O(1) elements per bucket.
                for (std::size_t a = i + 1; a < j; ++a) {
                    KeyedVertex tmp = sorted[a];
                    std::size_t b = a;
                    while (b > i) {
                        if (!vertex_before(tmp.v, sorted[b-1].v)) break;
                        sorted[b] = sorted[b-1];
                        --b;
                    }
                    sorted[b] = tmp;
                }
            }
            i = j;
        }
    }

    // [C91 §3.1]: a₀ and a_{m+1} are not in the canonical enumeration — add them.
    std::vector<FusionVertex> result;
    result.reserve(sorted.size() + 2);
    result.push_back(a_0);
    for (const auto& kv : sorted)
        result.push_back(kv.v);
    result.push_back(a_m1);

    state.sequence = std::move(result);

    // Map arc key → start of its block in `sequence`.  +1 for the prepended a₀.
    state.arc_starts.resize(num_arcs);
    for (std::size_t i = 0; i < num_arcs; ++i)
        state.arc_starts[i] = bucket_starts[i] + 1;
}

} // namespace chazelle
