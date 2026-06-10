// src/merge/fusion.cpp

#include "fusion.h"

#include <algorithm>

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
                    const RayShootingOracle& oracle,
                    bool require_hit) {
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
    if (require_hit)
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
        // `d2 <= d1` (rather than `d2 < d1`) selects case (i) in ties,
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
            // (set differently for case (i) leaving a₀ vs case (ii) leaving c₀).
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
            auto arc_region_above_chord = [&](std::size_t arc_idx) -> bool {
                const Arc& a = S2.arc(arc_idx);
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
                    // [C91 §2.4 tex 133]: single-edge arc on chord's edge —
                    // the arc that starts at the chord has start-y = chord.y
                    // (derived from the bounding chord, not stored).
                    assert(first_matches && last_matches &&
                           "adj arc must touch the chord endpoint via "
                           "first_edge or last_edge");
                    arc_starts_at_chord = symbolic_y_equal(
                        S2.arc_start_symbolic_y(arc_idx, C2),
                        ch.symbolic_y());
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
                bool arc0_above = arc_region_above_chord(adj->arcs[0]);
                bool arc1_above = arc_region_above_chord(adj->arcs[1]);
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
                arc_region_above_chord(ch.left_adj.arcs[0]);
            std::size_t r_arc = S2.arc(ch.left_adj.arcs[0]).region_node;
            std::size_t r_other = (ch.region[0] == r_arc)
                                ? ch.region[1] : ch.region[0];
            std::size_t above_r = left_arc_above ? r_arc : r_other;
            std::size_t below_r = left_arc_above ? r_other : r_arc;
            return leaving_downward ? below_r : above_r;
        }
        return initial_region;
    };

    // Record a₀c₀ — slot order = ascending x.  Per [C91 §2.1 tex 72]
    // `Side` labels the snake's geometric face (LEFT = walker-from-
    // start-to-end's left); two ∂C points of a horizontal chord can
    // share a Side (e.g. LEFT-face-of-ascending-edge + LEFT-face-of-
    // descending-edge), so slot name is independent of Side value.
    // [C91 §3.1 tex 224]: a₀ lives in C₁'s edge frame; c₀ lives in C₁
    // (case (ii)) or C₂ (case (i)).  rebuild_submap translates accordingly.
    const bool a0_on_c1 = true;
    const bool c0_on_c1 = !c0_on_c2;
    if (a0_point.x < c0.x)
        state.chords.push_back({a0.y, a0.edge, a0.side, c0.edge, c0.side,
                                a0_on_c1, c0_on_c1});
    else
        state.chords.push_back({a0.y, c0.edge, c0.side, a0.edge, a0.side,
                                c0_on_c1, a0_on_c1});

    if (c0_on_c2) {
        // [C91 §3.1 case (i)]: "set p = a₀, call the region of S₂ crossed by
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
        // [C91 §3.1 case (ii)]: "skip all the way to c₀ ... set p = c₀, call
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

        // [C91 §3.1 case (ii) / Lemma 2.1]: skipped vertices see ∂C₁ — paper-
        // guaranteed (their chord ends on ∂C₁ before c₀ in cw order).
        for (std::size_t skipped_i = 1; skipped_i < lo; ++skipped_i) {
            if (!state.sequence[skipped_i].is_companion) {
                assert(state.sequence[skipped_i].chord_idx != NONE &&
                       "[C91 §3.1 case (ii) (Lemma 2.1)]: skipped vertex must be a "
                       "chord endpoint of S₁ (sees ∂C₁ by construction)");
            }
        }

        return lo;
    }
}

// ── fuse_submaps ─────────────────────────────────────────────

void fuse_submaps(FusionState& state,
                      const Submap& S1, const Polygon& C1,
                      const Submap& S2, const Polygon& C2,
                      const RayShootingOracle& oracle1,
                      const RayShootingOracle& oracle2) {
    build_fusion_sequence(state, S1, C1);

    // [C91 §3.1 tex 224]: per-pass invalidation bitsets for rebuild_submap.
    state.invalidated_self.assign(S1.num_chords(), false);
    state.invalidated_other.assign(S2.num_chords(), false);

    // [C91 §3.1 tex 199]: k = index of arc A_k of S₁ containing p, tie-break
    // p ≠ a_k.  A_j runs from a_{j-1} to a_j cw on ∂C₁ (A_1 = a_0→a_1,
    // A_{m+1} = a_m→a_{m+1}).
    std::size_t k = fusion_startup(state, S1, C1, S2, C2, oracle1, oracle2);

    while (true) {
        // [C91 §3.1 tex 199]: p == a_{m+1} ⟹ terminate (no A_k defined).
        if (k >= state.sequence.size()) return;

        // [C91 §3.1 tex 194]: invariant (A) — "the points of ∂C seen by
        // a_0..a_{k-1} have all been determined already."  Proven by
        // induction over case (i)/(ii) action, not directly checkable
        // without auxiliary state; we trust the per-action chord pushes.
        //
        // [C91 §3.1 tex 195]: invariant (B) — p sees ∂C₂ from R.
        assert(state.s2_region != NONE &&
               state.s2_region < S2.num_nodes() &&
               !S2.node(state.s2_region).dead &&
               "[C91 §3.1 invariant (B)]: current S₂ region must be valid");
        assert(local_shoot(state.p,
                   shooting_direction(state.p_edge, state.p_side, C1),
                   state.s2_region, S2, C2, oracle2).hit &&
               "[C91 §3.1 invariant (B)]: p must see ∂C₂");

        // [C91 §3.1 tex 199]: R = current region; snapshot fixed for the
        // inner walk; updates take effect next outer iteration.
        const std::size_t R = state.s2_region;

        // Position of a fusion vertex: edge-interp at chord y, or junction
        // for companions.
        auto fv_point = [&](const FusionVertex& v) -> Point {
            if (v.is_companion)
                return C1.vertex(C1.num_vertices() - 1);
            const auto& e = C1.edge(v.edge);
            const Point& vs = C1.vertex(e.start_idx);
            const Point& ve = C1.vertex(e.end_idx);
            double y = v.y.y;
            double t = (y - vs.y) / (ve.y - vs.y);
            return Point{vs.x + t * (ve.x - vs.x), y, NONE};
        };

        // [C91 §3.1 tex 220]: case (i) test in O(f(γ₂)).
        struct CaseIResult { bool fires; RayHit s_hit; };
        auto case_i_test = [&](std::size_t j) -> CaseIResult {
            const FusionVertex& aj_v = state.sequence[j];
            Point aj_point = fv_point(aj_v);
            Side dir = shooting_direction(aj_v.edge, aj_v.side, C1);

            // shoot from a_j; no hit ⟹ a_j ∉ R.
            RayHit s_hit = local_shoot(aj_point, dir, R, S2, C2, oracle2,
                                       /*require_hit=*/false);
            if (!s_hit.hit) return {false, {}};

            // hit.side matches arc's stored Side ⟺ a_j ∈ R.
            const Arc& s_arc = S2.arc(s_hit.hit_arc_idx);
            if (s_hit.side != s_arc.first_side &&
                s_hit.side != s_arc.last_side)
                return {false, {}};

            // t = ∂C₁ point a_j sees.  O(1) at S₁ chord endpoints
            // (other endpoint); local_shoot in S₁ for companions.
            double t_x;
            if (aj_v.chord_idx != NONE) {
                const Chord& ch = S1.chord(aj_v.chord_idx);
                std::size_t other_edge = aj_v.is_left_endpoint
                    ? ch.right_edge : ch.left_edge;
                const auto& e = C1.edge(other_edge);
                const Point& vs = C1.vertex(e.start_idx);
                const Point& ve = C1.vertex(e.end_idx);
                double t_param = (ch.y - vs.y) / (ve.y - vs.y);
                t_x = vs.x + t_param * (ve.x - vs.x);
            } else {
                std::size_t aj_region_s1 = S1.arc(S1.end_arc).region_node;
                t_x = local_shoot(aj_point, dir, aj_region_s1,
                                  S1, C1, oracle1).x;
            }

            // (i) fires iff s closer than t in dir.
            double s_dist = (dir == LEFT)
                ? (aj_point.x - s_hit.x) : (s_hit.x - aj_point.x);
            double t_dist = (dir == LEFT)
                ? (aj_point.x - t_x) : (t_x - aj_point.x);
            return {s_dist < t_dist, s_hit};
        };

        // Position of an S₂ chord endpoint: edge-interp at chord y.
        auto s2_endpoint_point = [&](std::size_t edge, double y) -> Point {
            const auto& e = C2.edge(edge);
            const Point& vs = C2.vertex(e.start_idx);
            const Point& ve = C2.vertex(e.end_idx);
            double t = (y - vs.y) / (ve.y - vs.y);
            return Point{vs.x + t * (ve.x - vs.x), y, NONE};
        };

        // O(1) lookup of A_j's S₁ arc index (precomputed inverse map).
        auto Aj_arc_idx = [&](std::size_t j) -> std::size_t {
            return (j < state.arc_for_seq_pos.size())
                ? state.arc_for_seq_pos[j] : NONE;
        };

        // CW position on ∂C₁ as (trav_pos, within_edge); lex compare
        // gives the cw walk.  LEFT side: within = edge param t; RIGHT
        // side: 1−t (cw traversal is end→start).
        auto cw_position = [&](double y, std::size_t edge, Side side)
                           -> std::pair<std::size_t, double> {
            std::size_t junction_edge = C1.num_edges() - 1;
            std::size_t right_half_len = junction_edge + 1;
            std::size_t tp = (side == RIGHT)
                ? (junction_edge - edge)
                : (right_half_len + edge);
            const auto& e = C1.edge(edge);
            const Point& vs = C1.vertex(e.start_idx);
            const Point& ve = C1.vertex(e.end_idx);
            double t = (y - vs.y) / (ve.y - vs.y);
            return {tp, (side == LEFT) ? t : (1.0 - t)};
        };

        // [C91 §3.1 tex 222]: case (ii) test in O(f(γ₁)).
        struct CaseIIResult {
            bool fires;
            RayHit p_prime_hit;
            std::size_t chord_idx;
            bool q_is_left;
        };
        auto case_ii_test = [&](std::size_t j) -> CaseIIResult {
            CaseIIResult best{false, {}, NONE, false};
            auto best_cw = std::make_pair<std::size_t, double>(0, -1.0);

            std::size_t aj_arc = Aj_arc_idx(j);
            if (aj_arc == NONE) return best;
            Subarc aj_sub{state.sequence[j-1].edge,
                          state.sequence[j-1].side,
                          state.sequence[j].edge,
                          state.sequence[j].side};

            Side aj_left_side = state.sequence[j-1].side;
            Side aj_right_side = state.sequence[j].side;
            bool aj_non_wrapping = (aj_left_side == aj_right_side);

            auto p_cw = cw_position(state.p.y, state.p_edge, state.p_side);

            for (std::size_t ci : S2.node(R).incident_chords) {
                const Chord& chord_ab = S2.chord(ci);
                if (chord_ab.dead || chord_ab.is_null_length) continue;

                Point a_pt = s2_endpoint_point(chord_ab.left_edge,
                                                chord_ab.y);
                Point b_pt = s2_endpoint_point(chord_ab.right_edge,
                                                chord_ab.y);

                for (bool is_left : {true, false}) {
                    std::size_t q_edge = is_left ? chord_ab.left_edge
                                                  : chord_ab.right_edge;
                    Side q_side = is_left ? chord_ab.left_side
                                          : chord_ab.right_side;
                    Point q_point = is_left ? a_pt : b_pt;
                    Point other_point = is_left ? b_pt : a_pt;
                    Side shoot_dir = shooting_direction(q_edge, q_side, C2);

                    // "hit does not lie on ab": shoot must traverse ab —
                    // i.e., shoot_dir points toward the other endpoint.
                    Side dir_to_other = (other_point.x < q_point.x)
                                          ? LEFT : RIGHT;
                    if (shoot_dir != dir_to_other) continue;

                    // shoot from a toward A_j; no hit ⟹ disqualify.
                    RayHit hit = oracle1.shoot(q_point, shoot_dir,
                                                aj_arc, aj_sub);
                    if (!hit.hit) continue;

                    // "proper orientation": hit.Side matches A_j's Side
                    // at hit.edge.  Wrapping A_j (only the RIGHT→LEFT
                    // transition) splits into edge-range zones.
                    bool hit_side_ok = aj_non_wrapping
                        ? (hit.side == aj_left_side)
                        : (hit.side == RIGHT
                            ? hit.edge <= state.sequence[j-1].edge
                            : hit.edge <= state.sequence[j].edge);
                    if (!hit_side_ok) continue;

                    // "occurs before p along A_j": strict cw compare.
                    auto hit_cw = cw_position(hit.y, hit.edge, hit.side);
                    if (hit_cw <= p_cw) continue;

                    // back-shot from s in its natural direction; q' must
                    // lie between s and t for s↔a visibility.
                    Side s_back_dir = shooting_direction(hit.edge, hit.side,
                                                         C1);
                    Point s_point{hit.x, hit.y, NONE};
                    RayHit t_hit = local_shoot(
                        s_point, s_back_dir,
                        S1.arc(aj_arc).region_node,
                        S1, C1, oracle1);
                    double s_to_q = (s_back_dir == LEFT)
                        ? (s_point.x - q_point.x)
                        : (q_point.x - s_point.x);
                    double s_to_t = (s_back_dir == LEFT)
                        ? (s_point.x - t_hit.x)
                        : (t_hit.x - s_point.x);
                    if (s_to_q < 0.0 || s_to_q > s_to_t) continue;

                    // [C91 §3.1 tex 206]: pick LAST p' cw from p.
                    if (hit_cw > best_cw) {
                        best = {true, hit, ci, is_left};
                        best_cw = hit_cw;
                    }
                }
            }
            return best;
        };

        for (std::size_t j = k; ; ++j) {
            // [C91 §3.1 tex 203(iii) + tex 206 action]: j == m+2 ⟹ stop;
            // remaining a_k..a_{m+1} see ∂C₁ (informational — no chords
            // recorded, only S₁'s existing chords carry over in rebuild).
            if (j == state.sequence.size()) return;

            state.current_stop = j;

            // [C91 §3.1 tex 206 case (i)]: record a_j → s_hit; p = a_j;
            // R unchanged; k = j+1 (A_k tie-break).
            if (auto r = case_i_test(j); r.fires) {
                const FusionVertex& aj_v = state.sequence[j];
                Point aj_point = fv_point(aj_v);

                // a_j on C₁ (S₁ walker); s_hit on C₂ (oracle2 target).
                if (aj_point.x < r.s_hit.x)
                    state.chords.push_back({aj_v.y, aj_v.edge, aj_v.side,
                        r.s_hit.edge, r.s_hit.side,
                        /*left_on_c1=*/true, /*right_on_c1=*/false});
                else
                    state.chords.push_back({aj_v.y, r.s_hit.edge,
                        r.s_hit.side, aj_v.edge, aj_v.side,
                        /*left_on_c1=*/false, /*right_on_c1=*/true});

                // [C91 §3.1 tex 224]: a_j sees ∂C₂ now → mark its S₁
                // chord superseded (junction companions have NONE).
                if (aj_v.chord_idx != NONE) {
                    assert(aj_v.chord_idx < state.invalidated_self.size());
                    state.invalidated_self[aj_v.chord_idx] = true;
                }

                state.p = aj_point;
                state.p_edge = aj_v.edge;
                state.p_side = aj_v.side;
                k = j + 1;
                break;
            }
            // [C91 §3.1 tex 206 case (ii)]: record q' → p'; p = p';
            // current = R's neighbor across chord ab; k = j.
            if (auto r = case_ii_test(j); r.fires) {
                const Chord& chord_ab = S2.chord(r.chord_idx);
                std::size_t q_edge = r.q_is_left ? chord_ab.left_edge
                                                  : chord_ab.right_edge;
                Side q_side = r.q_is_left ? chord_ab.left_side
                                          : chord_ab.right_side;
                Point q_point = s2_endpoint_point(q_edge, chord_ab.y);
                Point p_prime{r.p_prime_hit.x, r.p_prime_hit.y, NONE};
                SymbolicY chord_y = chord_ab.symbolic_y();

                // q on C₂ (S₂ exit chord endpoint); p' on C₁ (A_j of S₁).
                if (q_point.x < p_prime.x)
                    state.chords.push_back({chord_y, q_edge, q_side,
                        r.p_prime_hit.edge, r.p_prime_hit.side,
                        /*left_on_c1=*/false, /*right_on_c1=*/true});
                else
                    state.chords.push_back({chord_y, r.p_prime_hit.edge,
                        r.p_prime_hit.side, q_edge, q_side,
                        /*left_on_c1=*/true, /*right_on_c1=*/false});

                // [C91 §3.1 tex 224]: q sees p' on A_j now → mark S₂'s
                // chord_ab superseded.
                assert(r.chord_idx < state.invalidated_other.size());
                state.invalidated_other[r.chord_idx] = true;

                state.p = p_prime;
                state.p_edge = r.p_prime_hit.edge;
                state.p_side = r.p_prime_hit.side;
                state.s2_region = (chord_ab.region[0] == R)
                                    ? chord_ab.region[1]
                                    : chord_ab.region[0];
                // [C91 §3.1 tex 199 + §2 tex 47 SoS]: k = j is correct
                // only if p' is not an S₁ chord endpoint (otherwise the
                // p ≠ a_k tie-break would set k = j+1 instead).  p' is
                // a back-shot hit on A_j's interior, which under SoS
                // is disjoint from the {a_1..a_m} enumeration.
#ifdef CHAZELLE_EXPENSIVE_ASSERTS
                {
                    SymbolicY pp_y{state.p.y, NONE};
                    for (std::size_t l = 1; l + 1 < state.sequence.size(); ++l) {
                        const auto& v = state.sequence[l];
                        if (v.edge == state.p_edge && v.side == state.p_side)
                            assert(!symbolic_y_equal(v.y, pp_y) &&
                                   "[C91 §3.1 tex 199 + §2 tex 47]: "
                                   "p' must not coincide with any a_l");
                    }
                }
#endif
                k = j;
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
    // endpoints in cw ∂C₁ order.  [C91 §2.4(iii) tex 138]: arc-sequence is in
    // ∂C order (LEFT ascending, then RIGHT descending).  cw ∂C₁ from a₀
    // is RIGHT then LEFT; we counting-sort endpoints by cw-arc position
    // in O(m).  Each endpoint is associated with one arc via adj_arcs:
    // vertex endpoint ⟹ the single adj arc; mid-edge endpoint ⟹ the
    // "starting" adj arc whose start-y (derived from its bounding chord
    // per [C91 §2.4 tex 133]) matches chord.y.

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

    // [C91 §2.4 tex 133]: mid-edge endpoint — of the two adj arcs,
    // exactly one starts at the chord (its start-y, derived from the
    // bounding chord, matches the chord's y).
    auto starting_arc = [&](const Chord::AdjArcs& adj,
                             SymbolicY chord_y) -> std::size_t {
        assert(adj.count == 2);
        if (symbolic_y_equal(S.arc_start_symbolic_y(adj.arcs[0], C), chord_y))
            return adj.arcs[0];
        assert(symbolic_y_equal(S.arc_start_symbolic_y(adj.arcs[1], C), chord_y) &&
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

    // Inverse: arc_for_seq_pos[j] = S₁ arc index for A_j.  O(num_arcs)
    // one-time cost so the case (ii) per-test lookup is O(1).
    state.arc_for_seq_pos.assign(state.sequence.size(), NONE);
    for (std::size_t cw_pos = 0; cw_pos < num_arcs; ++cw_pos) {
        std::size_t j = state.arc_starts[cw_pos];
        if (j >= state.arc_for_seq_pos.size()) continue;
        std::size_t arc_idx = (cw_pos < num_right)
            ? (lrb + cw_pos) : (cw_pos - num_right);
        state.arc_for_seq_pos[j] = arc_idx;
    }

    // [C91 §3.1 tex 209]: sequence is "a_0, a_1, ..., a_m, a_{m+1}" where
    // a_1..a_m are S's EXIT-chord endpoints (2 per non-null chord, since
    // null-length chords are skipped above per tex 224) and a_0/a_{m+1}
    // are the junction companions.  Lemma 2.3 bounds m at O(n/γ + 1) —
    // γ isn't in scope here, but the structural identity is checkable.
    {
        std::size_t exit_chord_count = 0;
        for (std::size_t ci = 0; ci < S.num_chords(); ++ci)
            if (!S.chord(ci).dead && !S.chord(ci).is_null_length)
                ++exit_chord_count;
        assert(state.sequence.size() == 2 * exit_chord_count + 2 &&
               "[C91 §3.1 tex 209]: |sequence| = 2·(#exit chords) + 2");
    }
}

// ── rebuild_submap ──────────────────────────────────────────────

// [C91 §3.1 tex 226]: set up S in normal form from the [tex 224] chord
// inventory — sort endpoints along ∂C (edge name then y), build the
// arc-sequence table + region tree (parenthesis sweep), fill in each
// chord's adjacent-arc pointers per [§2.4(ii) tex 137], skip the tree
// decomposition (§3.2's job).
//
// Time: O((n₁/γ₁ + n₂/γ₂ + 1) log(n₁+n₂)) — endpoint sort dominates.
void rebuild_submap(Submap& out_S,
                     const Polygon& C,
                     const Submap& S1, const Polygon& C1,
                     const Submap& S2, const Polygon& C2,
                     const FusionState& state1,
                     const FusionState& state2) {
    assert(out_S.num_nodes() == 0 && out_S.num_arcs() == 0 &&
           out_S.num_chords() == 0 &&
           "[C91 §3.1 tex 226]: rebuild_submap requires a fresh output submap");
    assert(C.num_edges() == C1.num_edges() + C2.num_edges() &&
           "[C91 §3 tex 160]: C = C₁ ∪ C₂ shares one vertex");

    const std::size_t n_c1_edges = C1.num_edges();
    const std::size_t n_edges    = C.num_edges();

    // ── Step 1: Chord inventory (per [C91 §3.1 tex 224]) ────────

    struct Pending {
        SymbolicY y;
        std::size_t left_edge_c;
        Side        left_side;
        std::size_t right_edge_c;
        Side        right_side;
        bool        is_null_length;
    };
    std::vector<Pending> pending;
    pending.reserve(state1.chords.size() + state2.chords.size() +
                    S1.num_live_chords() + S2.num_live_chords());

    auto edge_in_c = [&](std::size_t edge, bool on_c1) {
        return on_c1 ? edge : (n_c1_edges + edge);
    };

    auto ingest_discovered = [&](const FusionState& st) {
        for (const auto& dc : st.chords) {
            pending.push_back({dc.y,
                edge_in_c(dc.left_edge,  dc.left_on_c1),  dc.left_side,
                edge_in_c(dc.right_edge, dc.right_on_c1), dc.right_side,
                /*is_null_length=*/false});
        }
    };
    ingest_discovered(state1);
    ingest_discovered(state2);

    // [C91 §3.1 tex 224]: drop old chords whose endpoint sat on the C₁/C₂
    // junction wrap.  Per [C91 §2.1 tex 72 case 3], the junction is a
    // C-endpoint in V(Cᵢ) with two duplicate companions; in V(C) the
    // junction is interior and the wrap is gone, so the chord is no
    // longer maximal there.
    const std::size_t junction_vidx_in_c = C1.num_vertices() - 1;
    const std::size_t junction_tag = C1.vertex(junction_vidx_in_c).index;
    assert(junction_tag == C2.vertex(0).index &&
           "[C91 §3 tex 160]: junction shared between C₁ and C₂");

    assert(junction_vidx_in_c >= 1 && junction_vidx_in_c + 1 < C.num_vertices() &&
           "[C91 §3 tex 160]: junction is interior to C");

    // [C91 §3.1 tex 224] visibility filter: union the per-pass bitsets.
    // state₁'s self = S₁, other = S₂; state₂ swaps.
    auto bit_at = [](const std::vector<bool>& b, std::size_t i) {
        return i < b.size() && b[i];
    };
    auto s1_invalid = [&](std::size_t i) {
        return bit_at(state1.invalidated_self, i) ||
               bit_at(state2.invalidated_other, i);
    };
    auto s2_invalid = [&](std::size_t i) {
        return bit_at(state1.invalidated_other, i) ||
               bit_at(state2.invalidated_self, i);
    };

    auto ingest_old = [&](const Submap& S, bool on_c1,
                           std::size_t junction_edge, auto is_invalid) {
        for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
            const Chord& c = S.chord(ci);
            if (c.dead) continue;
            // [C91 §3.1 tex 224]: drop junction-wrap chord ([§2.1 tex 72
            // case 3] gone in V(C)) or one superseded by case (i)/(ii).
            if (!c.is_null_length && c.y_tag == junction_tag &&
                (c.left_edge == junction_edge ||
                 c.right_edge == junction_edge)) continue;
            if (is_invalid(ci)) continue;
            pending.push_back({c.symbolic_y(),
                edge_in_c(c.left_edge,  on_c1), c.left_side,
                edge_in_c(c.right_edge, on_c1), c.right_side,
                c.is_null_length});
        }
    };
    ingest_old(S1, /*on_c1=*/true,  n_c1_edges - 1, s1_invalid);
    ingest_old(S2, /*on_c1=*/false, 0,             s2_invalid);

    // [C91 §2.1 tex 72 case 2 + §3.1 tex 224]: y-extremum junction —
    // synthesize the inside-pair null-length chord.  Convention: store
    // both endpoints at (next-edge, LEFT) per the s2.4 layout; auxiliary
    // SoS tag = C.num_vertices() (distinct from any polygon-vertex tag).
    if (is_local_y_extremum(C.vertex(junction_vidx_in_c - 1),
                             C.vertex(junction_vidx_in_c),
                             C.vertex(junction_vidx_in_c + 1))) {
        Pending p;
        p.y.y = C.vertex(junction_vidx_in_c).y;
        p.y.tag = C.num_vertices();
        p.left_edge_c = p.right_edge_c = n_c1_edges;
        p.left_side = p.right_side = LEFT;
        p.is_null_length = true;
        pending.push_back(p);
    }

    // ── Step 2: Sort endpoints along ∂C ────────────────────────

    // [C91 §3.1 tex 226]: sort by edge name, then by y within an edge.
    // [C91 §2.4(iii) tex 138] canonical ∂C order: LEFT ascending (0..n-1)
    // then RIGHT descending; within an edge the y-tiebreak follows the
    // cw traversal direction.
    auto trav_pos = [&](std::size_t edge, Side side) -> std::size_t {
        return (side == LEFT) ? edge : (2 * n_edges - 1 - edge);
    };
    auto edge_ascends = [&](std::size_t edge) -> bool {
        const auto& e = C.edge(edge);
        return symbolic_y_less(symbolic_y_of(C.vertex(e.start_idx)),
                                symbolic_y_of(C.vertex(e.end_idx)));
    };

    struct Endpoint {
        std::size_t edge_c;
        Side        side;
        SymbolicY   y;
        std::size_t pending_idx;
        bool        is_left_slot;
    };
    std::vector<Endpoint> eps;
    eps.reserve(2 * pending.size());
    for (std::size_t i = 0; i < pending.size(); ++i) {
        const auto& p = pending[i];
        eps.push_back({p.left_edge_c,  p.left_side,  p.y, i, true});
        eps.push_back({p.right_edge_c, p.right_side, p.y, i, false});
    }
    std::sort(eps.begin(), eps.end(),
        [&](const Endpoint& a, const Endpoint& b) {
            std::size_t pa = trav_pos(a.edge_c, a.side);
            std::size_t pb = trav_pos(b.edge_c, b.side);
            if (pa != pb) return pa < pb;
            bool trav_asc = (a.side == LEFT) == edge_ascends(a.edge_c);
            return trav_asc ? symbolic_y_less(a.y, b.y)
                            : symbolic_y_greater(a.y, b.y);
        });

    // ── Step 3: Walk endpoints → arcs + regions ─────────────────

    // [C91 §2.4(ii) tex 137 + §2.2 tex 94]: mid-edge endpoint → 2 adj
    // (arc splits); polygon-vertex endpoint → 1 adj (no split).
    // Identification via SoS y match per [§2 tex 47].
    auto is_vertex_endpoint = [&](std::size_t edge_c, SymbolicY y) -> bool {
        const auto& e = C.edge(edge_c);
        return symbolic_y_equal(y, symbolic_y_of(C.vertex(e.start_idx))) ||
               symbolic_y_equal(y, symbolic_y_of(C.vertex(e.end_idx)));
    };

    struct ChordOut {
        std::size_t r_outer = NONE;
        std::size_t r_inner = NONE;
        Chord::AdjArcs left_adj{};
        Chord::AdjArcs right_adj{};
        bool first_seen = false;
        // Mid-edge endpoint: the second adjacent arc is the one we emit
        // next.  Set here at the first visit, filled in once that next
        // arc exists.  Vertex endpoint: unused.
        std::size_t pending_after_left  = NONE;
        std::size_t pending_after_right = NONE;
        // [C91 §2.4 tex 133]: for null-length chords, the null-arc we
        // emit in R_inner.  Recorded at the opening visit so that nested
        // null-length chords at the same y-extremum don't accidentally
        // reuse each other's null-arc as their closing adjacent arc.
        std::size_t r_inner_arc = NONE;
    };
    std::vector<ChordOut> chord_out(pending.size());

    // List of chords waiting for their second adjacent arc.  Cleared
    // every time a new arc is emitted, so its size is at most the size
    // of one endpoint group — total cost across the walk is O(N).
    std::vector<std::size_t> pending_chords;

    std::vector<Arc> arcs;
    arcs.reserve(2 * pending.size() + 2);

    const std::size_t r_start = out_S.add_node();
    std::size_t current_region = r_start;

    struct Cursor { std::size_t edge; Side side; };
    Cursor cursor{0, LEFT};

    auto emit_arc = [&](std::size_t end_edge, Side end_side,
                         std::size_t override_edge_count = NONE) {
        assert(cursor.side == end_side &&
               "[C91 §2.4(iii)]: arc must not straddle the L/R wrap");
        Arc a;
        a.first_edge = cursor.edge; a.first_side = cursor.side;
        a.last_edge  = end_edge;     a.last_side  = end_side;
        a.region_node = current_region;
        if (override_edge_count != NONE) {
            a.edge_count = override_edge_count;
        } else {
            std::size_t lo = std::min(a.first_edge, a.last_edge);
            std::size_t hi = std::max(a.first_edge, a.last_edge);
            a.edge_count = C.count_nonnull_edges(lo, hi);
        }
        std::size_t idx = arcs.size();
        arcs.push_back(a);
        cursor = {end_edge, end_side};
        return idx;
    };

    // Fill in the waiting mid-edge endpoints' second adjacent arc with
    // the arc we just emitted.
    auto patch_after_arcs = [&](std::size_t after_arc) {
        for (std::size_t pi : pending_chords) {
            ChordOut& co = chord_out[pi];
            if (co.pending_after_left != NONE) {
                assert(co.left_adj.count == 1);
                co.left_adj.arcs[1] = after_arc;
                co.left_adj.count   = 2;
                co.pending_after_left = NONE;
            }
            if (co.pending_after_right != NONE) {
                assert(co.right_adj.count == 1);
                co.right_adj.arcs[1] = after_arc;
                co.right_adj.count   = 2;
                co.pending_after_right = NONE;
            }
        }
        pending_chords.clear();
    };

    const std::size_t end_v_edge = n_edges - 1;

    auto process_group = [&](std::size_t i, std::size_t j) {
        // [tex 226]: emit the arc closing at this group's position.
        std::size_t before_arc = emit_arc(eps[i].edge_c, eps[i].side);
        patch_after_arcs(before_arc);

        // [C91 §2.2] tree + cw walk → proper parenthesis nesting.
        for (std::size_t k = i; k < j; ++k) {
            const Endpoint& e = eps[k];
            ChordOut& co = chord_out[e.pending_idx];
            bool is_vert = is_vertex_endpoint(e.edge_c, e.y);

            Chord::AdjArcs& adj = e.is_left_slot ? co.left_adj : co.right_adj;
            std::size_t& pending_after = e.is_left_slot
                ? co.pending_after_left : co.pending_after_right;

            if (!co.first_seen) {
                // Opening: most-recent arc is in R_outer (pre-toggle).
                adj.arcs[0] = arcs.size() - 1;
                adj.count   = 1;
                co.first_seen = true;
                co.r_outer = current_region;
                co.r_inner = out_S.add_node();
                current_region = co.r_inner;

                // [C91 §2.4 tex 133]: synthesize R_inner-bounding null-arc.
                if (pending[e.pending_idx].is_null_length) {
                    co.r_inner_arc = emit_arc(
                        e.edge_c, e.side, /*edge_count=*/0);
                    assert(is_vert &&
                           "[C91 §2.1 tex 72]: null-length endpoints are "
                           "polygon-vertex companions");
                    continue;
                }
            } else {
                // Closing: null-length uses own R_inner-arc; else most-recent.
                adj.arcs[0] = pending[e.pending_idx].is_null_length
                    ? co.r_inner_arc : arcs.size() - 1;
                assert(adj.arcs[0] != NONE);
                adj.count   = 1;
                current_region = co.r_outer;
            }

            if (!is_vert) {
                pending_after = adj.arcs[0];
                pending_chords.push_back(e.pending_idx);
            }
        }
    };

    std::size_t i = 0;
    bool wrap_done = false;
    while (i < eps.size()) {
        // L→R wrap at end_vertex on first RIGHT endpoint.
        if (!wrap_done && eps[i].side == RIGHT) {
            std::size_t a_idx = emit_arc(end_v_edge, LEFT);
            patch_after_arcs(a_idx);
            cursor = {end_v_edge, RIGHT};
            wrap_done = true;
            continue;
        }

        std::size_t j = i + 1;
        while (j < eps.size() &&
               trav_pos(eps[j].edge_c, eps[j].side) ==
                   trav_pos(eps[i].edge_c, eps[i].side) &&
               symbolic_y_equal(eps[j].y, eps[i].y)) ++j;
        process_group(i, j);
        i = j;
    }

    // Trailing wrap (if no RIGHT endpoint) + final arc to vertex(0).
    if (!wrap_done) {
        std::size_t a_idx = emit_arc(end_v_edge, LEFT);
        patch_after_arcs(a_idx);
        cursor = {end_v_edge, RIGHT};
    }
    std::size_t tail_arc = emit_arc(0, RIGHT);
    patch_after_arcs(tail_arc);

    // [C91 §3.1]: closed-loop ∂C + Jordan curve ⟹ sweep balances.
    assert(current_region == r_start &&
           "[C91 §3.1]: parenthesis sweep must return to the starting region");

    // ── Step 4: Hand to Submap in canonical order ──────────────
    // [C91 §2.4(iii) tex 138]: LEFT ascending then RIGHT descending.
    // Submap::add_arc asserts the order at insertion.

    std::vector<std::size_t> left_order, right_order;
    left_order.reserve(arcs.size());
    right_order.reserve(arcs.size());
    for (std::size_t ai = 0; ai < arcs.size(); ++ai) {
        (arcs[ai].first_side == LEFT ? left_order : right_order)
            .push_back(ai);
    }
    std::stable_sort(left_order.begin(), left_order.end(),
        [&](std::size_t a, std::size_t b) {
            return arcs[a].first_edge < arcs[b].first_edge;
        });
    std::stable_sort(right_order.begin(), right_order.end(),
        [&](std::size_t a, std::size_t b) {
            return arcs[a].first_edge > arcs[b].first_edge;
        });

    std::vector<std::size_t> arc_remap(arcs.size(), NONE);
    auto add_in_order = [&](const std::vector<std::size_t>& order) {
        for (std::size_t old_ai : order)
            arc_remap[old_ai] = out_S.add_arc(arcs[old_ai]);
    };
    add_in_order(left_order);
    add_in_order(right_order);

    out_S.start_vertex = 0;
    out_S.end_vertex   = C.num_vertices() - 1;

    assert(!left_order.empty() &&
           "[C91 §2.4(iii) tex 138]: normal-form requires ≥1 LEFT arc");
    out_S.start_arc = arc_remap[left_order.front()];
    out_S.end_arc   = arc_remap[left_order.back()];

    auto remap_adj = [&](Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            assert(adj.arcs[k] != NONE);
            adj.arcs[k] = arc_remap[adj.arcs[k]];
        }
    };
    for (std::size_t pi = 0; pi < pending.size(); ++pi) {
        const Pending& p = pending[pi];
        ChordOut& co = chord_out[pi];
        assert(co.first_seen && co.r_outer != NONE && co.r_inner != NONE &&
               "[C91 §3.1]: every chord must be visited by the walk");
        assert(co.left_adj.count + co.right_adj.count >= 2 &&
               co.left_adj.count + co.right_adj.count <= 4 &&
               "[C91 §2.4(ii) tex 137]: chord has 2, 3, or 4 adj arcs total");

        Chord c;
        c.region[0] = co.r_outer; c.region[1] = co.r_inner;
        c.left_edge = p.left_edge_c;   c.left_side  = p.left_side;
        c.right_edge = p.right_edge_c; c.right_side = p.right_side;
        c.y = p.y.y; c.y_tag = p.y.tag;
        c.is_null_length = p.is_null_length;
        c.left_adj = co.left_adj; c.right_adj = co.right_adj;
        remap_adj(c.left_adj); remap_adj(c.right_adj);
        out_S.add_chord(c);
    }
}

} // namespace chazelle
