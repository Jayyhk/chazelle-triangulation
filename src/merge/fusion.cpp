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

    // [C91 §2.2 tex 96]: every live arc ends at a live chord endpoint of
    // its region and is that chord's recorded before-arc, so the slot
    // walk above is exhaustive — wrap-spanning arcs included, since they
    // are single chord-bounded structures ([C91 §2.4 tex 142]).  The one
    // exception is the chordless submap's single closed arc, reachable
    // through the C-endpoint pointers ([C91 §2.4(iii) tex 138]).
    if (S.node(region).incident_chords.empty()) {
        assert(S.num_live_chords() == 0 &&
               "[C91 §2.2 tex 102]: a chord-free region exists only in "
               "the chordless (single-region) submap");
        assert(S.start_arc != NONE && S.start_arc == S.end_arc &&
               S.start_arc < S.num_arcs() && !S.arc(S.start_arc).dead &&
               S.arc(S.start_arc).region_node == region &&
               "[C91 §2.4(iii) tex 138]: the chordless submap's closed "
               "arc is the endpoint arc");
        out.push(S.start_arc);
    }

    // [C91 §3.1 tex 181 + §2.2 tex 96]: the boundary alternates arcs and
    // exit chords, so the arc-structure count equals the region's degree
    // — at most 4 ([C91 §2.3 tex 114]); 1 for the chordless region.
    assert(out.count == std::max<std::size_t>(S.node(region).degree(), 1) &&
           "[C91 §2.2 tex 96]: arc count == region degree (≥ 1)");

    return out;
}

// ── local_shoot ─────────────────────────────────────────────────

RayHit local_shoot(Point p, Side direction,
                    std::size_t region,
                    const Submap& S, const Polygon& C,
                    const RayShootingOracle& oracle,
                    bool require_hit) {
    // [C91 §3.1 tex 181]: check each region arc, take the nearest hit.
    // Conformality bounds the loop to ≤ 4 arcs ([C91 §2.3 tex 114];
    // wrap-spanning arcs are single structures, [C91 §2.4 tex 142]).
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
        assert_subarc_clockwise(sub);

        RayHit hit = oracle.shoot(p, direction, ai, sub);
        if (!hit.hit) continue;

        // [C91 §2.1 tex 70]: a ray that misses everything "wraps around
        // in the spherical plane until it hits C again" — hits carry the
        // wrap flag, and a through-infinity hit lies BEHIND the source
        // (signed distance ≤ 0).  Assert the convention rather than a
        // forward-only report.
        double hit_signed_dist = (direction == LEFT)
            ? (p.x - hit.x) : (hit.x - p.x);
        assert(((hit_signed_dist <= 0.0) == hit.wrapped) &&
               "[C91 §2.1 tex 70]: wrapped ⟺ the hit lies at or behind "
               "the source in the travel direction");

        hit.hit_arc_idx = ai;  // propagate O(1) context

        if (!best.hit) {
            best = hit;
        } else if (hit.wrapped != best.wrapped) {
            // [C91 §2.1 tex 70]: every direct hit precedes every
            // through-infinity hit in ray order.
            if (!hit.wrapped) best = hit;
        } else {
            // Nearest in the shooting direction (within the same wrap
            // class the ray order is ascending signed distance).
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
                // travel direction and the edge's orientation.  Under SoS +
                // exact arithmetic, true ties only arise from double-boundary
                // alternation on the SAME edge; the heuristic here also
                // covers floating-point coincidences across edges by giving
                // a deterministic answer keyed on hit.edge's orientation.
                assert(hit.edge < C.num_edges());
                const auto& e = C.edge(hit.edge);
                bool edge_ascending = symbolic_y_less(
                    symbolic_y_of(C.vertex(e.start_idx)),
                    symbolic_y_of(C.vertex(e.end_idx)));

                Side minus_x_face = edge_ascending ? LEFT : RIGHT;
                Side plus_x_face  = edge_ascending ? RIGHT : LEFT;

                // [C91 §2.1 tex 72]: a ray travels INTO the free region
                // (shooting_direction) and strikes the wall whose free
                // side faces it: traveling RIGHT it arrives from the −X
                // side and strikes the −X face; traveling LEFT the +X
                // face.  This holds at distance 0 too (e.g. the
                // null-length inside pair, whose siblings face each
                // other): the source's own wall faces AWAY from the ray
                // and is never struck, so no self-hit is possible.
                Side struck_first_face =
                    (direction == RIGHT) ? minus_x_face : plus_x_face;

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
    // [C91 §2.1 tex 72]: "any point of ∂C has a unique horizontal 'chord
    // direction' ... this direction points to the left of an observer
    // walking clockwise around ∂C."  [C91 §2.2 tex 96]: the clockwise
    // ∂C traversal walks each region's boundary COUNTERCLOCKWISE with
    // respect to the region — the free region is on the observer's
    // LEFT, so the chord direction points INTO the region the ∂C point
    // bounds (away from the curve).
    //
    // The canonical order ([C91 §2.4(iii) tex 138]) walks Side::LEFT
    // ascending — the observer climbs an ascending edge on its west
    // wall with the snake on the right — so LEFT of an ascending edge
    // is the west wall and shoots WEST (LEFT); the other three cases
    // follow by symmetry.
    assert(edge < C.num_edges());
    const auto& e = C.edge(edge);
    SymbolicY start_y = symbolic_y_of(C.vertex(e.start_idx));
    SymbolicY end_y = symbolic_y_of(C.vertex(e.end_idx));
    bool edge_ascending = symbolic_y_less(start_y, end_y);

    if (side == LEFT)
        return edge_ascending ? LEFT : RIGHT;
    else
        return edge_ascending ? RIGHT : LEFT;
}

// ── chord_runs_through_infinity ─────────────────────────────────

bool chord_runs_through_infinity(const Polygon& C, const Chord& c) {
    assert(!c.is_null_length &&
           "[C91 §2.2 tex 96]: null-length chords occupy a single point "
           "and run through nothing");
    double xl = edge_x_at_y(C, c.left_edge, c.symbolic_y());
    double xr = edge_x_at_y(C, c.right_edge, c.symbolic_y());
    assert(xl <= xr && "chord slots are ordered by ascending x");
    if (xl == xr) return true;      // outside duplicate pair at an extremum
    return shooting_direction(c.left_edge, c.left_side, C) == LEFT;
}

// ── arc/slot identification ─────────────────────────────────────

// [s2-adjacency-convention]: does `arc_idx` START at the chord's slot
// endpoint?  Tested on the FULL start position (edge + side + derived
// start-y, [C91 §2.4 tex 133]) — the y alone is ambiguous when the
// other adj arc begins at the chord's source vertex (same symbolic y,
// different ∂C position — e.g. a pocket arc running from the source
// vertex to the mid-edge endpoint).
static bool arc_starts_at_chord_slot(const Submap& S, const Polygon& C,
                                     const Chord& c, bool left_slot,
                                     std::size_t arc_idx) {
    const Arc& a = S.arc(arc_idx);
    std::size_t se = left_slot ? c.left_edge : c.right_edge;
    Side ss = left_slot ? c.left_side : c.right_side;
    return a.first_edge == se && a.first_side == ss &&
           symbolic_y_equal(S.arc_start_symbolic_y(arc_idx, C),
                            c.symbolic_y());
}

// ── fusion_startup ──────────────────────────────────────────────

std::size_t fusion_startup(FusionState& state,
                            const Submap& S1, const Polygon& C1,
                            const Submap& S2, const Polygon& C2,
                            const RayShootingOracle& oracle1,
                            const RayShootingOracle& oracle2) {
    // [C91 §3.1 tex 179]: sequence starts at a₀ and ends at a_{m+1} — the
    // junction companions delimiting the cw tour of the walked ∂C.
    assert(state.sequence.size() >= 2 &&
           "[C91 §3.1]: fusion sequence needs at least a₀ and a_{m+1}");
    const bool at_end = state.junction_at_end;
    const std::size_t junction_edge = at_end ? C1.num_edges() - 1 : 0;
    const Side a0_side  = at_end ? RIGHT : LEFT;
    const Side am1_side = at_end ? LEFT : RIGHT;
    const auto& a0 = state.sequence[0];
    assert(a0.is_companion && a0.side == a0_side &&
           "[C91 §3.1 tex 179]: sequence[0] = a₀ (tour-start companion)");
    assert(a0.edge == junction_edge &&
           "[C91 §3.1 tex 179]: a₀ duplicates C₁ ∩ C₂ = the walked curve's "
           "endpoint, so a₀.edge = the junction-incident boundary edge");
    const auto& a_m1 = state.sequence.back();
    assert(a_m1.is_companion && a_m1.side == am1_side &&
           "[C91 §3.1 tex 179]: sequence.back() = a_{m+1} (tour-end companion)");
    assert(a_m1.edge == junction_edge &&
           "[C91 §3.1 tex 179]: a_{m+1} duplicates C₁ ∩ C₂, so "
           "a_{m+1}.edge = the junction-incident boundary edge");

    // [C91 §3.1]: "Using local shooting, we find the point of ∂C₁ that a₀
    // sees with respect to C₁."  a₀ is the companion vertex just past
    // the walked curve's junction-endpoint turnaround in clockwise order:
    // junction at the end → the RIGHT companion of C₁'s end vertex, on
    // the end-turn arc (end_arc, [C91 §2.4 tex 142]); at the start → the
    // LEFT companion of vertex 0, on the start-turn arc (start_arc).
    // If a chord endpoint sits exactly at that companion, the wrap arc
    // ENDS there and a₀'s clockwise-entered region is the SUCCESSOR
    // arc's ([C91 §3.1 tex 191]: elect the region "that we locally enter
    // as we leave p in a clockwise traversal of ∂C₁").
    std::size_t a0_arc;
    {
        std::size_t turn_arc = at_end ? S1.end_arc : S1.start_arc;
        assert(turn_arc != NONE && turn_arc < S1.num_arcs() &&
               !S1.arc(turn_arc).dead &&
               "[C91 §2.4(iii) tex 138]: S₁'s endpoint arcs must be set "
               "(normal form)");
        const Arc& ta = S1.arc(turn_arc);
        bool ends_at_a0 =
            ta.last_edge == a0.edge && ta.last_side == a0.side &&
            symbolic_y_equal(S1.arc_end_symbolic_y(turn_arc, C1), a0.y);
        a0_arc = ends_at_a0 ? (turn_arc + 1) % S1.num_arcs() : turn_arc;
    }
    std::size_t a0_region_s1 = S1.arc(a0_arc).region_node;
    Side a0_dir = shooting_direction(a0.edge, a0.side, C1);

    // a₀ is the junction vertex; use its coords directly.
    std::size_t junction_v = at_end ? C1.num_vertices() - 1 : 0;
    Point a0_point = C1.vertex(junction_v);

    // [C91 §3.1 tex 191]: "Technically, it is not quite true that a₀ is
    // always a point of ∂C₂.  It coincides with one most often, but
    // when it sits at a local extremum (in the y-direction) it is not
    // one because of duplication. ... when c₀ cannot see a point of
    // ∂C₂, an infinitesimal deformation of ∂C₂ locally around a₀ can
    // make c₀ see one. ... for simplicity we still go on saying that c₀
    // sees a point of ∂C₂ with respect to C."  When the junction is a
    // y-extremum of C and a₀ is one of its INSIDE duplicates, the point
    // of ∂C₂ that a₀ sees is its facing duplicate at distance zero —
    // the junction's null-length chord, which enters S structurally
    // ([C91 §3.1 tex 224]: "this last category is where null-length
    // chords originate"), not through ray-shooting.  Case (i) applies
    // with c₀ = the coincident duplicate and no new chord is recorded.
    bool junction_inside_pair = false;
    {
        // Neighbors of the junction in the MERGED curve's order.
        const Point& walker_nb = at_end ? C1.vertex(C1.num_vertices() - 2)
                                        : C1.vertex(1);
        const Point& target_nb = at_end ? C2.vertex(1)
                                        : C2.vertex(C2.num_vertices() - 2);
        const Point& prev_pt = at_end ? walker_nb : target_nb;
        const Point& next_pt = at_end ? target_nb : walker_nb;
        const Point& v_pt = a0_point;
        if (is_local_y_extremum(prev_pt, v_pt, next_pt)) {
            // Inside ("wedge") face of the walker's junction-incident
            // edge, via the branch x-slopes at the vertex ([C91 §2.1
            // tex 72]; same computation as rebuild_submap's junction
            // synthesis).
            const bool is_max = point_y_above(v_pt, prev_pt) &&
                                point_y_above(v_pt, next_pt);
            const double t_prev = is_max
                ? (prev_pt.x - v_pt.x) * (v_pt.y - next_pt.y)
                : (prev_pt.x - v_pt.x) * (next_pt.y - v_pt.y);
            const double t_next = is_max
                ? (next_pt.x - v_pt.x) * (v_pt.y - prev_pt.y)
                : (next_pt.x - v_pt.x) * (prev_pt.y - v_pt.y);
            assert(t_prev != t_next &&
                   "[C91 §2 tex 47]: junction branches must have distinct "
                   "x-slopes (equal slopes ⟹ overlapping edges)");
            const bool prev_left = t_prev < t_next;
            Side inside_wall;
            if (at_end) {
                // Walker edge = prev→v: ascends into a max, descends
                // into a min; its −x face is LEFT iff ascending.
                const Side minus_x = is_max ? LEFT : RIGHT;
                const Side plus_x  = is_max ? RIGHT : LEFT;
                inside_wall = prev_left ? plus_x : minus_x;
            } else {
                // Walker edge = v→next: descends from a max, ascends
                // from a min.
                const Side minus_x = is_max ? RIGHT : LEFT;
                const Side plus_x  = is_max ? LEFT : RIGHT;
                inside_wall = prev_left ? minus_x : plus_x;
            }
            junction_inside_pair = (a0.side == inside_wall);
        }
    }

    RayHit hit_c1{};
    if (!junction_inside_pair)
        hit_c1 = local_shoot(a0_point, a0_dir, a0_region_s1,
                             S1, C1, oracle1);

    // [C91 §3.1]: a₀ "most often" coincides with a ∂C₂ point but not at
    // y-extrema due to duplication.  SoS ([C91 §2 tex 47]) gives every vertex
    // a unique symbolic y, so visibility is well-defined either way.
    //
    // [C91 §3.1]: "the normal-form representation of S₂ ... lets us find,
    // in constant time, in which region of S₂ the point a₀ lies."  The
    // junction is at the OPPOSITE end of the target curve ([C91 §3 tex
    // 160]): the target's FIRST vertex when the walker's junction is at
    // its end, and its LAST vertex otherwise.
    std::size_t target_junction_arc = at_end ? S2.start_arc : S2.end_arc;
    assert(target_junction_arc != NONE &&
           "[C91 §2.4(iii) tex 138]: S₂'s C-endpoint arc pointer must be "
           "set (normal form)");
    std::size_t a0_region_s2 = S2.arc(target_junction_arc).region_node;

    RayHit hit_c2{};
    if (!junction_inside_pair)
        hit_c2 = local_shoot(a0_point, a0_dir, a0_region_s2,
                             S2, C2, oracle2);

    // [C91 §3.1 tex 185]: c₀ = the closer of hit_c1 and hit_c2.  Both are
    // guaranteed to hit (Lemma 2.1).  In the tex 191 inside-pair case
    // the shots are skipped: c₀ is the coincident junction duplicate on
    // ∂C₂ and case (i) applies directly.
    bool c0_on_c2 = true;
    RayHit c0{};
    if (!junction_inside_pair) {
    assert(hit_c1.hit && "[C91 §3.1]: a₀ must see ∂C₁ (Lemma 2.1)");
    assert(hit_c2.hit && "[C91 §3.1]: a₀ must see ∂C₂ (Lemma 2.1)");
    {
        double d1 = (a0_dir == LEFT)
            ? (a0_point.x - hit_c1.x) : (hit_c1.x - a0_point.x);
        double d2 = (a0_dir == LEFT)
            ? (a0_point.x - hit_c2.x) : (hit_c2.x - a0_point.x);
        // [C91 §2.1 tex 70]: ray order is lexicographic in
        // (wrapped, signed distance) — direct hits precede all
        // through-infinity ones.  A junction at a global y-extremum
        // makes a₀'s shot wrap on one or both curves.
        // [C91 §3.1 tex 191]: "for simplicity we still go on saying that
        // c₀ sees a point of ∂C₂ with respect to C."  The paper allows
        // edge cases (a₀ at a y-extremum where c₀ technically cannot see
        // ∂C₂) and explicitly mandates defaulting to "c₀ on ∂C₂"; ties
        // therefore select case (i).
        bool c2_at_or_before =
            (hit_c2.wrapped != hit_c1.wrapped) ? !hit_c2.wrapped
                                               : (d2 <= d1);
        if (c2_at_or_before) {
            c0 = hit_c2; c0_on_c2 = true;
        } else {
            c0 = hit_c1; c0_on_c2 = false;
        }
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

            // [C91 §2.1 tex 72]: junction has TWO target-frame companions,
            // each incident on one chord → up to two distinct S₂ chords
            // share a₀'s symbolic y.  Pick the SPECIFIC chord that a₀c₀
            // lies on: one endpoint at a₀'s position on a₀'s ∂C side.  In
            // the target's frame the junction sits at the OPPOSITE end
            // ([C91 §3 tex 160]) — edge 0 when the walker's junction is at
            // its end, the target's last edge otherwise.  (a₀'s side label
            // carries across the junction: the snake's left/right faces
            // are continuous through an interior vertex of C.)
            const std::size_t tj_edge = at_end ? 0 : C2.num_edges() - 1;
            const bool incident_on_a0 =
                (ch.left_edge  == tj_edge && ch.left_side  == a0.side) ||
                (ch.right_edge == tj_edge && ch.right_side == a0.side);
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
    // [C91 §3.1 tex 224]: flags are walker-frame — a₀ lives on the walked
    // curve; c₀ on the walked curve (case (ii)) or the target (case (i)).
    // In the tex 191 inside-pair case there is nothing to record: a₀c₀
    // is the junction's null-length chord, which rebuild_submap
    // synthesizes structurally (tex 224).
    if (!junction_inside_pair) {
        const bool a0_on_walker = true;
        const bool c0_on_walker = !c0_on_c2;
        if (a0_point.x < c0.x)
            state.chords.push_back({a0.y, a0.edge, a0.side, c0.edge,
                                    c0.side, a0_on_walker, c0_on_walker});
        else
            state.chords.push_back({a0.y, c0.edge, c0.side, a0.edge,
                                    a0.side, c0_on_walker, a0_on_walker});
    }

    if (c0_on_c2) {
        // [C91 §3.1 case (i)]: "set p = a₀, call the region of S₂ crossed by
        // a₀c₀ current."  Main loop starts at k=1 (A_k = arc a₀→a₁).
        state.p = a0_point;
        state.p_edge = a0.edge;
        state.p_side = a0.side;
        state.p_y = a0.y;

        // The cw step from a₀ walks the junction-incident edge away from
        // the junction: to vertex junction_v − 1 (RIGHT half first) when
        // the junction is at the end, to vertex 1 (LEFT half first) when
        // it is at the start.  Its y-direction is the topological
        // "leaving p" vector.
        std::size_t next_v = at_end ? junction_v - 1 : junction_v + 1;
        SymbolicY y_here = symbolic_y_of(C1.vertex(junction_v));
        SymbolicY y_next = symbolic_y_of(C1.vertex(next_v));
        bool a0_leaving_downward = symbolic_y_less(y_next, y_here);

        state.s2_region = resolve_s2_region(a0_region_s2, a0_leaving_downward);
        return 1;
    } else {
        // [C91 §3.1 case (ii)]: "skip all the way to c₀ ... set p = c₀, call
        // the region of S₂ containing a₀ current."  c₀ is on an edge
        // interior, not a named vertex.
        // [C91 §2 tex 47]: c₀ shares its perturbed y with a₀ (the
        // horizontal visibility chord); it has no vertex tag of its own,
        // so it borrows a₀'s tag for SoS comparisons.  The tag MUST ride
        // on Point.index (the oracle reads its SymbolicY as {p.y,
        // p.index}); leaving it NONE shifts the ray infinitesimally off
        // a₀'s y and makes a hit that lands exactly on a vertex at that
        // y (e.g. a local-minimum junction) slip past — tripping the
        // [C91 §3.1 tex 181] "local shoot must hit" invariant.
        state.p = Point{c0.x, c0.y, a0.y.tag};
        state.p_edge = c0.edge;
        state.p_side = c0.side;
        state.p_y = a0.y;
        assert(state.p.index == state.p_y.tag &&
               "[C91 §2 tex 47]: p's SoS tag must match its symbolic y");

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
        std::size_t n_edges = C1.num_edges();
        auto trav_pos = [&](std::size_t edge, Side side) -> std::size_t {
            // CW tour from the junction ([C91 §3.1 tex 179]): junction at
            // the end → RIGHT half (descending edges) then LEFT half
            // (ascending); junction at the start → LEFT then RIGHT.
            if (at_end)
                return (side == RIGHT) ? (n_edges - 1) - edge
                                       : n_edges + edge;
            return (side == LEFT) ? edge
                                  : 2 * n_edges - 1 - edge;
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

        // Arc-table index → cw tour key (must mirror
        // build_fusion_sequence's cw_pos for the same orientation): the
        // origin arc — the one double-backing around the junction
        // endpoint ([C91 §2.4 tex 142]) — is split by the tour origin
        // into a leading piece (key 0) and a trailing piece (last key);
        // c₀'s hit position selects the piece.
        std::size_t cw_key = 0;
        {
            const std::size_t N = S1.num_arcs();
            const std::size_t origin_arc = at_end ? S1.end_arc
                                                  : S1.start_arc;
            if (c0_arc_idx == origin_arc) {
                const Arc& oa = S1.arc(origin_arc);
                bool leading;
                if (at_end) {
                    if (oa.first_side == LEFT && oa.last_side == RIGHT)
                        leading = hit_c1.side == RIGHT;
                    else if (oa.first_side == LEFT)   // L→L double wrap
                        leading = hit_c1.side == RIGHT ||
                                  (hit_c1.side == LEFT &&
                                   hit_c1.edge <= oa.last_edge);
                    else                              // R→R double wrap
                        leading = hit_c1.side == RIGHT &&
                                  hit_c1.edge >= oa.last_edge;
                } else {
                    if (oa.first_side == RIGHT && oa.last_side == LEFT)
                        leading = hit_c1.side == LEFT;
                    else if (oa.first_side == LEFT)   // L→L double wrap
                        leading = hit_c1.side == LEFT &&
                                  hit_c1.edge <= oa.last_edge;
                    else                              // R→R double wrap
                        leading = hit_c1.side == LEFT ||
                                  (hit_c1.side == RIGHT &&
                                   hit_c1.edge >= oa.last_edge);
                }
                cw_key = leading ? 0 : N;
            } else if (at_end) {
                std::size_t lrb = S1.left_right_boundary();
                std::size_t num_right = N - lrb;
                cw_key = ((c0_arc_idx >= lrb)
                              ? (c0_arc_idx - lrb)
                              : (num_right + c0_arc_idx)) + 1;
            } else {
                cw_key = c0_arc_idx + 1;
            }
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

        // [C91 §3.1 tex 188]: skipped a_l's "exit chord endpoints see all
        // belong to ∂C₁ and thus are available directly from S₁."  No new
        // chord is recorded for them — their S₁ chords carry over via
        // rebuild_submap's `ingest_old` pass.
        for (std::size_t skipped_i = 1; skipped_i < lo; ++skipped_i) {
            if (!state.sequence[skipped_i].is_companion) {
                assert(state.sequence[skipped_i].chord_idx != NONE &&
                       "[C91 §3.1 tex 188]: skipped vertex must be an S₁ "
                       "chord endpoint (visibility available from S₁)");
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

    // [C91 §3.1 tex 179]: tour orientation (see FusionState).
    const bool at_end = state.junction_at_end;
    const std::size_t junction_v = at_end ? C1.num_vertices() - 1 : 0;

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
                return C1.vertex(junction_v);
            return Point{edge_x_at_y(C1, v.edge, v.y), v.y.y, NONE};
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

            // [C91 §3.1 tex 220]: "Whether a_j lies in R can be directly
            // inferred from the local orientation of the hit at s and
            // which side of the double boundary is hit."  hit.side must
            // match the arc's coverage at hit.edge — per-leg for wrapped
            // arcs ([C91 §2.4 tex 142]; "arc-structures encode on which
            // side(s) of the double boundary the arcs lie").
            const Arc& s_arc = S2.arc(s_hit.hit_arc_idx);
            assert(S2.start_vertex != NONE && S2.end_vertex != NONE &&
                   "[C91 §2.4(iii)]: S₂'s C endpoints must be set");
            if (!s_arc.covers(s_hit.edge, s_hit.side,
                              S2.start_vertex, S2.end_vertex))
                return {false, {}};

            // t = ∂C₁ point a_j sees.  O(1) at S₁ chord endpoints
            // (other endpoint); local_shoot in S₁ for companions.
            double t_x;
            bool t_wrapped;
            if (aj_v.chord_idx != NONE) {
                const Chord& ch = S1.chord(aj_v.chord_idx);
                std::size_t other_edge = aj_v.is_left_endpoint
                    ? ch.right_edge : ch.left_edge;
                t_x = edge_x_at_y(C1, other_edge, ch.symbolic_y());
                // [C91 §2.1 tex 70]: the chord itself may run through
                // infinity — then its other endpoint lies at or behind
                // a_j in the shooting direction.
                double t_signed = (dir == LEFT) ? (aj_point.x - t_x)
                                                : (t_x - aj_point.x);
                t_wrapped = (t_signed <= 0.0);
            } else {
                // [C91 §3.1 tex 220] companion case: j == m+1 only.
                // (j == 0 is consumed by fusion_startup, which always
                // returns k ≥ 1; main-loop j starts at k.)  a_{m+1} is
                // the companion just BEFORE the junction-endpoint
                // turnaround in clockwise order, so it lies on the
                // turn-crossing arc (end_arc / start_arc, [C91 §2.4 tex
                // 142]) — unless a chord endpoint sits exactly at the
                // companion, where that arc STARTS and a_{m+1} belongs
                // to the predecessor.  Per [C91 §2.1 tex 72] a_{m+1} and
                // a₀ can be in different S₁ regions whenever S₁'s
                // junction chord is live.
                assert(aj_v.is_companion);
                assert(aj_v.side == (at_end ? LEFT : RIGHT) &&
                       "[C91 §3.1]: only a_{m+1} reaches the companion branch "
                       "(a₀ is consumed by fusion_startup)");
                std::size_t s1_arc = at_end ? S1.end_arc : S1.start_arc;
                {
                    const Arc& sa = S1.arc(s1_arc);
                    bool starts_at_am1 =
                        sa.first_edge == aj_v.edge &&
                        sa.first_side == aj_v.side &&
                        symbolic_y_equal(
                            S1.arc_start_symbolic_y(s1_arc, C1), aj_v.y);
                    if (starts_at_am1)
                        s1_arc = (s1_arc + S1.num_arcs() - 1)
                                 % S1.num_arcs();
                }
                std::size_t aj_region_s1 = S1.arc(s1_arc).region_node;
                RayHit t_hit = local_shoot(aj_point, dir, aj_region_s1,
                                           S1, C1, oracle1);
                t_x = t_hit.x;
                t_wrapped = t_hit.wrapped;
            }

            // (i) fires iff s strictly precedes t in ray order —
            // lexicographic (wrapped, distance) per [C91 §2.1 tex 70].
            double s_dist = (dir == LEFT)
                ? (aj_point.x - s_hit.x) : (s_hit.x - aj_point.x);
            double t_dist = (dir == LEFT)
                ? (aj_point.x - t_x) : (t_x - aj_point.x);
            bool s_first = (s_hit.wrapped != t_wrapped)
                ? !s_hit.wrapped : (s_dist < t_dist);
            return {s_first, s_hit};
        };

        // Position of an S₂ chord endpoint: edge-interp at chord y.
        auto s2_endpoint_point = [&](std::size_t edge, SymbolicY y) -> Point {
            return Point{edge_x_at_y(C2, edge, y), y.y, NONE};
        };

        // CW position on ∂C₁ as (trav_pos, within_edge); lex compare
        // gives the cw walk.  LEFT side: within = edge param t; RIGHT
        // side: 1−t (cw traversal is end→start).  trav_pos mirrors
        // build_fusion_sequence's tour for the pass's orientation.
        auto cw_position = [&](SymbolicY y, std::size_t edge, Side side)
                           -> std::pair<std::size_t, double> {
            std::size_t n_edges = C1.num_edges();
            std::size_t tp;
            if (at_end)
                tp = (side == RIGHT) ? (n_edges - 1) - edge
                                     : n_edges + edge;
            else
                tp = (side == LEFT) ? edge
                                    : 2 * n_edges - 1 - edge;
            double t = edge_t_at_y(C1, edge, y);
            return {tp, (side == LEFT) ? t : (1.0 - t)};
        };

        // [C91 §3.1 tex 222]: case (ii) test in O(f(γ₁)).
        struct CaseIIResult {
            bool fires;
            RayHit p_prime_hit;
            std::size_t chord_idx;
            bool q_is_left;
        };
        // [C91 §3.1 tex 199]: A_j is the segment of ∂C₁ from a_{j-1} to
        // a_j cw.  A wrap-spanning A_j is ONE arc-structure that
        // double-backs ([C91 §2.4 tex 142]) — wraps never split it — but
        // a null-length chord strictly inside A_j (its endpoints are not
        // enumeration stops, [C91 §3.1 tex 224]) still separates two
        // structures.  Compute the up-to-two spanning nonzero structures
        // (arc-after a_{j-1}, arc-before a_j); they coincide unless a
        // null chord intervenes.
        auto Aj_arcs = [&](std::size_t j) -> std::array<std::size_t, 2> {
            // Clockwise successor of an arc in S₁'s arc-sequence table:
            // arcs tile ∂C in table order ([C91 §2.4(iii) tex 138]), the
            // circle closing from the last entry (the start-turn arc,
            // [C91 §2.4 tex 142]) back to index 0.
            auto cw_successor = [&](std::size_t ai) -> std::size_t {
                return (ai + 1) % S1.num_arcs();
            };
            // Tour-start companion side (see FusionState).
            const Side a0_side  = at_end ? RIGHT : LEFT;
            auto arc_after = [&](const FusionVertex& v) -> std::size_t {
                if (v.is_companion) {
                    if (v.side != a0_side) return NONE;  // a_{m+1}
                    // a₀ sits just past the junction turnaround, on the
                    // turn-crossing arc — or on its successor when a
                    // chord endpoint coincides with a₀.
                    std::size_t turn_arc = at_end ? S1.end_arc
                                                  : S1.start_arc;
                    const Arc& ta = S1.arc(turn_arc);
                    bool ends_at_v =
                        ta.last_edge == v.edge && ta.last_side == v.side &&
                        symbolic_y_equal(
                            S1.arc_end_symbolic_y(turn_arc, C1), v.y);
                    return ends_at_v ? cw_successor(turn_arc) : turn_arc;
                }
                const Chord& c = S1.chord(v.chord_idx);
                const Chord::AdjArcs& adj = v.is_left_endpoint
                    ? c.left_adj : c.right_adj;
                // [C91 §2.4(ii) tex 137]: a polygon-vertex endpoint records
                // ONE adj arc — the arc ENDING there (before-arc).  The arc
                // starting at the endpoint is its clockwise successor.
                if (adj.count == 1) return cw_successor(adj.arcs[0]);
                // count==2: the starting arc, by full start position.
                return arc_starts_at_chord_slot(S1, C1, c,
                                                v.is_left_endpoint,
                                                adj.arcs[0])
                    ? adj.arcs[0] : adj.arcs[1];
            };
            auto arc_before = [&](const FusionVertex& v) -> std::size_t {
                if (v.is_companion) {
                    if (v.side == a0_side) return NONE;  // a₀
                    // a_{m+1} sits just before the turnaround, on the
                    // turn-crossing arc — or on its predecessor when a
                    // chord endpoint coincides with a_{m+1}.
                    std::size_t turn_arc = at_end ? S1.end_arc
                                                  : S1.start_arc;
                    const Arc& ta = S1.arc(turn_arc);
                    bool starts_at_v =
                        ta.first_edge == v.edge &&
                        ta.first_side == v.side &&
                        symbolic_y_equal(
                            S1.arc_start_symbolic_y(turn_arc, C1), v.y);
                    return starts_at_v
                        ? (turn_arc + S1.num_arcs() - 1) % S1.num_arcs()
                        : turn_arc;
                }
                const Chord& c = S1.chord(v.chord_idx);
                const Chord::AdjArcs& adj = v.is_left_endpoint
                    ? c.left_adj : c.right_adj;
                // [C91 §2.4(ii) tex 137]: vertex endpoint's single adj arc
                // IS the arc ending there — exactly the arc-before.
                if (adj.count == 1) return adj.arcs[0];
                // count==2: arc-before is the one NOT starting at the slot.
                return arc_starts_at_chord_slot(S1, C1, c,
                                                v.is_left_endpoint,
                                                adj.arcs[0])
                    ? adj.arcs[1] : adj.arcs[0];
            };
            return std::array<std::size_t, 2>{
                arc_after(state.sequence[j - 1]),
                arc_before(state.sequence[j])};
        };

        auto case_ii_test = [&](std::size_t j) -> CaseIIResult {
            CaseIIResult best{false, {}, NONE, false};
            auto best_cw = std::make_pair<std::size_t, double>(0, -1.0);

            auto aj_arcs_pair = Aj_arcs(j);
            if (aj_arcs_pair[0] == NONE && aj_arcs_pair[1] == NONE)
                return best;

            auto p_cw = cw_position(state.p_y, state.p_edge, state.p_side);

            // [C91 §3.1 tex 200]: case (ii) iterates the exit chords of R.
            // Null-length chords are iterated too but can never fire: the
            // tex 222 "hit must lie on ab" disqualification below rejects
            // every hit against a zero-length ab.  That is exactly why
            // [C91 §3.1 tex 224] carries null-length chords over to the
            // fusion unconditionally — their endpoints' visibility is
            // settled at distance 0 ([C91 §2.1 tex 72] inside pair).
            for (std::size_t ci : S2.node(R).incident_chords) {
                const Chord& chord_ab = S2.chord(ci);
                if (chord_ab.dead) continue;

                // [C91 §2.1 tex 70]: whether ab itself runs through
                // infinity decides the shape of the on-ab test below.
                const bool ab_wraps =
                    !chord_ab.is_null_length &&
                    chord_runs_through_infinity(C2, chord_ab);

                SymbolicY chord_ab_y = chord_ab.symbolic_y();
                Point a_pt = s2_endpoint_point(chord_ab.left_edge,
                                                chord_ab_y);
                Point b_pt = s2_endpoint_point(chord_ab.right_edge,
                                                chord_ab_y);

                for (bool is_left : {true, false}) {
                    std::size_t q_edge = is_left ? chord_ab.left_edge
                                                  : chord_ab.right_edge;
                    Side q_side = is_left ? chord_ab.left_side
                                          : chord_ab.right_side;
                    Point q_point = is_left ? a_pt : b_pt;
                    Point other_point = is_left ? b_pt : a_pt;
                    Side shoot_dir = shooting_direction(q_edge, q_side, C2);

                    // [C91 §3.1 tex 222]: shoot toward A_j.  A wrapping
                    // A_j is ONE arc-structure and ONE oracle call with
                    // a double-backing subarc ([C91 §2.4 tex 142] /
                    // [C91 §3.0(i) tex 169]: subarcs carry side flags).
                    // Only a null-length chord strictly inside A_j
                    // splits it into two structures, each shot with its
                    // own restricted subarc — ≤ 2 calls stay O(f(γ₁)).
                    for (std::size_t arc_slot = 0; arc_slot < 2; ++arc_slot) {
                        std::size_t aj_arc = aj_arcs_pair[arc_slot];
                        if (aj_arc == NONE) continue;
                        if (arc_slot == 1 &&
                            aj_arcs_pair[0] == aj_arcs_pair[1])
                            continue;  // single-arc: no second call

                        const Arc& aj_arc_struct = S1.arc(aj_arc);
                        Subarc aj_sub;
                        if (aj_arcs_pair[0] == aj_arcs_pair[1]) {
                            // Single-arc A_j: full span from a_{j-1} to a_j.
                            aj_sub = Subarc{state.sequence[j-1].edge,
                                            state.sequence[j-1].side,
                                            state.sequence[j].edge,
                                            state.sequence[j].side};
                        } else if (arc_slot == 0) {
                            // First structure: a_{j-1} → its end (the
                            // interior null-chord position).
                            aj_sub = Subarc{state.sequence[j-1].edge,
                                            state.sequence[j-1].side,
                                            aj_arc_struct.last_edge,
                                            aj_arc_struct.last_side};
                        } else {
                            // Second structure: its start (just past the
                            // interior null-chord position) → a_j.
                            aj_sub = Subarc{aj_arc_struct.first_edge,
                                            aj_arc_struct.first_side,
                                            state.sequence[j].edge,
                                            state.sequence[j].side};
                        }

                        assert_subarc_clockwise(aj_sub);
                        RayHit hit = oracle1.shoot(q_point, shoot_dir,
                                                    aj_arc, aj_sub);
                        if (!hit.hit) continue;

                        // [C91 §2.1 tex 70]: hits are ordered in the
                        // lexicographic (wrapped, distance) ray metric;
                        // wrapped ⟺ at-or-behind the source.
                        double q_to_hit = (shoot_dir == LEFT)
                            ? (q_point.x - hit.x) : (hit.x - q_point.x);
                        assert(((q_to_hit <= 0.0) == hit.wrapped) &&
                               "[C91 §2.1 tex 70]: wrapped ⟺ the hit "
                               "lies at or behind the source");

                        // [C91 §3.1 tex 222]: "or if the hit does not lie
                        // on ab" — the ray from q travels along the chord
                        // ab's line; a hit beyond the other endpoint would
                        // make the sightline pass through that endpoint,
                        // a point of ∂C₂, so q cannot see it with respect
                        // to C.  For a null-length ab no point of ∂C₁
                        // lies on the zero-length segment, so every hit
                        // is rejected ([C91 §2.1 tex 72] inside pair).
                        // A through-infinity ab ([C91 §2.1 tex 70] —
                        // e.g. an outside duplicate pair at a global
                        // y-extremum, where q_to_other ≤ 0) contains
                        // every ray position up to the far endpoint in
                        // the wrap metric.
                        double q_to_other = (shoot_dir == LEFT)
                            ? (q_point.x - other_point.x)
                            : (other_point.x - q_point.x);
                        if (chord_ab.is_null_length) continue;
                        if (ab_wraps) {
                            // ab through infinity: its free span covers
                            // every ray position up to the far endpoint
                            // in the wrap metric — on-ab ⟺
                            // lex(hit) ≤ (wrapped, q_to_other).
                            if (hit.wrapped && q_to_hit > q_to_other)
                                continue;
                        } else {
                            // ab direct: on-ab ⟺ direct and not beyond.
                            if (hit.wrapped || q_to_hit > q_to_other)
                                continue;
                        }

                        // "proper orientation": hit.Side matches A_j's
                        // Side at hit.edge ([C91 §3.1 tex 222]).  A
                        // wrapping A_j covers per-leg side/edge zones
                        // ([C91 §2.4 tex 142]) — test the shot subarc's
                        // leg coverage directly.
                        assert(S1.start_vertex != NONE &&
                               S1.end_vertex != NONE &&
                               "[C91 §2.4(iii)]: S₁'s C endpoints must "
                               "be set");
                        if (!subarc_covers_position(aj_sub, hit.edge,
                                                    hit.side,
                                                    S1.start_vertex,
                                                    S1.end_vertex))
                            continue;

                        // "occurs before p along A_j": strict cw compare.
                        // hit.y inherits the ray's perturbed source y =
                        // chord_ab's SymbolicY (RayHit carries only raw y).
                        auto hit_cw = cw_position(chord_ab_y, hit.edge, hit.side);
                        if (hit_cw <= p_cw) continue;

                        // back-shot from s in its natural direction; q
                        // must lie at or before t in the ray order for
                        // s↔a visibility ([C91 §3.1 tex 222]: "If
                        // shooting from s to t passes through a").
                        // [C91 §2.1 tex 70] wrap metric: q at distance 0
                        // (a junction duplicate at s's own position) is
                        // seen immediately; otherwise compare
                        // lexicographically by (behind-source, distance).
                        Side s_back_dir = shooting_direction(hit.edge,
                                                             hit.side, C1);
                        Point s_point{hit.x, hit.y, NONE};
                        RayHit t_hit = local_shoot(
                            s_point, s_back_dir,
                            aj_arc_struct.region_node,
                            S1, C1, oracle1);
                        double s_to_q = (s_back_dir == LEFT)
                            ? (s_point.x - q_point.x)
                            : (q_point.x - s_point.x);
                        double s_to_t = (s_back_dir == LEFT)
                            ? (s_point.x - t_hit.x)
                            : (t_hit.x - s_point.x);
                        if (s_to_q != 0.0) {
                            bool q_behind = (s_to_q < 0.0);
                            bool q_first = (q_behind != t_hit.wrapped)
                                ? !q_behind : (s_to_q <= s_to_t);
                            if (!q_first) continue;
                        }

                        // [C91 §3.1 tex 206]: pick LAST p' cw from p.
                        if (hit_cw > best_cw) {
                            best = {true, hit, ci, is_left};
                            best_cw = hit_cw;
                        }
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

                // a_j on the walked curve; s_hit on the target (walker-frame flags).
                if (aj_point.x < r.s_hit.x)
                    state.chords.push_back({aj_v.y, aj_v.edge, aj_v.side,
                        r.s_hit.edge, r.s_hit.side,
                        /*left_on_walker=*/true, /*right_on_walker=*/false});
                else
                    state.chords.push_back({aj_v.y, r.s_hit.edge,
                        r.s_hit.side, aj_v.edge, aj_v.side,
                        /*left_on_walker=*/false, /*right_on_walker=*/true});

                // [C91 §3.1 tex 224]: a_j sees ∂C₂ now → mark its S₁
                // chord superseded (junction companions have NONE).
                if (aj_v.chord_idx != NONE) {
                    assert(aj_v.chord_idx < state.invalidated_self.size());
                    state.invalidated_self[aj_v.chord_idx] = true;
                }

                state.p = aj_point;
                state.p_edge = aj_v.edge;
                state.p_side = aj_v.side;
                state.p_y = aj_v.y;
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
                SymbolicY chord_y = chord_ab.symbolic_y();
                Point q_point = s2_endpoint_point(q_edge, chord_y);
                // p' lies mid-edge on A_j but its perturbed y is
                // chord_ab's (the back-shot's source); carry that SoS tag
                // on Point.index so oracle shots from p keep the correct
                // y (see the case (ii) startup note above).
                Point p_prime{r.p_prime_hit.x, r.p_prime_hit.y, chord_y.tag};

                // q on C₂ (S₂ exit chord endpoint); p' on C₁ (A_j of S₁).
                if (q_point.x < p_prime.x)
                    state.chords.push_back({chord_y, q_edge, q_side,
                        r.p_prime_hit.edge, r.p_prime_hit.side,
                        /*left_on_walker=*/false, /*right_on_walker=*/true});
                else
                    state.chords.push_back({chord_y, r.p_prime_hit.edge,
                        r.p_prime_hit.side, q_edge, q_side,
                        /*left_on_walker=*/true, /*right_on_walker=*/false});

                // [C91 §3.1 tex 224]: q sees p' on A_j now → mark S₂'s
                // chord_ab superseded.
                assert(r.chord_idx < state.invalidated_other.size());
                state.invalidated_other[r.chord_idx] = true;

                state.p = p_prime;
                state.p_edge = r.p_prime_hit.edge;
                state.p_side = r.p_prime_hit.side;
                // p' = back-shot hit on ∂C₁ from a chord_ab endpoint;
                // its perturbed y is chord_ab's (the ray's source).
                state.p_y = chord_y;
                // [C91 §3.1 tex 222]: a firing chord has nonzero length —
                // the on-ab disqualification rejects every hit against a
                // zero-length ab (see the comment at the candidate loop).
                assert(!chord_ab.is_null_length &&
                       "[C91 §3.1 tex 222]: null-length chords cannot fire "
                       "case (ii) — no hit lies on a zero-length ab");
                // [C91 §3.1 tex 206]: "make current the region of S₂
                // which we enter as we locally cross the exit chord" —
                // p' lies on ab, so crossing it at p' enters the chord's
                // other side.
                state.s2_region = (chord_ab.region[0] == R)
                    ? chord_ab.region[1] : chord_ab.region[0];
                // [C91 §3.1 tex 199]: update k per its definition — A_k is
                // the arc containing p with the tie-break p ≠ a_k ("we
                // must choose the one starting (and not ending) at p").
                // p' lies on A_j strictly after p, so among the
                // enumeration stops it can coincide only with a_j (the
                // arc's terminal stop); in that case A_k = A_{j+1}.
                // Coincidence requires equal symbolic y (identical SoS
                // tags), possible only when both chords are sourced at
                // the shared junction vertex ([C91 §3 tex 160]).
                {
                    const FusionVertex& aj_end = state.sequence[j];
                    bool p_at_aj =
                        state.p_edge == aj_end.edge &&
                        state.p_side == aj_end.side &&
                        symbolic_y_equal(state.p_y, aj_end.y);
                    k = p_at_aj ? j + 1 : j;
                }
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
    const bool at_end = state.junction_at_end;

    // [C91 §3.1 tex 179]: junction = C₁ ∩ C₂ in the walked curve's frame —
    // its LAST vertex for pass 1, its FIRST for the symmetric pass.  The
    // junction's incident edge is the walked curve's last / first edge.
    std::size_t junction_v    = at_end ? n - 1 : 0;
    std::size_t junction_edge = at_end ? n - 2 : 0;

    // [C91 §3.1 tex 179]: a₀ / a_{m+1} = junction companions at the start
    // / end of the cw ∂C tour.  Junction at the end: tour = RIGHT half
    // then LEFT half, so a₀ is the RIGHT companion.  Junction at the
    // start: tour = LEFT half then RIGHT half, a₀ is the LEFT companion.
    SymbolicY junction_y = symbolic_y_of(C.vertex(junction_v));

    FusionVertex a_0;
    a_0.y = junction_y;
    a_0.edge = junction_edge;
    a_0.side = at_end ? RIGHT : LEFT;
    a_0.chord_idx = NONE;
    a_0.is_left_endpoint = false;
    a_0.is_companion = true;

    FusionVertex a_m1;
    a_m1.y = junction_y;
    a_m1.edge = junction_edge;
    a_m1.side = at_end ? LEFT : RIGHT;
    a_m1.chord_idx = NONE;
    a_m1.is_left_endpoint = false;
    a_m1.is_companion = true;

    // [C91 §3.1]: a₁, ..., aₘ = canonical enumeration of the walked
    // submap's exit chord endpoints in cw ∂C order from the junction.
    // [C91 §2.4(iii) tex 138]: the arc-sequence table is already in ∂C
    // order (LEFT ascending, then RIGHT descending); we counting-sort
    // endpoints by cw-arc position in O(m).  Each endpoint is associated
    // with one arc via adj_arcs: vertex endpoint ⟹ the single adj arc;
    // mid-edge endpoint ⟹ the "starting" adj arc whose start-y (derived
    // from its bounding chord per [C91 §2.4 tex 133]) matches chord.y.

    std::size_t num_arcs = S.num_arcs();
    std::size_t lrb = S.left_right_boundary();
    std::size_t num_right = num_arcs - lrb;

    // Arc-table index → cw tour position.  Exactly ONE arc is split by
    // the tour's starting turnaround ([C91 §2.4 tex 142]: the arc
    // double-backing around the junction endpoint — end_arc when the
    // junction is at the walked curve's end, start_arc otherwise): its
    // LEADING piece (just past the turn) opens the tour at key 0 and
    // its TRAILING piece closes it at key num_arcs.  Everything else
    // keeps one key.  Key domain: [0, num_arcs].
    const std::size_t origin_arc = at_end ? S.end_arc : S.start_arc;
    assert(origin_arc != NONE && origin_arc < num_arcs &&
           "[C91 §2.4(iii) tex 138]: the walked submap's endpoint arcs "
           "must be set (normal form)");

    // Is a ∂C position on the origin arc part of its LEADING piece —
    // between the tour's turnaround and the arc's end, clockwise?
    // Enumerated per wrap class ([C91 §2.4 tex 142]).
    auto origin_leading = [&](std::size_t edge, Side side) -> bool {
        const Arc& a = S.arc(origin_arc);
        if (at_end) {
            // Tour origin = C's END turn.
            if (a.first_side == LEFT && a.last_side == RIGHT)
                return side == RIGHT;                       // end wrap
            assert(a.first_side == a.last_side && a.wraps() &&
                   "[C91 §2.4 tex 142]: at_end origin arc must wrap C's "
                   "end vertex");
            if (a.first_side == LEFT)                       // L→L double
                return side == RIGHT ||
                       (side == LEFT && edge <= a.last_edge);
            return side == RIGHT && edge >= a.last_edge;    // R→R double
        }
        // Tour origin = C's START turn.
        if (a.first_side == RIGHT && a.last_side == LEFT)
            return side == LEFT;                            // start wrap
        assert(a.first_side == a.last_side && a.wraps() &&
               "[C91 §2.4 tex 142]: at_start origin arc must wrap C's "
               "start vertex");
        if (a.first_side == LEFT)                           // L→L double
            return side == LEFT && edge <= a.last_edge;
        return side == LEFT ||
               (side == RIGHT && edge >= a.last_edge);      // R→R double
    };

    auto cw_pos = [&](std::size_t arc_idx, std::size_t edge,
                      Side side) -> std::size_t {
        if (arc_idx == origin_arc)
            return origin_leading(edge, side) ? 0 : num_arcs;
        std::size_t base = at_end
            ? ((arc_idx >= lrb) ? arc_idx - lrb : num_right + arc_idx)
            : arc_idx;
        return base + 1;
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
    // exactly one STARTS at the chord endpoint, identified by its full
    // start position (arc_starts_at_chord_slot; the y alone can tie).
    auto starting_arc = [&](const Chord& c, bool left_slot) -> std::size_t {
        const Chord::AdjArcs& adj = left_slot ? c.left_adj : c.right_adj;
        assert(adj.count == 2);
        bool s0 = arc_starts_at_chord_slot(S, C, c, left_slot, adj.arcs[0]);
        assert(s0 != arc_starts_at_chord_slot(S, C, c, left_slot,
                                              adj.arcs[1]) &&
               "[C91 §2.4]: exactly one adj arc starts at a mid-edge "
               "chord endpoint");
        return s0 ? adj.arcs[0] : adj.arcs[1];
    };

    for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
        assert(!S.chord(ci).dead && "[C91 §2.4]: normal-form submap has no dead chords");
        const auto& c = S.chord(ci);
        // [C91 §3.1 tex 179]: only EXIT chords enter the enumeration.
        // Null-length chords ([C91 §2.2 tex 96]) are carried over to the fused
        // submap directly (tex 224); including them would blow Lemma 3.1's
        // O(n/γ + 1) bound (tex 209–210).
        if (c.is_null_length) continue;

        // LEFT endpoint → associated arc.
        std::size_t left_arc = (c.left_adj.count == 2)
            ? starting_arc(c, true)
            : c.left_adj.arcs[0];
        endpoints.push_back({make_vertex(c, ci, true),
                             cw_pos(left_arc, c.left_edge, c.left_side)});

        // RIGHT endpoint → associated arc.
        std::size_t right_arc = (c.right_adj.count == 2)
            ? starting_arc(c, false)
            : c.right_adj.arcs[0];
        endpoints.push_back({make_vertex(c, ci, false),
                             cw_pos(right_arc, c.right_edge,
                                    c.right_side)});
    }

    // Counting sort by clockwise arc key.  O(m) time, O(m) space.
    // [C91 §2.4 tex 142]: "canonical vertex enumerations in optimal
    // time" — O(m) where m = number of chord endpoints.
    const std::size_t num_keys = num_arcs + 1;   // origin arc: 2 pieces
    std::vector<KeyedVertex> sorted(endpoints.size());
    std::vector<std::size_t> bucket_starts;
    if (!endpoints.empty()) {
        std::vector<std::size_t> counts(num_keys + 1, 0);
        for (const auto& ep : endpoints)
            ++counts[ep.key + 1];
        for (std::size_t i = 1; i <= num_keys; ++i)
            counts[i] += counts[i - 1];

        bucket_starts.assign(counts.begin(), counts.begin() + static_cast<std::ptrdiff_t>(num_keys));

        for (const auto& ep : endpoints)
            sorted[counts[ep.key]++] = ep;
    } else {
        bucket_starts.assign(num_keys, 0);
    }

    // Within each bucket (same arc): O(1) elements (conformal degree ≤ 4).
    // Insertion-sort by y along the arc's traversal direction; O(m) total.
    std::size_t n_edges = C.num_edges();
    auto trav_pos = [&](std::size_t edge, Side side) -> std::size_t {
        if (at_end)
            return (side == RIGHT) ? junction_edge - edge
                                   : n_edges + edge;
        return (side == LEFT) ? edge
                              : 2 * n_edges - 1 - edge;
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

    // Map cw key → start of its block in `sequence` (+1 for the
    // prepended a₀).  Key domain [0, num_arcs]: the origin arc's two
    // pieces bracket the tour ([C91 §2.4 tex 142]).
    state.arc_starts.resize(num_keys);
    for (std::size_t i = 0; i < num_keys; ++i)
        state.arc_starts[i] = bucket_starts[i] + 1;

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
// chord's adjacent-arc pointers per [C91 §2.4(ii) tex 137], skip the tree
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

    // [C91 §3.1 tex 224]: DiscoveredChord flags are WALKER-frame ("true"
    // = endpoint on the curve the pass walked).  state1 walked C₁,
    // state2 walked C₂ — translate to C's frame accordingly.
    auto ingest_discovered = [&](const FusionState& st, bool walker_is_c1) {
        for (const auto& dc : st.chords) {
            pending.push_back({dc.y,
                edge_in_c(dc.left_edge,  dc.left_on_walker == walker_is_c1),
                dc.left_side,
                edge_in_c(dc.right_edge, dc.right_on_walker == walker_is_c1),
                dc.right_side,
                /*is_null_length=*/false});
        }
    };
    ingest_discovered(state1, /*walker_is_c1=*/true);
    ingest_discovered(state2, /*walker_is_c1=*/false);

    // Junction identity, used by the asserts below and the null-chord
    // synthesis.  No dedicated "drop junction-wrap chords" step exists
    // or is needed: an old chord whose endpoint sat on the junction
    // wrap of V(Cᵢ) has endpoints that are enumeration stops of the
    // fusion passes ([C91 §3.1 tex 179]), so the general tex-224
    // visibility filter below invalidates it exactly when the other
    // curve blocks it in V(C).
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

    auto ingest_old = [&](const Submap& S, bool on_c1, auto is_invalid) {
        for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
            const Chord& c = S.chord(ci);
            if (c.dead) continue;
            // [C91 §2.1 tex 72 case 3]: the junction is a C-endpoint of
            // Cᵢ, and C-endpoints are never y-extrema, so V(Cᵢ) has NO
            // null-length chord sourced at the junction vertex.  (The
            // junction becomes case 2 only in the merged C; that chord is
            // synthesized below.)
            assert(!(c.is_null_length && c.y_tag == junction_tag) &&
                   "[C91 §2.1 tex 72 case 3]: Sᵢ cannot contain a "
                   "null-length chord sourced at its own C-endpoint");
            if (is_invalid(ci)) continue;
            pending.push_back({c.symbolic_y(),
                edge_in_c(c.left_edge,  on_c1), c.left_side,
                edge_in_c(c.right_edge, on_c1), c.right_side,
                c.is_null_length});
        }
    };
    ingest_old(S1, /*on_c1=*/true,  s1_invalid);
    ingest_old(S2, /*on_c1=*/false, s2_invalid);

    // [C91 §2.1 tex 72 case 2 + §3.1 tex 224]: y-extremum junction —
    // synthesize the inside-pair null-length chord.  The inside pair
    // duplicates the junction vertex itself, so its symbolic y is the
    // vertex's own (y, SoS index) ([C91 §2 tex 47] — no auxiliary
    // perturbation source exists).  Both endpoints stored at C's edge
    // `junction_vidx_in_c` (the first C₂ edge) on the wedge ("inside of
    // the turn") side, computed from the two incident branches'
    // x-slopes at the vertex.
    if (is_local_y_extremum(C.vertex(junction_vidx_in_c - 1),
                             C.vertex(junction_vidx_in_c),
                             C.vertex(junction_vidx_in_c + 1))) {
        const Point& u = C.vertex(junction_vidx_in_c - 1);
        const Point& v = C.vertex(junction_vidx_in_c);
        const Point& w = C.vertex(junction_vidx_in_c + 1);
        const bool is_max = point_y_above(v, u) && point_y_above(v, w);

        // Branch x-positions infinitesimally off v (below for a max,
        // above for a min): x = v.x + ε·(other.x − v.x)/Δy with Δy > 0.
        // prev branch left of next branch ⟺ slope_prev < slope_next,
        // cross-multiplied with positive denominators (exact).  A raw-y
        // tie on one branch makes its term the decisive sign; a tie on
        // both would mean two raw-horizontal edges at the junction,
        // excluded by the simplicity of P (they would overlap).
        const double t_prev = is_max ? (u.x - v.x) * (v.y - w.y)
                                     : (u.x - v.x) * (w.y - v.y);
        const double t_next = is_max ? (w.x - v.x) * (v.y - u.y)
                                     : (w.x - v.x) * (u.y - v.y);
        assert(t_prev != t_next &&
               "[C91 §2 tex 47]: junction branches must have distinct "
               "x-slopes (equal slopes ⟹ overlapping edges, non-simple P)");
        const bool prev_left = t_prev < t_next;

        // The inside face of C's edge junction_vidx_in_c (v→w) is the
        // face toward the prev branch.  The edge descends for a max and
        // ascends for a min; its −x face is LEFT iff ascending.
        const Side minus_x_face = is_max ? RIGHT : LEFT;
        const Side plus_x_face  = is_max ? LEFT  : RIGHT;
        const Side inside = prev_left ? minus_x_face : plus_x_face;

        Pending p;
        p.y = symbolic_y_of(v);
        p.left_edge_c = p.right_edge_c = n_c1_edges;
        p.left_side = p.right_side = inside;
        p.is_null_length = true;
        pending.push_back(p);
    }

    // ── Step 1b: Deduplicate the inventory ──────────────────────
    //
    // [C91 §3.1 tex 224]: each visible pair must appear ONCE in S —
    // Lemma 2.2's tree structure admits no duplicate chords.  The same
    // chord can legitimately enter the inventory twice:
    //   - the junction-companion chords are shared between the two
    //     passes (pass 1's a_{m+1} is the same ∂C point as pass 2's a₀,
    //     so one pass's startup record can coincide with the other
    //     pass's case (i) record);
    //   - a startup chord can coincide with a surviving old chord when
    //     a₀ / a_{m+1} "are already part of the sequence"
    //     ([C91 §3.1 tex 179]).
    // Sort + adjacent-unique on the full chord key: O(M log M), within
    // the [tex 226-227] sorting budget.
    //
    // Canonical slot order first: slots are "ascending x"; a chord of
    // zero geometric length (coincident endpoint x's — e.g. the
    // junction duplicate pair of a straight-through junction, seen once
    // from each pass with the roles of a₀ and c₀ exchanged) leaves the
    // x-order vacuous, so ties are refined deterministically by
    // clockwise ∂C position.  Without this, the same chord can enter
    // the inventory twice with swapped slots and defeat the dedup.
    {
        auto tie_pos = [&](std::size_t edge, Side side) {
            return (side == LEFT) ? edge : (2 * n_edges - 1 - edge);
        };
        for (Pending& p : pending) {
            double xl = edge_x_at_y(C, p.left_edge_c, p.y);
            double xr = edge_x_at_y(C, p.right_edge_c, p.y);
            bool swap_slots =
                xl > xr ||
                (xl == xr && tie_pos(p.left_edge_c, p.left_side) >
                                 tie_pos(p.right_edge_c, p.right_side));
            if (swap_slots) {
                std::swap(p.left_edge_c, p.right_edge_c);
                std::swap(p.left_side, p.right_side);
            }
        }
    }
    {
        auto key_less = [](const Pending& a, const Pending& b) {
            if (a.y.y != b.y.y) return a.y.y < b.y.y;
            if (a.y.tag != b.y.tag) return a.y.tag < b.y.tag;
            if (a.left_edge_c != b.left_edge_c)
                return a.left_edge_c < b.left_edge_c;
            if (a.left_side != b.left_side)
                return a.left_side < b.left_side;
            if (a.right_edge_c != b.right_edge_c)
                return a.right_edge_c < b.right_edge_c;
            if (a.right_side != b.right_side)
                return a.right_side < b.right_side;
            return a.is_null_length < b.is_null_length;
        };
        auto key_eq = [&](const Pending& a, const Pending& b) {
            return !key_less(a, b) && !key_less(b, a);
        };
        std::sort(pending.begin(), pending.end(), key_less);
        pending.erase(std::unique(pending.begin(), pending.end(), key_eq),
                      pending.end());
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
    // Identification via SoS y match per [C91 §2 tex 47].
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

    // [C91 §2.4 tex 142]: an arc crossing one of C's endpoint
    // turnarounds is emitted as ONE double-backing structure
    // (cursor.side != end_side, or a same-side span around both turns).
    auto emit_arc = [&](std::size_t end_edge, Side end_side,
                         std::size_t override_edge_count = NONE) {
        Arc a;
        a.first_edge = cursor.edge; a.first_side = cursor.side;
        a.last_edge  = end_edge;     a.last_side  = end_side;
        a.region_node = current_region;
        if (override_edge_count != NONE) {
            a.edge_count = override_edge_count;
        } else {
            auto [lo, hi] =
                a.underlying_edge_range(0, C.num_vertices() - 1);
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

    // [C91 §2 tex 47]: SoS-exact test — does a group sit exactly at a
    // polygon vertex of its edge?  Within-edge ∂C traversal: LEFT follows
    // the edge (start→end); RIGHT reverses it (end→start).
    auto group_at_vertex = [&](const Endpoint& e, std::size_t vidx) -> bool {
        return symbolic_y_equal(e.y, symbolic_y_of(C.vertex(vidx)));
    };
    auto trav_start_vertex = [&](std::size_t edge, Side side) -> std::size_t {
        return (side == LEFT) ? C.edge(edge).start_idx
                              : C.edge(edge).end_idx;
    };
    auto trav_end_vertex = [&](std::size_t edge, Side side) -> std::size_t {
        return (side == LEFT) ? C.edge(edge).end_idx
                              : C.edge(edge).start_idx;
    };

    auto process_group = [&](std::size_t i, std::size_t j) {
        const Endpoint& g = eps[i];
        // [tex 226]: emit the arc closing at this group's position.
        //
        // [C91 §2.4 tex 133]: e₁..eₜ are the edges OF the arc.  A group
        // sitting at its edge's traversal-START vertex touches edge_c at
        // a single point, so the closing arc's true extent ends at the
        // PREVIOUS traversal edge; when the cursor already sits on edge_c
        // (a preceding group at the same vertex, or the position is C's
        // start/end vertex) the closing arc has zero length
        // ([C91 §2.2 tex 96]: "some arcs may be of zero length").
        std::size_t before_arc;
        std::size_t tsv = trav_start_vertex(g.edge_c, g.side);
        if (group_at_vertex(g, tsv)) {
            if (cursor.edge == g.edge_c && cursor.side == g.side) {
                before_arc = emit_arc(g.edge_c, g.side, /*edge_count=*/0);
            } else if (tsv == 0 || tsv + 1 == C.num_vertices()) {
                // [C91 §2.1 tex 72] case 3 + [C91 §2.4 tex 142]: the
                // group sits at a C-endpoint COMPANION just past the
                // turnaround; the closing arc passes THROUGH the turn
                // and ends at the companion — one double-backing
                // structure, its last pointer on the turnaround's
                // single incident edge.
                before_arc = emit_arc(g.edge_c, g.side);
            } else {
                std::size_t prev_edge = (g.side == LEFT) ? g.edge_c - 1
                                                         : g.edge_c + 1;
                before_arc = emit_arc(prev_edge, g.side);
                cursor = {g.edge_c, g.side};
            }
        } else {
            before_arc = emit_arc(g.edge_c, g.side);
        }
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
    bool tail_is_zero = false;
    while (i < eps.size()) {
        std::size_t j = i + 1;
        while (j < eps.size() &&
               trav_pos(eps[j].edge_c, eps[j].side) ==
                   trav_pos(eps[i].edge_c, eps[i].side) &&
               symbolic_y_equal(eps[j].y, eps[i].y)) ++j;
        process_group(i, j);

        // [C91 §2.4 tex 133]: a group at its edge's traversal-END vertex —
        // the next arc begins on the FOLLOWING traversal edge.  Leaving
        // the cursor on edge_c would put a zero-coverage edge into the
        // next arc's cache range and, worse, break the input-table
        // vertex derivation of the start-y (the polygon vertex read from
        // first_edge must be the arc's true starting vertex).  No advance
        // at C's own endpoints — the next arc passes THROUGH the
        // turnaround from the companion, so its pointer stays on the
        // turnaround's single incident edge ([C91 §2.4 tex 142]).
        {
            const Endpoint& g = eps[i];
            if (group_at_vertex(g, trav_end_vertex(g.edge_c, g.side))) {
                if (g.side == LEFT && g.edge_c + 1 <= end_v_edge)
                    cursor = {g.edge_c + 1, LEFT};
                else if (g.side == RIGHT && g.edge_c > 0)
                    cursor = {g.edge_c - 1, RIGHT};
                else if (g.side == RIGHT && g.edge_c == 0)
                    // [C91 §2.1 tex 72] case 3: a chord endpoint at the
                    // RIGHT companion of C's start — only the zero
                    // passage to the turnaround remains after it.
                    tail_is_zero = true;
            }
        }
        i = j;
    }

    // Final arc: from the cursor through C's end turnaround (if not yet
    // crossed), the whole remaining RIGHT side, ending just before C's
    // start turnaround — merged below with the sweep's first arc into
    // the ONE structure double-backing around C's start
    // ([C91 §2.4 tex 142]).  For a chordless inventory this single
    // emission IS the closed arc covering all of ∂C.
    std::size_t tail_piece = emit_arc(0, RIGHT,
                                      tail_is_zero ? 0 : NONE);
    patch_after_arcs(tail_piece);

    // [C91 §3.1]: closed-loop ∂C + Jordan curve ⟹ sweep balances.
    assert(current_region == r_start &&
           "[C91 §3.1]: parenthesis sweep must return to the starting region");

    // ── Step 3b: glue the sweep's first and last pieces across C's
    // start turnaround ([C91 §2.4 tex 142]: the turnaround is never a
    // chord endpoint, so exactly one arc passes through it and it is
    // stored as ONE double-backing structure).
    std::size_t head_piece = 0;
    bool head_merged = false;
    if (arcs.size() >= 2) {
        Arc& tail = arcs[tail_piece];
        const Arc& head = arcs[head_piece];
        assert(tail.region_node == head.region_node &&
               "[C91 §2.2 tex 96]: the pieces flanking C's start "
               "turnaround bound the same region");
        assert(head.first_edge == 0 && head.first_side == LEFT &&
               tail.last_edge == 0 && tail.last_side == RIGHT &&
               "[C91 §2.4(iii) tex 138]: the sweep starts and ends at "
               "C's start turnaround");
        tail.last_edge = head.last_edge;
        tail.last_side = head.last_side;
        if (tail.edge_count == 0 && head.edge_count == 0) {
            // Zero passage around the turnaround (chord endpoints at
            // both companions): stays a zero-length wrap arc.
        } else {
            auto [lo, hi] =
                tail.underlying_edge_range(0, C.num_vertices() - 1);
            tail.edge_count = C.count_nonnull_edges(lo, hi);
        }
        head_merged = true;
    }

    // ── Step 4: Hand to Submap in canonical order ──────────────
    // [C91 §2.4(iii) tex 138]: LEFT-starting arcs ascending first_edge,
    // then RIGHT-starting descending; ties keep sweep (traversal) order,
    // so the start-turn arc lands last.  Submap::add_arc asserts the
    // order at insertion.

    std::vector<std::size_t> left_order, right_order;
    left_order.reserve(arcs.size());
    right_order.reserve(arcs.size());
    for (std::size_t ai = 0; ai < arcs.size(); ++ai) {
        if (head_merged && ai == head_piece) continue;
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

    out_S.start_vertex = 0;
    out_S.end_vertex   = C.num_vertices() - 1;

    std::vector<std::size_t> arc_remap(arcs.size(), NONE);
    auto add_in_order = [&](const std::vector<std::size_t>& order) {
        for (std::size_t old_ai : order)
            arc_remap[old_ai] = out_S.add_arc(arcs[old_ai]);
    };
    add_in_order(left_order);
    add_in_order(right_order);
    if (head_merged) arc_remap[head_piece] = arc_remap[tail_piece];

    // [C91 §2.4(iii) tex 138]: the endpoint pointers — the start-turn
    // arc is the merged tail (last table entry); the end-turn arc is the
    // last LEFT-starting arc when one exists ([C91 §2.4 tex 144]), else
    // the merged tail double-backs around both endpoints.
    out_S.start_arc = arc_remap[tail_piece];
    out_S.end_arc = left_order.empty()
        ? arc_remap[tail_piece]
        : arc_remap[left_order.back()];
    assert(out_S.start_arc == out_S.num_arcs() - 1 &&
           "[C91 §2.4(iii) tex 138]: the start-turn arc sorts last");
    assert((out_S.arc(out_S.end_arc).wraps_end() ||
            out_S.num_live_chords() == 0) &&
           "[C91 §2.4 tex 142]: end_arc double-backs around C's end");

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

    // [C91 §2.2 tex 106]: finalize the edge_count caches — nonnull ∂C
    // edges per leg (arc.h::arc_boundary_edge_count).  The sweep's
    // provisional counts used underlying C-edge ranges, which
    // single-count the doubled-over edges of wrap arcs; the exact
    // per-leg count needs the arcs' endpoint ys, derivable only now
    // that every chord and endpoint pointer is in place.  O(#arcs),
    // within the [C91 §3.1 tex 226] rebuild budget.
    out_S.refresh_arc_edge_counts(C);
}

} // namespace chazelle
