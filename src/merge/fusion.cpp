// src/merge/fusion.cpp

#include "fusion.h"
#include "../visibility/naive_visibility.h"

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

// [C91 §3.1 tex 195]: the ∂C₁ tour's y-direction as the channel
// leaves position (edge, side, y) clockwise — along the face's
// traversal, or across the corner/turnaround when the position sits
// at the face's traversal end (the next face reverses the direction
// at an extremum corner or a turnaround, keeps it otherwise; both are
// captured by re-deriving the next face's own traversal direction).
static bool leaving_clockwise_downward(const Polygon& C,
                                       std::size_t edge, Side side,
                                       SymbolicY y) {
    auto trav_asc = [&](std::size_t e, Side s_) {
        const auto& ed = C.edge(e);
        const bool asc = symbolic_y_less(
            symbolic_y_of(C.vertex(ed.start_idx)),
            symbolic_y_of(C.vertex(ed.end_idx)));
        return (s_ == LEFT) ? asc : !asc;
    };
    const auto& ed = C.edge(edge);
    const std::size_t trav_end =
        (side == LEFT) ? ed.end_idx : ed.start_idx;
    if (!symbolic_y_equal(y, symbolic_y_of(C.vertex(trav_end))))
        return !trav_asc(edge, side);        // mid-face: follow the face
    // At the traversal end: the clockwise tour continues on the next
    // face — same side, next edge (LEFT ascends the table, RIGHT
    // descends), or the same edge's other side at C's turnarounds.
    if (side == LEFT) {
        if (edge + 1 < C.num_edges())
            return !trav_asc(edge + 1, LEFT);
        return !trav_asc(edge, RIGHT);       // end turnaround
    }
    if (edge > 0)
        return !trav_asc(edge - 1, RIGHT);
    return !trav_asc(edge, LEFT);            // start turnaround
}

// ── local_shoot ─────────────────────────────────────────────────

RayHit local_shoot(Point p, Side direction,
                    std::size_t region,
                    const Submap& S, const Polygon& C,
                    const RayShootingOracle& oracle,
                    bool require_hit, double source_x_offset) {
    // [C91 §3.1 tex 181]: check each region arc, take the nearest hit.
    // Conformality bounds the loop to ≤ 4 arcs ([C91 §2.3 tex 114];
    // wrap-spanning arcs are single structures, [C91 §2.4 tex 142]).
    RegionArcs arcs = collect_region_arcs(S, region);

    RayHit best;
    best.hit = false;

    // The per-arc subarcs, retained for the tex-246 local checking of
    // the winner below.
    Subarc subs[RegionArcs::MAX];

    for (std::size_t k = 0; k < arcs.count; ++k) {
        const std::size_t ai = arcs.arcs[k];
        const auto& a = S.arc(ai);

        // Build a Subarc covering the full arc.
        Subarc& sub = subs[k];
        sub.first_edge = a.first_edge;
        sub.first_side = a.first_side;
        sub.last_edge = a.last_edge;
        sub.last_side = a.last_side;
        // [C91 §3.0(i) tex 169]: α' is "specified by its two endpoints"
        // — carry the arc's exact endpoint ys.
        sub.first_y = S.arc_start_symbolic_y(ai, C);
        sub.last_y = S.arc_end_symbolic_y(ai, C);
        assert_subarc_clockwise(sub);

        RayHit hit = oracle.shoot(p, direction, ai, sub, source_x_offset);
        if (!hit.hit) continue;

        // [C91 §2.1 tex 70]: a ray that misses everything "wraps around
        // in the spherical plane until it hits C again" — hits carry the
        // wrap flag, and a through-infinity hit lies BEHIND the source
        // (signed distance ≤ 0).  Assert the convention rather than a
        // forward-only report.
        double hit_signed_dist = (direction == LEFT)
            ? (p.x - hit.x) : (hit.x - p.x);
        // [C91 §2.1 tex 70/§2 tex 47]: d < 0 ⟹ behind; d > 0 ⟹
        // forward; raw d == 0 shares the source's raw point — direct
        // when the perturbed x-offsets put the contact strictly
        // forward (perturbed_hit_forward), else met after a full
        // wrap.
        assert((hit_signed_dist < 0.0
                    ? hit.wrapped
                    : (hit_signed_dist > 0.0 ? !hit.wrapped : true)) &&
               "[C91 §2.1 tex 70]: wrap flag consistent with the signed "
               "travel distance");

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
                // Same raw travel position — two walls at one tag-matched
                // vertex (the outside wall shields the wedge, [C91 §2.1
                // tex 72]) or a strict crossing tied with the padding
                // cluster ([C91 §4 tex 316]); resolved by the SoS order
                // of [C91 §2 tex 47] (polygon.h::ray_contact_precedes).
                // The source's own wall faces AWAY from the ray and is
                // never struck, so no self-hit is possible at distance 0.
                if (ray_contact_precedes(C, SymbolicY{p.y, p.index},
                                         direction, hit.edge, hit.side,
                                         best.edge, best.side)) {
                    best = hit;
                }
            }
        }
    }

    // [C91 §3.2 tex 246]: "local shooting reports edges of P and does
    // not tell us if the point hit is on the desired arc or is the
    // companion of a point of the arc ... local checking can decide
    // which way it is in constant time."  Attribute the winner to the
    // region arc whose exact span contains it — the reporting arc when
    // it does (a shared-endpoint contact belongs to both, keep the
    // reporter), else any containing arc, else NONE (a companion
    // contact; [C91 §3.1 tex 220]'s case (i) reads NONE as a_j ∉ R).
    if (best.hit) {
        assert(S.start_vertex != NONE && S.end_vertex != NONE &&
               "[C91 §2.4(iii)]: C endpoints must be set");
        const SymbolicY ray_y{p.y, p.index};
        std::size_t attributed = NONE;
        for (std::size_t k = 0; k < arcs.count; ++k) {
            if (!subarc_contains_point(subs[k], C, best.edge, best.side,
                                       ray_y, S.start_vertex,
                                       S.end_vertex))
                continue;
            if (arcs.arcs[k] == best.hit_arc_idx) {
                attributed = arcs.arcs[k];
                break;
            }
            if (attributed == NONE) attributed = arcs.arcs[k];
        }
        best.hit_arc_idx = attributed;
    }

    // [C91 §3.1 tex 181]: regions are closed under visibility (Lemma 2.1
    // corollary) ⟹ shooting from inside a region always hits.
    if (require_hit) {
        assert(best.hit &&
               "[C91 §3.1 tex 181]: local shoot inside a region must hit (Lemma 2.1)");
        assert(best.hit_arc_idx != NONE &&
               "[C91 §2.2 Lemma 2.1]: an in-region shot's first contact "
               "lies ON the region's boundary arcs");
    }
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
    if (xl == xr) {
        // Raw-coincident endpoints: either an extremum's outside
        // duplicate pair (both AT the vertex — neither moves under
        // the perturbation ⟹ the through-infinity chord), or two
        // strict crossings hugging one point from opposite sides under
        // a foreign raw-tied level ([C91 §2 tex 47]).  The perturbed
        // x-offsets (polygon.h::perturbed_x_offset) decide which
        // endpoint is really west; the chord is direct iff that
        // endpoint's free side faces east.
        // (At a raw x-tie the recorded slot order is arbitrary, so it
        // must not be consulted.)
        const double rl =
            perturbed_x_offset(C, c.symbolic_y(), c.left_edge);
        const double rr =
            perturbed_x_offset(C, c.symbolic_y(), c.right_edge);
        if (rl == rr) return true;      // outside duplicate pair
        const bool left_is_west = rl < rr;
        const std::size_t we = left_is_west ? c.left_edge : c.right_edge;
        const Side ws = left_is_west ? c.left_side : c.right_side;
        return shooting_direction(we, ws, C) == LEFT;
    }
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

// [C91 §2.1 tex 72] + [C91 §3.1 tex 191]: when the junction C₁ ∩ C₂ is
// a y-extremum of the MERGED curve, the walker's junction-incident edge
// has a wedge ("inside of the turn") face, and the junction companion
// sitting on it is an INSIDE duplicate: the point of ∂C it sees is its
// distance-0 sibling across the junction's null-length chord — which
// enters S structurally ([C91 §3.1 tex 224]: "this last category is
// where null-length chords originate"), never through ray-shooting.
// The wedge face is computed from the branch x-slopes at the vertex
// (same computation as rebuild_submap's junction synthesis).
struct JunctionWedge {
    bool extremum = false;
    Side inside_wall = LEFT;    // of the WALKER's junction-incident edge
};
static JunctionWedge junction_wedge(const Polygon& C1, const Polygon& C2,
                                    bool at_end) {
    JunctionWedge out;
    // Neighbors of the junction in the MERGED curve's order.
    const Point& walker_nb = at_end ? C1.vertex(C1.num_vertices() - 2)
                                    : C1.vertex(1);
    const Point& target_nb = at_end ? C2.vertex(1)
                                    : C2.vertex(C2.num_vertices() - 2);
    const Point& prev_pt = at_end ? walker_nb : target_nb;
    const Point& next_pt = at_end ? target_nb : walker_nb;
    const Point& v_pt = at_end ? C1.vertex(C1.num_vertices() - 1)
                               : C1.vertex(0);
    if (!is_local_y_extremum(prev_pt, v_pt, next_pt)) return out;
    out.extremum = true;
    const bool is_max = point_y_above(v_pt, prev_pt) &&
                        point_y_above(v_pt, next_pt);
    // Branch x-order infinitesimally off v ([C91 §2 tex 47]; null
    // padding branches included, polygon.h).
    const bool prev_left = extremum_prev_branch_left(prev_pt, v_pt,
                                                     next_pt);
    if (at_end) {
        // Walker edge = prev→v: ascends into a max, descends into a
        // min; its −x face is LEFT iff ascending.
        const Side minus_x = is_max ? LEFT : RIGHT;
        const Side plus_x  = is_max ? RIGHT : LEFT;
        out.inside_wall = prev_left ? plus_x : minus_x;
    } else {
        // Walker edge = v→next: descends from a max, ascends from a
        // min.
        const Side minus_x = is_max ? RIGHT : LEFT;
        const Side plus_x  = is_max ? LEFT : RIGHT;
        out.inside_wall = prev_left ? minus_x : plus_x;
    }
    return out;
}

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
    const JunctionWedge jwedge = junction_wedge(C1, C2, at_end);
    const bool junction_inside_pair =
        jwedge.extremum && a0.side == jwedge.inside_wall;

    RayHit hit_c1{};
    if (!junction_inside_pair)
        hit_c1 = local_shoot(a0_point, a0_dir, a0_region_s1,
                             S1, C1, oracle1, /*require_hit=*/true,
                             perturbed_x_offset(C1, a0.y, a0.edge));

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
        // Cross-frame shot: a₀'s perturbed x-offset is computed in
        // ITS frame (C₁) — identical in every frame ([C91 §2 tex 47]).
        hit_c2 = local_shoot(a0_point, a0_dir, a0_region_s2,
                             S2, C2, oracle2, /*require_hit=*/true,
                             perturbed_x_offset(C1, a0.y, a0.edge));
#ifdef CHAZELLE_TRACE_FUSION
    if (!junction_inside_pair)
        std::fprintf(stderr,
            "[startup] a0=(%zu,%c,x=%.6g)@%g@%zu dir=%c region2=%zu "
            "c1=(%zu,%c,x=%.6g,w=%d) c2=(%zu,%c,x=%.6g,w=%d)\n",
            a0.edge, a0.side == LEFT ? 'L' : 'R', a0_point.x, a0.y.y,
            a0.y.tag, a0_dir == LEFT ? 'L' : 'R', a0_region_s2,
            hit_c1.edge, hit_c1.side == LEFT ? 'L' : 'R', hit_c1.x,
            (int)hit_c1.wrapped, hit_c2.edge,
            hit_c2.side == LEFT ? 'L' : 'R', hit_c2.x,
            (int)hit_c2.wrapped);
#endif

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
        bool c2_at_or_before;
        if (hit_c2.wrapped != hit_c1.wrapped) {
            c2_at_or_before = !hit_c2.wrapped;
        } else if (d2 != d1) {
            c2_at_or_before = d2 < d1;
        } else {
            // Same raw travel position (only at the junction point or
            // its padding cluster, [C91 §3 tex 160] + [C91 §4 tex
            // 316]): resolve by the wall order of [C91 §2 tex 47] in
            // the merged frame; wall-equivalent ties keep the tex-191
            // default ("we still go on saying that c₀ sees a point of
            // ∂C₂").
            const Polygon Cm = at_end ? Polygon(C1, C2)
                                      : Polygon(C2, C1);
            const std::size_t off_w = at_end ? 0 : C2.num_edges();
            const std::size_t off_t = at_end ? C1.num_edges() : 0;
            c2_at_or_before = !ray_contact_precedes(
                Cm, a0.y, a0_dir, off_w + hit_c1.edge, hit_c1.side,
                off_t + hit_c2.edge, hit_c2.side);
        }
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
            // lies on — the one at the target companion the sightline
            // passes through.  In the target's frame the junction sits
            // at the OPPOSITE end ([C91 §3 tex 160]) — edge 0 when the
            // walker's junction is at its end, the target's last edge
            // otherwise.  Side labels are per-edge-orientation notions
            // and FLIP across an extremum junction (the two junction
            // edges ascend/descend oppositely there), so geometric
            // face coincidence is a shooting-direction test, not a
            // side-label match; at a non-extremum junction the two
            // tests agree.  For a real sightline the coincident
            // companion's ray is COLLINEAR with a₀'s (same direction);
            // in the [C91 §3.1 tex 191] inside-pair case a₀c₀ is the
            // degenerate segment to the FACING sibling, whose ray
            // points the opposite way.
            const std::size_t tj_edge = at_end ? 0 : C2.num_edges() - 1;
            const Side coin_dir =
                junction_inside_pair ? (a0_dir == LEFT ? RIGHT : LEFT)
                                     : a0_dir;
            const bool incident_on_a0 =
                (ch.left_edge == tj_edge &&
                 shooting_direction(tj_edge, ch.left_side, C2) ==
                     coin_dir) ||
                (ch.right_edge == tj_edge &&
                 shooting_direction(tj_edge, ch.right_side, C2) ==
                     coin_dir);
            if (!incident_on_a0) continue;

            // [C91 §3.1 tex 191]: a₀c₀ lies on this S₂ chord.  Elect
            // "the region we locally enter as we leave p clockwise on
            // ∂C₁."  `leaving_downward` encodes the y-direction of
            // motion from p (set differently for case (i) leaving a₀
            // vs case (ii) leaving c₀).  The chord's above/below
            // regions are constant along its whole (spherical) length —
            // wrapped chords and the degenerate polar-cap configuration
            // included — per Submap::chord_regions_below_above.
            assert(!ch.is_null_length &&
                   "[C91 §2.1 tex 72]: chords at the junction's symbolic "
                   "y with an endpoint at the target's junction "
                   "companion are the companion-pair chords (non-null)");
            std::size_t below_r = NONE, above_r = NONE;
            S2.chord_regions_below_above(ci, C2, &below_r, &above_r);
            assert((below_r == ch.region[0] || below_r == ch.region[1]) &&
                   (above_r == ch.region[0] || above_r == ch.region[1]));
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
    state.invalidated_walker_chords.assign(S1.num_chords(), false);
    state.invalidated_target_chords.assign(S2.num_chords(), false);

    // [C91 §3.1 tex 199]: k = index of arc A_k of S₁ containing p, tie-break
    // p ≠ a_k.  A_j runs from a_{j-1} to a_j cw on ∂C₁ (A_1 = a_0→a_1,
    // A_{m+1} = a_m→a_{m+1}).
    std::size_t k = fusion_startup(state, S1, C1, S2, C2, oracle1, oracle2);

    // [C91 §3.1 tex 179]: tour orientation (see FusionState).
    const bool at_end = state.junction_at_end;
    const std::size_t junction_v = at_end ? C1.num_vertices() - 1 : 0;

    // [C91 §2.1 tex 72] + [C91 §3.1 tex 191/224]: the junction wedge,
    // fixed for the whole pass.  A stop AT the junction's inside-
    // duplicate position — the added companion (tex 179) or a real
    // endpoint of the walked submap's own junction-companion chord —
    // sees its facing sibling at distance 0: case (i) holds with the
    // junction's null-length chord, which enters S structurally
    // ("this last category is where null-length chords originate"),
    // and nothing new is recorded.
    const JunctionWedge jwedge =
        junction_wedge(C1, C2, state.junction_at_end);
    const std::size_t jedge =
        state.junction_at_end ? C1.num_edges() - 1 : 0;
    const SymbolicY jy = symbolic_y_of(
        C1.vertex(state.junction_at_end ? C1.num_vertices() - 1 : 0));
    auto stop_is_junction_inside = [&](const FusionVertex& v) -> bool {
        return jwedge.extremum && v.side == jwedge.inside_wall &&
               v.edge == jedge && symbolic_y_equal(v.y, jy);
    };
    // The TARGET's junction-incident edge and wedge ([C91 §3 tex 160]:
    // the junction sits at the target's opposite end; swapping the
    // curves flips the junction end).
    const std::size_t jedge_t =
        state.junction_at_end ? 0 : C2.num_edges() - 1;
    const JunctionWedge jwedge_t =
        junction_wedge(C2, C1, !state.junction_at_end);
    // The merged curve C = C₁ ∪ C₂ in P's order, as an O(1) view
    // ([C91 §2.4 tex 133]) — the frame for wall-order tie-breaks at the
    // junction point ([C91 §2 tex 47]).
    const Polygon Cm = state.junction_at_end ? Polygon(C1, C2)
                                             : Polygon(C2, C1);
    const std::size_t off_walker = state.junction_at_end
                                       ? 0 : C2.num_edges();
    const std::size_t off_target = state.junction_at_end
                                       ? C1.num_edges() : 0;

#ifdef CHAZELLE_TRACE_FUSION
    std::fprintf(stderr, "[fuse] walk C1 off=%zu len=%zu -> C2 off=%zu len=%zu at_end=%d k0=%zu seq=%zu\n",
                 C1.table_offset(), C1.num_vertices(), C2.table_offset(),
                 C2.num_vertices(), (int)state.junction_at_end, k,
                 state.sequence.size());
#endif
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
#ifndef NDEBUG
        {
            // [C91 §3.1 tex 195] invariant (B): "The point q of ∂C that
            // is seen by p belongs to ∂C₂ and the chord pq lies in the
            // region of S₂ called current."  Normally q lies on an arc
            // of the current region ([C91 §3.1 tex 181]); but when p
            // lies ON an exit chord of the region (the tex-191/tex-206
            // election places it there), its sight along the chord's
            // line may land exactly on that chord's ENDPOINT — a
            // boundary point of the region that no arc of the region
            // contains — with pq running along the exit chord inside
            // the region's closure.
            RayHit b = local_shoot(
                state.p,
                shooting_direction(state.p_edge, state.p_side, C1),
                state.s2_region, S2, C2, oracle2, /*require_hit=*/false,
                perturbed_x_offset(C1, state.p_y, state.p_edge));
            bool sees_region_arc = b.hit && b.hit_arc_idx != NONE;
            bool on_chord_level = false;
            for (std::size_t ci :
                 S2.node(state.s2_region).incident_chords) {
                const Chord& ch = S2.chord(ci);
                if (!ch.dead &&
                    symbolic_y_equal(ch.symbolic_y(), state.p_y))
                    on_chord_level = true;
            }
            // Third disjunct: p at the junction's inside-duplicate
            // position sees its distance-0 sibling ([C91 §2.1 tex 72]).
            const bool p_at_junction_inside =
                jwedge.extremum && state.p_side == jwedge.inside_wall &&
                state.p_edge == jedge && symbolic_y_equal(state.p_y, jy);
            assert((sees_region_arc || on_chord_level ||
                    p_at_junction_inside) &&
                   "[C91 §3.1 invariant (B)]: p must see ∂C₂ (through "
                   "the region's arcs, along its own exit chord, or as "
                   "the junction inside duplicate)");
        }
#endif

        // [C91 §3.1 tex 199]: R = current region; snapshot fixed for the
        // inner walk; updates take effect next outer iteration.
        const std::size_t R = state.s2_region;

        // Position of a fusion vertex: edge-interp at chord y, or junction
        // for companions.
        auto fv_point = [&](const FusionVertex& v) -> Point {
            if (v.is_companion)
                return C1.vertex(junction_v);
            // [C91 §2 tex 47]: the point rides its chord's symbolic
            // level — Point.index carries the SoS tag for the oracles
            // (SOS_NONE would sort BELOW every real tag at a raw-y tie,
            // inverting the perturbed order).
            return Point{edge_x_at_y(C1, v.edge, v.y), v.y.y, v.y.tag};
        };

        // [C91 §3.1 tex 220]: case (i) test in O(f(γ₂)).
        struct CaseIResult { bool fires; RayHit s_hit; };
        auto case_i_test = [&](std::size_t j) -> CaseIResult {
            const FusionVertex& aj_v = state.sequence[j];

            // [C91 §2.1 tex 72]: a null-length chord's endpoints are an
            // inside duplicate pair — each sees exactly its distance-0
            // sibling on ∂C₁ itself, so case (i) ("the point of ∂C that
            // a_j sees belongs to ∂C₂") never holds at such a stop, and
            // no ray is shot from it (ray_contact_precedes
            // precondition).
            if (aj_v.chord_idx != NONE &&
                S1.chord(aj_v.chord_idx).is_null_length)
                return {false, {}};

            Point aj_point = fv_point(aj_v);
            Side dir = shooting_direction(aj_v.edge, aj_v.side, C1);

            // shoot from a_j; no hit ⟹ a_j ∉ R.
            RayHit s_hit = local_shoot(aj_point, dir, R, S2, C2, oracle2,
                                       /*require_hit=*/false,
                                       perturbed_x_offset(
                                           C1, aj_v.y, aj_v.edge));
#ifdef CHAZELLE_TRACE_FUSION
            std::fprintf(stderr,
                         "[case-i] j=%zu from=(%zu,%c,x=%.6g) dir=%c "
                         "s_hit=%d (%zu,%c,x=%.6g w=%d arc=%zd)\n",
                         j, aj_v.edge, aj_v.side == LEFT ? 'L' : 'R',
                         aj_point.x, dir == LEFT ? 'L' : 'R',
                         (int)s_hit.hit, s_hit.edge,
                         s_hit.side == LEFT ? 'L' : 'R', s_hit.x,
                         (int)s_hit.wrapped, (ptrdiff_t)s_hit.hit_arc_idx);
#endif
            if (!s_hit.hit) return {false, {}};

            // [C91 §3.1 tex 220]: "Whether a_j lies in R can be directly
            // inferred from the local orientation of the hit at s and
            // which side of the double boundary is hit."  The tex-246
            // local checking in local_shoot attributed s to the region
            // arc whose exact span contains it — NONE means s is a
            // companion contact (the ray was absorbed by the far side
            // of an arc's underlying curve before reaching the region),
            // so a_j does not lie in R.
            if (s_hit.hit_arc_idx == NONE) return {false, {}};
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
            std::size_t t_edge;
            Side t_side;
            if (aj_v.chord_idx != NONE) {
                const Chord& ch = S1.chord(aj_v.chord_idx);
                std::size_t other_edge = aj_v.is_left_endpoint
                    ? ch.right_edge : ch.left_edge;
                t_edge = other_edge;
                t_side = aj_v.is_left_endpoint ? ch.right_side
                                               : ch.left_side;
                t_x = edge_x_at_y(C1, other_edge, ch.symbolic_y());
                // [C91 §2.1 tex 70]: the chord itself may run through
                // infinity — then its other endpoint lies at or behind
                // a_j in the shooting direction.
                double t_signed = (dir == LEFT) ? (aj_point.x - t_x)
                                                : (t_x - aj_point.x);
                if (t_signed != 0.0) {
                    t_wrapped = (t_signed < 0.0);
                } else {
                    // [C91 §2 tex 47]: t at raw distance 0 shares
                    // a_j's raw point (both endpoints of a_j's own
                    // chord hugging one raw-tied vertex, e.g. a direct
                    // wedge cross-chord) — the perturbed x-offsets
                    // decide strictly-forward (direct) vs met after a
                    // full wrap.
                    const double aj_offset = perturbed_x_offset(
                        C1, ch.symbolic_y(), aj_v.edge);
                    const double t_offset = perturbed_x_offset(
                        C1, ch.symbolic_y(), t_edge);
                    t_wrapped =
                        (aj_offset == t_offset) ||
                        ((dir == RIGHT) ? !(t_offset > aj_offset)
                                        : !(t_offset < aj_offset));
                }
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
                                           S1, C1, oracle1,
                                           /*require_hit=*/true,
                                           perturbed_x_offset(
                                               C1, aj_v.y, aj_v.edge));
                t_x = t_hit.x;
                t_wrapped = t_hit.wrapped;
                t_edge = t_hit.edge;
                t_side = t_hit.side;
            }

            // (i) fires iff s strictly precedes t in ray order —
            // lexicographic (wrapped, distance) per [C91 §2.1 tex 70];
            // a raw-position tie (the junction point / its padding
            // cluster) resolves by the [C91 §2 tex 47] wall order in
            // the merged frame, with the no-fire default on
            // wall-equivalence ("strictly precedes").
            double s_dist = (dir == LEFT)
                ? (aj_point.x - s_hit.x) : (s_hit.x - aj_point.x);
            double t_dist = (dir == LEFT)
                ? (aj_point.x - t_x) : (t_x - aj_point.x);
            bool s_first;
            if (s_hit.wrapped != t_wrapped)
                s_first = !s_hit.wrapped;
            else if (s_dist != t_dist)
                s_first = s_dist < t_dist;
            else
                s_first = ray_contact_precedes(
                    Cm, aj_v.y, dir, off_target + s_hit.edge, s_hit.side,
                    off_walker + t_edge, t_side);
            return {s_first, s_hit};
        };

        // Position of an S₂ chord endpoint: edge-interp at chord y.
        // [C91 §2 tex 47]: rides the chord's symbolic level (tag on
        // Point.index for the oracles' SoS comparisons).
        auto s2_endpoint_point = [&](std::size_t edge, SymbolicY y) -> Point {
            return Point{edge_x_at_y(C2, edge, y), y.y, y.tag};
        };

        // CW position on ∂C₁ as (trav_pos, within_edge, symbolic y);
        // trav_pos mirrors build_fusion_sequence's tour for the pass's
        // orientation.  LEFT side: within = edge param t; RIGHT side:
        // 1−t (cw traversal is end→start).
        struct CwPos {
            std::size_t tp;
            double t;
            SymbolicY y;
            std::size_t edge;
            Side side;
        };
        auto cw_position = [&](SymbolicY y, std::size_t edge, Side side)
                           -> CwPos {
            std::size_t n_edges = C1.num_edges();
            std::size_t tp;
            if (at_end)
                tp = (side == RIGHT) ? (n_edges - 1) - edge
                                     : n_edges + edge;
            else
                tp = (side == LEFT) ? edge
                                    : 2 * n_edges - 1 - edge;
            double t = edge_t_at_y(C1, edge, y);
            return {tp, (side == LEFT) ? t : (1.0 - t), y, edge, side};
        };
        // Strict cw order: (trav_pos, within-edge t), with raw-t ties
        // between SoS-distinct ys broken symbolically along the edge's
        // traversal direction — same rule as rebuild's vertex_before.
        // [C91 §2 tex 47]: "strictly follows p" (tex 200) and "the last
        // one encountered" (tex 206) are symbolic notions; collapsing
        // raw-y ties would drop candidates between same-raw-y chords.
        auto cw_less = [&](const CwPos& u, const CwPos& v) -> bool {
            if (u.tp != v.tp) return u.tp < v.tp;
            if (u.t != v.t) return u.t < v.t;
            if (symbolic_y_equal(u.y, v.y)) return false;
            assert(u.edge == v.edge && u.side == v.side &&
                   "trav_pos is injective on (edge, side)");
            const auto& e = C1.edge(u.edge);
            bool edge_ascending = symbolic_y_less(
                symbolic_y_of(C1.vertex(e.start_idx)),
                symbolic_y_of(C1.vertex(e.end_idx)));
            bool trav_asc = (u.side == LEFT) ? edge_ascending
                                             : !edge_ascending;
            return trav_asc ? symbolic_y_less(u.y, v.y)
                            : symbolic_y_greater(u.y, v.y);
        };
        // [C91 §3.1 tex 200/206]: "strictly follows p" and "the last
        // one encountered as we traverse ∂C₁ clockwise starting from p"
        // are CLOCKWISE-FROM-p notions — cyclic, anchored at the pass's
        // origin a₀.  When the junction inside wall lies on the
        // linearization's closing face (an at_end=0 pass whose a₀ is a
        // RIGHT-side wall, or symmetrically), A₁ legitimately crosses
        // the linear seam at the junction turnaround, so a plain linear
        // comparison would reject every candidate beyond the seam.
        auto cw_less_from_a0 = [&](const CwPos& u, const CwPos& v) -> bool {
            const CwPos a0_cw = cw_position(state.sequence[0].y,
                                            state.sequence[0].edge,
                                            state.sequence[0].side);
            const bool u_ge = !cw_less(u, a0_cw);
            const bool v_ge = !cw_less(v, a0_cw);
            if (u_ge != v_ge) return u_ge;  // [a₀, seam) precedes the rest
            return cw_less(u, v);
        };

        // [C91 §3.1 tex 222]: case (ii) test in O(f(γ₁)).
        struct CaseIIResult {
            bool fires;
            RayHit p_prime_hit;
            std::size_t chord_idx;
            bool q_is_left;
            // [C91 §3.1 tex 224]: true when the discovered pair is the
            // junction inside pair — its chord is the synthesized
            // null-length chord, so the action records nothing new.
            bool suppress_record = false;
        };
        // [C91 §3.1 tex 199]: A_j is the segment of ∂C₁ from a_{j-1} to
        // a_j cw.  A wrap-spanning A_j is ONE arc-structure that
        // double-backs ([C91 §2.4 tex 142]) — wraps never split it — but
        // a null-length chord strictly inside A_j (its endpoints are not
        // enumeration stops, [C91 §3.1 tex 224]) still separates two
        // structures — L interior null chords leave L + 1 structures.
        // Enumerate them ALL, each with its clipped subarc: the end
        // structures clipped at a_{j-1} / a_j, interior ones whole.
        // [C91 §3.1 tex 222] shoots "toward A_j" — skipping an interior
        // structure would silently drop candidates whose first contact
        // lies on it.  Σ(1+L) over all A_j is within tex 224's O(m)
        // shot budget (null chords count toward the submap's region
        // total, Lemma 2.3).
        struct AjSpan { std::size_t arc; Subarc sub; };
        auto Aj_spans = [&](std::size_t j) -> std::vector<AjSpan> {
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
            // [C91 §3.1 tex 179]: "it could happen that a₀ and
            // a_{m+1} are already part of the sequence" — a real chord
            // endpoint coincides with an added junction companion, so
            // consecutive stops can share one ∂C position and A_j is a
            // single point.  Its only point can never STRICTLY follow
            // p (tex 200), so case (ii) cannot fire on it; emit the
            // degenerate span for the oracle contract.  (Walking the
            // table from arc-after to arc-before would instead wrap
            // the whole boundary.)
            {
                const FusionVertex& va = state.sequence[j - 1];
                const FusionVertex& vb = state.sequence[j];
                if (va.edge == vb.edge && va.side == vb.side &&
                    symbolic_y_equal(va.y, vb.y)) {
                    std::vector<AjSpan> spans;
                    spans.push_back({arc_after(va),
                                     Subarc{va.edge, va.side, vb.edge,
                                            vb.side, va.y, vb.y}});
                    return spans;
                }
            }

            const std::size_t first = arc_after(state.sequence[j - 1]);
            const std::size_t last  = arc_before(state.sequence[j]);
            assert(first != NONE && last != NONE &&
                   "[C91 §3.1 tex 199]: A_j is delimited by real stops "
                   "(a₀ opens and a_{m+1} closes the tour, so neither "
                   "companion NONE case is reachable here)");

            // [C91 §3.0(i) tex 169]: each span's Subarc is "specified
            // by its two endpoints" — exact ys from the delimiting
            // stops (fusion-sequence positions) or the structure's own
            // endpoints (its bounding null-chord positions).
            std::vector<AjSpan> spans;
            if (first == last) {
                // Single-structure A_j: full span from a_{j-1} to a_j.
                spans.push_back({first,
                                 Subarc{state.sequence[j - 1].edge,
                                        state.sequence[j - 1].side,
                                        state.sequence[j].edge,
                                        state.sequence[j].side,
                                        state.sequence[j - 1].y,
                                        state.sequence[j].y}});
                return spans;
            }
            [[maybe_unused]] std::size_t guard = 0;
            for (std::size_t ai = first;; ai = cw_successor(ai)) {
                assert(++guard <= S1.num_arcs() &&
                       "[C91 §2.4(iii) tex 138]: the table walk from "
                       "arc-after(a_{j-1}) must reach arc-before(a_j)");
                const Arc& a = S1.arc(ai);
                assert(!a.dead &&
                       "[C91 §2.4]: normal-form S₁ has no dead arcs");
                if (ai == first) {
                    // First structure: a_{j-1} → its end (the first
                    // interior null-chord position).
                    spans.push_back({ai,
                                     Subarc{state.sequence[j - 1].edge,
                                            state.sequence[j - 1].side,
                                            a.last_edge, a.last_side,
                                            state.sequence[j - 1].y,
                                            S1.arc_end_symbolic_y(ai, C1)}});
                } else if (ai == last) {
                    // Last structure: its start (just past the final
                    // interior null-chord position) → a_j.
                    spans.push_back({ai,
                                     Subarc{a.first_edge, a.first_side,
                                            state.sequence[j].edge,
                                            state.sequence[j].side,
                                            S1.arc_start_symbolic_y(ai, C1),
                                            state.sequence[j].y}});
                } else if (a.edge_count != 0) {
                    // Interior structure between two null chords: whole.
                    spans.push_back({ai,
                                     Subarc{a.first_edge, a.first_side,
                                            a.last_edge, a.last_side,
                                            S1.arc_start_symbolic_y(ai, C1),
                                            S1.arc_end_symbolic_y(ai, C1)}});
                }
                // Zero-length interior structures are skipped: their
                // only points are inside-pair duplicates whose
                // visibility is settled at distance 0 ([C91 §2.1
                // tex 72]).
                if (ai == last) break;
            }
            return spans;
        };

        auto case_ii_test = [&](std::size_t j) -> CaseIIResult {
            CaseIIResult best{false, {}, NONE, false};
            CwPos best_cw{};
            bool best_cw_valid = false;

            auto aj_spans = Aj_spans(j);
            assert(!aj_spans.empty() &&
                   "[C91 §3.1 tex 199]: A_j spans at least one structure");

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
                    // [C91 §2.1 tex 72]: is q the TARGET's junction
                    // INSIDE duplicate?  Its visibility is exactly its
                    // facing sibling — the WALKER's junction inside
                    // duplicate — at distance 0 (the junction's
                    // null-length chord, which enters S structurally
                    // per [C91 §3.1 tex 224]).  Case (ii) may still
                    // fire on it (the paper's region re-entry mechanics
                    // rely on junction-chord endpoints), but only with
                    // p' = that sibling, and nothing new is recorded.
                    const bool q_is_junction_inside =
                        jwedge_t.extremum && q_edge == jedge_t &&
                        q_side == jwedge_t.inside_wall &&
                        symbolic_y_equal(chord_ab_y, jy);
                    Point q_point = is_left ? a_pt : b_pt;
                    Point other_point = is_left ? b_pt : a_pt;
                    Side shoot_dir = shooting_direction(q_edge, q_side, C2);

                    // [C91 §3.1 tex 222]: shoot toward A_j.  A wrapping
                    // A_j is ONE arc-structure and ONE oracle call with
                    // a double-backing subarc ([C91 §2.4 tex 142] /
                    // [C91 §3.0(i) tex 169]: subarcs carry side flags).
                    // Null-length chords strictly inside A_j split it
                    // into further structures, each shot with its own
                    // clipped subarc (see Aj_spans).
                    for (const AjSpan& span : aj_spans) {
                        const std::size_t aj_arc = span.arc;
                        const Arc& aj_arc_struct = S1.arc(aj_arc);
                        const Subarc& aj_sub = span.sub;

                        assert_subarc_clockwise(aj_sub);
                        // Cross-frame shot: q's perturbed x-offset in
                        // ITS frame (C₂) — the junction point is the
                        // one cross-curve raw coincidence
                        // ([C91 §3 tex 160]).
                        RayHit hit = oracle1.shoot(
                            q_point, shoot_dir, aj_arc, aj_sub,
                            perturbed_x_offset(C2, chord_ab_y, q_edge));
#ifdef CHAZELLE_TRACE_FUSION
                        std::fprintf(stderr,
                            "[case-ii] ci=%zu q=(%zu,%c,x=%.6g)@%g@%zu "
                            "dir=%c hit=%d (%zu,%c,x=%.6g,w=%d)\n",
                            ci, q_edge, q_side == LEFT ? 'L' : 'R',
                            q_point.x, chord_ab_y.y, chord_ab_y.tag,
                            shoot_dir == LEFT ? 'L' : 'R',
                            (int)hit.hit, hit.edge,
                            hit.side == LEFT ? 'L' : 'R', hit.x,
                            (int)hit.wrapped);
#endif
                        if (!hit.hit) continue;
#ifdef CHAZELLE_TRACE_FUSION
                        auto cii = [&](const char* stage) {
                            std::fprintf(stderr,
                                         "[case-ii]   reject at %s\n",
                                         stage);
                        };
#else
                        auto cii = [&](const char*) {};
#endif

                        // [C91 §2.1 tex 70]: hits are ordered in the
                        // lexicographic (wrapped, distance) ray metric;
                        // wrapped ⟺ at-or-behind the source.
                        double q_to_hit = (shoot_dir == LEFT)
                            ? (q_point.x - hit.x) : (hit.x - q_point.x);
                        // [C91 §2.1 tex 70/§2 tex 47]: d < 0 ⟹ behind;
                        // d > 0 ⟹ forward; raw d == 0 goes by the
                        // perturbed x-offsets (perturbed_hit_forward).
                        assert((q_to_hit < 0.0
                                    ? hit.wrapped
                                    : (q_to_hit > 0.0 ? !hit.wrapped
                                                      : true)) &&
                               "[C91 §2.1 tex 70]: wrap flag consistent "
                               "with the signed travel distance");

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
                        // A raw tie with the far endpoint (its padding
                        // cluster / the junction point) resolves by the
                        // [C91 §2 tex 47] wall order: on-ab ⟺ the hit
                        // is not strictly beyond the far endpoint.
                        const std::size_t oth_e = is_left
                            ? chord_ab.right_edge : chord_ab.left_edge;
                        const Side oth_s = is_left ? chord_ab.right_side
                                                   : chord_ab.left_side;
                        auto beyond_other = [&]() -> bool {
                            return ray_contact_precedes(
                                Cm, chord_ab_y, shoot_dir,
                                off_target + oth_e, oth_s,
                                off_walker + hit.edge, hit.side);
                        };
                        if (ab_wraps) {
                            // ab through infinity: its free span covers
                            // every ray position up to the far endpoint
                            // in the wrap metric — on-ab ⟺
                            // lex(hit) ≤ (wrapped, q_to_other).
                            if (hit.wrapped &&
                                (q_to_hit > q_to_other ||
                                 (q_to_hit == q_to_other &&
                                  beyond_other()))) {
                                cii("on-ab-wrap");
                                continue;
                            }
                        } else {
                            // ab direct: on-ab ⟺ direct and not beyond.
                            if (hit.wrapped || q_to_hit > q_to_other ||
                                (q_to_hit == q_to_other &&
                                 beyond_other())) {
                                cii("on-ab-direct");
                                continue;
                            }
                        }

                        // "proper orientation": hit.Side matches A_j's
                        // Side at hit.edge ([C91 §3.1 tex 222]).  A
                        // wrapping A_j covers per-leg side/edge zones
                        // ([C91 §2.4 tex 142]).  The [C91 §3.0(i)
                        // tex 169] oracle reports "the single point of
                        // α'" — an endpoint-exact contract, so the
                        // report lies on the shot subarc by definition.
                        // The hit rides ab's symbolic level chord_ab_y
                        // (RayHit carries only raw y).
                        assert(S1.start_vertex != NONE &&
                               S1.end_vertex != NONE &&
                               "[C91 §2.4(iii)]: S₁'s C endpoints must "
                               "be set");
                        // [C91 §3.1 tex 222]: "if it does not have the
                        // proper orientation which lets a see A_j
                        // without the other side of the double boundary
                        // of A_j interfering (a constant-time test),
                        // then the endpoint can be disqualified" — the
                        // [C91 §3.0(i) tex 169] report is the first
                        // ᾱ'-contact, companion side included
                        // ([C91 §3.2 tex 246]).
                        if (!subarc_contains_point(aj_sub, C1,
                                                   hit.edge, hit.side,
                                                   chord_ab_y,
                                                   S1.start_vertex,
                                                   S1.end_vertex)) {
                            cii("orientation");
                            continue;
                        }

                        // "occurs before p along A_j": strict cw compare.
                        // hit.y inherits the ray's perturbed source y =
                        // chord_ab's SymbolicY (RayHit carries only raw y).
                        auto hit_cw = cw_position(chord_ab_y, hit.edge, hit.side);
                        if (!cw_less_from_a0(p_cw, hit_cw)) {
                            cii("follows-p");
                            continue;
                        }

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
                        // [C91 §2 tex 47]: s sits on ab's sightline, so
                        // the back-shot rides ab's symbolic level.
                        Point s_point{hit.x, hit.y, chord_ab_y.tag};
                        // s is a point of ∂C₁ on hit.edge — the same
                        // frame as oracle1, so its perturbed x-offset
                        // classifies raw-distance-0 contacts around it
                        // ([C91 §2 tex 47]).
                        const double s_offset =
                            perturbed_x_offset(C1, chord_ab_y,
                                                hit.edge);
                        RayHit t_hit = local_shoot(
                            s_point, s_back_dir,
                            aj_arc_struct.region_node,
                            S1, C1, oracle1, /*require_hit=*/true,
                            s_offset);
                        double s_to_q = (s_back_dir == LEFT)
                            ? (s_point.x - q_point.x)
                            : (q_point.x - s_point.x);
                        double s_to_t = (s_back_dir == LEFT)
                            ? (s_point.x - t_hit.x)
                            : (t_hit.x - s_point.x);
                        bool suppress = false;
                        if (q_is_junction_inside) {
                            // [C91 §2.1 tex 72]: q sees ONLY its facing
                            // sibling — the walker's junction inside
                            // duplicate — at distance 0.  Case (ii)
                            // fires iff p' IS that sibling (which the
                            // strictly-follows-p test above admitted);
                            // the pair's chord is the junction's
                            // null-length chord, synthesized by the
                            // rebuild, so nothing is recorded.
                            if (!(hit.edge == jedge &&
                                  hit.side == jwedge.inside_wall &&
                                  jwedge.extremum))
                                continue;
                            suppress = true;
                        } else {
                            // [C91 §2.1 tex 70/§2 tex 47]: q at raw
                            // distance 0 shares s's raw point (the
                            // junction, [C91 §3 tex 160]) — the
                            // perturbed x-offsets decide whether q is
                            // strictly FORWARD of s (offsets are
                            // frame-independent; q's is computed in
                            // its own frame C₂).  Equal distances mean q
                            // and t sit at one raw point — resolve by
                            // the wall order of [C91 §2 tex 47] in the
                            // merged frame.
                            bool q_behind;
                            if (s_to_q != 0.0) {
                                q_behind = (s_to_q < 0.0);
                            } else {
                                const double q_offset = perturbed_x_offset(
                                    C2, chord_ab_y, q_edge);
                                q_behind =
                                    (s_offset == q_offset) ||
                                    ((s_back_dir == RIGHT)
                                         ? !(q_offset > s_offset)
                                         : !(q_offset < s_offset));
                            }
                            bool q_first;
                            if (q_behind != t_hit.wrapped)
                                q_first = !q_behind;
                            else if (s_to_q != s_to_t)
                                q_first = s_to_q < s_to_t;
                            else
                                q_first = !ray_contact_precedes(
                                    Cm, chord_ab_y, s_back_dir,
                                    off_walker + t_hit.edge, t_hit.side,
                                    off_target + q_edge, q_side);
                            if (!q_first) { cii("back-shot"); continue; }
                        }

                        // [C91 §3.1 tex 206]: pick LAST p' cw from p.
                        if (!best_cw_valid ||
                            cw_less_from_a0(best_cw, hit_cw)) {
                            best = {true, hit, ci, is_left, suppress};
                            best_cw = hit_cw;
                            best_cw_valid = true;
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
            if (j == state.sequence.size()) {
#ifdef CHAZELLE_TRACE_FUSION
                std::fprintf(stderr, "[fuse]  j=%zu (iii) terminate\n", j);
#endif
                return;
            }

            state.current_stop = j;
#ifdef CHAZELLE_TRACE_FUSION
            if (j < state.sequence.size())
                std::fprintf(stderr, "[fuse]   probe j=%zu stop=(%zu,%s,y=%g@%zu comp=%d ch=%zd) R=%zu\n",
                             j, state.sequence[j].edge,
                             state.sequence[j].side==LEFT?"L":"R",
                             state.sequence[j].y.y, state.sequence[j].y.tag,
                             (int)state.sequence[j].is_companion,
                             (ptrdiff_t)state.sequence[j].chord_idx, R);
#endif

            // [C91 §3.1 tex 191/199/224]: a stop at the junction's
            // INSIDE-duplicate position — the tour-end companion, or a
            // real endpoint of the walked submap's junction-companion
            // chord coinciding with it (tex 179: "it could happen that
            // a₀ and a_{m+1} are already part of the sequence").  Case
            // (i) holds with the seen point = the distance-0 sibling
            // (∂C₂): nothing new is recorded (the null chord is
            // synthesized structurally), the stop's own S₁ chord is
            // superseded ([C91 §3.1 tex 224] invalidation), p moves to
            // the stop, and R is unchanged; at a_{m+1} the pass simply
            // terminates (tex 199).
            if (stop_is_junction_inside(state.sequence[j])) {
                const FusionVertex& aj_v = state.sequence[j];
#ifdef CHAZELLE_TRACE_FUSION
                std::fprintf(stderr, "[fuse]  j=%zu inside-dup stop (%zu,%s) chord=%zd\n",
                             j, aj_v.edge, aj_v.side==LEFT?"L":"R",
                             (ptrdiff_t)aj_v.chord_idx);
#endif
                if (aj_v.chord_idx != NONE) {
                    assert(aj_v.chord_idx <
                           state.invalidated_walker_chords.size());
                    state.invalidated_walker_chords[aj_v.chord_idx] = true;
                }
                if (aj_v.is_companion) return;      // p = a_{m+1}
                state.p = fv_point(aj_v);
                state.p_edge = aj_v.edge;
                state.p_side = aj_v.side;
                state.p_y = aj_v.y;
                k = j + 1;
                break;
            }

            // [C91 §3.1 tex 206 case (i)]: record a_j → s_hit; p = a_j;
            // R unchanged; k = j+1 (A_k tie-break).
            if (auto r = case_i_test(j); r.fires) {
                const FusionVertex& aj_v = state.sequence[j];
                Point aj_point = fv_point(aj_v);
#ifdef CHAZELLE_TRACE_FUSION
                std::fprintf(stderr, "[fuse]  j=%zu case(i) stop=(%zu,%s,y=%g@%zu) hit=(%zu,%s,x=%g)\n",
                             j, aj_v.edge, aj_v.side==LEFT?"L":"R", aj_v.y.y,
                             aj_v.y.tag, r.s_hit.edge,
                             r.s_hit.side==LEFT?"L":"R", r.s_hit.x);
#endif

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
                    assert(aj_v.chord_idx < state.invalidated_walker_chords.size());
                    state.invalidated_walker_chords[aj_v.chord_idx] = true;
                }

                state.p = aj_point;
                state.p_edge = aj_v.edge;
                state.p_side = aj_v.side;
                state.p_y = aj_v.y;

                // [C91 §3.1 tex 195] invariant (B) tie rule: "If p lies
                // on a chord between two regions of S₂, then the
                // current region should be the one that we enter as we
                // locally leave p clockwise around ∂C₁."  When s lands
                // EXACTLY on an S₂ chord endpoint at the sightline's
                // level, the new p's sight rides that chord — tex 206's
                // "let the current region still be R" is the generic
                // action; the tie rule refines it here, exactly as the
                // startup's election does (tex 191).  Without it the
                // walk can keep R on the wrong side of the chord and
                // the Lemma 3.1 trichotomy breaks at a later stop (its
                // proof's "a' cannot be equal to p" presumes the
                // election).
                for (std::size_t ci :
                     S2.node(R).incident_chords) {
                    const Chord& ch = S2.chord(ci);
                    if (ch.dead || ch.is_null_length) continue;
                    if (!symbolic_y_equal(ch.symbolic_y(), aj_v.y))
                        continue;
                    const bool at_endpoint =
                        (ch.left_edge == r.s_hit.edge &&
                         ch.left_side == r.s_hit.side) ||
                        (ch.right_edge == r.s_hit.edge &&
                         ch.right_side == r.s_hit.side);
                    if (!at_endpoint) continue;
                    std::size_t below_r = NONE, above_r = NONE;
                    S2.chord_regions_below_above(ci, C2, &below_r,
                                                 &above_r);
                    // R is a per-iteration snapshot; the walk re-reads
                    // state.s2_region right after this action's break.
                    state.s2_region =
                        leaving_clockwise_downward(C1, aj_v.edge,
                                                   aj_v.side, aj_v.y)
                            ? below_r
                            : above_r;
#ifdef CHAZELLE_TRACE_FUSION
                    std::fprintf(stderr,
                                 "[fuse]   invariant-B tie: p on S2 "
                                 "chord %zu -> R=%zu\n", ci,
                                 state.s2_region);
#endif
                    break;
                }
                k = j + 1;
                break;
            }
            // [C91 §3.1 tex 206 case (ii)]: record q' → p'; p = p';
            // current = R's neighbor across chord ab; k = j.
            if (auto r = case_ii_test(j); r.fires) {
                const Chord& chord_ab = S2.chord(r.chord_idx);
#ifdef CHAZELLE_TRACE_FUSION
                std::fprintf(stderr, "[fuse]  j=%zu case(ii) chord=%zu p'=(%zu,%s,x=%g)\n",
                             j, r.chord_idx, r.p_prime_hit.edge,
                             r.p_prime_hit.side==LEFT?"L":"R",
                             r.p_prime_hit.x);
#endif
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
                // [C91 §3.1 tex 224]: the junction inside pair's chord
                // is the synthesized null-length chord — not recorded.
                if (!r.suppress_record) {
                    if (q_point.x < p_prime.x)
                        state.chords.push_back({chord_y, q_edge, q_side,
                            r.p_prime_hit.edge, r.p_prime_hit.side,
                            /*left_on_walker=*/false,
                            /*right_on_walker=*/true});
                    else
                        state.chords.push_back({chord_y,
                            r.p_prime_hit.edge, r.p_prime_hit.side,
                            q_edge, q_side,
                            /*left_on_walker=*/true,
                            /*right_on_walker=*/false});
                }

                // [C91 §3.1 tex 224]: q sees p' on A_j now → mark S₂'s
                // chord_ab superseded.
                assert(r.chord_idx < state.invalidated_target_chords.size());
                state.invalidated_target_chords[r.chord_idx] = true;

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
                // which we enter as we locally cross the exit chord at
                // p' along ∂C₁."  p' lies ON ab (the spherical chord —
                // for a wrapping ab the on-ab test above admits hits on
                // its through-infinity span, [C91 §2.1 tex 70]), whose
                // two sides are its above- and below-region along its
                // WHOLE length ([C91 §2.2 tex 96]: a chord separates
                // exactly its two regions).  The walker's continuation
                // past p' moves with its face's clockwise traversal
                // y-direction ([C91 §2.4(iii) tex 138]: LEFT follows
                // the edge, RIGHT reverses it), entering the above-
                // region when ascending.  (Flipping from R would be
                // wrong when the walker approached p' from the far
                // side of the chord's line.)
                {
                    const auto& pe = C1.edge(r.p_prime_hit.edge);
                    const bool e_asc = symbolic_y_less(
                        symbolic_y_of(C1.vertex(pe.start_idx)),
                        symbolic_y_of(C1.vertex(pe.end_idx)));
                    const bool leaving_up =
                        (r.p_prime_hit.side == LEFT) == e_asc;
                    std::size_t below = NONE, above = NONE;
                    S2.chord_regions_below_above(r.chord_idx, C2,
                                              &below, &above);
                    state.s2_region = leaving_up ? above : below;
                    assert((state.s2_region == chord_ab.region[0] ||
                            state.s2_region == chord_ab.region[1]) &&
                           "[C91 §3.1 tex 206]: the entered region is "
                           "one of the exit chord's two sides");
                }
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

    // ── Step 1: Chord inventory (per [C91 §3.1 tex 224]) ────────

    using Pending = PendingChord;
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
        return bit_at(state1.invalidated_walker_chords, i) ||
               bit_at(state2.invalidated_target_chords, i);
    };
    auto s2_invalid = [&](std::size_t i) {
        return bit_at(state1.invalidated_target_chords, i) ||
               bit_at(state2.invalidated_walker_chords, i);
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

        // Branch x-positions infinitesimally off v ([C91 §2 tex 47];
        // null padding branches included, polygon.h).
        const bool prev_left = extremum_prev_branch_left(u, v, w);

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

    // ── Steps 1b–4: shared normal-form realization ──────────────
    build_submap_from_chords(out_S, C, std::move(pending));
}

// [C91 §3.1 tex 226] + [C91 §2.2 tex 90–96]: realize the normal-form
// submap of a chord inventory over C (see fusion.h).
void build_submap_from_chords(Submap& out_S, const Polygon& C,
                              std::vector<PendingChord> pending) {
    assert(out_S.num_nodes() == 0 && out_S.num_arcs() == 0 &&
           out_S.num_chords() == 0 &&
           "[C91 §3.1 tex 226]: builder requires a fresh output submap");

    using Pending = PendingChord;
    const std::size_t n_edges = C.num_edges();

#ifdef CHAZELLE_EXPENSIVE_ASSERTS
    // [C91 §2.2 tex 90–96]: a submap's chords are visible pairs with
    // respect to C — each non-null chord's endpoints must see each
    // other ([C91 §2.1 tex 74]).  Gated: O(k·m) brute-force shots.
    for (const Pending& p : pending) {
        if (p.is_null_length) continue;
        const double xl = edge_x_at_y(C, p.left_edge_c, p.y);
        const double xr = edge_x_at_y(C, p.right_edge_c, p.y);
        const Side dl = shooting_direction(p.left_edge_c, p.left_side, C);
        Point pl{xl, p.y.y, p.y.tag};
        RayHit h = naive_first_contact(C, pl, p.y, dl,
                                       /*source_edge=*/p.left_edge_c);
        [[maybe_unused]] bool ok =
            h.hit && h.x == xr &&
            (h.edge == p.right_edge_c ||
             (symbolic_y_equal(p.y, symbolic_y_of(C.vertex(std::max(
                  h.edge, p.right_edge_c)))) &&
              (h.edge + 1 == p.right_edge_c ||
               p.right_edge_c + 1 == h.edge)));
        if (!ok) {
            std::fprintf(stderr,
                         "BAD CHORD y=%g@%zu L=(%zu,%s,x=%g) "
                         "R=(%zu,%s,x=%g) first-contact=(%zu,%s,x=%g)\n",
                         p.y.y, p.y.tag, p.left_edge_c,
                         p.left_side == LEFT ? "L" : "R", xl,
                         p.right_edge_c,
                         p.right_side == LEFT ? "L" : "R", xr,
                         h.edge, h.side == LEFT ? "L" : "R", h.x);
        }
        assert(ok &&
               "[C91 §2.2 tex 90]: every inventory chord must join a "
               "mutually visible pair with respect to C");
    }
#endif

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
    // Canonical endpoint labels first.  [C91 §3.0(i) tex 169]: a report
    // names "the edge of P that contains" the point hit, and a chord
    // endpoint at a NON-extremum vertex of C is contained in BOTH
    // incident edges — the same ∂C point can enter the inventory under
    // two names: pass 1's junction records label C₁∩C₂ with C₁'s last
    // edge while pass 2's label it with C₂'s first ([C91 §3.1 tex 179]),
    // and an oracle hit at a vertex-endpoint stop can disagree with the
    // stop's own Sᵢ slot label.  Rewrite to the LOWER incident edge so
    // the dedup key is label-free.  At a y-extremum the same-side
    // points on the two incident edges are DISTINCT duplication
    // companions ([C91 §2.1 tex 72]) and are left alone.
    {
        auto canonicalize_vertex_label = [&](std::size_t& edge_c,
                                             SymbolicY y) {
            if (edge_c == 0) return;
            const auto& e = C.edge(edge_c);
            // [C91 §2 tex 47]: SoS-exact — y matching the edge's
            // traversal-lower endpoint means the point IS that vertex
            // (the vertex shared with edge_c − 1).
            if (!symbolic_y_equal(y, symbolic_y_of(C.vertex(e.start_idx))))
                return;
            const std::size_t vidx = e.start_idx;
            assert(vidx == edge_c &&
                   "[C91 §2.4 tex 133]: edge k spans vertices k → k+1");
            assert(vidx + 1 < C.num_vertices());
            if (is_local_y_extremum(C.vertex(vidx - 1), C.vertex(vidx),
                                    C.vertex(vidx + 1)))
                return;
            edge_c = edge_c - 1;
        };
        for (Pending& p : pending) {
            canonicalize_vertex_label(p.left_edge_c, p.y);
            canonicalize_vertex_label(p.right_edge_c, p.y);
        }
    }

    // Canonical slot order next: slots are "ascending x"; a chord of
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
    std::vector<Endpoint> endpoints;
    endpoints.reserve(2 * pending.size());
    for (std::size_t i = 0; i < pending.size(); ++i) {
        const auto& p = pending[i];
        endpoints.push_back({p.left_edge_c,  p.left_side,  p.y, i, true});
        endpoints.push_back({p.right_edge_c, p.right_side, p.y, i, false});
    }
    // Linear ∂C order of a position: face rank, then symbolic y along
    // the face's traversal direction.
    auto lin_less = [&](std::size_t pa, SymbolicY ya, std::size_t pb,
                        SymbolicY yb) -> bool {
        if (pa != pb) return pa < pb;
        if (symbolic_y_equal(ya, yb)) return false;
        // trav rank pa encodes (edge, side); recover traversal sense.
        const std::size_t e = (pa < n_edges) ? pa : (2 * n_edges - 1 - pa);
        const Side sd = (pa < n_edges) ? LEFT : RIGHT;
        bool trav_asc = (sd == LEFT) == edge_ascends(e);
        return trav_asc ? symbolic_y_less(ya, yb)
                        : symbolic_y_greater(ya, yb);
    };
    auto partner_pos = [&](const Endpoint& e, std::size_t* pp,
                           SymbolicY* py) {
        const auto& p = pending[e.pending_idx];
        const std::size_t pe = e.is_left_slot ? p.right_edge_c
                                              : p.left_edge_c;
        const Side ps = e.is_left_slot ? p.right_side : p.left_side;
        *pp = trav_pos(pe, ps);
        *py = p.y;
    };
    std::sort(endpoints.begin(), endpoints.end(),
        [&](const Endpoint& a, const Endpoint& b) {
            std::size_t pa = trav_pos(a.edge_c, a.side);
            std::size_t pb = trav_pos(b.edge_c, b.side);
            if (lin_less(pa, a.y, pb, b.y)) return true;
            if (lin_less(pb, b.y, pa, a.y)) return false;
            // Same ∂C position ([C91 §2 tex 47]: a padded curve's
            // endpoint-turnaround cluster stacks several chord
            // endpoints at one wall point, [C91 §4 tex 316]).  Lemma
            // 2.2's parenthesis system ([C91 §2.2 tex 102]) decides
            // the local order by the partners: an endpoint whose
            // partner lies linearly BEFORE this position closes its
            // chord and precedes every opener; within closers the
            // inner (later-opened) chord closes first, within openers
            // the outer (farther-reaching) chord opens first — both
            // are descending partner order.
            if (a.pending_idx == b.pending_idx) return a.is_left_slot;
            std::size_t qa, qb;
            SymbolicY qya, qyb;
            partner_pos(a, &qa, &qya);
            partner_pos(b, &qb, &qyb);
            const bool a_closer = lin_less(qa, qya, pa, a.y);
            const bool b_closer = lin_less(qb, qyb, pb, b.y);
            if (a_closer != b_closer) return a_closer;
            return lin_less(qb, qyb, qa, qya);
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
    std::vector<std::size_t> awaiting_after_arc;

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
        for (std::size_t pi : awaiting_after_arc) {
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
        awaiting_after_arc.clear();
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
        const Endpoint& g = endpoints[i];
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
        //
        // [C91 §2.4 tex 133 + §2.2 tex 96]: chord endpoints COINCIDENT
        // at one ∂C position (a padded curve's turnaround cluster,
        // [C91 §4 tex 316]) are separated along ∂C by zero-length
        // arcs — "some arcs may be of zero length", and an arc never
        // contains a chord attachment, even a coincident one
        // ([C91 §2.3 tex 106]).  `separator_ready` tracks whether the
        // last emitted arc already bounds the current region at this
        // position (the group's closing arc for the first endpoint;
        // a null chord's own R_inner arc after its opening).
        bool separator_ready = true;
        for (std::size_t k = i; k < j; ++k) {
            const Endpoint& e = endpoints[k];
            if (!separator_ready) {
                std::size_t sep = emit_arc(e.edge_c, e.side,
                                           /*edge_count=*/0);
                patch_after_arcs(sep);
            }
            separator_ready = false;
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
                    separator_ready = true;   // the R_inner arc itself
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
                awaiting_after_arc.push_back(e.pending_idx);
            }
        }
    };

    std::size_t i = 0;
    bool tail_is_zero = false;
    while (i < endpoints.size()) {
        std::size_t j = i + 1;
        while (j < endpoints.size() &&
               trav_pos(endpoints[j].edge_c, endpoints[j].side) ==
                   trav_pos(endpoints[i].edge_c, endpoints[i].side) &&
               symbolic_y_equal(endpoints[j].y, endpoints[i].y)) ++j;
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
            const Endpoint& g = endpoints[i];
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
