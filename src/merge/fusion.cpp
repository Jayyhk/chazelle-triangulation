/// src/merge/fusion.cpp

#include "fusion.h"

#include <cmath>

namespace chazelle {

// ── collect_region_arcs ─────────────────────────────────────────

RegionArcs collect_region_arcs(const Submap& S, std::size_t region) {
    assert(region < S.num_nodes() && !S.node(region).dead);

    RegionArcs out;
    // [C91 §3.1] (tex 181): "The claim on the time follows from the
    // conformality of Sᵢ."  Assert local conformality before iteration
    // to strictly enforce O(1) evaluation bound.
    assert(S.node(region).degree() <= 4 &&
           "§2.3: conformal regions MUST have degree ≤ 4");
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
               "§2.4: normal-form (compacted) submap must have no "
               "dead chords in incident_chords");
        check_adj(S.chord(ci).left_adj);
        check_adj(S.chord(ci).right_adj);
    }

    // [C91 §2.4(iii)] (tex 138): "the endpoints of C are identified
    // by... pointers to the arc-structures whose corresponding arcs
    // pass through the endpoints."  Normal-form submaps must have
    // start_arc and end_arc set.  Assert instead of silently skipping.
    auto check_endpoint_arc = [&](std::size_t arc_idx) {
        assert(arc_idx != NONE && arc_idx < S.num_arcs() &&
               "§2.4(iii): endpoint arc must be set for normal-form submaps");
        assert(!S.arc(arc_idx).dead &&
               "§2.4(iii): endpoint arc must be live");
        if (S.arc(arc_idx).region_node != region) return;
        bool dup = false;
        for (std::size_t i = 0; i < out.count; ++i)
            if (out.arcs[i] == arc_idx) { dup = true; break; }
        if (!dup) out.push(arc_idx);
    };
    check_endpoint_arc(S.start_arc);
    check_endpoint_arc(S.end_arc);

    // [C91 §3.1] (tex 181): "the conformality of Sᵢ, which ensures that
    // at most four arcs need to be checked."  A degree-d region has exactly
    // one arc per chord-gap around its boundary, so degree ≤ 4 ⟹ ≤ 4 arcs.
    assert(out.count <= 4 &&
           "§3.1 (tex 181): conformal region has at most 4 arcs (degree ≤ 4)");

    return out;
}

// ── local_shoot ─────────────────────────────────────────────────

RayHit local_shoot(Point p, Side direction,
                    std::size_t region,
                    const Submap& S, const Polygon& C,
                    const RayShootingOracle& oracle) {
    // [C91 §3.1]: "Using the appropriate ray-shooters, we can find
    // that point by checking each arc in turn and finding the nearest
    // hit.  The claim on the time follows from the conformality of Sᵢ,
    // which ensures that at most four arcs need to be checked."
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

        // [C91 §3.0] (tex 166): "for each region arc α of Sᵢ
        // specified by a pointer to its arc-structure."
        RayHit hit = oracle.shoot(p, direction, ai, sub);
        if (!hit.hit) continue;

        // [C91 §3.0(i) tex 169]: "reports the single point of α' (if any)
        // that a ray of light shot from p in the given direction would hit."
        // Oracle contract is forward-only: hit.x must lie at or past p.x
        // along the shooting direction.  A backward hit would have negative
        // signed distance and "beat" every legitimate forward hit, silently
        // corrupting best-of-nearest selection.
        double hit_signed_dist = (direction == LEFT)
            ? (p.x - hit.x) : (hit.x - p.x);
        assert(hit_signed_dist >= 0.0 &&
               "§3.0(i) tex 169: ray-shooting oracle must report a forward "
               "hit (hit_signed_dist >= 0 in the shooting direction)");

        // [C91 §3.1]: Track the explicit arc index to propagate O(1) context
        hit.hit_arc_idx = ai;

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
                // [C91 §2.1] (tex 72): Double boundary disambiguation.
                // "We give each edge of C an infinitesimal width so as to make
                // the curve C into a very thin simple polygon… we refer to the
                // two sides of the double boundary as one would speak of the
                // left and right sides of a snake."  The two sides of the same
                // edge sit at the same geometric location but are topologically
                // distinct ∂C points; ray-shooting must select the side struck
                // first by the (infinitesimally thick) ray.  See also §2.4
                // (tex 142) on double-backing detection.
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

    // [C91 §3.1] (tex 181): "the point of ∂Cᵢ that it sees lies on
    // one of the region's arcs.  This is because regions are closed
    // under visibility, which is a corollary of Lemma 2.1."
    assert(best.hit &&
           "§3.1: local shooting inside a region must always hit "
           "one of the region's arcs (corollary of Lemma 2.1)");
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
    // [C91 §3.1]: The fusion sequence always contains at least a₀ and
    // a_{m+1} (the two companion vertices at the junction).
    assert(state.sequence.size() >= 2 &&
           "§3.1: fusion sequence must have at least a₀ and a_{m+1}");
    const auto& a0 = state.sequence[0];
    assert(a0.is_companion && a0.side == RIGHT);
    // [C91 §3.1 tex 179]: a_{m+1} is the LEFT companion at the junction,
    // end of the clockwise tour around ∂C₁.  Symmetric to a₀.
    const auto& a_m1 = state.sequence.back();
    assert(a_m1.is_companion && a_m1.side == LEFT &&
           "§3.1 tex 179: sequence.back() must be the a_{m+1} LEFT companion");

    // [C91 §3.1]: "Using local shooting, we find the point of ∂C₁
    // that a₀ sees with respect to C₁."
    //
    // a₀ is at c_end of C₁, RIGHT side.  The region in S₁ is the
    // one containing the end_arc (or its RIGHT-side neighbor).
    // [C91 §2.4(iii)]: end_arc must be set for normal-form submaps.
    assert(S1.end_arc != NONE &&
           "§2.4(iii): S₁ must have end_arc set (normal form)");
    std::size_t a0_region_s1 = S1.arc(S1.end_arc).region_node;
    Side a0_dir = shooting_direction(a0.edge, a0.side, C1);

    // a₀ IS the junction vertex (last vertex of C₁).  Use its
    // coordinates directly — no interpolation needed.
    std::size_t junction_v = C1.num_vertices() - 1;
    Point a0_point = C1.vertex(junction_v);

    RayHit hit_c1 = local_shoot(a0_point, a0_dir, a0_region_s1,
                                 S1, C1, oracle1);

    // [C91 §3.1]: "it is not quite true that a₀ is always a point
    // of ∂C₂.  It coincides with one most often, but when it sits at
    // a local extremum (in the y-direction) it is not one because of
    // duplication."  SoS perturbation (§2.0) handles this: distinct
    // symbolic y-coordinates prevent exact coincidence at extrema,
    // so the visibility relation is always well-defined.

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
    // [C91 §3.1] (tex 185): Both local_shoot calls are guaranteed to
    // hit (Lemma 2.1 corollary — regions are closed under visibility).
    // c₀ is the nearer of hit_c1 and hit_c2 in the shooting direction.
    assert(hit_c1.hit &&
           "§3.1: a₀ must see a point on ∂C₁ (Lemma 2.1)");
    assert(hit_c2.hit &&
           "§3.1: a₀ must see a point on ∂C₂ (Lemma 2.1)");
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
        // [C91 §3.1] (tex 191): If a₀ lies on an S₂ chord, that chord
        // must be incident on initial_region.  Conformal ⟹ degree ≤ 4,
        // so this scan is O(1) — matching the paper's O(1) claim (tex 185).
        SymbolicY a0y = a0.y;
        for (std::size_t ci : S2.node(initial_region).incident_chords) {
            assert(!S2.chord(ci).dead &&
                   "§2.4: normal-form (compacted) submap must have no "
                   "dead chords in incident_chords");
            if (!symbolic_y_equal(S2.chord(ci).symbolic_y(), a0y))
                continue;

            // a₀c₀ lies on this S₂ chord.  Determine which of its
            // two regions we "locally enter as we leave p in a
            // clockwise traversal of ∂C₁."
            //
            // We parameterize `leaving_downward` because the topological
            // vector leaving `p` depends entirely on whether we are in
            // Case 1 (leaving a₀) or Case 2 (leaving c₀).

            // Determine which region is above/below the chord by
            // comparing adj arcs' key_y to chord.y.  Check both
            // endpoints — left_adj or right_adj — whichever has 2
            // adj arcs (non-vertex endpoint per §2.2 tex 94).
            const auto& ch = S2.chord(ci);
            const Chord::AdjArcs* adj = nullptr;
            if (ch.left_adj.count == 2)
                adj = &ch.left_adj;
            else if (ch.right_adj.count == 2)
                adj = &ch.right_adj;

            if (adj) {
                // [C91 §2.4] (tex 137): At a non-vertex endpoint exactly one
                // arc STARTS here (key_y == chord.y) and one ENDS here
                // (key_y != chord.y).  Always use the ending arc for
                // above/below determination — its key_y is unambiguously
                // above or below the chord on both LEFT and RIGHT sides of ∂C.
                //
                // Using the starting arc is WRONG for RIGHT-side endpoints:
                // its key_y == chord.y yields symbolic_y_less = false ("not
                // below" → classified above), but a RIGHT-side starting arc
                // descends downward, so its region is actually BELOW.
                SymbolicY ky0 = S2.arc(adj->arcs[0]).key_symbolic_y();
                SymbolicY ky1 = S2.arc(adj->arcs[1]).key_symbolic_y();
                std::size_t r0 = S2.arc(adj->arcs[0]).region_node;
                std::size_t r1 = S2.arc(adj->arcs[1]).region_node;
                // Exactly one arc starts here (key_y == chord.y).
                bool arc0_is_starting = symbolic_y_equal(ky0, ch.symbolic_y());
                bool arc1_is_starting = symbolic_y_equal(ky1, ch.symbolic_y());
                assert(arc0_is_starting != arc1_is_starting &&
                       "§2.4: exactly one adj arc must start at the chord "
                       "endpoint (key_y == chord.y); one starts, one ends");
                // Use the ending arc (key_y != chord.y) for direction.
                bool use_arc1 = arc0_is_starting;
                SymbolicY ky_ref = use_arc1 ? ky1 : ky0;
                std::size_t r_ref   = use_arc1 ? r1  : r0;
                std::size_t r_other = use_arc1 ? r0  : r1;
                bool ref_below = symbolic_y_less(ky_ref, ch.symbolic_y());
                std::size_t below_r = ref_below ? r_ref   : r_other;
                std::size_t above_r = ref_below ? r_other : r_ref;
                return leaving_downward ? below_r : above_r;
            }
            // Both endpoints are C₂ vertices (count==1 on both sides).
            // Under SoS (§2, tex 47), the junction vertex has a unique
            // symbolic y.  The junction is C₂'s start vertex (an
            // endpoint), so §2.1 guarantees no NLC there.  No S₂ chord
            // can match the junction's symbolic y — this path is
            // unreachable.
            assert(false &&
                   "§2.1/SoS: S₂ chord at junction's symbolic y is "
                   "impossible (junction is C₂ endpoint, no NLC)");
        }
        return initial_region;
    };

    if (c0_on_c2) {
        // [C91 §3.1 Case 1]: "If c₀ belongs to ∂C₂, then we set
        // p = a₀ and we call the region of S₂ crossed by a₀c₀
        // current: the start-up phase is over."
        //
        // Main loop starts at k where A_k contains p = a₀ and
        // p ≠ a_k.  A_k = arc from a₀ to a₁, so k=1.
        // [C91 §3.1 Case 1]: "we set p = a₀"
        state.p = a0_point;
        state.p_edge = a0.edge;
        state.p_side = a0.side;

        // From a₀ (RIGHT companion at c_end), clockwise ∂C₁
        // goes from vertex junction_v toward vertex junction_v-1.
        // The y-direction of that step mathematically defines leaving p.
        std::size_t junction_v = C1.num_vertices() - 1;
        SymbolicY y_here = symbolic_y_of(C1.vertex(junction_v));
        SymbolicY y_next = symbolic_y_of(C1.vertex(junction_v - 1));
        bool a0_leaving_downward = symbolic_y_less(y_next, y_here);

        state.s2_region = resolve_s2_region(a0_region_s2, a0_leaving_downward);
        // [C91 §2.2] (tex 96): exit chords connect opposite sides of the double boundary.
        // a0.side == RIGHT (asserted above); assert guarantees c0.side == LEFT.
        assert(c0.side != a0.side &&
               "§2.2: chord must connect opposite sides of the double boundary");
        assert(c0.side == LEFT && a0.side == RIGHT);
        state.chords.push_back({a0.y, c0.edge, c0.side, a0.edge, a0.side});
        return 1;
    } else {
        // [C91 §3.1 Case 2]: "If c₀ belongs to ∂C₁... We can
        // therefore skip all the way to c₀.  Now, however, c₀ sees
        // a point of ∂C₂, namely a₀, so we set p = c₀ and call the
        // region of S₂ containing a₀ current."
        // [C91 §3.1 Case 2]: "we set p = c₀"
        // c₀ is on an edge interior (not necessarily a vertex), so
        // Point.index = NONE (no vertex identity for SoS).
        state.p = Point{c0.x, c0.y, NONE};
        state.p_edge = c0.edge;
        state.p_side = c0.side;

        // [C91 §3.1 Case 2]: Determine topological slope leaving c₀
        assert(c0.edge < C1.num_edges());
        const auto& e = C1.edge(c0.edge);
        bool edge_ascending = symbolic_y_less(
            symbolic_y_of(C1.vertex(e.start_idx)),
            symbolic_y_of(C1.vertex(e.end_idx)));
        
        bool c0_leaving_downward = (c0.side == LEFT) ? !edge_ascending : edge_ascending;

        // [C91 §3.1 Case 2]: "call the region of S₂ containing a₀ current."
        // We unconditionally resolve the region, because while a₀ is a start vertex,
        // it may be a reflex vertex generated by an arbitrary interior cut splitting,
        // allowing a horizontal S₂ chord to originate from the exact same Y coordinate.
        state.s2_region = resolve_s2_region(a0_region_s2, c0_leaving_downward);
        // [C91 §2.2] (tex 96): exit chords connect opposite sides of the double boundary.
        // a0.side == RIGHT (asserted above); assert guarantees c0.side == LEFT.
        assert(c0.side != a0.side &&
               "§2.2: chord must connect opposite sides of the double boundary");
        assert(c0.side == LEFT && a0.side == RIGHT);
        state.chords.push_back({a0.y, c0.edge, c0.side, a0.edge, a0.side});

        // [C91 §3.1 Case 2]: "We can therefore skip all the way to c₀."
        // Find the first fusion vertex at or past c₀'s position on ∂C₁.
        // All vertices between a₀ and c₀ (clockwise) see points on ∂C₁
        // (by Lemma 2.1) and are already known from S₁.
        //
        // [C91 §3.1 Case 2]: "We can therefore skip all the way to c₀."
        // O(1) via arc_starts lookup (the hit arc from local_shoot indexes
        // directly into the fusion sequence) + O(1) within-bucket scan
        // (conformal degree ≤ 4 → at most O(1) endpoints per arc bucket).
        std::size_t junction_edge = C1.num_edges() - 1;
        std::size_t right_half_len = junction_edge + 1;
        auto trav_pos = [&](std::size_t edge, Side side) -> std::size_t {
            if (side == RIGHT)
                return junction_edge - edge;
            else
                return right_half_len + edge;
        };

        // c₀ lies on ∂C₁ at the same horizontal level as a₀ (they share a
        // horizontal visibility chord).  c₀ is a point on an edge interior,
        // not a named vertex, so it has no vertex index of its own.
        // We assign it the SoS tag of the junction vertex (the vertex that
        // defines this y-level: a₀ is the companion vertex at the junction,
        // so a₀.y.tag IS the junction vertex's index).  This is correct
        // because every comparison downstream (vertex_past_c0, symbolic_y_geq/
        // symbolic_y_leq) pairs c₀'s symbolic-y against fusion vertices whose
        // tags are all distinct from the junction vertex's tag — they are
        // non-companion chord endpoints.  The junction tag is therefore
        // never involved in a tie-break within the comparison, and the
        // raw-y equality of the horizontal chord is never ambiguous under SoS.
        SymbolicY c0_y = a0.y;
        std::size_t c0_pos = trav_pos(c0.edge, c0.side);

        // Determine y-ordering within an edge on ∂C₁ traversal.
        // RIGHT half traverses edges in descending order; within an
        // edge, ∂C₁ goes from end_vertex to start_vertex (reversed).
        // LEFT half follows edge direction (start to end).
        auto vertex_past_c0 = [&](const FusionVertex& v) -> bool {
            std::size_t v_pos = trav_pos(v.edge, v.side);
            if (v_pos != c0_pos) return v_pos > c0_pos;
            // Same edge and side: compare y along traversal direction.
            // On RIGHT side, ∂C₁ traverses the edge in reverse, so
            // "past c₀" means lower y (if edge ascends) or higher y
            // (if edge descends).  On LEFT side, it follows the edge.
            assert(v.edge < C1.num_edges());
            const auto& e = C1.edge(v.edge);
            bool edge_ascending = symbolic_y_less(
                symbolic_y_of(C1.vertex(e.start_idx)),
                symbolic_y_of(C1.vertex(e.end_idx)));
            // LEFT: traversal follows edge direction.
            //   ascending edge → later in traversal = higher y.
            //   v is past c₀ if v.y ≥ c₀.y (ascending) or v.y ≤ c₀.y (descending).
            // RIGHT: traversal is reversed.
            bool traversal_ascending =
                (v.side == LEFT) ? edge_ascending : !edge_ascending;
            if (traversal_ascending)
                return symbolic_y_geq(v.y, c0_y);
            else
                return symbolic_y_leq(v.y, c0_y);
        };

        // [C91 §3.1 Case 2]: "We can therefore skip all the way to c₀."
        // We evaluate this natively in O(1) time. c0 was struck on an explicitly
        // tracked arc in S1 via local_shoot. We know precisely which block it occupies.
        assert(hit_c1.hit_arc_idx != NONE && "§3.1: c0 must carry hit context");
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

        // [C91 §3.1 tex 179, 188]: c₀ is the visibility hit from a₀ on ∂C₁ and
        // hence strictly precedes a_{m+1} on the clockwise tour.  a_{m+1} sits
        // at the maximal trav_pos (LEFT side at the junction's SymbolicY), so
        // vertex_past_c0(a_{m+1}) is always true and the while loop must exit
        // with lo ≤ sequence.size() − 1.
        assert(lo < state.sequence.size() &&
               "§3.1 tex 179, 188: a_{m+1} sits at maximal trav_pos so "
               "vertex_past_c0(a_{m+1}) holds; the loop above must terminate "
               "with lo ≤ sequence.size() − 1");

        // §3.1 Case 2 (Lemma 2.1): All skipped vertices see ∂C₁.
        // This is a paper-guaranteed structural property; only verify
        // under CHAZELLE_EXPENSIVE_ASSERTS to avoid inflating startup
        // from O(f(γ₂)) to O(m) in normal debug builds.
#ifdef CHAZELLE_EXPENSIVE_ASSERTS
        for (std::size_t skipped_i = 1; skipped_i < lo; ++skipped_i) {
            if (!state.sequence[skipped_i].is_companion) {
                assert(state.sequence[skipped_i].chord_idx != NONE &&
                       "§3.1 Case 2 (Lemma 2.1): skipped vertex must be a "
                       "chord endpoint of S₁ (sees ∂C₁ by construction)");
            }
        }
#endif

        return lo;
    }
}

// ── fuse_s1_into_s2 ─────────────────────────────────────────────

void fuse_s1_into_s2(FusionState& state,
                      const Submap& S1, const Polygon& C1,
                      const Submap& S2, const Polygon& C2,
                      const RayShootingOracle& oracle1,
                      const RayShootingOracle& oracle2) {
    // [C91 §3.1]: Build the fusion vertex sequence.
    build_fusion_sequence(state, S1, C1);
    state.current_stop = 0;

    // [C91 §3.1 Start-Up]: Initialize p and current S₂ region.
    std::size_t start_idx = fusion_startup(state, S1, C1, S2, C2,
                                            oracle1, oracle2);

    // [C91 §3.1 invariant (A)]: "The points of ∂C that are seen by the exit
    // chord endpoints of S₁ ... have all been determined already."
    assert(!state.chords.empty() &&
           "§3.1 invariant (A): startup must have produced at least one chord");

    // [C91 §3.1 invariant (B)]: "The point q of ∂C that is seen by p belongs
    // to ∂C₂ and the chord pq lies in the region of S₂ called current."
    // local_shoot is O(f(γ₂)) — gated to avoid doubling startup's own bound.
#ifdef CHAZELLE_EXPENSIVE_ASSERTS
    assert(local_shoot(state.p,
               shooting_direction(state.p_edge, state.p_side, C1),
               state.s2_region, S2, C2, oracle2).hit &&
           "§3.1 invariant (B): p must see ∂C₂ after startup");
#endif

    // [C91 §3.1]: "We let a variable p run through ∂C₁ in clockwise
    // order, stopping at a₀, ..., a_{m+1}."
    for (std::size_t i = start_idx; i < state.sequence.size(); ++i) {
        state.current_stop = i;

        // [C91 §3.1 invariant (B)]: "The point q of ∂C that is seen
        // by p belongs to ∂C₂ and the chord pq lies in the region
        // of S₂ called current."
        assert(state.s2_region != NONE &&
               state.s2_region < S2.num_nodes() &&
               !S2.node(state.s2_region).dead &&
               "§3.1 invariant (B): current S₂ region must be valid");
        // [C91 §3.1 invariant (B)] (tex 195): "The point q of ∂C
        // that is seen by p belongs to ∂C₂ and the chord pq lies in
        // the region of S₂ called current."
        // local_shoot is O(f(γ₂)) per call — running it every iteration
        // would inflate the main loop from O(m + …) to O(m·f(γ₂)) in
        // debug builds, violating the paper's per-merge bound.  Gate
        // under CHAZELLE_EXPENSIVE_ASSERTS (matches the §3.1 startup
        // skipped-vertex pattern at line ~502).
#ifdef CHAZELLE_EXPENSIVE_ASSERTS
        assert(local_shoot(
                   state.p,
                   shooting_direction(state.p_edge, state.p_side, C1),
                   state.s2_region, S2, C2, oracle2).hit &&
               "§3.1 invariant (B): p must see a point q on ∂C₂");
#endif

        // [C91 §3.1 invariant (A)]: "The points of ∂C that are seen
        // by the exit chord endpoints of S₁ on the portion of ∂C₁
        // running clockwise from a₀ to the point p in its current
        // position have all been determined already."
        //
        // All fusion vertices before index i have been processed
        // (either by startup or previous iterations).  Their chords
        // are in state.chords.
        // Count non-companion vertices before i that should have
        // been resolved.
        // At minimum, the startup produced at least one chord.
        if (i > 0)
            assert(!state.chords.empty() &&
                   "§3.1 invariant (A): startup must have "
                   "produced at least one chord");

        // TODO: (§3.1) At each stop:
        //   1. Determine what p sees on ∂C₂ (local_shoot into S₂ taking O(f(γ₂))).
        //   2. Determine what p sees on ∂C₁ (CRITICAL: if p has a chord_idx, do NOT use
        //      local_shoot into S₁! Just return the opposite chord endpoint in O(1) time
        //      per §3.1 deviation rules. Only fallback to local_shoot for companion vertices).
        //   3. Record discovered chord and update s2_region for the next stop.
    }
}

// ── build_fusion_sequence ───────────────────────────────────────

void build_fusion_sequence(FusionState& state, const Submap& S,
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
    // Clockwise ∂C₁: RIGHT descending (junction → c_start) → LEFT ascending (c_start → junction).
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
    // [C91 §2.4(iii)] (tex 142): "our representation is powerful enough
    // to let us perform canonical vertex/region enumerations in optimal
    // time."  The arc-sequence table is in ∂C order (LEFT ascending,
    // then RIGHT descending).  Clockwise ∂C₁ from a₀ = RIGHT section
    // then LEFT section.  We enumerate chord endpoints in O(m) via
    // counting sort keyed on the clockwise arc-table position.
    //
    // Each chord endpoint is associated with an arc via adj_arcs:
    //   count=2 (non-vertex endpoint): the "starting" arc whose key_y
    //     equals the chord's y.
    //   count=1 (vertex endpoint): the single containing arc.

    std::size_t num_arcs = S.num_arcs();
    std::size_t lrb = S.left_right_boundary();
    std::size_t num_right = num_arcs - lrb;

    // Map arc-table index → clockwise position.
    // Clockwise: RIGHT arcs (table[lrb..end]) then LEFT (table[0..lrb-1]).
    auto cw_pos = [&](std::size_t arc_idx) -> std::size_t {
        if (arc_idx >= lrb)
            return arc_idx - lrb;             // RIGHT: [0, num_right)
        else
            return num_right + arc_idx;       // LEFT:  [num_right, num_arcs)
    };

    // Collect chord endpoints with their clockwise arc key.
    struct KeyedVertex {
        FusionVertex v;
        std::size_t key;   // clockwise arc position
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

    // For a count=2 endpoint, the "starting" arc is the one whose
    // key_y equals the chord's y (it begins at the chord endpoint).
    auto starting_arc = [&](const Chord::AdjArcs& adj,
                             SymbolicY chord_y) -> std::size_t {
        assert(adj.count == 2);
        if (symbolic_y_equal(S.arc(adj.arcs[0]).key_symbolic_y(), chord_y))
            return adj.arcs[0];
        // [C91 §2.4]: For a non-vertex endpoint, exactly one of the two
        // adjacent arcs starts at the chord endpoint (key_y == chord_y).
        assert(symbolic_y_equal(S.arc(adj.arcs[1]).key_symbolic_y(), chord_y) &&
               "§2.4: one of the two adj arcs must start at the chord endpoint");
        return adj.arcs[1];
    };

    for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
        assert(!S.chord(ci).dead &&
               "§2.4: normal-form (compacted) submap must have no "
               "dead chords");
        const auto& c = S.chord(ci);
        // [C91 §3.1] (tex 179, tex 142): "the canonical vertex enumeration
        // enumerates the exit chord endpoints."  Exit chords are explicitly
        // distinguished from null-length chords (tex 96).  NLCs are carried
        // over automatically to the fused submap (tex 224) and must not enter
        // the shooting traversal.  Including them would inflate the sequence
        // from O(n/γ + 1) to O(n), violating Lemma 3.1 (tex 209-210).
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
    // [C91 §2.4] (tex 142): "canonical vertex enumerations in optimal
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

    // Within each bucket (same arc), O(1) elements for conformal
    // submaps (degree ≤ 4).  Sort by y along the arc's traversal
    // direction.  Use insertion sort — O(1) per bucket, O(m) total.
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

    // [C91 §3.1]: "we assume that [a₀ and a_{m+1}] are not [in the
    // sequence] and therefore add them to the sequence."
    // Prepend a₀, append a_{m+1}.
    std::vector<FusionVertex> result;
    result.reserve(sorted.size() + 2);
    result.push_back(a_0);
    for (const auto& kv : sorted)
        result.push_back(kv.v);
    result.push_back(a_m1);

    state.sequence = std::move(result);

    // [C91 §3.1]: Map arc keys directly to their contiguous bounds in sequence.
    // bucket_starts.size() == num_arcs by construction (both branches above),
    // so every slot is assigned — no NONE sentinels remain.
    state.arc_starts.resize(num_arcs);
    for (std::size_t i = 0; i < num_arcs; ++i) {
        state.arc_starts[i] = bucket_starts[i] + 1; // +1 because a_0 is prepended
    }
}

} // namespace chazelle
