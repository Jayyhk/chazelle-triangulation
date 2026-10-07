#include "fusion.h"
#include "../trace.h"

#include <algorithm>

namespace chazelle::animation {

static void record_discovery(const FusionState& state, const Polygon& first,
                             const Polygon& second) {
    if (auto* trace = AnimationTrace::current()) {
        const Polygon combined =
            state.junction_at_end ? Polygon(first, second) : Polygon(second, first);
        const std::size_t first_offset = state.junction_at_end ? 0 : second.num_edges();
        const std::size_t second_offset = state.junction_at_end ? first.num_edges() : 0;
        const auto& discovered = state.chords.back();
        Chord chord;
        chord.y = discovered.y.y;
        chord.y_tag = discovered.y.tag;
        chord.left_edge =
            discovered.left_edge + (discovered.left_on_first_curve ? first_offset : second_offset);
        chord.left_side = discovered.left_side;
        chord.right_edge = discovered.right_edge +
                           (discovered.right_on_first_curve ? first_offset : second_offset);
        chord.right_side = discovered.right_side;
        trace->chord("fusion_chord", combined, chord);
    }
}

static bool leaving_clockwise_downward(const Polygon& curve, std::size_t edge, Side side,
                                       const SymbolicY& y) {
    auto trav_asc = [&](std::size_t e, Side s_) {
        const auto& ed = curve.edge(e);
        const bool asc = symbolic_y_less(symbolic_y_of(curve.vertex(ed.start_idx)),
                                         symbolic_y_of(curve.vertex(ed.end_idx)));
        return (s_ == LEFT) ? asc : !asc;
    };
    const auto& ed = curve.edge(edge);
    const std::size_t trav_end = (side == LEFT) ? ed.end_idx : ed.start_idx;
    if (!symbolic_y_equal(y, symbolic_y_of(curve.vertex(trav_end))))
        return !trav_asc(edge, side);

    if (side == LEFT) {
        if (edge + 1 < curve.num_edges())
            return !trav_asc(edge + 1, LEFT);
        return !trav_asc(edge, RIGHT);
    }
    if (edge > 0)
        return !trav_asc(edge - 1, RIGHT);
    return !trav_asc(edge, LEFT);
}

RayHit local_shoot(const Point& p, Side direction, std::size_t region, const Submap& submap,
                   const Polygon& curve, const RayShootingOracle& oracle, bool require_hit,
                   const SourceOffset& source_x_offset, bool record) {
    const AnimationTrace::QueryRecording recording(record);
    RegionArcs arcs = collect_region_arcs(submap, region);

    RayHit best;
    best.hit = false;

    Subarc subs[RegionArcs::MAX];

    for (std::size_t k = 0; k < arcs.count; ++k) {
        const std::size_t ai = arcs.arcs[k];
        const auto& a = submap.arc(ai);

        Subarc& sub = subs[k];
        sub.first_edge = a.first_edge;
        sub.first_side = a.first_side;
        sub.last_edge = a.last_edge;
        sub.last_side = a.last_side;

        sub.first_y = submap.arc_start_symbolic_y(ai, curve);
        sub.last_y = submap.arc_end_symbolic_y(ai, curve);
        assert_subarc_clockwise(sub);

        RayHit hit = oracle.shoot(p, direction, ai, sub, source_x_offset);
        if (!hit.hit)
            continue;

        Exact hit_signed_dist = (direction == LEFT) ? (p.x - hit.x) : (hit.x - p.x);

        assert(
            (hit_signed_dist < 0.0 ? hit.wrapped : (hit_signed_dist > 0.0 ? !hit.wrapped : true)) &&
            "[C91 §2.1 tex 70]: wrap flag consistent with the signed "
            "travel distance");

        hit.hit_arc_idx = ai;

        if (!best.hit) {
            best = hit;
        } else if (hit.wrapped != best.wrapped) {
            if (!hit.wrapped)
                best = hit;
        } else {
            Exact best_dist = (direction == LEFT) ? (p.x - best.x) : (best.x - p.x);
            Exact hit_dist = (direction == LEFT) ? (p.x - hit.x) : (hit.x - p.x);

            if (hit_dist < best_dist) {
                best = hit;
            } else if (hit_dist == best_dist) {
                if (ray_contact_precedes(curve, SymbolicY{p.y, p.index}, direction, hit.edge,
                                         hit.side, best.edge, best.side)) {
                    best = hit;
                }
            }
        }
    }

    if (best.hit) {
        assert(submap.start_vertex != NONE && submap.end_vertex != NONE &&
               "[C91 §2.4(iii)]: C endpoints must be set");
        const SymbolicY ray_y{p.y, p.index};
        std::size_t attributed = NONE;
        for (std::size_t k = 0; k < arcs.count; ++k) {
            if (!subarc_contains_point(subs[k], curve, best.edge, best.side, ray_y,
                                       submap.start_vertex, submap.end_vertex))
                continue;
            if (arcs.arcs[k] == best.hit_arc_idx) {
                attributed = arcs.arcs[k];
                break;
            }
            if (attributed == NONE)
                attributed = arcs.arcs[k];
        }
        best.hit_arc_idx = attributed;
    }

    if (require_hit) {
        assert(best.hit && "[C91 §3.1 tex 181]: local shoot inside a region must hit (Lemma 2.1)");
        assert(best.hit_arc_idx != NONE &&
               "[C91 §2.2 Lemma 2.1]: an in-region shot's first contact "
               "lies ON the region's boundary arcs");
    }
    if (record)
        if (auto* trace = AnimationTrace::current())
            trace->ray(curve, p, direction, best, region);
    return best;
}

struct SharedVertexCompanion {
    bool at_y_extremum = false;
    Side inside_side = LEFT;
};
static SharedVertexCompanion shared_vertex_companion(const Polygon& first_curve,
                                                     const Polygon& second_curve, bool at_end) {
    SharedVertexCompanion out;

    const Point& first_curve_neighbor =
        at_end ? first_curve.vertex(first_curve.num_vertices() - 2) : first_curve.vertex(1);
    const Point& second_curve_neighbor =
        at_end ? second_curve.vertex(1) : second_curve.vertex(second_curve.num_vertices() - 2);
    Point prev_pt = at_end ? first_curve_neighbor : second_curve_neighbor;
    Point next_pt = at_end ? second_curve_neighbor : first_curve_neighbor;
    prev_pt.x -= at_end ? first_curve.edge_horizontal_shift(first_curve.num_edges() - 1)
                        : second_curve.edge_horizontal_shift(second_curve.num_edges() - 1);
    next_pt.x +=
        at_end ? second_curve.edge_horizontal_shift(0) : first_curve.edge_horizontal_shift(0);
    const Point& v_pt =
        at_end ? first_curve.vertex(first_curve.num_vertices() - 1) : first_curve.vertex(0);
    if (!is_local_y_extremum(prev_pt, v_pt, next_pt))
        return out;
    out.at_y_extremum = true;
    const bool is_max = point_y_above(v_pt, prev_pt) && point_y_above(v_pt, next_pt);

    const bool prev_left = extremum_prev_branch_left(prev_pt, v_pt, next_pt);
    if (at_end) {
        const Side minus_x = is_max ? LEFT : RIGHT;
        const Side plus_x = is_max ? RIGHT : LEFT;
        out.inside_side = prev_left ? plus_x : minus_x;
    } else {
        const Side minus_x = is_max ? RIGHT : LEFT;
        const Side plus_x = is_max ? LEFT : RIGHT;
        out.inside_side = prev_left ? minus_x : plus_x;
    }
    return out;
}

std::size_t fusion_startup(FusionState& state, const Submap& first_submap,
                           const Polygon& first_curve, const Submap& second_submap,
                           const Polygon& second_curve, const RayShootingOracle& oracle1,
                           const RayShootingOracle& oracle2) {
    assert(state.sequence.size() >= 2 &&
           "[C91 §3.1]: fusion sequence needs at least a₀ and a_{m+1}");
    const bool at_end = state.junction_at_end;
    if (auto* trace = AnimationTrace::current()) {
        trace->checkpoint("fusion", at_end ? 0 : 1);
        trace->merge_inputs(first_submap, second_submap);
    }
    [[maybe_unused]] const std::size_t junction_edge = at_end ? first_curve.num_edges() - 1 : 0;
    [[maybe_unused]] const Side a0_side = at_end ? RIGHT : LEFT;
    [[maybe_unused]] const Side am1_side = at_end ? LEFT : RIGHT;
    const auto& a0 = state.sequence[0];
    assert(a0.is_companion && a0.side == a0_side &&
           "[C91 §3.1 tex 179]: sequence[0] = a₀ (tour-start companion)");
    assert(a0.edge == junction_edge &&
           "[C91 §3.1 tex 179]: a₀ duplicates C₁ ∩ C₂ = the walked curve's "
           "endpoint, so a₀.edge = the junction-incident boundary edge");
    [[maybe_unused]] const auto& a_m1 = state.sequence.back();
    assert(a_m1.is_companion && a_m1.side == am1_side &&
           "[C91 §3.1 tex 179]: sequence.back() = a_{m+1} (tour-end companion)");
    assert(a_m1.edge == junction_edge && "[C91 §3.1 tex 179]: a_{m+1} duplicates C₁ ∩ C₂, so "
                                         "a_{m+1}.edge = the junction-incident boundary edge");

    std::size_t a0_arc;
    {
        std::size_t turn_arc = at_end ? first_submap.end_arc : first_submap.start_arc;
        assert(turn_arc != NONE && turn_arc < first_submap.num_arcs() &&
               !first_submap.arc(turn_arc).dead &&
               "[C91 §2.4(iii) tex 138]: S₁'s endpoint arcs must be set "
               "(normal form)");
        const Arc& ta = first_submap.arc(turn_arc);
        bool ends_at_a0 =
            ta.last_edge == a0.edge && ta.last_side == a0.side &&
            symbolic_y_equal(first_submap.arc_end_symbolic_y(turn_arc, first_curve), a0.y);
        a0_arc = ends_at_a0 ? (turn_arc + 1) % first_submap.num_arcs() : turn_arc;
    }
    std::size_t a0_region_s1 = first_submap.arc(a0_arc).region_node;
    Side a0_dir = shooting_direction(a0.edge, a0.side, first_curve);

    std::size_t junction_v = at_end ? first_curve.num_vertices() - 1 : 0;
    const Point& a0_point = first_curve.vertex(junction_v);

    const SharedVertexCompanion first_curve_companion =
        shared_vertex_companion(first_curve, second_curve, at_end);
    const bool junction_inside_pair =
        first_curve_companion.at_y_extremum && a0.side == first_curve_companion.inside_side;

    RayHit hit_c1{};
    if (!junction_inside_pair)
        hit_c1 = local_shoot(a0_point, a0_dir, a0_region_s1, first_submap, first_curve, oracle1,
                             true, perturbed_x_offset(first_curve, a0.y, a0.edge));

    std::size_t target_junction_arc = at_end ? second_submap.start_arc : second_submap.end_arc;
    assert(target_junction_arc != NONE &&
           "[C91 §2.4(iii) tex 138]: S₂'s C-endpoint arc pointer must be "
           "set (normal form)");
    std::size_t a0_region_s2 = second_submap.arc(target_junction_arc).region_node;

    RayHit hit_c2{};
    if (!junction_inside_pair)

        hit_c2 = local_shoot(a0_point, a0_dir, a0_region_s2, second_submap, second_curve, oracle2,
                             true, perturbed_x_offset(first_curve, a0.y, a0.edge));

    bool c0_on_c2 = true;
    RayHit c0{};
    if (!junction_inside_pair) {
        assert(hit_c1.hit && "[C91 §3.1]: a₀ must see ∂C₁ (Lemma 2.1)");
        assert(hit_c2.hit && "[C91 §3.1]: a₀ must see ∂C₂ (Lemma 2.1)");
        {
            Exact d1 = (a0_dir == LEFT) ? (a0_point.x - hit_c1.x) : (hit_c1.x - a0_point.x);
            Exact d2 = (a0_dir == LEFT) ? (a0_point.x - hit_c2.x) : (hit_c2.x - a0_point.x);

            bool c2_at_or_before;
            if (hit_c2.wrapped != hit_c1.wrapped) {
                c2_at_or_before = !hit_c2.wrapped;
            } else if (d2 != d1) {
                c2_at_or_before = d2 < d1;
            } else {
                const Polygon Cm = at_end ? Polygon(first_curve, second_curve)
                                          : Polygon(second_curve, first_curve);
                const std::size_t off_w = at_end ? 0 : second_curve.num_edges();
                const std::size_t off_t = at_end ? first_curve.num_edges() : 0;
                c2_at_or_before =
                    !ray_contact_precedes(Cm, a0.y, a0_dir, off_w + hit_c1.edge, hit_c1.side,
                                          off_t + hit_c2.edge, hit_c2.side);
            }
            if (c2_at_or_before) {
                c0 = hit_c2;
                c0_on_c2 = true;
            } else {
                c0 = hit_c1;
                c0_on_c2 = false;
            }
        }
    }

    auto resolve_s2_region = [&](std::size_t initial_region, bool leaving_downward) -> std::size_t {
        SymbolicY a0y = a0.y;
        for (std::size_t ci : second_submap.node(initial_region).incident_chords) {
            assert(!second_submap.chord(ci).dead &&
                   "[C91 §2.4]: normal-form (compacted) submap must have no "
                   "dead chords in incident_chords");
            const auto& ch = second_submap.chord(ci);
            if (!symbolic_y_equal(ch.symbolic_y(), a0y))
                continue;

            const std::size_t tj_edge = at_end ? 0 : second_curve.num_edges() - 1;
            const Side coin_dir = junction_inside_pair ? (a0_dir == LEFT ? RIGHT : LEFT) : a0_dir;
            const bool incident_on_a0 =
                (ch.left_edge == tj_edge &&
                 shooting_direction(tj_edge, ch.left_side, second_curve) == coin_dir) ||
                (ch.right_edge == tj_edge &&
                 shooting_direction(tj_edge, ch.right_side, second_curve) == coin_dir);
            if (!incident_on_a0)
                continue;

            assert(!ch.is_null_length && "[C91 §2.1 tex 72]: chords at the junction's symbolic "
                                         "y with an endpoint at the target's junction "
                                         "companion are the companion-pair chords (non-null)");
            std::size_t below_r = NONE, above_r = NONE;
            second_submap.chord_regions_below_above(ci, second_curve, &below_r, &above_r);
            assert((below_r == ch.region[0] || below_r == ch.region[1]) &&
                   (above_r == ch.region[0] || above_r == ch.region[1]));
            return leaving_downward ? below_r : above_r;
        }
        return initial_region;
    };

    if (!junction_inside_pair) {
        const bool a0_on_first_curve = true;
        const bool c0_on_first_curve = !c0_on_c2;
        if (a0_point.x < c0.x)
            state.chords.push_back(
                {a0.y, a0.edge, a0.side, c0.edge, c0.side, a0_on_first_curve, c0_on_first_curve});
        else
            state.chords.push_back(
                {a0.y, c0.edge, c0.side, a0.edge, a0.side, c0_on_first_curve, a0_on_first_curve});
        record_discovery(state, first_curve, second_curve);
    }

    if (c0_on_c2) {
        state.p = a0_point;
        state.p_edge = a0.edge;
        state.p_side = a0.side;
        state.p_y = a0.y;

        std::size_t next_v = at_end ? junction_v - 1 : junction_v + 1;
        SymbolicY y_here = symbolic_y_of(first_curve.vertex(junction_v));
        SymbolicY y_next = symbolic_y_of(first_curve.vertex(next_v));
        bool a0_leaving_downward = symbolic_y_less(y_next, y_here);

        state.current_region = resolve_s2_region(a0_region_s2, a0_leaving_downward);
        return 1;
    } else {
        state.p = Point{c0.x, c0.y, a0.y.tag};
        state.p_edge = c0.edge;
        state.p_side = c0.side;
        state.p_y = a0.y;
        assert(state.p.index == state.p_y.tag &&
               "[C91 §2 tex 47]: p's SoS tag must match its symbolic y");

        assert(c0.edge < first_curve.num_edges());
        const auto& e = first_curve.edge(c0.edge);
        bool edge_ascending = symbolic_y_less(symbolic_y_of(first_curve.vertex(e.start_idx)),
                                              symbolic_y_of(first_curve.vertex(e.end_idx)));
        bool c0_leaving_downward = (c0.side == LEFT) ? !edge_ascending : edge_ascending;

        state.current_region = resolve_s2_region(a0_region_s2, c0_leaving_downward);

        std::size_t n_edges = first_curve.num_edges();
        auto trav_pos = [&](std::size_t edge, Side side) -> std::size_t {
            if (at_end)
                return (side == RIGHT) ? (n_edges - 1) - edge : n_edges + edge;
            return (side == LEFT) ? edge : 2 * n_edges - 1 - edge;
        };

        SymbolicY c0_y = a0.y;
        std::size_t c0_pos = trav_pos(c0.edge, c0.side);

        auto vertex_past_c0 = [&](const FusionVertex& v) -> bool {
            std::size_t v_pos = trav_pos(v.edge, v.side);
            if (v_pos != c0_pos)
                return v_pos > c0_pos;
            assert(v.edge < first_curve.num_edges());
            const auto& e = first_curve.edge(v.edge);
            bool edge_ascending = symbolic_y_less(symbolic_y_of(first_curve.vertex(e.start_idx)),
                                                  symbolic_y_of(first_curve.vertex(e.end_idx)));
            bool traversal_ascending = (v.side == LEFT) ? edge_ascending : !edge_ascending;
            return traversal_ascending ? symbolic_y_geq(v.y, c0_y) : symbolic_y_leq(v.y, c0_y);
        };

        assert(hit_c1.hit_arc_idx != NONE && "[C91 §3.1]: c0 must carry hit context");
        std::size_t c0_arc_idx = hit_c1.hit_arc_idx;

        std::size_t cw_key = 0;
        {
            const std::size_t N = first_submap.num_arcs();
            const std::size_t origin_arc = at_end ? first_submap.end_arc : first_submap.start_arc;
            if (c0_arc_idx == origin_arc) {
                const Arc& oa = first_submap.arc(origin_arc);
                bool leading;
                if (at_end) {
                    if (oa.first_side == LEFT && oa.last_side == RIGHT)
                        leading = hit_c1.side == RIGHT;
                    else if (oa.first_side == LEFT)
                        leading = hit_c1.side == RIGHT ||
                                  (hit_c1.side == LEFT && hit_c1.edge <= oa.last_edge);
                    else
                        leading = hit_c1.side == RIGHT && hit_c1.edge >= oa.last_edge;
                } else {
                    if (oa.first_side == RIGHT && oa.last_side == LEFT)
                        leading = hit_c1.side == LEFT;
                    else if (oa.first_side == LEFT)
                        leading = hit_c1.side == LEFT && hit_c1.edge <= oa.last_edge;
                    else
                        leading = hit_c1.side == LEFT ||
                                  (hit_c1.side == RIGHT && hit_c1.edge >= oa.last_edge);
                }
                cw_key = leading ? 0 : N;
            } else if (at_end) {
                std::size_t lrb = first_submap.left_right_boundary();
                std::size_t num_right = N - lrb;
                cw_key = ((c0_arc_idx >= lrb) ? (c0_arc_idx - lrb) : (num_right + c0_arc_idx)) + 1;
            } else {
                cw_key = c0_arc_idx + 1;
            }
        }

        std::size_t lo = state.arc_starts[cw_key];

        while (lo < state.sequence.size() && !vertex_past_c0(state.sequence[lo])) {
            lo++;
        }

        assert(lo < state.sequence.size() &&
               "[C91 §3.1 tex 179, 188]: skip-to-c₀ must land within sequence");

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

void fuse_submaps(FusionState& state, const Submap& first_submap, const Polygon& first_curve,
                  const Submap& second_submap, const Polygon& second_curve,
                  const RayShootingOracle& oracle1, const RayShootingOracle& oracle2) {
    build_fusion_sequence(state, first_submap, first_curve);

    state.invalidated_first_chords.assign(first_submap.num_chords(), false);
    state.invalidated_second_chords.assign(second_submap.num_chords(), false);

    std::size_t k = fusion_startup(state, first_submap, first_curve, second_submap, second_curve,
                                   oracle1, oracle2);

    const bool at_end = state.junction_at_end;
    const std::size_t junction_v = at_end ? first_curve.num_vertices() - 1 : 0;

    const SharedVertexCompanion first_curve_companion =
        shared_vertex_companion(first_curve, second_curve, state.junction_at_end);
    const std::size_t jedge = state.junction_at_end ? first_curve.num_edges() - 1 : 0;
    const SymbolicY jy = symbolic_y_of(
        first_curve.vertex(state.junction_at_end ? first_curve.num_vertices() - 1 : 0));
    auto stop_is_junction_inside = [&](const FusionVertex& v) -> bool {
        return first_curve_companion.at_y_extremum && v.side == first_curve_companion.inside_side &&
               v.edge == jedge && symbolic_y_equal(v.y, jy);
    };

    const std::size_t jedge_t = state.junction_at_end ? 0 : second_curve.num_edges() - 1;
    const SharedVertexCompanion second_curve_companion =
        shared_vertex_companion(second_curve, first_curve, !state.junction_at_end);

    const Polygon Cm = state.junction_at_end ? Polygon(first_curve, second_curve)
                                             : Polygon(second_curve, first_curve);
    const std::size_t first_curve_edge_offset =
        state.junction_at_end ? 0 : second_curve.num_edges();
    const std::size_t off_target = state.junction_at_end ? first_curve.num_edges() : 0;

    while (true) {
        if (k >= state.sequence.size())
            return;

        assert(state.current_region != NONE && state.current_region < second_submap.num_nodes() &&
               !second_submap.node(state.current_region).dead &&
               "[C91 §3.1 invariant (B)]: current S₂ region must be valid");
#ifndef NDEBUG
        {
            RayHit b =
                local_shoot(state.p, shooting_direction(state.p_edge, state.p_side, first_curve),
                            state.current_region, second_submap, second_curve, oracle2, false,
                            perturbed_x_offset(first_curve, state.p_y, state.p_edge), false);
            bool sees_region_arc = b.hit && b.hit_arc_idx != NONE;
            bool on_chord_level = false;
            for (std::size_t ci : second_submap.node(state.current_region).incident_chords) {
                const Chord& ch = second_submap.chord(ci);
                if (!ch.dead && symbolic_y_equal(ch.symbolic_y(), state.p_y))
                    on_chord_level = true;
            }

            const bool p_at_junction_inside = first_curve_companion.at_y_extremum &&
                                              state.p_side == first_curve_companion.inside_side &&
                                              state.p_edge == jedge &&
                                              symbolic_y_equal(state.p_y, jy);
            assert((sees_region_arc || on_chord_level || p_at_junction_inside) &&
                   "[C91 §3.1 invariant (B)]: p must see ∂C₂ (through "
                   "the region's arcs, along its own exit chord, or as "
                   "the junction inside duplicate)");
        }
#endif

        const std::size_t R = state.current_region;

        auto fv_point = [&](const FusionVertex& v) -> Point {
            if (v.is_companion)
                return first_curve.vertex(junction_v);

            return Point{edge_x_at_y(first_curve, v.edge, v.y), v.y.y, v.y.tag};
        };

        struct CaseIResult {
            bool fires;
            RayHit s_hit;
        };
        auto case_i_test = [&](std::size_t j) -> CaseIResult {
            const FusionVertex& aj_v = state.sequence[j];

            if (aj_v.chord_idx != NONE && first_submap.chord(aj_v.chord_idx).is_null_length)
                return {false, {}};

            Point aj_point = fv_point(aj_v);
            Side dir = shooting_direction(aj_v.edge, aj_v.side, first_curve);

            RayHit s_hit = local_shoot(aj_point, dir, R, second_submap, second_curve, oracle2,
                                       false, perturbed_x_offset(first_curve, aj_v.y, aj_v.edge));
            if (!s_hit.hit)
                return {false, {}};

            if (s_hit.hit_arc_idx == NONE)
                return {false, {}};
            const Arc& s_arc = second_submap.arc(s_hit.hit_arc_idx);
            assert(second_submap.start_vertex != NONE && second_submap.end_vertex != NONE &&
                   "[C91 §2.4(iii)]: S₂'s C endpoints must be set");
            if (!s_arc.covers(s_hit.edge, s_hit.side, second_submap.start_vertex,
                              second_submap.end_vertex))
                return {false, {}};

            Exact t_x;
            bool t_wrapped;
            std::size_t t_edge;
            Side t_side;
            if (aj_v.chord_idx != NONE) {
                const Chord& ch = first_submap.chord(aj_v.chord_idx);
                std::size_t other_edge = aj_v.is_left_endpoint ? ch.right_edge : ch.left_edge;
                t_edge = other_edge;
                t_side = aj_v.is_left_endpoint ? ch.right_side : ch.left_side;
                t_x = edge_x_at_y(first_curve, other_edge, ch.symbolic_y());

                Exact t_signed = (dir == LEFT) ? (aj_point.x - t_x) : (t_x - aj_point.x);
                if (t_signed != 0.0) {
                    t_wrapped = (t_signed < 0.0);
                } else {
                    const Exact aj_offset =
                        perturbed_x_offset(first_curve, ch.symbolic_y(), aj_v.edge);
                    const Exact t_offset = perturbed_x_offset(first_curve, ch.symbolic_y(), t_edge);
                    t_wrapped =
                        (aj_offset == t_offset) ||
                        ((dir == RIGHT) ? !(t_offset > aj_offset) : !(t_offset < aj_offset));
                }
            } else {
                assert(aj_v.is_companion);
                assert(aj_v.side == (at_end ? LEFT : RIGHT) &&
                       "[C91 §3.1]: only a_{m+1} reaches the companion branch "
                       "(a₀ is consumed by fusion_startup)");
                std::size_t s1_arc = at_end ? first_submap.end_arc : first_submap.start_arc;
                {
                    const Arc& sa = first_submap.arc(s1_arc);
                    bool starts_at_am1 =
                        sa.first_edge == aj_v.edge && sa.first_side == aj_v.side &&
                        symbolic_y_equal(first_submap.arc_start_symbolic_y(s1_arc, first_curve),
                                         aj_v.y);
                    if (starts_at_am1)
                        s1_arc = (s1_arc + first_submap.num_arcs() - 1) % first_submap.num_arcs();
                }
                std::size_t aj_region_s1 = first_submap.arc(s1_arc).region_node;
                RayHit t_hit =
                    local_shoot(aj_point, dir, aj_region_s1, first_submap, first_curve, oracle1,
                                true, perturbed_x_offset(first_curve, aj_v.y, aj_v.edge));
                t_x = t_hit.x;
                t_wrapped = t_hit.wrapped;
                t_edge = t_hit.edge;
                t_side = t_hit.side;
            }

            Exact s_dist = (dir == LEFT) ? (aj_point.x - s_hit.x) : (s_hit.x - aj_point.x);
            Exact t_dist = (dir == LEFT) ? (aj_point.x - t_x) : (t_x - aj_point.x);
            bool s_first;
            if (s_hit.wrapped != t_wrapped)
                s_first = !s_hit.wrapped;
            else if (s_dist != t_dist)
                s_first = s_dist < t_dist;
            else
                s_first = ray_contact_precedes(Cm, aj_v.y, dir, off_target + s_hit.edge, s_hit.side,
                                               first_curve_edge_offset + t_edge, t_side);
            return {s_first, s_hit};
        };

        auto s2_endpoint_point = [&](std::size_t edge, const SymbolicY& y) -> Point {
            return Point{edge_x_at_y(second_curve, edge, y), y.y, y.tag};
        };

        struct ClockwiseBoundaryPosition {
            std::size_t tp;
            Exact t;
            SymbolicY y;
            std::size_t edge;
            Side side;
        };
        auto cw_position = [&](const SymbolicY& y, std::size_t edge,
                               Side side) -> ClockwiseBoundaryPosition {
            std::size_t n_edges = first_curve.num_edges();
            std::size_t tp;
            if (at_end)
                tp = (side == RIGHT) ? (n_edges - 1) - edge : n_edges + edge;
            else
                tp = (side == LEFT) ? edge : 2 * n_edges - 1 - edge;
            Exact t = edge_t_at_y(first_curve, edge, y);
            return {tp, (side == LEFT) ? t : (1.0 - t), y, edge, side};
        };

        auto cw_less = [&](const ClockwiseBoundaryPosition& u,
                           const ClockwiseBoundaryPosition& v) -> bool {
            if (u.tp != v.tp)
                return u.tp < v.tp;
            if (u.t != v.t)
                return u.t < v.t;
            if (symbolic_y_equal(u.y, v.y))
                return false;
            assert(u.edge == v.edge && u.side == v.side && "trav_pos is injective on (edge, side)");
            const auto& e = first_curve.edge(u.edge);
            bool edge_ascending = symbolic_y_less(symbolic_y_of(first_curve.vertex(e.start_idx)),
                                                  symbolic_y_of(first_curve.vertex(e.end_idx)));
            bool trav_asc = (u.side == LEFT) ? edge_ascending : !edge_ascending;
            return trav_asc ? symbolic_y_less(u.y, v.y) : symbolic_y_greater(u.y, v.y);
        };

        auto cw_less_from_a0 = [&](const ClockwiseBoundaryPosition& u,
                                   const ClockwiseBoundaryPosition& v) -> bool {
            const ClockwiseBoundaryPosition a0_cw =
                cw_position(state.sequence[0].y, state.sequence[0].edge, state.sequence[0].side);
            const bool u_ge = !cw_less(u, a0_cw);
            const bool v_ge = !cw_less(v, a0_cw);
            if (u_ge != v_ge)
                return u_ge;
            return cw_less(u, v);
        };

        struct CaseIIResult {
            bool fires;
            RayHit p_prime_hit;
            std::size_t chord_idx;
            bool q_is_left;

            bool suppress_record = false;
        };

        struct FusionArcInterval {
            std::size_t arc;
            Subarc sub;
        };
        auto Aj_spans = [&](std::size_t j) -> std::vector<FusionArcInterval> {
            auto cw_successor = [&](std::size_t ai) -> std::size_t {
                return (ai + 1) % first_submap.num_arcs();
            };

            [[maybe_unused]] const Side a0_side = at_end ? RIGHT : LEFT;
            auto arc_after = [&](const FusionVertex& v) -> std::size_t {
                if (v.is_companion) {
                    if (v.side != a0_side)
                        return NONE;

                    std::size_t turn_arc = at_end ? first_submap.end_arc : first_submap.start_arc;
                    const Arc& ta = first_submap.arc(turn_arc);
                    bool ends_at_v =
                        ta.last_edge == v.edge && ta.last_side == v.side &&
                        symbolic_y_equal(first_submap.arc_end_symbolic_y(turn_arc, first_curve),
                                         v.y);
                    return ends_at_v ? cw_successor(turn_arc) : turn_arc;
                }
                const Chord& c = first_submap.chord(v.chord_idx);
                const Chord::AdjArcs& adj = v.is_left_endpoint ? c.left_adj : c.right_adj;

                if (adj.count == 1)
                    return cw_successor(adj.arcs[0]);

                return arc_starts_at_chord_slot(first_submap, first_curve, c, v.is_left_endpoint,
                                                adj.arcs[0])
                           ? adj.arcs[0]
                           : adj.arcs[1];
            };
            auto arc_before = [&](const FusionVertex& v) -> std::size_t {
                if (v.is_companion) {
                    if (v.side == a0_side)
                        return NONE;

                    std::size_t turn_arc = at_end ? first_submap.end_arc : first_submap.start_arc;
                    const Arc& ta = first_submap.arc(turn_arc);
                    bool starts_at_v =
                        ta.first_edge == v.edge && ta.first_side == v.side &&
                        symbolic_y_equal(first_submap.arc_start_symbolic_y(turn_arc, first_curve),
                                         v.y);
                    return starts_at_v
                               ? (turn_arc + first_submap.num_arcs() - 1) % first_submap.num_arcs()
                               : turn_arc;
                }
                const Chord& c = first_submap.chord(v.chord_idx);
                const Chord::AdjArcs& adj = v.is_left_endpoint ? c.left_adj : c.right_adj;

                if (adj.count == 1)
                    return adj.arcs[0];

                return arc_starts_at_chord_slot(first_submap, first_curve, c, v.is_left_endpoint,
                                                adj.arcs[0])
                           ? adj.arcs[1]
                           : adj.arcs[0];
            };

            {
                const FusionVertex& va = state.sequence[j - 1];
                const FusionVertex& vb = state.sequence[j];
                if (va.edge == vb.edge && va.side == vb.side && symbolic_y_equal(va.y, vb.y)) {
                    std::vector<FusionArcInterval> spans;
                    spans.push_back(
                        {arc_after(va), Subarc{va.edge, va.side, vb.edge, vb.side, va.y, vb.y}});
                    return spans;
                }
            }

            const std::size_t first = arc_after(state.sequence[j - 1]);
            const std::size_t last = arc_before(state.sequence[j]);
            assert(first != NONE && last != NONE &&
                   "[C91 §3.1 tex 199]: A_j is delimited by real stops "
                   "(a₀ opens and a_{m+1} closes the tour, so neither "
                   "companion NONE case is reachable here)");

            std::vector<FusionArcInterval> spans;
            if (first == last) {
                spans.push_back(
                    {first, Subarc{state.sequence[j - 1].edge, state.sequence[j - 1].side,
                                   state.sequence[j].edge, state.sequence[j].side,
                                   state.sequence[j - 1].y, state.sequence[j].y}});
                return spans;
            }
            [[maybe_unused]] std::size_t guard = 0;
            for (std::size_t ai = first;; ai = cw_successor(ai)) {
                ++guard;
                assert(guard <= first_submap.num_arcs() &&
                       "[C91 §2.4(iii) tex 138]: the table walk from "
                       "arc-after(a_{j-1}) must reach arc-before(a_j)");
                const Arc& a = first_submap.arc(ai);
                assert(!a.dead && "[C91 §2.4]: normal-form S₁ has no dead arcs");
                if (ai == first) {
                    spans.push_back(
                        {ai, Subarc{state.sequence[j - 1].edge, state.sequence[j - 1].side,
                                    a.last_edge, a.last_side, state.sequence[j - 1].y,
                                    first_submap.arc_end_symbolic_y(ai, first_curve)}});
                } else if (ai == last) {
                    spans.push_back({ai, Subarc{a.first_edge, a.first_side, state.sequence[j].edge,
                                                state.sequence[j].side,
                                                first_submap.arc_start_symbolic_y(ai, first_curve),
                                                state.sequence[j].y}});
                } else if (a.edge_count != 0) {
                    spans.push_back(
                        {ai, Subarc{a.first_edge, a.first_side, a.last_edge, a.last_side,
                                    first_submap.arc_start_symbolic_y(ai, first_curve),
                                    first_submap.arc_end_symbolic_y(ai, first_curve)}});
                }

                if (ai == last)
                    break;
            }
            return spans;
        };

        auto case_ii_test = [&](std::size_t j) -> CaseIIResult {
            CaseIIResult best{false, {}, NONE, false};
            ClockwiseBoundaryPosition best_cw{};
            bool best_cw_valid = false;

            auto aj_spans = Aj_spans(j);
            assert(!aj_spans.empty() && "[C91 §3.1 tex 199]: A_j spans at least one structure");

            auto p_cw = cw_position(state.p_y, state.p_edge, state.p_side);

            for (std::size_t ci : second_submap.node(R).incident_chords) {
                const Chord& chord_ab = second_submap.chord(ci);
                if (chord_ab.dead)
                    continue;

                const bool ab_wraps =
                    !chord_ab.is_null_length && chord_runs_through_infinity(second_curve, chord_ab);

                SymbolicY chord_ab_y = chord_ab.symbolic_y();
                Point a_pt = s2_endpoint_point(chord_ab.left_edge, chord_ab_y);
                Point b_pt = s2_endpoint_point(chord_ab.right_edge, chord_ab_y);

                for (bool is_left : {true, false}) {
                    std::size_t q_edge = is_left ? chord_ab.left_edge : chord_ab.right_edge;
                    Side q_side = is_left ? chord_ab.left_side : chord_ab.right_side;

                    const bool q_is_junction_inside =
                        second_curve_companion.at_y_extremum && q_edge == jedge_t &&
                        q_side == second_curve_companion.inside_side &&
                        symbolic_y_equal(chord_ab_y, jy);
                    Point q_point = is_left ? a_pt : b_pt;
                    Point other_point = is_left ? b_pt : a_pt;
                    Side shoot_dir = shooting_direction(q_edge, q_side, second_curve);

                    for (const FusionArcInterval& span : aj_spans) {
                        const std::size_t aj_arc = span.arc;
                        const Arc& aj_arc_struct = first_submap.arc(aj_arc);
                        const Subarc& aj_sub = span.sub;

                        assert_subarc_clockwise(aj_sub);

                        RayHit hit =
                            oracle1.shoot(q_point, shoot_dir, aj_arc, aj_sub,
                                          perturbed_x_offset(second_curve, chord_ab_y, q_edge));
                        if (!hit.hit)
                            continue;

                        Exact q_to_hit =
                            (shoot_dir == LEFT) ? (q_point.x - hit.x) : (hit.x - q_point.x);

                        assert((q_to_hit < 0.0 ? hit.wrapped
                                               : (q_to_hit > 0.0 ? !hit.wrapped : true)) &&
                               "[C91 §2.1 tex 70]: wrap flag consistent "
                               "with the signed travel distance");

                        Exact q_to_other = (shoot_dir == LEFT) ? (q_point.x - other_point.x)
                                                               : (other_point.x - q_point.x);
                        if (chord_ab.is_null_length)
                            continue;

                        const std::size_t oth_e =
                            is_left ? chord_ab.right_edge : chord_ab.left_edge;
                        const Side oth_s = is_left ? chord_ab.right_side : chord_ab.left_side;
                        auto beyond_other = [&]() -> bool {
                            return ray_contact_precedes(
                                Cm, chord_ab_y, shoot_dir, off_target + oth_e, oth_s,
                                first_curve_edge_offset + hit.edge, hit.side);
                        };
                        if (ab_wraps) {
                            if (hit.wrapped && (q_to_hit > q_to_other ||
                                                (q_to_hit == q_to_other && beyond_other()))) {
                                continue;
                            }
                        } else {
                            if (hit.wrapped || q_to_hit > q_to_other ||
                                (q_to_hit == q_to_other && beyond_other())) {
                                continue;
                            }
                        }

                        assert(first_submap.start_vertex != NONE &&
                               first_submap.end_vertex != NONE &&
                               "[C91 §2.4(iii)]: S₁'s C endpoints must "
                               "be set");

                        if (!subarc_contains_point(aj_sub, first_curve, hit.edge, hit.side,
                                                   chord_ab_y, first_submap.start_vertex,
                                                   first_submap.end_vertex)) {
                            continue;
                        }

                        auto hit_cw = cw_position(chord_ab_y, hit.edge, hit.side);
                        if (!cw_less_from_a0(p_cw, hit_cw)) {
                            continue;
                        }

                        Side s_back_dir = shooting_direction(hit.edge, hit.side, first_curve);

                        Point s_point{hit.x, hit.y, chord_ab_y.tag};

                        const Exact s_offset =
                            perturbed_x_offset(first_curve, chord_ab_y, hit.edge);
                        RayHit t_hit =
                            local_shoot(s_point, s_back_dir, aj_arc_struct.region_node,
                                        first_submap, first_curve, oracle1, true, s_offset);
                        Exact s_to_q = (s_back_dir == LEFT) ? (s_point.x - q_point.x)
                                                            : (q_point.x - s_point.x);
                        Exact s_to_t =
                            (s_back_dir == LEFT) ? (s_point.x - t_hit.x) : (t_hit.x - s_point.x);
                        bool suppress = false;
                        if (q_is_junction_inside) {
                            if (!(hit.edge == jedge &&
                                  hit.side == first_curve_companion.inside_side &&
                                  first_curve_companion.at_y_extremum))
                                continue;
                            suppress = true;
                        } else {
                            bool q_behind;
                            if (s_to_q != 0.0) {
                                q_behind = (s_to_q < 0.0);
                            } else {
                                const Exact q_offset =
                                    perturbed_x_offset(second_curve, chord_ab_y, q_edge);
                                q_behind = (s_offset == q_offset) ||
                                           ((s_back_dir == RIGHT) ? !(q_offset > s_offset)
                                                                  : !(q_offset < s_offset));
                            }
                            bool q_first;
                            if (q_behind != t_hit.wrapped)
                                q_first = !q_behind;
                            else if (s_to_q != s_to_t)
                                q_first = s_to_q < s_to_t;
                            else
                                q_first =
                                    !ray_contact_precedes(Cm, chord_ab_y, s_back_dir,
                                                          first_curve_edge_offset + t_hit.edge,
                                                          t_hit.side, off_target + q_edge, q_side);
                            if (!q_first)
                                continue;
                        }

                        if (!best_cw_valid || cw_less_from_a0(best_cw, hit_cw)) {
                            best = {true, hit, ci, is_left, suppress};
                            best_cw = hit_cw;
                            best_cw_valid = true;
                        }
                    }
                }
            }
            return best;
        };

        for (std::size_t j = k;; ++j) {
            if (j == state.sequence.size()) {
                return;
            }

            state.current_stop = j;
            if (auto* trace = AnimationTrace::current())
                trace->cursor(first_curve, fv_point(state.sequence[j]));

            if (stop_is_junction_inside(state.sequence[j])) {
                const FusionVertex& aj_v = state.sequence[j];
                if (aj_v.chord_idx != NONE) {
                    assert(aj_v.chord_idx < state.invalidated_first_chords.size());
                    if (!state.invalidated_first_chords[aj_v.chord_idx])
                        if (auto* trace = AnimationTrace::current())
                            trace->invalidate(first_submap, first_curve, aj_v.chord_idx);
                    state.invalidated_first_chords[aj_v.chord_idx] = true;
                }
                if (aj_v.is_companion)
                    return;
                state.p = fv_point(aj_v);
                state.p_edge = aj_v.edge;
                state.p_side = aj_v.side;
                state.p_y = aj_v.y;
                k = j + 1;
                break;
            }

            if (auto r = case_i_test(j); r.fires) {
                const FusionVertex& aj_v = state.sequence[j];
                Point aj_point = fv_point(aj_v);

                if (aj_point.x < r.s_hit.x)
                    state.chords.push_back(
                        {aj_v.y, aj_v.edge, aj_v.side, r.s_hit.edge, r.s_hit.side, true, false});
                else
                    state.chords.push_back(
                        {aj_v.y, r.s_hit.edge, r.s_hit.side, aj_v.edge, aj_v.side, false, true});

                record_discovery(state, first_curve, second_curve);

                if (aj_v.chord_idx != NONE) {
                    assert(aj_v.chord_idx < state.invalidated_first_chords.size());
                    if (!state.invalidated_first_chords[aj_v.chord_idx])
                        if (auto* trace = AnimationTrace::current())
                            trace->invalidate(first_submap, first_curve, aj_v.chord_idx);
                    state.invalidated_first_chords[aj_v.chord_idx] = true;
                }

                state.p = aj_point;
                state.p_edge = aj_v.edge;
                state.p_side = aj_v.side;
                state.p_y = aj_v.y;

                for (std::size_t ci : second_submap.node(R).incident_chords) {
                    const Chord& ch = second_submap.chord(ci);
                    if (ch.dead || ch.is_null_length)
                        continue;
                    if (!symbolic_y_equal(ch.symbolic_y(), aj_v.y))
                        continue;

                    auto endpoint_matches = [&](std::size_t ce, Side cs) {
                        if (cs != r.s_hit.side)
                            return false;
                        if (ce == r.s_hit.edge)
                            return true;
                        const std::size_t lo_e = std::min(ce, r.s_hit.edge);
                        const std::size_t hi_e = std::max(ce, r.s_hit.edge);
                        if (lo_e + 1 != hi_e)
                            return false;
                        const std::size_t sv = hi_e;
                        return symbolic_y_of(second_curve.vertex(sv)).tag == aj_v.y.tag &&
                               !second_curve.is_y_extremum(sv);
                    };
                    const bool at_endpoint = endpoint_matches(ch.left_edge, ch.left_side) ||
                                             endpoint_matches(ch.right_edge, ch.right_side);
                    if (!at_endpoint)
                        continue;
                    std::size_t below_r = NONE, above_r = NONE;
                    second_submap.chord_regions_below_above(ci, second_curve, &below_r, &above_r);

                    state.current_region =
                        leaving_clockwise_downward(first_curve, aj_v.edge, aj_v.side, aj_v.y)
                            ? below_r
                            : above_r;
                    break;
                }
                k = j + 1;
                break;
            }

            if (auto r = case_ii_test(j); r.fires) {
                const Chord& chord_ab = second_submap.chord(r.chord_idx);
                std::size_t q_edge = r.q_is_left ? chord_ab.left_edge : chord_ab.right_edge;
                Side q_side = r.q_is_left ? chord_ab.left_side : chord_ab.right_side;
                SymbolicY chord_y = chord_ab.symbolic_y();
                Point q_point = s2_endpoint_point(q_edge, chord_y);

                Point p_prime{r.p_prime_hit.x, r.p_prime_hit.y, chord_y.tag};

                if (!r.suppress_record) {
                    if (q_point.x < p_prime.x)
                        state.chords.push_back({chord_y, q_edge, q_side, r.p_prime_hit.edge,
                                                r.p_prime_hit.side, false, true});
                    else
                        state.chords.push_back({chord_y, r.p_prime_hit.edge, r.p_prime_hit.side,
                                                q_edge, q_side, true, false});
                    record_discovery(state, first_curve, second_curve);
                }

                assert(r.chord_idx < state.invalidated_second_chords.size());
                if (!state.invalidated_second_chords[r.chord_idx])
                    if (auto* trace = AnimationTrace::current())
                        trace->invalidate(second_submap, second_curve, r.chord_idx);
                state.invalidated_second_chords[r.chord_idx] = true;

                state.p = p_prime;
                state.p_edge = r.p_prime_hit.edge;
                state.p_side = r.p_prime_hit.side;

                state.p_y = chord_y;

                assert(!chord_ab.is_null_length &&
                       "[C91 §3.1 tex 222]: null-length chords cannot fire "
                       "case (ii) — no hit lies on a zero-length ab");

                {
                    const auto& pe = first_curve.edge(r.p_prime_hit.edge);
                    const bool e_asc =
                        symbolic_y_less(symbolic_y_of(first_curve.vertex(pe.start_idx)),
                                        symbolic_y_of(first_curve.vertex(pe.end_idx)));
                    const bool leaving_up = (r.p_prime_hit.side == LEFT) == e_asc;
                    std::size_t below = NONE, above = NONE;
                    second_submap.chord_regions_below_above(r.chord_idx, second_curve, &below,
                                                            &above);
                    state.current_region = leaving_up ? above : below;
                    assert((state.current_region == chord_ab.region[0] ||
                            state.current_region == chord_ab.region[1]) &&
                           "[C91 §3.1 tex 206]: the entered region is "
                           "one of the exit chord's two sides");
                }

                {
                    const FusionVertex& aj_end = state.sequence[j];
                    bool p_at_aj = state.p_edge == aj_end.edge && state.p_side == aj_end.side &&
                                   symbolic_y_equal(state.p_y, aj_end.y);
                    k = p_at_aj ? j + 1 : j;
                }
                break;
            }
        }
    }
}

void build_fusion_sequence(FusionState& state, const Submap& submap, const Polygon& curve) {
    std::size_t n = curve.num_vertices();
    assert(n >= 2);
    const bool at_end = state.junction_at_end;

    std::size_t junction_v = at_end ? n - 1 : 0;
    std::size_t junction_edge = at_end ? n - 2 : 0;

    SymbolicY junction_y = symbolic_y_of(curve.vertex(junction_v));

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

    std::size_t num_arcs = submap.num_arcs();
    std::size_t lrb = submap.left_right_boundary();
    std::size_t num_right = num_arcs - lrb;

    const std::size_t origin_arc = at_end ? submap.end_arc : submap.start_arc;
    assert(origin_arc != NONE && origin_arc < num_arcs &&
           "[C91 §2.4(iii) tex 138]: the walked submap's endpoint arcs "
           "must be set (normal form)");

    auto origin_leading = [&](std::size_t edge, Side side) -> bool {
        const Arc& a = submap.arc(origin_arc);
        if (at_end) {
            if (a.first_side == LEFT && a.last_side == RIGHT)
                return side == RIGHT;
            assert(a.first_side == a.last_side && a.wraps() &&
                   "[C91 §2.4 tex 142]: at_end origin arc must wrap C's "
                   "end vertex");
            if (a.first_side == LEFT)
                return side == RIGHT || (side == LEFT && edge <= a.last_edge);
            return side == RIGHT && edge >= a.last_edge;
        }

        if (a.first_side == RIGHT && a.last_side == LEFT)
            return side == LEFT;
        assert(a.first_side == a.last_side && a.wraps() &&
               "[C91 §2.4 tex 142]: at_start origin arc must wrap C's "
               "start vertex");
        if (a.first_side == LEFT)
            return side == LEFT && edge <= a.last_edge;
        return side == LEFT || (side == RIGHT && edge >= a.last_edge);
    };

    auto cw_pos = [&](std::size_t arc_idx, std::size_t edge, Side side) -> std::size_t {
        if (arc_idx == origin_arc)
            return origin_leading(edge, side) ? 0 : num_arcs;
        std::size_t base =
            at_end ? ((arc_idx >= lrb) ? arc_idx - lrb : num_right + arc_idx) : arc_idx;
        return base + 1;
    };

    struct KeyedVertex {
        FusionVertex v;
        std::size_t key;
    };
    std::vector<KeyedVertex> endpoints;
    endpoints.reserve(2 * submap.num_live_chords());

    auto make_vertex = [](const Chord& c, std::size_t ci, bool is_left) -> FusionVertex {
        FusionVertex fv;
        fv.y = c.symbolic_y();
        fv.edge = is_left ? c.left_edge : c.right_edge;
        fv.side = is_left ? c.left_side : c.right_side;
        fv.chord_idx = ci;
        fv.is_left_endpoint = is_left;
        fv.is_companion = false;
        return fv;
    };

    auto starting_arc = [&](const Chord& c, bool left_slot) -> std::size_t {
        const Chord::AdjArcs& adj = left_slot ? c.left_adj : c.right_adj;
        assert(adj.count == 2);
        bool s0 = arc_starts_at_chord_slot(submap, curve, c, left_slot, adj.arcs[0]);
        assert(s0 != arc_starts_at_chord_slot(submap, curve, c, left_slot, adj.arcs[1]) &&
               "[C91 §2.4]: exactly one adj arc starts at a mid-edge "
               "chord endpoint");
        return s0 ? adj.arcs[0] : adj.arcs[1];
    };

    for (std::size_t ci = 0; ci < submap.num_chords(); ++ci) {
        assert(!submap.chord(ci).dead && "[C91 §2.4]: normal-form submap has no dead chords");
        const auto& c = submap.chord(ci);

        if (c.is_null_length)
            continue;

        std::size_t left_arc = (c.left_adj.count == 2) ? starting_arc(c, true) : c.left_adj.arcs[0];
        endpoints.push_back({make_vertex(c, ci, true), cw_pos(left_arc, c.left_edge, c.left_side)});

        std::size_t right_arc =
            (c.right_adj.count == 2) ? starting_arc(c, false) : c.right_adj.arcs[0];
        endpoints.push_back(
            {make_vertex(c, ci, false), cw_pos(right_arc, c.right_edge, c.right_side)});
    }

    const std::size_t num_keys = num_arcs + 1;
    std::vector<KeyedVertex> sorted(endpoints.size());
    std::vector<std::size_t> bucket_starts;
    if (!endpoints.empty()) {
        std::vector<std::size_t> counts(num_keys + 1, 0);
        for (const auto& ep : endpoints)
            ++counts[ep.key + 1];
        for (std::size_t i = 1; i <= num_keys; ++i)
            counts[i] += counts[i - 1];

        bucket_starts.assign(counts.begin(),
                             counts.begin() + static_cast<std::ptrdiff_t>(num_keys));

        for (const auto& ep : endpoints)
            sorted[counts[ep.key]++] = ep;
    } else {
        bucket_starts.assign(num_keys, 0);
    }

    std::size_t n_edges = curve.num_edges();
    auto trav_pos = [&](std::size_t edge, Side side) -> std::size_t {
        if (at_end)
            return (side == RIGHT) ? junction_edge - edge : n_edges + edge;
        return (side == LEFT) ? edge : 2 * n_edges - 1 - edge;
    };
    auto vertex_before = [&](const FusionVertex& u, const FusionVertex& v) -> bool {
        std::size_t u_pos = trav_pos(u.edge, u.side);
        std::size_t v_pos = trav_pos(v.edge, v.side);
        if (u_pos != v_pos)
            return u_pos < v_pos;
        assert(u.edge < curve.num_edges());
        const auto& e = curve.edge(u.edge);
        bool edge_ascending = symbolic_y_less(symbolic_y_of(curve.vertex(e.start_idx)),
                                              symbolic_y_of(curve.vertex(e.end_idx)));
        bool trav_asc = (u.side == LEFT) ? edge_ascending : !edge_ascending;
        return trav_asc ? symbolic_y_less(u.y, v.y) : symbolic_y_greater(u.y, v.y);
    };

    {
        std::size_t i = 0;
        while (i < sorted.size()) {
            std::size_t j = i + 1;
            while (j < sorted.size() && sorted[j].key == sorted[i].key)
                ++j;

            assert(j - i <= 8 && "[C91 §2.3 tex 114 / §2.4 tex 137]: bounded arc bucket");

            if (j - i > 1) {
                for (std::size_t a = i + 1; a < j; ++a) {
                    KeyedVertex tmp = sorted[a];
                    std::size_t b = a;
                    while (b > i) {
                        if (!vertex_before(tmp.v, sorted[b - 1].v))
                            break;
                        sorted[b] = sorted[b - 1];
                        --b;
                    }
                    sorted[b] = tmp;
                }
            }
            i = j;
        }
    }

    std::vector<FusionVertex> result;
    result.reserve(sorted.size() + 2);
    result.push_back(a_0);
    for (const auto& kv : sorted)
        result.push_back(kv.v);
    result.push_back(a_m1);

    state.sequence = std::move(result);

    state.arc_starts.resize(num_keys);
    for (std::size_t i = 0; i < num_keys; ++i)
        state.arc_starts[i] = bucket_starts[i] + 1;

    {
        std::size_t exit_chord_count = 0;
        for (std::size_t ci = 0; ci < submap.num_chords(); ++ci)
            if (!submap.chord(ci).dead && !submap.chord(ci).is_null_length)
                ++exit_chord_count;
        assert(state.sequence.size() == 2 * exit_chord_count + 2 &&
               "[C91 §3.1 tex 209]: |sequence| = 2·(#exit chords) + 2");
    }
}

}
