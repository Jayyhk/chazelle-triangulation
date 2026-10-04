#include "conformality.h"
#include "../submap/shielding.h"
#include "../visibility/naive_visibility.h"
#include "fusion.h"

#include <algorithm>

namespace chazelle {

namespace {

struct ClockwiseBoundaryPosition {
    std::size_t clockwise_edge_index = 0;
    SymbolicY y{};
    bool y_ascending = true;
};

ClockwiseBoundaryPosition boundary_position(const Polygon& curve, std::size_t edge, Side side,
                                            const SymbolicY& y) {
    assert(edge < curve.num_edges());
    ClockwiseBoundaryPosition t;
    t.clockwise_edge_index = (side == LEFT) ? edge : (2 * curve.num_edges() - 1 - edge);
    const auto& e = curve.edge(edge);
    bool asc = symbolic_y_less(symbolic_y_of(curve.vertex(e.start_idx)),
                               symbolic_y_of(curve.vertex(e.end_idx)));
    t.y_ascending = (side == LEFT) ? asc : !asc;
    t.y = y;
    return t;
}

int compare_boundary_positions(const ClockwiseBoundaryPosition& a,
                               const ClockwiseBoundaryPosition& b) {
    if (a.clockwise_edge_index != b.clockwise_edge_index)
        return a.clockwise_edge_index < b.clockwise_edge_index ? -1 : 1;
    int c = symbolic_y_compare(a.y, b.y);
    if (c == 0)
        return 0;
    return a.y_ascending ? c : -c;
}

bool precedes_from_clockwise_origin(const ClockwiseBoundaryPosition& base,
                                    const ClockwiseBoundaryPosition& x,
                                    const ClockwiseBoundaryPosition& y) {
    bool first_wraps_origin = compare_boundary_positions(x, base) < 0;
    bool second_wraps_origin = compare_boundary_positions(y, base) < 0;
    if (first_wraps_origin != second_wraps_origin)
        return second_wraps_origin;
    return compare_boundary_positions(x, y) <= 0;
}

bool in_closed_clockwise_interval(const ClockwiseBoundaryPosition& x,
                                  const ClockwiseBoundaryPosition& lo,
                                  const ClockwiseBoundaryPosition& hi) {
    return precedes_from_clockwise_origin(lo, x, hi);
}

bool in_half_open_clockwise_interval(const ClockwiseBoundaryPosition& x,
                                     const ClockwiseBoundaryPosition& lo,
                                     const ClockwiseBoundaryPosition& hi) {
    return precedes_from_clockwise_origin(lo, x, hi) && compare_boundary_positions(x, hi) != 0;
}

ClockwiseBoundaryPosition fused_arc_start(const Submap& submap, const Polygon& curve,
                                          std::size_t arc_idx) {
    const Arc& a = submap.arc(arc_idx);
    return boundary_position(curve, a.first_edge, a.first_side,
                             submap.arc_start_symbolic_y(arc_idx, curve));
}

ClockwiseBoundaryPosition fused_arc_end(const Submap& submap, const Polygon& curve,
                                        std::size_t arc_idx) {
    const Arc& a = submap.arc(arc_idx);
    return boundary_position(curve, a.last_edge, a.last_side,
                             submap.arc_end_symbolic_y(arc_idx, curve));
}

bool fused_arc_contains(const Submap& submap, const Polygon& curve, std::size_t arc_idx,
                        std::size_t edge, Side side, const SymbolicY& y) {
    const Arc& a = submap.arc(arc_idx);
    assert(submap.start_vertex != NONE && submap.end_vertex != NONE &&
           "[C91 §2.4(iii)]: S's C endpoints must be set");
    if (!a.covers(edge, side, submap.start_vertex, submap.end_vertex))
        return false;
    ClockwiseBoundaryPosition x = boundary_position(curve, edge, side, y);
    ClockwiseBoundaryPosition s = fused_arc_start(submap, curve, arc_idx);
    ClockwiseBoundaryPosition e = fused_arc_end(submap, curve, arc_idx);
    if (compare_boundary_positions(s, e) <= 0)
        return compare_boundary_positions(s, x) <= 0 && compare_boundary_positions(x, e) <= 0;

    return compare_boundary_positions(s, x) <= 0 || compare_boundary_positions(x, e) <= 0;
}

}

std::vector<ArcSource> identify_arc_sources(const Submap& submap, const Polygon& curve,
                                            const Submap& first_submap, const Polygon& first_curve,
                                            const Submap& second_submap,
                                            const Polygon& second_curve) {
    const std::size_t n1e = first_curve.num_edges();
    assert(curve.num_edges() == first_curve.num_edges() + second_curve.num_edges() &&
           "[C91 §3 tex 160]: C = C₁ ∪ C₂ shares one vertex");

    std::vector<ArcSource> out(submap.num_arcs());
    for (std::size_t ai = 0; ai < submap.num_arcs(); ++ai) {
        const Arc& a = submap.arc(ai);
        assert(!a.dead && "[C91 §3.1 tex 226]: fused submap is freshly built (no dead arcs)");

        assert(!(a.first_side == a.last_side && a.wraps()) &&
               "[C91 §3.2 tex 244]: no fused arc double-wraps (junction "
               "chords break both sides)");

        const bool on_first_curve = a.first_edge < n1e;
        assert((a.last_edge < n1e) == on_first_curve &&
               "[C91 §3.2 tex 244]: fused arc must lie on one ∂Cᵢ only");

        const Submap& input_submap = on_first_curve ? first_submap : second_submap;
        const Polygon& input_curve = on_first_curve ? first_curve : second_curve;
        const std::size_t off = on_first_curve ? 0 : n1e;

        const std::size_t edge_i = a.first_edge - off;
        const SymbolicY sy = submap.arc_start_symbolic_y(ai, curve);

        auto cands = input_submap.double_identify(edge_i, sy, input_curve);
        assert(cands.count >= 1 && "[C91 §2.4 tex 144]: ∂Cᵢ point must lie on some Sᵢ arc");

        assert(input_submap.start_vertex != NONE && input_submap.end_vertex != NONE &&
               "[C91 §2.4(iii)]: Sᵢ's C endpoints must be set");
        std::size_t chosen = NONE;
        std::size_t covering_ender = NONE;
        for (std::size_t k = 0; k < cands.count; ++k) {
            std::size_t si_arc = cands.arcs[k];
            const Arc& sa = input_submap.arc(si_arc);

            if (!fused_arc_contains(input_submap, input_curve, si_arc, edge_i, a.first_side, sy))
                continue;
            bool ends_here =
                sa.last_edge == edge_i && sa.last_side == a.first_side &&
                symbolic_y_equal(input_submap.arc_end_symbolic_y(si_arc, input_curve), sy);
            if (ends_here) {
                covering_ender = si_arc;
                continue;
            }
            assert(chosen == NONE && "[C91 §2.4(iii)]: arcs tile ∂Cᵢ — exactly one covers "
                                     "clockwise-forward from any point");
            chosen = si_arc;
        }

        if (chosen == NONE && a.edge_count == 0)
            chosen = covering_ender;
        assert(chosen != NONE && "[C91 §3.2 tex 244]: every fused arc lies inside an Sᵢ arc");
        out[ai] = ArcSource{on_first_curve, chosen};
    }
    return out;
}

FusedRegionCycle fused_region_cycle(const Submap& submap, const Polygon& curve,
                                    [[maybe_unused]] std::size_t region,
                                    const std::vector<std::size_t>& arcs_of_region) {
    assert(region < submap.num_nodes() && !submap.node(region).dead);

    struct Keyed {
        std::size_t arc;
        ClockwiseBoundaryPosition pos;
    };
    std::vector<Keyed> sorted;
    sorted.reserve(arcs_of_region.size());
    for (std::size_t ai : arcs_of_region) {
        assert(ai < submap.num_arcs() && !submap.arc(ai).dead &&
               submap.arc(ai).region_node == region &&
               "[C91 §3.2]: inventory entries must be live arcs of the region");
        sorted.push_back({ai, fused_arc_start(submap, curve, ai)});
    }
    assert(!sorted.empty() && "[C91 §2.2 tex 96]: region has ≥ 1 arc");
    std::sort(sorted.begin(), sorted.end(), [](const Keyed& a, const Keyed& b) {
        return compare_boundary_positions(a.pos, b.pos) < 0;
    });
#ifndef NDEBUG
    for (std::size_t i = 1; i < sorted.size(); ++i)
        assert(compare_boundary_positions(sorted[i - 1].pos, sorted[i].pos) != 0 &&
               "[C91 §2.2 tex 96]: region arcs must start at distinct "
               "∂C positions");
#endif

    FusedRegionCycle cycle;
    for (const Keyed& k : sorted) {
        assert(cycle.count < FusedRegionCycle::MAX_ARCS &&
               "[C91 §3.2 tex 238]: fused region arc count is bounded "
               "(≤ 2 runs of constant length)");

        cycle.arcs[cycle.count++] = CycleArc{k.arc, submap.arc(k.arc).edge_count == 0};
    }
    return cycle;
}

RayHit local_shoot_fused(Point p, const SymbolicY& p_y, Side direction,
                         const FusedRegionCycle& cycle, const FusedShootContext& ctx,
                         bool require_hit, std::size_t source_edge_c) {
    const SourceOffset source_x_offset =
        (source_edge_c == NONE) ? SourceOffset{}
                                : SourceOffset{perturbed_x_offset(*ctx.curve, p_y, source_edge_c)};
    assert(ctx.submap && ctx.curve && ctx.first_curve && ctx.first_ray_shooter && ctx.arc_sources);
    const Submap& submap = *ctx.submap;
    const Polygon& curve = *ctx.curve;
    const std::size_t n1e = ctx.first_curve->num_edges();

    p.index = p_y.tag;

    RayHit best;
    best.hit = false;

    for (std::size_t li = 0; li < cycle.count; ++li) {
        {
            const std::size_t ai = cycle.arcs[li].arc;
            const Arc& a = submap.arc(ai);
            const ArcSource& pr = (*ctx.arc_sources)[ai];
            const std::size_t off = pr.on_first_curve ? 0 : n1e;
            const RayShootingOracle& oracle =
                pr.on_first_curve ? *ctx.first_ray_shooter : *ctx.second_ray_shooter;

            Subarc target;
            target.first_edge = a.first_edge - off;
            target.first_side = a.first_side;
            target.last_edge = a.last_edge - off;
            target.last_side = a.last_side;

            target.first_y = submap.arc_start_symbolic_y(ai, curve);
            target.last_y = submap.arc_end_symbolic_y(ai, curve);
            assert_subarc_clockwise(target);

            RayHit hit = oracle.shoot(p, direction, pr.input_arc, target, source_x_offset);
            if (!hit.hit)
                continue;

            Exact d = (direction == LEFT) ? (p.x - hit.x) : (hit.x - p.x);
            if (hit.wrapped)
                assert(d <= 0.0 && "[C91 §2.1 tex 70]: a wrapped hit lies at or behind "
                                   "the source in the travel direction (d == 0: a "
                                   "raw-coincident wall not strictly forward, met "
                                   "after a full wrap)");
            else
                assert((d > 0.0 || (d == 0.0 && source_edge_c != NONE)) &&
                       "[C91 §3.0(i) tex 169/§2 tex 47]: a direct hit lies "
                       "strictly forward — at raw distance 0 only by the "
                       "source's own perturbed x-offset order "
                       "(perturbed_hit_forward)");

            hit.edge += off;

            hit.hit_arc_idx = ai;

            if (!best.hit) {
                best = hit;
                continue;
            }

            Exact bd = (direction == LEFT) ? (p.x - best.x) : (best.x - p.x);
            if (hit.wrapped != best.wrapped) {
                if (!hit.wrapped)
                    best = hit;
            } else if (d < bd) {
                best = hit;
            } else if (d == bd) {
                [[maybe_unused]] auto opposes_ray = [&](const RayHit& h) {
                    assert(h.edge < curve.num_edges());
                    const auto& e = curve.edge(h.edge);
                    bool asc = symbolic_y_less(symbolic_y_of(curve.vertex(e.start_idx)),
                                               symbolic_y_of(curve.vertex(e.end_idx)));
                    Side minus_x_face = asc ? LEFT : RIGHT;
                    Side plus_x_face = asc ? RIGHT : LEFT;
                    return h.side == ((direction == RIGHT) ? minus_x_face : plus_x_face);
                };
                assert(opposes_ray(hit) && opposes_ray(best) &&
                       "[C91 §3.0(i) tex 169]: a reported hit is a wall "
                       "opposing the ray's travel");

                if (ray_contact_precedes(curve, p_y, direction, hit.edge, hit.side, best.edge,
                                         best.side))
                    best = hit;
            }
        }
    }

    if (best.hit) {
        std::size_t attributed = NONE;
        for (std::size_t li = 0; li < cycle.count; ++li) {
            const std::size_t ai = cycle.arcs[li].arc;
            if (!fused_arc_contains(submap, curve, ai, best.edge, best.side, p_y))
                continue;
            if (ai == best.hit_arc_idx) {
                attributed = ai;
                break;
            }
            if (attributed == NONE)
                attributed = ai;
        }
        best.hit_arc_idx = attributed;
    }

    if (require_hit) {
        assert(best.hit && "[C91 §3.2 tex 244]: local shoot within a fused region must hit");
        assert(best.hit_arc_idx != NONE &&
               "[C91 §2.2 Lemma 2.1]: an in-region shot's first contact "
               "lies ON the region's boundary arcs");
    }
    return best;
}

namespace {

bool duplicates_region_chord(const Submap& submap, std::size_t region, const SymbolicY& y,
                             std::size_t e1, Side s1, std::size_t e2, Side s2) {
    for (std::size_t ci : submap.node(region).incident_chords) {
        const Chord& c = submap.chord(ci);
        if (c.dead)
            continue;
        if (!symbolic_y_equal(c.symbolic_y(), y))
            continue;
        bool fwd =
            c.left_edge == e1 && c.left_side == s1 && c.right_edge == e2 && c.right_side == s2;
        bool rev =
            c.left_edge == e2 && c.left_side == s2 && c.right_edge == e1 && c.right_side == s1;
        if (fwd || rev)
            return true;
    }
    return false;
}

struct CandidateShot {
    bool success = false;
    RayHit hit{};
};

struct ExistingChord {
    bool exists = false;
    std::size_t other_edge = NONE;
    Side other_side = LEFT;
};

ExistingChord existing_region_chord_at(const Submap& submap, std::size_t region, std::size_t edge,
                                       Side side, const SymbolicY& y) {
    for (std::size_t ci : submap.node(region).incident_chords) {
        const Chord& c = submap.chord(ci);
        if (c.dead)
            continue;
        if (!symbolic_y_equal(c.symbolic_y(), y))
            continue;
        if (c.left_edge == edge && c.left_side == side)
            return {true, c.right_edge, c.right_side};
        if (c.right_edge == edge && c.right_side == side)
            return {true, c.left_edge, c.left_side};
    }
    return {};
}

CandidateShot try_candidate_vertex(std::size_t edge_c, Side side, const SymbolicY& y,
                                   const Submap& submap, const Polygon& curve, std::size_t region,
                                   std::size_t A2, const FusedRegionCycle& cycle,
                                   const FusedShootContext& fctx) {
    {
        const auto& ed = curve.edge(edge_c);
        for (std::size_t vi : {ed.start_idx, ed.end_idx}) {
            if (symbolic_y_of(curve.vertex(vi)).tag != y.tag)
                continue;
            if (is_inside_companion(curve, edge_c, side, vi))
                return CandidateShot{};
            break;
        }
    }
    if (existing_region_chord_at(submap, region, edge_c, side, y).exists)
        return CandidateShot{};

    Point p{edge_x_at_y(curve, edge_c, y), y.y, y.tag};
    Side dir = shooting_direction(edge_c, side, curve);
    RayHit hit = local_shoot_fused(p, y, dir, cycle, fctx, true, edge_c);

    CandidateShot out;
    out.hit = hit;
    out.success = hit.hit_arc_idx == A2 &&
                  !duplicates_region_chord(submap, region, y, edge_c, side, hit.edge, hit.side);
    return out;
}

struct AlphaPoint {
    std::size_t edge_a = NONE;
    Side side = LEFT;
    std::size_t edge_c = NONE;
    Exact x = 0.0;
    ClockwiseBoundaryPosition pos_a{};
};

struct PieceSearchContext {
    const Submap* submap;
    const Polygon* curve;
    std::size_t region;
    std::size_t A1;
    std::size_t A2;
    const FusedRegionCycle* cycle;
    const FusedShootContext* fctx;

    const ArcPiece* piece;
    const Polygon* input_curve;
    std::size_t edge_off;
    std::size_t lo;
    Side s;

    ClockwiseBoundaryPosition a2_start{};
    ClockwiseBoundaryPosition a2_end{};

    ClockwiseBoundaryPosition c_pos_a{}, d_pos_a{};
    ClockwiseBoundaryPosition c_pos_c{}, d_pos_c{};
};

VisiblePoint make_success(const PieceSearchContext& ctx, std::size_t p_edge_c, Side p_side,
                          const SymbolicY& y, const RayHit& hit) {
    VisiblePoint vp;
    vp.found = true;
    assert(fused_arc_contains(*ctx.submap, *ctx.curve, ctx.A1, p_edge_c, p_side, y) &&
           "[C91 §3.2]: the candidate vertex must lie on A₁");
    vp.p_table_arc = ctx.A1;
    vp.p_edge = p_edge_c;
    vp.p_side = p_side;
    vp.p_x = edge_x_at_y(*ctx.curve, p_edge_c, y);
    vp.y = y;
    vp.q_table_arc = hit.hit_arc_idx;
    vp.q_edge = hit.edge;
    vp.q_side = hit.side;
    vp.q_x = hit.x;
    return vp;
}

VisiblePoint descend_step(const PieceSearchContext& ctx, const TreeDecomposition& td,
                          std::size_t node_idx, std::size_t* next_node) {
    const Submap& Sa = *ctx.piece->submap;
    const Polygon& Ca = *ctx.piece->curve;
    const TreeDecompositionNode& node = td.node(node_idx);
    const Chord& ab = Sa.chord(node.chord_idx);

    if (ab.is_null_length) {
        std::size_t lens = NONE;
        for (std::size_t ri = 0; ri < 2 && lens == NONE; ++ri) {
            const std::size_t r = ab.region[ri];
            RegionArcs ra = collect_region_arcs(Sa, r);
            bool empty = true;
            for (std::size_t k = 0; k < ra.count && empty; ++k)
                empty = Sa.arc(ra.arcs[k]).edge_count == 0;
            if (empty)
                lens = r;
        }
        assert(lens != NONE && "[C91 §2.3 tex 105]: a null-length chord bounds the empty "
                               "region hugging the apex turn on its corner side");

        *next_node = (lens == ab.region[0]) ? node.right_child : node.left_child;
        assert(*next_node != NONE && "[C91 §2.3]: internal TD node has two children");
        return VisiblePoint{};
    }

    const SymbolicY y_ab = ab.symbolic_y();

    auto alpha_point = [&](bool left) -> AlphaPoint {
        AlphaPoint ap;
        ap.edge_a = left ? ab.left_edge : ab.right_edge;
        ap.side = left ? ab.left_side : ab.right_side;
        ap.edge_c = ctx.edge_off + ctx.lo + ap.edge_a;
        ap.x = edge_x_at_y(Ca, ap.edge_a, y_ab);
        ap.pos_a = boundary_position(Ca, ap.edge_a, ap.side, y_ab);
        return ap;
    };
    AlphaPoint el = alpha_point(true);
    AlphaPoint er = alpha_point(false);
    assert(compare_boundary_positions(el.pos_a, er.pos_a) != 0 &&
           "[C91 §2.1 tex 70]: non-null chord endpoints are distinct ∂ᾱ points");

    const bool ab_wraps = chord_runs_through_infinity(Ca, ab);

    bool left_is_a = precedes_from_clockwise_origin(ctx.c_pos_a, el.pos_a, er.pos_a);
    const AlphaPoint& A = left_is_a ? el : er;
    const AlphaPoint& B = left_is_a ? er : el;
    bool a_is_left_endpoint = left_is_a;

    bool a_in = (A.side == ctx.s);
    bool b_in = (B.side == ctx.s);

    assert(!(b_in && !a_in) && "[C91 §2.5]: b ∈ α requires a ∈ α (circle order)");

    VisiblePoint fail;
    bool a_prime_exists = false, b_prime_exists = false;
    ClockwiseBoundaryPosition a_prime_pos{}, b_prime_pos{};
    std::size_t ap_edge = NONE, bp_edge = NONE;
    Side ap_side = LEFT, bp_side = LEFT;
    Exact ap_x = 0.0, bp_x = 0.0;

    auto endpoint_shot = [&](const AlphaPoint& from, const AlphaPoint& other, bool other_in,
                             bool* prime_exists, ClockwiseBoundaryPosition* prime_pos,
                             std::size_t* pe, Side* ps, Exact* px, VisiblePoint* success) -> bool {
        Side dir = shooting_direction(from.edge_c, from.side, *ctx.curve);

        if (!ab_wraps)
            assert(((dir == RIGHT) == (other.x > from.x)) &&
                   "[C91 §2.1 tex 70]: a direct chord's shooting "
                   "direction points toward the other endpoint");
        else if (other.x != from.x)
            assert(((dir == RIGHT) == (other.x < from.x)) &&
                   "[C91 §2.1 tex 70]: a wrapping chord's shooting "
                   "direction points away from the other endpoint");

        {
            const auto& ed = ctx.curve->edge(from.edge_c);
            for (std::size_t vi : {ed.start_idx, ed.end_idx}) {
                if (symbolic_y_of(ctx.curve->vertex(vi)).tag != y_ab.tag)
                    continue;
                if (!is_inside_companion(*ctx.curve, from.edge_c, from.side, vi))
                    break;
                const std::size_t sib_edge = (from.edge_c == vi) ? vi - 1 : vi;
                const Side sib_side =
                    is_inside_companion(*ctx.curve, sib_edge, LEFT, vi) ? LEFT : RIGHT;
                assert(is_inside_companion(*ctx.curve, sib_edge, sib_side, vi) &&
                       "[C91 §2.1 tex 72]: exactly one face of the other "
                       "edge at an interior extremum is the side facing the inside companion");
                if (fused_arc_contains(*ctx.submap, *ctx.curve, ctx.A2, sib_edge, sib_side, y_ab)) {
                    RayHit h;
                    h.hit = true;
                    h.x = from.x;
                    h.y = y_ab.y;
                    h.edge = sib_edge;
                    h.side = sib_side;
                    h.wrapped = false;
                    h.hit_arc_idx = ctx.A2;
                    *success = make_success(ctx, from.edge_c, from.side, y_ab, h);
                    return true;
                }
                *prime_exists = true;
                *prime_pos = boundary_position(*ctx.curve, sib_edge, sib_side, y_ab);
                *pe = sib_edge;
                *ps = sib_side;
                *px = from.x;
                return false;
            }
        }

        {
            ExistingChord rv =
                existing_region_chord_at(*ctx.submap, ctx.region, from.edge_c, from.side, y_ab);
            if (rv.exists) {
                *prime_exists = true;
                *prime_pos = boundary_position(*ctx.curve, rv.other_edge, rv.other_side, y_ab);
                *pe = rv.other_edge;
                *ps = rv.other_side;
                *px = edge_x_at_y(*ctx.curve, rv.other_edge, y_ab);
                return false;
            }
        }
        CandidateShot shot =
            try_candidate_vertex(from.edge_c, from.side, y_ab, *ctx.submap, *ctx.curve, ctx.region,
                                 ctx.A2, *ctx.cycle, *ctx.fctx);

        Exact d_hit = (dir == LEFT) ? (from.x - shot.hit.x) : (shot.hit.x - from.x);
        Exact d_other = (dir == LEFT) ? (from.x - other.x) : (other.x - from.x);
        [[maybe_unused]] const bool hit_at_other = shot.hit.wrapped == ab_wraps && d_hit == d_other;
        const bool hit_before_other =
            (!shot.hit.wrapped && ab_wraps) || (shot.hit.wrapped == ab_wraps && d_hit < d_other);
        assert((hit_at_other || hit_before_other) &&
               "[C91 §2.2 Lemma 2.1]: the other chord endpoint is on ∂C, "
               "so the first hit cannot lie beyond it in the wrap metric");
        if (shot.success) {
            *success = make_success(ctx, from.edge_c, from.side, y_ab, shot.hit);
            return true;
        }
        if (hit_before_other) {
            *prime_exists = true;
            *prime_pos = boundary_position(*ctx.curve, shot.hit.edge, shot.hit.side, y_ab);
            *pe = shot.hit.edge;
            *ps = shot.hit.side;
            *px = shot.hit.x;
        } else if (!other_in) {
            if (fused_arc_contains(*ctx.submap, *ctx.curve, ctx.A2, other.edge_c, other.side,
                                   y_ab)) {
                RayHit h;
                h.hit = true;
                h.x = other.x;
                h.y = y_ab.y;
                h.edge = other.edge_c;
                h.side = other.side;
                h.wrapped = ab_wraps;
                h.hit_arc_idx = ctx.A2;
                *success = make_success(ctx, from.edge_c, from.side, y_ab, h);
                return true;
            }
            *prime_exists = true;
            *prime_pos = boundary_position(*ctx.curve, other.edge_c, other.side, y_ab);
            *pe = other.edge_c;
            *ps = other.side;
            *px = other.x;
        }

        return false;
    };

    if (a_in) {
        VisiblePoint success;
        if (endpoint_shot(A, B, b_in, &a_prime_exists, &a_prime_pos, &ap_edge, &ap_side, &ap_x,
                          &success))
            return success;
    }
    if (b_in) {
        VisiblePoint success;
        if (endpoint_shot(B, A, a_in, &b_prime_exists, &b_prime_pos, &bp_edge, &bp_side, &bp_x,
                          &success))
            return success;

        assert(a_prime_exists == b_prime_exists &&
               "[C91 §2.5 Lemma 2.4]: a' and b' exist together when both "
               "a, b ∈ B₁ ∪ B₂");
    }

    auto alpha_avoids_open = [&](const ClockwiseBoundaryPosition& x,
                                 const ClockwiseBoundaryPosition& y) {
        return in_closed_clockwise_interval(ctx.c_pos_a, y, x) &&
               in_closed_clockwise_interval(ctx.d_pos_a, y, x) &&
               precedes_from_clockwise_origin(y, ctx.c_pos_a, ctx.d_pos_a);
    };
    bool piece1_empty = alpha_avoids_open(A.pos_a, B.pos_a);
    bool piece2_empty = alpha_avoids_open(B.pos_a, A.pos_a);
    assert(!(piece1_empty && piece2_empty) &&
           "[C91 §3.2]: α has interior points, so some piece meets it");

    bool c_in_piece1 = in_half_open_clockwise_interval(ctx.c_pos_a, A.pos_a, B.pos_a);

    bool reject_piece1;
    if (piece1_empty) {
        reject_piece1 = true;
    } else if (piece2_empty) {
        reject_piece1 = false;
    } else {
        assert(a_in && "[C91 §2.5]: both B₁ and B₂ nonempty requires a ∈ B₁ ∪ B₂");
        bool a_eq_b = a_prime_exists && b_prime_exists && ap_edge == bp_edge &&
                      ap_side == bp_side && ap_x == bp_x;
        [[maybe_unused]] const std::size_t num_pieces =
            shielding_piece_count(a_prime_exists, b_prime_exists, a_eq_b);

        auto face_trav_up = [&](std::size_t e_, Side s_) {
            const auto& ed_ = ctx.curve->edge(e_);
            bool asc_ = symbolic_y_less(symbolic_y_of(ctx.curve->vertex(ed_.start_idx)),
                                        symbolic_y_of(ctx.curve->vertex(ed_.end_idx)));
            return (s_ == LEFT) ? asc_ : !asc_;
        };
        const bool piece1_side_up = face_trav_up(A.edge_c, A.side);
        bool first_on_piece1_side;
        if (a_prime_exists) {
            first_on_piece1_side = (face_trav_up(ap_edge, ap_side) == piece1_side_up);
        } else if (compare_boundary_positions(ctx.c_pos_a, A.pos_a) == 0) {
            first_on_piece1_side = !c_in_piece1;
        } else {
            first_on_piece1_side = c_in_piece1;
        }

        auto at_or_before_d = [&](const ClockwiseBoundaryPosition& u,
                                  const ClockwiseBoundaryPosition& v) {
            return precedes_from_clockwise_origin(ctx.d_pos_c, u, v);
        };
        std::size_t piece_idx;
        if (!a_prime_exists) {
            assert(!b_prime_exists);
            piece_idx = 0;
        } else if (!b_prime_exists) {
            piece_idx = at_or_before_d(a_prime_pos, ctx.a2_start) ? 0 : 1;

            if (piece_idx == 1)
                assert(at_or_before_d(ctx.a2_end, a_prime_pos) &&
                       "[C91 §2.5 Lemma 2.4]: A₂ must lie in one piece of A");
        } else if (a_eq_b) {
            piece_idx = at_or_before_d(a_prime_pos, ctx.a2_start) ? 0 : 1;
            if (piece_idx == 1)
                assert(at_or_before_d(ctx.a2_end, a_prime_pos));
        } else {
            assert(at_or_before_d(b_prime_pos, a_prime_pos) &&
                   "[C91 §2.5 Lemma 2.4]: piece order along A is c, a', b', d");
            if (at_or_before_d(a_prime_pos, ctx.a2_start)) {
                piece_idx = 0;
            } else if (at_or_before_d(b_prime_pos, ctx.a2_start)) {
                piece_idx = 1;
                assert(at_or_before_d(ctx.a2_end, a_prime_pos));
            } else {
                piece_idx = 2;
                assert(at_or_before_d(ctx.a2_end, b_prime_pos));
            }
        }
        assert(piece_idx < num_pieces &&
               "[C91 §2.5 Lemma 2.4]: A₂'s piece must be one of the 1–3 pieces");

        const std::size_t flips = a_eq_b ? 0 : piece_idx;
        const bool a2_on_piece1_side = first_on_piece1_side == (flips % 2 == 0);
        reject_piece1 = !a2_on_piece1_side;
    }

    const Chord::AdjArcs& a_adj = a_is_left_endpoint ? ab.left_adj : ab.right_adj;
    std::size_t probe_arc;
    bool probe_in_piece1;
    if (a_adj.count == 2) {
        auto starts_at_a = [&](std::size_t ai) {
            return arc_starts_at_chord_slot(Sa, Ca, ab, a_is_left_endpoint, ai);
        };
        bool first_starts = starts_at_a(a_adj.arcs[0]);
        probe_arc = first_starts ? a_adj.arcs[0] : a_adj.arcs[1];
        assert((first_starts || starts_at_a(a_adj.arcs[1])) &&
               "[C91 §2.4(ii)]: one adj slot must hold the arc starting "
               "at the chord endpoint");
        probe_in_piece1 = true;
    } else {
        probe_arc = a_adj.arcs[0];
        probe_in_piece1 = false;
    }
    std::size_t probe_region = Sa.arc(probe_arc).region_node;
    assert((probe_region == ab.region[0] || probe_region == ab.region[1]) &&
           "[C91 §2.4(ii)]: adj arc belongs to one of the chord's regions");
    bool probe_child_is_left = (probe_region == ab.region[0]);

    bool keep_piece1 = !reject_piece1;
    bool keep_probe_side = (keep_piece1 == probe_in_piece1);
    bool go_left = (keep_probe_side == probe_child_is_left);
    *next_node = go_left ? node.left_child : node.right_child;
    assert(*next_node != NONE && "[C91 §2.3]: internal TD node has two children");
    return fail;
}

VisiblePoint search_piece(const PieceSearchContext& ctx) {
    const Submap& Sa = *ctx.piece->submap;
    const Polygon& Ca = *ctx.piece->curve;

    const TreeDecomposition& td = Sa.tree_decomposition();
    assert(!td.empty() && "[C91 §3.0(ii)(3) tex 170]: normal-form S_α carries its tree "
                          "decomposition");

    std::size_t node_idx = td.root();
    while (td.node(node_idx).is_internal()) {
        std::size_t next = NONE;
        VisiblePoint vp = descend_step(ctx, td, node_idx, &next);
        if (vp.found)
            return vp;
        node_idx = next;
    }

    std::size_t leaf_region = td.node(node_idx).region_idx;
    RegionArcs arcs = collect_region_arcs(Sa, leaf_region);
    for (std::size_t k = 0; k < arcs.count; ++k) {
        const Arc& ra = Sa.arc(arcs.arcs[k]);

        assert(Sa.start_vertex != NONE && Sa.end_vertex != NONE &&
               "[C91 §2.4(iii)]: S_α's endpoints must be set");
        ArcSideRange ranges[3];
        std::size_t range_count = ra.side_ranges(Sa.start_vertex, Sa.end_vertex, ranges);
        for (std::size_t g = 0; g < range_count; ++g) {
            if (ranges[g].side != ctx.s)
                continue;
            for (std::size_t ee = ranges[g].first_edge; ee <= ranges[g].last_edge; ++ee) {
                for (std::size_t vv = ee; vv <= ee + 1; ++vv) {
                    SymbolicY vy = symbolic_y_of(Ca.vertex(vv));

                    {
                        std::size_t edge_c = ctx.edge_off + ctx.lo + ee;
                        std::size_t vidx_c =
                            symbolic_y_equal(vy, symbolic_y_of(ctx.curve->vertex(edge_c)))
                                ? edge_c
                                : edge_c + 1;
                        if (is_inside_companion(*ctx.curve, edge_c, ctx.s, vidx_c))
                            continue;
                    }

                    ClockwiseBoundaryPosition cand = boundary_position(Ca, ee, ctx.s, vy);
                    ClockwiseBoundaryPosition a_start =
                        boundary_position(Ca, ra.first_edge, ra.first_side,
                                          Sa.arc_start_symbolic_y(arcs.arcs[k], Ca));
                    ClockwiseBoundaryPosition a_end = boundary_position(
                        Ca, ra.last_edge, ra.last_side, Sa.arc_end_symbolic_y(arcs.arcs[k], Ca));
                    if (!in_closed_clockwise_interval(cand, a_start, a_end))
                        continue;

                    std::size_t edge_c = ctx.edge_off + ctx.lo + ee;
                    CandidateShot shot =
                        try_candidate_vertex(edge_c, ctx.s, vy, *ctx.submap, *ctx.curve, ctx.region,
                                             ctx.A2, *ctx.cycle, *ctx.fctx);
                    if (shot.success)
                        return make_success(ctx, edge_c, ctx.s, vy, shot.hit);
                }
            }
        }
    }

    return VisiblePoint{};
}

}

VisiblePoint find_visible_point(const Submap& submap, const Polygon& curve, std::size_t region,
                                std::size_t A1, std::size_t A2, const FusedRegionCycle& cycle,
                                const std::vector<ArcSource>& arc_sources,
                                const ConformalityOracles& oracles) {
    assert(submap.arc(A1).edge_count > 0 && submap.arc(A2).edge_count > 0 &&
           "[C91 §2.1 tex 70/72]: zero-length arcs are single ∂C points "
           "whose visibility is already realized — never candidates");

    const std::size_t n1e = oracles.first_curve->num_edges();

    FusedShootContext fctx;
    fctx.submap = &submap;
    fctx.curve = &curve;
    fctx.first_curve = oracles.first_curve;
    fctx.second_curve = oracles.second_curve;
    fctx.first_ray_shooter = oracles.first_ray_shooter;
    fctx.second_ray_shooter = oracles.second_ray_shooter;
    fctx.arc_sources = &arc_sources;

    ClockwiseBoundaryPosition a2_start = fused_arc_start(submap, curve, A2);
    ClockwiseBoundaryPosition a2_end = fused_arc_end(submap, curve, A2);

    auto on_A1 = [&](std::size_t edge_c, Side side, const SymbolicY& y) -> bool {
        return fused_arc_contains(submap, curve, A1, edge_c, side, y);
    };

    {
        const ArcSource& pr = arc_sources[A1];

        const std::size_t off = pr.on_first_curve ? 0 : n1e;
        const Polygon& input_curve =
            pr.on_first_curve ? *oracles.first_curve : *oracles.second_curve;
        const ArcCuttingOracle& cutter =
            pr.on_first_curve ? *oracles.first_arc_cutter : *oracles.second_arc_cutter;
        const std::size_t piece_count_bound =
            pr.on_first_curve ? oracles.first_piece_count_bound : oracles.second_piece_count_bound;
        const std::size_t piece_granularity_bound = pr.on_first_curve
                                                        ? oracles.first_piece_granularity_bound
                                                        : oracles.second_piece_granularity_bound;

        const Arc& a1_struct = submap.arc(A1);
        Subarc target;
        target.first_edge = a1_struct.first_edge - off;
        target.first_side = a1_struct.first_side;
        target.last_edge = a1_struct.last_edge - off;
        target.last_side = a1_struct.last_side;

        target.first_y = submap.arc_start_symbolic_y(A1, curve);
        target.last_y = submap.arc_end_symbolic_y(A1, curve);
        assert_subarc_clockwise(target);

        std::vector<ArcPiece> pieces = cutter.cut(pr.input_arc, target);

        assert_cut_postconditions(input_curve, target, pieces.data(), pieces.size(),
                                  piece_count_bound, piece_granularity_bound);

        for (const ArcPiece& piece : pieces) {
            if (piece.is_boundary_piece) {
                assert(piece.subarc.first_edge == piece.subarc.last_edge);
                const std::size_t e = piece.subarc.first_edge;
                const Side s = piece.subarc.first_side;
                for (std::size_t vv = e; vv <= e + 1; ++vv) {
                    SymbolicY vy = symbolic_y_of(input_curve.vertex(vv));
                    if (!on_A1(off + e, s, vy))
                        continue;

                    {
                        std::size_t vidx_c =
                            symbolic_y_equal(vy, symbolic_y_of(curve.vertex(off + e)))
                                ? off + e
                                : off + e + 1;
                        if (is_inside_companion(curve, off + e, s, vidx_c))
                            continue;
                    }
                    CandidateShot shot = try_candidate_vertex(off + e, s, vy, submap, curve, region,
                                                              A2, cycle, fctx);
                    if (shot.success) {
                        VisiblePoint vp;
                        vp.found = true;
                        vp.p_table_arc = A1;
                        vp.p_edge = off + e;
                        vp.p_side = s;
                        vp.p_x = edge_x_at_y(curve, off + e, vy);
                        vp.y = vy;
                        vp.q_table_arc = shot.hit.hit_arc_idx;
                        vp.q_edge = shot.hit.edge;
                        vp.q_side = shot.hit.side;
                        vp.q_x = shot.hit.x;
                        return vp;
                    }
                }
                continue;
            }

            const std::size_t lo = std::min(piece.subarc.first_edge, piece.subarc.last_edge);
            [[maybe_unused]] const std::size_t hi =
                std::max(piece.subarc.first_edge, piece.subarc.last_edge);
            assert(piece.subarc.first_side == piece.subarc.last_side &&
                   "[C91 §3.0(ii)(2) tex 170]: pieces do not double-back");
            const Side s = piece.subarc.first_side;

            assert(piece.curve->vertex(0).index == input_curve.vertex(lo).index &&
                   "[C91 §3.0(ii)(3) tex 170]: ᾱ must be the ascending "
                   "vertex-to-vertex subchain of Cᵢ over the piece's range");
            assert(piece.curve->num_edges() == hi - lo + 1);

            PieceSearchContext pctx;
            pctx.submap = &submap;
            pctx.curve = &curve;
            pctx.region = region;
            pctx.A1 = A1;
            pctx.A2 = A2;
            pctx.cycle = &cycle;
            pctx.fctx = &fctx;
            pctx.piece = &piece;
            pctx.input_curve = &input_curve;
            pctx.edge_off = off;
            pctx.lo = lo;
            pctx.s = s;
            pctx.a2_start = a2_start;
            pctx.a2_end = a2_end;

            const Polygon& Ca = *piece.curve;
            const std::size_t na = Ca.num_edges();
            if (s == LEFT) {
                pctx.c_pos_a = boundary_position(Ca, 0, LEFT, symbolic_y_of(Ca.vertex(0)));
                pctx.d_pos_a = boundary_position(Ca, na - 1, LEFT, symbolic_y_of(Ca.vertex(na)));
                pctx.c_pos_c =
                    boundary_position(curve, off + lo, LEFT, symbolic_y_of(Ca.vertex(0)));
                pctx.d_pos_c =
                    boundary_position(curve, off + lo + na - 1, LEFT, symbolic_y_of(Ca.vertex(na)));
            } else {
                pctx.c_pos_a = boundary_position(Ca, na - 1, RIGHT, symbolic_y_of(Ca.vertex(na)));
                pctx.d_pos_a = boundary_position(Ca, 0, RIGHT, symbolic_y_of(Ca.vertex(0)));
                pctx.c_pos_c = boundary_position(curve, off + lo + na - 1, RIGHT,
                                                 symbolic_y_of(Ca.vertex(na)));
                pctx.d_pos_c =
                    boundary_position(curve, off + lo, RIGHT, symbolic_y_of(Ca.vertex(0)));
            }

            VisiblePoint vp = search_piece(pctx);
            if (vp.found)
                return vp;
        }
    }

    return VisiblePoint{};
}

static void restore_regions(Submap& submap, const Polygon& curve,
                            const ConformalityOracles& oracles,
                            std::vector<ArcSource> arc_sources) {
    std::vector<std::vector<std::size_t>> region_arcs(submap.num_nodes());
    for (std::size_t ai = 0; ai < submap.num_arcs(); ++ai) {
        const Arc& a = submap.arc(ai);
        if (a.dead)
            continue;
        region_arcs[a.region_node].push_back(ai);
    }

    std::vector<std::size_t> work;
    for (std::size_t r = 0; r < submap.num_nodes(); ++r)
        if (!submap.node(r).dead && !region_arcs[r].empty())
            work.push_back(r);

    while (!work.empty()) {
        std::size_t r = work.back();
        work.pop_back();
        assert(!submap.node(r).dead && "[C91 §3.2]: regions are never removed here");

        FusedRegionCycle cycle = fused_region_cycle(submap, curve, r, region_arcs[r]);

#ifndef NDEBUG
        {
            std::size_t runs = 0;
            for (std::size_t li = 0; li < cycle.count; ++li) {
                std::size_t prev = (li + cycle.count - 1) % cycle.count;
                if (cycle.count == 1 || arc_sources[cycle.arcs[li].arc].on_first_curve !=
                                            arc_sources[cycle.arcs[prev].arc].on_first_curve)
                    ++runs;
            }
            assert(runs <= 2 && "[C91 §3.2 tex 238]: a fused region's arcs form at "
                                "most two runs, one per operand");
        }
#endif
        if (cycle.count <= 4)
            continue;

        bool found = false;
        const std::size_t k = cycle.count;
        for (std::size_t i = 0; i < k && !found; ++i) {
            if (cycle.arcs[i].is_zero_length)
                continue;
            for (std::size_t j = 0; j < k && !found; ++j) {
                if (j == i || j == (i + 1) % k || j == (i + k - 1) % k)
                    continue;
                if (cycle.arcs[j].is_zero_length)
                    continue;

                VisiblePoint vp =
                    find_visible_point(submap, curve, r, cycle.arcs[i].arc, cycle.arcs[j].arc,
                                       cycle, arc_sources, oracles);
                if (!vp.found)
                    continue;

#ifdef CHAZELLE_EXPENSIVE_ASSERTS
                {
                    Point pp{vp.p_x, vp.y.y, vp.y.tag};
                    RayHit h = naive_first_contact(curve, pp, vp.y,
                                                   shooting_direction(vp.p_edge, vp.p_side, curve),
                                                   vp.p_edge);
                    [[maybe_unused]] bool ok = h.hit && h.x == vp.q_x;
                    assert(ok && "[C91 §3.2 tex 264]: Lemma 3.2's point must "
                                 "actually see its partner w.r.t. C");
                }
#endif

                std::vector<std::size_t> flat;
                flat.reserve(k);
                for (std::size_t li = 0; li < k; ++li)
                    flat.push_back(cycle.arcs[li].arc);

                Submap::ChordPointSpec p{vp.p_table_arc, vp.p_edge, vp.p_side, vp.p_x};
                Submap::ChordPointSpec q{vp.q_table_arc, vp.q_edge, vp.q_side, vp.q_x};
                auto res = submap.insert_chord(p, q, vp.y, r, flat.data(), flat.size(), curve);

                arc_sources.resize(submap.num_arcs());
                arc_sources[res.p_after_arc] = arc_sources[vp.p_table_arc];
                arc_sources[res.q_after_arc] = arc_sources[vp.q_table_arc];

                region_arcs.resize(submap.num_nodes());
                region_arcs[r].clear();
                region_arcs[res.new_region].clear();
                flat.push_back(res.p_after_arc);
                flat.push_back(res.q_after_arc);
                for (std::size_t ai : flat)
                    region_arcs[submap.arc(ai).region_node].push_back(ai);

                work.push_back(r);
                work.push_back(res.new_region);
                found = true;
            }
        }

        assert(found && "[C91 §3.2 Lemma 3.3]: a region with more than four arcs "
                        "must admit a chord between nonconsecutive arcs");
    }

#ifndef NDEBUG
    for (std::size_t r = 0; r < submap.num_nodes(); ++r) {
        if (submap.node(r).dead || region_arcs.size() <= r)
            continue;
        if (region_arcs[r].empty())
            continue;
        FusedRegionCycle cycle = fused_region_cycle(submap, curve, r, region_arcs[r]);
        assert(cycle.count <= 4 && "[C91 §3.2 tex 264]: every region must end with ≤ 4 arcs");
    }
#endif
    assert(submap.is_conformal() && "[C91 §2.3 tex 114]: S must be conformal after §3.2");
    submap.assert_tree_property();
}

void restore_conformality(Submap& submap, const Polygon& curve,
                          const ConformalityOracles& oracles) {
    assert(oracles.first_submap && oracles.second_submap && oracles.first_curve &&
           oracles.second_curve && oracles.first_ray_shooter && oracles.second_ray_shooter &&
           oracles.first_arc_cutter && oracles.second_arc_cutter &&
           "[C91 §3.0 tex 166–170]: merging requires both input oracles");
    auto sources = identify_arc_sources(submap, curve, *oracles.first_submap, *oracles.first_curve,
                                        *oracles.second_submap, *oracles.second_curve);
    restore_regions(submap, curve, oracles, std::move(sources));
}

void restore_conformality(Submap& submap, const Polygon& curve,
                          const RayShootingOracle& ray_shooter, const ArcCuttingOracle& arc_cutter,
                          std::size_t piece_count_bound, std::size_t piece_granularity_bound) {
    std::vector<ArcSource> sources;
    sources.reserve(submap.num_arcs());
    for (std::size_t arc = 0; arc < submap.num_arcs(); ++arc)
        sources.push_back({true, arc});
    const ConformalityOracles oracles{
        .first_submap = &submap,
        .first_curve = &curve,
        .first_ray_shooter = &ray_shooter,
        .first_arc_cutter = &arc_cutter,
        .first_piece_count_bound = piece_count_bound,
        .first_piece_granularity_bound = piece_granularity_bound,
    };
    restore_regions(submap, curve, oracles, std::move(sources));
}

}
