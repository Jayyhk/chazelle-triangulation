#include "ray_shooting.h"
#include "../submap/boundary_geometry.h"
#include "../trace.h"

#include <algorithm>
#include <array>

namespace chazelle::animation {

namespace {

struct BoundaryPosition {
    std::size_t edge = NONE;
    SymbolicY y{};
};

bool edge_ascends(const Polygon& curve, std::size_t e) {
    const auto& ed = curve.edge(e);
    return symbolic_y_less(symbolic_y_of(curve.vertex(ed.start_idx)),
                           symbolic_y_of(curve.vertex(ed.end_idx)));
}

int compare_boundary_positions(const Polygon& curve, const BoundaryPosition& a,
                               const BoundaryPosition& b) {
    if (a.edge != b.edge)
        return a.edge < b.edge ? -1 : 1;
    if (symbolic_y_equal(a.y, b.y))
        return 0;
    bool asc = edge_ascends(curve, a.edge);
    bool less = asc ? symbolic_y_less(a.y, b.y) : symbolic_y_greater(a.y, b.y);
    return less ? -1 : 1;
}

struct BoundarySideInterval {
    Side side = LEFT;
    std::size_t lo_edge = NONE, hi_edge = NONE;
    BoundaryPosition lower_position, upper_position;
};

std::size_t decompose_arc_boundary(const Submap& submap, const Polygon& curve, std::size_t ai,
                                   BoundarySideInterval out[3]) {
    const Arc& a = submap.arc(ai);
    assert(submap.start_vertex != NONE && submap.end_vertex != NONE &&
           "[C91 §2.4(iii)]: C endpoints must be set");
    ArcSideRange ranges[3];
    std::size_t n = a.side_ranges(submap.start_vertex, submap.end_vertex, ranges);

    SymbolicY sy_start = submap.arc_start_symbolic_y(ai, curve);
    SymbolicY sy_end = submap.arc_end_symbolic_y(ai, curve);
    SymbolicY start_wrap_y = symbolic_y_of(curve.vertex(submap.start_vertex));
    SymbolicY end_wrap_y = symbolic_y_of(curve.vertex(submap.end_vertex));

    for (std::size_t i = 0; i < n; ++i) {
        BoundarySideInterval l;
        l.side = ranges[i].side;
        l.lo_edge = ranges[i].first_edge;
        l.hi_edge = ranges[i].last_edge;

        bool first_interval = (i == 0);
        bool last_interval = (i + 1 == n);
        if (l.side == LEFT) {
            l.lower_position = {l.lo_edge, first_interval ? sy_start : start_wrap_y};
            l.upper_position = {l.hi_edge, last_interval ? sy_end : end_wrap_y};
        } else {
            l.lower_position = {l.lo_edge, last_interval ? sy_end : start_wrap_y};
            l.upper_position = {l.hi_edge, first_interval ? sy_start : end_wrap_y};
        }
        out[i] = l;
    }
    return n;
}

bool arc_covers(const Submap& submap, std::size_t ai, std::size_t edge, Side side) {
    assert(submap.start_vertex != NONE && submap.end_vertex != NONE);
    return submap.arc(ai).covers(edge, side, submap.start_vertex, submap.end_vertex);
}

Side struck_side(const Polygon& curve, std::size_t e, Side dir) {
    Side minus_x = edge_ascends(curve, e) ? LEFT : RIGHT;
    return (dir == RIGHT) ? minus_x : (minus_x == LEFT ? RIGHT : LEFT);
}

struct NearestRayHit {
    bool hit = false;
    Exact x = 0.0;
    std::size_t edge = NONE;
    Side side = LEFT;
    bool wrapped = false;
    Exact signed_distance = 0.0;
    SourceOffset source_x_offset = SOURCE_OFFSET_NONE;

    void offer(const Polygon& curve, const SymbolicY& sy, Side dir, const Exact& candidate_x,
               std::size_t candidate_edge, Side candidate_side, const Exact& candidate_distance,
               std::size_t query) {
        bool candidate_wraps =
            (candidate_distance < 0.0) ||
            (candidate_distance == 0.0 &&
             !perturbed_hit_forward(curve, sy, dir, source_x_offset, candidate_edge));
        bool better;
        if (!hit)
            better = true;
        else if (candidate_wraps != wrapped)
            better = !candidate_wraps;
        else if (candidate_distance != signed_distance)
            better = candidate_distance < signed_distance;
        else
            better =
                ray_contact_precedes(curve, sy, dir, candidate_edge, candidate_side, edge, side);
        if (auto* trace = AnimationTrace::current(); trace && query != NONE)
            trace->point("search_candidate", curve, {candidate_x, sy.y, sy.tag},
                         {{"query", query},
                          {"edge", candidate_edge},
                          {"accepted", better},
                          {"wrapped", candidate_wraps}});
        if (better) {
            hit = true;
            x = candidate_x;
            edge = candidate_edge;
            side = candidate_side;
            wrapped = candidate_wraps;
            signed_distance = candidate_distance;
        }
    }
};

void scan_edge_range(const Polygon& curve, std::size_t lo, std::size_t hi, const Point& p,
                     const SymbolicY& sy, Side dir, NearestRayHit& best,
                     const BoundarySideInterval* clip = nullptr, std::size_t query = NONE) {
    std::size_t e = lo;
    while (e <= hi) {
        if (auto* trace = AnimationTrace::current(); trace && query != NONE)
            trace->record("search_edge", {{"query", query}, {"edge", e}});
        const std::size_t nn = curve.next_nonnull_edge(e);
        if (nn > e) {
            const std::size_t run_hi = std::min(hi, nn - 1);
            const Point& va = curve.vertex(e);
            const std::size_t vi = curve.local_index_of_tag(sy.tag);
            if (sy.y == va.y && vi != NONE && vi >= e && vi <= run_hi + 1) {
                assert(symbolic_y_equal(sy, symbolic_y_of(curve.vertex(vi))) &&
                       "[C91 §2.4 tex 133]: a symbolic tag names its input vertex");
                const Exact x = curve.vertex(vi).x;
                const Exact signed_distance = (dir == RIGHT) ? (x - p.x) : (p.x - x);
                auto offer_if_inside = [&](std::size_t candidate_edge) {
                    const BoundaryPosition contact{candidate_edge, sy};
                    if (!clip ||
                        (compare_boundary_positions(curve, clip->lower_position, contact) <= 0 &&
                         compare_boundary_positions(curve, contact, clip->upper_position) <= 0))
                        best.offer(curve, sy, dir, x, candidate_edge,
                                   struck_side(curve, candidate_edge, dir), signed_distance, query);
                };
                if (vi > e)
                    offer_if_inside(vi - 1);
                if (vi <= run_hi)
                    offer_if_inside(vi);
            }
            e = run_hi + 1;
            continue;
        }
        Exact x;
        const BoundaryPosition contact{e, sy};
        if ((!clip || (compare_boundary_positions(curve, clip->lower_position, contact) <= 0 &&
                       compare_boundary_positions(curve, contact, clip->upper_position) <= 0)) &&
            edge_crossing_x(curve, e, sy, &x)) {
            Exact signed_distance = (dir == RIGHT) ? (x - p.x) : (p.x - x);
            best.offer(curve, sy, dir, x, e, struck_side(curve, e, dir), signed_distance, query);
        }
        ++e;
    }
}

void scan_region(const Submap& submap, const Polygon& curve, const std::vector<std::size_t>& arcs,
                 const Point& p, const SymbolicY& sy, Side dir, NearestRayHit& best,
                 std::size_t query) {
    for (std::size_t ai : arcs) {
        if (auto* trace = AnimationTrace::current())
            trace->record("search_arc", {{"query", query}, {"arc", ai}});
        BoundarySideInterval ranges[3];
        std::size_t n = decompose_arc_boundary(submap, curve, ai, ranges);
        for (std::size_t i = 0; i < n; ++i)
            scan_edge_range(curve, ranges[i].lo_edge, ranges[i].hi_edge, p, sy, dir, best,
                            &ranges[i], query);
    }
}

std::size_t chord_before_arc(const Submap& submap, const Polygon& curve, const Chord& c,
                             bool left_slot) {
    const Chord::AdjArcs& adj = left_slot ? c.left_adj : c.right_adj;
    if (adj.count == 1)
        return adj.arcs[0];
    assert(adj.count == 2);
    bool s0 = arc_starts_at_chord_slot(submap, curve, c, left_slot, adj.arcs[0]);
    assert(s0 != arc_starts_at_chord_slot(submap, curve, c, left_slot, adj.arcs[1]) &&
           "[C91 §2.2 tex 94/96, §2.4(ii)]: exactly one adj arc starts at a "
           "mid-edge chord endpoint");
    return s0 ? adj.arcs[1] : adj.arcs[0];
}
std::size_t chord_after_arc(const Submap& submap, const Polygon& curve, const Chord& c,
                            bool left_slot) {
    const Chord::AdjArcs& adj = left_slot ? c.left_adj : c.right_adj;
    if (adj.count == 1)
        return (adj.arcs[0] + 1) % submap.num_arcs();
    assert(adj.count == 2);
    bool s0 = arc_starts_at_chord_slot(submap, curve, c, left_slot, adj.arcs[0]);
    assert(s0 != arc_starts_at_chord_slot(submap, curve, c, left_slot, adj.arcs[1]) &&
           "[C91 §2.2 tex 94/96, §2.4(ii)]: exactly one adj arc starts at a "
           "mid-edge chord endpoint");
    return s0 ? adj.arcs[0] : adj.arcs[1];
}

}

RayShootingStructure::RayShootingStructure(const Submap& submap, const Polygon& curve,
                                           std::size_t granularity)
    : submap_(&submap), curve_(&curve), granularity_(granularity) {
#ifdef CHAZELLE_EXPENSIVE_ASSERTS
    submap.check_invariants(curve);
    assert(submap.is_granular(granularity, curve) &&
           "[C91 Lemma 3.6 tex 311]: S must be γ-granular");
#endif
    assert(submap.is_conformal() && "[C91 §3.4 tex 286]: S must be conformal");
    assert(submap.is_semigranular(granularity) &&
           "[C91 §3.4 tex 286]: S must be γ-granular — the naive "
           "region scans rely on O(γ) edges per region (tex 306)");

    if (auto* trace = AnimationTrace::current())
        animation_structure_ = trace->record("search_structure", {{"owner", trace->map_id(submap)},
                                                                  {"curve", trace->curve(curve)},
                                                                  {"granularity", granularity}});
    build_faces();
    if (auto* trace = AnimationTrace::current())
        trace->indices("search_faces", {{"structure", animation_structure_}}, "regions",
                       region_of_face_);
    if (face_count_ > 1) {
        build_dual_graph_and_decomposition();
        build_vertical_line();
    }
    if (auto* trace = AnimationTrace::current()) {
        for (std::size_t index = 0; index < vertical_line_crossings_.size(); ++index) {
            const auto& crossing = vertical_line_crossings_[index];
            trace->search_crossing(animation_structure_, curve, submap.chord(crossing.chord), index,
                                   crossing.region_below, crossing.region_above);
        }
        trace->indices(
            "search_ready",
            {{"structure", animation_structure_}, {"infinity_region", region_at_infinity_}},
            "subsets", separator_decomposition_.subset);
    }
}

void RayShootingStructure::build_faces() {
    const Submap& submap = *submap_;
    const Polygon& curve = *curve_;

    arcs_of_region_.assign(submap.num_nodes(), {});
    for (std::size_t ai = 0; ai < submap.num_arcs(); ++ai) {
        assert(!submap.arc(ai).dead && "[C91 §2.4]: normal form has no dead arcs");
        arcs_of_region_[submap.arc(ai).region_node].push_back(ai);
    }

    face_of_region_.assign(submap.num_nodes(), NONE);
    region_of_face_.clear();
    for (std::size_t r = 0; r < submap.num_nodes(); ++r) {
        if (submap.node(r).dead)
            continue;
        bool nonempty = false;
        for (std::size_t ai : arcs_of_region_[r]) {
            BoundarySideInterval ranges[3];
            std::size_t n = decompose_arc_boundary(submap, curve, ai, ranges);
            for (std::size_t i = 0; i < n && !nonempty; ++i)
                if (compare_boundary_positions(curve, ranges[i].lower_position,
                                               ranges[i].upper_position) != 0)
                    nonempty = true;
        }
        for (std::size_t ci : submap.node(r).incident_chords)
            if (!submap.chord(ci).dead && !submap.chord(ci).is_null_length)
                nonempty = true;
        if (nonempty) {
            face_of_region_[r] = region_of_face_.size();
            region_of_face_.push_back(r);
        }

        assert(arcs_of_region_[r].size() <= 4 &&
               "[C91 §2.3 tex 114]: conformal region has ≤ 4 arc-structures");
    }
    face_count_ = region_of_face_.size();
    assert(face_count_ >= 1 && "V(C) has at least one nonempty region");
}

void RayShootingStructure::build_dual_graph_and_decomposition() {
    const Submap& submap = *submap_;
    const Polygon& curve = *curve_;

    struct DualEdge {
        std::size_t fa, fb;
        bool dead = false;
        std::size_t first_arc = NONE, second_arc = NONE, chord = NONE;
    };
    std::vector<DualEdge> edges;

    struct Interval {
        BoundaryPosition lo, hi;
        std::size_t arc;
        std::size_t interval;
        std::size_t region;
    };
    std::vector<Interval> left_iv, right_iv;
    for (std::size_t ai = 0; ai < submap.num_arcs(); ++ai) {
        BoundarySideInterval ranges[3];
        std::size_t n = decompose_arc_boundary(submap, curve, ai, ranges);
        for (std::size_t i = 0; i < n; ++i) {
            if (compare_boundary_positions(curve, ranges[i].lower_position,
                                           ranges[i].upper_position) == 0)
                continue;
            Interval iv{ranges[i].lower_position, ranges[i].upper_position, ai, i,
                        submap.arc(ai).region_node};
            (ranges[i].side == LEFT ? left_iv : right_iv).push_back(iv);
        }
    }
    auto by_lo = [&](const Interval& a, const Interval& b) {
        int o = compare_boundary_positions(curve, a.lo, b.lo);
        if (o != 0)
            return o < 0;
        return compare_boundary_positions(curve, a.hi, b.hi) < 0;
    };
    std::sort(left_iv.begin(), left_iv.end(), by_lo);
    std::sort(right_iv.begin(), right_iv.end(), by_lo);

    auto retain = [&](const std::vector<Interval>& src, std::vector<BoundaryInterval>& dst) {
        dst.clear();
        dst.reserve(src.size());
        for (const Interval& iv : src)
            dst.push_back(BoundaryInterval{iv.lo.edge, iv.hi.edge, iv.lo.y, iv.hi.y, iv.region});
    };
    retain(left_iv, left_intervals_);
    retain(right_iv, right_intervals_);

    std::vector<std::array<std::vector<std::size_t>, 3>> interval_events(submap.num_arcs());
    {
        std::size_t i = 0, j = 0;
        while (i < left_iv.size() && j < right_iv.size()) {
            const Interval& a = left_iv[i];
            const Interval& b = right_iv[j];
            const BoundaryPosition& lo =
                (compare_boundary_positions(curve, a.lo, b.lo) >= 0) ? a.lo : b.lo;
            const BoundaryPosition& hi =
                (compare_boundary_positions(curve, a.hi, b.hi) <= 0) ? a.hi : b.hi;
            if (compare_boundary_positions(curve, lo, hi) < 0) {
                std::size_t fa = face_of_region_[a.region];
                std::size_t fb = face_of_region_[b.region];
                assert(fa != NONE && fb != NONE &&
                       "positive-length boundary implies nonempty regions");
                if (fa != fb) {
                    std::size_t id = edges.size();
                    edges.push_back({fa, fb, false, a.arc, b.arc, NONE});
                    interval_events[a.arc][a.interval].push_back(id);
                    interval_events[b.arc][b.interval].push_back(id);
                }
            }
            if (compare_boundary_positions(curve, a.hi, b.hi) <= 0)
                ++i;
            else
                ++j;
        }
    }

    std::vector<std::size_t> chord_edge(submap.num_chords(), NONE);
    for (std::size_t ci = 0; ci < submap.num_chords(); ++ci) {
        const Chord& c = submap.chord(ci);
        assert(!c.dead && "[C91 §2.4]: normal form has no dead chords");
        if (c.is_null_length)
            continue;
        std::size_t fa = face_of_region_[c.region[0]];
        std::size_t fb = face_of_region_[c.region[1]];
        assert(fa != NONE && fb != NONE && "a positive-length chord bounds two nonempty regions");
        assert(fa != fb && "[C91 §2.2]: a chord separates two regions");
        chord_edge[ci] = edges.size();
        edges.push_back({fa, fb, false, NONE, NONE, ci});
    }

    std::vector<std::vector<std::size_t>> rot(face_count_);
    for (std::size_t f = 0; f < face_count_; ++f) {
        std::size_t r = region_of_face_[f];
        const auto& arcs = arcs_of_region_[r];
        assert(!arcs.empty() && "a face's region has boundary arcs");
        std::vector<std::size_t> sorted_arcs(arcs);
        std::sort(sorted_arcs.begin(), sorted_arcs.end());

        std::size_t start = sorted_arcs[0];
        std::size_t cur = start;
        std::size_t emitted = 0;
        do {
            ++emitted;
            assert(emitted <= 2 * arcs.size() + 2 && "region boundary cycle must close");

            {
                ArcSideRange ranges[3];
                std::size_t n =
                    submap.arc(cur).side_ranges(submap.start_vertex, submap.end_vertex, ranges);
                for (std::size_t li = 0; li < n; ++li) {
                    const auto& evs = interval_events[cur][li];
                    if (ranges[li].side == LEFT) {
                        for (std::size_t id : evs)
                            rot[f].push_back(id);
                    } else {
                        for (std::size_t k = evs.size(); k-- > 0;)
                            rot[f].push_back(evs[k]);
                    }
                }
            }

            std::size_t via_chord = NONE;
            bool via_left_slot = false;
            for (std::size_t ci : submap.node(r).incident_chords) {
                const Chord& c = submap.chord(ci);
                for (bool left_slot : {true, false}) {
                    if (chord_before_arc(submap, curve, c, left_slot) != cur)
                        continue;

                    assert(via_chord == NONE && "an arc ends at exactly one chord endpoint");
                    via_chord = ci;
                    via_left_slot = left_slot;
                }
            }

            assert(via_chord != NONE && "[C91 §2.2 tex 96]: every arc ends at a chord "
                                        "endpoint of its region");
            if (chord_edge[via_chord] != NONE)
                rot[f].push_back(chord_edge[via_chord]);
            std::size_t next =
                chord_after_arc(submap, curve, submap.chord(via_chord), !via_left_slot);
            assert(submap.arc(next).region_node == r &&
                   "[C91 §2.2 tex 96]: the boundary cycle stays in the "
                   "region");
            cur = next;
        } while (cur != start);
        assert(emitted == arcs.size() && "the boundary cycle visits every arc of the region once");
    }

    {
        std::vector<std::size_t> order(edges.size());
        for (std::size_t i = 0; i < order.size(); ++i)
            order[i] = i;
        auto key = [&](std::size_t e) {
            auto [x, y] = std::minmax(edges[e].fa, edges[e].fb);
            return std::pair<std::size_t, std::size_t>(x, y);
        };
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return key(a) < key(b) || (key(a) == key(b) && a < b);
        });
        for (std::size_t i = 1; i < order.size(); ++i)
            if (key(order[i]) == key(order[i - 1]))
                edges[order[i]].dead = true;
    }
    std::vector<std::size_t> edge_map(edges.size(), NONE);
    std::vector<std::pair<std::size_t, std::size_t>> live_edges;
    for (std::size_t e = 0; e < edges.size(); ++e) {
        if (edges[e].dead)
            continue;
        edge_map[e] = live_edges.size();
        live_edges.emplace_back(edges[e].fa, edges[e].fb);
        dual_edges_.emplace_back(edges[e].fa, edges[e].fb);
    }
    std::vector<std::vector<std::size_t>> live_rot(face_count_);
    for (std::size_t f = 0; f < face_count_; ++f)
        for (std::size_t id : rot[f])
            if (edge_map[id] != NONE)
                live_rot[f].push_back(edge_map[id]);

    assert((face_count_ < 3 || live_edges.size() <= 3 * face_count_ - 6) &&
           "[C91 §3.4 tex 297]: |E(G)| ≤ 3μ − 6");

    EmbeddedPlanarGraph G(face_count_, live_edges, live_rot);

    {
        std::vector<bool> vis(face_count_, false);
        std::vector<std::size_t> q{0};
        vis[0] = true;
        for (std::size_t qi = 0; qi < q.size(); ++qi)
            for (std::size_t h : G.incident_halves(q[qi])) {
                std::size_t w = G.half_to(h);
                if (!vis[w]) {
                    vis[w] = true;
                    q.push_back(w);
                }
            }
        assert(q.size() == face_count_ && "[C91 §3.4 tex 295]: the dual graph G is connected");
    }

    if (auto* trace = AnimationTrace::current()) {
        for (const auto& edge : edges) {
            if (!edge.dead)
                trace->record("search_graph_edge", {{"structure", animation_structure_},
                                                    {"first", edge.fa},
                                                    {"second", edge.fb},
                                                    {"first_arc", edge.first_arc},
                                                    {"second_arc", edge.second_arc},
                                                    {"chord", edge.chord}});
        }
    }
    separator_decomposition_ = build_separator_decomposition(G, animation_structure_);

    subset_faces_.assign(separator_decomposition_.num_subsets, {});
    for (std::size_t f = 0; f < face_count_; ++f) {
        if (separator_decomposition_.subset[f] == NONE)
            separator_faces_.push_back(f);
        else
            subset_faces_[separator_decomposition_.subset[f]].push_back(f);
    }
}

void RayShootingStructure::build_vertical_line() {
    const Submap& submap = *submap_;
    const Polygon& curve = *curve_;

    std::vector<std::size_t> wrapped;
    for (std::size_t ci = 0; ci < submap.num_chords(); ++ci) {
        const Chord& c = submap.chord(ci);
        if (c.is_null_length)
            continue;
        if (chord_runs_through_infinity(curve, c))
            wrapped.push_back(ci);
    }

    if (wrapped.empty()) {
        std::size_t vt = curve.max_y_vertex();
        if (curve.is_endpoint(vt)) {
            region_at_infinity_ = (vt == 0) ? submap.arc(submap.start_arc).region_node
                                            : submap.arc(submap.end_arc).region_node;
        } else {
            assert(curve.is_y_extremum(vt) && "the interior global y-max is a local extremum");

            const Point& v = curve.vertex(vt);
            const bool prev_left = curve.previous_branch_left(vt);

            const Side outside_in = prev_left ? LEFT : RIGHT;
            auto ids = submap.double_identify(vt - 1, symbolic_y_of(v), curve);
            std::size_t found = NONE;
            for (std::size_t ai : ids)
                if (arc_covers(submap, ai, vt - 1, outside_in)) {
                    found = ai;
                    break;
                }
            assert(found != NONE && "an arc covers the outside face at the global maximum");
            region_at_infinity_ = submap.arc(found).region_node;
        }
        assert(face_of_region_[region_at_infinity_] != NONE && "the polar region is nonempty");
        return;
    }

    struct Cross {
        std::size_t chord, below, above;
    };
    std::vector<Cross> candidate_x;
    std::vector<std::size_t> at_region(submap_->num_nodes(), 0);
    for (std::size_t ci : wrapped) {
        Cross c{ci, NONE, NONE};
        submap.chord_regions_below_above(ci, curve, &c.below, &c.above);
        candidate_x.push_back(c);
        ++at_region[c.below];
        ++at_region[c.above];
        assert(at_region[c.below] <= 2 && at_region[c.above] <= 2 &&
               "[C91 §3.4 tex 306]: the cut regions lie on a path");
    }

    std::vector<std::size_t> next_by_below(submap_->num_nodes(), NONE);
    for (std::size_t k = 0; k < candidate_x.size(); ++k) {
        assert(next_by_below[candidate_x[k].below] == NONE);
        next_by_below[candidate_x[k].below] = k;
    }
    std::size_t first = NONE;
    for (std::size_t k = 0; k < candidate_x.size(); ++k)
        if (at_region[candidate_x[k].below] == 1) {
            assert(first == NONE && "exactly one south-cap crossing exists");
            first = k;
        }
    assert(first != NONE);
    std::size_t curk = first;
    while (true) {
        const Cross& c = candidate_x[curk];
        LineCrossing lc;
        lc.chord = c.chord;
        lc.y = submap.chord(c.chord).symbolic_y();
        lc.region_below = c.below;
        lc.region_above = c.above;
        if (!vertical_line_crossings_.empty()) {
            assert(vertical_line_crossings_.back().region_above == c.below &&
                   "[C91 §3.4 tex 306]: consecutive segments share their "
                   "crossing chord's region");
            assert(symbolic_y_less(vertical_line_crossings_.back().y, lc.y) &&
                   "[C91 §3.4 tex 306]: sorting the intersections comes "
                   "for free");
        }
        vertical_line_crossings_.push_back(lc);
        std::size_t nk = next_by_below[c.above];
        if (nk == NONE)
            break;
        curk = nk;
    }
    assert(vertical_line_crossings_.size() == candidate_x.size() &&
           "every crossing lies on the single bottom-to-top path");
}

void RayShootingStructure::regions_at_boundary(std::size_t edge, Side side, const SymbolicY& y,
                                               std::vector<std::size_t>& out,
                                               std::size_t query) const {
    const Polygon& curve = *curve_;
    const auto& list = (side == LEFT) ? left_intervals_ : right_intervals_;
    assert(!list.empty());
    BoundaryPosition pos{edge, y};

    std::size_t lo = 0, hi = list.size();
    while (lo < hi) {
        std::size_t mid = (lo + hi) / 2;
        BoundaryPosition mlo{list[mid].lo_edge, list[mid].lo_y};
        const bool follows = compare_boundary_positions(curve, mlo, pos) <= 0;
        if (auto* trace = AnimationTrace::current())
            trace->record("boundary_search", {{"query", query},
                                              {"region", list[mid].region},
                                              {"lo", lo},
                                              {"hi", hi},
                                              {"mid", mid},
                                              {"follows", follows}});
        if (follows)
            lo = mid + 1;
        else
            hi = mid;
    }
    auto contains = [&](std::size_t i) {
        BoundaryPosition ilo{list[i].lo_edge, list[i].lo_y};
        BoundaryPosition ihi{list[i].hi_edge, list[i].hi_y};
        return compare_boundary_positions(curve, ilo, pos) <= 0 &&
               compare_boundary_positions(curve, pos, ihi) <= 0;
    };
    auto push = [&](std::size_t r) {
        for (std::size_t x : out)
            if (x == r)
                return;
        out.push_back(r);
        if (auto* trace = AnimationTrace::current())
            trace->record("boundary_identify", {{"query", query}, {"region", r}});
    };

    for (std::size_t i = lo; i-- > 0 && contains(i);)
        push(list[i].region);
    assert(!out.empty() && "[C91 §2.4 tex 144]: every ∂C contact identifies a region");
}

RayHit RayShootingStructure::shoot_toward_boundary(const Point& p, Side dir,
                                                   const SourceOffset& source_x_offset) const {
    const Submap& submap = *submap_;
    const Polygon& curve = *curve_;
    SymbolicY sy{p.y, p.index};
    auto* trace = AnimationTrace::current();
    const std::size_t query = trace ? trace->point("search_begin", curve, p,
                                                   {{"structure", animation_structure_},
                                                    {"direction", static_cast<std::size_t>(dir)}})
                                    : NONE;

    auto to_rayhit = [&](const NearestRayHit& c) {
        RayHit h;
        if (!c.hit) {
            if (trace) {
                const auto result = trace->event_count();
                trace->ray(curve, p, dir, h);
                trace->record("search_end", {{"query", query}, {"result", result}});
            }
            return h;
        }
        h.hit = true;
        h.x = c.x;
        h.y = p.y;
        h.edge = c.edge;
        h.side = c.side;
        h.wrapped = c.wrapped;
        if (trace) {
            const auto result = trace->event_count();
            trace->ray(curve, p, dir, h);
            trace->record("search_end", {{"query", query}, {"result", result}});
        }
        return h;
    };

    if (face_count_ <= 1) {
        NearestRayHit best;
        best.source_x_offset = source_x_offset;
        scan_edge_range(curve, 0, curve.num_edges() - 1, p, sy, dir, best, nullptr, query);
        return to_rayhit(best);
    }

    std::vector<std::size_t> subsets;
    auto add_subset = [&](std::size_t region) {
        std::size_t f = face_of_region_[region];
        assert(f != NONE);
        std::size_t sub = separator_decomposition_.subset[f];

        if (sub == NONE)
            return;
        for (std::size_t s : subsets)
            if (s == sub)
                return;
        subsets.push_back(sub);
        if (trace)
            trace->record("search_subset", {{"query", query}, {"subset", sub}, {"region", region}});
    };

    NearestRayHit dstar_best;
    dstar_best.source_x_offset = source_x_offset;
    for (std::size_t f : separator_faces_) {
        if (trace)
            trace->record("search_scan", {{"query", query}, {"face", f}, {"separator", 1}});
        scan_region(submap, curve, arcs_of_region_[region_of_face_[f]], p, sy, dir, dstar_best,
                    query);
    }

    if (dstar_best.hit) {
        std::vector<std::size_t> near;
        regions_at_boundary(dstar_best.edge, dstar_best.side, sy, near, query);
        bool all_dstar = true;
        for (std::size_t r : near) {
            std::size_t f = face_of_region_[r];
            assert(f != NONE && "a struck arc bounds a nonempty region");
            if (separator_decomposition_.subset[f] == NONE)
                continue;
            all_dstar = false;
            add_subset(r);
        }
        if (all_dstar) {
            return to_rayhit(dstar_best);
        }
    } else {
        if (vertical_line_crossings_.empty()) {
            add_subset(region_at_infinity_);
        } else {
            std::size_t lo = 0, hi = vertical_line_crossings_.size();
            while (lo < hi) {
                std::size_t mid = (lo + hi) / 2;
                const bool below = symbolic_y_less(sy, vertical_line_crossings_[mid].y);
                if (trace)
                    trace->record(
                        "vertical_search",
                        {{"query", query}, {"lo", lo}, {"hi", hi}, {"mid", mid}, {"below", below}});
                if (below)
                    hi = mid;
                else
                    lo = mid + 1;
            }

            if (lo > 0 && symbolic_y_equal(sy, vertical_line_crossings_[lo - 1].y)) {
                add_subset(vertical_line_crossings_[lo - 1].region_below);
                add_subset(vertical_line_crossings_[lo - 1].region_above);
            } else if (lo == 0) {
                add_subset(vertical_line_crossings_[0].region_below);
            } else if (lo == vertical_line_crossings_.size()) {
                add_subset(vertical_line_crossings_.back().region_above);
            } else {
                assert(vertical_line_crossings_[lo - 1].region_above ==
                       vertical_line_crossings_[lo].region_below);
                add_subset(vertical_line_crossings_[lo].region_below);
            }
        }
    }

    NearestRayHit best = dstar_best;
    for (std::size_t sub : subsets) {
        const auto& members = subset_faces_[sub];
#ifndef NDEBUG
        __extension__ typedef unsigned __int128 u128;
        assert((u128)members.size() * members.size() * members.size() <=
               (u128)face_count_ * face_count_);
#endif
        for (std::size_t f : members) {
            if (trace)
                trace->record("search_scan", {{"query", query}, {"face", f}, {"separator", 0}});
            scan_region(submap, curve, arcs_of_region_[region_of_face_[f]], p, sy, dir, best,
                        query);
        }
    }

    if (!best.hit) {
        assert((symbolic_y_greater(sy, symbolic_y_of(curve.vertex(curve.max_y_vertex()))) ||
                symbolic_y_less(sy, symbolic_y_of(curve.vertex(curve.min_y_vertex())))) &&
               "[C91 §2.1 tex 70]: a wrapping ray inside C's y-range "
               "must hit C");
    }
    return to_rayhit(best);
}

}
