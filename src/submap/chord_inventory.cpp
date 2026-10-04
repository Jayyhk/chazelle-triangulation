#include "chord_inventory.h"
#include "boundary_geometry.h"

#ifdef CHAZELLE_EXPENSIVE_ASSERTS
#include "../visibility/naive_visibility.h"
#endif

#include <algorithm>
#include <utility>

namespace chazelle {

void canonicalize_chord(PendingChord& p, const Polygon& curve) {
    const std::size_t n_edges = curve.num_edges();
    auto canonicalize_vertex_label = [&](std::size_t& edge_c, const SymbolicY& y) {
        if (edge_c == 0)
            return;
        const auto& e = curve.edge(edge_c);

        if (!symbolic_y_equal(y, symbolic_y_of(curve.vertex(e.start_idx))))
            return;
        const std::size_t vidx = e.start_idx;
        assert(vidx == edge_c && "[C91 §2.4 tex 133]: edge k spans vertices k → k+1");
        assert(vidx + 1 < curve.num_vertices());
        if (is_local_y_extremum(curve.vertex(vidx - 1), curve.vertex(vidx), curve.vertex(vidx + 1)))
            return;
        edge_c = edge_c - 1;
    };
    canonicalize_vertex_label(p.left_edge_c, p.y);
    canonicalize_vertex_label(p.right_edge_c, p.y);

    auto tie_pos = [&](std::size_t edge, Side side) {
        return (side == LEFT) ? edge : (2 * n_edges - 1 - edge);
    };
    Exact xl = edge_x_at_y(curve, p.left_edge_c, p.y);
    Exact xr = edge_x_at_y(curve, p.right_edge_c, p.y);
    bool swap_slots = xl > xr || (xl == xr && tie_pos(p.left_edge_c, p.left_side) >
                                                  tie_pos(p.right_edge_c, p.right_side));
    if (swap_slots) {
        std::swap(p.left_edge_c, p.right_edge_c);
        std::swap(p.left_side, p.right_side);
    }
}

static bool chord_less(const PendingChord& a, const PendingChord& b) {
    if (a.y.y != b.y.y)
        return a.y.y < b.y.y;
    if (a.y.tag != b.y.tag)
        return a.y.tag < b.y.tag;
    if (a.left_edge_c != b.left_edge_c)
        return a.left_edge_c < b.left_edge_c;
    if (a.left_side != b.left_side)
        return a.left_side < b.left_side;
    if (a.right_edge_c != b.right_edge_c)
        return a.right_edge_c < b.right_edge_c;
    if (a.right_side != b.right_side)
        return a.right_side < b.right_side;
    return a.is_null_length < b.is_null_length;
}

bool chord_endpoint_precedes(const Polygon& curve, const std::vector<PendingChord>& pending,
                             const ChordEndpoint& first, const ChordEndpoint& second) {
    const std::size_t n_edges = curve.num_edges();
    auto trav_pos = [&](std::size_t edge, Side side) -> std::size_t {
        return (side == LEFT) ? edge : (2 * n_edges - 1 - edge);
    };
    auto edge_ascends = [&](std::size_t edge) -> bool {
        const auto& e = curve.edge(edge);
        return symbolic_y_less(symbolic_y_of(curve.vertex(e.start_idx)),
                               symbolic_y_of(curve.vertex(e.end_idx)));
    };

    auto lin_less = [&](std::size_t pa, const SymbolicY& ya, std::size_t pb,
                        const SymbolicY& yb) -> bool {
        if (pa != pb)
            return pa < pb;
        if (symbolic_y_equal(ya, yb))
            return false;

        const std::size_t e = (pa < n_edges) ? pa : (2 * n_edges - 1 - pa);
        const Side sd = (pa < n_edges) ? LEFT : RIGHT;
        bool trav_asc = (sd == LEFT) == edge_ascends(e);
        return trav_asc ? symbolic_y_less(ya, yb) : symbolic_y_greater(ya, yb);
    };
    auto partner_pos = [&](const ChordEndpoint& e, std::size_t* pp, SymbolicY* py) {
        const auto& p = pending[e.pending_idx];
        const std::size_t pe = e.is_left_slot ? p.right_edge_c : p.left_edge_c;
        const Side ps = e.is_left_slot ? p.right_side : p.left_side;
        *pp = trav_pos(pe, ps);
        *py = p.y;
    };
    const ChordEndpoint& a = first;
    const ChordEndpoint& b = second;

    std::size_t pa = trav_pos(a.edge_c, a.side);
    std::size_t pb = trav_pos(b.edge_c, b.side);
    if (lin_less(pa, a.y, pb, b.y))
        return true;
    if (lin_less(pb, b.y, pa, a.y))
        return false;

    if (a.pending_idx == b.pending_idx)
        return a.is_left_slot && !b.is_left_slot;
    std::size_t qa, qb;
    SymbolicY qya, qyb;
    partner_pos(a, &qa, &qya);
    partner_pos(b, &qb, &qyb);
    if (qa == qb && symbolic_y_equal(qya, qyb))
        return a.is_left_slot && !b.is_left_slot;
    const bool a_closer = lin_less(qa, qya, pa, a.y);
    const bool b_closer = lin_less(qb, qyb, pb, b.y);
    if (a_closer != b_closer)
        return a_closer;
    return lin_less(qb, qyb, qa, qya);
}

void build_submap_from_chords(Submap& out_S, const Polygon& curve,
                              std::vector<PendingChord> pending) {
    for (auto& chord : pending)
        canonicalize_chord(chord, curve);
    std::sort(pending.begin(), pending.end(), chord_less);
    pending.erase(std::unique(pending.begin(), pending.end(),
                              [](const auto& a, const auto& b) {
                                  return !chord_less(a, b) && !chord_less(b, a);
                              }),
                  pending.end());
    std::vector<ChordEndpoint> endpoints;
    endpoints.reserve(2 * pending.size());
    for (std::size_t i = 0; i < pending.size(); ++i) {
        const auto& p = pending[i];
        endpoints.push_back({p.left_edge_c, p.left_side, p.y, i, true});
        endpoints.push_back({p.right_edge_c, p.right_side, p.y, i, false});
    }
    std::sort(endpoints.begin(), endpoints.end(), [&](const auto& a, const auto& b) {
        return chord_endpoint_precedes(curve, pending, a, b);
    });
    build_submap_from_ordered_chords(out_S, curve, pending, endpoints);
}

void build_submap_from_ordered_chords(Submap& out_S, const Polygon& curve,
                                      const std::vector<PendingChord>& pending,
                                      const std::vector<ChordEndpoint>& endpoints) {
    assert(out_S.num_nodes() == 0 && out_S.num_arcs() == 0 && out_S.num_chords() == 0 &&
           "[C91 §3.1 tex 226]: builder requires a fresh output submap");

    using Pending = PendingChord;
    const std::size_t n_edges = curve.num_edges();

#ifdef CHAZELLE_EXPENSIVE_ASSERTS
    for (const Pending& p : pending) {
        if (p.is_null_length)
            continue;
        const Exact xl = edge_x_at_y(curve, p.left_edge_c, p.y);
        const Exact xr = edge_x_at_y(curve, p.right_edge_c, p.y);
        const Side dl = shooting_direction(p.left_edge_c, p.left_side, curve);
        Point pl{xl, p.y.y, p.y.tag};
        RayHit h = naive_first_contact(curve, pl, p.y, dl, p.left_edge_c);
        [[maybe_unused]] bool ok =
            h.hit && h.x == xr && h.side == p.right_side &&
            (h.edge == p.right_edge_c ||
             (symbolic_y_equal(p.y,
                               symbolic_y_of(curve.vertex(std::max(h.edge, p.right_edge_c)))) &&
              !curve.is_y_extremum(std::max(h.edge, p.right_edge_c)) &&
              (h.edge + 1 == p.right_edge_c || p.right_edge_c + 1 == h.edge)));
        assert(ok && "[C91 §2.2 tex 92]: every inventory chord must join a "
                     "mutually visible pair with respect to C");
    }
#endif

    auto trav_pos = [&](std::size_t edge, Side side) -> std::size_t {
        return (side == LEFT) ? edge : (2 * n_edges - 1 - edge);
    };

    assert(endpoints.size() == 2 * pending.size() &&
           "[C91 §2.4 tex 138]: each chord has two boundary endpoints");
#ifndef NDEBUG
    std::vector<bool> visited(endpoints.size(), false);
    for (std::size_t i = 0; i < endpoints.size(); ++i) {
        const auto& e = endpoints[i];
        assert(e.pending_idx < pending.size());
        const auto& p = pending[e.pending_idx];
        const std::size_t slot = 2 * e.pending_idx + (e.is_left_slot ? 0 : 1);
        assert(!visited[slot] &&
               "[C91 §2.4 tex 138]: the boundary traversal visits each endpoint once");
        visited[slot] = true;
        assert(e.edge_c == (e.is_left_slot ? p.left_edge_c : p.right_edge_c) &&
               e.side == (e.is_left_slot ? p.left_side : p.right_side) &&
               symbolic_y_equal(e.y, p.y));
        assert((i == 0 || !chord_endpoint_precedes(curve, pending, e, endpoints[i - 1])) &&
               "[C91 §2.4 tex 138]: endpoints follow the canonical boundary traversal");
    }
#endif

    auto is_vertex_endpoint = [&](std::size_t edge_c, const SymbolicY& y) -> bool {
        const auto& e = curve.edge(edge_c);
        return symbolic_y_equal(y, symbolic_y_of(curve.vertex(e.start_idx))) ||
               symbolic_y_equal(y, symbolic_y_of(curve.vertex(e.end_idx)));
    };

    struct ChordOut {
        std::size_t r_outer = NONE;
        std::size_t r_inner = NONE;
        Chord::AdjArcs left_adj{};
        Chord::AdjArcs right_adj{};
        bool first_seen = false;

        std::size_t pending_after_left = NONE;
        std::size_t pending_after_right = NONE;

        std::size_t r_inner_arc = NONE;
    };
    std::vector<ChordOut> chord_out(pending.size());

    std::vector<std::size_t> awaiting_after_arc;

    std::vector<Arc> arcs;
    arcs.reserve(2 * pending.size() + 2);

    const std::size_t r_start = out_S.add_node();
    std::size_t current_region = r_start;

    struct Cursor {
        std::size_t edge;
        Side side;
    };
    Cursor cursor{0, LEFT};

    auto emit_arc = [&](std::size_t end_edge, Side end_side,
                        std::size_t override_edge_count = NONE) {
        Arc a;
        a.first_edge = cursor.edge;
        a.first_side = cursor.side;
        a.last_edge = end_edge;
        a.last_side = end_side;
        a.region_node = current_region;
        if (override_edge_count != NONE) {
            a.edge_count = override_edge_count;
        } else {
            auto [lo, hi] = a.underlying_edge_range(0, curve.num_vertices() - 1);
            a.edge_count = curve.count_nonnull_edges(lo, hi);
        }
        std::size_t idx = arcs.size();
        arcs.push_back(a);
        cursor = {end_edge, end_side};
        return idx;
    };

    auto patch_after_arcs = [&](std::size_t after_arc) {
        for (std::size_t pi : awaiting_after_arc) {
            ChordOut& co = chord_out[pi];
            if (co.pending_after_left != NONE) {
                assert(co.left_adj.count == 1);
                co.left_adj.arcs[1] = after_arc;
                co.left_adj.count = 2;
                co.pending_after_left = NONE;
            }
            if (co.pending_after_right != NONE) {
                assert(co.right_adj.count == 1);
                co.right_adj.arcs[1] = after_arc;
                co.right_adj.count = 2;
                co.pending_after_right = NONE;
            }
        }
        awaiting_after_arc.clear();
    };

    const std::size_t end_v_edge = n_edges - 1;

    auto group_at_vertex = [&](const ChordEndpoint& e, std::size_t vidx) -> bool {
        return symbolic_y_equal(e.y, symbolic_y_of(curve.vertex(vidx)));
    };
    auto trav_start_vertex = [&](std::size_t edge, Side side) -> std::size_t {
        return (side == LEFT) ? curve.edge(edge).start_idx : curve.edge(edge).end_idx;
    };
    auto trav_end_vertex = [&](std::size_t edge, Side side) -> std::size_t {
        return (side == LEFT) ? curve.edge(edge).end_idx : curve.edge(edge).start_idx;
    };

    auto process_group = [&](std::size_t i, std::size_t j) {
        const ChordEndpoint& g = endpoints[i];

        std::size_t before_arc;
        std::size_t tsv = trav_start_vertex(g.edge_c, g.side);
        if (group_at_vertex(g, tsv)) {
            if (cursor.edge == g.edge_c && cursor.side == g.side) {
                before_arc = emit_arc(g.edge_c, g.side, 0);
            } else if (tsv == 0 || tsv + 1 == curve.num_vertices()) {
                before_arc = emit_arc(g.edge_c, g.side);
            } else {
                std::size_t prev_edge = (g.side == LEFT) ? g.edge_c - 1 : g.edge_c + 1;
                before_arc = emit_arc(prev_edge, g.side);
                cursor = {g.edge_c, g.side};
            }
        } else {
            before_arc = emit_arc(g.edge_c, g.side);
        }
        patch_after_arcs(before_arc);

        bool separator_ready = true;
        for (std::size_t k = i; k < j; ++k) {
            const ChordEndpoint& e = endpoints[k];
            if (!separator_ready) {
                std::size_t sep = emit_arc(e.edge_c, e.side, 0);
                patch_after_arcs(sep);
            }
            separator_ready = false;
            ChordOut& co = chord_out[e.pending_idx];
            bool is_vert = is_vertex_endpoint(e.edge_c, e.y);

            Chord::AdjArcs& adj = e.is_left_slot ? co.left_adj : co.right_adj;
            std::size_t& pending_after =
                e.is_left_slot ? co.pending_after_left : co.pending_after_right;

            if (!co.first_seen) {
                adj.arcs[0] = arcs.size() - 1;
                adj.count = 1;
                co.first_seen = true;
                co.r_outer = current_region;
                co.r_inner = out_S.add_node();
                current_region = co.r_inner;

                if (pending[e.pending_idx].is_null_length) {
                    co.r_inner_arc = emit_arc(e.edge_c, e.side, 0);
                    assert(is_vert && "[C91 §2.1 tex 72]: null-length endpoints are "
                                      "polygon-vertex companions");
                    separator_ready = true;
                    continue;
                }
            } else {
                adj.arcs[0] =
                    pending[e.pending_idx].is_null_length ? co.r_inner_arc : arcs.size() - 1;
                assert(adj.arcs[0] != NONE);
                adj.count = 1;
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
               symbolic_y_equal(endpoints[j].y, endpoints[i].y))
            ++j;
        process_group(i, j);

        {
            const ChordEndpoint& g = endpoints[i];
            if (group_at_vertex(g, trav_end_vertex(g.edge_c, g.side))) {
                if (g.side == LEFT && g.edge_c + 1 <= end_v_edge)
                    cursor = {g.edge_c + 1, LEFT};
                else if (g.side == RIGHT && g.edge_c > 0)
                    cursor = {g.edge_c - 1, RIGHT};
                else if (g.side == RIGHT && g.edge_c == 0)

                    tail_is_zero = true;
            }
        }
        i = j;
    }

    std::size_t tail_piece = emit_arc(0, RIGHT, tail_is_zero ? 0 : NONE);
    patch_after_arcs(tail_piece);

    assert(current_region == r_start &&
           "[C91 §3.1]: parenthesis sweep must return to the starting region");

    std::size_t head_piece = 0;
    bool head_merged = false;
    if (arcs.size() >= 2) {
        Arc& tail = arcs[tail_piece];
        const Arc& head = arcs[head_piece];
        assert(tail.region_node == head.region_node &&
               "[C91 §2.2 tex 96]: the pieces flanking C's start "
               "turnaround bound the same region");
        assert(head.first_edge == 0 && head.first_side == LEFT && tail.last_edge == 0 &&
               tail.last_side == RIGHT &&
               "[C91 §2.4(iii) tex 138]: the sweep starts and ends at "
               "C's start turnaround");
        tail.last_edge = head.last_edge;
        tail.last_side = head.last_side;
        if (tail.edge_count == 0 && head.edge_count == 0) {
        } else {
            auto [lo, hi] = tail.underlying_edge_range(0, curve.num_vertices() - 1);
            tail.edge_count = curve.count_nonnull_edges(lo, hi);
        }
        head_merged = true;
    }

    std::vector<std::size_t> left_order, right_order;
    left_order.reserve(arcs.size());
    right_order.reserve(arcs.size());
    for (std::size_t ai = 0; ai < arcs.size(); ++ai) {
        if (head_merged && ai == head_piece)
            continue;
        (arcs[ai].first_side == LEFT ? left_order : right_order).push_back(ai);
    }
    for (std::size_t k = 1; k < left_order.size(); ++k)
        assert(arcs[left_order[k - 1]].first_edge <= arcs[left_order[k]].first_edge &&
               "[C91 §2.4 tex 138]: the boundary sweep already orders LEFT arcs");
    for (std::size_t k = 1; k < right_order.size(); ++k)
        assert(arcs[right_order[k - 1]].first_edge >= arcs[right_order[k]].first_edge &&
               "[C91 §2.4 tex 138]: the boundary sweep already orders RIGHT arcs");

    out_S.start_vertex = 0;
    out_S.end_vertex = curve.num_vertices() - 1;

    std::vector<std::size_t> arc_remap(arcs.size(), NONE);
    auto add_in_order = [&](const std::vector<std::size_t>& order) {
        for (std::size_t old_ai : order)
            arc_remap[old_ai] = out_S.add_arc(arcs[old_ai]);
    };
    add_in_order(left_order);
    add_in_order(right_order);
    if (head_merged)
        arc_remap[head_piece] = arc_remap[tail_piece];

    out_S.start_arc = arc_remap[tail_piece];
    out_S.end_arc = left_order.empty() ? arc_remap[tail_piece] : arc_remap[left_order.back()];
    assert(out_S.start_arc == out_S.num_arcs() - 1 &&
           "[C91 §2.4(iii) tex 138]: the start-turn arc sorts last");
    assert((out_S.arc(out_S.end_arc).wraps_end() || out_S.num_live_chords() == 0) &&
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
        c.region[0] = co.r_outer;
        c.region[1] = co.r_inner;
        c.left_edge = p.left_edge_c;
        c.left_side = p.left_side;
        c.right_edge = p.right_edge_c;
        c.right_side = p.right_side;
        c.y = p.y.y;
        c.y_tag = p.y.tag;
        c.is_null_length = p.is_null_length;
        c.left_adj = co.left_adj;
        c.right_adj = co.right_adj;
        remap_adj(c.left_adj);
        remap_adj(c.right_adj);
        out_S.add_chord(c);
    }

    out_S.refresh_arc_edge_counts(curve);
}

}
