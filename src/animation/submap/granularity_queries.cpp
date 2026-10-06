#include "../polygon/polygon.h"
#include "submap.h"

#include <algorithm>

namespace chazelle::animation {

std::size_t Submap::region_weight(std::size_t node_idx) const noexcept {
    assert(node_idx < nodes_.size() && !nodes_[node_idx].dead);

    const auto& nd = nodes_[node_idx];
    std::size_t max_count = 0;

    auto check_adj = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            std::size_t ai = adj.arcs[k];
            assert(ai != NONE && ai < arc_sequence_.size() && !arc_sequence_[ai].dead &&
                   "[C91 §2.4(ii)]: adj_arc must be valid + live");
            const auto& a = arc_sequence_[ai];
            if (a.region_node == node_idx && a.edge_count > max_count)
                max_count = a.edge_count;
        }
    };
    for (std::size_t ci : nd.incident_chords) {
        assert(ci < chords_.size());
        const auto& ch = chords_[ci];
        if (ch.dead)
            continue;
        check_adj(ch.left_adj);
        check_adj(ch.right_adj);
    }

    if (nd.incident_chords.empty()) {
        assert(num_live_chords() == 0 && "[C91 §2.2 tex 102]: a chord-free region exists only in "
                                         "the chordless (single-region) submap");
        assert(start_arc != NONE && start_arc == end_arc && start_arc < arc_sequence_.size() &&
               !arc_sequence_[start_arc].dead &&
               "[C91 §2.4(iii) tex 138]: the chordless submap's closed "
               "arc is the endpoint arc");
        const auto& a = arc_sequence_[start_arc];
        if (a.region_node == node_idx && a.edge_count > max_count)
            max_count = a.edge_count;
    }

    return max_count;
}

bool Submap::is_conformal() const noexcept {
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i].dead)
            continue;
        if (nodes_[i].degree() > 4)
            return false;
    }
    return true;
}

bool Submap::is_semigranular(std::size_t granularity) const noexcept {
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i].dead)
            continue;
        if (region_weight(i) > granularity)
            return false;
    }
    return true;
}

std::size_t Submap::simulated_contraction_weight(std::size_t chord_idx,
                                                 const Polygon& polygon) const noexcept {
    assert(chord_idx < chords_.size() && !chords_[chord_idx].dead);
    const auto& c = chords_[chord_idx];
    return simulated_contraction_weight(chord_idx, polygon, region_weight(c.region[0]),
                                        region_weight(c.region[1]));
}

std::size_t Submap::simulated_contraction_weight(std::size_t chord_idx, const Polygon& polygon,
                                                 std::size_t w0, std::size_t w1) const noexcept {
    assert(chord_idx < chords_.size() && !chords_[chord_idx].dead);
    const auto& c = chords_[chord_idx];
    assert(c.region[0] != NONE && c.region[1] != NONE && c.region[0] != c.region[1] &&
           "[C91 §2.4(i)/§2.2 tex 102]: chord regions valid + distinct (tree)");

    std::size_t max_count = std::max(w0, w1);

    if (num_live_chords() == 1) {
        assert(start_vertex != NONE && end_vertex != NONE && end_vertex > start_vertex &&
               "[C91 §2.4(iii) tex 138]: C endpoints must be identified");

        std::size_t merged = 2 * polygon.count_nonnull_edges(start_vertex, end_vertex - 1);
        return std::max(max_count, merged);
    }

    auto endpoint_vertex = [&](std::size_t edge) -> std::size_t {
        assert(edge < polygon.num_edges() && "[C91 §2.2]: invalid edge");
        const auto& e = polygon.edge(edge);
        SymbolicY chord_y{c.y, c.y_tag};
        if (symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.start_idx))))
            return e.start_idx;
        if (symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.end_idx))))
            return e.end_idx;
        return NONE;
    };

    struct AdjacentArcPair {
        std::size_t before, after;
    };
    AdjacentArcPair pairs[2];
    std::size_t np = 0;

    if (c.is_null_length) {
        std::size_t sl = c.left_adj.arcs[0];
        std::size_t sr = c.right_adj.arcs[0];
        bool sl_inner = arc_sequence_[sl].region_node == c.region[1];
        [[maybe_unused]] bool sr_inner = arc_sequence_[sr].region_node == c.region[1];
        assert(sl_inner != sr_inner &&
               "[C91 §2.1 tex 72]: exactly one slot holds the inner null arc");
        std::size_t inner = sl_inner ? sl : sr;
        std::size_t outer = sl_inner ? sr : sl;
        std::size_t v = endpoint_vertex(c.left_edge);
        assert(v != NONE && "[C91 §2.1 tex 72]: null-length chord endpoints are polygon vertices");
        const Arc& oa = arc_sequence_[outer];
        bool outer_is_before = (oa.last_side == c.left_side) &&
                               (oa.last_edge == c.left_edge || oa.last_edge + 1 == c.left_edge ||
                                oa.last_edge == c.left_edge + 1) &&
                               symbolic_y_equal(arc_end_symbolic_y(outer, polygon), c.symbolic_y());
        std::size_t missing = find_junction_arc(c, true, c.left_edge, c.left_side, v,
                                                outer_is_before, outer, inner, polygon);
        std::size_t before = outer_is_before ? outer : missing;
        std::size_t after = outer_is_before ? missing : outer;
        pairs[np++] = {before, inner};
        pairs[np++] = {inner, after};
    } else {
        auto collect = [&](const Chord::AdjArcs& adj, bool ql, std::size_t edge, Side side) {
            std::size_t v = endpoint_vertex(edge);
            if (v == NONE) {
                assert(adj.count == 2 && "[C91 §2.2 tex 94]: non-vertex endpoint needs 2 adj arcs");
                pairs[np++] = {adj.arcs[0], adj.arcs[1]};
                return;
            }
            assert(adj.count == 1 && "[C91 §2.2 tex 94]: vertex endpoint records one adj arc");
            std::size_t after =
                find_junction_arc(c, ql, edge, side, v, true, adj.arcs[0], NONE, polygon);
            pairs[np++] = {adj.arcs[0], after};
        };
        collect(c.left_adj, true, c.left_edge, c.left_side);
        collect(c.right_adj, false, c.right_edge, c.right_side);
    }
    assert(np == 2 && "[C91 §2.2 tex 94]: chord removal glues arcs at its two endpoints");

    auto chain_count = [&](const std::size_t* members, std::size_t n) {
        for (std::size_t i = 0; i < n; ++i)
            assert(members[i] < arc_sequence_.size() && !arc_sequence_[members[i]].dead &&
                   "[C91 §2.4(ii)]: glue mates must be valid + live");

        Arc acc = arc_sequence_[members[0]];
        bool acc_point = arc_is_point(members[0], polygon);
        for (std::size_t i = 1; i < n; ++i) {
            const Arc& d = arc_sequence_[members[i]];
            if (arc_is_point(members[i], polygon)) {
            } else if (acc_point) {
                acc = d;
                acc_point = false;
            } else {
                acc.last_edge = d.last_edge;
                acc.last_side = d.last_side;
            }
        }
        if (acc_point)
            return;
        std::size_t merged = arc_boundary_edge_count(acc, polygon, start_vertex, end_vertex,
                                                     arc_start_symbolic_y(members[0], polygon),
                                                     arc_end_symbolic_y(members[n - 1], polygon));
        if (merged > max_count)
            max_count = merged;
    };

    if (np == 2) {
        bool shared_lr = (pairs[0].after == pairs[1].before);
        bool shared_rl = (pairs[1].after == pairs[0].before);

        assert(pairs[0].before != pairs[1].before && pairs[0].after != pairs[1].after &&
               "[C91 §2.2 tex 94]: an arc cannot take the same junction "
               "role at both chord endpoints");

        assert(!(shared_lr && shared_rl) &&
               "[C91 §2.2 tex 102]: doubly-shared glue pairs occur only "
               "for the last chord (closure early-return)");
        if (shared_lr) {
            std::size_t chain[3] = {pairs[0].before, pairs[0].after, pairs[1].after};
            chain_count(chain, 3);
            return max_count;
        }
        if (shared_rl) {
            std::size_t chain[3] = {pairs[1].before, pairs[1].after, pairs[0].after};
            chain_count(chain, 3);
            return max_count;
        }
    }
    for (std::size_t i = 0; i < np; ++i) {
        std::size_t chain[2] = {pairs[i].before, pairs[i].after};
        chain_count(chain, 2);
    }

    return max_count;
}

bool Submap::is_granular(std::size_t granularity, const Polygon& polygon) const noexcept {
    if (!is_semigranular(granularity))
        return false;

    if (num_live_chords() == 0)
        return true;

    for (std::size_t ci = 0; ci < chords_.size(); ++ci) {
        const auto& c = chords_[ci];
        if (c.dead)
            continue;
        std::size_t d0 = nodes_[c.region[0]].degree();
        std::size_t d1 = nodes_[c.region[1]].degree();
        if (d0 >= 3 && d1 >= 3)
            continue;
        if (simulated_contraction_weight(ci, polygon) <= granularity)
            return false;
    }

#ifndef NDEBUG
    if (is_conformal()) {
        assert(start_vertex != NONE && end_vertex != NONE &&
               "[C91 §2.4(iii) tex 138]: submap with chords needs start/end_vertex");
        std::size_t n_c = end_vertex - start_vertex + 1;
        std::size_t bound = 2 * (8 * (n_c - 1) / (granularity + 1));
        assert(num_live_nodes() <= bound && "[C91 §2.3 Lemma 2.3]: V ≤ 2·⌊8(|C|−1)/(γ+1)⌋");
    }
#endif

    return true;
}

}
