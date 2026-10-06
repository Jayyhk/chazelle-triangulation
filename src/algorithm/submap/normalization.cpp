#include "../polygon/polygon.h"
#include "submap.h"

#include <algorithm>
#include <utility>

namespace chazelle {

void Submap::compact() {
    for (const auto& a : arc_sequence_) {
        if (a.dead)
            continue;
        assert(a.region_node < nodes_.size() && !nodes_[a.region_node].dead &&
               "[C91 §2.2 tex 96]: live arc must point to a live region "
               "(remove_chord's slot walk is exhaustive)");
    }

    std::vector<std::size_t> region_remap(nodes_.size(), NONE);
    {
        std::size_t j = 0;
        for (std::size_t i = 0; i < nodes_.size(); ++i) {
            if (!nodes_[i].dead) {
                region_remap[i] = j;
                if (j != i)
                    nodes_[j] = std::move(nodes_[i]);
                ++j;
            }
        }
        nodes_.resize(j);
    }

    std::vector<std::size_t> chord_remap(chords_.size(), NONE);
    {
        std::size_t j = 0;
        for (std::size_t i = 0; i < chords_.size(); ++i) {
            if (!chords_[i].dead) {
                chord_remap[i] = j;
                if (j != i)
                    chords_[j] = std::move(chords_[i]);
                ++j;
            }
        }
        chords_.resize(j);
    }

    std::vector<std::size_t> arc_remap(arc_sequence_.size(), NONE);
    {
        std::size_t j = 0;
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            if (!arc_sequence_[i].dead) {
                arc_remap[i] = j;
                if (j != i)
                    arc_sequence_[j] = arc_sequence_[i];
                ++j;
            }
        }
        arc_sequence_.resize(j);
    }

    for (auto& a : arc_sequence_) {
        assert(a.region_node != NONE && region_remap[a.region_node] != NONE);
        a.region_node = region_remap[a.region_node];
    }

    auto remap_adjacent_arcs = [&](Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            if (adj.arcs[k] != NONE) {
                assert(arc_remap[adj.arcs[k]] != NONE);
                adj.arcs[k] = arc_remap[adj.arcs[k]];
            }
        }
    };
    for (auto& ch : chords_) {
        for (auto& r : ch.region) {
            assert(r != NONE && region_remap[r] != NONE);
            r = region_remap[r];
        }
        remap_adjacent_arcs(ch.left_adj);
        remap_adjacent_arcs(ch.right_adj);
    }

    for (auto& nd : nodes_) {
        std::size_t retained_chord_count = 0;
        for (std::size_t ci : nd.incident_chords) {
            if (chord_remap[ci] != NONE)
                nd.incident_chords[retained_chord_count++] = chord_remap[ci];
        }
        nd.incident_chords.resize(retained_chord_count);
    }

    if (start_arc != NONE) {
        assert(arc_remap[start_arc] != NONE);
        start_arc = arc_remap[start_arc];
    }
    if (end_arc != NONE) {
        assert(arc_remap[end_arc] != NONE);
        end_arc = arc_remap[end_arc];
    }

    left_right_boundary_ = 0;
    for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
        if (arc_sequence_[i].first_side == LEFT)
            left_right_boundary_ = i + 1;
    }

    compacted_ = true;
    tree_decomp_dirty_ = true;
}

void Submap::normalize(const Polygon& polygon) {
    assert(is_conformal() && "[C91 §2.4(iv)]: normal form requires a conformal submap");

    compact();

    const std::size_t M = arc_sequence_.size();
    if (M > 0) {
        const std::size_t n_edges = polygon.num_edges();
        auto edge_ascends = [&](std::size_t e) -> bool {
            const auto& ed = polygon.edge(e);
            return symbolic_y_less(symbolic_y_of(polygon.vertex(ed.start_idx)),
                                   symbolic_y_of(polygon.vertex(ed.end_idx)));
        };

        struct Key {
            std::size_t trav;
            SymbolicY y;
            bool asc;
            bool zero;
            std::size_t idx;
        };
        std::vector<Key> keys(M);
        for (std::size_t i = 0; i < M; ++i) {
            const Arc& a = arc_sequence_[i];
            assert(!a.dead && "compact() leaves no dead arcs");
            std::size_t trav =
                (a.first_side == LEFT) ? a.first_edge : 2 * n_edges - 1 - a.first_edge;
            bool asc = (a.first_side == LEFT) == edge_ascends(a.first_edge);
            keys[i] = Key{trav, arc_start_symbolic_y(i, polygon), asc, a.edge_count == 0, i};
        }
        std::sort(keys.begin(), keys.end(), [](const Key& a, const Key& b) {
            if (a.trav != b.trav)
                return a.trav < b.trav;
            if (!symbolic_y_equal(a.y, b.y))
                return a.asc ? symbolic_y_less(a.y, b.y) : symbolic_y_greater(a.y, b.y);
            return a.zero && !b.zero;
        });

        for (std::size_t i = 0; i + 1 < M; ++i) {
            assert((keys[i].trav != keys[i + 1].trav ||
                    !symbolic_y_equal(keys[i].y, keys[i + 1].y) ||
                    (keys[i].zero && !keys[i + 1].zero)) &&
                   "[C91 §2.4(iii)]: distinct arcs cannot share a ∂C "
                   "start point and length class");
        }

        std::vector<std::size_t> arc_remap(M, NONE);
        {
            std::vector<Arc> sorted;
            sorted.reserve(M);
            for (std::size_t p = 0; p < M; ++p) {
                arc_remap[keys[p].idx] = p;
                sorted.push_back(arc_sequence_[keys[p].idx]);
            }
            arc_sequence_ = std::move(sorted);
        }
        auto remap_adjacent_arcs = [&](Chord::AdjArcs& adj) {
            for (std::size_t k = 0; k < adj.count; ++k) {
                assert(adj.arcs[k] != NONE && arc_remap[adj.arcs[k]] != NONE);
                adj.arcs[k] = arc_remap[adj.arcs[k]];
            }
        };
        for (auto& ch : chords_) {
            remap_adjacent_arcs(ch.left_adj);
            remap_adjacent_arcs(ch.right_adj);
        }
        if (start_arc != NONE)
            start_arc = arc_remap[start_arc];
        if (end_arc != NONE)
            end_arc = arc_remap[end_arc];

        left_right_boundary_ = 0;
        for (std::size_t i = 0; i < M; ++i)
            if (arc_sequence_[i].first_side == LEFT)
                left_right_boundary_ = i + 1;
    }

    check_invariants(polygon);

    build_tree_decomposition();
}

}
