#include "../polygon/polygon.h"
#include "submap.h"

#include <algorithm>

namespace chazelle::animation {

void Submap::assert_tree_property() const {
#ifndef NDEBUG
    assert(num_live_nodes() >= 1 && "[C91 §2.2]: submap must have at least one region");

    assert(num_live_nodes() == num_live_chords() + 1 &&
           "[C91 §2.2]: submap tree property: num_regions = num_chords + 1");
#endif
}

void Submap::check_invariants() const {
#ifndef NDEBUG
    assert_tree_property();

    for (std::size_t i = 0; i < chords_.size(); ++i) {
        const auto& c = chords_[i];
        if (c.dead)
            continue;
        assert(c.region[0] < nodes_.size() && !nodes_[c.region[0]].dead &&
               "[C91 §2.2]: chord region[0] invalid or dead");
        assert(c.region[1] < nodes_.size() && !nodes_[c.region[1]].dead &&
               "[C91 §2.2]: chord region[1] invalid or dead");
    }

    for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
        if (arc_sequence_[i].dead)
            continue;
        assert(arc_sequence_[i].region_node < nodes_.size() &&
               !nodes_[arc_sequence_[i].region_node].dead &&
               "[C91 §2.2]: arc region_node invalid or dead");
    }

    for (std::size_t i = 0; i < chords_.size(); ++i) {
        const auto& c = chords_[i];
        if (c.dead)
            continue;
        auto check_adj = [&](const Chord::AdjArcs& adj) {
            for (std::size_t j = 0; j < adj.count; ++j) {
                std::size_t ai = adj.arcs[j];
                assert(ai < arc_sequence_.size() && !arc_sequence_[ai].dead &&
                       "[C91 §2.4(ii)]: adj_arc index invalid or dead");
                const auto& a = arc_sequence_[ai];
                assert((a.region_node == c.region[0] || a.region_node == c.region[1]) &&
                       "[C91 §2.4(ii)]: adj_arc must belong to one of "
                       "the chord's endpoint regions");
            }
        };
        check_adj(c.left_adj);
        check_adj(c.right_adj);
    }

    {
        bool seen_right = false;
        std::size_t prev_live = NONE;
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            if (arc_sequence_[i].dead)
                continue;
            if (arc_sequence_[i].first_side == RIGHT)
                seen_right = true;
            if (seen_right) {
                assert(arc_sequence_[i].first_side == RIGHT &&
                       "[C91 §2.4(iii)]: LEFT arc after RIGHT arc in "
                       "arc-sequence table violates ∂C order");
            }
            if (prev_live != NONE) {
                if (arc_sequence_[prev_live].first_side == LEFT &&
                    arc_sequence_[i].first_side == LEFT) {
                    assert(arc_sequence_[i].first_edge >= arc_sequence_[prev_live].first_edge &&
                           "[C91 §2.4(iii)]: LEFT arcs must have ascending "
                           "first_edge for double_identify binary search");
                }
                if (arc_sequence_[prev_live].first_side == RIGHT &&
                    arc_sequence_[i].first_side == RIGHT) {
                    assert(arc_sequence_[i].first_edge <= arc_sequence_[prev_live].first_edge &&
                           "[C91 §2.4(iii)]: RIGHT arcs must have descending "
                           "first_edge for double_identify binary search");
                }
            }
            prev_live = i;
        }
    }

    if (num_live_arcs() > 0) {
        assert(start_arc != NONE && "[C91 §2.4(iii)]: start_arc must be set when arcs exist");
        assert(end_arc != NONE && "[C91 §2.4(iii)]: end_arc must be set when arcs exist");
        assert(start_arc < arc_sequence_.size() && !arc_sequence_[start_arc].dead &&
               "[C91 §2.4(iii)]: start_arc out of range or dead");
        assert(end_arc < arc_sequence_.size() && !arc_sequence_[end_arc].dead &&
               "[C91 §2.4(iii)]: end_arc out of range or dead");
        assert(start_vertex != NONE && "[C91 §2.4(iii)]: start_vertex required when arcs exist");
        assert(end_vertex != NONE && end_vertex > start_vertex &&
               "[C91 §2.4(iii)]: end_vertex required when arcs exist");

        std::size_t last_live = NONE;
        std::size_t last_live_left = NONE;
        std::size_t wrap_count = 0;
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            const Arc& a = arc_sequence_[i];
            if (a.dead)
                continue;
            last_live = i;
            if (a.first_side == LEFT)
                last_live_left = i;
            if (a.wraps())
                ++wrap_count;
        }

        assert(start_arc == last_live && "[C91 §2.4(iii) tex 138]: start_arc must be the last "
                                         "table entry (its cw start position is maximal)");

        if (last_live_left != NONE) {
            assert(end_arc == last_live_left && "[C91 §2.4 tex 144]: end_arc must be the last "
                                                "LEFT-starting arc");
        } else {
            assert(end_arc == last_live && "[C91 §2.4 tex 142]: with no LEFT-starting arc the "
                                           "double-wrap arc covers both turnarounds");
        }

        const Arc& sa = arc_sequence_[start_arc];
        const Arc& ea = arc_sequence_[end_arc];

        bool closed = (num_live_chords() == 0);
        if (closed) {
            assert(num_live_arcs() == 1 && start_arc == end_arc &&
                   "[C91 §2.2 tex 96]: a chordless submap is one region "
                   "bounded by the single closed arc (all of ∂C)");
            assert(sa.first_side == LEFT && sa.last_side == RIGHT &&
                   sa.first_edge == start_vertex && sa.last_edge == start_vertex &&
                   "[C91 §2.4 tex 142/138]: the closed arc is stored cut "
                   "at C's start turnaround");
        } else {
            assert(ea.wraps_end() && "[C91 §2.4 tex 142]: end_arc must double-back around "
                                     "C's end vertex (one arc-structure, never split)");
            assert(sa.wraps_start() && "[C91 §2.4 tex 142]: start_arc must double-back around "
                                       "C's start vertex (one arc-structure, never split)");
        }

        std::size_t expected_wraps = (start_arc == end_arc) ? 1 : 2;
        assert(wrap_count == expected_wraps &&
               "[C91 §2.4 tex 142]: only the arcs passing through C's "
               "endpoints double-back");
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            const Arc& a = arc_sequence_[i];
            if (a.dead || !a.wraps())
                continue;
            assert((i == start_arc || i == end_arc) &&
                   "[C91 §2.4 tex 142]: a wrapping arc must be an "
                   "endpoint arc");
        }
    }
#endif
}

void Submap::check_invariants([[maybe_unused]] const Polygon& polygon) const {
#ifndef NDEBUG
    check_invariants();

    {
        std::size_t i = 0;
        while (i < arc_sequence_.size()) {
            if (arc_sequence_[i].dead) {
                ++i;
                continue;
            }

            std::size_t j = i;
            while (true) {
                std::size_t next = j + 1;
                while (next < arc_sequence_.size() && arc_sequence_[next].dead)
                    ++next;
                if (next >= arc_sequence_.size())
                    break;
                if (arc_sequence_[next].first_side != arc_sequence_[i].first_side ||
                    arc_sequence_[next].first_edge != arc_sequence_[i].first_edge)
                    break;
                j = next;
            }

            if (j > i) {
                const auto& e = polygon.edge(arc_sequence_[i].first_edge);
                bool edge_ascending = symbolic_y_less(symbolic_y_of(polygon.vertex(e.start_idx)),
                                                      symbolic_y_of(polygon.vertex(e.end_idx)));
                bool asc = (arc_sequence_[i].first_side == LEFT) ? edge_ascending : !edge_ascending;
                std::size_t prev_live = i;
                for (std::size_t k = i + 1; k <= j; ++k) {
                    if (arc_sequence_[k].dead)
                        continue;
                    if (asc) {
                        assert(symbolic_y_leq(arc_start_symbolic_y(prev_live, polygon),
                                              arc_start_symbolic_y(k, polygon)) &&
                               "[C91 §2.4 tex 144]: same-first_edge run must "
                               "be start-y-monotonic (ascending)");
                    } else {
                        assert(symbolic_y_geq(arc_start_symbolic_y(prev_live, polygon),
                                              arc_start_symbolic_y(k, polygon)) &&
                               "[C91 §2.4 tex 144]: same-first_edge run must "
                               "be start-y-monotonic (descending)");
                    }
                    prev_live = k;
                }
            }
            i = j + 1;
        }
    }

    if (num_live_arcs() > 0) {
        assert(start_vertex != NONE && end_vertex != NONE && end_vertex > 0 &&
               "[C91 §2.4(iii)]: start/end_vertex must be set when arcs exist");
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            const auto& a = arc_sequence_[i];
            if (a.dead)
                continue;

            assert(arc_start_symbolic_y(i, polygon).tag != SOS_NONE &&
                   "[C91 §2 tex 47]: live arc's start position needs SoS tag");
            if (a.edge_count == 0)
                continue;

            assert(a.first_edge < polygon.num_edges() && a.last_edge < polygon.num_edges() &&
                   "[C91 §2.4(iii)]: arc edges must be valid input-table indices");
            std::size_t actual = arc_boundary_edge_count(a, polygon, start_vertex, end_vertex,
                                                         arc_start_symbolic_y(i, polygon),
                                                         arc_end_symbolic_y(i, polygon));
            assert(a.edge_count == actual && "[C91 §2.2 tex 106]: arc.edge_count cache must match "
                                             "the arc's on each side nonnull ∂C edge count "
                                             "(arc_boundary_edge_count)");
        }
    }

    auto matches_an_endpoint = [&](std::size_t edge_idx, const SymbolicY& chord_y) -> bool {
        assert(edge_idx < polygon.num_edges());
        const auto& e = polygon.edge(edge_idx);
        return symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.start_idx))) ||
               symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.end_idx)));
    };
    for (std::size_t ci = 0; ci < chords_.size(); ++ci) {
        const auto& c = chords_[ci];
        if (c.dead)
            continue;

        SymbolicY chord_y{c.y, c.y_tag};

        if (c.is_null_length) {
            const std::size_t v = polygon.local_index_of_tag(c.y_tag);
            assert(v != NONE && "[C91 §2.1 tex 72]: null chord level names its extremum");
            assert(symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(v))) &&
                   "[C91 §2.1 tex 72]: null chord lies at its extremum's y");
            assert(polygon.is_y_extremum(v) &&
                   "[C91 §2.1 tex 72]: null-length chord source must be an "
                   "interior local y-extremum of C");
            assert((c.left_edge == v || (v >= 1 && c.left_edge == v - 1)) &&
                   "[C91 §2.1 tex 72]: null-length chord endpoints sit at "
                   "the extremum vertex (edge v-1 or v)");
            continue;
        }

        assert((c.left_adj.count == 1) == matches_an_endpoint(c.left_edge, chord_y) &&
               "[C91 §2.2 tex 94]: LEFT endpoint count == 1 ⟺ endpoint is a polygon vertex");
        assert((c.right_adj.count == 1) == matches_an_endpoint(c.right_edge, chord_y) &&
               "[C91 §2.2 tex 94]: RIGHT endpoint count == 1 ⟺ endpoint is a polygon vertex");

        auto arc_spans_vertex = [&](std::size_t arc_idx, std::size_t edge_idx) -> bool {
            const auto& a = arc_sequence_[arc_idx];
            auto [elo, ehi] = a.underlying_edge_range(start_vertex, end_vertex);
            const auto e = polygon.edge(edge_idx);
            const std::size_t v =
                symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.start_idx))) ? e.start_idx
                                                                                      : e.end_idx;
            return v >= elo && v <= ehi + 1;
        };
        if (c.left_adj.count == 1) {
            assert(
                arc_spans_vertex(c.left_adj.arcs[0], c.left_edge) &&
                "[C91 §2.2 tex 94]: LEFT vertex endpoint's adj arc must span the polygon vertex");
        }
        if (c.right_adj.count == 1) {
            assert(
                arc_spans_vertex(c.right_adj.arcs[0], c.right_edge) &&
                "[C91 §2.2 tex 94]: RIGHT vertex endpoint's adj arc must span the polygon vertex");
        }
    }
#endif
}

}
