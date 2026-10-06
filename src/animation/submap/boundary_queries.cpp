#include "../polygon/polygon.h"
#include "submap.h"

#include <algorithm>

namespace chazelle::animation {

namespace {

bool arc_region_above_chord(const Submap& submap, const Polygon& curve, const Chord& ch,
                            std::size_t edge_at, Side side_at, std::size_t arc_idx) {
    const Arc& a = submap.arc(arc_idx);
    bool first_matches = (a.first_edge == edge_at && a.first_side == side_at);
    bool last_matches = (a.last_edge == edge_at && a.last_side == side_at);
    bool starts_at_chord;
    if (first_matches && !last_matches) {
        starts_at_chord = true;
    } else if (!first_matches && last_matches) {
        starts_at_chord = false;
    } else {
        assert(first_matches && last_matches &&
               "adj arc must touch the chord endpoint via first or last");
        starts_at_chord =
            symbolic_y_equal(submap.arc_start_symbolic_y(arc_idx, curve), ch.symbolic_y());
    }
    std::size_t adj_v;
    if (starts_at_chord) {
        adj_v = (a.first_side == LEFT) ? curve.edge(a.first_edge).end_idx
                                       : curve.edge(a.first_edge).start_idx;
    } else {
        adj_v = (a.last_side == LEFT) ? curve.edge(a.last_edge).start_idx
                                      : curve.edge(a.last_edge).end_idx;
    }
    const SymbolicY vy = symbolic_y_of(curve.vertex(adj_v));
    if (symbolic_y_equal(vy, ch.symbolic_y())) {
        const SymbolicY maxy = symbolic_y_of(curve.vertex(curve.max_y_vertex()));
        if (symbolic_y_equal(ch.symbolic_y(), maxy))
            return true;
        assert(
            symbolic_y_equal(ch.symbolic_y(), symbolic_y_of(curve.vertex(curve.min_y_vertex()))) &&
            "[C91 §2.1 tex 70]: a polar-cap chord sits at C's "
            "global y-extremum");
        return false;
    }
    return symbolic_y_greater(vy, ch.symbolic_y());
}

}

void Submap::chord_regions_below_above(std::size_t ci, const Polygon& curve, std::size_t* below,
                                       std::size_t* above) const {
    const Submap& submap = *this;
    const Chord& c = submap.chord(ci);
    const Chord::AdjArcs* adj = nullptr;
    std::size_t edge_at = NONE;
    Side side_at = LEFT;
    if (c.left_adj.count == 2) {
        adj = &c.left_adj;
        edge_at = c.left_edge;
        side_at = c.left_side;
    } else if (c.right_adj.count == 2) {
        adj = &c.right_adj;
        edge_at = c.right_edge;
        side_at = c.right_side;
    }
    if (adj) {
        bool a0 = arc_region_above_chord(submap, curve, c, edge_at, side_at, adj->arcs[0]);
        [[maybe_unused]] bool a1 =
            arc_region_above_chord(submap, curve, c, edge_at, side_at, adj->arcs[1]);
        assert(a0 != a1 && "[C91 §2.2 tex 96]: the two adj arcs at a mid-edge chord "
                           "endpoint lie on opposite sides of the chord");
        std::size_t r0 = submap.arc(adj->arcs[0]).region_node;
        std::size_t r1 = submap.arc(adj->arcs[1]).region_node;
        *above = a0 ? r0 : r1;
        *below = a0 ? r1 : r0;
        return;
    }

    bool left_above =
        arc_region_above_chord(submap, curve, c, c.left_edge, c.left_side, c.left_adj.arcs[0]);
    std::size_t r_arc = submap.arc(c.left_adj.arcs[0]).region_node;
    std::size_t r_other = (c.region[0] == r_arc) ? c.region[1] : c.region[0];
    *above = left_above ? r_arc : r_other;
    *below = left_above ? r_other : r_arc;
}

Submap::DoubleIdentifyResult Submap::double_identify(std::size_t edge_idx, SymbolicY y,
                                                     const Polygon& polygon) const {
    DoubleIdentifyResult result;
    assert(!arc_sequence_.empty() && "[C91 §2.4 tex 144]: double identification requires a submap");
    assert(start_vertex != NONE && end_vertex != NONE && edge_idx >= start_vertex &&
           edge_idx < end_vertex && edge_idx < polygon.num_edges() &&
           "[C91 §2.4 tex 144]: e is an edge of C");
    [[maybe_unused]] Exact query_x = 0.0;
    assert(edge_crossing_x(polygon, edge_idx, y, &query_x) &&
           "[C91 §2.4 tex 144]: q lies on the supplied edge e");

    assert(compacted_ && "[C91 §2.4]: double_identify requires compacted arc-sequence");

    assert(y.tag != SOS_NONE && "[C91 §2.4]: double_identify requires a valid SoS y-tag");

    assert(start_arc == arc_sequence_.size() - 1 &&
           "[C91 §2.4(iii) tex 138]: start_arc must be the last table entry");
    assert((left_right_boundary_ == 0 || end_arc == left_right_boundary_ - 1) &&
           "[C91 §2.4 tex 144]: end_arc must be the last LEFT-starting arc");

    std::size_t left_begin = 0;
    std::size_t left_end = left_right_boundary_;
    std::size_t right_begin = left_right_boundary_;
    std::size_t right_end = arc_sequence_.size();

    auto search_half = [&](std::size_t lo, std::size_t hi, bool ascending, std::size_t seam) {
        const Side half_side = ascending ? LEFT : RIGHT;

        std::size_t blo = lo, bhi = hi;
        while (blo < bhi) {
            std::size_t mid = blo + (bhi - blo) / 2;
            bool advance = ascending ? (arc_sequence_[mid].first_edge < edge_idx)
                                     : (arc_sequence_[mid].first_edge > edge_idx);
            if (advance)
                blo = mid + 1;
            else
                bhi = mid;
        }

        std::size_t bend = blo;
        if (bend < hi && arc_sequence_[bend].first_edge == edge_idx) {
            std::size_t slo = blo, shi = hi;
            while (slo < shi) {
                std::size_t mid = slo + (shi - slo) / 2;
                if (arc_sequence_[mid].first_edge == edge_idx)
                    slo = mid + 1;
                else
                    shi = mid;
            }
            bend = slo;
        }

        auto arc_contains_edge = [&](std::size_t ai) -> bool {
            assert(ai < arc_sequence_.size() && "[C91 §2.4]: invalid arc index");
            const auto& a = arc_sequence_[ai];
            assert(!a.dead && "[C91 §2.4]: arc_contains_edge on dead arc");
            assert(!a.wraps() || ai == start_arc || ai == end_arc);
            return a.covers(edge_idx, half_side, start_vertex, end_vertex);
        };

        std::size_t boundary_arc = NONE;
        if (blo > lo) {
            if (arc_contains_edge(blo - 1))
                boundary_arc = blo - 1;
        } else if (seam != NONE && !arc_sequence_[seam].dead && !(seam >= blo && seam < bend) &&
                   arc_contains_edge(seam)) {
            boundary_arc = seam;
        }

        std::size_t interval_len = bend - blo;

        if (interval_len == 0 && boundary_arc == NONE)
            return;

        if (interval_len == 0) {
            result.push(boundary_arc);
            return;
        }
        if (interval_len == 1 && boundary_arc == NONE) {
            result.push(blo);
            return;
        }

        if (interval_len == 1 && boundary_arc != NONE) {
            SymbolicY junction_y = arc_start_symbolic_y(blo, polygon);
            if (symbolic_y_equal(junction_y, y)) {
                result.push(blo);
                result.push(boundary_arc);
            } else {
                assert(edge_idx < polygon.num_edges());
                const auto& e = polygon.edge(edge_idx);
                bool edge_ascending = symbolic_y_less(symbolic_y_of(polygon.vertex(e.start_idx)),
                                                      symbolic_y_of(polygon.vertex(e.end_idx)));
                bool traversal_ascending = ascending ? edge_ascending : !edge_ascending;
                bool y_in_boundary = traversal_ascending ? symbolic_y_less(y, junction_y)
                                                         : symbolic_y_greater(y, junction_y);
                result.push(y_in_boundary ? boundary_arc : blo);
            }
            return;
        }

        bool keys_ascending;
        {
            assert(edge_idx < polygon.num_edges());
            const auto& e = polygon.edge(edge_idx);
            bool edge_ascending = symbolic_y_less(symbolic_y_of(polygon.vertex(e.start_idx)),
                                                  symbolic_y_of(polygon.vertex(e.end_idx)));
            keys_ascending = ascending ? edge_ascending : !edge_ascending;
        }

        std::size_t ylo = blo, yhi = bend;
        while (ylo < yhi) {
            std::size_t mid = ylo + (yhi - ylo) / 2;
            SymbolicY mid_y = arc_start_symbolic_y(mid, polygon);
            if (keys_ascending) {
                if (symbolic_y_leq(mid_y, y))
                    ylo = mid + 1;
                else
                    yhi = mid;
            } else {
                if (symbolic_y_geq(mid_y, y))
                    ylo = mid + 1;
                else
                    yhi = mid;
            }
        }

        if (ylo > blo) {
            std::size_t p = ylo - 1;
            result.push(p);
            if (symbolic_y_equal(arc_start_symbolic_y(p, polygon), y) && p > blo) {
                result.push(p - 1);
                for (std::size_t i = p - 1; i > blo; --i) {
                    if (!symbolic_y_equal(arc_start_symbolic_y(i, polygon), y))
                        break;
                    result.push(i - 1);
                }
            }
        } else {
            assert(boundary_arc != NONE && "[C91 §2.4 tex 144]: query precedes the run's first "
                                           "start-y with no boundary arc — off-edge query or "
                                           "corrupt arc-sequence table");
        }

        if (boundary_arc != NONE) {
            const Arc& ba = arc_sequence_[boundary_arc];

            bool starts_here = ba.first_edge == edge_idx && ba.first_side == half_side;
            if (starts_here && symbolic_y_equal(arc_start_symbolic_y(boundary_arc, polygon), y)) {
                result.push(boundary_arc);
            } else {
                SymbolicY first_y = arc_start_symbolic_y(blo, polygon);
                bool in_boundary =
                    keys_ascending ? symbolic_y_leq(y, first_y) : symbolic_y_geq(y, first_y);
                if (in_boundary)
                    result.push(boundary_arc);
            }
        }
    };

    search_half(left_begin, left_end, true, arc_sequence_.size() - 1);
    assert(result.count <= 3 && "[C91 §2.4 tex 144]: ≤ 3 arcs per ∂C half");

    [[maybe_unused]] std::size_t left_count = result.count;
    search_half(right_begin, right_end, false, end_arc);
    assert(result.count - left_count <= 3 && "[C91 §2.4 tex 144]: ≤ 3 arcs per ∂C half");
    assert(result.count <= DoubleIdentifyResult::MAX &&
           "[C91 §2.4 tex 144]: ≤ 6 arcs at any point");

    return result;
}

SymbolicY Submap::arc_start_symbolic_y(std::size_t arc_idx, const Polygon& polygon) const {
    assert(arc_idx < arc_sequence_.size() && !arc_sequence_[arc_idx].dead &&
           "[C91 §2.4]: arc index must be valid + live");
    const Arc& a = arc_sequence_[arc_idx];

    for (std::size_t ci : nodes_[a.region_node].incident_chords) {
        const Chord& c = chords_[ci];
        if (c.dead)
            continue;
        if (c.is_null_length) {
            if (c.right_adj.count == 1 && c.right_adj.arcs[0] == arc_idx &&
                a.region_node == c.region[1])
                return c.symbolic_y();
            if (c.left_adj.count == 1 && c.left_adj.arcs[0] == arc_idx &&
                a.region_node == c.region[1])
                return c.symbolic_y();
            continue;
        }

        if (c.left_adj.count == 2 && c.left_adj.arcs[1] == arc_idx)
            return c.symbolic_y();
        if (c.right_adj.count == 2 && c.right_adj.arcs[1] == arc_idx)
            return c.symbolic_y();
    }

    if (start_vertex != NONE && end_vertex != NONE) {
        auto companion_chord_at = [&](std::size_t edge, Side side, const SymbolicY& vy) -> bool {
            for (std::size_t ci : nodes_[a.region_node].incident_chords) {
                const Chord& c = chords_[ci];
                if (c.dead)
                    continue;
                if (!symbolic_y_equal(c.symbolic_y(), vy))
                    continue;
                if ((c.left_edge == edge && c.left_side == side) ||
                    (c.right_edge == edge && c.right_side == side))
                    return true;
            }
            return false;
        };
        if (arc_idx == start_arc && a.first_side == RIGHT && a.first_edge == start_vertex) {
            SymbolicY vy = symbolic_y_of(polygon.vertex(start_vertex));
            if (companion_chord_at(a.first_edge, RIGHT, vy))
                return vy;
        }
        if (arc_idx == end_arc && a.first_side == LEFT && a.first_edge + 1 == end_vertex) {
            SymbolicY vy = symbolic_y_of(polygon.vertex(end_vertex));
            if (companion_chord_at(a.first_edge, LEFT, vy))
                return vy;
        }
    }

    {
        const std::size_t vend = (a.first_side == LEFT) ? a.first_edge + 1 : a.first_edge;
        const SymbolicY vy = symbolic_y_of(polygon.vertex(vend));
        const std::size_t other_e = (a.first_edge == vend) ? vend - 1 : vend;
        auto corner_label_matches = [&](std::size_t ce, Side cs) {
            if (ce == a.first_edge && cs == a.first_side)
                return true;
            if (ce != other_e || other_e >= polygon.num_edges())
                return false;
            if (polygon.is_y_extremum(vend))
                return is_inside_companion(polygon, ce, cs, vend) ==
                       is_inside_companion(polygon, a.first_edge, a.first_side, vend);
            return cs == a.first_side;
        };
        for (std::size_t ci : nodes_[a.region_node].incident_chords) {
            const Chord& c = chords_[ci];
            if (c.dead || c.is_null_length)
                continue;
            if (!symbolic_y_equal(c.symbolic_y(), vy))
                continue;
            const Chord::AdjArcs* adj = nullptr;
            if (corner_label_matches(c.left_edge, c.left_side))
                adj = &c.left_adj;
            else if (corner_label_matches(c.right_edge, c.right_side))
                adj = &c.right_adj;
            if (!adj || adj->count != 1)
                continue;
            if (adj->arcs[0] == arc_idx)
                continue;
            return c.symbolic_y();
        }
    }

    assert(a.first_edge < polygon.num_edges() &&
           "[C91 §2.4(iii)]: first_edge must be a valid edge index");
    std::size_t vidx = (a.first_side == LEFT) ? a.first_edge : a.first_edge + 1;
    assert(vidx < polygon.num_vertices() &&
           "[C91 §2.4(iii)]: arc-start polygon vertex must be valid");
    return symbolic_y_of(polygon.vertex(vidx));
}

SymbolicY Submap::arc_end_symbolic_y(std::size_t arc_idx, const Polygon& polygon) const {
    assert(arc_idx < arc_sequence_.size() && !arc_sequence_[arc_idx].dead &&
           "[C91 §2.4]: arc index must be valid + live");
    const Arc& a = arc_sequence_[arc_idx];

    for (std::size_t ci : nodes_[a.region_node].incident_chords) {
        const Chord& c = chords_[ci];
        if (c.dead)
            continue;
        if (c.is_null_length) {
            if (c.right_adj.count == 1 && c.right_adj.arcs[0] == arc_idx &&
                a.region_node == c.region[1])
                return c.symbolic_y();
            if (c.left_adj.count == 1 && c.left_adj.arcs[0] == arc_idx &&
                a.region_node == c.region[1])
                return c.symbolic_y();
            continue;
        }

        if (c.left_adj.arcs[0] == arc_idx)
            return c.symbolic_y();
        if (c.right_adj.arcs[0] == arc_idx)
            return c.symbolic_y();
    }

    assert(a.last_edge < polygon.num_edges() &&
           "[C91 §2.4(iii)]: last_edge must be a valid edge index");
    std::size_t vidx = (a.last_side == LEFT) ? a.last_edge + 1 : a.last_edge;
    assert(vidx < polygon.num_vertices() &&
           "[C91 §2.4(iii)]: arc-end polygon vertex must be valid");
    return symbolic_y_of(polygon.vertex(vidx));
}

}
