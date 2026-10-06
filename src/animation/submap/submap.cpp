#include "submap.h"
#include "../polygon/polygon.h"
#include "../trace.h"

#include <algorithm>
#include <utility>

namespace chazelle::animation {

std::size_t SubmapNode::degree() const noexcept {
    return incident_chords.size();
}

std::size_t Submap::add_node() {
    std::size_t idx = nodes_.size();
    nodes_.push_back(SubmapNode{});
    tree_decomp_dirty_ = true;
    return idx;
}

std::size_t Submap::add_arc(Arc arc) {
    assert(arc.region_node != NONE && arc.region_node < nodes_.size() &&
           !nodes_[arc.region_node].dead && "[C91 §2.4(ii) tex 137]: arc.region_node must be live");

    std::size_t idx = arc_sequence_.size();

    if (!arc_sequence_.empty()) {
        Side prev_side = arc_sequence_.back().first_side;
        if (arc.first_side == LEFT) {
            assert(prev_side == LEFT && "[C91 §2.4(iii)]: LEFT arc after RIGHT violates ∂C order");
            assert(arc.first_edge >= arc_sequence_.back().first_edge &&
                   "[C91 §2.4(iii)]: LEFT arcs must ascend in first_edge");
        } else if (prev_side == RIGHT) {
            assert(arc.first_edge <= arc_sequence_.back().first_edge &&
                   "[C91 §2.4(iii)]: RIGHT arcs must descend in first_edge");
        }
    }

    if (arc.first_side == LEFT)
        left_right_boundary_ = idx + 1;
    arc_sequence_.push_back(arc);

    if (arc.wraps_end())
        end_arc = idx;
    bool closed_encoding = arc.first_side == LEFT && arc.last_side == RIGHT &&
                           arc.first_edge == arc.last_edge && start_vertex != NONE &&
                           arc.first_edge == start_vertex;
    if (arc.wraps_start() || closed_encoding)
        start_arc = idx;

    tree_decomp_dirty_ = true;
    return idx;
}

std::size_t Submap::add_chord(Chord chord) {
    std::size_t idx = chords_.size();

    assert(chord.left_adj.count >= 1 && chord.left_adj.count <= 2 &&
           "[C91 §2.4(ii)]: LEFT endpoint adj count ∈ [1,2]");
    assert(chord.right_adj.count >= 1 && chord.right_adj.count <= 2 &&
           "[C91 §2.4(ii)]: RIGHT endpoint adj count ∈ [1,2]");
    assert(chord.left_adj.count + chord.right_adj.count >= 2 &&
           chord.left_adj.count + chord.right_adj.count <= 4 &&
           "[C91 §2.4(ii) tex 137]: total adj count ∈ [2,4]");
    auto check_adj = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            assert(adj.arcs[k] != NONE && adj.arcs[k] < arc_sequence_.size() &&
                   "[C91 §2.4(ii)]: adj_arc index must be valid");
            assert(!arc_sequence_[adj.arcs[k]].dead && "[C91 §2.4(ii)]: adj_arc must be live");
            assert((arc_sequence_[adj.arcs[k]].region_node == chord.region[0] ||
                    arc_sequence_[adj.arcs[k]].region_node == chord.region[1]) &&
                   "[C91 §2.4(ii)]: adj_arc must belong to one of the "
                   "chord's two regions");
        }

        if (adj.count == 2) {
            assert(arc_sequence_[adj.arcs[0]].region_node !=
                       arc_sequence_[adj.arcs[1]].region_node &&
                   "[C91 §2.2 tex 96]: mid-edge endpoint's two adj arcs "
                   "must lie in the chord's two distinct regions");
        }
    };
    check_adj(chord.left_adj);
    check_adj(chord.right_adj);

    assert(chord.region[0] != chord.region[1] &&
           "[C91 §2.2 tex 102]: chord must connect two distinct regions");

    assert(chord.y_tag != SOS_NONE &&
           "[C91 §2 tex 47 (SoS)]: chord must carry its symbolic level tag");

    if (chord.is_null_length) {
        assert(chord.left_edge == chord.right_edge && chord.left_side == chord.right_side &&
               "[C91 §2.1 tex 72]: null-length chord endpoints must coincide");
        assert(chord.left_adj.count == 1 && chord.right_adj.count == 1 &&
               "[C91 §2.1 tex 72]: null-length chord has 1 adj arc per ∂C side");
    }

    chords_.push_back(std::move(chord));
    ++live_chords_;

    for (std::size_t r : chords_[idx].region) {
        assert(r != NONE && r < nodes_.size() && !nodes_[r].dead &&
               "[C91 §2.4(i)]: chord must connect two LIVE regions");
        nodes_[r].incident_chords.push_back(idx);
    }

    tree_decomp_dirty_ = true;
    return idx;
}

bool Submap::arc_is_point(std::size_t arc_idx, const Polygon& polygon) const {
    assert(arc_idx < arc_sequence_.size() && !arc_sequence_[arc_idx].dead);
    const Arc& a = arc_sequence_[arc_idx];

    return a.edge_count == 0 && a.first_side == a.last_side && a.first_edge == a.last_edge &&
           symbolic_y_equal(arc_start_symbolic_y(arc_idx, polygon),
                            arc_end_symbolic_y(arc_idx, polygon));
}

std::size_t Submap::find_junction_arc(const Chord& c, bool query_left, std::size_t edge, Side side,
                                      std::size_t vertex_idx, bool want_after, std::size_t exclude,
                                      std::size_t exclude2, const Polygon& polygon) const {
    assert(vertex_idx != NONE && vertex_idx < polygon.num_vertices() &&
           "[C91 §2.2 tex 96]: junction lookup requires a polygon vertex");
    assert(
        edge < polygon.num_edges() &&
        (vertex_idx == polygon.edge(edge).start_idx || vertex_idx == polygon.edge(edge).end_idx) &&
        "[C91 §2.2 tex 96]: junction vertex must bound the endpoint edge");

    SymbolicY jy = c.symbolic_y();

    std::size_t expected_region = NONE;
    if (!c.is_null_length) {
        const Chord::AdjArcs& other = query_left ? c.right_adj : c.left_adj;
        std::size_t other_before;
        if (other.count == 1) {
            other_before = other.arcs[0];
        } else {
            const std::size_t oe = query_left ? c.right_edge : c.left_edge;
            const Side os = query_left ? c.right_side : c.left_side;
            auto starts_there = [&](std::size_t ai) {
                const Arc& a = arc_sequence_[ai];
                return a.first_edge == oe && a.first_side == os &&
                       symbolic_y_equal(arc_start_symbolic_y(ai, polygon), jy);
            };
            other_before = starts_there(other.arcs[0]) ? other.arcs[1] : other.arcs[0];
        }
        const std::size_t rb = arc_sequence_[other_before].region_node;
        assert((rb == c.region[0] || rb == c.region[1]) &&
               "[C91 §2.4(ii)]: slot arcs bound the chord's regions");
        expected_region = want_after ? rb : (c.region[0] == rb ? c.region[1] : c.region[0]);
    }

    std::size_t found_zero = NONE;
    std::size_t found_nonzero = NONE;

    auto consider = [&](std::size_t ai) {
        if (ai == NONE || ai == exclude || ai == exclude2)
            return;
        assert(ai < arc_sequence_.size());
        const Arc& a = arc_sequence_[ai];
        if (a.dead)
            return;
        if (a.region_node != c.region[0] && a.region_node != c.region[1])
            return;
        if (expected_region != NONE && a.region_node != expected_region)
            return;
        Side a_side = want_after ? a.first_side : a.last_side;
        std::size_t a_edge = want_after ? a.first_edge : a.last_edge;
        if (a_edge + 1 < vertex_idx || a_edge > vertex_idx)
            return;
        if (a_side != side) {
            if (!polygon.is_y_extremum(vertex_idx))
                return;
            const std::size_t other_e = (edge == vertex_idx) ? vertex_idx - 1 : vertex_idx;
            if (a_edge != other_e)
                return;
            if (is_inside_companion(polygon, a_edge, a_side, vertex_idx) !=
                is_inside_companion(polygon, edge, side, vertex_idx))
                return;
        }
        SymbolicY ay =
            want_after ? arc_start_symbolic_y(ai, polygon) : arc_end_symbolic_y(ai, polygon);
        if (!symbolic_y_equal(ay, jy))
            return;
        bool zero = (a.edge_count == 0);
        std::size_t& slot = zero ? found_zero : found_nonzero;
        assert((slot == NONE || slot == ai) && "[C91 §2.2 tex 96]: at most one mate per junction "
                                               "(distinct arcs cannot share a ∂C start/end point)");
        slot = ai;
    };

    for (std::size_t r : c.region) {
        for (std::size_t ci : nodes_[r].incident_chords) {
            const Chord& ch = chords_[ci];
            if (ch.dead)
                continue;
            for (std::size_t k = 0; k < ch.left_adj.count; ++k)
                consider(ch.left_adj.arcs[k]);
            for (std::size_t k = 0; k < ch.right_adj.count; ++k)
                consider(ch.right_adj.arcs[k]);
        }
    }

    std::size_t result = (found_zero != NONE) ? found_zero : found_nonzero;
    assert(result != NONE && "[C91 §2.2 tex 96]: an interior junction always has both a "
                             "before- and an after-arc");
    return result;
}

std::size_t Submap::remove_chord(std::size_t chord_idx, const Polygon& polygon) {
    assert(chord_idx < chords_.size());
    auto& c = chords_[chord_idx];
    assert(!c.dead && "removing an already-dead chord");
    if (auto* trace = AnimationTrace::current())
        trace->remove(*this, polygon, chord_idx, c);
    assert(c.region[0] != NONE && c.region[1] != NONE);

    std::size_t r0 = c.region[0];
    std::size_t r1 = c.region[1];
    assert(!nodes_[r0].dead && !nodes_[r1].dead);

    if (num_live_chords() == 1) {
        assert(start_vertex != NONE && end_vertex != NONE && end_vertex > start_vertex &&
               "[C91 §2.4(iii) tex 138]: C endpoints must be identified");

        std::size_t keep = NONE;
        auto sweep_slots = [&](const Chord::AdjArcs& adj) {
            for (std::size_t k = 0; k < adj.count; ++k) {
                std::size_t ai = adj.arcs[k];
                assert(ai != NONE && ai < arc_sequence_.size() &&
                       "[C91 §2.4(ii)]: adj_arc must be valid");
                if (keep == NONE) {
                    assert(!arc_sequence_[ai].dead && "[C91 §2.4(ii)]: adj_arc must be live");
                    keep = ai;
                } else if (ai != keep && !arc_sequence_[ai].dead) {
                    arc_sequence_[ai].dead = true;
                    if (auto* trace = AnimationTrace::current())
                        trace->delete_arc(*this, ai);
                    compacted_ = false;
                }
            }
        };
        sweep_slots(c.left_adj);
        sweep_slots(c.right_adj);
        assert(keep != NONE && "[C91 §2.4(ii)]: a chord references its regions' arcs");

        Arc& a = arc_sequence_[keep];
        a.first_edge = start_vertex;
        a.first_side = LEFT;
        a.last_edge = start_vertex;
        a.last_side = RIGHT;

        a.edge_count = 2 * polygon.count_nonnull_edges(start_vertex, end_vertex - 1);
        a.region_node = r0;
        start_arc = keep;
        end_arc = keep;

        c.dead = true;
        --live_chords_;
        nodes_[r1].dead = true;
        nodes_[r0].incident_chords.clear();
        tree_decomp_dirty_ = true;
        if (auto* trace = AnimationTrace::current()) {
            trace->build_arc(*this, polygon, keep, a, symbolic_y_of(polygon.vertex(0)),
                             symbolic_y_of(polygon.vertex(0)));
            trace->settled("contract_end", *this);
        }
        return r0;
    }

    auto endpoint_vertex = [&](std::size_t edge, const Exact& ey,
                               std::size_t ey_tag) -> std::size_t {
        assert(edge < polygon.num_edges() && "[C91 §2.2]: invalid edge index");
        const auto& e = polygon.edge(edge);
        SymbolicY chord_y{ey, ey_tag};
        if (symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.start_idx))))
            return e.start_idx;
        if (symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.end_idx))))
            return e.end_idx;
        return NONE;
    };

    std::size_t left_vertex = endpoint_vertex(c.left_edge, c.y, c.y_tag);
    std::size_t right_vertex = endpoint_vertex(c.right_edge, c.y, c.y_tag);
    [[maybe_unused]] bool left_is_vertex = left_vertex != NONE;
    [[maybe_unused]] bool right_is_vertex = right_vertex != NONE;

    struct GluedArcEndpoints {
        std::size_t arc;
        SymbolicY start, end;
    };
    GluedArcEndpoints glued[2];
    std::size_t n_glued = 0;
    auto start_of = [&](std::size_t ai) -> SymbolicY {
        for (std::size_t k = 0; k < n_glued; ++k)
            if (glued[k].arc == ai)
                return glued[k].start;
        return arc_start_symbolic_y(ai, polygon);
    };
    auto end_of = [&](std::size_t ai) -> SymbolicY {
        for (std::size_t k = 0; k < n_glued; ++k)
            if (glued[k].arc == ai)
                return glued[k].end;
        return arc_end_symbolic_y(ai, polygon);
    };
    auto record_glued = [&](std::size_t ai, const SymbolicY& sy, const SymbolicY& ey) {
        for (std::size_t k = 0; k < n_glued; ++k)
            if (glued[k].arc == ai) {
                glued[k].start = sy;
                glued[k].end = ey;
                return;
            }
        assert(n_glued < 2 && "[C91 §2.2 tex 94]: at most 2 glues");
        glued[n_glued++] = {ai, sy, ey};
    };

    auto glue_arcs = [&](std::size_t ai, std::size_t aj) {
        assert(ai != NONE && ai < arc_sequence_.size() && aj != NONE && aj < arc_sequence_.size() &&
               !arc_sequence_[ai].dead && !arc_sequence_[aj].dead &&
               "[C91 §2.4(ii)]: adj_arc must be valid + live");

        assert(ai != aj && "[C91 §2.2 tex 94]: glued arc pair must be distinct "
                           "(the last-chord closure path handles ai == aj)");

        auto& a_keep = arc_sequence_[ai];
        auto& a_dead = arc_sequence_[aj];

        assert(a_keep.last_side == a_dead.first_side &&
               "[C91 §2.2 tex 96]: glue mates share the junction's ∂C side");
        assert((a_keep.last_edge == a_dead.first_edge ||
                (a_keep.last_side == LEFT && a_keep.last_edge + 1 == a_dead.first_edge) ||
                (a_keep.last_side == RIGHT && a_keep.last_edge == a_dead.first_edge + 1)) &&
               "[C91 §2.2 tex 96]: glue mates' edges must coincide or be "
               "traversal-consecutive at the junction");

        auto is_point_here = [&](std::size_t idx) {
            const Arc& a = arc_sequence_[idx];
            return a.edge_count == 0 && a.first_side == a.last_side &&
                   a.first_edge == a.last_edge && symbolic_y_equal(start_of(idx), end_of(idx));
        };
        const bool dead_is_point = is_point_here(aj);
        const bool keep_is_point = is_point_here(ai);
        SymbolicY merged_start_y{}, merged_end_y{};
        if (dead_is_point) {
            merged_start_y = start_of(ai);
            merged_end_y = end_of(ai);
        } else if (keep_is_point) {
            merged_start_y = start_of(aj);
            merged_end_y = end_of(aj);
            a_keep.first_edge = a_dead.first_edge;
            a_keep.first_side = a_dead.first_side;
            a_keep.last_edge = a_dead.last_edge;
            a_keep.last_side = a_dead.last_side;
            a_keep.edge_count = a_dead.edge_count;
        } else {
            merged_start_y = start_of(ai);
            merged_end_y = end_of(aj);
            a_keep.last_edge = a_dead.last_edge;
            a_keep.last_side = a_dead.last_side;
            a_keep.edge_count = arc_boundary_edge_count(a_keep, polygon, start_vertex, end_vertex,
                                                        merged_start_y, merged_end_y);
        }
        record_glued(ai, merged_start_y, merged_end_y);

        a_dead.dead = true;
        if (auto* trace = AnimationTrace::current())
            trace->delete_arc(*this, aj);
        compacted_ = false;

        if (start_arc == aj)
            start_arc = ai;
        if (end_arc == aj)
            end_arc = ai;

        auto replace_arc = [&](Chord::AdjArcs& adj) {
            for (std::size_t k = 0; k < adj.count; ++k)
                if (adj.arcs[k] == aj)
                    adj.arcs[k] = ai;
        };
        for (std::size_t ri : {r0, r1}) {
            for (std::size_t ci : nodes_[ri].incident_chords) {
                auto& other = chords_[ci];
                if (other.dead)
                    continue;
                replace_arc(other.left_adj);
                replace_arc(other.right_adj);
            }
        }
    };

    if (c.is_null_length) {
        assert(c.left_adj.count == 1 && c.right_adj.count == 1 &&
               "[C91 §2.2]: null-length chord must have 1 adj arc per side");

        assert(left_is_vertex && right_is_vertex &&
               "[C91 §2.1 tex 72]: null-length chord endpoints are polygon vertices");
        assert(left_vertex == right_vertex &&
               "[C91 §2.1 tex 72]: null-length chord endpoints share the vertex");
        assert(left_vertex != start_vertex && left_vertex != end_vertex &&
               "[C91 §2.1 tex 72]: null-length chords arise only at "
               "non-endpoint local extrema");

        std::size_t sl = c.left_adj.arcs[0];
        std::size_t sr = c.right_adj.arcs[0];
        bool sl_inner = arc_sequence_[sl].region_node == c.region[1];
        [[maybe_unused]] bool sr_inner = arc_sequence_[sr].region_node == c.region[1];
        assert(sl_inner != sr_inner &&
               "[C91 §2.1 tex 72]: exactly one slot holds the inner null arc");
        std::size_t inner = sl_inner ? sl : sr;
        std::size_t outer = sl_inner ? sr : sl;
        assert(arc_sequence_[inner].edge_count == 0 &&
               "[C91 §2.1 tex 72]: the inner region of a null-length chord "
               "is bounded by a null arc");

        const Arc& oa = arc_sequence_[outer];
        bool outer_is_before = (oa.last_side == c.left_side) &&
                               (oa.last_edge == c.left_edge || oa.last_edge + 1 == c.left_edge ||
                                oa.last_edge == c.left_edge + 1) &&
                               symbolic_y_equal(arc_end_symbolic_y(outer, polygon), c.symbolic_y());
        if (outer_is_before) {
            glue_arcs(outer, inner);
            std::size_t after = find_junction_arc(c, true, c.left_edge, c.left_side, left_vertex,
                                                  true, outer, inner, polygon);
            glue_arcs(outer, after);
        } else {
            std::size_t before = find_junction_arc(c, true, c.left_edge, c.left_side, left_vertex,
                                                   false, outer, inner, polygon);
            glue_arcs(before, inner);
            glue_arcs(before, outer);
        }
    } else {
        auto glue_endpoint = [&](const Chord::AdjArcs& adj, bool ql, std::size_t edge, Side side,
                                 std::size_t vtx) {
            if (vtx == NONE) {
                assert(adj.count == 2 && "[C91 §2.2 tex 94]: non-vertex endpoint needs 2 adj arcs");
                glue_arcs(adj.arcs[0], adj.arcs[1]);
                return;
            }
            assert(adj.count == 1 && "[C91 §2.2 tex 94]: vertex endpoint records one adj arc");
            std::size_t before = adj.arcs[0];
            std::size_t after =
                find_junction_arc(c, ql, edge, side, vtx, true, before, NONE, polygon);
            glue_arcs(before, after);
        };
        glue_endpoint(c.left_adj, true, c.left_edge, c.left_side, left_vertex);
        glue_endpoint(c.right_adj, false, c.right_edge, c.right_side, right_vertex);
    }

    auto reassign_live = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            std::size_t ai = adj.arcs[k];
            assert(ai != NONE && ai < arc_sequence_.size() && !arc_sequence_[ai].dead &&
                   "[C91 §2.4(ii)]: adj arcs must be live after glueing");
            if (arc_sequence_[ai].region_node == r1)
                arc_sequence_[ai].region_node = r0;
        }
    };
    for (std::size_t ci : nodes_[r1].incident_chords) {
        if (ci == chord_idx)
            continue;
        const auto& ch = chords_[ci];
        if (ch.dead)
            continue;
        reassign_live(ch.left_adj);
        reassign_live(ch.right_adj);
    }
    reassign_live(c.left_adj);
    reassign_live(c.right_adj);

    for (std::size_t ci : nodes_[r1].incident_chords) {
        if (ci == chord_idx)
            continue;
        auto& ch = chords_[ci];
        if (ch.dead)
            continue;
        nodes_[r0].incident_chords.push_back(ci);
        if (ch.region[0] == r1)
            ch.region[0] = r0;
        if (ch.region[1] == r1)
            ch.region[1] = r0;
    }

    {
        auto& ic = nodes_[r0].incident_chords;
        ic.erase(std::remove(ic.begin(), ic.end(), chord_idx), ic.end());
    }

    c.dead = true;
    --live_chords_;
    nodes_[r1].dead = true;
    tree_decomp_dirty_ = true;
    if (auto* trace = AnimationTrace::current()) {
        for (std::size_t i = 0; i < n_glued; ++i) {
            const auto& endpoint = glued[i];
            if (arc_sequence_[endpoint.arc].dead)
                continue;
            trace->build_arc(*this, polygon, endpoint.arc, arc_sequence_[endpoint.arc],
                             endpoint.start, endpoint.end);
        }
        for (std::size_t ci : nodes_[r0].incident_chords)
            trace->build_chord(*this, polygon, ci, chords_[ci]);
        trace->settled("contract_end", *this);
    }
    return r0;
}

std::size_t Submap::num_live_nodes() const noexcept {
    std::size_t n = 0;
    for (const auto& nd : nodes_)
        if (!nd.dead)
            ++n;
    return n;
}

std::size_t Submap::num_live_chords() const noexcept {
#ifdef CHAZELLE_EXPENSIVE_ASSERTS
    std::size_t n = 0;
    for (const auto& ch : chords_)
        if (!ch.dead)
            ++n;
    assert(n == live_chords_ && "live chord counter must match the table");
#endif
    return live_chords_;
}

std::size_t Submap::num_live_arcs() const noexcept {
    std::size_t n = 0;
    for (const auto& a : arc_sequence_)
        if (!a.dead)
            ++n;
    return n;
}

void Submap::refresh_arc_edge_counts(const Polygon& polygon) {
    assert(start_vertex != NONE && end_vertex != NONE &&
           "[C91 §2.4(iii)]: C endpoints must be identified");
    for (std::size_t ai = 0; ai < arc_sequence_.size(); ++ai) {
        Arc& a = arc_sequence_[ai];
        if (a.dead || a.edge_count == 0)
            continue;
        a.edge_count = arc_boundary_edge_count(a, polygon, start_vertex, end_vertex,
                                               arc_start_symbolic_y(ai, polygon),
                                               arc_end_symbolic_y(ai, polygon));
    }
}

Submap::InsertChordResult Submap::insert_chord(const ChordPointSpec& p, const ChordPointSpec& q,
                                               SymbolicY y, std::size_t region,
                                               const std::size_t* cycle, std::size_t cycle_len,
                                               const Polygon& polygon) {
    assert(region < nodes_.size() && !nodes_[region].dead &&
           "[C91 §3.2]: insert_chord requires a live region");
    assert(y.tag != SOS_NONE && "[C91 §2 tex 47 (SoS)]: chord must carry its source SoS tag");
    assert(p.arc != q.arc && "[C91 §3.2 Lemma 3.3]: the chord connects two distinct "
                             "(nonconsecutive) arcs of the region");
    assert((p.edge != q.edge || p.side != q.side) &&
           "[C91 §2.1 tex 70]: a chord connects two distinct ∂C points. "
           "Equal x IS legal: the outside duplicate pair at a y-extremum "
           "([C91 §2.1 tex 72]) sees itself through infinity ([C91 §2.1 "
           "tex 70]), giving a zero-geometric-length wrap chord between "
           "distinct ∂C points (chord_runs_through_infinity)");
    assert(cycle_len >= 2 && "[C91 §3.2]: region cycle must hold both arcs");

    std::size_t ip = NONE, iq = NONE;
    for (std::size_t i = 0; i < cycle_len; ++i) {
        std::size_t ai = cycle[i];
        assert(ai < arc_sequence_.size() && !arc_sequence_[ai].dead &&
               arc_sequence_[ai].region_node == region &&
               "[C91 §2.2 tex 96]: cycle entries must be live arcs of the region");
        if (ai == p.arc)
            ip = i;
        if (ai == q.arc)
            iq = i;
    }
    assert(ip != NONE && iq != NONE &&
           "[C91 §3.2]: p.arc and q.arc must appear in the region cycle");

#ifndef NDEBUG
    auto assert_not_chord_endpoint = [&](std::size_t edge, Side side) {
        for (std::size_t ci : nodes_[region].incident_chords) {
            const Chord& c = chords_[ci];
            if (c.dead)
                continue;
            if (!symbolic_y_equal(c.symbolic_y(), y))
                continue;
            assert(!((c.left_edge == edge && c.left_side == side) ||
                     (c.right_edge == edge && c.right_side == side)) &&
                   "[C91 §2.1 tex 70]: new chord endpoint coincides with an "
                   "existing chord endpoint (visibility already realized)");
        }
    };
    assert_not_chord_endpoint(p.edge, p.side);
    assert_not_chord_endpoint(q.edge, q.side);
#endif

    auto split_arc = [&](const ChordPointSpec& sp) -> std::size_t {
        SymbolicY start_y = arc_start_symbolic_y(sp.arc, polygon);
        SymbolicY end_y = arc_end_symbolic_y(sp.arc, polygon);

        Arc before = arc_sequence_[sp.arc];
        assert(before.covers(sp.edge, sp.side, start_vertex, end_vertex) &&
               "[C91 §3.2]: split point must lie on the arc "
               "([C91 §2.4 tex 142]: wrap arcs cover on each side)");

        bool at_start = (sp.edge == before.first_edge && sp.side == before.first_side &&
                         symbolic_y_equal(y, start_y));
        bool at_end = (sp.edge == before.last_edge && sp.side == before.last_side &&
                       symbolic_y_equal(y, end_y));

        assert(!(at_start && at_end) && "[C91 §3.2]: cannot split a zero-length arc");

        Arc after = before;
        after.first_edge = sp.edge;
        after.first_side = sp.side;
        before.last_edge = sp.edge;
        before.last_side = sp.side;

        if (!at_start && !at_end) {
            const auto& pe = polygon.edge(sp.edge);
            SymbolicY y_lo = symbolic_y_of(polygon.vertex(pe.start_idx));
            SymbolicY y_hi = symbolic_y_of(polygon.vertex(pe.end_idx));
            bool at_trav_end =
                (sp.side == LEFT) ? symbolic_y_equal(y, y_hi) : symbolic_y_equal(y, y_lo);
            bool at_trav_start =
                (sp.side == LEFT) ? symbolic_y_equal(y, y_lo) : symbolic_y_equal(y, y_hi);
            assert(!(at_trav_start && at_trav_end) &&
                   "[C91 §2 tex 47 (SoS)]: an edge's endpoints have "
                   "distinct symbolic ys");
            if (at_trav_end) {
                std::size_t trav_end_v = (sp.side == LEFT) ? pe.end_idx : pe.start_idx;
                if (trav_end_v != start_vertex && trav_end_v != end_vertex) {
                    after.first_edge = (sp.side == LEFT) ? sp.edge + 1 : sp.edge - 1;
                    assert(after.first_edge < polygon.num_edges() &&
                           "[C91 §3.2]: vertex split's after-half must have "
                           "a following edge (its span is nonempty)");
                }
            } else if (at_trav_start) {
                std::size_t trav_start_v = (sp.side == LEFT) ? pe.start_idx : pe.end_idx;
                if (trav_start_v != start_vertex && trav_start_v != end_vertex) {
                    before.last_edge = (sp.side == LEFT) ? sp.edge - 1 : sp.edge + 1;
                    assert(before.last_edge < polygon.num_edges() &&
                           "[C91 §3.2]: vertex split's before-half must have "
                           "a preceding edge (its span is nonempty)");
                }
            }
        }

        if (at_start) {
            before.edge_count = 0;
        } else {
            before.edge_count =
                arc_boundary_edge_count(before, polygon, start_vertex, end_vertex, start_y, y);
        }
        if (at_end) {
            after.edge_count = 0;
        } else {
            after.edge_count =
                arc_boundary_edge_count(after, polygon, start_vertex, end_vertex, y, end_y);
        }

        arc_sequence_[sp.arc] = before;
        std::size_t after_idx = arc_sequence_.size();
        arc_sequence_.push_back(after);
        if (auto* trace = AnimationTrace::current())
            trace->split_arc(*this, polygon, sp.arc, after_idx, start_y, y, end_y);
        return after_idx;
    };
    std::size_t p_after = split_arc(p);
    std::size_t q_after = split_arc(q);

    std::size_t r_new = add_node();
    arc_sequence_[p_after].region_node = r_new;
    if (auto* trace = AnimationTrace::current())
        trace->arc_owner(*this, p_after, r_new);
    for (std::size_t i = (ip + 1) % cycle_len;; i = (i + 1) % cycle_len) {
        assert(i != ip && "[C91 §3.2]: chain walk must reach q before p");
        arc_sequence_[cycle[i]].region_node = r_new;
        if (auto* trace = AnimationTrace::current())
            trace->arc_owner(*this, cycle[i], r_new);
        if (i == iq)
            break;
    }

    auto repoint = [&](Chord::AdjArcs& adj, std::size_t old_arc, std::size_t after_idx) {
        assert(adj.count >= 1 && adj.count <= 2 &&
               "[C91 §2.4(ii) tex 137]: each chord endpoint has one or two adjacent arcs");
        for (std::size_t k = 0; k < adj.count; ++k) {
            if (adj.arcs[k] != old_arc)
                continue;
            if (!(adj.count == 2 && k == 1))
                adj.arcs[k] = after_idx;
        }
    };
    for (std::size_t ci : nodes_[region].incident_chords) {
        Chord& c = chords_[ci];
        assert(!c.dead && "[C91 §2.4]: incident_chords holds live chords only");
        repoint(c.left_adj, p.arc, p_after);
        repoint(c.left_adj, q.arc, q_after);
        repoint(c.right_adj, p.arc, p_after);
        repoint(c.right_adj, q.arc, q_after);
    }

    {
        std::vector<std::size_t> stay;
        stay.reserve(nodes_[region].incident_chords.size());
        for (std::size_t ci : nodes_[region].incident_chords) {
            Chord& c = chords_[ci];
            std::size_t side_region = NONE;
            auto scan = [&](const Chord::AdjArcs& adj) {
                for (std::size_t k = 0; k < adj.count; ++k) {
                    std::size_t rn = arc_sequence_[adj.arcs[k]].region_node;
                    if (rn != region && rn != r_new)
                        continue;
                    assert((side_region == NONE || side_region == rn) &&
                           "[C91 §3.2]: a chord's boundary footprint must "
                           "lie entirely in one chain (pq crosses no chord)");
                    side_region = rn;
                }
            };
            scan(c.left_adj);
            scan(c.right_adj);
            assert(side_region != NONE &&
                   "[C91 §2.4(ii)]: chord must reference an arc of its region");
            if (side_region == r_new) {
                if (c.region[0] == region)
                    c.region[0] = r_new;
                else {
                    assert(c.region[1] == region);
                    c.region[1] = r_new;
                }
                nodes_[r_new].incident_chords.push_back(ci);
            } else {
                stay.push_back(ci);
            }
        }
        nodes_[region].incident_chords = std::move(stay);
    }

    auto endpoint_is_vertex = [&](std::size_t edge) -> bool {
        assert(edge < polygon.num_edges());
        const auto& e = polygon.edge(edge);
        return symbolic_y_equal(y, symbolic_y_of(polygon.vertex(e.start_idx))) ||
               symbolic_y_equal(y, symbolic_y_of(polygon.vertex(e.end_idx)));
    };
    auto make_adj = [&](const ChordPointSpec& sp, std::size_t after_idx) -> Chord::AdjArcs {
        Chord::AdjArcs adj;
        adj.arcs[0] = sp.arc;
        if (endpoint_is_vertex(sp.edge)) {
            adj.count = 1;
        } else {
            adj.arcs[1] = after_idx;
            adj.count = 2;
        }
        return adj;
    };

    Chord nc;
    nc.region[0] = region;
    nc.region[1] = r_new;
    nc.y = y.y;
    nc.y_tag = y.tag;
    nc.is_null_length = false;

    auto tie_pos = [&](std::size_t edge, Side side) {
        return (side == LEFT) ? edge : (2 * polygon.num_edges() - 1 - edge);
    };
    const bool p_left =
        p.x < q.x || (p.x == q.x && tie_pos(p.edge, p.side) < tie_pos(q.edge, q.side));
    const ChordPointSpec& lp = p_left ? p : q;
    const ChordPointSpec& rp = p_left ? q : p;
    std::size_t l_after = p_left ? p_after : q_after;
    std::size_t r_after = p_left ? q_after : p_after;
    nc.left_edge = lp.edge;
    nc.left_side = lp.side;
    nc.right_edge = rp.edge;
    nc.right_side = rp.side;
    nc.left_adj = make_adj(lp, l_after);
    nc.right_adj = make_adj(rp, r_after);
    std::size_t chord_idx = add_chord(nc);
    if (auto* trace = AnimationTrace::current())
        trace->insert(*this, polygon, chord_idx, nc);

    auto repoint_wrap = [&](const ChordPointSpec& sp, std::size_t after_idx) {
        if (end_arc == sp.arc && !arc_sequence_[sp.arc].wraps_end() &&
            arc_sequence_[after_idx].wraps_end())
            end_arc = after_idx;
        if (start_arc == sp.arc && !arc_sequence_[sp.arc].wraps_start() &&
            arc_sequence_[after_idx].wraps_start())
            start_arc = after_idx;
    };
    repoint_wrap(p, p_after);
    repoint_wrap(q, q_after);

    compacted_ = false;
    tree_decomp_dirty_ = true;

    if (auto* trace = AnimationTrace::current()) {
        for (std::size_t ci : nodes_[region].incident_chords)
            trace->build_chord(*this, polygon, ci, chords_[ci]);
        for (std::size_t ci : nodes_[r_new].incident_chords)
            trace->build_chord(*this, polygon, ci, chords_[ci]);
        trace->settled("split_end", *this);
    }
    return InsertChordResult{chord_idx, r_new, p_after, q_after};
}

void Submap::build_tree_decomposition() {
    tree_decomp_.build(*this);
    tree_decomp_dirty_ = false;
}

}
