#include "bounded_regions.h"
#include "../merge/oracle.h"
#include "../submap/boundary_geometry.h"
#include "../submap/chord_inventory.h"
#include "../trace.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <iterator>
#include <utility>
#include <vector>

namespace chazelle::animation {

namespace {

struct VertexOccurrence {
    std::size_t vertex;
    std::size_t edge;
    Side side;
};

struct HeightRun {
    std::size_t first;
    std::size_t count;
    bool increasing;

    std::size_t next() const {
        return increasing ? first : first + count - 1;
    }

    void advance() {
        assert(count > 0);
        --count;
        if (increasing)
            ++first;
    }
};

struct VertexChords {
    std::array<std::size_t, 3> chords{};
    std::size_t count = 0;

    void add(std::size_t chord) {
        assert(count < chords.size() &&
               "[C91 §2.1 tex 72]: a vertex supplies two chords, or a null chord and two others");
        chords[count++] = chord;
    }
};

Subarc whole_arc(const Polygon& curve, const Submap& submap, std::size_t arc) {
    const Arc& a = submap.arc(arc);
    return {a.first_edge,
            a.first_side,
            a.last_edge,
            a.last_side,
            submap.arc_start_symbolic_y(arc, curve),
            submap.arc_end_symbolic_y(arc, curve)};
}

bool contains_endpoint(const Polygon& curve, const Submap& submap, std::size_t arc,
                       std::size_t edge, Side side, const SymbolicY& y) {
    return subarc_contains_point(whole_arc(curve, submap, arc), curve, edge, side, y, 0,
                                 curve.num_vertices() - 1);
}

void append_height_runs(const Polygon& curve, const std::vector<VertexOccurrence>& vertices,
                        std::size_t first, std::vector<HeightRun>& runs) {
    while (first < vertices.size()) {
        std::size_t last = first + 1;
        bool increasing = true;
        if (last < vertices.size())
            increasing = symbolic_y_leq(symbolic_y_of(curve.vertex(vertices[first].vertex)),
                                        symbolic_y_of(curve.vertex(vertices[last].vertex)));
        while (last < vertices.size()) {
            const SymbolicY previous = symbolic_y_of(curve.vertex(vertices[last - 1].vertex));
            const SymbolicY current = symbolic_y_of(curve.vertex(vertices[last].vertex));
            if (increasing ? symbolic_y_greater(previous, current)
                           : symbolic_y_less(previous, current))
                break;
            ++last;
        }
        runs.push_back({first, last - first, increasing});
        first = last;
    }
}

RayHit first_region_contact(const Polygon& curve, const Submap& submap, const RegionArcs& arcs,
                            std::size_t vertex, std::size_t source_edge, Side direction,
                            RegionCompletionWork& work) {
    const Point& origin = curve.vertex(vertex);
    const SymbolicY y = symbolic_y_of(origin);
    const SourceOffset offset = perturbed_x_offset(curve, y, source_edge);
    RayHit best;
    Exact best_distance = 0;
    auto consider = [&](std::size_t arc, std::size_t edge) {
        ++work.ray_edge_tests;
        Exact x;
        if (!edge_crossing_x(curve, edge, y, &x))
            return;
        const bool ascending = point_y_below(curve.vertex(edge), curve.vertex(edge + 1));
        const Side west = ascending ? LEFT : RIGHT;
        const Side struck = direction == RIGHT ? west : (west == LEFT ? RIGHT : LEFT);
        if (!contains_endpoint(curve, submap, arc, edge, struck, y))
            return;
        const Exact distance = direction == RIGHT ? x - origin.x : origin.x - x;
        const bool wrapped =
            distance < 0 ||
            (distance == 0 && !perturbed_hit_forward(curve, y, direction, offset, edge));
        if (!best.hit ||
            (wrapped != best.wrapped ? !wrapped
             : distance != best_distance
                 ? distance < best_distance
                 : ray_contact_precedes(curve, y, direction, edge, struck, best.edge, best.side))) {
            best = {true, x, y.y, edge, struck, wrapped, arc};
            best_distance = distance;
        }
    };
    for (const std::size_t arc : arcs) {
        ArcSideRange ranges[3];
        const std::size_t count = submap.arc(arc).side_ranges(0, curve.num_vertices() - 1, ranges);
        for (std::size_t r = 0; r < count; ++r) {
            for (std::size_t edge = curve.next_nonnull_edge(ranges[r].first_edge);
                 edge <= ranges[r].last_edge; edge = curve.next_nonnull_edge(edge + 1)) {
                consider(arc, edge);
                if (edge + 1 == curve.num_edges())
                    break;
            }
        }
        if (vertex > 0 && curve.edge_is_null(vertex - 1))
            consider(arc, vertex - 1);
        if (vertex < curve.num_edges() && curve.edge_is_null(vertex))
            consider(arc, vertex);
    }
    assert(best.hit &&
           "[C91 §4.2 tex 367, §2.1 tex 70]: a missing chord stays in its region and meets C");
    if (auto* trace = AnimationTrace::current())
        trace->ray(curve, origin, direction, best);
    return best;
}

std::size_t endpoint_arc(const Polygon& curve, const Submap& submap, const RegionArcs& arcs,
                         const ChordEndpoint& endpoint) {
    for (const std::size_t arc : arcs)
        if (contains_endpoint(curve, submap, arc, endpoint.edge_c, endpoint.side, endpoint.y))
            return arc;
    assert(false && "[C91 §4.2 tex 367]: both endpoints of a missing chord lie on its region arcs");
    return NONE;
}

bool source_is_present(const Polygon& curve, const std::vector<PendingChord>& chords,
                       const VertexChords& inventory, const VertexOccurrence& source) {
    const bool inside = is_inside_companion(curve, source.edge, source.side, source.vertex);
    const std::size_t edge =
        !curve.is_endpoint(source.vertex) && !curve.is_y_extremum(source.vertex) ? source.vertex - 1
                                                                                 : source.edge;
    for (std::size_t i = 0; i < inventory.count; ++i) {
        const PendingChord& chord = chords[inventory.chords[i]];
        if (inside ? chord.is_null_length
                   : ((chord.left_edge_c == edge && chord.left_side == source.side) ||
                      (chord.right_edge_c == edge && chord.right_side == source.side)))
            return true;
    }
    return false;
}

std::vector<ChordEndpoint> order_endpoints(const Polygon& curve, const Submap& submap,
                                           const std::vector<PendingChord>& chords,
                                           std::vector<std::vector<ChordEndpoint>>& new_endpoints,
                                           std::vector<std::vector<ChordEndpoint>>& old_endpoints,
                                           RegionCompletionWork& work) {
    const std::size_t n_edges = curve.num_edges();
    auto position = [&](std::size_t edge, Side side) {
        return side == LEFT ? edge : 2 * n_edges - 1 - edge;
    };
    auto increasing = [&](const ChordEndpoint& e) {
        return (e.side == LEFT) ==
               point_y_below(curve.vertex(e.edge_c), curve.vertex(e.edge_c + 1));
    };
    for (std::size_t arc = 0; arc < submap.num_arcs(); ++arc) {
        auto& old = old_endpoints[arc];
        assert(old.size() <= 4 &&
               "[C91 §2.4 tex 137]: an arc has two ends with bounded incidences");
        std::sort(old.begin(), old.end(),
                  [](const auto& a, const auto& b) { return symbolic_y_less(a.y, b.y); });
        auto& added = new_endpoints[arc];
        for (std::size_t i = 1; i < added.size(); ++i)
            assert(symbolic_y_leq(added[i - 1].y, added[i].y) &&
                   "[C91 §4.2 tex 367]: region vertices are processed by height");
        std::vector<ChordEndpoint> merged;
        merged.reserve(old.size() + added.size());
        std::merge(old.begin(), old.end(), added.begin(), added.end(), std::back_inserter(merged),
                   [](const auto& a, const auto& b) { return symbolic_y_less(a.y, b.y); });
        added = std::move(merged);
    }
    const Arc& start = submap.arc(submap.start_arc);
    const std::size_t start_position = position(start.first_edge, start.first_side);
    const SymbolicY start_y = submap.arc_start_symbolic_y(submap.start_arc, curve);
    auto in_head = [&](const ChordEndpoint& e) {
        const std::size_t p = position(e.edge_c, e.side);
        return p != start_position ? p < start_position
                                   : (increasing(e) ? symbolic_y_less(e.y, start_y)
                                                    : symbolic_y_greater(e.y, start_y));
    };
    std::vector<ChordEndpoint> sequence;
    sequence.reserve(2 * chords.size());
    auto append = [&](std::size_t arc, bool head) {
        const auto& list = new_endpoints[arc];
        auto belongs = [&](const auto& e) { return arc != submap.start_arc || in_head(e) == head; };
        for (const auto& e : list) {
            ++work.endpoint_visits;
            if (belongs(e) && increasing(e))
                sequence.push_back(e);
        }
        for (auto it = list.rbegin(); it != list.rend(); ++it) {
            ++work.endpoint_visits;
            if (belongs(*it) && !increasing(*it))
                sequence.push_back(*it);
        }
    };
    append(submap.start_arc, true);
    for (std::size_t arc = 0; arc < submap.num_arcs(); ++arc)
        append(arc, false);
    assert(sequence.size() == 2 * chords.size());
    std::vector<std::size_t> next(2 * n_edges + 1, 0);
    for (const auto& e : sequence)
        ++next[position(e.edge_c, e.side) + 1];
    for (std::size_t i = 1; i < next.size(); ++i)
        next[i] += next[i - 1];
    std::vector<ChordEndpoint> endpoints(sequence.size());
    for (const auto& e : sequence) {
        endpoints[next[position(e.edge_c, e.side)]++] = e;
        ++work.endpoint_visits;
    }
    for (std::size_t first = 0; first < endpoints.size();) {
        std::size_t last = first + 1;
        while (last < endpoints.size() && endpoints[first].edge_c == endpoints[last].edge_c &&
               endpoints[first].side == endpoints[last].side &&
               symbolic_y_equal(endpoints[first].y, endpoints[last].y))
            ++last;
        assert(last - first <= 4 &&
               "[C91 §2.1 tex 72, §2.4 tex 137]: a companion has bounded endpoint incidences");
        std::sort(endpoints.begin() + static_cast<std::ptrdiff_t>(first),
                  endpoints.begin() + static_cast<std::ptrdiff_t>(last),
                  [&](const auto& a, const auto& b) {
                      return chord_endpoint_precedes(curve, chords, a, b);
                  });
        work.endpoint_visits += last - first;
        first = last;
    }
    return endpoints;
}

}

RegionCompletion complete_bounded_regions(const Polygon& curve, const Submap& submap,
                                          [[maybe_unused]] std::size_t granularity) {
    assert(granularity >= 1 && granularity <= REGION_COMPLETION_GRANULARITY &&
           "[C91 §4.2 tex 367]: this is the constant-grade base case");
    assert(submap.is_conformal() && submap.is_semigranular(granularity) &&
           submap.is_granular(granularity, curve) &&
           "[C91 §4.2 tex 367]: bounded-region completion requires bounded conformal regions");
    assert(!submap.tree_decomposition().empty() &&
           "[C91 Lemma 4.2 tex 364]: the input submap is in normal form");
    RegionCompletion result;
    std::vector<PendingChord> chords;
    chords.reserve(3 * curve.num_vertices());
    std::vector<VertexChords> inventory(curve.num_vertices());
    std::vector<std::vector<ChordEndpoint>> new_endpoints(submap.num_arcs());
    std::vector<std::vector<ChordEndpoint>> old_endpoints(submap.num_arcs());
    const std::size_t first_tag = curve.vertex(0).index;
    auto endpoint = [&](std::size_t chord, bool left) {
        const auto& c = chords[chord];
        return ChordEndpoint{left ? c.left_edge_c : c.right_edge_c,
                             left ? c.left_side : c.right_side, c.y, chord, left};
    };
    for (std::size_t i = 0; i < submap.num_chords(); ++i) {
        const Chord& c = submap.chord(i);
        assert(!c.dead && "[C91 Lemma 4.2 tex 364]: normal form has no deleted records");
        const std::size_t index = chords.size();
        chords.push_back({c.symbolic_y(), c.left_edge, c.left_side, c.right_edge, c.right_side,
                          c.is_null_length});
        canonicalize_chord(chords.back(), curve);
        if (c.y_tag >= first_tag && c.y_tag - first_tag < inventory.size()) {
            const std::size_t v = c.y_tag - first_tag;
            auto at_vertex = [&](std::size_t edge) {
                return symbolic_y_equal(c.symbolic_y(), symbolic_y_of(curve.vertex(v))) &&
                       (edge == v || edge + 1 == v);
            };
            if (at_vertex(c.left_edge) || at_vertex(c.right_edge))
                inventory[v].add(index);
        }
        old_endpoints[c.left_adj.arcs[0]].push_back(endpoint(index, true));
        old_endpoints[c.right_adj.arcs[0]].push_back(endpoint(index, false));
    }
    for (std::size_t region = 0; region < submap.num_nodes(); ++region) {
        const RegionArcs arcs = collect_region_arcs(submap, region);
        std::vector<VertexOccurrence> vertices;
        std::vector<HeightRun> runs;
        for (const std::size_t arc : arcs) {
            const std::size_t first = vertices.size();
            ArcSideRange ranges[3];
            const std::size_t count =
                submap.arc(arc).side_ranges(0, curve.num_vertices() - 1, ranges);
            for (std::size_t r = 0; r < count; ++r) {
                auto visit = [&](std::size_t edge) {
                    for (std::size_t j = 0; j < 2; ++j) {
                        const std::size_t v = ranges[r].side == LEFT ? edge + j : edge + 1 - j;
                        if (contains_endpoint(curve, submap, arc, edge, ranges[r].side,
                                              symbolic_y_of(curve.vertex(v))))
                            vertices.push_back({v, edge, ranges[r].side});
                    }
                };
                if (ranges[r].side == LEFT) {
                    for (std::size_t edge = ranges[r].first_edge; edge <= ranges[r].last_edge;
                         ++edge)
                        visit(edge);
                } else {
                    for (std::size_t edge = ranges[r].last_edge + 1; edge-- > ranges[r].first_edge;)
                        visit(edge);
                }
            }
            append_height_runs(curve, vertices, first, runs);
        }
        assert(runs.size() <= 4 * (2 * REGION_COMPLETION_GRANULARITY + 9) &&
               "[C91 Lemma 2.3 tex 124]: bounded regions have a constant number of nonnull edges");
        result.work.vertex_occurrences += vertices.size();
        while (true) {
            std::size_t chosen = NONE;
            for (std::size_t i = 0; i < runs.size(); ++i) {
                if (runs[i].count == 0)
                    continue;
                ++result.work.height_comparisons;
                if (chosen == NONE ||
                    point_y_below(curve.vertex(vertices[runs[i].next()].vertex),
                                  curve.vertex(vertices[runs[chosen].next()].vertex)))
                    chosen = i;
            }
            if (chosen == NONE)
                break;
            const VertexOccurrence source = vertices[runs[chosen].next()];
            runs[chosen].advance();
            VertexChords& existing = inventory[source.vertex];
            if (source_is_present(curve, chords, existing, source))
                continue;
            PendingChord chord;
            chord.y = symbolic_y_of(curve.vertex(source.vertex));
            if (is_inside_companion(curve, source.edge, source.side, source.vertex)) {
                chord.left_edge_c = chord.right_edge_c = source.vertex;
                chord.left_side = chord.right_side =
                    is_inside_companion(curve, source.vertex, LEFT, source.vertex) ? LEFT : RIGHT;
                chord.is_null_length = true;
            } else {
                const RayHit hit = first_region_contact(
                    curve, submap, arcs, source.vertex, source.edge,
                    shooting_direction(source.edge, source.side, curve), result.work);
                chord.left_edge_c = source.edge;
                chord.left_side = source.side;
                chord.right_edge_c = hit.edge;
                chord.right_side = hit.side;
            }
            canonicalize_chord(chord, curve);
            const std::size_t index = chords.size();
            chords.push_back(std::move(chord));
            existing.add(index);
            for (const bool left : {true, false}) {
                const ChordEndpoint e = endpoint(index, left);
                new_endpoints[endpoint_arc(curve, submap, arcs, e)].push_back(e);
            }
        }
    }
    const auto endpoints =
        order_endpoints(curve, submap, chords, new_endpoints, old_endpoints, result.work);
#ifndef NDEBUG
    for (std::size_t v = 0; v < curve.num_vertices(); ++v) {
        for (const Side side : {LEFT, RIGHT}) {
            if (v > 0)
                assert(source_is_present(curve, chords, inventory[v], {v, v - 1, side}) &&
                       "[C91 §2.1 tex 72]: every companion is incident upon its visibility chord");
            if (v < curve.num_edges())
                assert(source_is_present(curve, chords, inventory[v], {v, v, side}) &&
                       "[C91 §2.1 tex 72]: every companion is incident upon its visibility chord");
        }
    }
#endif
    build_submap_from_ordered_chords(result.visibility_map, curve, chords, endpoints);
    assert(result.visibility_map.is_conformal() && result.visibility_map.is_semigranular(1) &&
           "[C91 §2.1 tex 70–72, §4.2 tex 367]: completing every companion produces V(C)");
    return result;
}

}
