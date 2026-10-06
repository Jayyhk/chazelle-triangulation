#include "algorithm/submap/chord_inventory.h"
#include "algorithm/visibility/down_phase.h"
#include "support/pipeline_fixtures.h"
#include "support/random.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <tuple>

using namespace chazelle;

namespace {

using ChordKey = std::tuple<Exact, std::size_t, std::size_t, Side, std::size_t, Side, bool>;

std::vector<ChordKey> chord_keys(const Submap& map, const Polygon& curve,
                                 bool vertex_chords_only = false) {
    std::vector<ChordKey> result;
    result.reserve(map.num_chords());
    for (std::size_t index = 0; index < map.num_chords(); ++index) {
        const Chord& chord = map.chord(index);
        assert(!chord.dead);
        PendingChord value{chord.symbolic_y(), chord.left_edge,  chord.left_side,
                           chord.right_edge,   chord.right_side, chord.is_null_length};
        canonicalize_chord(value, curve);
        const std::size_t source = curve.local_index_of_tag(value.y.tag);
        assert(source != NONE && symbolic_y_equal(value.y, symbolic_y_of(curve.vertex(source))));
        if (vertex_chords_only && value.left_edge_c != source && value.left_edge_c + 1 != source &&
            value.right_edge_c != source && value.right_edge_c + 1 != source)
            continue;
        if (value.is_null_length)
            value.left_edge_c = value.right_edge_c = source;
        result.emplace_back(value.y.y, value.y.tag, value.left_edge_c, value.left_side,
                            value.right_edge_c, value.right_side, value.is_null_length);
    }
    std::sort(result.begin(), result.end());
    return result;
}

struct SubmapState {
    std::vector<std::size_t> fields;
    std::vector<Exact> heights;
    bool operator==(const SubmapState&) const = default;
};

SubmapState submap_state(const Submap& map) {
    SubmapState state;
    const auto append = [&](std::initializer_list<std::size_t> fields) {
        state.fields.insert(state.fields.end(), fields);
    };
    append({map.start_vertex, map.end_vertex, map.start_arc, map.end_arc, map.left_right_boundary(),
            map.num_nodes(), map.num_arcs(), map.num_chords(), map.tree_decomposition().root()});
    for (std::size_t index = 0; index < map.num_nodes(); ++index) {
        const auto& node = map.node(index);
        append({node.dead, node.incident_chords.size()});
        state.fields.insert(state.fields.end(), node.incident_chords.begin(),
                            node.incident_chords.end());
    }
    for (std::size_t index = 0; index < map.num_arcs(); ++index) {
        const auto& arc = map.arc(index);
        append({arc.first_edge, arc.first_side, arc.last_edge, arc.last_side, arc.region_node,
                arc.edge_count, arc.dead});
    }
    for (std::size_t index = 0; index < map.num_chords(); ++index) {
        const auto& chord = map.chord(index);
        append({chord.region[0], chord.region[1], chord.left_adj.count, chord.left_adj.arcs[0],
                chord.left_adj.arcs[1], chord.right_adj.count, chord.right_adj.arcs[0],
                chord.right_adj.arcs[1], chord.left_edge, chord.left_side, chord.right_edge,
                chord.right_side, chord.y_tag, chord.is_null_length, chord.dead});
        state.heights.push_back(chord.y);
    }
    for (std::size_t index = 0; index < map.tree_decomposition().size(); ++index) {
        const auto& node = map.tree_decomposition().node(index);
        append({node.chord_idx, node.region_idx, node.parent, node.left_child, node.right_child});
    }
    return state;
}

std::vector<SubmapState> stored_submaps(const UpPhase& up) {
    std::vector<SubmapState> result;
    for (std::size_t grade = 0; grade < up.graded().num_grades(); ++grade)
        for (std::size_t chain = 0; chain < up.graded().num_chains(grade); ++chain)
            result.push_back(submap_state(up.chain_submap(grade, chain)));
    return result;
}

void check_result(const UpPhase& up, std::size_t grade, std::size_t chain) {
    const Polygon& curve = up.graded().chain(grade, chain);
    const auto& canonical = up.chain_submap(grade, chain);
    canonical.check_invariants(curve);
    assert(canonical.is_conformal() &&
           canonical.is_granular(UpPhase::grade_granularity(grade), curve));
    const auto result = compute_visibility_map(up, grade, chain);
    const Submap& map = result.visibility_map;
    map.check_invariants(curve);
    assert(map.is_conformal() && map.is_semigranular(1));
    assert(map.num_nodes() == map.num_live_nodes() && map.num_chords() == map.num_live_chords() &&
           map.num_arcs() == map.num_live_arcs() && !map.tree_decomposition().empty());
    assert(map.num_nodes() == map.num_chords() + 1 &&
           map.tree_decomposition().size() == 2 * map.num_nodes() - 1);
    const auto reference = build_full_visibility_map(curve);
    const auto expected = chord_keys(reference, curve);
    const auto actual = chord_keys(map, curve);
    assert(actual == expected &&
           "[C91 Theorem 4.3 tex 390]: up-phase followed by down-phase recovers V(C)");
    const auto inherited = chord_keys(canonical, curve, true);
    assert(std::includes(actual.begin(), actual.end(), inherited.begin(), inherited.end()) &&
           "[C91 §2.2 tex 90, §4.2]: retain vertex chords when discarding augmentation");
    std::size_t parameter = grade;
    std::size_t rounds = 0;
    while (UpPhase::grade_granularity(parameter) > REGION_COMPLETION_GRANULARITY) {
        parameter = ceil_beta(parameter);
        ++rounds;
    }
    assert(result.refinement_rounds == rounds && (rounds == 0 || result.region_boundaries > 0));
    const auto& work = result.completion_work;
    assert(work.vertex_occurrences <= 64 * curve.num_vertices() &&
           work.height_comparisons <= 68 * work.vertex_occurrences &&
           work.ray_edge_tests <= 2048 * curve.num_vertices() &&
           work.endpoint_visits <= 64 * curve.num_vertices());
}

void check_curve(const std::vector<Point>& vertices, bool all_chains) {
    const UpPhase up(vertices);
    const auto& graded = up.graded();
    const auto positions = graded.original_vertex_positions();
    assert(positions.size() == vertices.size() && positions.front() == 0 &&
           positions.back() + 1 == graded.curve().num_vertices());
    assert(graded.curve().num_edges() == (std::size_t{1} << graded.maximum_grade()) &&
           graded.curve().num_vertices() < 2 * vertices.size());
    for (std::size_t vertex = 0; vertex < vertices.size(); ++vertex) {
        const Point& retained = graded.curve().vertex(positions[vertex]);
        assert(retained.x == vertices[vertex].x && retained.y == vertices[vertex].y);
    }
    const auto before = stored_submaps(up);
    for (std::size_t grade = 0; grade < graded.num_grades(); ++grade) {
        const std::size_t count = graded.num_chains(grade);
        if (all_chains) {
            for (std::size_t chain = 0; chain < count; ++chain)
                check_result(up, grade, chain);
        } else {
            check_result(up, grade, count - 1);
        }
    }
    assert(stored_submaps(up) == before &&
           "[C91 §4.2 tex 362-367]: refinement preserves the stored up-phase data");
    if (all_chains) {
        const auto first = compute_visibility_map(up);
        const auto second = compute_visibility_map(up);
        assert(chord_keys(first.visibility_map, graded.curve()) ==
               chord_keys(second.visibility_map, graded.curve()));
    }
    assert(stored_submaps(up) == before);
}

}

int main() {
    for (const auto& fixture : test::polygon_fixtures())
        for (bool reverse : {false, true})
            check_curve(test::boundary_order(fixture, reverse, fixture.size() / 2), true);
    test::DeterministicRandomGenerator random(41991);
    for (std::size_t sample = 0; sample < 12; ++sample) {
        std::vector<Point> vertices;
        const std::size_t count = 18 + random.next() % 24;
        vertices.reserve(count);
        for (std::size_t vertex = 0; vertex < count; ++vertex)
            vertices.push_back({Exact(vertex) / 3, Exact(random.next() % 7) / 5, vertex});
        for (Point& point : vertices) {
            if (sample % 3 == 0)
                std::swap(point.x, point.y);
            if (sample % 3 == 1)
                point.x = -point.x;
        }
        check_curve(test::boundary_order(vertices, sample % 2 == 0, 0, 101), true);
    }
    for (std::size_t count : {2U, 1025U, 1026U, 2050U}) {
        std::vector<Point> vertices;
        vertices.reserve(count);
        for (std::size_t vertex = 0; vertex < count; ++vertex)
            vertices.push_back({0, Exact(vertex), vertex + 23});
        check_curve(vertices, false);
    }
    std::puts("Section 4 end-to-end tests passed");
}
