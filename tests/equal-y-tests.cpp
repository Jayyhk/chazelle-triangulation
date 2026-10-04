#include "merge/ray_shooting.h"
#include "polygon/polygon.h"
#include "visibility/naive_visibility.h"
#include "visibility/up_phase.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <tuple>
#include <utility>
#include <vector>

using namespace chazelle;

namespace {

struct ExactPerturbedHit {
    bool hit = false;
    bool wrapped = false;
    Exact distance = 0;
    std::vector<std::size_t> edges;
};

Exact perturbed_height(const Point& point, const Exact& epsilon) {
    return point.y - Exact(point.index) * epsilon;
}

ExactPerturbedHit exact_perturbed_first_contact(const Polygon& curve, const Point& origin,
                                                Side direction, const Exact& epsilon) {
    const Exact query_y = perturbed_height(origin, epsilon);
    ExactPerturbedHit best;
    for (std::size_t edge = 0; edge < curve.num_edges(); ++edge) {
        const Point& start = curve.vertex(edge);
        const Point& end = curve.vertex(edge + 1);
        const Exact start_y = perturbed_height(start, epsilon);
        const Exact end_y = perturbed_height(end, epsilon);
        const bool crosses =
            (start_y <= query_y && query_y <= end_y) || (end_y <= query_y && query_y <= start_y);
        if (!crosses)
            continue;
        assert(start_y != end_y);
        const Exact x = start.x + (query_y - start_y) * (end.x - start.x) / (end_y - start_y);
        const Exact distance = direction == RIGHT ? x - origin.x : origin.x - x;
        const bool wrapped = distance <= 0;
        if (!best.hit || (wrapped != best.wrapped ? !wrapped : distance < best.distance)) {
            best.hit = true;
            best.wrapped = wrapped;
            best.distance = distance;
            best.edges.clear();
            best.edges.push_back(edge);
        } else if (wrapped == best.wrapped && distance == best.distance) {
            best.edges.push_back(edge);
        }
    }
    return best;
}

void check_exact_perturbed_hit(const Polygon& curve, const Point& origin, Side direction,
                               const RayHit& actual, const Exact& epsilon) {
    const ExactPerturbedHit expected =
        exact_perturbed_first_contact(curve, origin, direction, epsilon);
    assert(actual.hit == expected.hit && "[C91 §2 tex 47]: symbolic perturbation preserves hits");
    if (!actual.hit)
        return;
    assert(actual.wrapped == expected.wrapped &&
           "[C91 §2.1 tex 70]: symbolic perturbation preserves wrapping through infinity");
    assert(std::find(expected.edges.begin(), expected.edges.end(), actual.edge) !=
               expected.edges.end() &&
           "[C91 §2 tex 47, §3.4 Lemma 3.6]: the symbolic ray must find the first perturbed edge");
    const Point& start = curve.vertex(actual.edge);
    const Point& end = curve.vertex(actual.edge + 1);
    const Exact expected_x =
        origin.index == start.index
            ? start.x
            : (origin.index == end.index
                   ? end.x
                   : start.x + (origin.y - start.y) * (end.x - start.x) / (end.y - start.y));
    assert(actual.x == expected_x &&
           "[C91 §2 tex 47]: reported coordinates retain their exact unperturbed values");
    const bool ascending = perturbed_height(start, epsilon) < perturbed_height(end, epsilon);
    const Side face_to_left = ascending ? LEFT : RIGHT;
    const Side expected_side =
        direction == RIGHT ? face_to_left : (face_to_left == LEFT ? RIGHT : LEFT);
    assert(actual.side == expected_side &&
           "[C91 §2.1 tex 72]: the ray strikes the appropriate double-boundary side");
}

using ChordKey = std::tuple<std::size_t, std::size_t, Side, std::size_t, Side, bool>;

std::vector<ChordKey> chord_keys(const Submap& submap) {
    std::vector<ChordKey> keys;
    keys.reserve(submap.num_live_chords());
    for (std::size_t index = 0; index < submap.num_chords(); ++index) {
        const Chord& chord = submap.chord(index);
        if (chord.dead)
            continue;
        auto first = std::pair{chord.left_edge, chord.left_side};
        auto second = std::pair{chord.right_edge, chord.right_side};
        if (second < first)
            std::swap(first, second);
        keys.emplace_back(chord.y_tag, first.first, first.second, second.first, second.second,
                          chord.is_null_length);
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

void check_full_map_perturbation(const Polygon& curve, const Exact& epsilon) {
    std::vector<Point> vertices;
    vertices.reserve(curve.num_vertices());
    for (std::size_t index = 0; index < curve.num_vertices(); ++index) {
        Point vertex = curve.vertex(index);
        vertex.y = perturbed_height(vertex, epsilon);
        vertices.push_back(std::move(vertex));
    }
    const Polygon perturbed_curve(std::move(vertices));
    Submap original = build_full_visibility_map(curve);
    Submap perturbed = build_full_visibility_map(perturbed_curve);
    original.check_invariants(curve);
    perturbed.check_invariants(perturbed_curve);
    assert(
        chord_keys(original) == chord_keys(perturbed) &&
        "[C91 §2 tex 47, §2.1]: symbolic and exact perturbations produce the same visibility chords");
}

void check_double_identification(const Submap& submap, const Polygon& curve) {
    for (std::size_t edge = 0; edge < curve.num_edges(); ++edge) {
        for (std::size_t vertex = 0; vertex < curve.num_vertices(); ++vertex) {
            const SymbolicY query_y = symbolic_y_of(curve.vertex(vertex));
            Exact query_x;
            if (!edge_crossing_x(curve, edge, query_y, &query_x))
                continue;
            std::vector<std::size_t> expected;
            for (std::size_t arc = 0; arc < submap.num_arcs(); ++arc) {
                const Arc& entry = submap.arc(arc);
                const Subarc boundary{entry.first_edge,
                                      entry.first_side,
                                      entry.last_edge,
                                      entry.last_side,
                                      submap.arc_start_symbolic_y(arc, curve),
                                      submap.arc_end_symbolic_y(arc, curve)};
                if (subarc_contains_point(boundary, curve, edge, LEFT, query_y, 0,
                                          curve.num_vertices() - 1) ||
                    subarc_contains_point(boundary, curve, edge, RIGHT, query_y, 0,
                                          curve.num_vertices() - 1))
                    expected.push_back(arc);
            }
            const auto actual = submap.double_identify(edge, query_y, curve);
            assert(
                actual.count == expected.size() &&
                "[C91 §2.4 tex 144]: equal raw heights preserve exhaustive double identification");
            for (std::size_t arc : expected)
                assert(std::find(actual.arcs.data(), actual.arcs.data() + actual.count, arc) !=
                       actual.arcs.data() + actual.count);
        }
    }
}

void check_curve(std::vector<Point> vertices) {
    const Polygon original(vertices);
    for (const Exact& epsilon : {Exact(0x1p-160), Exact(0x1p-240)})
        check_full_map_perturbation(original, epsilon);

    const UpPhase up(std::move(vertices));
    for (std::size_t grade = 0; grade <= up.graded().maximum_grade(); ++grade) {
        for (std::size_t chain = 0; chain < up.graded().num_chains(grade); ++chain) {
            const Polygon& curve = up.graded().chain(grade, chain);
            const Submap& submap = up.chain_submap(grade, chain);
            submap.check_invariants(curve);
            assert(submap.is_conformal() &&
                   submap.is_granular(UpPhase::grade_granularity(grade), curve) &&
                   "[C91 §4.1 tex 327]: equal raw heights preserve canonical submaps");
            const RayShootingStructure& shooter = up.chain_structure(grade, chain);
            for (std::size_t vertex = 0; vertex < curve.num_vertices(); ++vertex) {
                const Point& origin = curve.vertex(vertex);
                for (Side direction : {LEFT, RIGHT}) {
                    const RayHit actual = shooter.shoot_toward_boundary(origin, direction);
                    for (const Exact& epsilon : {Exact(0x1p-160), Exact(0x1p-240)})
                        check_exact_perturbed_hit(curve, origin, direction, actual, epsilon);
                }
            }
        }
    }
    const std::size_t last = up.graded().curve().num_vertices() - 1;
    const std::size_t first = last >= 8 ? 1 : 0;
    const std::size_t end = last >= 8 ? last - 1 : last;
    auto portion = up.compute_canonical_portion(first, end);
    portion.submap.check_invariants(portion.curve);
    assert(portion.submap.is_conformal() &&
           portion.submap.is_granular(canonical_granularity(end - first), portion.curve));
    check_double_identification(portion.submap, portion.curve);
}

std::vector<std::vector<Point>> curves() {
    std::vector<std::vector<Point>> fixtures{
        {{0, 0, 0}, {12, 0, 1}, {12, 12, 2}, {0, 12, 3}},
        {{0, 0, 0},
         {12, 0, 1},
         {12, 12, 2},
         {0, 12, 3},
         {0, 2, 4},
         {10, 2, 5},
         {10, 10, 6},
         {2, 10, 7},
         {2, 4, 8},
         {8, 4, 9},
         {8, 8, 10},
         {4, 8, 11},
         {4, 6, 12},
         {6, 6, 13}},
    };
    std::vector<Point> horizontal;
    horizontal.reserve(66);
    for (std::size_t index = 0; index < 66; ++index)
        horizontal.push_back({Exact(index), 0, index});
    fixtures.push_back(std::move(horizontal));
    std::vector<Point> comb{{0, 0, 0}};
    comb.reserve(69);
    for (std::size_t tooth = 1; tooth < 18; ++tooth) {
        const Exact x = Exact(2 * tooth);
        for (const auto& coordinate : {std::pair{x - 1, Exact(0)}, std::pair{x - 1, Exact(8)},
                                       std::pair{x, Exact(8)}, std::pair{x, Exact(0)}})
            comb.push_back({coordinate.first, coordinate.second, comb.size()});
    }
    fixtures.push_back(std::move(comb));
    std::vector<Point> snake;
    snake.reserve(72);
    for (std::size_t row = 0; row < 9; ++row) {
        for (std::size_t column = 0; column < 8; ++column) {
            const std::size_t x = row % 2 == 0 ? column : 7 - column;
            snake.push_back({Exact(x), Exact(2 * row), snake.size()});
        }
    }
    fixtures.push_back(std::move(snake));
    return fixtures;
}

void check_five_vertex_height_sequences() {
    for (std::size_t sequence = 0; sequence < 243; ++sequence) {
        std::size_t remaining_heights = sequence;
        std::vector<Point> vertices;
        vertices.reserve(5);
        for (std::size_t index = 0; index < 5; ++index) {
            vertices.push_back({Exact(2 * index), Exact(remaining_heights % 3), index});
            remaining_heights /= 3;
        }
        check_curve(std::move(vertices));
    }
}

void check_large_vertex_tags() {
    std::vector<Point> vertices;
    vertices.reserve(33);
    for (std::size_t index = 0; index < 33; ++index)
        vertices.push_back({Exact(index), 0, SOS_NONE - 257 + index});
    check_curve(std::move(vertices));
}

}

int main() {
    for (const auto& fixture : curves()) {
        for (std::size_t rotation = 0; rotation < 4; ++rotation) {
            for (bool reverse : {false, true}) {
                std::vector<Point> vertices = fixture;
                if (reverse)
                    std::reverse(vertices.begin(), vertices.end());
                for (std::size_t index = 0; index < vertices.size(); ++index) {
                    Point& vertex = vertices[index];
                    for (std::size_t turn = 0; turn < rotation; ++turn) {
                        const Exact x = vertex.x;
                        vertex.x = -vertex.y;
                        vertex.y = x;
                    }
                    vertex.index = index + 17;
                }
                check_curve(std::move(vertices));
            }
        }
    }
    check_five_vertex_height_sequences();
    check_large_vertex_tags();
    std::printf("[C91 §2/§3/§4.1 equal-y tests]: all passed\n");
}
