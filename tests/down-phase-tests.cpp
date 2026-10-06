#include "algorithm/merge/granularity.h"
#include "algorithm/submap/chord_inventory.h"
#include "algorithm/visibility/down_phase.h"
#include "support/assertions.h"
#include "support/random.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <tuple>

using namespace chazelle;

namespace {

using ChordKey = std::tuple<std::size_t, std::size_t, Side, std::size_t, Side, bool>;

std::vector<ChordKey> keys(const Submap& map, const Polygon& curve) {
    std::vector<ChordKey> result;
    for (std::size_t i = 0; i < map.num_chords(); ++i) {
        const Chord& chord = map.chord(i);
        PendingChord value{chord.symbolic_y(), chord.left_edge,  chord.left_side,
                           chord.right_edge,   chord.right_side, chord.is_null_length};
        canonicalize_chord(value, curve);
        if (value.is_null_length) {
            const std::size_t vertex = curve.local_index_of_tag(value.y.tag);
            value.left_edge_c = value.right_edge_c = vertex;
        }
        result.emplace_back(value.y.tag, value.left_edge_c, value.left_side, value.right_edge_c,
                            value.right_side, value.is_null_length);
    }
    std::sort(result.begin(), result.end());
    return result;
}

void check_map(const Polygon& curve, const Submap& result, bool augmented = false) {
    result.check_invariants(curve);
    assert(result.is_conformal() && result.is_semigranular(1));
    assert(!result.tree_decomposition().empty());
    const auto reference = build_full_visibility_map(curve);
    const auto actual = keys(result, curve);
    const auto expected = keys(reference, curve);
    assert((augmented
                ? std::includes(actual.begin(), actual.end(), expected.begin(), expected.end())
                : actual == expected) &&
           "[C91 Lemma 4.2 tex 364–381]: the down-phase recovers exactly V(C)");
}

void check_visibility(std::vector<Point> vertices) {
    const UpPhase up(std::move(vertices));
    const Polygon& curve = up.graded().curve();
    const auto result = compute_visibility_map(up);
    check_map(curve, result.visibility_map);
}

void check_refinement(std::vector<Point> vertices) {
    const UpPhase up(std::move(vertices));
    const Polygon& curve = up.graded().curve();
    assert(up.graded().maximum_grade() >= 6);
    Submap coarse = up.chain_submap(up.graded().maximum_grade(), 0);
    enforce_granularity(coarse, curve, 4);
    coarse.normalize(curve);
    Submap refined = refine_visibility_submap(up, curve, coarse, 6);
    assert(refined.is_conformal() && refined.is_granular(2, curve));
    auto full = complete_bounded_regions(curve, refined, 2).visibility_map;
    full.build_tree_decomposition();
    check_map(curve, full, true);
}

std::vector<Point> subdivide_curve(const std::vector<Point>& vertices, std::size_t edge_count) {
    const std::size_t original_edges = vertices.size() - 1;
    assert(edge_count >= original_edges);
    std::vector<Point> subdivided;
    subdivided.reserve(edge_count + 1);
    for (std::size_t edge = 0; edge < original_edges; ++edge) {
        const Point& first = vertices[edge];
        const Point& last = vertices[edge + 1];
        const std::size_t count =
            edge_count / original_edges + (edge < edge_count % original_edges);
        for (std::size_t part = 0; part < count; ++part) {
            const Exact position = Exact{part} / count;
            subdivided.push_back({first.x + position * (last.x - first.x),
                                  first.y + position * (last.y - first.y),
                                  vertices.front().index + subdivided.size()});
        }
    }
    Point last = vertices.back();
    last.index = vertices.front().index + subdivided.size();
    subdivided.push_back(std::move(last));
    return subdivided;
}

void check_infinitesimals() {
    const Exact first = Exact::infinitesimal(3);
    const Exact second = Exact::infinitesimal(2);
    const Exact third = Exact::infinitesimal(1);
    const Exact fourth = Exact::infinitesimal(0);
    assert(first > 0 && first < Exact{1} / Exact{2});
    assert(second < first * first && third < second * second && fourth < third * third);
    assert(1 / first > Exact{std::numeric_limits<std::uint64_t>::max()});
    assert((1 + first) / (1 - first) > 1 && (1 + first) * (1 - first) == 1 - first * first);
    Exact value = (first + second) / (first - second);
    value += value;
    assert(value > 2 && value < 3);
    value /= value;
    assert(value == 1);
    const Exact quotient = (first * second + third) / (first * second);
    assert(quotient > 1 && quotient < 1 + second);
    const Exact finite = (2 + first) / (1 - second);
    assert(finite > 1 && finite < 3 && -finite < -1 && -finite > -3);
    assert(finite * first > 0 && finite * first < 1 && finite / first > 100);
    assert((1 + second) / (1 - first) > 1 + first);
}

void check_boundary_views() {
    const Exact period = 4 / Exact::infinitesimal(3);
    const Polygon west = Polygon::symbolic_curve({{2, 0, 71}, {1, 1, 73}}, -period);
    const Polygon east = Polygon::symbolic_curve({{1, 1, 73}, {3, 0, 75}}, period);
    const Polygon joined(west, east);
    assert(joined.num_edges() == 2 && joined.vertex(1).index == 73);
    assert(&joined.vertex(0) == &west.vertex(0) && &joined.vertex(2) == &east.vertex(1));
    assert(joined.edge_horizontal_delta(0) == -1 - period &&
           joined.edge_horizontal_delta(1) == 2 + period);
    const SymbolicY low{Exact{1} / 4, 79};
    const SymbolicY high{Exact{3} / 4, 81};
    const Exact west_x = edge_x_at_y(joined, 0, low);
    const Exact east_x = edge_x_at_y(joined, 1, high);
    assert(west_x == Exact{7} / 4 - period / 4);
    assert(east_x == Exact{3} / 2 + period / 4);
    const Polygon reversed = joined.reversed();
    assert(edge_x_at_y(reversed, 1, low) == west_x && edge_x_at_y(reversed, 0, high) == east_x);
    assert(reversed.local_index_of_tag(71) == 2 && reversed.local_index_of_tag(75) == 0);
    assert(reversed.subchain(0, 2).edge_horizontal_shift(0) == -period);
    assert(joined.previous_branch_left(1) && !reversed.previous_branch_left(1));
    const Polygon equal_x = Polygon::symbolic_curve({{1, 0, 71}, {1, 1, 73}}, period);
    const SymbolicY below_end{1, 75};
    const SymbolicY above_start{0, 69};
    assert(perturbed_x_offset(equal_x, below_end, 0) == -period &&
           perturbed_x_offset(equal_x, above_start, 0) == period);
    assert(perturbed_x_offset(equal_x.reversed(), below_end, 0) == -period);
    assert(perturbed_hit_forward(equal_x, below_end, LEFT, Exact{0}, 0) &&
           !perturbed_hit_forward(equal_x, below_end, RIGHT, Exact{0}, 0));
}

void check_random_refinements() {
    chazelle::test::DeterministicRandomGenerator random(83);
    for (std::size_t sample = 0; sample < 24; ++sample) {
        std::vector<Point> vertices;
        const std::size_t count = 9 + random.next() % 20;
        vertices.reserve(count);
        for (std::size_t i = 0; i < count; ++i)
            vertices.push_back({Exact{i}, Exact{random.next() % 7}, i});
        if (sample % 2 == 0)
            std::reverse(vertices.begin(), vertices.end());
        for (std::size_t i = 0; i < count; ++i) {
            if (sample % 3 == 0) {
                const Exact x = vertices[i].x;
                vertices[i].x = -vertices[i].y;
                vertices[i].y = x;
            }
            vertices[i].index = i + 17;
        }
        check_visibility(std::move(vertices));
    }
}

void check_valid_refinements() {
    chazelle::test::DeterministicRandomGenerator random(337);
    for (std::size_t sample = 0; sample < 8; ++sample) {
        const std::size_t count = 65 + random.next() % 16;
        std::vector<Point> vertices;
        vertices.reserve(count);
        for (std::size_t i = 0; i < count; ++i)
            vertices.push_back({Exact{i} / 3, Exact{random.next() % 17} / 5, i});
        if (sample % 2 == 0)
            std::reverse(vertices.begin(), vertices.end());
        for (std::size_t i = 0; i < count; ++i) {
            if (sample % 4 == 0) {
                const Exact x = vertices[i].x;
                vertices[i].x = -vertices[i].y;
                vertices[i].y = x;
            } else if (sample % 4 == 1) {
                vertices[i].x = -vertices[i].x;
            } else if (sample % 4 == 2) {
                vertices[i].y = -vertices[i].y;
            }
            vertices[i].index = i + 101;
        }
        check_refinement(std::move(vertices));
    }
}

void check_exhaustive_small_curves() {
    std::size_t cases = 9;
    for (std::size_t count = 2; count <= 6; ++count) {
        for (std::size_t sample = 0; sample < cases; ++sample) {
            std::vector<Point> vertices;
            std::size_t heights = sample;
            for (std::size_t i = 0; i < count; ++i) {
                vertices.push_back({Exact{i}, Exact{heights % 3}, i});
                heights /= 3;
            }
            check_visibility(std::move(vertices));
        }
        cases *= 3;
    }
}

void check_refinement_grade_contract() {
    const UpPhase up({{0, 0, 0}, {1, 1, 1}, {2, 0, 2}, {3, 1, 3}, {4, 0, 4}});
    const Polygon& curve = up.graded().curve();
    Submap coarse = up.chain_submap(2, 0);
    enforce_granularity(coarse, curve, 8);
    coarse.normalize(curve);
    chazelle::test::require_assertion_abort(
        [&] { (void)refine_visibility_submap(up, curve, coarse, 11); });
}

void check_exact_height_perturbation() {
    std::vector<Point> vertices;
    vertices.reserve(65);
    for (std::size_t i = 0; i < 65; ++i)
        vertices.push_back({Exact{i} / 3, Exact{(13 * i) % 17} / 5, i + 101});
    const UpPhase up(std::move(vertices));
    const Polygon& curve = up.graded().curve();
    vertices.clear();
    const Exact epsilon = Exact{1} / (Exact{1000000} * curve.num_vertices() * curve.num_vertices());
    for (std::size_t i = 0; i < curve.num_vertices(); ++i) {
        Point point = curve.vertex(i);
        point.y -= Exact{point.index} * epsilon;
        vertices.push_back(std::move(point));
    }
    const UpPhase perturbed(std::move(vertices));
    const Polygon& perturbed_curve = perturbed.graded().curve();
    for (std::size_t first = 0; first < curve.num_vertices(); ++first)
        for (std::size_t second = 0; second < curve.num_vertices(); ++second)
            assert(point_y_below(curve.vertex(first), curve.vertex(second)) ==
                   (perturbed_curve.vertex(first).y < perturbed_curve.vertex(second).y));
    const auto original = compute_visibility_map(up);
    const auto explicit_perturbation = compute_visibility_map(perturbed);
    check_map(curve, original.visibility_map);
    check_map(perturbed_curve, explicit_perturbation.visibility_map);
    assert(keys(original.visibility_map, curve) ==
               keys(explicit_perturbation.visibility_map, perturbed_curve) &&
           "[C91 §2 tex 47]: exact height perturbation agrees with the symbolic visibility map");
}

void check_tree(std::size_t count, std::size_t shape) {
    Submap map;
    for (std::size_t i = 0; i < count; ++i) {
        map.add_node();
        Arc arc;
        arc.region_node = i;
        arc.edge_count = 1;
        map.add_arc(arc);
    }
    std::vector<std::size_t> available{0};
    chazelle::test::DeterministicRandomGenerator random(47);
    for (std::size_t i = 1; i < count; ++i) {
        const std::size_t candidate = random.next() % available.size();
        const std::size_t parent = shape == 0   ? i - 1
                                   : shape == 1 ? (i - 1) / 3
                                                : available[candidate];
        Chord chord;
        chord.y_tag = i;
        chord.region[0] = i % 2 == 0 ? parent : i;
        chord.region[1] = i % 2 == 0 ? i : parent;
        chord.left_adj = {{chord.region[0]}, 1};
        chord.right_adj = {{chord.region[1]}, 1};
        map.add_chord(chord);
        if (shape == 2) {
            if (map.node(parent).degree() == 4) {
                available[candidate] = available.back();
                available.pop_back();
            }
            available.push_back(i);
        }
    }
    map.build_tree_decomposition();
    const auto& tree = map.tree_decomposition();
    assert(tree.size() == 2 * count - 1);
    std::vector<std::size_t> first(tree.size());
    std::vector<std::size_t> last(tree.size());
    std::vector<std::size_t> leaves(count, NONE);
    std::size_t clock = 0;
    auto visit = [&](auto&& self, std::size_t node) -> std::size_t {
        first[node] = clock++;
        const auto& current = tree.node(node);
        std::size_t chords = 0;
        if (current.is_leaf()) {
            assert(leaves[current.region_idx] == NONE);
            leaves[current.region_idx] = node;
        } else {
            assert(tree.node(current.left_child).parent == node &&
                   tree.node(current.right_child).parent == node);
            const std::size_t left = self(self, current.left_child);
            const std::size_t right = self(self, current.right_child);
            chords = left + right + 1;
            const std::size_t bound = chords - (chords + 3) / 4;
            assert(left <= bound && right <= bound);
        }
        last[node] = clock;
        return chords;
    };
    assert(visit(visit, tree.root()) == count - 1);
    for (std::size_t i = 0; i < tree.size(); ++i) {
        const auto& node = tree.node(i);
        if (node.is_leaf())
            continue;
        const auto& chord = map.chord(node.chord_idx);
        for (std::size_t side = 0; side < 2; ++side) {
            const std::size_t child = side == 0 ? node.left_child : node.right_child;
            const std::size_t leaf = leaves[chord.region[side]];
            assert(first[child] <= first[leaf] && first[leaf] < last[child]);
        }
    }
}

void check_large_grade() {
    std::vector<Point> vertices;
    vertices.reserve(2049);
    for (std::size_t i = 0; i < 2049; ++i)
        vertices.push_back({0, Exact{i}, i + 23});
    std::puts("Preparing grade-11 up-phase");
    const UpPhase up(std::move(vertices));
    std::puts("Refining grade-11 submap");
    const auto result = compute_visibility_map(up);
    std::puts("Checking grade-11 result");
    assert(result.refinement_rounds == 1 && result.region_boundaries > 0);
    check_map(up.graded().curve(), result.visibility_map);
}

void check_nonroot_chains() {
    std::vector<Point> vertices;
    vertices.reserve(17);
    for (std::size_t i = 0; i < 17; ++i)
        vertices.push_back({Exact{i}, Exact{(13 * i) % 7}, i + 41});
    const UpPhase up(std::move(vertices));
    for (std::size_t grade = 0; grade < up.graded().maximum_grade(); ++grade)
        for (std::size_t index = 0; index < up.graded().num_chains(grade); ++index) {
            const auto result = compute_visibility_map(up, grade, index);
            check_map(up.graded().chain(grade, index), result.visibility_map);
        }
    vertices.clear();
    for (std::size_t i = 0; i < 129; ++i)
        vertices.push_back({Exact{i}, Exact{(13 * i) % 7}, i + 41});
    const UpPhase larger_up(std::move(vertices));
    const Polygon& curve = larger_up.graded().chain(6, 1);
    Submap coarse = larger_up.chain_submap(6, 1);
    enforce_granularity(coarse, curve, 4);
    coarse.normalize(curve);
    Submap refined = refine_visibility_submap(larger_up, curve, coarse, 6);
    assert(refined.is_conformal() && refined.is_granular(2, curve));
    auto full = complete_bounded_regions(curve, refined, 2).visibility_map;
    full.build_tree_decomposition();
    check_map(curve, full, true);
}

}

int main() {
    check_infinitesimals();
    check_boundary_views();
    for (std::size_t count : {std::size_t{1}, std::size_t{2}, std::size_t{5}, std::size_t{31},
                              std::size_t{2049}, std::size_t{32769}})
        for (std::size_t shape = 0; shape < 3; ++shape)
            check_tree(count, shape);
    const std::vector<std::vector<Point>> fixtures{
        {{0, 0, 0},
         {1, 2, 1},
         {2, 1, 2},
         {3, 3, 3},
         {4, 0, 4},
         {5, 2, 5},
         {6, -1, 6},
         {7, 3, 7},
         {8, 0, 8}},
        {{0, 0, 0},
         {1, 0, 1},
         {3, 0, 2},
         {4, 0, 3},
         {7, 0, 4},
         {8, 0, 5},
         {10, 0, 6},
         {11, 0, 7},
         {14, 0, 8}},
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
         {4, 8, 11}},
    };
    for (const auto& input : fixtures) {
        check_refinement(subdivide_curve(input, 64));
        for (std::size_t rotation = 0; rotation < 4; ++rotation) {
            for (bool reverse : {false, true}) {
                auto vertices = input;
                if (reverse)
                    std::reverse(vertices.begin(), vertices.end());
                for (std::size_t i = 0; i < vertices.size(); ++i) {
                    for (std::size_t turn = 0; turn < rotation; ++turn) {
                        const Exact x = vertices[i].x;
                        vertices[i].x = -vertices[i].y;
                        vertices[i].y = x;
                    }
                    vertices[i].index = i + 31;
                }
                check_visibility(std::move(vertices));
            }
        }
    }
    std::puts("Transformed down-phase fixtures passed");
    check_random_refinements();
    std::puts("Random full visibility maps passed");
    check_valid_refinements();
    std::puts("Refinements satisfying Lemma 4.2 passed");
    check_exhaustive_small_curves();
    std::puts("Exhaustive small visibility maps passed");
    check_refinement_grade_contract();
    check_exact_height_perturbation();
    std::puts("Exact height perturbation passed");
    check_nonroot_chains();
    std::puts("Nonroot down-phase chains passed");
    check_large_grade();
    std::puts("Section 4.2 down-phase tests passed");
}
