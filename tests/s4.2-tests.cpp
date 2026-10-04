#include "merge/granularity.h"
#include "submap/chord_inventory.h"
#include "support/random.h"
#include "visibility/bounded_regions.h"
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

using ChordKey = std::tuple<std::size_t, std::size_t, Side, std::size_t, Side, bool>;

std::vector<ChordKey> chord_keys(const Submap& submap, const Polygon& curve) {
    std::vector<ChordKey> keys;
    for (std::size_t i = 0; i < submap.num_chords(); ++i) {
        const Chord& c = submap.chord(i);
        assert(!c.dead);
        PendingChord canonical{c.symbolic_y(), c.left_edge,  c.left_side,
                               c.right_edge,   c.right_side, c.is_null_length};
        canonicalize_chord(canonical, curve);
        auto first = std::pair{canonical.left_edge_c, canonical.left_side};
        auto second = std::pair{canonical.right_edge_c, canonical.right_side};
        if (second < first)
            std::swap(first, second);
        keys.emplace_back(c.y_tag, first.first, first.second, second.first, second.second,
                          c.is_null_length);
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

void check_completion(const Polygon& curve, const Submap& coarse, std::size_t granularity) {
    const auto completed = complete_bounded_regions(curve, coarse, granularity);
    const Submap& full = completed.visibility_map;
    full.check_invariants(curve);
    const Submap reference = build_full_visibility_map(curve);
    const auto actual_keys = chord_keys(full, curve);
    auto reference_keys = chord_keys(reference, curve);
    const auto coarse_keys = chord_keys(coarse, curve);
    reference_keys.insert(reference_keys.end(), coarse_keys.begin(), coarse_keys.end());
    std::sort(reference_keys.begin(), reference_keys.end());
    reference_keys.erase(std::unique(reference_keys.begin(), reference_keys.end()),
                         reference_keys.end());
    if (actual_keys != reference_keys) {
        std::fprintf(stderr, "Chord mismatch: %zu vertices, first tag %zu\n", curve.num_vertices(),
                     curve.vertex(0).index);
        std::fprintf(stderr, "Actual %zu, expected %zu\n", actual_keys.size(),
                     reference_keys.size());
        for (const auto& key : reference_keys)
            if (std::find(actual_keys.begin(), actual_keys.end(), key) == actual_keys.end())
                std::fprintf(stderr, "Expected %zu: %zu/%u -- %zu/%u null %u\n", std::get<0>(key),
                             std::get<1>(key), std::get<2>(key), std::get<3>(key), std::get<4>(key),
                             std::get<5>(key));
        for (const auto& key : actual_keys)
            if (std::count(actual_keys.begin(), actual_keys.end(), key) !=
                std::count(reference_keys.begin(), reference_keys.end(), key))
                std::fprintf(stderr, "Extra %zu: %zu/%u -- %zu/%u null %u\n", std::get<0>(key),
                             std::get<1>(key), std::get<2>(key), std::get<3>(key), std::get<4>(key),
                             std::get<5>(key));
    }
    assert(actual_keys == reference_keys &&
           "[C91 Lemma 4.2 tex 367]: bounded-region completion recovers every full-map chord");
    assert(full.is_conformal() && full.is_semigranular(1));
    const std::size_t n = curve.num_vertices();
    const std::size_t m = coarse.num_arcs();
    assert(completed.work.vertex_occurrences <= 4 * n + 8 * m);
    assert(completed.work.height_comparisons <= 68 * completed.work.vertex_occurrences);
    assert(completed.work.ray_edge_tests <= 64 * n);
    assert(completed.work.endpoint_visits <= 12 * full.num_chords());
}

void check_curve(std::vector<Point> vertices) {
    const UpPhase up(std::move(vertices));
    for (std::size_t grade = 0; grade <= up.graded().maximum_grade(); ++grade) {
        assert(UpPhase::grade_granularity(grade) <= REGION_COMPLETION_GRANULARITY);
        for (std::size_t chain = 0; chain < up.graded().num_chains(grade); ++chain)
            check_completion(up.graded().chain(grade, chain), up.chain_submap(grade, chain),
                             UpPhase::grade_granularity(grade));
    }
    const std::size_t last = up.graded().curve().num_vertices() - 1;
    if (last >= 8) {
        const auto portion = up.compute_canonical_portion(1, last - 1);
        check_completion(portion.curve, portion.submap,
                         canonical_granularity(portion.curve.num_edges()));
    }
}

std::vector<std::vector<Point>> fixtures() {
    std::vector<std::vector<Point>> curves{
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
         {4, 8, 11}},
    };
    for (std::size_t length : {std::size_t{2}, std::size_t{18}, std::size_t{66}}) {
        std::vector<Point> horizontal;
        std::vector<Point> comb;
        for (std::size_t v = 0; v < length; ++v) {
            horizontal.push_back({Exact(v), 0, v});
            comb.push_back({Exact(v), Exact(v % 2 == 0 ? 0 : 9), v});
        }
        curves.push_back(std::move(horizontal));
        curves.push_back(std::move(comb));
    }
    return curves;
}

void check_height_sequences() {
    for (std::size_t sequence = 0; sequence < 243; ++sequence) {
        std::size_t heights = sequence;
        std::vector<Point> vertices;
        for (std::size_t v = 0; v < 5; ++v) {
            vertices.push_back({Exact(v), Exact(heights % 3), v + 31});
            heights /= 3;
        }
        check_curve(std::move(vertices));
    }
}

void check_endpoint_comparison() {
    const Polygon curve({{0, 0, 7}, {4, 3, 8}});
    std::vector<PendingChord> chords{{symbolic_y_of(curve.vertex(0)), 0, LEFT, 0, RIGHT, false}};
    const ChordEndpoint left{0, LEFT, chords[0].y, 0, true};
    const ChordEndpoint right{0, RIGHT, chords[0].y, 0, false};
    assert(!chord_endpoint_precedes(curve, chords, left, left));
    assert(!chord_endpoint_precedes(curve, chords, right, right));
    assert(chord_endpoint_precedes(curve, chords, left, right));
    assert(!chord_endpoint_precedes(curve, chords, right, left));

    const Polygon turn({{0, 0, 7}, {1, 3, 8}, {2, 0, 9}});
    chords = {{symbolic_y_of(turn.vertex(1)), 1, RIGHT, 1, RIGHT, true}};
    const ChordEndpoint null_left{1, RIGHT, chords[0].y, 0, true};
    const ChordEndpoint null_right{1, RIGHT, chords[0].y, 0, false};
    assert(!chord_endpoint_precedes(turn, chords, null_left, null_left));
    assert(!chord_endpoint_precedes(turn, chords, null_right, null_right));
    assert(chord_endpoint_precedes(turn, chords, null_left, null_right));
    assert(!chord_endpoint_precedes(turn, chords, null_right, null_left));
    chords.push_back(chords[0]);
    const std::vector<ChordEndpoint> copies{
        null_left, null_right, {1, RIGHT, chords[1].y, 1, true}, {1, RIGHT, chords[1].y, 1, false}};
    auto precedes = [&](const auto& a, const auto& b) {
        return chord_endpoint_precedes(turn, chords, a, b);
    };
    auto equivalent = [&](const auto& a, const auto& b) {
        return !precedes(a, b) && !precedes(b, a);
    };
    for (const auto& a : copies)
        for (const auto& b : copies)
            for (const auto& c : copies) {
                assert(!(precedes(a, b) && precedes(b, c)) || precedes(a, c));
                assert(!(equivalent(a, b) && equivalent(b, c)) || equivalent(a, c));
            }
}

void check_chordless_input() {
    for (const auto& vertices : std::vector<std::vector<Point>>{
             {{0, 0, 0}, {3, 5, 1}}, {{0, 0, 0}, {3, 5, 1}, {6, 0, 2}}}) {
        const Polygon curve(vertices);
        Submap coarse = build_full_visibility_map(curve);
        enforce_granularity(coarse, curve, REGION_COMPLETION_GRANULARITY);
        coarse.normalize(curve);
        assert(coarse.num_chords() == 0);
        check_completion(curve, coarse, REGION_COMPLETION_GRANULARITY);
    }
}

void check_random_curves() {
    chazelle::test::DeterministicRandomGenerator random(42);
    for (std::size_t sample = 0; sample < 30; ++sample) {
        std::vector<Point> vertices;
        const std::size_t n = 3 + random.next() % 96;
        vertices.reserve(n);
        for (std::size_t v = 0; v < n; ++v)
            vertices.push_back({Exact(v), Exact(random.next() % 21) / Exact(3), v});
        check_curve(std::move(vertices));
    }
}

void check_grade_ten() {
    std::vector<Point> vertices;
    vertices.reserve(1025);
    for (std::size_t v = 0; v < 1025; ++v)
        vertices.push_back({Exact(v), 0, v});
    const UpPhase up(std::move(vertices));
    assert(up.graded().maximum_grade() == 10);
    check_completion(up.graded().curve(), up.chain_submap(10, 0), UpPhase::grade_granularity(10));
}

}

int main() {
    check_endpoint_comparison();
    for (const auto& fixture : fixtures()) {
        for (std::size_t rotation = 0; rotation < 4; ++rotation) {
            for (const bool reverse : {false, true}) {
                auto vertices = fixture;
                if (reverse)
                    std::reverse(vertices.begin(), vertices.end());
                for (std::size_t i = 0; i < vertices.size(); ++i) {
                    auto& vertex = vertices[i];
                    for (std::size_t turn = 0; turn < rotation; ++turn) {
                        const Exact x = vertex.x;
                        vertex.x = -vertex.y;
                        vertex.y = x;
                    }
                    vertex.index = i + 17;
                }
                check_curve(std::move(vertices));
            }
        }
    }
    check_height_sequences();
    check_chordless_input();
    check_random_curves();
    check_grade_ten();
    std::puts("§4.2 bounded-region completion tests passed");
}
