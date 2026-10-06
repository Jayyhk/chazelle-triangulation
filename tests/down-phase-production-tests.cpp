#include "algorithm/submap/chord_inventory.h"
#include "algorithm/visibility/down_phase.h"
#include "support/random.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <tuple>
#include <vector>

using namespace chazelle;

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "%s\n", message);
        std::abort();
    }
}

using ChordKey = std::tuple<std::size_t, std::size_t, Side, std::size_t, Side, bool>;

std::vector<ChordKey> keys(const Submap& map, const Polygon& curve) {
    std::vector<ChordKey> result;
    result.reserve(map.num_chords());
    for (std::size_t i = 0; i < map.num_chords(); ++i) {
        const Chord& chord = map.chord(i);
        PendingChord value{chord.symbolic_y(), chord.left_edge,  chord.left_side,
                           chord.right_edge,   chord.right_side, chord.is_null_length};
        canonicalize_chord(value, curve);
        if (value.is_null_length)
            value.left_edge_c = value.right_edge_c = curve.local_index_of_tag(value.y.tag);
        const std::size_t source = curve.local_index_of_tag(value.y.tag);
        require(source != NONE && symbolic_y_equal(value.y, symbolic_y_of(curve.vertex(source))),
                "[C91 §2.1 tex 72]: every full-map chord retains its original source level");
        result.emplace_back(value.y.tag, value.left_edge_c, value.left_side, value.right_edge_c,
                            value.right_side, value.is_null_length);
    }
    std::sort(result.begin(), result.end());
    return result;
}

void check_large_curve(std::vector<Point> vertices) {
    const UpPhase up(std::move(vertices));
    const Polygon& curve = up.graded().curve();
    const auto result = compute_visibility_map(up);
    require(result.refinement_rounds == 1 && result.region_boundaries > 0,
            "[C91 Lemma 4.2 tex 364–381]: grade 11 must exercise inductive refinement");
    require(result.visibility_map.is_conformal() && result.visibility_map.is_semigranular(1) &&
                !result.visibility_map.tree_decomposition().empty(),
            "[C91 §§2.1/2.4]: V(C) is conformal and represented in normal form");
    const Submap reference = build_full_visibility_map(curve);
    require(keys(result.visibility_map, curve) == keys(reference, curve),
            "[C91 Theorem 4.3 tex 390]: the production algorithm must recover exactly V(C)");
    const auto& work = result.completion_work;
    require(work.vertex_occurrences <= 64 * curve.num_vertices() &&
                work.ray_edge_tests <= 2048 * curve.num_vertices() &&
                work.endpoint_visits <= 64 * curve.num_vertices(),
            "[C91 §4.2 tex 367]: bounded completion performs linear work");
}

}

int main() {
    chazelle::test::DeterministicRandomGenerator random(337);
    for (std::size_t sample = 0; sample < 4; ++sample) {
        std::vector<Point> vertices;
        vertices.reserve(2049);
        for (std::size_t i = 0; i < 2049; ++i)
            vertices.push_back({Exact{i} / 3, Exact{random.next() % 17} / 5, i});
        if (sample % 2 == 0)
            std::reverse(vertices.begin(), vertices.end());
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            if (sample % 2 == 0) {
                const Exact x = vertices[i].x;
                vertices[i].x = -vertices[i].y;
                vertices[i].y = x;
            } else if (sample == 1) {
                vertices[i].x = -vertices[i].x;
            }
            vertices[i].index = i + 101;
        }
        std::printf("Checking production grade-11 curve %zu\n", sample);
        std::fflush(stdout);
        check_large_curve(std::move(vertices));
    }
    std::puts("Production down-phase tests passed");
}
