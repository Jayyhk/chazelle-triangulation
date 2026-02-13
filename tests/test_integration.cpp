/// Submap + visibility data structure tests, plus end-to-end integration
/// tests feeding polygons through the full Chazelle pipeline.
///
///   - Submap: add_node, add_chord, remove_chord, weight computation
///   - Convert: convert_submap_to_fm
///   - End-to-end: polygon → up_phase → down_phase → convert → triangulate

#include "visibility/submap.h"
#include "visibility/arc_structure.h"
#include "visibility/chord.h"
#include "geometry/polygon.h"
#include "geometry/point.h"
#include "triangulate/convert.h"
#include "triangulate/triangulate.h"
#include "triangulate/vertex_node.h"
#include "chazelle/grade_storage.h"
#include "chazelle/up_phase.h"
#include "chazelle/down_phase.h"
#include "common.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <unordered_set>
#include <vector>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════════
//  Submap unit tests
// ════════════════════════════════════════════════════════════════════

static void test_submap_add_nodes() {
    std::printf("Test 1: Submap add_node\n");

    Submap sm;
    auto r0 = sm.add_node();
    auto r1 = sm.add_node();
    auto r2 = sm.add_node();

    assert(r0 == 0 && r1 == 1 && r2 == 2);
    assert(sm.num_nodes() == 3);
    assert(sm.node(0).degree() == 0);
    assert(sm.node(0).is_empty());
    assert(!sm.node(0).deleted);

    std::printf("  PASS\n\n");
}

static void test_submap_add_chord() {
    std::printf("Test 2: Submap add_chord\n");

    Submap sm;
    sm.add_node();
    sm.add_node();

    Chord c;
    c.region[0] = 0;
    c.region[1] = 1;
    c.y = 5.0;
    c.left_vertex = 3;
    c.right_vertex = 7;
    auto ci = sm.add_chord(c);

    assert(ci == 0);
    assert(sm.num_chords() == 1);
    assert(sm.node(0).degree() == 1);
    assert(sm.node(1).degree() == 1);
    assert(sm.node(0).incident_chords[0] == 0);
    assert(sm.chord(0).y == 5.0);

    std::printf("  PASS\n\n");
}

static void test_submap_remove_chord() {
    std::printf("Test 3: Submap remove_chord (merge regions)\n");

    Submap sm;
    sm.add_node();  // region 0
    sm.add_node();  // region 1
    sm.add_node();  // region 2

    // Chain: 0 —c0— 1 —c1— 2
    Chord c0;
    c0.region[0] = 0; c0.region[1] = 1;
    sm.add_chord(c0);

    Chord c1;
    c1.region[0] = 1; c1.region[1] = 2;
    sm.add_chord(c1);

    assert(sm.node(1).degree() == 2);

    // Remove chord 0 (merge regions 0 and 1 → survivor = 0).
    auto survivor = sm.remove_chord(0);
    assert(survivor == 0);
    assert(sm.node(1).deleted);

    // Chord 1 should now connect region 0 to region 2.
    assert(sm.chord(1).region[0] == 0 || sm.chord(1).region[1] == 0);
    assert(sm.chord(1).region[0] == 2 || sm.chord(1).region[1] == 2);

    // Region 0 should have degree 1 (only chord 1 remains).
    assert(sm.node(0).degree() == 1);

    std::printf("  PASS\n\n");
}

static void test_submap_add_arc_and_weight() {
    std::printf("Test 4: Submap arc + weight computation\n");

    Submap sm;
    sm.add_node();  // region 0

    ArcStructure a;
    a.first_edge = 3;
    a.last_edge = 7;
    a.region_node = 0;
    a.edge_count = 5;
    auto ai = sm.add_arc(a);
    sm.node(0).arcs.push_back(ai);

    ArcStructure b;
    b.first_edge = 10;
    b.last_edge = 12;
    b.region_node = 0;
    b.edge_count = 3;
    auto bi = sm.add_arc(b);
    sm.node(0).arcs.push_back(bi);

    sm.recompute_weight(0);
    // Weight = max edge_count = max(5, 3) = 5.
    assert(sm.node(0).weight == 5);
    assert(!sm.node(0).is_empty());

    std::printf("  PASS\n\n");
}

static void test_submap_is_conformal() {
    std::printf("Test 5: Submap conformality check\n");

    Submap sm;
    sm.add_node();  // degree will be 5 → not conformal

    // Add 5 neighbors, each connected by a chord.
    for (int i = 0; i < 5; ++i) {
        sm.add_node();
        Chord c;
        c.region[0] = 0;
        c.region[1] = static_cast<std::size_t>(i + 1);
        sm.add_chord(c);
    }

    assert(sm.node(0).degree() == 5);
    assert(!sm.is_conformal()); // degree > 4

    // Remove one chord to bring degree to 4.
    sm.remove_chord(4); // merges region 5 into region 0
    // Now region 0 has chords 0,1,2,3 → degree 4.
    assert(sm.is_conformal());

    std::printf("  PASS\n\n");
}

static void test_submap_semigranularity() {
    std::printf("Test 6: Submap semigranularity check\n");

    Submap sm;
    sm.add_node();

    ArcStructure a;
    a.edge_count = 10;
    a.region_node = 0;
    auto ai = sm.add_arc(a);
    sm.node(0).arcs.push_back(ai);
    sm.recompute_weight(0);

    assert(sm.is_semigranular(10));   // weight ≤ γ=10
    assert(sm.is_semigranular(100));  // weight ≤ γ=100
    assert(!sm.is_semigranular(9));   // weight > γ=9

    std::printf("  PASS\n\n");
}

static void test_submap_remove_transfers_arcs() {
    std::printf("Test 7: remove_chord transfers arcs to survivor\n");

    Submap sm;
    sm.add_node();  // 0
    sm.add_node();  // 1

    // Region 0 has arc A, region 1 has arc B.
    ArcStructure arcA;
    arcA.edge_count = 3; arcA.region_node = 0;
    auto aiA = sm.add_arc(arcA);
    sm.node(0).arcs.push_back(aiA);

    ArcStructure arcB;
    arcB.edge_count = 7; arcB.region_node = 1;
    auto aiB = sm.add_arc(arcB);
    sm.node(1).arcs.push_back(aiB);

    Chord c;
    c.region[0] = 0; c.region[1] = 1;
    sm.add_chord(c);

    sm.recompute_weight(0);
    sm.recompute_weight(1);
    assert(sm.node(0).weight == 3);
    assert(sm.node(1).weight == 7);

    // Remove the chord → merge into region 0.
    sm.remove_chord(0);

    // Region 0 should now have both arcs.
    assert(sm.node(0).arcs.size() == 2);
    // Weight should be max(3, 7) = 7.
    assert(sm.node(0).weight == 7);
    // Arc B's region_node should be updated.
    assert(sm.arc(aiB).region_node == 0);

    std::printf("  PASS\n\n");
}

// ════════════════════════════════════════════════════════════════════
//  ArcStructure tests
// ════════════════════════════════════════════════════════════════════

static void test_arc_structure_queries() {
    std::printf("Test 8: ArcStructure queries\n");

    ArcStructure null_arc;
    null_arc.first_edge = 5;
    null_arc.last_edge = 5;
    null_arc.edge_count = 0;
    assert(null_arc.is_null());
    assert(null_arc.is_single_edge());

    ArcStructure single;
    single.first_edge = 3;
    single.last_edge = 3;
    single.edge_count = 1;
    assert(!single.is_null());
    assert(single.is_single_edge());

    ArcStructure multi;
    multi.first_edge = 3;
    multi.last_edge = 7;
    multi.edge_count = 5;
    assert(!multi.is_null());
    assert(!multi.is_single_edge());

    std::printf("  PASS\n\n");
}

// ════════════════════════════════════════════════════════════════════
//  End-to-end integration: full Chazelle pipeline
// ════════════════════════════════════════════════════════════════════

/// Helper: run the full pipeline on a polygon and verify triangle count.
static void run_full_pipeline(const char* name,
                              std::vector<Point> vertices,
                              bool expect_success = true) {
    std::printf("  %s (%zu vertices)...", name, vertices.size());
    std::fflush(stdout);

    std::size_t n = vertices.size();

    try {
        Polygon polygon(std::move(vertices));

        GradeStorage storage;
        fprintf(stderr, "  [%s] up_phase...\n", name);
        up_phase(polygon, storage);
        fprintf(stderr, "  [%s] down_phase...\n", name);
        Submap vp = down_phase(polygon, storage);
        fprintf(stderr, "  [%s] convert_submap_to_fm...\n", name);
        auto fm = convert_submap_to_fm(vp, polygon);
        fprintf(stderr, "  [%s] fm_triangulate...\n", name);
        auto triangles = fm_triangulate(fm.nodes, fm.traps, fm.first, fm.last);
        fprintf(stderr, "  [%s] done (%zu triangles)\n", name, triangles.size());

        // A simple polygon with n vertices → exactly n-2 triangles.
        if (triangles.size() != n - 2) {
            std::printf(" FAIL (expected %zu triangles, got %zu)\n",
                        n - 2, triangles.size());
            // Debug: print visibility map info.
            std::printf("    VP: %zu nodes, %zu chords, %zu arcs\n",
                        vp.num_nodes(), vp.num_chords(), vp.num_arcs());
            std::printf("    FM: %zu nodes, %zu traps\n",
                        fm.nodes.size(), fm.traps.size());
            for (std::size_t i = 0; i < vp.num_chords(); ++i) {
                const auto& c = vp.chord(i);
                std::printf("    chord %zu: y=%.3f lv=%zu rv=%zu r0=%zu r1=%zu null=%d\n",
                            i, c.y, c.left_vertex, c.right_vertex,
                            c.region[0], c.region[1], (int)c.is_null_length);
            }
            for (std::size_t i = 0; i < vp.num_nodes(); ++i) {
                const auto& nd = vp.node(i);
                if (nd.deleted) continue;
                std::printf("    node %zu: %zu arcs, %zu ichords, w=%zu\n",
                            i, nd.arcs.size(), nd.incident_chords.size(), nd.weight);
                for (std::size_t ai : nd.arcs) {
                    const auto& a = vp.arc(ai);
                    std::printf("      arc %zu: edges [%zu,%zu]\n",
                                ai, a.first_edge, a.last_edge);
                }
            }
            for (auto& [a, b, c] : triangles) {
                std::printf("    tri: %zu %zu %zu\n", a, b, c);
            }
            assert(false);
        }

        // All triangle vertex indices should be in [0, n).
        for (auto& [a, b, c] : triangles) {
            assert(a < n && b < n && c < n);
            assert(a != b && b != c && a != c);
        }

        // All original vertices should appear in at least one triangle.
        std::unordered_set<std::size_t> used;
        for (auto& [a, b, c] : triangles) {
            used.insert(a);
            used.insert(b);
            used.insert(c);
        }
        for (std::size_t i = 0; i < n; ++i) {
            assert(used.count(i) > 0);
        }

        std::printf(" PASS (%zu triangles)\n", triangles.size());
    } catch (const std::exception& ex) {
        if (expect_success) {
            std::printf(" FAIL (exception: %s)\n", ex.what());
            assert(false);
        } else {
            std::printf(" OK (expected failure: %s)\n", ex.what());
        }
    }
}

static void test_e2e_triangle() {
    std::printf("Test 9: end-to-end triangle\n");
    run_full_pipeline("triangle", {
        {0, 0, 0}, {4, 0, 1}, {2, 3, 2}
    });
    std::printf("\n");
}

static void test_e2e_square() {
    std::printf("Test 10: end-to-end square\n");
    run_full_pipeline("square", {
        {0, 0, 0}, {4, 0, 1}, {4, 4, 2}, {0, 4, 3}
    });
    std::printf("\n");
}

static void test_e2e_pentagon() {
    std::printf("Test 11: end-to-end pentagon\n");
    run_full_pipeline("pentagon", {
        {4, 0, 0}, {5, 2, 1}, {2.5, 4, 2}, {0, 2, 3}, {1, 0, 4}
    });
    std::printf("\n");
}

static void test_e2e_hexagon() {
    std::printf("Test 12: end-to-end hexagon\n");
    run_full_pipeline("hexagon", {
        {2, 0, 0}, {4, 1, 1}, {4, 3, 2}, {2, 4, 3}, {0, 3, 4}, {0, 1, 5}
    });
    std::printf("\n");
}

static void test_e2e_concave_L() {
    std::printf("Test 13: end-to-end L-shape (concave)\n");

    // An L-shaped polygon (6 vertices, one reflex corner):
    //
    //   5 ──── 4
    //   |      |
    //   |  3 ──2
    //   |  |
    //   0 ──1
    run_full_pipeline("L-shape", {
        {0, 0, 0}, {2, 0, 1}, {2, 2, 2}, {4, 2, 3}, {4, 4, 4}, {0, 4, 5}
    });
    std::printf("\n");
}

static void test_e2e_regular_polygons() {
    std::printf("Test 14: end-to-end regular convex polygons\n");

    for (std::size_t n : {7u, 8u, 9u, 10u, 12u, 16u}) {
        std::vector<Point> pts(n);
        for (std::size_t i = 0; i < n; ++i) {
            double angle = 2.0 * M_PI * static_cast<double>(i) /
                           static_cast<double>(n);
            pts[i] = {10.0 * std::cos(angle), 10.0 * std::sin(angle), i};
        }
        char label[64];
        std::snprintf(label, sizeof(label), "regular %zu-gon", n);
        run_full_pipeline(label, pts);
    }
    std::printf("\n");
}

static void test_e2e_star() {
    std::printf("Test 15: end-to-end star polygon\n");

    const std::size_t N = 10;
    std::vector<Point> pts(N);
    for (std::size_t i = 0; i < N; ++i) {
        double angle = 2.0 * M_PI * static_cast<double>(i) /
                       static_cast<double>(N);
        double r = (i % 2 == 0) ? 3.0 : 1.2;
        pts[i] = {r * std::cos(angle), r * std::sin(angle), i};
    }
    run_full_pipeline("star-10", pts);
    std::printf("\n");
}

static void test_e2e_comb() {
    std::printf("Test 16: end-to-end comb polygon\n");

    // A comb: a base rectangle with teeth extending upward.
    //
    //   *  *  *  *        ← teeth tips
    //   |  |  |  |
    //   *--*  *--*  *--*  ← teeth bases
    //   |              |
    //   *──────────────*  ← base
    //
    // 4 teeth, each 2 units wide, 2 units tall, gaps of 1.
    // Slightly jitter gap-base y-coordinates to avoid collinearity.
    std::vector<Point> pts;
    std::size_t idx = 0;

    // Bottom-left to bottom-right (base).
    pts.push_back({0, 0, idx++});
    pts.push_back({13, 0, idx++});

    // Right side up.
    pts.push_back({13, 2, idx++});

    // Teeth (right to left).
    double x = 12;
    for (int tooth = 0; tooth < 4; ++tooth) {
        // Tooth right edge top.
        pts.push_back({x, 4, idx++});
        // Tooth left edge top.
        pts.push_back({x - 2, 4, idx++});
        // Tooth left edge base.
        pts.push_back({x - 2, 2.01 + 0.001 * tooth, idx++});

        if (tooth < 3) {
            // Gap base.
            pts.push_back({x - 3, 2.02 + 0.001 * tooth, idx++});
        }
        x -= 3;
    }

    // Left side down.
    pts.push_back({0, 2, idx++});

    run_full_pipeline("comb", pts);
    std::printf("\n");
}

static void test_e2e_zigzag_monotone() {
    std::printf("Test 17: end-to-end tall zigzag (y-monotone)\n");

    // A y-monotone zigzag.
    std::vector<Point> pts;
    std::size_t idx = 0;
    const std::size_t teeth = 6;

    // Right chain going up.
    for (std::size_t i = 0; i <= teeth; ++i) {
        double y = static_cast<double>(i) * 2.0;
        double x = (i % 2 == 0) ? 4.0 : 6.0;
        pts.push_back({x, y, idx++});
    }
    // Left chain going down.
    for (std::size_t i = teeth; i > 0; --i) {
        double y = static_cast<double>(i) * 2.0 - 1.0;
        double x = (i % 2 == 0) ? 0.0 : 2.0;
        pts.push_back({x, y, idx++});
    }

    run_full_pipeline("zigzag-monotone", pts);
    std::printf("\n");
}

// ════════════════════════════════════════════════════════════════════
//  Granularity constant check
// ════════════════════════════════════════════════════════════════════

static void test_compute_granularity() {
    std::printf("Test 18: compute_granularity values\n");

    // γ = 2^⌈λ/5⌉
    assert(compute_granularity(0) == 1);   // 2^⌈0/5⌉ = 2^0 = 1
    assert(compute_granularity(1) == 2);   // 2^⌈1/5⌉ = 2^1 = 2
    assert(compute_granularity(5) == 2);   // 2^⌈5/5⌉ = 2^1 = 2
    assert(compute_granularity(6) == 4);   // 2^⌈6/5⌉ = 2^2 = 4
    assert(compute_granularity(10) == 4);  // 2^⌈10/5⌉ = 2^2 = 4
    assert(compute_granularity(11) == 8);  // 2^⌈11/5⌉ = 2^3 = 8
    assert(compute_granularity(20) == 16); // 2^⌈20/5⌉ = 2^4 = 16

    std::printf("  PASS\n\n");
}

// ════════════════════════════════════════════════════════════════════
//  main
// ════════════════════════════════════════════════════════════════════

int main() {
    std::printf("=== Submap + integration tests ===\n\n");

    test_submap_add_nodes();
    test_submap_add_chord();
    test_submap_remove_chord();
    test_submap_add_arc_and_weight();
    test_submap_is_conformal();
    test_submap_semigranularity();
    test_submap_remove_transfers_arcs();
    test_arc_structure_queries();

    std::printf("--- End-to-end pipeline ---\n\n");
    test_e2e_triangle();
    test_e2e_square();
    test_e2e_pentagon();
    test_e2e_hexagon();
    test_e2e_concave_L();
    test_e2e_regular_polygons();
    test_e2e_star();
    test_e2e_comb();
    test_e2e_zigzag_monotone();
    test_compute_granularity();

    std::printf("All %d tests passed.\n", 18);
    return 0;
}
