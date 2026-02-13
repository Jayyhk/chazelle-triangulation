/// Comprehensive tests for the Fournier-Montuno triangulation module.
///
/// Tests Algorithm 3 (triangulate_monotone) and Algorithm 2+3 (fm_triangulate)
/// on a variety of polygon shapes and configurations.

#include "triangulate/trapezoid.h"
#include "triangulate/vertex_node.h"
#include "triangulate/monotone.h"
#include "triangulate/triangulate.h"
#include "geometry/point.h"
#include "geometry/perturbation.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace chazelle;

// ── Helpers ─────────────────────────────────────────────────────────

/// Print a triangle list.
static void print_triangles(const char* label,
                            const std::vector<Triangle>& tris) {
    std::printf("  %s (%zu triangles):\n", label, tris.size());
    for (auto& [a, b, c] : tris)
        std::printf("    (%zu, %zu, %zu)\n", a, b, c);
}

/// Validate that every triangle has three distinct vertex indices
/// drawn from [0, n).
static void check_valid_triangles(const std::vector<Triangle>& tris,
                                  std::size_t n) {
    for (auto& [a, b, c] : tris) {
        assert(a < n && b < n && c < n);
        assert(a != b && b != c && a != c);
    }
}

/// Validate that no edge (undirected) appears in more than 2 triangles
/// (each internal diagonal is shared by exactly 2, boundary edges by 1).
static void check_no_overlapping_triangles(const std::vector<Triangle>& tris) {
    // Count how many triangles each edge appears in.
    struct PairHash {
        std::size_t operator()(std::pair<std::size_t, std::size_t> p) const {
            return std::hash<std::size_t>()(p.first) * 31 +
                   std::hash<std::size_t>()(p.second);
        }
    };
    std::unordered_map<std::pair<std::size_t, std::size_t>, int, PairHash> edge_count;
    for (auto& [a, b, c] : tris) {
        auto e = [](std::size_t u, std::size_t v) {
            return std::make_pair(std::min(u, v), std::max(u, v));
        };
        edge_count[e(a, b)]++;
        edge_count[e(b, c)]++;
        edge_count[e(a, c)]++;
    }
    for (auto& [edge, count] : edge_count) {
        assert(count <= 2);
    }
}

/// Check that all original vertices appear in at least one triangle.
static void check_all_vertices_used(const std::vector<Triangle>& tris,
                                    std::size_t n) {
    std::unordered_set<std::size_t> used;
    for (auto& [a, b, c] : tris) {
        used.insert(a);
        used.insert(b);
        used.insert(c);
    }
    for (std::size_t i = 0; i < n; ++i) {
        assert(used.count(i) > 0);
    }
}

/// Build a regular convex polygon with n vertices centered at origin.
static std::vector<Point> regular_polygon(std::size_t n) {
    std::vector<Point> pts(n);
    for (std::size_t i = 0; i < n; ++i) {
        double angle = 2.0 * M_PI * static_cast<double>(i) /
                       static_cast<double>(n);
        pts[i] = {std::cos(angle), std::sin(angle), i};
    }
    return pts;
}

// ════════════════════════════════════════════════════════════════════
//  Algorithm 3: triangulate_monotone
// ════════════════════════════════════════════════════════════════════

static void test_triangle() {
    std::printf("Test 1: triangle (3 vertices)\n");

    std::vector<Point> pts = {
        {0.0, 0.0, 0}, {4.0, 0.0, 1}, {2.0, 3.0, 2}
    };
    auto nodes = build_vertex_list(pts);
    std::vector<Triangle> out;
    triangulate_monotone(nodes, 0, 2, out);
    assert(out.size() == 1);
    check_valid_triangles(out, 3);
    std::printf("  PASS\n\n");
}

static void test_monotone_quad() {
    std::printf("Test 2: monotone quadrilateral (4 vertices)\n");

    std::vector<Point> pts = {
        {2.0, 0.0, 0}, {4.0, 1.0, 1}, {2.0, 3.0, 2}, {0.0, 1.0, 3}
    };
    auto nodes = build_vertex_list(pts);
    std::vector<Triangle> out;
    triangulate_monotone(nodes, 0, 2, out);
    assert(out.size() == 2);
    check_valid_triangles(out, 4);
    check_no_overlapping_triangles(out);
    std::printf("  PASS\n\n");
}

static void test_monotone_tall_zigzag() {
    std::printf("Test 3: tall y-monotone zigzag (8 vertices)\n");

    // A y-monotone polygon shaped like a zigzag.
    // Left chain descends, right chain descends.
    //
    //     0 (top, y=8)
    //    / |
    //   7   1   (y=7 and y=7)
    //  /     |
    // 6       2  (y=5 and y=5)
    //  |     /
    //   5   3   (y=3 and y=3)
    //    | /
    //     4 (bottom, y=0)
    // CCW order: 0(top) → 1(right-top) → ... → 4(bottom) → ... → 7(left-top)
    // is CW.  Reverse to get CCW:
    //   0(top) → 7(left-top) → 6(left) → 5(left-bottom) →
    //   4(bottom) → 3(right-bottom) → 2(right) → 1(right-top)
    std::vector<Point> pts = {
        {4.0, 8.0, 0},   // top
        {2.0, 7.0, 1},   // left-top
        {1.0, 5.0, 2},   // left
        {2.0, 3.0, 3},   // left-bottom
        {4.0, 0.0, 4},   // bottom
        {6.0, 3.0, 5},   // right-bottom
        {7.0, 5.0, 6},   // right
        {6.0, 7.0, 7},   // right-top
    };

    auto nodes = build_vertex_list(pts);

    // This is a monotone polygon. Treat as unimonotone with base 0→4.
    // Need a trapezoid to create the base diagonal.
    std::vector<Trapezoid> traps;
    traps.push_back({0, 4, FM_NONE, FM_NONE});
    nodes[0].trapezoid_idx = 0;

    auto out = fm_triangulate(nodes, traps, 0, 7);
    print_triangles("zigzag", out);
    assert(out.size() == 6); // 8 vertices → 6 triangles
    check_valid_triangles(out, 8);
    check_no_overlapping_triangles(out);
    check_all_vertices_used(out, 8);
    std::printf("  PASS\n\n");
}

// ════════════════════════════════════════════════════════════════════
//  Algorithm 2+3: fm_triangulate (non-monotone with trapezoids)
// ════════════════════════════════════════════════════════════════════

static void test_fm_pentagon() {
    std::printf("Test 4: fm_triangulate on pentagon (5 vertices)\n");

    std::vector<Point> pts = {
        {4.0, 0.0, 0}, {5.0, 2.0, 1}, {2.5, 4.0, 2},
        {0.0, 2.0, 3}, {1.0, 0.0, 4}
    };
    auto nodes = build_vertex_list(pts);
    std::vector<Trapezoid> traps;
    traps.push_back({2, 0, FM_NONE, FM_NONE});
    nodes[2].trapezoid_idx = 0;

    auto out = fm_triangulate(nodes, traps, 0, 4);
    print_triangles("pentagon", out);
    assert(out.size() == 3);
    check_valid_triangles(out, 5);
    check_no_overlapping_triangles(out);
    check_all_vertices_used(out, 5);
    std::printf("  PASS\n\n");
}

static void test_fm_hexagon() {
    std::printf("Test 5: fm_triangulate on hexagon (6 vertices)\n");

    std::vector<Point> pts = {
        {2.0, 5.0, 0}, {0.0, 4.0, 1}, {0.0, 2.0, 2},
        {2.0, 0.0, 3}, {4.0, 2.0, 4}, {4.0, 4.0, 5}
    };
    auto nodes = build_vertex_list(pts);
    std::vector<Trapezoid> traps;
    traps.push_back({0, 3, FM_NONE, FM_NONE});
    nodes[0].trapezoid_idx = 0;

    auto out = fm_triangulate(nodes, traps, 0, 5);
    print_triangles("hexagon", out);
    assert(out.size() == 4);
    check_valid_triangles(out, 6);
    check_no_overlapping_triangles(out);
    check_all_vertices_used(out, 6);
    std::printf("  PASS\n\n");
}

static void test_fm_two_trapezoids() {
    std::printf("Test 6: fm_triangulate with 2 Class-B trapezoids (7 vertices)\n");

    // Heptagon-like shape with two non-adjacent trapezoid diagonals.
    //
    //       1 (y=5)
    //      / |
    //     2   0 (y=4 and y=3)
    //    /     |
    //   3       6 (y=2 and y=2)
    //    |     /
    //     4   5 (y=1 and y=1)
    //      | /
    //   (connect 4->5)
    std::vector<Point> pts = {
        {5.0, 3.0, 0},
        {3.0, 5.0, 1},
        {1.0, 4.0, 2},
        {0.0, 2.0, 3},
        {1.0, 0.0, 4},
        {4.0, 0.0, 5},
        {5.0, 2.0, 6},
    };

    auto nodes = build_vertex_list(pts);
    std::vector<Trapezoid> traps;

    // Trapezoid 0: connects vertex 1 (top) to vertex 4 (not adjacent → Class B)
    traps.push_back({1, 4, FM_NONE, FM_NONE});
    nodes[1].trapezoid_idx = 0;

    // Trapezoid 1: connects vertex 1 (top) to vertex 5 (not adjacent → Class B)
    traps.push_back({1, 5, FM_NONE, FM_NONE});
    nodes[5].trapezoid_idx = 1;

    auto out = fm_triangulate(nodes, traps, 0, 6);
    print_triangles("two-trap heptagon", out);
    assert(out.size() == 5); // 7 vertices → 5 triangles
    check_valid_triangles(out, 7);
    check_no_overlapping_triangles(out);
    check_all_vertices_used(out, 7);
    std::printf("  PASS\n\n");
}

static void test_fm_concave_arrow() {
    std::printf("Test 7: concave arrow shape (5 vertices, concavity)\n");

    // Arrow pointing right with a concave dent at vertex 2:
    //
    //   0 (0,4) ------ 1 (4,4)
    //     |                /
    //       2 (2,2)       /    <- concave vertex
    //     /              /
    //   4 (0,0) ------ 3 (4,0)
    //
    // CCW order
    std::vector<Point> pts = {
        {0.0, 4.0, 0},
        {0.0, 0.0, 1},
        {4.0, 0.0, 2},
        {2.0, 2.0, 3},  // reflex vertex
        {4.0, 4.0, 4},
    };

    auto nodes = build_vertex_list(pts);
    std::vector<Trapezoid> traps;
    // Trapezoid connecting top-level vertices across the concavity.
    traps.push_back({0, 3, FM_NONE, FM_NONE});
    nodes[0].trapezoid_idx = 0;

    auto out = fm_triangulate(nodes, traps, 0, 4);
    print_triangles("concave arrow", out);
    assert(out.size() == 3);
    check_valid_triangles(out, 5);
    check_no_overlapping_triangles(out);
    check_all_vertices_used(out, 5);
    std::printf("  PASS\n\n");
}

static void test_fm_class_a_only() {
    std::printf("Test 8: all Class-A trapezoids (no diagonals needed)\n");

    // A unimonotone polygon where all trapezoids are Class A
    // (adjacent vertices).  fm_triangulate should fall straight through
    // to triangulate_monotone.
    //
    //   0 (top)
    //   |
    //   1
    //   |
    //   2 (bottom)
    // CCW triangle
    std::vector<Point> pts = {
        {0.0, 4.0, 0},
        {0.0, 0.0, 1},
        {2.0, 2.0, 2},
    };
    auto nodes = build_vertex_list(pts);
    std::vector<Trapezoid> traps;
    // Class A: top=0, bottom=1 — they ARE adjacent (0→1).
    traps.push_back({0, 1, FM_NONE, FM_NONE});
    nodes[0].trapezoid_idx = 0;

    auto out = fm_triangulate(nodes, traps, 0, 2);
    assert(out.size() == 1);
    check_valid_triangles(out, 3);
    std::printf("  PASS\n\n");
}

static void test_fm_octagon_multi_split() {
    std::printf("Test 9: regular octagon with multiple trapezoid splits\n");

    auto pts = regular_polygon(8);
    auto nodes = build_vertex_list(pts);
    std::vector<Trapezoid> traps;

    // Create non-crossing Class-B trapezoid diagonals.
    // For a regular octagon (CCW, starting at angle 0):
    //   0:(1,0)  1:(0.7,0.7)  2:(0,1)  3:(-0.7,0.7)
    //   4:(-1,0) 5:(-0.7,-0.7) 6:(0,-1) 7:(0.7,-0.7)
    //
    // Use a fan from vertex 0: diagonals 0→3 and 0→5.
    // These are non-crossing and split the octagon into 3 pieces:
    //   0→1→2→3→0,  0→3→4→5→0,  0→5→6→7→0
    traps.push_back({0, 3, FM_NONE, FM_NONE});
    nodes[0].trapezoid_idx = 0;

    traps.push_back({0, 5, FM_NONE, FM_NONE});
    nodes[5].trapezoid_idx = 1;

    auto out = fm_triangulate(nodes, traps, 0, 7);
    print_triangles("octagon", out);
    assert(out.size() == 6); // 8 vertices → 6 triangles
    check_valid_triangles(out, 8);
    check_no_overlapping_triangles(out);
    check_all_vertices_used(out, 8);
    std::printf("  PASS\n\n");
}

static void test_fm_large_polygon() {
    std::printf("Test 10: regular 20-gon stress test\n");

    const std::size_t N = 20;
    auto pts = regular_polygon(N);
    auto nodes = build_vertex_list(pts);
    std::vector<Trapezoid> traps;

    // Single non-adjacent trapezoid diagonal: 0→10.
    traps.push_back({0, 10, FM_NONE, FM_NONE});
    nodes[0].trapezoid_idx = 0;

    auto out = fm_triangulate(nodes, traps, 0, N - 1);
    print_triangles("20-gon", out);
    assert(out.size() == N - 2);
    check_valid_triangles(out, N);
    check_no_overlapping_triangles(out);
    check_all_vertices_used(out, N);
    std::printf("  PASS\n\n");
}

// ════════════════════════════════════════════════════════════════════
//  Algorithm 3 edge cases
// ════════════════════════════════════════════════════════════════════

static void test_monotone_thin_spike() {
    std::printf("Test 11: thin spike monotone polygon\n");

    // Very thin polygon (nearly degenerate):
    //   0 (0,0)  1 (10,0)  2 (5, 0.001)
    // This is a valid (very flat) triangle.
    std::vector<Point> pts = {
        {0.0, 0.0, 0}, {10.0, 0.0, 1}, {5.0, 0.001, 2}
    };
    auto nodes = build_vertex_list(pts);
    std::vector<Triangle> out;
    triangulate_monotone(nodes, 0, 2, out);
    assert(out.size() == 1);
    check_valid_triangles(out, 3);
    std::printf("  PASS\n\n");
}

static void test_monotone_5_vertices() {
    std::printf("Test 12: monotone pentagon (5 vertices via triangulate_monotone)\n");

    // A y-monotone polygon where one chain is 3 vertices and the other
    // is the base edge.
    //
    //    2 (top, y=4)
    //   / |
    //  1   3
    //   | /
    //    0 (bottom, y=0)  <- base edge 2->0
    //   /
    //  4 (y=1, left)
    //
    // Wait — that creates a crossing.  Let's use a proper monotone shape:
    //
    //    2 (2, 4)        <- top
    //   / |
    //  3   1             <- left and right (y=2)
    //   | /
    //    0 (2, 0)        <- bottom
    //    |
    //    4 (2, -1)       ← below bottom on same x
    //
    // Actually let's just make a simple 5-vertex monotone:
    //   top at vertex 0, bottom at vertex 2
    //   Right chain: 0→1→2
    //   Left chain:  0→4→3→2
    // CCW order
    std::vector<Point> pts = {
        {2.0, 5.0, 0},   // top
        {0.0, 4.0, 1},   // left-high
        {0.0, 1.0, 2},   // left-low
        {3.0, 0.0, 3},   // bottom
        {4.0, 3.0, 4},   // right
    };

    auto nodes = build_vertex_list(pts);
    // Make it unimonotone with base 0→3.
    std::vector<Trapezoid> traps;
    traps.push_back({0, 3, FM_NONE, FM_NONE});
    nodes[0].trapezoid_idx = 0;

    auto out = fm_triangulate(nodes, traps, 0, 4);
    print_triangles("monotone pentagon", out);
    assert(out.size() == 3);
    check_valid_triangles(out, 5);
    check_no_overlapping_triangles(out);
    check_all_vertices_used(out, 5);
    std::printf("  PASS\n\n");
}

static void test_no_trapezoids() {
    std::printf("Test 13: unimonotone polygon with no trapezoids at all\n");

    // A triangle with no trapezoid records — algo2 should immediately
    // fall through to triangulate_monotone.
    std::vector<Point> pts = {
        {0.0, 0.0, 0}, {4.0, 0.0, 1}, {2.0, 3.0, 2}
    };
    auto nodes = build_vertex_list(pts);
    std::vector<Trapezoid> traps; // empty

    auto out = fm_triangulate(nodes, traps, 0, 2);
    assert(out.size() == 1);
    check_valid_triangles(out, 3);
    std::printf("  PASS\n\n");
}

// ════════════════════════════════════════════════════════════════════
//  VertexNode / build_vertex_list
// ════════════════════════════════════════════════════════════════════

static void test_build_vertex_list() {
    std::printf("Test 14: build_vertex_list circularity\n");

    std::vector<Point> pts = {
        {0, 0, 10}, {1, 0, 11}, {2, 0, 12}, {3, 0, 13}
    };
    auto nodes = build_vertex_list(pts);

    assert(nodes.size() == 4);

    // Check circular links.
    assert(nodes[0].prev == 3);
    assert(nodes[0].next == 1);
    assert(nodes[3].prev == 2);
    assert(nodes[3].next == 0);

    // Check coordinate/index copying.
    assert(nodes[0].index == 10);
    assert(nodes[0].x == 0.0);
    assert(nodes[2].index == 12);

    // All trapezoid_idx should be FM_NONE.
    for (auto& n : nodes)
        assert(n.trapezoid_idx == FM_NONE);

    std::printf("  PASS\n\n");
}

// ════════════════════════════════════════════════════════════════════
//  Larger stress tests
// ════════════════════════════════════════════════════════════════════

static void test_fm_50gon() {
    std::printf("Test 15: regular 50-gon with fan of trapezoid diagonals\n");

    const std::size_t N = 50;
    auto pts = regular_polygon(N);
    auto nodes = build_vertex_list(pts);
    std::vector<Trapezoid> traps;

    // Single non-adjacent diagonal: 0→25 (opposite side).
    traps.push_back({0, 25, FM_NONE, FM_NONE});
    nodes[0].trapezoid_idx = 0;

    auto out = fm_triangulate(nodes, traps, 0, N - 1);
    assert(out.size() == N - 2);
    check_valid_triangles(out, N);
    check_no_overlapping_triangles(out);
    check_all_vertices_used(out, N);
    std::printf("  PASS\n\n");
}

static void test_fm_star_polygon() {
    std::printf("Test 16: star-shaped polygon (10 vertices, alternating radii)\n");

    // Star polygon: vertices alternate between inner and outer radius.
    const std::size_t N = 10;
    std::vector<Point> pts(N);
    for (std::size_t i = 0; i < N; ++i) {
        double angle = 2.0 * M_PI * static_cast<double>(i) /
                       static_cast<double>(N);
        double r = (i % 2 == 0) ? 2.0 : 0.8;
        pts[i] = {r * std::cos(angle), r * std::sin(angle), i};
    }

    auto nodes = build_vertex_list(pts);
    std::vector<Trapezoid> traps;

    // Find y-max and y-min vertices for a single horizontal cut.
    std::size_t y_max_idx = 0, y_min_idx = 0;
    for (std::size_t i = 1; i < N; ++i) {
        if (pts[i].y > pts[y_max_idx].y) y_max_idx = i;
        if (pts[i].y < pts[y_min_idx].y) y_min_idx = i;
    }

    if (nodes[y_max_idx].next != y_min_idx &&
        nodes[y_min_idx].next != y_max_idx) {
        // Not adjacent → Class B
        traps.push_back({y_max_idx, y_min_idx, FM_NONE, FM_NONE});
        nodes[y_max_idx].trapezoid_idx = 0;
    }

    auto out = fm_triangulate(nodes, traps, 0, N - 1);
    assert(out.size() == N - 2);
    check_valid_triangles(out, N);
    check_no_overlapping_triangles(out);
    check_all_vertices_used(out, N);
    std::printf("  PASS\n\n");
}

// ════════════════════════════════════════════════════════════════════
//  main
// ════════════════════════════════════════════════════════════════════

int main() {
    std::printf("=== Fournier-Montuno comprehensive tests ===\n\n");

    test_triangle();
    test_monotone_quad();
    test_monotone_tall_zigzag();
    test_fm_pentagon();
    test_fm_hexagon();
    test_fm_two_trapezoids();
    test_fm_concave_arrow();
    test_fm_class_a_only();
    test_fm_octagon_multi_split();
    test_fm_large_polygon();
    test_monotone_thin_spike();
    test_monotone_5_vertices();
    test_no_trapezoids();
    test_build_vertex_list();
    test_fm_50gon();
    test_fm_star_polygon();

    std::printf("All %d tests passed.\n", 16);
    return 0;
}
