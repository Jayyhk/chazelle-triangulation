/// Smoke-test for the Fournier-Montuno triangulation module.
///
/// Tests:
///   1. Algorithm 3 (triangulate_monotone) on a simple monotone polygon.
///   2. Algorithm 2+3 (fm_triangulate) on a trapezoidized polygon.

#include "triangulate/trapezoid.h"
#include "triangulate/vertex_node.h"
#include "triangulate/monotone.h"
#include "triangulate/triangulate.h"
#include "geometry/point.h"

#include <cassert>
#include <cstddef>
#include <cstdio>
#include <vector>

using namespace chazelle;

/// Print a triangle list.
static void print_triangles(const char* label,
                            const std::vector<Triangle>& tris) {
    std::printf("  %s (%zu triangles):\n", label, tris.size());
    for (auto& [a, b, c] : tris)
        std::printf("    (%zu, %zu, %zu)\n", a, b, c);
}

/// Test 1: unimonotone triangle (3 vertices — trivial).
static void test_triangle() {
    std::printf("Test 1: triangle (3 vertices)\n");

    // CCW triangle:  (0,0) (4,0) (2,3)
    std::vector<Point> pts = {
        {0.0, 0.0, 0},
        {4.0, 0.0, 1},
        {2.0, 3.0, 2},
    };
    auto nodes = build_vertex_list(pts);
    std::vector<Triangle> out;
    triangulate_monotone(nodes, 0, 2, out);
    print_triangles("monotone", out);
    assert(out.size() == 1);
    std::printf("  PASS\n\n");
}

/// Test 2: unimonotone quadrilateral (4 vertices).
///
///     2
///    / \
///   3   1
///    \ /
///     0
///
/// Base edge: 2→0.  Monotone chain: 0→1→2 (right side).
static void test_monotone_quad() {
    std::printf("Test 2: monotone quadrilateral (4 vertices)\n");

    std::vector<Point> pts = {
        {2.0, 0.0, 0},  // bottom
        {4.0, 1.0, 1},  // right
        {2.0, 3.0, 2},  // top
        {0.0, 1.0, 3},  // left
    };
    auto nodes = build_vertex_list(pts);
    std::vector<Triangle> out;
    triangulate_monotone(nodes, 0, 2, out);
    print_triangles("monotone", out);
    assert(out.size() == 2);
    std::printf("  PASS\n\n");
}

/// Test 3: pentagon with one Class B trapezoid, fed through fm_triangulate.
///
///       2
///      / \
///     3   1
///    / \ / \
///   4---+---0
///
/// We'll manually construct one trapezoid connecting vertex 2 (top)
/// and vertex 0 (bottom) with a Class B diagonal (non-adjacent), and
/// one trapezoid connecting vertex 2 and vertex 4 that is Class A (adjacent).
static void test_fm_pentagon() {
    std::printf("Test 3: fm_triangulate on pentagon (5 vertices)\n");

    std::vector<Point> pts = {
        {4.0, 0.0, 0},
        {5.0, 2.0, 1},
        {2.5, 4.0, 2},
        {0.0, 2.0, 3},
        {1.0, 0.0, 4},
    };
    auto nodes = build_vertex_list(pts);

    // Create trapezoids.
    std::vector<Trapezoid> traps;

    // Trapezoid 0: top = vertex 2, bottom = vertex 0  (Class B — not adjacent).
    traps.push_back({2, 0, FM_NONE, FM_NONE});
    nodes[2].trapezoid_idx = 0;

    auto out = fm_triangulate(nodes, traps, 0, 4);
    print_triangles("fm_triangulate", out);
    // Pentagon → 3 triangles.
    assert(out.size() == 3);
    std::printf("  PASS\n\n");
}

/// Test 4: larger monotone polygon (6 vertices).
///
///     0 (top)
///    / \
///   5   1
///   |   |
///   4   2
///    \ /
///     3 (bottom)
///
/// y-monotone: left chain 0→5→4→3, right chain 0→1→2→3.
/// Base edge: 0→3 does not exist directly → this is NOT unimonotone
/// but it IS monotone.  Algorithm 3 handles it if we set up the
/// linked list correctly as a unimonotone (base = 0→3 via a diagonal).
static void test_monotone_hexagon() {
    std::printf("Test 4: monotone hexagon (6 vertices)\n");

    // CCW hexagon:
    //     0 (top, y=5)
    //    / \
    //   5   1
    //   |   |
    //   4   2
    //    \ /
    //     3 (bottom, y=0)
    //
    // CCW order: 0→5→4→3→2→1→0
    std::vector<Point> pts = {
        {2.0, 5.0, 0},  // top
        {0.0, 4.0, 1},  // left-top
        {0.0, 2.0, 2},  // left-bottom
        {2.0, 0.0, 3},  // bottom
        {4.0, 2.0, 4},  // right-bottom
        {4.0, 4.0, 5},  // right-top
    };
    auto nodes = build_vertex_list(pts);
    std::vector<Trapezoid> traps;

    // Trapezoid: top = 0, bottom = 3 (Class B: vertices 0 and 3 are not adjacent)
    traps.push_back({0, 3, FM_NONE, FM_NONE});
    nodes[0].trapezoid_idx = 0;

    auto out = fm_triangulate(nodes, traps, 0, 5);
    print_triangles("fm_triangulate", out);
    assert(out.size() == 4); // hexagon -> 4 triangles
    std::printf("  PASS\n\n");
}

int main() {
    std::printf("=== Fournier-Montuno smoke tests ===\n\n");

    test_triangle();
    test_monotone_quad();
    test_fm_pentagon();
    test_monotone_hexagon();

    std::printf("All tests passed.\n");
    return 0;
}
