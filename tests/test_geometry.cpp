/// Unit tests for the geometry module:
///   - Point / Side types
///   - Perturbation predicates (perturbed_y_compare, orient2d, extrema)
///   - Polygon construction, padding, chain ranges, vertex classification

#include "geometry/point.h"
#include "geometry/edge.h"
#include "geometry/perturbation.h"
#include "geometry/polygon.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════════
//  Perturbation predicates
// ════════════════════════════════════════════════════════════════════

static void test_perturbed_y_compare_basic() {
    std::printf("  perturbed_y_compare: basic ordering\n");

    Point a{0.0, 1.0, 0};
    Point b{0.0, 2.0, 1};
    assert(perturbed_y_compare(a, b) < 0); // a below b
    assert(perturbed_y_compare(b, a) > 0); // b above a
    assert(perturbed_y_compare(a, a) == 0); // same point
}

static void test_perturbed_y_compare_tiebreak_x() {
    std::printf("  perturbed_y_compare: tie-break on x\n");

    Point a{1.0, 5.0, 0};
    Point b{3.0, 5.0, 1};
    // Same y → break on x.  a.x < b.x → a < b.
    assert(perturbed_y_compare(a, b) < 0);
    assert(perturbed_y_compare(b, a) > 0);
}

static void test_perturbed_y_compare_tiebreak_index() {
    std::printf("  perturbed_y_compare: tie-break on index\n");

    Point a{5.0, 5.0, 2};
    Point b{5.0, 5.0, 7};
    // Same y, same x → break on index.  a.index < b.index → a < b.
    assert(perturbed_y_compare(a, b) < 0);
    assert(perturbed_y_compare(b, a) > 0);
}

static void test_perturbed_y_less() {
    std::printf("  perturbed_y_less\n");

    Point a{0.0, 1.0, 0};
    Point b{0.0, 2.0, 1};
    assert(perturbed_y_less(a, b));
    assert(!perturbed_y_less(b, a));
    assert(!perturbed_y_less(a, a)); // not strictly less
}

static void test_perturbed_y_order_functor() {
    std::printf("  PerturbedYOrder functor\n");

    PerturbedYOrder cmp;
    Point a{0.0, 1.0, 0};
    Point b{0.0, 2.0, 1};
    Point c{0.0, 1.0, 5}; // same y as a, same x, higher index

    assert(cmp(a, b));
    assert(!cmp(b, a));
    assert(cmp(a, c));  // a.index < c.index
    assert(!cmp(a, a)); // strict weak ordering: !(a < a)
}

// ── orient2d ────────────────────────────────────────────────────────

static void test_orient2d_ccw() {
    std::printf("  orient2d: CCW triangle\n");

    Point a{0.0, 0.0, 0};
    Point b{4.0, 0.0, 1};
    Point c{2.0, 3.0, 2};
    // (b-a)×(c-a) = 4·3 - 0·2 = 12 > 0 → CCW
    assert(orient2d(a, b, c) > 0.0);
    assert(is_convex(a, b, c));
    assert(!is_reflex(a, b, c));
}

static void test_orient2d_cw() {
    std::printf("  orient2d: CW triangle\n");

    Point a{0.0, 0.0, 0};
    Point b{2.0, 3.0, 1};
    Point c{4.0, 0.0, 2};
    // (b-a)×(c-a) = 2·0 - 3·4 = -12 < 0 → CW
    assert(orient2d(a, b, c) < 0.0);
    assert(is_reflex(a, b, c));
    assert(!is_convex(a, b, c));
}

static void test_orient2d_collinear() {
    std::printf("  orient2d: collinear\n");

    Point a{0.0, 0.0, 0};
    Point b{2.0, 2.0, 1};
    Point c{4.0, 4.0, 2};
    assert(orient2d(a, b, c) == 0.0);
    assert(!is_convex(a, b, c));
    assert(!is_reflex(a, b, c));
}

// ── y-extrema ───────────────────────────────────────────────────────

static void test_local_y_minimum() {
    std::printf("  is_local_y_minimum\n");

    Point prev{0.0, 3.0, 0};
    Point curr{1.0, 1.0, 1};
    Point next{2.0, 4.0, 2};
    assert(is_local_y_minimum(prev, curr, next));
    assert(!is_local_y_maximum(prev, curr, next));
    assert(is_local_y_extremum(prev, curr, next));
}

static void test_local_y_maximum() {
    std::printf("  is_local_y_maximum\n");

    Point prev{0.0, 1.0, 0};
    Point curr{1.0, 5.0, 1};
    Point next{2.0, 2.0, 2};
    assert(is_local_y_maximum(prev, curr, next));
    assert(!is_local_y_minimum(prev, curr, next));
    assert(is_local_y_extremum(prev, curr, next));
}

static void test_not_extremum() {
    std::printf("  not an extremum (monotone passage)\n");

    Point prev{0.0, 1.0, 0};
    Point curr{1.0, 3.0, 1};
    Point next{2.0, 5.0, 2};
    // Strictly increasing → not an extremum.
    assert(!is_local_y_minimum(prev, curr, next));
    assert(!is_local_y_maximum(prev, curr, next));
    assert(!is_local_y_extremum(prev, curr, next));
}

// ── horizontal ray intersection ─────────────────────────────────────

static void test_horizontal_ray_x_intercept() {
    std::printf("  horizontal_ray_x_intercept\n");

    Point start{0.0, 0.0, 0};
    Point end{4.0, 4.0, 1};
    // Line from (0,0) to (4,4).  At y=2 → t=0.5, x = 0 + 0.5*4 = 2.
    double x = horizontal_ray_x_intercept(start, end, 2.0);
    assert(std::fabs(x - 2.0) < 1e-12);

    // At y=0 → x=0.
    x = horizontal_ray_x_intercept(start, end, 0.0);
    assert(std::fabs(x - 0.0) < 1e-12);

    // At y=4 → x=4.
    x = horizontal_ray_x_intercept(start, end, 4.0);
    assert(std::fabs(x - 4.0) < 1e-12);
}

// ════════════════════════════════════════════════════════════════════
//  Polygon construction and padding
// ════════════════════════════════════════════════════════════════════

static void test_polygon_triangle() {
    std::printf("  Polygon: triangle (3 vertices → pad to 3)\n");

    // 3 vertices → 2 edges.  2 = 2^1 → need 2^1 + 1 = 3 vertices.
    // Already exactly right.
    std::vector<Point> pts = {{0, 0, 0}, {4, 0, 1}, {2, 3, 2}};
    Polygon poly(pts);

    assert(poly.original_size() == 3);
    assert(poly.num_vertices() == 3);   // 2^1 + 1 = 3
    assert(poly.num_edges() == 2);      // 2^1 = 2
    assert(poly.num_grades() == 1);     // p = 1
}

static void test_polygon_quad_padding() {
    std::printf("  Polygon: quadrilateral (4 vertices → pad to 5)\n");

    // 4 vertices → 3 edges.  bit_ceil(3) = 4 = 2^2 → need 2^2 + 1 = 5.
    std::vector<Point> pts = {{0, 0, 0}, {4, 0, 1}, {4, 4, 2}, {0, 4, 3}};
    Polygon poly(pts);

    assert(poly.original_size() == 4);
    assert(poly.num_vertices() == 5);   // 2^2 + 1 = 5
    assert(poly.num_edges() == 4);      // 2^2 = 4
    assert(poly.num_grades() == 2);     // p = 2
}

static void test_polygon_pentagon_padding() {
    std::printf("  Polygon: pentagon (5 vertices → pad to 5)\n");

    // 5 vertices → 4 edges.  4 = 2^2 → already a power of 2.
    // Need 2^2 + 1 = 5.  Already correct.
    std::vector<Point> pts = {
        {4, 0, 0}, {5, 2, 1}, {2.5, 4, 2}, {0, 2, 3}, {1, 0, 4}
    };
    Polygon poly(pts);

    assert(poly.original_size() == 5);
    assert(poly.num_vertices() == 5);   // 2^2 + 1 = 5
    assert(poly.num_edges() == 4);      // 2^2 = 4
    assert(poly.num_grades() == 2);     // p = 2
}

static void test_polygon_hexagon_padding() {
    std::printf("  Polygon: hexagon (6 vertices → pad to 9)\n");

    // 6 vertices → 5 edges.  bit_ceil(5) = 8 = 2^3 → need 2^3 + 1 = 9.
    std::vector<Point> pts = {
        {2, 0, 0}, {4, 1, 1}, {4, 3, 2}, {2, 4, 3}, {0, 3, 4}, {0, 1, 5}
    };
    Polygon poly(pts);

    assert(poly.original_size() == 6);
    assert(poly.num_vertices() == 9);   // 2^3 + 1 = 9
    assert(poly.num_edges() == 8);      // 2^3 = 8
    assert(poly.num_grades() == 3);     // p = 3
}

static void test_polygon_exact_power_of_two() {
    std::printf("  Polygon: 9 vertices (no padding needed)\n");

    // 9 vertices → 8 edges.  8 = 2^3.  Need 2^3 + 1 = 9.  Exact!
    std::vector<Point> pts;
    for (std::size_t i = 0; i < 9; ++i) {
        double angle = 2.0 * M_PI * static_cast<double>(i) / 9.0;
        pts.push_back({std::cos(angle), std::sin(angle), i});
    }
    Polygon poly(pts);

    assert(poly.original_size() == 9);
    assert(poly.num_vertices() == 9);
    assert(poly.num_edges() == 8);
    assert(poly.num_grades() == 3);
}

static void test_polygon_padding_collinear_vertices() {
    std::printf("  Polygon: padded vertices are collinear on last edge\n");

    // 4 vertices → padded to 5.  The 5th vertex should lie on the
    // segment from vertex 2 to vertex 3 (the last edge before padding).
    std::vector<Point> pts = {
        {0, 0, 0}, {10, 0, 1}, {10, 10, 2}, {0, 10, 3}
    };
    Polygon poly(pts);

    assert(poly.num_vertices() == 5);
    // The inserted vertex (index 3 after re-indexing) should be collinear
    // between original vertex 2 and original vertex 3.
    // Original vertex 2 = (10,10), original vertex 3 = (0,10).
    // Midpoint at t=0.5 → (5, 10).
    const Point& inserted = poly.vertex(3);
    assert(std::fabs(inserted.y - 10.0) < 1e-12);
    assert(inserted.x > 0.0 && inserted.x < 10.0);
}

// ── Chain ranges ────────────────────────────────────────────────────

static void test_chain_ranges_grade0() {
    std::printf("  chain_range: grade 0 (individual edges)\n");

    // 9 vertices → 8 edges, p = 3.
    std::vector<Point> pts;
    for (std::size_t i = 0; i < 9; ++i) {
        pts.push_back({static_cast<double>(i), 0.0, i});
    }
    Polygon poly(pts);

    // Grade 0: 2^(3-0) = 8 chains, each with 2^0 = 1 edge.
    assert(poly.num_chains(0) == 8);
    for (std::size_t j = 0; j < 8; ++j) {
        auto cr = poly.chain_range(0, j);
        assert(cr.start_vertex == j);
        assert(cr.end_vertex == j + 1);
        assert(cr.num_edges() == 1);
        assert(cr.num_vertices() == 2);
    }
}

static void test_chain_ranges_top_grade() {
    std::printf("  chain_range: top grade (entire curve)\n");

    std::vector<Point> pts;
    for (std::size_t i = 0; i < 9; ++i) {
        pts.push_back({static_cast<double>(i), 0.0, i});
    }
    Polygon poly(pts);

    // Grade p = 3: 2^(3-3) = 1 chain, covering all 8 edges.
    assert(poly.num_chains(3) == 1);
    auto cr = poly.chain_range(3, 0);
    assert(cr.start_vertex == 0);
    assert(cr.end_vertex == 8);
    assert(cr.num_edges() == 8);
    assert(cr.num_vertices() == 9);
}

static void test_chain_ranges_middle_grade() {
    std::printf("  chain_range: middle grade\n");

    std::vector<Point> pts;
    for (std::size_t i = 0; i < 9; ++i) {
        pts.push_back({static_cast<double>(i), 0.0, i});
    }
    Polygon poly(pts);

    // Grade 1: 2^(3-1) = 4 chains, each with 2^1 = 2 edges.
    assert(poly.num_chains(1) == 4);
    auto cr0 = poly.chain_range(1, 0);
    assert(cr0.start_vertex == 0 && cr0.end_vertex == 2);
    auto cr3 = poly.chain_range(1, 3);
    assert(cr3.start_vertex == 6 && cr3.end_vertex == 8);
}

// ── Vertex classification ───────────────────────────────────────────

static void test_is_endpoint() {
    std::printf("  is_endpoint\n");

    std::vector<Point> pts = {
        {0, 0, 0}, {4, 0, 1}, {2, 3, 2}
    };
    Polygon poly(pts);

    assert(poly.is_endpoint(0));
    assert(poly.is_endpoint(poly.num_vertices() - 1));
    assert(!poly.is_endpoint(1));
}

static void test_is_y_extremum() {
    std::printf("  is_y_extremum\n");

    // Diamond shape:  top(y=4), right(y=2), bottom(y=0), left(y=2)
    // Padded to 5 vertices (vertex 3 inserted on edge bottom→left).
    // After padding: 0:(0,4) 1:(2,2) 2:(0,0) 3:(interpolated) 4:(-2,2)
    // Vertex 1 is not an extremum (monotone passage 4→2→0).
    // Vertex 2 is a local min if neighbors are above it.
    std::vector<Point> pts = {
        {0, 4, 0}, {2, 2, 1}, {0, 0, 2}, {-2, 2, 3}
    };
    Polygon poly(pts);

    // Endpoints are never extrema.
    assert(!poly.is_y_extremum(0));
    assert(!poly.is_y_extremum(poly.num_vertices() - 1));

    // Vertex 1: prev=(0,4), curr=(2,2), next=(0,0) → decreasing y → not extremum.
    assert(!poly.is_y_extremum(1));

    // Vertex 2: prev=(2,2), curr=(0,0), next=(interpolated between (0,0) and (-2,2))
    // The next vertex after padding is between vertex 2 and vertex 3.
    // Its y > 0, so vertex 2 is a local minimum.
    assert(poly.is_y_extremum(2));
}

// ── Edge table ──────────────────────────────────────────────────────

static void test_edge_table() {
    std::printf("  Edge table\n");

    std::vector<Point> pts = {{0, 0, 0}, {4, 0, 1}, {2, 3, 2}};
    Polygon poly(pts);

    assert(poly.num_edges() == 2);
    for (std::size_t i = 0; i < poly.num_edges(); ++i) {
        const auto& e = poly.edge(i);
        assert(e.index == i);
        assert(e.start_idx == i);
        assert(e.end_idx == i + 1);
    }
}

// ════════════════════════════════════════════════════════════════════
//  main
// ════════════════════════════════════════════════════════════════════

int main() {
    std::printf("=== Geometry tests ===\n\n");

    std::printf("Perturbation predicates:\n");
    test_perturbed_y_compare_basic();
    test_perturbed_y_compare_tiebreak_x();
    test_perturbed_y_compare_tiebreak_index();
    test_perturbed_y_less();
    test_perturbed_y_order_functor();
    test_orient2d_ccw();
    test_orient2d_cw();
    test_orient2d_collinear();
    test_local_y_minimum();
    test_local_y_maximum();
    test_not_extremum();
    test_horizontal_ray_x_intercept();

    std::printf("\nPolygon construction:\n");
    test_polygon_triangle();
    test_polygon_quad_padding();
    test_polygon_pentagon_padding();
    test_polygon_hexagon_padding();
    test_polygon_exact_power_of_two();
    test_polygon_padding_collinear_vertices();

    std::printf("\nChain ranges:\n");
    test_chain_ranges_grade0();
    test_chain_ranges_top_grade();
    test_chain_ranges_middle_grade();

    std::printf("\nVertex classification:\n");
    test_is_endpoint();
    test_is_y_extremum();

    std::printf("\nEdge table:\n");
    test_edge_table();

    std::printf("\nAll geometry tests passed.\n");
    return 0;
}
