/// tests/2.1-tests.cpp — Tests for §2.1: double boundary ∂C classification.

#include "polygon/polygon.h"
#include "polygon/perturbation.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════
//  1. is_endpoint
// ════════════════════════════════════════════════════════════════

static void test_is_endpoint() {
    // Triangle: v0, v1, v2.  Endpoints are v0 and v2.
    Polygon tri({
        {0.0, 0.0, 0}, {4.0, 1.0, 1}, {2.0, 3.0, 2}
    });
    assert( tri.is_endpoint(0));
    assert(!tri.is_endpoint(1));
    assert( tri.is_endpoint(2));

    // 5-vertex curve: endpoints are v0 and v4.
    Polygon curve({
        {0,0,0}, {1,2,1}, {3,1,2}, {5,4,3}, {6,0,4}
    });
    assert( curve.is_endpoint(0));
    assert(!curve.is_endpoint(1));
    assert(!curve.is_endpoint(2));
    assert(!curve.is_endpoint(3));
    assert( curve.is_endpoint(4));

    std::printf("  [PASS] is_endpoint\n");
}

// ════════════════════════════════════════════════════════════════
//  2. is_local_y_minimum / maximum / extremum (free functions)
// ════════════════════════════════════════════════════════════════

static void test_local_extremum_fns() {
    // Clear y-minimum: middle vertex below both neighbors.
    Point a{0, 5, 0}, b{2, 1, 1}, c{4, 5, 2};
    assert( is_local_y_minimum(a, b, c));
    assert(!is_local_y_maximum(a, b, c));
    assert( is_local_y_extremum(a, b, c));

    // Clear y-maximum: middle vertex above both neighbors.
    Point d{0, 1, 0}, e{2, 5, 1}, f{4, 1, 2};
    assert(!is_local_y_minimum(d, e, f));
    assert( is_local_y_maximum(d, e, f));
    assert( is_local_y_extremum(d, e, f));

    // Monotone: ascending.  Not an extremum.
    Point g{0, 0, 0}, h{1, 1, 1}, i{2, 2, 2};
    assert(!is_local_y_minimum(g, h, i));
    assert(!is_local_y_maximum(g, h, i));
    assert(!is_local_y_extremum(g, h, i));

    // Monotone: descending.  Not an extremum.
    Point j{0, 2, 0}, k{1, 1, 1}, l{2, 0, 2};
    assert(!is_local_y_minimum(j, k, l));
    assert(!is_local_y_maximum(j, k, l));
    assert(!is_local_y_extremum(j, k, l));

    std::printf("  [PASS] local_extremum_fns\n");
}

// ════════════════════════════════════════════════════════════════
//  3. SoS-induced extrema (same raw y, resolved by tag)
// ════════════════════════════════════════════════════════════════

static void test_sos_extremum() {
    // All three at raw y=0.  Under SoS: idx 0 highest, idx 2 lowest.
    // So v1 (idx=1) is between v0 (above) and v2 (below) → NOT extremum.
    Point a{0, 0, 0}, b{1, 0, 1}, c{2, 0, 2};
    assert(!is_local_y_extremum(a, b, c));

    // v0 (idx=0) has highest SoS y.  v1 (idx=1) is middle.
    // As prev=v2(idx=2, lowest), curr=v0(idx=0, highest), next=v1(idx=1, middle):
    // v0 is above both → y-maximum.
    assert(is_local_y_maximum(c, a, b));

    // prev=v0(idx=0, highest), curr=v2(idx=2, lowest), next=v1(idx=1, middle):
    // v2 is below both → y-minimum.
    assert(is_local_y_minimum(a, c, b));

    // Two at same y, one different.
    // v0={y=5, idx=0}, v1={y=5, idx=1}, v2={y=0, idx=2}.
    // SoS: v0(y=5+ε₀) > v1(y=5+ε₁) > v2(y=0).
    // v1 is between v0 and v2 → NOT extremum.
    Point d{0, 5, 0}, e{1, 5, 1}, f{2, 0, 2};
    assert(!is_local_y_extremum(d, e, f));

    // v0(y=5+ε₀) above v2(y=0) and above v1(y=5+ε₁) → y-maximum.
    assert(is_local_y_maximum(e, d, f));

    std::printf("  [PASS] sos_extremum\n");
}

// ════════════════════════════════════════════════════════════════
//  4. Polygon::is_y_extremum (combines endpoint + extremum)
// ════════════════════════════════════════════════════════════════

static void test_polygon_is_y_extremum() {
    // V-shape: v0(0,0) → v1(2,3) → v2(4,0).
    // v1 is a y-maximum (above both neighbors).
    // v0, v2 are endpoints → never extrema.
    Polygon v_shape({
        {0, 0, 0}, {2, 3, 1}, {4, 0, 2}
    });
    assert(!v_shape.is_y_extremum(0)); // endpoint
    assert( v_shape.is_y_extremum(1)); // y-maximum
    assert(!v_shape.is_y_extremum(2)); // endpoint

    // W-shape: v0(0,4) → v1(1,0) → v2(2,3) → v3(3,0) → v4(4,4).
    // v1: y-minimum (below v0 and v2).
    // v2: y-maximum (above v1 and v3).
    // v3: y-minimum (below v2 and v4).
    Polygon w_shape({
        {0, 4, 0}, {1, 0, 1}, {2, 3, 2}, {3, 0, 3}, {4, 4, 4}
    });
    assert(!w_shape.is_y_extremum(0)); // endpoint
    assert( w_shape.is_y_extremum(1)); // y-minimum
    assert( w_shape.is_y_extremum(2)); // y-maximum
    assert( w_shape.is_y_extremum(3)); // y-minimum
    assert(!w_shape.is_y_extremum(4)); // endpoint

    // Monotone ascending: no interior extrema.
    Polygon mono({
        {0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}
    });
    assert(!mono.is_y_extremum(0)); // endpoint
    assert(!mono.is_y_extremum(1)); // monotone
    assert(!mono.is_y_extremum(2)); // monotone
    assert(!mono.is_y_extremum(3)); // endpoint

    std::printf("  [PASS] polygon_is_y_extremum\n");
}

// ════════════════════════════════════════════════════════════════
//  5. Fig 2.2 cases: classification matches the three cases
// ════════════════════════════════════════════════════════════════

static void test_fig22_cases() {
    // Curve with all three cases:
    // v0(endpoint) → v1(non-extremum) → v2(y-max) → v3(non-extremum) → v4(endpoint)
    Polygon curve({
        {0, 0, 0}, {1, 2, 1}, {2, 5, 2}, {3, 3, 3}, {4, 1, 4}
    });

    // Case 3 (Fig 2.2.3): endpoints → 2 ∂C vertices (duplicates)
    assert(curve.is_endpoint(0) && !curve.is_y_extremum(0));
    assert(curve.is_endpoint(4) && !curve.is_y_extremum(4));

    // Case 2 (Fig 2.2.2): y-extremum, non-endpoint → 4 ∂C vertices + NLC
    assert(!curve.is_endpoint(2) && curve.is_y_extremum(2));

    // Case 1 (Fig 2.2.1): non-extremum, non-endpoint → 2 ∂C vertices
    assert(!curve.is_endpoint(1) && !curve.is_y_extremum(1));
    assert(!curve.is_endpoint(3) && !curve.is_y_extremum(3));

    std::printf("  [PASS] fig22_cases\n");
}

// ════════════════════════════════════════════════════════════════
//  6. Side enum values
// ════════════════════════════════════════════════════════════════

static void test_side_enum() {
    assert(LEFT == 0);
    assert(RIGHT == 1);
    assert(LEFT != RIGHT);

    Side s = LEFT;
    assert(s == LEFT);
    s = RIGHT;
    assert(s == RIGHT);

    std::printf("  [PASS] side_enum\n");
}

// ════════════════════════════════════════════════════════════════
//  7. Polygon basic construction
// ════════════════════════════════════════════════════════════════

static void test_polygon_construction() {
    Polygon p({
        {0, 0, 0}, {4, 0, 1}, {4, 4, 2}, {0, 4, 3}
    });
    assert(p.num_vertices() == 4);
    assert(p.num_edges() == 3);

    // Edges are consecutive directed pairs.
    assert(p.edge(0).start_idx == 0 && p.edge(0).end_idx == 1);
    assert(p.edge(1).start_idx == 1 && p.edge(1).end_idx == 2);
    assert(p.edge(2).start_idx == 2 && p.edge(2).end_idx == 3);

    // Vertex data preserved.
    assert(p.vertex(0).x == 0.0 && p.vertex(0).y == 0.0);
    assert(p.vertex(2).x == 4.0 && p.vertex(2).y == 4.0);

    std::printf("  [PASS] polygon_construction\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::printf("§2.1 tests:\n");
    test_is_endpoint();
    test_local_extremum_fns();
    test_sos_extremum();
    test_polygon_is_y_extremum();
    test_fig22_cases();
    test_side_enum();
    test_polygon_construction();
    std::printf("All §2.1 tests passed.\n");
    return 0;
}
