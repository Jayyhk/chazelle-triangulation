/// tests/test_section2_0.cpp — Tests for §2.0: SoS y-perturbation.

#include "polygon/point.h"
#include "polygon/perturbation.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <vector>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════
//  1. Distinct y-coordinates — no SoS needed
// ════════════════════════════════════════════════════════════════

static void test_distinct_y() {
    assert(symbolic_y_less({1.0, 0}, {2.0, 0}));
    assert(!symbolic_y_less({2.0, 0}, {1.0, 0}));
    assert(!symbolic_y_less({1.0, 0}, {1.0, 0}));

    assert(symbolic_y_compare({1.0, 0}, {2.0, 0}) == -1);
    assert(symbolic_y_compare({2.0, 0}, {1.0, 0}) == +1);
    assert(symbolic_y_compare({5.0, 3}, {5.0, 3}) ==  0);

    std::printf("  [PASS] distinct_y\n");
}

// ════════════════════════════════════════════════════════════════
//  2. SoS tie-breaking — same raw y, different tags
// ════════════════════════════════════════════════════════════════

static void test_sos_tiebreak() {
    // Lower index → higher perturbed y (per Edelsbrunner [10]).
    // tag 0 is above tag 1 when raw y is equal.
    assert(!symbolic_y_less({4.0, 0}, {4.0, 1})); // tag 0 is NOT below tag 1
    assert( symbolic_y_less({4.0, 1}, {4.0, 0})); // tag 1 IS below tag 0

    assert(symbolic_y_compare({4.0, 0}, {4.0, 1}) == +1); // tag 0 above tag 1
    assert(symbolic_y_compare({4.0, 1}, {4.0, 0}) == -1); // tag 1 below tag 0

    // Three points at same raw y: ascending order is tag 2, 1, 0.
    std::vector<SymbolicY> v = {{4.0, 0}, {4.0, 2}, {4.0, 1}};
    std::sort(v.begin(), v.end(), symbolic_y_less);
    assert(v[0].tag == 2); // lowest
    assert(v[1].tag == 1);
    assert(v[2].tag == 0); // highest

    std::printf("  [PASS] sos_tiebreak\n");
}

// ════════════════════════════════════════════════════════════════
//  3. SoS transitivity
// ════════════════════════════════════════════════════════════════

static void test_transitivity() {
    Point a{0.0, 0.0, 0};
    Point b{0.0, 0.0, 1};
    Point c{0.0, 0.0, 2};

    // a above b, b above c, a above c — no cycles.
    assert(point_y_order(a, b) == +1);
    assert(point_y_order(b, c) == +1);
    assert(point_y_order(a, c) == +1);

    // Sort ascending: [c, b, a].
    std::vector<Point> pts = {a, b, c};
    std::sort(pts.begin(), pts.end(), point_y_below);
    assert(pts[0].index == 2); // c lowest
    assert(pts[1].index == 1); // b middle
    assert(pts[2].index == 0); // a highest

    std::printf("  [PASS] transitivity\n");
}

// ════════════════════════════════════════════════════════════════
//  4. Equality and identity
// ════════════════════════════════════════════════════════════════

static void test_equality() {
    assert( symbolic_y_equal({3.0, 5}, {3.0, 5}));
    assert(!symbolic_y_equal({3.0, 5}, {3.0, 6})); // different tag
    assert(!symbolic_y_equal({3.0, 5}, {3.1, 5})); // different y

    // compare returns 0 only for identical (y, tag).
    assert(symbolic_y_compare({3.0, 5}, {3.0, 5}) == 0);
    assert(symbolic_y_compare({3.0, 5}, {3.0, 6}) != 0);

    std::printf("  [PASS] equality\n");
}

// ════════════════════════════════════════════════════════════════
//  5. leq semantics
// ════════════════════════════════════════════════════════════════

static void test_leq() {
    SymbolicY a{1.0, 0};
    SymbolicY b{2.0, 0};

    assert( symbolic_y_leq(a, a)); // reflexive
    assert( symbolic_y_leq(a, b)); // a < b → a ≤ b
    assert(!symbolic_y_leq(b, a)); // b > a → b ≤ a is false

    // Same raw y, different tags.
    SymbolicY c{1.0, 0}; // higher perturbed y
    SymbolicY d{1.0, 1}; // lower perturbed y
    assert( symbolic_y_leq(d, c)); // d below c → d ≤ c
    assert(!symbolic_y_leq(c, d)); // c above d → c ≤ d is false

    // geq: mirror of leq.
    assert( symbolic_y_geq(a, a)); // reflexive
    assert( symbolic_y_geq(b, a)); // b > a → b ≥ a
    assert(!symbolic_y_geq(a, b)); // a < b → a ≥ b is false
    assert( symbolic_y_geq(c, d)); // c above d → c ≥ d
    assert(!symbolic_y_geq(d, c)); // d below c → d ≥ c is false

    // greater: strict.
    assert(!symbolic_y_greater(a, a)); // irreflexive
    assert( symbolic_y_greater(b, a)); // b > a
    assert(!symbolic_y_greater(a, b)); // a < b
    assert( symbolic_y_greater(c, d)); // c above d
    assert(!symbolic_y_greater(d, c)); // d below c

    std::printf("  [PASS] leq_geq_greater\n");
}

// ════════════════════════════════════════════════════════════════
//  6. Point-level helpers
// ════════════════════════════════════════════════════════════════

static void test_point_helpers() {
    Point p{3.0, 5.0, 7};
    SymbolicY sy = symbolic_y_of(p);
    assert(sy.y == 5.0);
    assert(sy.tag == 7);

    Point a{0.0, 1.0, 0};
    Point b{0.0, 2.0, 0};
    assert(point_y_below(a, b));
    assert(!point_y_below(b, a));
    assert(point_y_above(b, a));
    assert(!point_y_above(a, b));
    assert(!point_y_above(a, a)); // irreflexive
    assert(point_y_order(a, b) == -1);
    assert(point_y_order(b, a) == +1);
    assert(point_y_order(a, a) ==  0);

    std::printf("  [PASS] point_helpers\n");
}

// ════════════════════════════════════════════════════════════════
//  7. SOS_NONE sentinel
// ════════════════════════════════════════════════════════════════

static void test_sos_none() {
    // SOS_NONE is max size_t → largest tag → lowest perturbed y.
    SymbolicY none_tag{0.0, SOS_NONE};
    SymbolicY zero_tag{0.0, 0};

    assert(symbolic_y_less(none_tag, zero_tag)); // none is below everything
    assert(!symbolic_y_less(zero_tag, none_tag));

    // Two SOS_NONE at same y are equal.
    SymbolicY none2{0.0, SOS_NONE};
    assert(symbolic_y_equal(none_tag, none2));

    std::printf("  [PASS] sos_none\n");
}

// ════════════════════════════════════════════════════════════════
//  8. Square polygon vertices (regression)
// ════════════════════════════════════════════════════════════════

static void test_square_vertices() {
    Point v0{0.0, 0.0, 0};
    Point v1{4.0, 0.0, 1};
    Point v2{4.0, 4.0, 2};
    Point v3{0.0, 4.0, 3};

    // Same raw y=0: v0 (idx 0) above v1 (idx 1).
    assert(point_y_order(v0, v1) == +1);
    assert(point_y_below(v1, v0));

    // Same raw y=4: v2 (idx 2) above v3 (idx 3).
    assert(point_y_order(v2, v3) == +1);
    assert(point_y_below(v3, v2));

    // Different raw y: v2 (y=4) above v0 (y=0).
    assert(point_y_order(v2, v0) == +1);
    assert(point_y_below(v0, v2));

    // Full ascending sort: v1, v0, v3, v2.
    std::vector<Point> pts = {v0, v1, v2, v3};
    std::sort(pts.begin(), pts.end(), point_y_below);
    assert(pts[0].index == 1); // v1 lowest
    assert(pts[1].index == 0); // v0
    assert(pts[2].index == 3); // v3
    assert(pts[3].index == 2); // v2 highest

    std::printf("  [PASS] square_vertices\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::printf("§2.0 tests:\n");
    test_distinct_y();
    test_sos_tiebreak();
    test_transitivity();
    test_equality();
    test_leq();
    test_point_helpers();
    test_sos_none();
    test_square_vertices();
    std::printf("All §2.0 tests passed.\n");
    return 0;
}
