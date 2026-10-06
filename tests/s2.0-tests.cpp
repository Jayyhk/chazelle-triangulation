#include "algorithm/polygon/perturbation.h"
#include "algorithm/polygon/point.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <vector>

using namespace chazelle;

static void test_distinct_y() {
    assert(symbolic_y_less({1.0, 0}, {2.0, 0}));
    assert(!symbolic_y_less({2.0, 0}, {1.0, 0}));
    assert(!symbolic_y_less({1.0, 0}, {1.0, 0}));

    assert(symbolic_y_compare({1.0, 0}, {2.0, 0}) == -1);
    assert(symbolic_y_compare({2.0, 0}, {1.0, 0}) == +1);
    assert(symbolic_y_compare({5.0, 3}, {5.0, 3}) == 0);

    std::printf("  [PASS] distinct_y\n");
}

static void test_sos_tiebreak() {
    assert(!symbolic_y_less({4.0, 0}, {4.0, 1}));
    assert(symbolic_y_less({4.0, 1}, {4.0, 0}));

    assert(symbolic_y_compare({4.0, 0}, {4.0, 1}) == +1);
    assert(symbolic_y_compare({4.0, 1}, {4.0, 0}) == -1);

    std::vector<SymbolicY> v = {{4.0, 0}, {4.0, 2}, {4.0, 1}};
    std::sort(v.begin(), v.end(), symbolic_y_less);
    assert(v[0].tag == 2);
    assert(v[1].tag == 1);
    assert(v[2].tag == 0);

    std::printf("  [PASS] sos_tiebreak\n");
}

static void test_transitivity() {
    Point a{0.0, 0.0, 0};
    Point b{0.0, 0.0, 1};
    Point c{0.0, 0.0, 2};

    assert(point_y_order(a, b) == +1);
    assert(point_y_order(b, c) == +1);
    assert(point_y_order(a, c) == +1);

    std::vector<Point> pts = {a, b, c};
    std::sort(pts.begin(), pts.end(), point_y_below);
    assert(pts[0].index == 2);
    assert(pts[1].index == 1);
    assert(pts[2].index == 0);

    std::printf("  [PASS] transitivity\n");
}

static void test_equality() {
    assert(symbolic_y_equal({3.0, 5}, {3.0, 5}));
    assert(!symbolic_y_equal({3.0, 5}, {3.0, 6}));
    assert(!symbolic_y_equal({3.0, 5}, {3.1, 5}));

    assert(symbolic_y_compare({3.0, 5}, {3.0, 5}) == 0);
    assert(symbolic_y_compare({3.0, 5}, {3.0, 6}) != 0);

    std::printf("  [PASS] equality\n");
}

static void test_leq() {
    SymbolicY a{1.0, 0};
    SymbolicY b{2.0, 0};

    assert(symbolic_y_leq(a, a));
    assert(symbolic_y_leq(a, b));
    assert(!symbolic_y_leq(b, a));

    SymbolicY c{1.0, 0};
    SymbolicY d{1.0, 1};
    assert(symbolic_y_leq(d, c));
    assert(!symbolic_y_leq(c, d));

    assert(symbolic_y_geq(a, a));
    assert(symbolic_y_geq(b, a));
    assert(!symbolic_y_geq(a, b));
    assert(symbolic_y_geq(c, d));
    assert(!symbolic_y_geq(d, c));

    assert(!symbolic_y_greater(a, a));
    assert(symbolic_y_greater(b, a));
    assert(!symbolic_y_greater(a, b));
    assert(symbolic_y_greater(c, d));
    assert(!symbolic_y_greater(d, c));

    std::printf("  [PASS] leq_geq_greater\n");
}

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
    assert(!point_y_above(a, a));
    assert(point_y_order(a, b) == -1);
    assert(point_y_order(b, a) == +1);
    assert(point_y_order(a, a) == 0);

    std::printf("  [PASS] point_helpers\n");
}

static void test_sos_none() {
    SymbolicY none_tag{0.0, SOS_NONE};
    SymbolicY zero_tag{0.0, 0};

    assert(symbolic_y_less(none_tag, zero_tag));
    assert(!symbolic_y_less(zero_tag, none_tag));

    SymbolicY none2{0.0, SOS_NONE};
    assert(symbolic_y_equal(none_tag, none2));

    std::printf("  [PASS] sos_none\n");
}

static void test_square_vertices() {
    Point v0{0.0, 0.0, 0};
    Point v1{4.0, 0.0, 1};
    Point v2{4.0, 4.0, 2};
    Point v3{0.0, 4.0, 3};

    assert(point_y_order(v0, v1) == +1);
    assert(point_y_below(v1, v0));

    assert(point_y_order(v2, v3) == +1);
    assert(point_y_below(v3, v2));

    assert(point_y_order(v2, v0) == +1);
    assert(point_y_below(v0, v2));

    std::vector<Point> pts = {v0, v1, v2, v3};
    std::sort(pts.begin(), pts.end(), point_y_below);
    assert(pts[0].index == 1);
    assert(pts[1].index == 0);
    assert(pts[2].index == 3);
    assert(pts[3].index == 2);

    std::printf("  [PASS] square_vertices\n");
}

int main() {
    std::printf("[C91 §2.0 tests]:\n");
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
