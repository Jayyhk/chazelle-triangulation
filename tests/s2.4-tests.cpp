#include "algorithm/polygon/polygon.h"
#include "algorithm/submap/submap.h"
#include "support/assertions.h"

#include <cassert>
#include <cstdio>

using chazelle::test::require_assertion_abort;

using namespace chazelle;

static Polygon test_polygon() {
    return Polygon({{0, 0, 0}, {1, 3, 1}, {2, 1, 2}, {3, 4, 3}, {4, 2, 4}});
}

static Submap make_chordless(const Polygon& poly, std::size_t end_vertex) {
    Submap s;
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = end_vertex;
    Arc a;
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2 * poly.count_nonnull_edges(0, end_vertex - 1);
    std::size_t ai = s.add_arc(a);
    assert(s.start_arc == ai && s.end_arc == ai &&
           "[C91 §2.4(iii) tex 138]: the closed arc is both endpoint arcs");
    return s;
}

static void test_double_identify_simple() {
    auto poly = test_polygon();
    Submap s = make_chordless(poly, 1);

    auto result = s.double_identify(0, {0.0, 0}, poly);
    assert(result.count == 1);

    std::printf("  [PASS] double_identify_simple\n");
}

static void test_double_identify_multi_edge() {
    auto poly = test_polygon();
    Submap s = make_chordless(poly, 4);

    auto r = s.double_identify(2, {1.5, 0}, poly);
    assert(r.count == 1);

    r = s.double_identify(0, {0.0, 0}, poly);
    assert(r.count == 1);

    require_assertion_abort([&] { (void)s.double_identify(5, {0.0, 0}, poly); });

    std::printf("  [PASS] double_identify_multi_edge\n");
}

static void test_double_identify_same_edge() {
    auto poly = test_polygon();

    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 2;

    Arc a;
    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = LEFT;
    a.last_side = LEFT;
    a.region_node = r1;
    a.edge_count = 0;
    std::size_t N = s.add_arc(a);

    a = {};
    a.first_edge = 1;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = LEFT;
    a.region_node = r0;
    a.edge_count = 2 * poly.count_nonnull_edges(0, 1);
    std::size_t W = s.add_arc(a);

    Chord c;
    c = {};
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_edge = 1;
    c.right_edge = 1;
    c.left_side = LEFT;
    c.right_side = LEFT;
    c.y = poly.vertex(1).y;
    c.y_tag = 1;
    c.is_null_length = true;
    c.left_adj = {{W}, 1};
    c.right_adj = {{N}, 1};
    s.add_chord(c);

    assert(s.start_arc == W && s.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

    auto r = s.double_identify(1, {poly.vertex(1).y, 1}, poly);
    assert(r.count == 2);

    std::printf("  [PASS] double_identify_same_edge\n");
}

static void test_endpoint_pointers() {
    auto poly = test_polygon();
    Submap s = make_chordless(poly, 2);

    assert(s.start_arc == 0);
    assert(s.end_arc == 0);
    assert(s.start_vertex == 0);
    assert(s.end_vertex == 2);

    s.check_invariants();

    std::printf("  [PASS] endpoint_pointers\n");
}

static void test_arc_sequence_ordering() {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a;
    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = LEFT;
    a.last_side = LEFT;
    a.region_node = r1;
    a.edge_count = 1;
    std::size_t a1 = s.add_arc(a);
    a = {};
    a.first_edge = 2;
    a.last_edge = 2;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = r2;
    a.edge_count = 2;
    std::size_t aE = s.add_arc(a);
    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r1;
    a.edge_count = 1;
    std::size_t a4 = s.add_arc(a);
    a = {};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = r0;
    a.edge_count = 1;
    std::size_t aS = s.add_arc(a);

    Chord c;
    c = {};
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_adj = {{aS}, 1};
    c.right_adj = {{a4}, 1};
    c.left_edge = 0;
    c.right_edge = 0;
    c.y = 3.0;
    c.y_tag = 1;
    s.add_chord(c);
    c = {};
    c.region[0] = r1;
    c.region[1] = r2;
    c.left_adj = {{a1}, 1};
    c.right_adj = {{aE}, 1};
    c.left_edge = 1;
    c.right_edge = 1;
    c.y = 1.0;
    c.y_tag = 2;
    s.add_chord(c);

    assert(s.start_arc == aS && s.end_arc == aE);
    s.check_invariants();

    std::printf("  [PASS] arc_sequence_ordering\n");
}

static void test_double_identify_extremum() {
    auto poly = test_polygon();

    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();

    s.start_vertex = 0;
    s.end_vertex = 2;

    Arc a;

    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = LEFT;
    a.last_side = LEFT;
    a.region_node = r1;
    a.edge_count = 0;
    std::size_t N = s.add_arc(a);

    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = r0;
    a.edge_count = 2 * poly.count_nonnull_edges(1, 1);
    std::size_t arc2 = s.add_arc(a);

    a = {};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r2;
    a.edge_count = 0;
    std::size_t Z = s.add_arc(a);

    a = {};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = r0;
    a.edge_count = 2 * poly.count_nonnull_edges(0, 0);
    std::size_t arc1 = s.add_arc(a);

    Chord c;
    c = {};
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_edge = 1;
    c.right_edge = 1;
    c.left_side = LEFT;
    c.right_side = LEFT;
    c.y = poly.vertex(1).y;
    c.y_tag = 1;
    c.is_null_length = true;
    c.left_adj = {{arc1}, 1};
    c.right_adj = {{N}, 1};
    s.add_chord(c);

    c = {};
    c.region[0] = r0;
    c.region[1] = r2;
    c.left_edge = 1;
    c.left_side = RIGHT;
    c.left_adj = {{arc2}, 1};
    c.right_edge = 0;
    c.right_side = RIGHT;
    c.right_adj = {{Z}, 1};
    c.y = poly.vertex(1).y;
    c.y_tag = 1;
    s.add_chord(c);

    assert(s.start_arc == arc1 && s.end_arc == arc2);
    s.check_invariants(poly);

    auto rl = s.double_identify(1, {poly.vertex(1).y, 1}, poly);
    assert(rl.count == 2);
    auto rr = s.double_identify(0, {poly.vertex(1).y, 1}, poly);
    assert(rr.count == 2);
    assert(rl.count + rr.count <= Submap::DoubleIdentifyResult::MAX);

    std::printf("  [PASS] double_identify_extremum\n");
}

static void test_double_identify_miss() {
    auto poly = test_polygon();
    Submap s = make_chordless(poly, 1);

    require_assertion_abort([&] { (void)s.double_identify(5, {0.0, 0}, poly); });

    std::printf("  [PASS] double_identify_miss\n");
}

static void test_check_invariants_polygon_positive() {
    auto poly = test_polygon();
    Submap s = make_chordless(poly, 4);
    s.check_invariants(poly);

    std::printf("  [PASS] check_invariants_polygon_positive\n");
}

static Polygon run_polygon() {
    return Polygon({{0, 0, 0}, {1, 10, 1}, {2, 1, 2}, {3, 8, 3}, {4, 4, 4}, {5, 2, 5}});
}

static Submap build_run_submap(const Polygon& poly, bool swap_p2_p3) {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();
    std::size_t r3 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 5;

    Arc a;
    auto mk = [&](std::size_t fe, Side fs, std::size_t le, Side ls, std::size_t region,
                  std::size_t ec) {
        a = {};
        a.first_edge = fe;
        a.first_side = fs;
        a.last_edge = le;
        a.last_side = ls;
        a.region_node = region;
        a.edge_count = ec;
        return s.add_arc(a);
    };

    std::size_t P1 = mk(1, LEFT, 1, LEFT, r1, 1);
    std::size_t P2 = NONE, P3 = NONE;
    if (swap_p2_p3) {
        P3 = mk(1, LEFT, 2, LEFT, r3, 2);
        P2 = mk(1, LEFT, 1, LEFT, r2, 1);
    } else {
        P2 = mk(1, LEFT, 1, LEFT, r2, 1);
        P3 = mk(1, LEFT, 2, LEFT, r3, 2);
    }
    std::size_t P2b = mk(2, LEFT, 2, LEFT, r2, 1);
    std::size_t P1b = mk(2, LEFT, 2, LEFT, r1, 1);

    std::size_t W = mk(3, LEFT, 1, LEFT, r0,
                       poly.count_nonnull_edges(3, 4) + poly.count_nonnull_edges(0, 4) +
                           poly.count_nonnull_edges(0, 1));

    Chord c;

    c = {};
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_edge = 1;
    c.left_side = LEFT;
    c.left_adj = {{W, P1}, 2};
    c.right_edge = 2;
    c.right_side = LEFT;
    c.right_adj = {{P1b}, 1};
    c.y = 8.0;
    c.y_tag = 3;
    s.add_chord(c);

    c = {};
    c.region[0] = r1;
    c.region[1] = r2;
    c.left_edge = 1;
    c.left_side = LEFT;
    c.left_adj = {{P1, P2}, 2};
    c.right_edge = 2;
    c.right_side = LEFT;
    c.right_adj = {{P2b, P1b}, 2};
    c.y = 4.0;
    c.y_tag = 4;
    s.add_chord(c);

    c = {};
    c.region[0] = r2;
    c.region[1] = r3;
    c.left_edge = 1;
    c.left_side = LEFT;
    c.left_adj = {{P2, P3}, 2};
    c.right_edge = 2;
    c.right_side = LEFT;
    c.right_adj = {{P3, P2b}, 2};
    c.y = 2.0;
    c.y_tag = 5;
    s.add_chord(c);

    assert(s.start_arc == W && s.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");
    return s;
}

static void test_check_invariants_polygon_size_3_run() {
    auto poly = run_polygon();
    Submap s = build_run_submap(poly, false);
    s.check_invariants(poly);

    std::printf("  [PASS] check_invariants_polygon_size_3_run\n");
}

static void test_check_invariants_polygon_wrapped_arc() {
    Polygon poly({{0, 0, 0}, {1, 1, 1}, {2, 0, 2}});

    Submap s;
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 2;

    Arc a;
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2 * poly.count_nonnull_edges(0, 1);
    std::size_t ai0 = s.add_arc(a);
    assert(s.start_arc == ai0 && s.end_arc == ai0);

    s.check_invariants(poly);

    std::printf("  [PASS] check_invariants_polygon_wrapped_arc\n");
}

static void test_non_monotonic_run_fires() {
    require_assertion_abort([] {
        auto poly = run_polygon();
        Submap s = build_run_submap(poly, true);
        s.check_invariants(poly);
    });
    std::printf("  [PASS] non_monotonic_run_fires\n");
}

int main() {
    std::printf("[C91 §2.4 tests]:\n");
    test_double_identify_simple();
    test_double_identify_multi_edge();
    test_double_identify_same_edge();
    test_endpoint_pointers();
    test_arc_sequence_ordering();
    test_double_identify_extremum();
    test_double_identify_miss();
    test_check_invariants_polygon_positive();
    test_check_invariants_polygon_size_3_run();
    test_check_invariants_polygon_wrapped_arc();
    test_non_monotonic_run_fires();
    std::printf("All §2.4 tests passed.\n");
    return 0;
}
