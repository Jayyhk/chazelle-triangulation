#include "polygon/polygon.h"
#include "submap/submap.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

static Polygon test_polygon() {
    return Polygon({{0, 0, 0}, {1, 1, 1}, {2, 0.5, 2}, {3, 1.5, 3}});
}

static Submap build_conformal_submap() {
    Submap s;
    s.add_node();
    s.add_node();
    s.add_node();

    s.start_vertex = 0;
    s.end_vertex = 3;

    Arc a;

    a = {};
    a.first_edge = 0;
    a.last_edge = 1;
    a.first_side = LEFT;
    a.last_side = LEFT;
    a.region_node = 1;
    a.edge_count = 2;
    s.add_arc(a);

    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 2;
    a.edge_count = 4;
    s.add_arc(a);

    a = {};
    a.first_edge = 1;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = 1;
    a.edge_count = 2;
    s.add_arc(a);

    a = {};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = 0;
    a.edge_count = 2;
    s.add_arc(a);

    Chord c;

    c = {};
    c.region[0] = 0;
    c.region[1] = 1;
    c.left_edge = 0;
    c.right_edge = 0;
    c.y = 0.5;
    c.y_tag = 2;
    c.left_adj = {{3, 0}, 2};
    c.right_adj = {{2, 3}, 2};
    s.add_chord(c);

    c = {};
    c.region[0] = 1;
    c.region[1] = 2;
    c.left_edge = 1;
    c.right_edge = 1;
    c.y = 1.5;
    c.y_tag = 3;
    c.left_adj = {{0, 1}, 2};
    c.right_adj = {{1, 2}, 2};
    s.add_chord(c);

    assert(s.start_arc == 3 && s.end_arc == 1);
    s.check_invariants(test_polygon());

    return s;
}

static void test_is_conformal() {
    Submap s = build_conformal_submap();
    assert(s.is_conformal());

    s.add_node();
    s.add_node();
    s.add_node();
    Chord c;
    c = {};
    c.region[0] = 1;
    c.region[1] = 3;
    c.y_tag = 0;
    c.left_adj = {{0}, 1};
    c.right_adj = {{0}, 1};
    s.add_chord(c);
    c = {};
    c.region[0] = 1;
    c.region[1] = 4;
    c.y_tag = 0;
    c.left_adj = {{0}, 1};
    c.right_adj = {{0}, 1};
    s.add_chord(c);
    c = {};
    c.region[0] = 1;
    c.region[1] = 5;
    c.y_tag = 0;
    c.left_adj = {{0}, 1};
    c.right_adj = {{0}, 1};
    s.add_chord(c);

    assert(s.node(1).degree() == 5);
    assert(!s.is_conformal());

    std::printf("  [PASS] is_conformal\n");
}

static void test_is_semigranular() {
    Submap s = build_conformal_submap();

    assert(s.is_semigranular(4));
    assert(s.is_semigranular(5));
    assert(!s.is_semigranular(3));
    assert(!s.is_semigranular(1));

    std::printf("  [PASS] is_semigranular\n");
}

static void test_simulated_contraction_weight() {
    Submap s = build_conformal_submap();

    assert(s.simulated_contraction_weight(0, test_polygon()) == 4);

    assert(s.simulated_contraction_weight(1, test_polygon()) == 6);

    std::printf("  [PASS] simulated_contraction_weight\n");
}

static void test_is_granular() {
    Submap s = build_conformal_submap();

    assert(!s.is_granular(3, test_polygon()));

    assert(!s.is_granular(4, test_polygon()));

    assert(!s.is_granular(6, test_polygon()));

    Submap single;
    single.add_node();
    single.start_vertex = 0;
    single.end_vertex = 3;
    Arc a;
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 6;
    single.add_arc(a);
    assert(single.is_granular(6, test_polygon()));
    assert(single.is_granular(8, test_polygon()));

    std::printf("  [PASS] is_granular\n");
}

static void test_is_granular_true() {
    Submap s;
    s.add_node();
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 3;

    Arc a;

    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 1;
    a.edge_count = 4;
    std::size_t E = s.add_arc(a);

    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = 0;
    a.edge_count = 4;
    std::size_t submap = s.add_arc(a);

    Chord c;
    c.region[0] = 0;
    c.region[1] = 1;

    c.left_edge = 1;
    c.right_edge = 1;
    c.y = 1.5;
    c.y_tag = 3;
    c.left_adj = {{submap, E}, 2};
    c.right_adj = {{E, submap}, 2};
    s.add_chord(c);

    assert(s.start_arc == submap && s.end_arc == E &&
           "[C91 §2.4(iii) tex 138]: endpoint arcs auto-registered");
    s.check_invariants(test_polygon());

    assert(s.simulated_contraction_weight(0, test_polygon()) == 6);

    assert(s.is_semigranular(4));
    assert(s.is_granular(4, test_polygon()));

    assert(s.is_semigranular(6));
    assert(!s.is_granular(6, test_polygon()));

    assert(!s.is_granular(3, test_polygon()));

    Submap single;
    single.add_node();
    single.start_vertex = 0;
    single.end_vertex = 2;
    a = {};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 4;
    single.add_arc(a);
    assert(single.is_granular(4, test_polygon()));

    std::printf("  [PASS] is_granular_true\n");
}

static void test_tree_decomposition() {
    Submap s = build_conformal_submap();
    s.build_tree_decomposition();

    const auto& td = s.tree_decomposition();
    assert(!td.empty());

    assert(td.size() == 5);

    assert(td.root() != NONE);
    assert(td.node(td.root()).is_internal());

    std::size_t internals = 0, leaves = 0;
    for (std::size_t i = 0; i < td.size(); ++i) {
        if (td.node(i).is_internal())
            ++internals;
        else
            ++leaves;
    }
    assert(internals == 2);
    assert(leaves == 3);

    std::printf("  [PASS] tree_decomposition\n");
}

static void test_td_single_region() {
    Submap s;
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 3;

    Arc a;
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 3;
    s.add_arc(a);

    s.build_tree_decomposition();
    const auto& td = s.tree_decomposition();

    assert(td.size() == 1);
    assert(td.node(0).is_leaf());
    assert(td.node(0).region_idx == 0);

    std::printf("  [PASS] td_single_region\n");
}

int main() {
    std::printf("[C91 §2.3 tests]:\n");
    test_is_conformal();
    test_is_semigranular();
    test_simulated_contraction_weight();
    test_is_granular();
    test_is_granular_true();
    test_tree_decomposition();
    test_td_single_region();
    std::printf("All §2.3 tests passed.\n");
    return 0;
}
