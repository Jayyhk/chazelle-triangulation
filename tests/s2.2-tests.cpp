#include "algorithm/merge/fusion.h"
#include "algorithm/polygon/polygon.h"
#include "algorithm/submap/chord_inventory.h"
#include "algorithm/submap/submap.h"
#include "support/assertions.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

using chazelle::test::require_assertion_abort;

static Polygon test_polygon() {
    return Polygon({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4}});
}

static Submap build_3region_submap() {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();

    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a1;
    a1.first_edge = 1;
    a1.last_edge = 1;
    a1.first_side = LEFT;
    a1.last_side = LEFT;
    a1.region_node = r1;
    a1.edge_count = 1;
    std::size_t ai1 = s.add_arc(a1);

    Arc aE;
    aE.first_edge = 2;
    aE.last_edge = 2;
    aE.first_side = LEFT;
    aE.last_side = RIGHT;
    aE.region_node = r2;
    aE.edge_count = 4;
    std::size_t aiE = s.add_arc(aE);

    Arc a4;
    a4.first_edge = 1;
    a4.last_edge = 1;
    a4.first_side = RIGHT;
    a4.last_side = RIGHT;
    a4.region_node = r1;
    a4.edge_count = 1;
    std::size_t ai4 = s.add_arc(a4);

    Arc aS;
    aS.first_edge = 0;
    aS.last_edge = 0;
    aS.first_side = RIGHT;
    aS.last_side = LEFT;
    aS.region_node = r0;
    aS.edge_count = 2;
    std::size_t aiS = s.add_arc(aS);

    Chord c0;
    c0.region[0] = r0;
    c0.region[1] = r1;
    c0.left_adj = {{aiS}, 1};
    c0.right_adj = {{ai4}, 1};
    c0.left_edge = 0;
    c0.right_edge = 0;
    c0.y = 1.0;
    c0.y_tag = 1;
    s.add_chord(c0);

    Chord c1;
    c1.region[0] = r1;
    c1.region[1] = r2;
    c1.left_adj = {{ai1}, 1};
    c1.right_adj = {{aiE}, 1};
    c1.left_edge = 1;
    c1.right_edge = 1;
    c1.y = 2.0;
    c1.y_tag = 2;
    s.add_chord(c1);

    assert(s.start_arc == aiS && s.end_arc == aiE);

    return s;
}

static void test_count_nonnull_edges() {
    Polygon tri({{0, 0, 0}, {4, 1, 1}, {2, 3, 2}});
    assert(tri.count_nonnull_edges(0, 0) == 1);
    assert(tri.count_nonnull_edges(1, 1) == 1);
    assert(tri.count_nonnull_edges(0, 1) == 2);

    Polygon p({{0, 0, 0}, {3, 3, 1}, {3, 3, 2}, {5, 1, 3}});
    assert(p.count_nonnull_edges(1, 1) == 0);
    assert(p.count_nonnull_edges(0, 2) == 2);

    std::printf("  [PASS] count_nonnull_edges\n");
}

static void test_submap_construction() {
    Submap s = build_3region_submap();

    assert(s.num_nodes() == 3);
    assert(s.num_chords() == 2);

    assert(s.num_arcs() == 4);

    s.assert_tree_property();

    assert(s.node(0).degree() == 1);
    assert(s.node(1).degree() == 2);
    assert(s.node(2).degree() == 1);

    std::printf("  [PASS] submap_construction\n");
}

static void test_check_invariants() {
    Submap s = build_3region_submap();
    s.check_invariants();

    std::printf("  [PASS] check_invariants\n");
}

static void test_region_weight() {
    Submap s = build_3region_submap();

    assert(s.region_weight(0) == 2);

    assert(s.region_weight(1) == 1);

    assert(s.region_weight(2) == 4);

    std::printf("  [PASS] region_weight\n");
}

static void test_remove_chord() {
    Submap s = build_3region_submap();

    auto poly = test_polygon();
    [[maybe_unused]] std::size_t survivor = s.remove_chord(1, poly);

    assert(s.num_live_nodes() == 2);
    assert(s.num_live_chords() == 1);
    s.assert_tree_property();

    std::printf("  [PASS] remove_chord\n");
}

static void test_remove_all_chords() {
    Submap s = build_3region_submap();

    auto poly = test_polygon();
    s.remove_chord(0, poly);
    s.remove_chord(1, poly);

    s.assert_tree_property();

    assert(s.num_live_nodes() == 1);
    assert(s.num_live_chords() == 0);

    s.compact();
    assert(s.num_nodes() == 1);
    assert(s.num_chords() == 0);

    std::printf("  [PASS] remove_all_chords\n");
}

static void test_chordless_region_weight() {
    auto poly = test_polygon();
    Submap s;
    std::size_t r0 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a;
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = r0;
    a.edge_count = 2 * poly.count_nonnull_edges(0, 3);
    std::size_t ai = s.add_arc(a);
    assert(s.start_arc == ai && s.end_arc == ai &&
           "[C91 §2.4(iii) tex 138]: the closed arc is both endpoint arcs");

    assert(s.region_weight(r0) == 8);
    s.check_invariants(poly);

    std::printf("  [PASS] chordless_region_weight\n");
}

static void test_remove_chord_4_adj_arcs() {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a;
    a = {};
    a.first_edge = 1;
    a.last_edge = 2;
    a.first_side = LEFT;
    a.last_side = LEFT;
    a.region_node = r1;
    a.edge_count = 2;
    std::size_t B = s.add_arc(a);
    a = {};
    a.first_edge = 3;
    a.last_edge = 3;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = r2;
    a.edge_count = 2;
    std::size_t E = s.add_arc(a);
    a = {};
    a.first_edge = 2;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r1;
    a.edge_count = 2;
    std::size_t D = s.add_arc(a);
    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = r0;
    a.edge_count = 4;
    std::size_t submap = s.add_arc(a);

    Chord c;
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_adj = {{submap, B}, 2};
    c.right_adj = {{D, submap}, 2};
    c.left_edge = 1;
    c.right_edge = 1;
    c.left_side = LEFT;
    c.right_side = RIGHT;
    c.y = 1.5;
    c.y_tag = 99;
    s.add_chord(c);

    Chord cv;
    cv.region[0] = r1;
    cv.region[1] = r2;
    cv.left_adj = {{B}, 1};
    cv.right_adj = {{E}, 1};
    cv.left_edge = 2;
    cv.right_edge = 2;
    cv.left_side = LEFT;
    cv.right_side = RIGHT;
    cv.y = 3.0;
    cv.y_tag = 3;
    s.add_chord(cv);

    auto poly = test_polygon();
    s.remove_chord(0, poly);

    assert(s.num_live_arcs() == 2 && "two glues chain S,B,D into one arc");
    bool found_merged = false;
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead)
            continue;
        if (s.arc(i).first_side == RIGHT && s.arc(i).last_side == LEFT) {
            assert(s.arc(i).first_edge == 2 && s.arc(i).last_edge == 2 &&
                   s.arc(i).edge_count == 6 && "merged start-wrap arc must span ᾱ = [0,2]");
            found_merged = true;
        }
    }
    assert(found_merged && "the glued chain stays one double-backing "
                           "arc-structure ([C91 §2.4 tex 142])");
    s.assert_tree_property();

    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] remove_chord_4_adj_arcs\n");
}

static void test_remove_chord_merge_at_vertex() {
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
    a.last_side = RIGHT;
    a.region_node = r1;
    a.edge_count = 2;
    std::size_t E = s.add_arc(a);
    a = {};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = r0;
    a.edge_count = 2;
    std::size_t submap = s.add_arc(a);

    Chord c;
    c.region[0] = r0;
    c.region[1] = r1;

    c.left_adj = {{submap}, 1};
    c.right_adj = {{E}, 1};

    c.left_edge = 0;
    c.right_edge = 0;
    c.left_side = LEFT;
    c.right_side = RIGHT;
    c.y = 1.0;
    c.y_tag = 1;
    s.add_chord(c);

    auto poly = test_polygon();

    s.check_invariants(poly);

    s.remove_chord(0, poly);

    assert(s.num_live_arcs() == 1 && "[C91 §2.2 tex 96]: vertex endpoints must glue their arcs");
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead)
            continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 &&
               "closed arc is cut at C's start turnaround "
               "([C91 §2.4(iii) tex 138])");
        assert(s.arc(i).edge_count == 4 && "closed arc spans all of C through the vertex");
    }

    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] remove_chord_merge_at_vertex\n");
}

static void test_check_invariants_offset_subchain() {
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
    a.last_side = RIGHT;
    a.region_node = r1;
    a.edge_count = 2;
    std::size_t E = s.add_arc(a);
    a = {};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = r0;
    a.edge_count = 2;
    std::size_t submap = s.add_arc(a);

    Chord c;
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_adj = {{submap}, 1};
    c.right_adj = {{E}, 1};
    c.left_edge = 0;
    c.right_edge = 0;
    c.left_side = LEFT;
    c.right_side = RIGHT;

    c.y = 2.0;
    c.y_tag = 2;
    s.add_chord(c);

    auto poly = test_polygon().subchain(1, 3);
    s.check_invariants(poly);

    s.remove_chord(0, poly);
    assert(s.num_live_arcs() == 1 && "[C91 §2.2 tex 96]: vertex endpoints must glue their arcs");
    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] check_invariants_offset_subchain\n");
}

static void test_remove_chord_2_adj_arcs() {
    Submap s;
    std::size_t r_main = s.add_node();
    std::size_t r_cap = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a;

    a = {};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = r_main;
    a.edge_count = 8;
    std::size_t M = s.add_arc(a);

    a = {};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = r_cap;
    a.edge_count = 0;
    std::size_t Z = s.add_arc(a);

    Chord c;
    c.region[0] = r_main;
    c.region[1] = r_cap;
    c.left_adj = {{Z}, 1};
    c.right_adj = {{M}, 1};
    c.left_edge = 0;
    c.right_edge = 0;
    c.left_side = LEFT;
    c.right_side = RIGHT;
    c.y = 0.0;
    c.y_tag = 0;
    s.add_chord(c);

    assert(s.start_arc == Z && s.end_arc == M &&
           "[C91 §2.4(iii) tex 138]: endpoint arcs auto-registered");

    auto poly = test_polygon();
    s.remove_chord(0, poly);

    assert(s.num_live_arcs() == 1 && "companion-chord removal closes ∂C into one arc");
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead)
            continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 && s.arc(i).edge_count == 8 &&
               "[C91 §2.4 tex 142]: closed arc covering all of C");
    }
    s.assert_tree_property();

    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] remove_chord_2_adj_arcs\n");
}

static Polygon non_monotone_polygon() {
    return Polygon({{0, 0, 0}, {10, 10, 1}, {20, -5, 2}, {30, 5, 3}, {40, 15, 4}});
}

static void test_remove_chord_3_adj_arcs() {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a;

    a = {};
    a.first_edge = 0;
    a.last_edge = 2;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = r0;
    a.edge_count = 6;
    std::size_t E = s.add_arc(a);

    a = {};
    a.first_edge = 2;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = r1;
    a.edge_count = 4;
    std::size_t submap = s.add_arc(a);

    Chord c;
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_adj = {{submap}, 1};
    c.right_adj = {{E, submap}, 2};
    c.left_edge = 0;
    c.left_side = LEFT;
    c.right_edge = 2;
    c.right_side = RIGHT;
    c.y = 0.0;
    c.y_tag = 0;
    s.add_chord(c);

    assert(s.start_arc == submap && s.end_arc == E);

    auto poly = non_monotone_polygon();
    s.remove_chord(0, poly);

    assert(s.num_live_arcs() == 1 && "last-chord removal closes ∂C");
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead)
            continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 && s.arc(i).edge_count == 8 &&
               "closed arc covering all of C");
    }
    s.assert_tree_property();

    std::printf("  [PASS] remove_chord_3_adj_arcs\n");
}

static void test_remove_chord_shared_arc_left() {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a;

    a = {};
    a.first_edge = 1;
    a.last_edge = 3;
    a.first_side = LEFT;
    a.last_side = LEFT;
    a.region_node = r1;
    a.edge_count = 3;
    std::size_t A = s.add_arc(a);

    a = {};
    a.first_edge = 3;
    a.last_edge = 1;
    a.first_side = LEFT;
    a.last_side = LEFT;
    a.region_node = r0;
    a.edge_count = 7;
    std::size_t W = s.add_arc(a);

    Chord c;
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_edge = 1;
    c.left_side = LEFT;
    c.right_edge = 3;
    c.right_side = LEFT;
    c.y = 1.5;
    c.y_tag = 99;
    c.left_adj = {{W, A}, 2};
    c.right_adj = {{A, W}, 2};
    s.add_chord(c);

    assert(s.start_arc == W && s.end_arc == W);

    auto poly = test_polygon();

    assert(s.simulated_contraction_weight(0, poly) == 8);

    s.remove_chord(0, poly);
    s.assert_tree_property();
    assert(s.num_live_arcs() == 1 && "shared-arc removal of the last chord closes ∂C into one arc");
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead)
            continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 && s.arc(i).edge_count == 8 &&
               "closed arc must span the full glued range");
    }

    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] remove_chord_shared_arc_left\n");
}

static void test_remove_chord_shared_arc_right() {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a;

    a = {};
    a.first_edge = 3;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r1;
    a.edge_count = 3;
    std::size_t A = s.add_arc(a);

    a = {};
    a.first_edge = 1;
    a.last_edge = 3;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r0;
    a.edge_count = 7;
    std::size_t W = s.add_arc(a);

    Chord c;
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_edge = 3;
    c.left_side = RIGHT;
    c.right_edge = 1;
    c.right_side = RIGHT;
    c.y = 1.5;
    c.y_tag = 99;
    c.left_adj = {{W, A}, 2};
    c.right_adj = {{A, W}, 2};
    s.add_chord(c);

    assert(s.start_arc == W && s.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

    auto poly = test_polygon();
    assert(s.simulated_contraction_weight(0, poly) == 8);

    s.remove_chord(0, poly);
    s.assert_tree_property();
    assert(s.num_live_arcs() == 1);
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead)
            continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 && s.arc(i).edge_count == 8 &&
               "closed arc must span the full glued range");
    }

    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] remove_chord_shared_arc_right\n");
}

static void test_remove_chord_fully_wrapped_closes() {
    auto poly = test_polygon();
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a;
    a = {};
    a.first_edge = 1;
    a.last_edge = 2;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = r1;

    a.edge_count = poly.count_nonnull_edges(1, 3) + poly.count_nonnull_edges(2, 3);
    std::size_t A = s.add_arc(a);

    a = {};
    a.first_edge = 2;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = r0;

    a.edge_count = poly.count_nonnull_edges(0, 2) + poly.count_nonnull_edges(0, 1);
    std::size_t B = s.add_arc(a);

    Chord c;
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_edge = 1;
    c.left_side = LEFT;
    c.right_edge = 2;
    c.right_side = RIGHT;
    c.y = 1.5;
    c.y_tag = 99;
    c.left_adj = {{B, A}, 2};
    c.right_adj = {{A, B}, 2};
    s.add_chord(c);

    assert(s.start_arc == B && s.end_arc == A &&
           "[C91 §2.4(iii) tex 138]: wrap arcs auto-registered");

    s.remove_chord(0, poly);

    assert(s.num_live_nodes() == 1 && s.num_live_chords() == 0);
    assert(s.num_live_arcs() == 1 && "[C91 §2.2 tex 94]: last-chord removal closes ∂C");
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead)
            continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 && s.arc(i).edge_count == 8 &&
               "[C91 §2.4 tex 142]: the closed arc covers all of C");
    }
    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] remove_chord_fully_wrapped_closes\n");
}

static void test_null_length_chord() {
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
    a.edge_count = 4;
    std::size_t W = s.add_arc(a);

    Chord c;
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_adj = {{W}, 1};
    c.right_adj = {{N}, 1};
    c.left_edge = 1;
    c.right_edge = 1;
    c.left_side = LEFT;
    c.right_side = LEFT;
    c.is_null_length = true;
    c.y = 1.0;
    c.y_tag = 1;
    s.add_chord(c);

    assert(s.start_arc == W && s.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

    assert(s.chord(0).is_null_length);
    assert(s.node(r1).degree() == 1);

    assert(s.region_weight(r1) == 0);

    s.assert_tree_property();

    SymbolicY sy = s.chord(0).symbolic_y();
    assert(sy.y == 1.0 && sy.tag == 1);

    auto poly = Polygon({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}});
    s.remove_chord(0, poly);
    assert(s.num_live_arcs() == 1 && s.num_live_nodes() == 1);
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead)
            continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 && s.arc(i).edge_count == 4);
    }
    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] null_length_chord\n");
}

static void test_null_length_chord_right_side() {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 2;

    Arc a;
    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r1;
    a.edge_count = 0;
    std::size_t N = s.add_arc(a);

    a = {};
    a.first_edge = 0;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r0;
    a.edge_count = 4;
    std::size_t W = s.add_arc(a);

    Chord c;
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_adj = {{W}, 1};
    c.right_adj = {{N}, 1};
    c.left_edge = 1;
    c.right_edge = 1;
    c.left_side = RIGHT;
    c.right_side = RIGHT;
    c.is_null_length = true;
    c.y = 1.0;
    c.y_tag = 1;
    s.add_chord(c);

    assert(s.start_arc == W && s.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");
    assert(s.chord(0).is_null_length);
    assert(s.node(r1).degree() == 1);
    assert(s.region_weight(r1) == 0);
    s.assert_tree_property();

    std::printf("  [PASS] null_length_chord_right_side\n");
}

static void test_empty_submap_fires() {
    require_assertion_abort([] {
        Submap s;
        s.assert_tree_property();
    });
    std::printf("  [PASS] empty_submap_fires\n");
}

static void test_null_length_chord_mismatched_edges_fires() {
    require_assertion_abort([] {
        Submap s;
        std::size_t r0 = s.add_node();
        std::size_t r1 = s.add_node();

        Arc a;
        a = {};
        a.first_edge = 0;
        a.last_edge = 0;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = r0;
        a.edge_count = 1;
        std::size_t ai0 = s.add_arc(a);

        a = {};
        a.first_edge = 1;
        a.last_edge = 1;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = r1;
        a.edge_count = 0;
        std::size_t ai1 = s.add_arc(a);

        Chord c;
        c.region[0] = r0;
        c.region[1] = r1;
        c.left_edge = 0;
        c.right_edge = 1;
        c.left_side = LEFT;
        c.right_side = LEFT;
        c.is_null_length = true;
        c.y = 1.0;
        c.y_tag = 1;
        c.left_adj = {{ai0}, 1};
        c.right_adj = {{ai1}, 1};
        s.add_chord(c);
    });
    std::printf("  [PASS] null_length_chord_mismatched_edges_fires\n");
}

static void test_null_length_chord_mismatched_sides_fires() {
    require_assertion_abort([] {
        Submap s;
        std::size_t r0 = s.add_node();
        std::size_t r1 = s.add_node();

        Arc a;
        a = {};
        a.first_edge = 1;
        a.last_edge = 1;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = r0;
        a.edge_count = 1;
        std::size_t ai0 = s.add_arc(a);

        a = {};
        a.first_edge = 1;
        a.last_edge = 1;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = r1;
        a.edge_count = 0;
        std::size_t ai1 = s.add_arc(a);

        Chord c;
        c.region[0] = r0;
        c.region[1] = r1;
        c.left_edge = 1;
        c.right_edge = 1;
        c.left_side = LEFT;
        c.right_side = RIGHT;
        c.is_null_length = true;
        c.y = 1.0;
        c.y_tag = 1;
        c.left_adj = {{ai0}, 1};
        c.right_adj = {{ai1}, 1};
        s.add_chord(c);
    });
    std::printf("  [PASS] null_length_chord_mismatched_sides_fires\n");
}

static void test_null_length_chord_missing_y_tag_fires() {
    require_assertion_abort([] {
        Submap s;
        std::size_t r0 = s.add_node();
        std::size_t r1 = s.add_node();

        Arc a;
        a = {};
        a.first_edge = 1;
        a.last_edge = 1;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = r0;
        a.edge_count = 1;
        std::size_t ai0 = s.add_arc(a);

        a = {};
        a.first_edge = 1;
        a.last_edge = 1;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = r1;
        a.edge_count = 0;
        std::size_t ai1 = s.add_arc(a);

        Chord c;
        c.region[0] = r0;
        c.region[1] = r1;
        c.left_edge = 1;
        c.right_edge = 1;
        c.left_side = LEFT;
        c.right_side = LEFT;
        c.is_null_length = true;
        c.y = 1.0;

        c.left_adj = {{ai0}, 1};
        c.right_adj = {{ai1}, 1};
        s.add_chord(c);
    });
    std::printf("  [PASS] null_length_chord_missing_y_tag_fires\n");
}

static void test_null_length_chord_wrong_adj_count_fires() {
    require_assertion_abort([] {
        Submap s;
        std::size_t r0 = s.add_node();
        std::size_t r1 = s.add_node();

        Arc a;
        a = {};
        a.first_edge = 1;
        a.last_edge = 1;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = r0;
        a.edge_count = 1;
        std::size_t ai0 = s.add_arc(a);

        a = {};
        a.first_edge = 1;
        a.last_edge = 1;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = r0;
        a.edge_count = 1;
        std::size_t ai0b = s.add_arc(a);

        a = {};
        a.first_edge = 1;
        a.last_edge = 1;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = r1;
        a.edge_count = 0;
        std::size_t ai1 = s.add_arc(a);

        Chord c;
        c.region[0] = r0;
        c.region[1] = r1;
        c.left_edge = 1;
        c.right_edge = 1;
        c.left_side = LEFT;
        c.right_side = LEFT;
        c.is_null_length = true;
        c.y = 1.0;
        c.y_tag = 1;
        c.left_adj = {{ai0, ai0b}, 2};
        c.right_adj = {{ai1}, 1};
        s.add_chord(c);
    });
    std::printf("  [PASS] null_length_chord_wrong_adj_count_fires\n");
}

static void test_wrong_edge_count_fires() {
    require_assertion_abort([] {
        Polygon poly({{0, 0, 0}, {1, 3, 1}, {2, 1, 2}, {3, 4, 3}, {4, 2, 4}});
        Submap s;
        s.add_node();
        s.start_vertex = 0;
        s.end_vertex = 4;

        Arc a;
        a = {};
        a.first_edge = 0;
        a.last_edge = 0;
        a.first_side = LEFT;
        a.last_side = RIGHT;
        a.region_node = 0;
        a.edge_count = 99;
        s.add_arc(a);

        s.check_invariants(poly);
    });
    std::printf("  [PASS] wrong_edge_count_fires\n");
}

static void test_augmented_nonvertex_level() {
    Polygon curve({{0, 0, 0}, {0, 4, 1}, {4, 5, 2}, {4, -1, 3}});
    PendingChord p;
    p.y = {2.0, 42};
    p.left_edge_c = 0;
    p.left_side = RIGHT;
    p.right_edge_c = 2;
    p.right_side = RIGHT;
    Submap submap;
    build_submap_from_chords(submap, curve, {p});
    submap.check_invariants(curve);
    assert(submap.num_chords() == 1 && submap.num_nodes() == 2);
    assert(submap.chord(0).left_adj.count == 2 && submap.chord(0).right_adj.count == 2);
    auto ids = submap.double_identify(0, p.y, curve);
    assert(ids.count >= 2 && "[C91 §2.4 tex 144]: identify augmented endpoints");
    const std::size_t weight = submap.simulated_contraction_weight(0, curve);
    const std::size_t survivor = submap.remove_chord(0, curve);
    assert(submap.region_weight(survivor) == weight);
    submap.normalize(curve);
    submap.check_invariants(curve);
    assert(submap.num_nodes() == 1 && submap.num_chords() == 0);
    std::printf("  [PASS] augmented_nonvertex_level\n");
}

int main() {
    std::printf("[C91 §2.2 tests]:\n");
    test_count_nonnull_edges();
    test_submap_construction();
    test_check_invariants();
    test_region_weight();
    test_remove_chord();
    test_remove_all_chords();
    test_chordless_region_weight();
    test_remove_chord_4_adj_arcs();
    test_remove_chord_merge_at_vertex();
    test_check_invariants_offset_subchain();
    test_remove_chord_2_adj_arcs();
    test_remove_chord_3_adj_arcs();
    test_remove_chord_shared_arc_left();
    test_remove_chord_shared_arc_right();
    test_remove_chord_fully_wrapped_closes();
    test_null_length_chord();
    test_null_length_chord_right_side();
    test_empty_submap_fires();
    test_null_length_chord_mismatched_edges_fires();
    test_null_length_chord_mismatched_sides_fires();
    test_null_length_chord_missing_y_tag_fires();
    test_null_length_chord_wrong_adj_count_fires();
    test_wrong_edge_count_fires();
    test_augmented_nonvertex_level();
    std::printf("All §2.2 tests passed.\n");
    return 0;
}
