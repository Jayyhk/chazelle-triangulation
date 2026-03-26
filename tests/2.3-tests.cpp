/// tests/2.3-tests.cpp — Tests for §2.3: conformality, granularity, tree decomposition.

#include "submap/submap.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════
//  Helper: build a conformal submap (degree ≤ 4).
//  Linear chain: r0 — c0 — r1 — c1 — r2
// ════════════════════════════════════════════════════════════════

static Submap build_conformal_submap() {
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1
    s.add_node(); // r2

    // Arcs: 2 per region, LEFT then RIGHT.
    Arc a;

    // r0 LEFT
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    s.add_arc(a);
    // r1 LEFT
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 1; a.edge_count = 3;
    s.add_arc(a);
    // r2 LEFT
    a = {}; a.first_edge = 2; a.last_edge = 2; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 2; a.edge_count = 1;
    s.add_arc(a);
    // r2 RIGHT
    a = {}; a.first_edge = 2; a.last_edge = 2; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 2; a.edge_count = 1;
    s.add_arc(a);
    // r1 RIGHT
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 3;
    s.add_arc(a);
    // r0 RIGHT
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    s.add_arc(a);

    Chord c;
    // c0: r0 — r1
    c = {}; c.region[0] = 0; c.region[1] = 1;
    c.adj_arcs = {0, 1, 5, 4}; c.num_adj_arcs = 4;
    s.add_chord(c);
    // c1: r1 — r2
    c = {}; c.region[0] = 1; c.region[1] = 2;
    c.adj_arcs = {1, 2, 4, 3}; c.num_adj_arcs = 4;
    s.add_chord(c);

    return s;
}

// ════════════════════════════════════════════════════════════════
//  1. is_conformal
// ════════════════════════════════════════════════════════════════

static void test_is_conformal() {
    Submap s = build_conformal_submap();
    assert(s.is_conformal());

    // Add a 5th chord to r1 → degree 5 → not conformal.
    // (Need 2 more nodes for the extra chords.)
    s.add_node(); // r3
    s.add_node(); // r4
    s.add_node(); // r5
    Chord c;
    c = {}; c.region[0] = 1; c.region[1] = 3;
    c.adj_arcs = {1, 1, NONE, NONE}; c.num_adj_arcs = 2;
    s.add_chord(c);
    c = {}; c.region[0] = 1; c.region[1] = 4;
    c.adj_arcs = {1, 1, NONE, NONE}; c.num_adj_arcs = 2;
    s.add_chord(c);
    c = {}; c.region[0] = 1; c.region[1] = 5;
    c.adj_arcs = {1, 1, NONE, NONE}; c.num_adj_arcs = 2;
    s.add_chord(c);

    // r1 now has degree 5.
    assert(s.node(1).degree() == 5);
    assert(!s.is_conformal());

    std::printf("  [PASS] is_conformal\n");
}

// ════════════════════════════════════════════════════════════════
//  2. is_semigranular
// ════════════════════════════════════════════════════════════════

static void test_is_semigranular() {
    Submap s = build_conformal_submap();
    // Weights: r0=2, r1=3, r2=1.

    assert(s.is_semigranular(3));  // all ≤ 3
    assert(s.is_semigranular(4));  // all ≤ 4
    assert(!s.is_semigranular(2)); // r1 has weight 3 > 2

    std::printf("  [PASS] is_semigranular\n");
}

// ════════════════════════════════════════════════════════════════
//  3. simulated_contraction_weight
// ════════════════════════════════════════════════════════════════

static void test_simulated_contraction_weight() {
    Submap s = build_conformal_submap();

    // c0 adj_arcs = {0,1,5,4}.
    // Pair [0]+[1]: arc 0 (ec=2) + arc 1 (ec=3) = 5.
    // Pair [5]+[4]: arc 5 (ec=2) + arc 4 (ec=3) = 5.
    // Individual max = 3.  Merged max = 5.
    assert(s.simulated_contraction_weight(0) == 5);

    // c1 adj_arcs = {1,2,4,3}.
    // Pair [1]+[2]: arc 1 (ec=3) + arc 2 (ec=1) = 4.
    // Pair [4]+[3]: arc 4 (ec=3) + arc 3 (ec=1) = 4.
    // Individual max = 3.  Merged max = 4.
    assert(s.simulated_contraction_weight(1) == 4);

    std::printf("  [PASS] simulated_contraction_weight\n");
}

// ════════════════════════════════════════════════════════════════
//  4. is_granular
// ════════════════════════════════════════════════════════════════

static void test_is_granular() {
    Submap s = build_conformal_submap();
    // Weights: r0=2, r1=3, r2=1.
    // Degrees: r0=1, r1=2, r2=1.
    // Contraction weights (with merging): c0→5, c1→4.
    //
    // For γ=3: semigranular (all ≤ 3). Condition (ii):
    // c0 contraction = 5 > 3 ✓. c1 contraction = 4 > 3 ✓.
    // All degree<3 contractions exceed γ → granular!
    assert(s.is_granular(3));

    // For γ=2: not semigranular (r1 weight 3 > 2).
    assert(!s.is_granular(2));

    // For γ=4: semigranular. c0→5 > 4 ✓, c1→4 > 4? No → NOT granular.
    assert(!s.is_granular(4));

    // Single-region submap (no chords) → granular by default.
    Submap single;
    single.add_node();
    Arc a; a.first_edge = 0; a.last_edge = 2; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 3;
    single.add_arc(a);
    assert(single.is_granular(3));
    assert(single.is_granular(5));

    std::printf("  [PASS] is_granular\n");
}

// ════════════════════════════════════════════════════════════════
//  4b. is_granular — true positive (with chords)
// ════════════════════════════════════════════════════════════════

static void test_is_granular_true() {
    // r0 — c0 — r1.  γ = 3.
    // r0 has a 2-edge arc adjacent to c0.
    // r1 has a 2-edge arc adjacent to c0.
    // Individual weights: r0=2, r1=2.  Both ≤ 3 → semigranular.
    // Contraction of c0: adj arcs merge → 2+2 = 4.  4 > 3 → condition (ii) holds.
    // → granular!
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1

    Arc a;
    // r0 LEFT arc (2 edges, adjacent to c0)
    a = {}; a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = s.add_arc(a);

    // r1 LEFT arc (2 edges, adjacent to c0)
    a = {}; a.first_edge = 2; a.last_edge = 3; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 1; a.edge_count = 2;
    std::size_t ai1 = s.add_arc(a);

    // r1 RIGHT arc
    a = {}; a.first_edge = 3; a.last_edge = 2; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 2;
    std::size_t ai2 = s.add_arc(a);

    // r0 RIGHT arc
    a = {}; a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai3 = s.add_arc(a);

    Chord c;
    c.region[0] = 0; c.region[1] = 1;
    c.adj_arcs = {ai0, ai1, ai3, ai2}; c.num_adj_arcs = 4;
    s.add_chord(c);

    // Verify contraction weight accounts for merging.
    assert(s.simulated_contraction_weight(0) == 4); // 2+2 from merged adj arcs

    // γ=3: semigranular (both weights 2 ≤ 3), contraction 4 > 3 → granular.
    assert(s.is_semigranular(3));
    assert(s.is_granular(3));

    // γ=4: semigranular, but contraction 4 > 4 is false → NOT granular.
    assert(s.is_semigranular(4));
    assert(!s.is_granular(4));

    // γ=1: not semigranular (weight 2 > 1).
    assert(!s.is_granular(1));

    // No chords → granular by default.
    Submap single;
    single.add_node();
    a = {}; a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    single.add_arc(a);
    assert(single.is_granular(2));

    std::printf("  [PASS] is_granular_true\n");
}

// ════════════════════════════════════════════════════════════════
//  5. Tree decomposition
// ════════════════════════════════════════════════════════════════

static void test_tree_decomposition() {
    Submap s = build_conformal_submap();
    s.build_tree_decomposition();

    const auto& td = s.tree_decomposition();
    assert(!td.empty());

    // 3 regions + 2 chords = 5 TD nodes (2 internal + 3 leaves).
    assert(td.size() == 5);

    // Root exists and is internal (a chord).
    assert(td.root() != NONE);
    assert(td.node(td.root()).is_internal());

    // Count internal vs leaves.
    std::size_t internals = 0, leaves = 0;
    for (std::size_t i = 0; i < td.size(); ++i) {
        if (td.node(i).is_internal()) ++internals;
        else ++leaves;
    }
    assert(internals == 2);
    assert(leaves == 3);

    std::printf("  [PASS] tree_decomposition\n");
}

// ════════════════════════════════════════════════════════════════
//  6. Tree decomposition of single-region submap
// ════════════════════════════════════════════════════════════════

static void test_td_single_region() {
    Submap s;
    s.add_node();
    Arc a; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1;
    s.add_arc(a);

    s.build_tree_decomposition();
    const auto& td = s.tree_decomposition();

    assert(td.size() == 1);          // just the leaf
    assert(td.node(0).is_leaf());
    assert(td.node(0).region_idx == 0);

    std::printf("  [PASS] td_single_region\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::printf("§2.3 tests:\n");
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
