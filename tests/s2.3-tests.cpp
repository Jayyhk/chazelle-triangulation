/// tests/2.3-tests.cpp — Tests for §2.3: conformality, granularity, tree decomposition.

#include "submap/submap.h"
#include "polygon/polygon.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

static Polygon test_polygon() {
    return Polygon({
        {0,0,0}, {1,1,1}, {2,2,2}, {3,3,3}, {4,4,4}
    });
}

// ════════════════════════════════════════════════════════════════
//  Helper: build a conformal submap (degree ≤ 4).
//  Linear chain: r0 — c0 — r1 — c1 — r2
// ════════════════════════════════════════════════════════════════

static Submap build_conformal_submap() {
    // Polygon: {0,0,0},{1,1,1},{2,2,2},{3,3,3},{4,4,4}
    // Edges 0-3, all nonnull.
    //
    // Submap: r0 — c0 — r1 — c1 — r2
    // c0 at y=0.5 on edge 0 (non-vertex), c1 at y=1.5 on edge 1 (non-vertex).
    // Adjacent arcs at each chord endpoint share the junction edge
    // (§2.2 tex 94: "glueing back ∂C at those points").
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1
    s.add_node(); // r2

    Arc a;

    // LEFT half (ascending first_edge in ∂C order):
    // a0: r0, edge 0 only
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    s.add_arc(a); // idx 0
    // a1: r1, starts at c0 on edge 0, ends at c1 on edge 1
    a = {}; a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 1; a.edge_count = 3; a.key_y = 0.5;
    s.add_arc(a); // idx 1
    // a2: r2, starts at c1 on edge 1
    a = {}; a.first_edge = 1; a.last_edge = 2; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 2; a.edge_count = 1; a.key_y = 1.5;
    s.add_arc(a); // idx 2

    // RIGHT half (descending first_edge in ∂C order):
    // a3: r2
    a = {}; a.first_edge = 2; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 2; a.edge_count = 1;
    s.add_arc(a); // idx 3
    // a4: r1, starts at c1 on edge 1, ends at c0 on edge 0
    a = {}; a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 3;
    s.add_arc(a); // idx 4
    // a5: r0, starts at c0 on edge 0
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    s.add_arc(a); // idx 5

    Chord c;
    // c0: r0 — r1, on edge 0 at y=0.5 (non-vertex)
    c = {}; c.region[0] = 0; c.region[1] = 1;
    c.left_edge = 0; c.right_edge = 0; c.y = 0.5;
    c.left_adj = {{0, 1}, 2}; c.right_adj = {{4, 5}, 2};
    s.add_chord(c);
    // c1: r1 — r2, on edge 1 at y=1.5 (non-vertex)
    c = {}; c.region[0] = 1; c.region[1] = 2;
    c.left_edge = 1; c.right_edge = 1; c.y = 1.5;
    c.left_adj = {{1, 2}, 2}; c.right_adj = {{3, 4}, 2};
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
    c.left_adj = {{1}, 1}; c.right_adj = {{1}, 1};
    s.add_chord(c);
    c = {}; c.region[0] = 1; c.region[1] = 4;
    c.left_adj = {{1}, 1}; c.right_adj = {{1}, 1};
    s.add_chord(c);
    c = {}; c.region[0] = 1; c.region[1] = 5;
    c.left_adj = {{1}, 1}; c.right_adj = {{1}, 1};
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
    // Pair [0]+[1]: arc 0 (ec=2) + arc 1 (ec=3) - shared(1) = 4.
    // Pair [5]+[4]: arc 5 (ec=2) + arc 4 (ec=3) - shared(1) = 4.
    // Individual max = 3.  Merged max = 4.
    assert(s.simulated_contraction_weight(0, test_polygon()) == 4);

    // c1 adj_arcs = {1,2,4,3}.
    // Pair [1]+[2]: arc 1 (ec=3) + arc 2 (ec=1) - shared(1) = 3.
    // Pair [4]+[3]: arc 4 (ec=3) + arc 3 (ec=1) - shared(1) = 3.
    // Individual max = 3.  Merged max = 3.
    assert(s.simulated_contraction_weight(1, test_polygon()) == 3);

    std::printf("  [PASS] simulated_contraction_weight\n");
}

// ════════════════════════════════════════════════════════════════
//  4. is_granular
// ════════════════════════════════════════════════════════════════

static void test_is_granular() {
    Submap s = build_conformal_submap();
    // Weights: r0=2, r1=3, r2=1.
    // Degrees: r0=1, r1=2, r2=1.
    // Contraction weights (with shared edge deduction): c0→4, c1→3.
    //
    // For γ=2: not semigranular (r1 weight 3 > 2).
    assert(!s.is_granular(2, test_polygon()));

    // For γ=3: semigranular (all ≤ 3). Condition (ii):
    // c0 contraction = 4 > 3 ✓. c1 contraction = 3 > 3? No.
    // NOT granular (c1 is removable).
    assert(!s.is_granular(3, test_polygon()));

    // For γ=4: semigranular. c0→4 > 4? No → NOT granular.
    assert(!s.is_granular(4, test_polygon()));

    // Single-region submap (no chords) → granular by default.
    Submap single;
    single.add_node();
    Arc a; a.first_edge = 0; a.last_edge = 2; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 3;
    single.add_arc(a);
    assert(single.is_granular(3, test_polygon()));
    assert(single.is_granular(5, test_polygon()));

    std::printf("  [PASS] is_granular\n");
}

// ════════════════════════════════════════════════════════════════
//  4b. is_granular — true positive (with chords)
// ════════════════════════════════════════════════════════════════

static void test_is_granular_true() {
    // r0 — c0 — r1.  γ = 3.
    // c0 at y=1.5 on edge 1 (non-vertex).
    // Adjacent arcs share junction edge 1 at both chord endpoints.
    // Individual weights: r0=2, r1=2.  Both ≤ 3 → semigranular.
    // Contraction of c0: adj arcs merge → 2+2-1 = 3.  3 > 2 → granular for γ=2.
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1

    Arc a;
    // r0 LEFT arc: edges [0,1]
    a = {}; a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = s.add_arc(a);

    // r1 LEFT arc: starts at c0 on edge 1, edges [1,2]
    a = {}; a.first_edge = 1; a.last_edge = 2; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 1; a.edge_count = 2; a.key_y = 1.5;
    std::size_t ai1 = s.add_arc(a);

    // r1 RIGHT arc: edges [2,1]
    a = {}; a.first_edge = 2; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 2;
    std::size_t ai2 = s.add_arc(a);

    // r0 RIGHT arc: starts at c0 on edge 1, edges [1,0]
    a = {}; a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai3 = s.add_arc(a);

    Chord c;
    c.region[0] = 0; c.region[1] = 1;
    c.left_edge = 1; c.right_edge = 1; c.y = 1.5;
    c.left_adj = {{ai0, ai1}, 2}; c.right_adj = {{ai2, ai3}, 2};
    s.add_chord(c);

    // Verify contraction weight accounts for merging with
    // shared edge deduction: 2+2-1 = 3.
    assert(s.simulated_contraction_weight(0, test_polygon()) == 3);

    // γ=2: semigranular (both weights 2 ≤ 2), contraction 3 > 2 → granular.
    assert(s.is_semigranular(2));
    assert(s.is_granular(2, test_polygon()));

    // γ=3: semigranular, but contraction 3 > 3 is false → NOT granular.
    assert(s.is_semigranular(3));
    assert(!s.is_granular(3, test_polygon()));

    // γ=1: not semigranular (weight 2 > 1).
    assert(!s.is_granular(1, test_polygon()));

    // No chords → granular by default.
    Submap single;
    single.add_node();
    a = {}; a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    single.add_arc(a);
    assert(single.is_granular(2, test_polygon()));

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
