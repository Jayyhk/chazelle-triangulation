// tests/2.3-tests.cpp — Tests for §2.3: conformality, granularity, tree decomposition.

#include "submap/submap.h"
#include "polygon/polygon.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

static Polygon test_polygon() {
    // 4-vertex non-monotone polygon with vertices at the chord y values
    // (y=0.5 at index 2, y=1.5 at index 3).  Chord SoS tags identify
    // these source vertices per §2.1 tex 70.  Three edges keep region
    // weights bounded so the γ=2 granularity demo below remains achievable.
    return Polygon({
        {0, 0,   0},
        {1, 1,   1},
        {2, 0.5, 2},
        {3, 1.5, 3}
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

    // Polygon has 4 nonnull edges (all consecutive vertex pairs differ in
    // both coords).  Arc edge_count = polygon.count_nonnull_edges(lo, hi)
    // over the arc's underlying edge range, per §2.2 tex 106 cache invariant.
    Arc a;

    // LEFT half (ascending first_edge in ∂C order):
    // a0: r0, edge 0 only — count_nonnull(0,0) = 1.
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1;
    a.key_y = 0.0; a.key_y_tag = 0; // C-start at vertex 0
    s.add_arc(a); // idx 0
    // a1: r1, edges 0..1 — count_nonnull(0,1) = 2.
    a = {}; a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 1; a.edge_count = 2;
    a.key_y = 0.5; a.key_y_tag = 2; // vertex 2 (y=0.5)
    s.add_arc(a); // idx 1
    // a2: r2, edges 1..2 (last LEFT arc reaches end_vertex=3) —
    //     count_nonnull(1,2) = 2.
    a = {}; a.first_edge = 1; a.last_edge = 2; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 2; a.edge_count = 2;
    a.key_y = 1.5; a.key_y_tag = 3; // vertex 3 (y=1.5)
    s.add_arc(a); // idx 2

    // RIGHT half (descending first_edge in ∂C order):
    // a3: r2, edges 1..2 — count_nonnull(1,2) = 2.
    a = {}; a.first_edge = 2; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 2; a.edge_count = 2;
    a.key_y = 1.5; a.key_y_tag = 3; // C-end at vertex 3
    s.add_arc(a); // idx 3
    // a4: r1, edges 0..1 — count_nonnull(0,1) = 2.
    a = {}; a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 2;
    a.key_y = 1.5; a.key_y_tag = 3; // r1 RIGHT arc starts at chord c1 (vertex 3 y)
    s.add_arc(a); // idx 4
    // a5: r0, edge 0 only — count_nonnull(0,0) = 1.
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 1;
    a.key_y = 0.5; a.key_y_tag = 2; // r0 RIGHT arc starts at chord c0 (vertex 2 y)
    s.add_arc(a); // idx 5

    Chord c;
    // c0: r0 — r1.  Horizontal at y=0.5 (vertex 2's y) — non-NLC chord
    // sourced at vertex 2 per §2.1 tex 70; tag identifies the source.
    c = {}; c.region[0] = 0; c.region[1] = 1;
    c.left_edge = 0; c.right_edge = 0; c.y = 0.5; c.y_tag = 2;
    c.left_adj = {{0, 1}, 2}; c.right_adj = {{4, 5}, 2};
    s.add_chord(c);
    // c1: r1 — r2.  Horizontal at y=1.5 (vertex 3's y).
    c = {}; c.region[0] = 1; c.region[1] = 2;
    c.left_edge = 1; c.right_edge = 1; c.y = 1.5; c.y_tag = 3;
    c.left_adj = {{1, 2}, 2}; c.right_adj = {{3, 4}, 2};
    s.add_chord(c);

    // §2.4(iii) tex 138: endpoint pointers (LEFT-side, at C's two ends).
    s.start_arc = 0; s.end_arc = 2;
    s.start_vertex = 0; s.end_vertex = 3;
    s.check_invariants(test_polygon());

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
    c = {}; c.region[0] = 1; c.region[1] = 3; c.y_tag = 0;
    c.left_adj = {{1}, 1}; c.right_adj = {{1}, 1};
    s.add_chord(c);
    c = {}; c.region[0] = 1; c.region[1] = 4; c.y_tag = 0;
    c.left_adj = {{1}, 1}; c.right_adj = {{1}, 1};
    s.add_chord(c);
    c = {}; c.region[0] = 1; c.region[1] = 5; c.y_tag = 0;
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
    // Weights (paper-strict): r0=1, r1=2, r2=2.

    assert(s.is_semigranular(2));  // all ≤ 2
    assert(s.is_semigranular(3));  // all ≤ 3
    assert(!s.is_semigranular(1)); // r1, r2 have weight 2 > 1

    std::printf("  [PASS] is_semigranular\n");
}

// ════════════════════════════════════════════════════════════════
//  3. simulated_contraction_weight
// ════════════════════════════════════════════════════════════════

static void test_simulated_contraction_weight() {
    Submap s = build_conformal_submap();

    // c0 contraction merges r0 (weight 1) and r1 (weight 2):
    //   left: a0(ec=1) + a1(ec=2) - shared_nonnull(edge 0)=1 = 2.
    //   right: a4(ec=2) + a5(ec=1) - 1 = 2.
    //   max(initial=2, 2, 2) = 2.
    assert(s.simulated_contraction_weight(0, test_polygon()) == 2);

    // c1 contraction merges r1 (weight 2) and r2 (weight 2):
    //   left: a1(ec=2) + a2(ec=2) - shared_nonnull(edge 1)=1 = 3.
    //   right: a3(ec=2) + a4(ec=2) - 1 = 3.
    //   max(initial=2, 3, 3) = 3.
    assert(s.simulated_contraction_weight(1, test_polygon()) == 3);

    std::printf("  [PASS] simulated_contraction_weight\n");
}

// ════════════════════════════════════════════════════════════════
//  4. is_granular
// ════════════════════════════════════════════════════════════════

static void test_is_granular() {
    Submap s = build_conformal_submap();
    // Weights (paper-strict): r0=1, r1=2, r2=2.
    // Degrees: r0=1, r1=2, r2=1.
    // Contraction weights: c0→2, c1→3.
    //
    // For γ=2: semigranular ✓.  c0=2 NOT > 2 → NOT granular.
    assert(!s.is_granular(2, test_polygon()));

    // For γ=3: semigranular ✓.  c0=2 NOT > 3, c1=3 NOT > 3 → NOT granular.
    assert(!s.is_granular(3, test_polygon()));

    // For γ=4: semigranular ✓.  c0=2 NOT > 4 → NOT granular.
    assert(!s.is_granular(4, test_polygon()));

    // Single-region submap (no chords) → granular by default.
    Submap single;
    single.add_node();
    Arc a; a.first_edge = 0; a.last_edge = 2; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 3;
    a.key_y_tag = 0;
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
    a.key_y = 0.0; a.key_y_tag = 0; // C-start at vertex 0
    std::size_t ai0 = s.add_arc(a);

    // r1 LEFT arc: starts at c0 on edge 1, edges [1,2]
    a = {}; a.first_edge = 1; a.last_edge = 2; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 1; a.edge_count = 2;
    a.key_y = 1.5; a.key_y_tag = 3; // chord c0 source at vertex 3
    std::size_t ai1 = s.add_arc(a);

    // r1 RIGHT arc: edges [2,1]
    a = {}; a.first_edge = 2; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 2;
    a.key_y = 1.5; a.key_y_tag = 3; // C-end at vertex 3
    std::size_t ai2 = s.add_arc(a);

    // r0 RIGHT arc: starts at c0 on edge 1, edges [1,0]
    a = {}; a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    a.key_y = 1.5; a.key_y_tag = 3; // chord c0 source at vertex 3
    std::size_t ai3 = s.add_arc(a);

    Chord c;
    c.region[0] = 0; c.region[1] = 1;
    // c0 at y=1.5 = vertex 3's y (per test_polygon); tag identifies source.
    c.left_edge = 1; c.right_edge = 1; c.y = 1.5; c.y_tag = 3;
    c.left_adj = {{ai0, ai1}, 2}; c.right_adj = {{ai2, ai3}, 2};
    s.add_chord(c);

    // §2.4(iii) tex 138: endpoint pointers (LEFT-side, at C's two ends).
    s.start_arc = ai0; s.end_arc = ai1;
    s.start_vertex = 0; s.end_vertex = 3;
    s.check_invariants(test_polygon());

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
    a.key_y_tag = 0;
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
    a.key_y_tag = 0;
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
