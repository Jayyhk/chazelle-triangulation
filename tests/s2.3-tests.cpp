// tests/2.3-tests.cpp — Tests for [C91 §2.3]: conformality, granularity, tree decomposition.

#include "submap/submap.h"
#include "polygon/polygon.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

static Polygon test_polygon() {
    // 4-vertex non-monotone polygon with vertices at the chord y values
    // (y=0.5 at index 2, y=1.5 at index 3).  Chord SoS tags identify
    // these source vertices per [C91 §2 tex 47].  Three edges keep region
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
    // Polygon: test_polygon() — 4 vertices, edges 0-2, all nonnull.
    //
    // Submap: r0 — c0 — r1 — c1 — r2
    // c0 at y=0.5 on edge 0 (non-vertex), c1 at y=1.5 on edge 1 (non-vertex).
    // Adjacent arcs at each chord endpoint share the junction edge
    // ([C91 §2.2 tex 94]: "glueing back ∂C at those points").
    //
    // [C91 §2.4 tex 142]: wrap-spanning arcs are single double-backing
    // structures — r2's arc passes around C's end vertex, r0's around
    // C's start vertex — so the table holds 4 arc-structures.
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1
    s.add_node(); // r2

    // [C91 §2.4(iii) tex 138]: endpoints first, so add_arc classifies
    // the wrap arcs.
    s.start_vertex = 0; s.end_vertex = 3;

    // Arc edge_count = nonnull ∂C edges per leg ([C91 §2.2 tex 106]:
    // both sides of a doubled-over C edge count; [C91 §2.4 tex 142]).
    Arc a;

    // a1: r1, LEFT edges 0..1 — count_nonnull(0,1) = 2.
    a = {}; a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 1; a.edge_count = 2;
    s.add_arc(a); // idx 0
    // E: r2's end-wrap arc — LEFT [1,2] around C's end vertex, RIGHT
    // [1,2] back down to c1's endpoint; legs 2 + 2 → count = 4.
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 2; a.edge_count = 4;
    s.add_arc(a); // idx 1
    // a4: r1, RIGHT edges 1..0 — count_nonnull(0,1) = 2.
    a = {}; a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 2;
    s.add_arc(a); // idx 2
    // S: r0's start-wrap arc — RIGHT [0,0] around C's start vertex,
    // LEFT [0,0] up to c0's endpoint; legs 1 + 1 → count = 2.
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    s.add_arc(a); // idx 3

    Chord c;
    // c0: r0 — r1.  Horizontal exit chord at y=0.5 (vertex 2's y)
    // sourced at vertex 2 per [C91 §2 tex 47]; tag identifies the source.
    // Mid-edge endpoints record {before, after} ([C91 §2.4(ii) tex 137]):
    // LEFT (0,LEFT): {S, a1}; RIGHT (0,RIGHT): {a4, S}.
    c = {}; c.region[0] = 0; c.region[1] = 1;
    c.left_edge = 0; c.right_edge = 0; c.y = 0.5; c.y_tag = 2;
    c.left_adj = {{3, 0}, 2}; c.right_adj = {{2, 3}, 2};
    s.add_chord(c);
    // c1: r1 — r2.  Horizontal at y=1.5 (vertex 3's y).
    // LEFT (1,LEFT): {a1, E}; RIGHT (1,RIGHT): {E, a4}.
    c = {}; c.region[0] = 1; c.region[1] = 2;
    c.left_edge = 1; c.right_edge = 1; c.y = 1.5; c.y_tag = 3;
    c.left_adj = {{0, 1}, 2}; c.right_adj = {{1, 2}, 2};
    s.add_chord(c);

    // [C91 §2.4(iii) tex 138]: add_arc auto-registered the endpoint
    // pointers (start_arc = S, the last table entry; end_arc = E).
    assert(s.start_arc == 3 && s.end_arc == 1);
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
    c.left_adj = {{0}, 1}; c.right_adj = {{0}, 1};
    s.add_chord(c);
    c = {}; c.region[0] = 1; c.region[1] = 4; c.y_tag = 0;
    c.left_adj = {{0}, 1}; c.right_adj = {{0}, 1};
    s.add_chord(c);
    c = {}; c.region[0] = 1; c.region[1] = 5; c.y_tag = 0;
    c.left_adj = {{0}, 1}; c.right_adj = {{0}, 1};
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
    // Weights (paper-strict, [C91 §2.2 tex 106]): r0=2, r1=2, r2=4.

    assert(s.is_semigranular(4));  // all ≤ 4
    assert(s.is_semigranular(5));  // all ≤ 5
    assert(!s.is_semigranular(3)); // r2 has weight 4 > 3
    assert(!s.is_semigranular(1)); // r0, r1 have weight 2 > 1

    std::printf("  [PASS] is_semigranular\n");
}

// ════════════════════════════════════════════════════════════════
//  3. simulated_contraction_weight
// ════════════════════════════════════════════════════════════════

static void test_simulated_contraction_weight() {
    Submap s = build_conformal_submap();

    // c0 contraction merges r0 (weight 2) and r1 (weight 2): the
    // glue chains run S+a1 and a4+S — shared arc S links them into one
    // chain {a4, S, a1} ([C91 §2.2 tex 94]) = a start-wrap arc
    // R(1)→L(1); legs RIGHT[0,1] + LEFT[0,1] → 4 ([C91 §2.2 tex 106]).
    assert(s.simulated_contraction_weight(0, test_polygon()) == 4);

    // c1 contraction merges r1 (weight 2) and r2 (weight 2): shared
    // arc E links the glue pairs into the chain {a1, E, a4}
    // ([C91 §2.2 tex 94]) = an end-wrap arc L(0)→R(0); legs
    // LEFT[0,2] + RIGHT[0,2] → 6 ([C91 §2.2 tex 106]).
    assert(s.simulated_contraction_weight(1, test_polygon()) == 6);

    std::printf("  [PASS] simulated_contraction_weight\n");
}

// ════════════════════════════════════════════════════════════════
//  4. is_granular
// ════════════════════════════════════════════════════════════════

static void test_is_granular() {
    Submap s = build_conformal_submap();
    // Weights (paper-strict, [C91 §2.2 tex 106]): r0=2, r1=2, r2=4.
    // Degrees: r0=1, r1=2, r2=1.
    // Contraction weights: c0→4, c1→6.
    //
    // For γ=3: NOT semigranular (r2=4 > 3) → NOT granular.
    assert(!s.is_granular(3, test_polygon()));

    // For γ=4: semigranular ✓.  c0=4 NOT > 4 → NOT granular.
    assert(!s.is_granular(4, test_polygon()));

    // For γ=6: semigranular ✓.  c1=6 NOT > 6 → NOT granular.
    assert(!s.is_granular(6, test_polygon()));

    // Chordless submap (single region bounded by the closed arc,
    // [C91 §2.4 tex 142/138]) → granular by default ([C91 §2.3 tex 123]).
    Submap single;
    single.add_node();
    single.start_vertex = 0; single.end_vertex = 3;
    Arc a; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 6;   // closed arc: both ∂C sides
    single.add_arc(a);
    assert(single.is_granular(6, test_polygon()));
    assert(single.is_granular(8, test_polygon()));

    std::printf("  [PASS] is_granular\n");
}

// ════════════════════════════════════════════════════════════════
//  4b. is_granular — true positive (with chords)
// ════════════════════════════════════════════════════════════════

static void test_is_granular_true() {
    // r0 — c0 — r1.
    // c0 at y=1.5 on edge 1 (non-vertex).
    // Adjacent arcs share junction edge 1 at both chord endpoints.
    // Individual weights ([C91 §2.2 tex 106]): r0=4, r1=4 (each wrap
    // arc covers both ∂C sides of two edges).  Contracting c0 closes
    // ∂C → the closed arc weighs 2×3 = 6.  6 > 4 → granular for γ=4.
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1
    s.start_vertex = 0; s.end_vertex = 3;

    // [C91 §2.4 tex 142]: r1's arc double-backs around C's end vertex,
    // r0's around C's start vertex — one structure each.
    Arc a;
    // E: r1's end-wrap arc — LEFT [1,2] + RIGHT [1,2]; legs 2+2 → 4.
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 4;
    std::size_t E = s.add_arc(a);

    // S: r0's start-wrap arc — RIGHT [0,1] + LEFT [0,1]; legs 2+2 → 4.
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 4;
    std::size_t S = s.add_arc(a);

    Chord c;
    c.region[0] = 0; c.region[1] = 1;
    // c0 at y=1.5 = vertex 3's y (per test_polygon); tag identifies
    // source.  Mid-edge endpoints ([C91 §2.4(ii) tex 137]): LEFT
    // (1,LEFT): {before=S, after=E}; RIGHT (1,RIGHT): {before=E, after=S}.
    c.left_edge = 1; c.right_edge = 1; c.y = 1.5; c.y_tag = 3;
    c.left_adj = {{S, E}, 2}; c.right_adj = {{E, S}, 2};
    s.add_chord(c);

    assert(s.start_arc == S && s.end_arc == E &&
           "[C91 §2.4(iii) tex 138]: endpoint arcs auto-registered");
    s.check_invariants(test_polygon());

    // [C91 §2.2 tex 94 / §2.4 tex 142]: contracting the LAST chord
    // closes ∂C — the merged weight is the nonnull ∂C count of BOTH
    // sides of C ([C91 §2.2 tex 106]).
    assert(s.simulated_contraction_weight(0, test_polygon()) == 6);

    // γ=4: semigranular (both weights 4 ≤ 4), contraction 6 > 4 → granular.
    assert(s.is_semigranular(4));
    assert(s.is_granular(4, test_polygon()));

    // γ=6: semigranular, but contraction 6 > 6 is false → NOT granular.
    assert(s.is_semigranular(6));
    assert(!s.is_granular(6, test_polygon()));

    // γ=3: not semigranular (weight 4 > 3).
    assert(!s.is_granular(3, test_polygon()));

    // No chords → granular by default (single closed arc covering all
    // of C = vertices 0..2, [C91 §2.4 tex 142/138]).
    Submap single;
    single.add_node();
    single.start_vertex = 0; single.end_vertex = 2;
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 4;   // closed arc: both ∂C sides
    single.add_arc(a);
    assert(single.is_granular(4, test_polygon()));

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

    // 3 regions + 2 chords = 5 tree decomposition nodes (2 internal + 3 leaves).
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
    s.start_vertex = 0; s.end_vertex = 3;
    // The chordless submap's single closed arc ([C91 §2.4 tex 142/138]).
    Arc a; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 3;
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
