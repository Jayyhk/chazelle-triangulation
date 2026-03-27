/// tests/2.2-tests.cpp — Tests for §2.2: submaps, weight, chord removal.

#include "polygon/polygon.h"
#include "submap/submap.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

// A polygon with 5 vertices (4 edges) for chord removal tests.
// Chord endpoints in the test submaps are at edges 0-3.
static Polygon test_polygon() {
    return Polygon({
        {0,0,0}, {1,1,1}, {2,2,2}, {3,3,3}, {4,4,4}
    });
}

// ════════════════════════════════════════════════════════════════
//  Helper: build a simple submap with 3 regions and 2 chords.
//
//  Region 0 — chord 0 — Region 1 — chord 1 — Region 2
//
//  Arcs (LEFT side): a0(reg 0), a1(reg 1), a2(reg 2)
//  Arcs (RIGHT side): a3(reg 2), a4(reg 1), a5(reg 0)
// ════════════════════════════════════════════════════════════════

static Submap build_3region_submap() {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();

    // LEFT arcs (∂C order: a0, a1, a2)
    Arc a0; a0.first_edge = 0; a0.last_edge = 0; a0.first_side = LEFT; a0.last_side = LEFT;
    a0.region_node = r0; a0.edge_count = 1; a0.key_y = 0.0; a0.key_y_tag = 0;
    std::size_t ai0 = s.add_arc(a0);

    Arc a1; a1.first_edge = 1; a1.last_edge = 1; a1.first_side = LEFT; a1.last_side = LEFT;
    a1.region_node = r1; a1.edge_count = 1; a1.key_y = 1.0; a1.key_y_tag = 1;
    std::size_t ai1 = s.add_arc(a1);

    Arc a2; a2.first_edge = 2; a2.last_edge = 3; a2.first_side = LEFT; a2.last_side = LEFT;
    a2.region_node = r2; a2.edge_count = 2; a2.key_y = 2.0; a2.key_y_tag = 2;
    std::size_t ai2 = s.add_arc(a2);

    // RIGHT arcs (∂C order: a3, a4, a5)
    Arc a3; a3.first_edge = 3; a3.last_edge = 2; a3.first_side = RIGHT; a3.last_side = RIGHT;
    a3.region_node = r2; a3.edge_count = 2; a3.key_y = 3.0; a3.key_y_tag = 3;
    std::size_t ai3 = s.add_arc(a3);

    Arc a4; a4.first_edge = 1; a4.last_edge = 1; a4.first_side = RIGHT; a4.last_side = RIGHT;
    a4.region_node = r1; a4.edge_count = 1; a4.key_y = 1.0; a4.key_y_tag = 1;
    std::size_t ai4 = s.add_arc(a4);

    Arc a5; a5.first_edge = 0; a5.last_edge = 0; a5.first_side = RIGHT; a5.last_side = RIGHT;
    a5.region_node = r0; a5.edge_count = 1; a5.key_y = 0.0; a5.key_y_tag = 0;
    std::size_t ai5 = s.add_arc(a5);

    // Chord 0: between r0 and r1, adj arcs a0,a1 (LEFT) and a5,a4 (RIGHT)
    Chord c0;
    c0.region[0] = r0; c0.region[1] = r1;
    c0.left_adj = {{ai0, ai1}, 2}; c0.right_adj = {{ai5, ai4}, 2};
    c0.left_edge = 0; c0.right_edge = 0; c0.y = 1.0; c0.y_tag = 1;
    s.add_chord(c0);

    // Chord 1: between r1 and r2, adj arcs a1,a2 (LEFT) and a4,a3 (RIGHT)
    Chord c1;
    c1.region[0] = r1; c1.region[1] = r2;
    c1.left_adj = {{ai1, ai2}, 2}; c1.right_adj = {{ai4, ai3}, 2};
    c1.left_edge = 1; c1.right_edge = 1; c1.y = 2.0; c1.y_tag = 2;
    s.add_chord(c1);

    s.start_arc = ai0;
    s.end_arc = ai2;
    s.start_vertex = 0;
    s.end_vertex = 4;

    return s;
}

// ════════════════════════════════════════════════════════════════
//  1. count_nonnull_edges (Polygon)
// ════════════════════════════════════════════════════════════════

static void test_count_nonnull_edges() {
    Polygon tri({{0,0,0}, {4,1,1}, {2,3,2}});
    assert(tri.count_nonnull_edges(0, 0) == 1);
    assert(tri.count_nonnull_edges(1, 1) == 1);
    assert(tri.count_nonnull_edges(0, 1) == 2);

    // Zero-length edge.
    Polygon p({{0,0,0}, {3,3,1}, {3,3,2}, {5,1,3}});
    assert(p.count_nonnull_edges(1, 1) == 0);
    assert(p.count_nonnull_edges(0, 2) == 2);

    std::printf("  [PASS] count_nonnull_edges\n");
}

// ════════════════════════════════════════════════════════════════
//  2. Submap construction + tree property
// ════════════════════════════════════════════════════════════════

static void test_submap_construction() {
    Submap s = build_3region_submap();

    assert(s.num_nodes() == 3);
    assert(s.num_chords() == 2);
    assert(s.num_arcs() == 6);

    // Tree property: 3 regions = 2 chords + 1.
    s.assert_tree_property();

    // Node degrees.
    assert(s.node(0).degree() == 1); // r0: chord 0
    assert(s.node(1).degree() == 2); // r1: chord 0, chord 1
    assert(s.node(2).degree() == 1); // r2: chord 1

    std::printf("  [PASS] submap_construction\n");
}

// ════════════════════════════════════════════════════════════════
//  3. check_invariants
// ════════════════════════════════════════════════════════════════

static void test_check_invariants() {
    Submap s = build_3region_submap();
    s.check_invariants(); // should not fire

    std::printf("  [PASS] check_invariants\n");
}

// ════════════════════════════════════════════════════════════════
//  4. region_weight
// ════════════════════════════════════════════════════════════════

static void test_region_weight() {
    Submap s = build_3region_submap();

    // r0 has arcs a0(ec=1) and a5(ec=1) → weight = 1
    assert(s.region_weight(0) == 1);
    // r1 has arcs a1(ec=1) and a4(ec=1) → weight = 1
    assert(s.region_weight(1) == 1);
    // r2 has arcs a2(ec=2) and a3(ec=2) → weight = 2
    assert(s.region_weight(2) == 2);

    std::printf("  [PASS] region_weight\n");
}

// ════════════════════════════════════════════════════════════════
//  5. remove_chord (region merge + arc merge)
// ════════════════════════════════════════════════════════════════

static void test_remove_chord() {
    Submap s = build_3region_submap();

    // Remove chord 1 (between r1 and r2).
    auto poly = test_polygon();
    [[maybe_unused]] std::size_t survivor = s.remove_chord(1, poly);

    // 2 live regions, 1 live chord remain (dead entries tombstoned).
    assert(s.num_live_nodes() == 2);
    assert(s.num_live_chords() == 1);
    s.assert_tree_property();

    std::printf("  [PASS] remove_chord\n");
}

// ════════════════════════════════════════════════════════════════
//  6. remove all chords → single region
// ════════════════════════════════════════════════════════════════

static void test_remove_all_chords() {
    Submap s = build_3region_submap();

    auto poly = test_polygon();
    s.remove_chord(0, poly);  // tombstones chord 0
    s.remove_chord(1, poly);  // tombstones chord 1 (stable index)

    // 1 live region, 0 live chords.
    s.assert_tree_property();

    assert(s.num_live_nodes() == 1);
    assert(s.num_live_chords() == 0);

    // After compact, raw counts match live counts.
    s.compact();
    assert(s.num_nodes() == 1);
    assert(s.num_chords() == 0);

    std::printf("  [PASS] remove_all_chords\n");
}

// ════════════════════════════════════════════════════════════════
//  7. Empty region weight = 0
// ════════════════════════════════════════════════════════════════

static void test_empty_region_weight() {
    // [C91 §2.2]: "weight of a region as 0 if the region is empty."
    Submap s;
    std::size_t r0 = s.add_node();
    // r0 has no arcs → weight = 0.
    assert(s.region_weight(r0) == 0);

    // Add a zero-length arc (edge_count = 0) → still weight 0.
    Arc a;
    a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 0;
    s.add_arc(a);
    assert(s.region_weight(r0) == 0);

    std::printf("  [PASS] empty_region_weight\n");
}

// ════════════════════════════════════════════════════════════════
//  8. remove_chord arc merging verification
// ════════════════════════════════════════════════════════════════

static void test_remove_chord_4_adj_arcs() {
    // Build a submap where the chord endpoint is NOT a polygon vertex.
    // The chord is at y=1.5 (between polygon vertices at y=1 and y=2).
    // Removing it should merge the adjacent arcs.
    //
    // [C91 §2.4(iii)]: Arcs partition ∂C contiguously.  At a non-vertex
    // chord endpoint interior to edge 1, both adjacent LEFT arcs share
    // edge 1 (a_keep.last_edge = 1, a_dead.first_edge = 1).  Similarly
    // for the RIGHT pair.
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1

    Arc a;
    // r0 LEFT arc: edges 0–1 (up to chord at y=1.5 on edge 1)
    a = {}; a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = s.add_arc(a);
    // r1 LEFT arc: edges 1–3 (from chord at y=1.5 on edge 1 onward)
    a = {}; a.first_edge = 1; a.last_edge = 3; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 1; a.edge_count = 3;
    std::size_t ai1 = s.add_arc(a);
    // r1 RIGHT arc: edges 3–1 (descending, from end down to chord)
    a = {}; a.first_edge = 3; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 3;
    std::size_t ai2 = s.add_arc(a);
    // r0 RIGHT arc: edges 1–0 (descending, from chord down to start)
    a = {}; a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai3 = s.add_arc(a);

    Chord c;
    c.region[0] = 0; c.region[1] = 1;
    // LEFT endpoint (edge 1): ai0 ends at chord, ai1 starts at chord.
    // RIGHT endpoint (edge 1): ai2 ends at chord, ai3 starts at chord.
    c.left_adj = {{ai0, ai1}, 2}; c.right_adj = {{ai2, ai3}, 2};
    // Chord at y=1.5 — NOT a polygon vertex (polygon vertices are at y=0,1,2,3,4).
    c.left_edge = 1; c.right_edge = 1;
    c.y = 1.5; c.y_tag = 99; // tag 99 doesn't match any vertex index
    s.add_chord(c);

    auto poly = test_polygon();
    s.remove_chord(0, poly);

    // After removal: adj arc pairs should merge.  Dead arcs are
    // tombstoned; 2 live arcs remain.
    assert(s.num_live_arcs() == 2 && "4 adj_arcs: 2 merges → 2 live arcs");
    // LEFT merge: ai0(ec=2) + ai1(ec=3) - 1(shared edge 1) = 4.
    // RIGHT merge: ai2(ec=3) + ai3(ec=2) - 1(shared edge 1) = 4.
    bool found_merged = false;
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        if (s.arc(i).edge_count == 4 && s.arc(i).first_side == LEFT)
            found_merged = true;
    }
    assert(found_merged && "arc merging should produce merged arc "
           "(2+3-1=4 nonnull edges)");

    std::printf("  [PASS] remove_chord_4_adj_arcs\n");
}

// ════════════════════════════════════════════════════════════════
//  8b. remove_chord — NO merge at vertex endpoints
// ════════════════════════════════════════════════════════════════

static void test_remove_chord_no_merge_at_vertex() {
    // Chord endpoint IS a polygon vertex → arcs should NOT merge.
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1

    Arc a;
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1;
    std::size_t ai0 = s.add_arc(a);
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 1; a.edge_count = 1;
    std::size_t ai1 = s.add_arc(a);
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 1;
    std::size_t ai2 = s.add_arc(a);
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 1;
    std::size_t ai3 = s.add_arc(a);

    Chord c;
    c.region[0] = 0; c.region[1] = 1;
    c.left_adj = {{ai0, ai1}, 2}; c.right_adj = {{ai3, ai2}, 2};
    // Chord at y=1.0, tag=1 — matches polygon vertex 1 at (1,1,1).
    c.left_edge = 0; c.right_edge = 0;
    c.y = 1.0; c.y_tag = 1;
    s.add_chord(c);

    auto poly = test_polygon();
    s.remove_chord(0, poly);

    // No arcs should have been killed (vertex endpoint → no merge).
    assert(s.num_live_arcs() == 4 &&
           "vertex endpoint should NOT merge arcs");
    // All live arcs still have edge_count = 1 (no merging happened).
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        assert(s.arc(i).edge_count == 1);
    }

    std::printf("  [PASS] remove_chord_no_merge_at_vertex\n");
}


// ════════════════════════════════════════════════════════════════
//  8c. remove_chord with 2 adj_arcs (chord at C's start vertex)
// ════════════════════════════════════════════════════════════════

static void test_remove_chord_2_adj_arcs() {
    // [C91 §2.4(ii)]: "two... arcs adjacent to it."
    // Total adj_arcs = 2 happens only when both endpoints are at
    // C's start or end vertex, where only 1 arc exists per side.
    // Both endpoints are polygon vertices → no merging.
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();

    Arc a;
    // Region r0: arcs near C's start (edges 0-1)
    a = {}; a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 2; a.key_y = 0.0; a.key_y_tag = 0;
    std::size_t ai0 = s.add_arc(a);

    // Region r1: arcs away from C's start (edges 2-3)
    a = {}; a.first_edge = 2; a.last_edge = 3; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 2; a.key_y = 2.0; a.key_y_tag = 2;
    std::size_t ai1 = s.add_arc(a);

    // RIGHT arcs (descending first_edge)
    a = {}; a.first_edge = 3; a.last_edge = 2; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r1; a.edge_count = 2; a.key_y = 3.0; a.key_y_tag = 3;
    s.add_arc(a);

    a = {}; a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r0; a.edge_count = 2; a.key_y = 1.0; a.key_y_tag = 1;
    std::size_t ai3 = s.add_arc(a);

    // Chord at vertex 0 (C's start): y=0, tag=0 matches vertex 0.
    // LEFT endpoint on edge 0: only a0 starts here → count=1.
    // RIGHT endpoint on edge 0: only a3 ends here → count=1.
    Chord c;
    c.region[0] = r0; c.region[1] = r1;
    c.left_adj = {{ai0}, 1}; c.right_adj = {{ai3}, 1};
    c.left_edge = 0; c.right_edge = 0;
    c.left_side = LEFT; c.right_side = RIGHT;
    c.y = 0.0; c.y_tag = 0;
    s.add_chord(c);

    s.start_arc = ai0;
    s.end_arc = ai1;
    s.start_vertex = 0;
    s.end_vertex = 4;

    auto poly = test_polygon();
    s.remove_chord(0, poly);

    // Both endpoints are polygon vertices → no arcs merged.
    assert(s.num_live_arcs() == 4 &&
           "2 adj_arcs at C's start: no arc should be killed");
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        assert(s.arc(i).edge_count == 2);
    }
    s.assert_tree_property();

    std::printf("  [PASS] remove_chord_2_adj_arcs\n");
}

// ════════════════════════════════════════════════════════════════
//  8d. remove_chord with 3 adj_arcs (1 at C's start, 2 at non-vertex)
// ════════════════════════════════════════════════════════════════

// Non-monotone polygon: y=0 at v0 also intersects edge 2 on the
// RIGHT side (interior), giving count=1 on LEFT + count=2 on RIGHT.
static Polygon non_monotone_polygon() {
    return Polygon({
        {0, 0, 0}, {10, 10, 1}, {20, -5, 2}, {30, 5, 3}, {40, 15, 4}
    });
}

static void test_remove_chord_3_adj_arcs() {
    // [C91 §2.4(ii)]: "three... arcs adjacent to it."
    // LEFT endpoint at vertex 0 (C's start) → count=1, vertex, no merge.
    // RIGHT endpoint at y=0 on edge 2 (non-vertex) → count=2, merge.
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();

    Arc a;
    // Region r0: LEFT arc spanning entire LEFT side (edges 0-3)
    a = {}; a.first_edge = 0; a.last_edge = 3; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 4; a.key_y = 0.0; a.key_y_tag = 0;
    std::size_t ai0 = s.add_arc(a);

    // Region r0: RIGHT arc from v4 to P (edges 3-2 on RIGHT)
    a = {}; a.first_edge = 3; a.last_edge = 2; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r0; a.edge_count = 2; a.key_y = 15.0; a.key_y_tag = 4;
    std::size_t ai1 = s.add_arc(a);

    // Region r1: RIGHT arc from P to v0 (edges 2-0 on RIGHT)
    a = {}; a.first_edge = 2; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r1; a.edge_count = 3; a.key_y = 0.0; a.key_y_tag = 100;
    std::size_t ai2 = s.add_arc(a);

    // Chord at y=0, tag=0.
    // LEFT endpoint on edge 0 → vertex 0 (match), count=1.
    // RIGHT endpoint on edge 2 → v2(y=-5,tag=2) and v3(y=5,tag=3),
    //   neither matches (0.0, 0), so non-vertex, count=2.
    Chord c;
    c.region[0] = r0; c.region[1] = r1;
    c.left_adj = {{ai0}, 1};
    c.right_adj = {{ai1, ai2}, 2};
    c.left_edge = 0; c.left_side = LEFT;
    c.right_edge = 2; c.right_side = RIGHT;
    c.y = 0.0; c.y_tag = 0;
    s.add_chord(c);

    s.start_arc = ai0;
    s.end_arc = ai0;
    s.start_vertex = 0;
    s.end_vertex = 4;

    auto poly = non_monotone_polygon();
    s.remove_chord(0, poly);

    // LEFT: vertex endpoint → no merge. a0 survives (ec=4).
    // RIGHT: non-vertex → merge a1+a2. shared_edge=2, nonnull=1.
    //   merged_count = 2 + 3 - 1 = 4.
    assert(s.num_live_arcs() == 2 && "3 adj_arcs: 2 live arcs remain");
    bool found_left = false, found_merged = false;
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        if (s.arc(i).first_side == LEFT && s.arc(i).edge_count == 4)
            found_left = true;
        if (s.arc(i).first_side == RIGHT && s.arc(i).edge_count == 4)
            found_merged = true;
    }
    assert(found_left && "LEFT arc should survive unmerged (ec=4)");
    assert(found_merged && "RIGHT arcs should merge (2+3-1=4)");

    s.assert_tree_property();

    std::printf("  [PASS] remove_chord_3_adj_arcs\n");
}

// ════════════════════════════════════════════════════════════════
//  9. Null-length chord
// ════════════════════════════════════════════════════════════════

static void test_null_length_chord() {
    // [C91 §2.1]: "one of these pairs... gives rise to a chord of
    // null length."  NLCs create empty regions.
    Submap s;
    std::size_t r0 = s.add_node(); // main region
    std::size_t r1 = s.add_node(); // empty NLC region

    // LEFT arc before NLC.
    Arc a;
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 1;
    std::size_t ai0 = s.add_arc(a);

    // Zero-length arc in the NLC empty region.
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 0;
    std::size_t ai1 = s.add_arc(a);

    // LEFT arc after NLC.
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 1;
    std::size_t ai2 = s.add_arc(a);

    // RIGHT arc.
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r0; a.edge_count = 1;
    s.add_arc(a);

    // NLC chord.
    Chord c;
    c.region[0] = r0; c.region[1] = r1;
    c.left_adj = {{ai0}, 1}; c.right_adj = {{ai1}, 1};
    c.left_edge = 1; c.right_edge = 1;
    c.left_side = LEFT; c.right_side = LEFT;
    c.is_null_length = true;
    c.y = 1.0; c.y_tag = 1;
    s.add_chord(c);

    assert(s.chord(0).is_null_length);
    assert(s.node(r1).degree() == 1); // pendant

    // NLC region is empty (weight 0).
    assert(s.region_weight(r1) == 0);

    // Tree property holds.
    s.assert_tree_property();

    // Symbolic y accessor.
    SymbolicY sy = s.chord(0).symbolic_y();
    assert(sy.y == 1.0 && sy.tag == 1);

    std::printf("  [PASS] null_length_chord\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::printf("§2.2 tests:\n");
    test_count_nonnull_edges();
    test_submap_construction();
    test_check_invariants();
    test_region_weight();
    test_remove_chord();
    test_remove_all_chords();
    test_empty_region_weight();
    test_remove_chord_4_adj_arcs();
    test_remove_chord_no_merge_at_vertex();
    test_remove_chord_2_adj_arcs();
    test_remove_chord_3_adj_arcs();
    test_null_length_chord();
    std::printf("All §2.2 tests passed.\n");
    return 0;
}
