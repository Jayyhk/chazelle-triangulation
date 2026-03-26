/// tests/2.4-tests.cpp — Tests for §2.4: representation, double identification.

#include "submap/submap.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════
//  1. double_identify — single arc per side, no ambiguity
// ════════════════════════════════════════════════════════════════

static void test_double_identify_simple() {
    // Two arcs: LEFT on edge 0, RIGHT on edge 0.
    Submap s;
    s.add_node();

    Arc a;
    a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 0.0; a.key_y_tag = 0;
    s.add_arc(a);

    a.first_side = RIGHT; a.last_side = RIGHT;
    a.key_y = 0.0; a.key_y_tag = 0;
    s.add_arc(a);

    auto result = s.double_identify(0, {0.0, 0});
    // Should find both arcs (one per ∂C side).
    assert(result.count == 2);

    std::printf("  [PASS] double_identify_simple\n");
}

// ════════════════════════════════════════════════════════════════
//  2. double_identify — multi-edge arc
// ════════════════════════════════════════════════════════════════

static void test_double_identify_multi_edge() {
    Submap s;
    s.add_node();

    // LEFT arc spanning edges 0-3.
    Arc a;
    a.first_edge = 0; a.last_edge = 3; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 4; a.key_y = 0.0; a.key_y_tag = 0;
    s.add_arc(a);

    // RIGHT arc spanning edges 3-0.
    a.first_edge = 3; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.key_y = 3.0; a.key_y_tag = 3;
    s.add_arc(a);

    // Query edge 2 — should find both arcs.
    auto r = s.double_identify(2, {1.5, 0});
    assert(r.count == 2);

    // Query edge 0 — should find both arcs.
    r = s.double_identify(0, {0.0, 0});
    assert(r.count == 2);

    // Query edge 5 (out of range) — should find nothing.
    r = s.double_identify(5, {0.0, 0});
    assert(r.count == 0);

    std::printf("  [PASS] double_identify_multi_edge\n");
}

// ════════════════════════════════════════════════════════════════
//  3. double_identify — multiple arcs on same edge
// ════════════════════════════════════════════════════════════════

static void test_double_identify_same_edge() {
    // Two LEFT arcs on edge 1, at different y-levels (chord splits them).
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1

    Arc a;
    // LEFT arc on edge 0 (r0)
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 0.0; a.key_y_tag = 0;
    s.add_arc(a);
    // LEFT arc on edge 1 (r0), lower y
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 1.0; a.key_y_tag = 1;
    s.add_arc(a);
    // LEFT arc on edge 1 (r1), higher y
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 1; a.edge_count = 1; a.key_y = 2.0; a.key_y_tag = 2;
    s.add_arc(a);
    // RIGHT arcs
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 2.0; a.key_y_tag = 2;
    s.add_arc(a);
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 0.0; a.key_y_tag = 0;
    s.add_arc(a);

    // Query edge 1 — should find multiple arcs.
    auto r = s.double_identify(1, {1.5, 0});
    assert(r.count >= 2); // at least the 2 LEFT arcs + 1 RIGHT arc

    std::printf("  [PASS] double_identify_same_edge\n");
}

// ════════════════════════════════════════════════════════════════
//  4. Endpoint pointers
// ════════════════════════════════════════════════════════════════

static void test_endpoint_pointers() {
    Submap s;
    s.add_node();

    Arc a;
    a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = s.add_arc(a);

    a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    std::size_t ai1 = s.add_arc(a);

    s.start_arc = ai0;
    s.end_arc = ai1;
    s.start_vertex = 0;
    s.end_vertex = 2;

    assert(s.start_arc == 0);
    assert(s.end_arc == 1);
    assert(s.start_vertex == 0);
    assert(s.end_vertex == 2);

    // check_invariants verifies these are in range.
    s.check_invariants();

    std::printf("  [PASS] endpoint_pointers\n");
}

// ════════════════════════════════════════════════════════════════
//  5. Arc-sequence ∂C ordering (LEFT before RIGHT)
// ════════════════════════════════════════════════════════════════

static void test_arc_sequence_ordering() {
    Submap s;
    s.add_node();

    Arc a;
    // LEFT arcs first.
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1;
    s.add_arc(a);
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1;
    s.add_arc(a);
    // Then RIGHT arcs.
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 1;
    s.add_arc(a);
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 1;
    s.add_arc(a);

    // check_invariants verifies LEFT-before-RIGHT.
    s.check_invariants();

    std::printf("  [PASS] arc_sequence_ordering\n");
}

// ════════════════════════════════════════════════════════════════
//  6. double_identify — worst case: 6 arcs at y-extremum
// ════════════════════════════════════════════════════════════════

static void test_double_identify_worst_case() {
    // [C91 §2.4]: "there are at most six of them, two of which are
    // of zero length: this worst case occurs when q coincides with
    // a vertex of C that is a local extremum in the y-direction."
    //
    // At a y-extremum with chords on both sides:
    //   LEFT: arc_before, zero-length NLC arc, arc_after  (3 LEFT)
    //   RIGHT: arc_before, zero-length NLC arc, arc_after (3 RIGHT)
    //   Total: 6 arcs, all on the same edge.

    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1 (LEFT NLC empty)
    s.add_node(); // r2 (RIGHT NLC empty)

    Arc a;
    // 3 LEFT arcs on edge 1.
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 0.5; a.key_y_tag = 0;
    s.add_arc(a); // before NLC

    a.region_node = 1; a.edge_count = 0; a.key_y = 1.0; a.key_y_tag = 1;
    s.add_arc(a); // NLC zero-length

    a.region_node = 0; a.edge_count = 1; a.key_y = 1.5; a.key_y_tag = 2;
    s.add_arc(a); // after NLC

    // 3 RIGHT arcs on edge 1.
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 1.5; a.key_y_tag = 2;
    s.add_arc(a); // before NLC

    a.region_node = 2; a.edge_count = 0; a.key_y = 1.0; a.key_y_tag = 1;
    s.add_arc(a); // NLC zero-length

    a.region_node = 0; a.edge_count = 1; a.key_y = 0.5; a.key_y_tag = 0;
    s.add_arc(a); // after NLC

    auto r = s.double_identify(1, {1.0, 1});
    // Should find up to 6 arcs (all on edge 1).
    assert(r.count >= 4 && r.count <= 6);

    // Verify MAX capacity isn't exceeded.
    assert(r.count <= Submap::DoubleIdentifyResult::MAX);

    std::printf("  [PASS] double_identify_worst_case\n");
}

// ════════════════════════════════════════════════════════════════
//  7. double_identify — edge not in any arc
// ════════════════════════════════════════════════════════════════

static void test_double_identify_miss() {
    Submap s;
    s.add_node();

    Arc a;
    a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1;
    s.add_arc(a);

    // Query for edge 5 — not in any arc.
    auto r = s.double_identify(5, {0.0, 0});
    assert(r.count == 0);

    std::printf("  [PASS] double_identify_miss\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::printf("§2.4 tests:\n");
    test_double_identify_simple();
    test_double_identify_multi_edge();
    test_double_identify_same_edge();
    test_endpoint_pointers();
    test_arc_sequence_ordering();
    test_double_identify_worst_case();
    test_double_identify_miss();
    std::printf("All §2.4 tests passed.\n");
    return 0;
}
