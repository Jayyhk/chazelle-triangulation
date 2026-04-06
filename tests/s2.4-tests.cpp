/// tests/2.4-tests.cpp — Tests for §2.4: representation, double identification.

#include "submap/submap.h"
#include "polygon/polygon.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

// Zigzag polygon with both ascending and descending edges.
// Edges: 0(up), 1(down), 2(up), 3(down).  5 vertices, 4 edges.
// Distinct y-coordinates and distinct indices for SoS.
static Polygon test_polygon() {
    return Polygon({
        {0,0,0}, {1,3,1}, {2,1,2}, {3,4,3}, {4,2,4}
    });
}

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
    std::size_t ai0 = s.add_arc(a);

    a.first_side = RIGHT; a.last_side = RIGHT;
    a.key_y = 0.0; a.key_y_tag = 0;
    std::size_t ai1 = s.add_arc(a);

    // §2.4 (tex 144): end_arc = last LEFT arc (the turnaround point).
    s.start_arc = ai0; s.end_arc = ai0;
    s.start_vertex = 0; s.end_vertex = 1;

    auto poly = test_polygon();
    auto result = s.double_identify(0, {0.0, 0}, poly);
    assert(result.count == 2);

    std::printf("  [PASS] double_identify_simple\n");
}

// ════════════════════════════════════════════════════════════════
//  2. double_identify — multi-edge arc
// ════════════════════════════════════════════════════════════════

static void test_double_identify_multi_edge() {
    Submap s;
    s.add_node();

    Arc a;
    a.first_edge = 0; a.last_edge = 3; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 4; a.key_y = 0.0; a.key_y_tag = 0;
    std::size_t ai0 = s.add_arc(a);

    a.first_edge = 3; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.key_y = 3.0; a.key_y_tag = 3;
    std::size_t ai1 = s.add_arc(a);

    // §2.4 (tex 144): end_arc = last LEFT arc (ai0 is the only LEFT arc).
    s.start_arc = ai0; s.end_arc = ai0;
    s.start_vertex = 0; s.end_vertex = 4;

    auto poly = test_polygon();
    auto r = s.double_identify(2, {1.5, 0}, poly);
    assert(r.count == 2);

    r = s.double_identify(0, {0.0, 0}, poly);
    assert(r.count == 2);

    r = s.double_identify(5, {0.0, 0}, poly);
    assert(r.count == 0);

    std::printf("  [PASS] double_identify_multi_edge\n");
}

// ════════════════════════════════════════════════════════════════
//  3. double_identify — multiple arcs on same edge
// ════════════════════════════════════════════════════════════════

static void test_double_identify_same_edge() {
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1

    Arc a;
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 0.0; a.key_y_tag = 0;
    std::size_t ai0 = s.add_arc(a);

    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 1.0; a.key_y_tag = 1;
    s.add_arc(a);

    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 1; a.edge_count = 1; a.key_y = 2.0; a.key_y_tag = 2;
    std::size_t ai_end_left = s.add_arc(a); // last LEFT arc — §2.4 (tex 144)

    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 2.0; a.key_y_tag = 2;
    s.add_arc(a);

    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 0.0; a.key_y_tag = 0;
    s.add_arc(a);

    // §2.4 (tex 144): end_arc = last LEFT arc (ai_end_left).
    s.start_arc = ai0; s.end_arc = ai_end_left;
    s.start_vertex = 0; s.end_vertex = 2;

    auto poly = test_polygon();
    auto r = s.double_identify(1, {1.5, 0}, poly);
    assert(r.count >= 2);

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
    // §2.4 (tex 144): end_arc = last LEFT arc (ai0 is the only LEFT arc).
    s.end_arc = ai0;
    s.start_vertex = 0;
    s.end_vertex = 2;

    assert(s.start_arc == 0);
    assert(s.end_arc == 0);
    assert(s.start_vertex == 0);
    assert(s.end_vertex == 2);

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
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1;
    std::size_t ai0 = s.add_arc(a);
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1;
    std::size_t ai_end = s.add_arc(a);  // contains c_end_edge=1
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 1;
    s.add_arc(a);
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 1;
    s.add_arc(a);

    s.start_arc = ai0; s.end_arc = ai_end;
    s.start_vertex = 0; s.end_vertex = 2;

    s.check_invariants();

    std::printf("  [PASS] arc_sequence_ordering\n");
}

// ════════════════════════════════════════════════════════════════
//  6. double_identify — worst case: 6 arcs at y-extremum
// ════════════════════════════════════════════════════════════════

static void test_double_identify_worst_case() {
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1 (LEFT NLC empty)
    s.add_node(); // r2 (RIGHT NLC empty)

    Arc a;
    // 3 LEFT arcs on edge 1.
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 0.5; a.key_y_tag = 0;
    std::size_t ai0 = s.add_arc(a);

    a.region_node = 1; a.edge_count = 0; a.key_y = 1.0; a.key_y_tag = 1;
    s.add_arc(a);

    a.region_node = 0; a.edge_count = 1; a.key_y = 1.5; a.key_y_tag = 2;
    std::size_t ai_end_left = s.add_arc(a); // last LEFT arc — §2.4 (tex 144)

    // 3 RIGHT arcs on edge 1.
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 1; a.key_y = 1.5; a.key_y_tag = 2;
    s.add_arc(a);

    a.region_node = 2; a.edge_count = 0; a.key_y = 1.0; a.key_y_tag = 1;
    s.add_arc(a);

    a.region_node = 0; a.edge_count = 1; a.key_y = 0.5; a.key_y_tag = 0;
    s.add_arc(a);

    // §2.4 (tex 144): end_arc = last LEFT arc (ai_end_left = arc 2).
    s.start_arc = ai0; s.end_arc = ai_end_left;
    s.start_vertex = 0; s.end_vertex = 2;

    auto poly = test_polygon();
    auto r = s.double_identify(1, {1.0, 1}, poly);
    assert(r.count >= 4 && r.count <= 6);
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
    std::size_t ai0 = s.add_arc(a);

    a.first_side = RIGHT; a.last_side = RIGHT;
    std::size_t ai1 = s.add_arc(a);

    // §2.4 (tex 144): end_arc = last LEFT arc (ai0 is the only LEFT arc).
    s.start_arc = ai0; s.end_arc = ai0;
    s.start_vertex = 0; s.end_vertex = 1;

    auto poly = test_polygon();
    auto r = s.double_identify(5, {0.0, 0}, poly);
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
