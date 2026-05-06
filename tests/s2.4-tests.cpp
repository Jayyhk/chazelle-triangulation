/// tests/2.4-tests.cpp — Tests for §2.4: representation, double identification.

#include "submap/submap.h"
#include "polygon/polygon.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <functional>
#include <sys/prctl.h>
#include <sys/wait.h>
#include <unistd.h>

namespace {
// Death-test helper: forks; child runs `fn` with stderr silenced and core
// dumps disabled (PR_SET_DUMPABLE=0 skips the kernel's coredump pipe — on
// Ubuntu/apport systems this avoids ~1s of crash-reporter latency per
// abort).  Returns true iff the child terminated by SIGABRT.
bool assert_fires(std::function<void()> fn) {
    pid_t pid = fork();
    if (pid == 0) {
        prctl(PR_SET_DUMPABLE, 0);
        if (freopen("/dev/null", "w", stderr) == nullptr) std::_Exit(2);
        fn();
        std::_Exit(0);
    }
    int status = 0;
    if (waitpid(pid, &status, 0) < 0) return false;
    return WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
}
}

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
//  8. check_invariants(polygon) — direct positive test
// ════════════════════════════════════════════════════════════════

static void test_check_invariants_polygon_positive() {
    // [C91 §2.4 (i)–(iii) + §2.2]: Polygon-aware overload accepts a
    // well-formed normal-form submap (key_y monotonic + edge_count
    // cache consistent with polygon).
    auto poly = test_polygon();
    Submap s;
    s.add_node();

    Arc a;
    a = {}; a.first_edge = 0; a.last_edge = 3; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0;
    a.edge_count = poly.count_nonnull_edges(0, 3);
    a.key_y = poly.vertex(0).y; a.key_y_tag = 0;
    std::size_t ai0 = s.add_arc(a);

    a = {}; a.first_edge = 3; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = poly.count_nonnull_edges(0, 3);
    a.key_y = poly.vertex(4).y; a.key_y_tag = 4;
    s.add_arc(a);

    s.start_arc = ai0; s.end_arc = ai0;
    s.start_vertex = 0; s.end_vertex = 4;

    s.check_invariants(poly);

    std::printf("  [PASS] check_invariants_polygon_positive\n");
}

// ════════════════════════════════════════════════════════════════
//  9. check_invariants(polygon) — size-3 same-first_edge run
//     (§2.4 tex 144 monotonicity beyond trivial size-2)
// ════════════════════════════════════════════════════════════════

static void test_check_invariants_polygon_size_3_run() {
    // [C91 §2.4 tex 144]: arcs sharing first_edge must be key_y
    // monotonic.  Build a 3-arc run on edge 1 (LEFT side) separated
    // by 2 NLCs; SoS tags 2/5/6 give descending perturbed y in array
    // order.  check_invariants(polygon) infers descending and verifies
    // all consecutive pairs.
    auto poly = test_polygon();

    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();
    double yv = poly.vertex(1).y;

    Arc a;
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = poly.count_nonnull_edges(1, 1);
    a.key_y = yv; a.key_y_tag = 2;
    std::size_t ai0 = s.add_arc(a);

    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 0;
    a.key_y = yv; a.key_y_tag = 5;
    std::size_t ai1 = s.add_arc(a);

    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r2; a.edge_count = 0;
    a.key_y = yv; a.key_y_tag = 6;
    std::size_t ai2 = s.add_arc(a);

    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r0; a.edge_count = poly.count_nonnull_edges(1, 1);
    a.key_y = yv; a.key_y_tag = 2;
    s.add_arc(a);

    Chord c;
    c = {}; c.region[0] = r0; c.region[1] = r1;
    c.left_edge = 1; c.right_edge = 1; c.left_side = LEFT; c.right_side = LEFT;
    c.y = yv; c.y_tag = 5; c.is_null_length = true;
    c.left_adj = {{ai0}, 1}; c.right_adj = {{ai1}, 1};
    s.add_chord(c);

    c = {}; c.region[0] = r1; c.region[1] = r2;
    c.left_edge = 1; c.right_edge = 1; c.left_side = LEFT; c.right_side = LEFT;
    c.y = yv; c.y_tag = 6; c.is_null_length = true;
    c.left_adj = {{ai1}, 1}; c.right_adj = {{ai2}, 1};
    s.add_chord(c);

    s.start_arc = ai0; s.end_arc = ai2;
    s.start_vertex = 1; s.end_vertex = 2;

    s.check_invariants(poly);

    std::printf("  [PASS] check_invariants_polygon_size_3_run\n");
}

// ════════════════════════════════════════════════════════════════
//  10. check_invariants(polygon) — wrapped (LEFT→RIGHT) endpoint arc
//      (§2.4 tex 142 double-backing)
// ════════════════════════════════════════════════════════════════

static void test_check_invariants_polygon_wrapped_arc() {
    // [C91 §2.4 tex 142]: "an arc might wrap around both sides of C,
    // something we call double-backing."  Single-arc submap whose only
    // arc wraps at c_end (LEFT→RIGHT); the edge_count cache check
    // (§2.2 tex 106) uses the union-range count.
    Polygon poly({{0,0,0}, {1,1,1}, {2,0,2}});

    Submap s;
    s.add_node();

    Arc a;
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = poly.count_nonnull_edges(0, 1);
    a.key_y = poly.vertex(0).y; a.key_y_tag = 0;
    std::size_t ai0 = s.add_arc(a);

    s.start_arc = ai0; s.end_arc = ai0;
    s.start_vertex = 0; s.end_vertex = 2;

    s.check_invariants(poly);

    std::printf("  [PASS] check_invariants_polygon_wrapped_arc\n");
}

// ════════════════════════════════════════════════════════════════
//  11. Death: non-monotonic key_y in a same-first_edge run
//      (§2.4 tex 144)
// ════════════════════════════════════════════════════════════════

static void test_non_monotonic_run_fires() {
    // [C91 §2.4 tex 144]: same-first_edge run must be key_y monotonic
    // in the inferred direction.  Endpoints span tags 2..6 (descending
    // inferred) but middle arc has tag 1 (HIGHER perturbed y than
    // first) — non-monotonic.
    assert(assert_fires([]{
        auto poly = test_polygon();
        Submap s;
        std::size_t r0 = s.add_node();
        std::size_t r1 = s.add_node();
        std::size_t r2 = s.add_node();
        double yv = poly.vertex(1).y;

        Arc a;
        a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r0; a.edge_count = poly.count_nonnull_edges(1, 1);
        a.key_y = yv; a.key_y_tag = 2;
        std::size_t ai0 = s.add_arc(a);

        a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r1; a.edge_count = 0;
        a.key_y = yv; a.key_y_tag = 1;            // ← OUT OF ORDER
        std::size_t ai1 = s.add_arc(a);

        a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r2; a.edge_count = 0;
        a.key_y = yv; a.key_y_tag = 6;
        std::size_t ai2 = s.add_arc(a);

        a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
        a.region_node = r0; a.edge_count = poly.count_nonnull_edges(1, 1);
        a.key_y = yv; a.key_y_tag = 2;
        s.add_arc(a);

        Chord c;
        c = {}; c.region[0] = r0; c.region[1] = r1;
        c.left_edge = 1; c.right_edge = 1; c.left_side = LEFT; c.right_side = LEFT;
        c.y = yv; c.y_tag = 1; c.is_null_length = true;
        c.left_adj = {{ai0}, 1}; c.right_adj = {{ai1}, 1};
        s.add_chord(c);

        c = {}; c.region[0] = r1; c.region[1] = r2;
        c.left_edge = 1; c.right_edge = 1; c.left_side = LEFT; c.right_side = LEFT;
        c.y = yv; c.y_tag = 6; c.is_null_length = true;
        c.left_adj = {{ai1}, 1}; c.right_adj = {{ai2}, 1};
        s.add_chord(c);

        s.start_arc = ai0; s.end_arc = ai2;
        s.start_vertex = 1; s.end_vertex = 2;

        s.check_invariants(poly);  // ← non-monotonic must fire
    }));
    std::printf("  [PASS] non_monotonic_run_fires\n");
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
    test_check_invariants_polygon_positive();
    test_check_invariants_polygon_size_3_run();
    test_check_invariants_polygon_wrapped_arc();
    test_non_monotonic_run_fires();
    std::printf("All §2.4 tests passed.\n");
    return 0;
}
