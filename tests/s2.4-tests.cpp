// tests/2.4-tests.cpp — Tests for [C91 §2.4]: representation, double identification.

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

// The chordless submap over C = vertices [0, end_vertex]: one region
// bounded by the single closed arc — all of ∂C, one arc-structure
// stored cut at C's start turnaround ([C91 §2.4 tex 142/138]).
static Submap make_chordless(const Polygon& poly, std::size_t end_vertex) {
    Submap s;
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = end_vertex;
    Arc a;
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2 * poly.count_nonnull_edges(0, end_vertex - 1);
    std::size_t ai = s.add_arc(a);
    assert(s.start_arc == ai && s.end_arc == ai &&
           "[C91 §2.4(iii) tex 138]: the closed arc is both endpoint arcs");
    return s;
}

// ════════════════════════════════════════════════════════════════
//  1. double_identify — single closed arc, no ambiguity
// ════════════════════════════════════════════════════════════════

static void test_double_identify_simple() {
    auto poly = test_polygon();
    Submap s = make_chordless(poly, 1);

    // [C91 §2.4 tex 144]: "at most two distinct arcs" pass through a
    // non-chord-endpoint point — here the closed arc covers the point
    // on BOTH ∂C sides but is ONE arc-structure ([C91 §2.4 tex 142]),
    // reported once.
    auto result = s.double_identify(0, {0.0, 0}, poly);
    assert(result.count == 1);

    std::printf("  [PASS] double_identify_simple\n");
}

// ════════════════════════════════════════════════════════════════
//  2. double_identify — multi-edge closed arc
// ════════════════════════════════════════════════════════════════

static void test_double_identify_multi_edge() {
    auto poly = test_polygon();
    Submap s = make_chordless(poly, 4);

    // Interior point: the closed arc covers it on both sides — one
    // structure ([C91 §2.4 tex 142]).
    auto r = s.double_identify(2, {1.5, 0}, poly);
    assert(r.count == 1);

    r = s.double_identify(0, {0.0, 0}, poly);
    assert(r.count == 1);

    r = s.double_identify(5, {0.0, 0}, poly);
    assert(r.count == 0);

    std::printf("  [PASS] double_identify_multi_edge\n");
}

// ════════════════════════════════════════════════════════════════
//  3. double_identify — null-length chord at a y-extremum vertex
// ════════════════════════════════════════════════════════════════

static void test_double_identify_same_edge() {
    // [C91 §2.1 tex 72]: the inside pair at the extremum vertex 1 (a
    // local maximum of the zigzag) carries a null-length chord; its
    // symbolic y is the vertex's own (y, SoS index).  The outer region's
    // single arc W runs from the pair the long way around ∂C — through
    // BOTH endpoint turnarounds — back to the pair: ONE double-wrap
    // structure ([C91 §2.4 tex 142]), first=(1,LEFT) > last=(0,LEFT).
    auto poly = test_polygon();

    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 2;

    Arc a;
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 0;
    std::size_t N = s.add_arc(a);

    a = {}; a.first_edge = 1; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 2 * poly.count_nonnull_edges(0, 1);
    std::size_t W = s.add_arc(a);

    Chord c;
    c = {}; c.region[0] = r0; c.region[1] = r1;
    c.left_edge = 1; c.right_edge = 1; c.left_side = LEFT; c.right_side = LEFT;
    c.y = poly.vertex(1).y; c.y_tag = 1; c.is_null_length = true;
    c.left_adj = {{W}, 1}; c.right_adj = {{N}, 1};
    s.add_chord(c);

    assert(s.start_arc == W && s.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

    // At the null position, exactly the two structures pass through:
    // the inner null arc and W (which both starts and ends there).
    auto r = s.double_identify(1, {poly.vertex(1).y, 1}, poly);
    assert(r.count == 2);

    std::printf("  [PASS] double_identify_same_edge\n");
}

// ════════════════════════════════════════════════════════════════
//  4. Endpoint pointers
// ════════════════════════════════════════════════════════════════

static void test_endpoint_pointers() {
    auto poly = test_polygon();
    Submap s = make_chordless(poly, 2);

    // [C91 §2.4(iii) tex 138]: one arc passes through each of C's
    // turnarounds; for the chordless submap it is the same closed arc.
    assert(s.start_arc == 0);
    assert(s.end_arc == 0);
    assert(s.start_vertex == 0);
    assert(s.end_vertex == 2);

    s.check_invariants();

    std::printf("  [PASS] endpoint_pointers\n");
}

// ════════════════════════════════════════════════════════════════
//  5. Arc-sequence ∂C ordering (LEFT-starting before RIGHT-starting)
// ════════════════════════════════════════════════════════════════

static void test_arc_sequence_ordering() {
    // Chain r0 — c0 — r1 — c1 — r2 over C = vertices 0..4 with both
    // wrap arcs stored as single double-backing structures
    // ([C91 §2.4 tex 142]): table order LEFT-starting ascending
    // first_edge, then RIGHT-starting descending, start-turn arc last
    // ([C91 §2.4(iii) tex 138]).
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();
    s.start_vertex = 0; s.end_vertex = 4;

    Arc a;
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 1;
    std::size_t a1 = s.add_arc(a);
    a = {}; a.first_edge = 2; a.last_edge = 2; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = r2; a.edge_count = 2;          // end wrap: ᾱ = [2,3]
    std::size_t aE = s.add_arc(a);
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r1; a.edge_count = 1;
    std::size_t a4 = s.add_arc(a);
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 1;          // start wrap: ᾱ = [0,0]
    std::size_t aS = s.add_arc(a);

    Chord c;
    c = {}; c.region[0] = r0; c.region[1] = r1;
    c.left_adj = {{aS}, 1}; c.right_adj = {{a4}, 1};
    c.left_edge = 0; c.right_edge = 0; c.y = 3.0; c.y_tag = 1;
    s.add_chord(c);
    c = {}; c.region[0] = r1; c.region[1] = r2;
    c.left_adj = {{a1}, 1}; c.right_adj = {{aE}, 1};
    c.left_edge = 1; c.right_edge = 1; c.y = 1.0; c.y_tag = 2;
    s.add_chord(c);

    assert(s.start_arc == aS && s.end_arc == aE);
    s.check_invariants();

    std::printf("  [PASS] arc_sequence_ordering\n");
}

// ════════════════════════════════════════════════════════════════
//  6. double_identify — the [C91 §2.4 tex 144] extremum worst case
// ════════════════════════════════════════════════════════════════

static void test_double_identify_extremum() {
    // [C91 §2.4 tex 144]: "there are at most six of them, two of which
    // are of zero length: this worst case occurs when q coincides with
    // a vertex of C that is a local extremum."  [C91 §2.1 tex 72]: the
    // extremum vertex 1 has an inside pair (null-length chord cN + null
    // arc N) and an outside pair whose two companions are joined by the
    // through-infinity chord cO with the zero arc Z between them.  On
    // this 2-edge C the four remaining arc slots collapse into the two
    // wrap arcs ([C91 §2.4 tex 142]): arc2 (end wrap) serves as both
    // after-inside and before-outside, arc1 (start wrap) as both
    // after-outside and before-inside.
    auto poly = test_polygon();

    Submap s;
    std::size_t r0 = s.add_node();  // main region
    std::size_t r1 = s.add_node();  // inside the null-length chord
    std::size_t r2 = s.add_node();  // polar cap beyond cO

    s.start_vertex = 0;
    s.end_vertex = 2;

    Arc a;
    // N: inside pair's null arc (LEFT side of vertex 1).
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 0;
    std::size_t N = s.add_arc(a);

    // arc2: from the inside pair down LEFT edge 1, around C's end
    // vertex, up RIGHT edge 1 to the first outside companion.
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = r0; a.edge_count = 2 * poly.count_nonnull_edges(1, 1);
    std::size_t arc2 = s.add_arc(a);

    // Z: the outside pair's zero arc (RIGHT side of vertex 1).
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r2; a.edge_count = 0;
    std::size_t Z = s.add_arc(a);

    // arc1: from the second outside companion down RIGHT edge 0, around
    // C's start vertex, up LEFT edge 0 to the inside pair.
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 2 * poly.count_nonnull_edges(0, 0);
    std::size_t arc1 = s.add_arc(a);

    // cN: the inside pair's null-length chord — symbolic y is the
    // vertex's own (y, SoS index) ([C91 §2 tex 47]).
    Chord c;
    c = {}; c.region[0] = r0; c.region[1] = r1;
    c.left_edge = 1; c.right_edge = 1; c.left_side = LEFT; c.right_side = LEFT;
    c.y = poly.vertex(1).y; c.y_tag = 1; c.is_null_length = true;
    c.left_adj = {{arc1}, 1}; c.right_adj = {{N}, 1};
    s.add_chord(c);

    // cO: the outside pair's chord — runs through infinity from one
    // companion to the other ([C91 §2.1 tex 70/72]); vertex endpoints,
    // one adj arc each ([C91 §2.2 tex 94]).
    c = {}; c.region[0] = r0; c.region[1] = r2;
    c.left_edge = 1; c.left_side = RIGHT; c.left_adj = {{arc2}, 1};
    c.right_edge = 0; c.right_side = RIGHT; c.right_adj = {{Z}, 1};
    c.y = poly.vertex(1).y; c.y_tag = 1;
    s.add_chord(c);

    assert(s.start_arc == arc1 && s.end_arc == arc2);
    s.check_invariants(poly);

    // Query at the extremum on edge 1: the structures through it are
    // arc2 and N; on edge 0: arc1 and Z.  All four structures at the
    // vertex, two of zero length ([C91 §2.4 tex 144]).
    auto rl = s.double_identify(1, {poly.vertex(1).y, 1}, poly);
    assert(rl.count == 2);
    auto rr = s.double_identify(0, {poly.vertex(1).y, 1}, poly);
    assert(rr.count == 2);
    assert(rl.count + rr.count <= Submap::DoubleIdentifyResult::MAX);

    std::printf("  [PASS] double_identify_extremum\n");
}

// ════════════════════════════════════════════════════════════════
//  7. double_identify — edge not in any arc
// ════════════════════════════════════════════════════════════════

static void test_double_identify_miss() {
    auto poly = test_polygon();
    Submap s = make_chordless(poly, 1);

    auto r = s.double_identify(5, {0.0, 0}, poly);
    assert(r.count == 0);

    std::printf("  [PASS] double_identify_miss\n");
}

// ════════════════════════════════════════════════════════════════
//  8. check_invariants(polygon) — direct positive test
// ════════════════════════════════════════════════════════════════

static void test_check_invariants_polygon_positive() {
    // [C91 §2.4(i)–(iii) + §2.2]: Polygon-aware overload accepts a
    // well-formed chordless submap (single closed arc, edge_count cache
    // consistent with the polygon).
    auto poly = test_polygon();
    Submap s = make_chordless(poly, 4);
    s.check_invariants(poly);

    std::printf("  [PASS] check_invariants_polygon_positive\n");
}

// ════════════════════════════════════════════════════════════════
//  9. check_invariants(polygon) — size-3 same-first_edge run
//     ([C91 §2.4 tex 144] monotonicity beyond trivial size-2)
// ════════════════════════════════════════════════════════════════

// Polygon for the nested-pocket run fixtures: edge 1 (y 10→1) carries
// three mid-edge chord endpoints at the vertex ys 8, 4, 2.
static Polygon run_polygon() {
    return Polygon({
        {0,0,0}, {1,10,1}, {2,1,2}, {3,8,3}, {4,4,4}, {5,2,5}
    });
}

// Nested pockets c1 ⊃ c2 ⊃ c3 between edges 1 and 2 ([C91 §2.2 tex
// 102]: proper parenthesis nesting).  Edge 1's LEFT face carries a
// 3-arc same-first_edge run (P1, P2, P3 with start-ys 8, 4, 2 derived
// from their bounding chords, [C91 §2.4 tex 133]); r0's arc W runs from
// c1's far (vertex-3) endpoint around BOTH of C's turnarounds — one
// double-wrap structure ([C91 §2.4 tex 142]).
//
// @param swap_p2_p3  insert P3 before P2 (breaks [C91 §2.4 tex 144]
//                    run monotonicity — for the death test).
static Submap build_run_submap(const Polygon& poly, bool swap_p2_p3) {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();
    std::size_t r3 = s.add_node();
    s.start_vertex = 0; s.end_vertex = 5;

    Arc a;
    auto mk = [&](std::size_t fe, Side fs, std::size_t le, Side ls,
                  std::size_t region, std::size_t ec) {
        a = {}; a.first_edge = fe; a.first_side = fs;
        a.last_edge = le; a.last_side = ls;
        a.region_node = region; a.edge_count = ec;
        return s.add_arc(a);
    };

    std::size_t P1 = mk(1, LEFT, 1, LEFT, r1, 1);
    std::size_t P2 = NONE, P3 = NONE;
    if (swap_p2_p3) {
        P3 = mk(1, LEFT, 2, LEFT, r3, 2);
        P2 = mk(1, LEFT, 1, LEFT, r2, 1);
    } else {
        P2 = mk(1, LEFT, 1, LEFT, r2, 1);
        P3 = mk(1, LEFT, 2, LEFT, r3, 2);
    }
    std::size_t P2b = mk(2, LEFT, 2, LEFT, r2, 1);
    std::size_t P1b = mk(2, LEFT, 2, LEFT, r1, 1);
    // W: double-wrap ([C91 §2.4 tex 142]) — starts at vertex 3 leaving
    // on edge 3, covers LEFT [3,4] + all of RIGHT + LEFT [0, part of 1];
    // ᾱ = all of C.  [C91 §2.2 tex 106]: nonnull ∂C edges per leg =
    // 2 + 5 + 2.
    std::size_t W = mk(3, LEFT, 1, LEFT, r0,
                       poly.count_nonnull_edges(3, 4) +
                       poly.count_nonnull_edges(0, 4) +
                       poly.count_nonnull_edges(0, 1));

    Chord c;
    // c1: y=8 (vertex 3's y).  LEFT endpoint mid-edge-1 → {W, P1};
    // RIGHT endpoint AT vertex 3 → 1 slot {P1b} ([C91 §2.2 tex 94]).
    c = {}; c.region[0] = r0; c.region[1] = r1;
    c.left_edge = 1; c.left_side = LEFT; c.left_adj = {{W, P1}, 2};
    c.right_edge = 2; c.right_side = LEFT; c.right_adj = {{P1b}, 1};
    c.y = 8.0; c.y_tag = 3;
    s.add_chord(c);
    // c2: y=4 (vertex 4's y), both endpoints mid-edge.
    c = {}; c.region[0] = r1; c.region[1] = r2;
    c.left_edge = 1; c.left_side = LEFT; c.left_adj = {{P1, P2}, 2};
    c.right_edge = 2; c.right_side = LEFT; c.right_adj = {{P2b, P1b}, 2};
    c.y = 4.0; c.y_tag = 4;
    s.add_chord(c);
    // c3: y=2 (vertex 5's y), both endpoints mid-edge.
    c = {}; c.region[0] = r2; c.region[1] = r3;
    c.left_edge = 1; c.left_side = LEFT; c.left_adj = {{P2, P3}, 2};
    c.right_edge = 2; c.right_side = LEFT; c.right_adj = {{P3, P2b}, 2};
    c.y = 2.0; c.y_tag = 5;
    s.add_chord(c);

    assert(s.start_arc == W && s.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");
    return s;
}

static void test_check_invariants_polygon_size_3_run() {
    // [C91 §2.4 tex 144]: arcs sharing first_edge must be start-y
    // monotonic.  P1(8), P2(4), P3(2) on descending edge 1 —
    // check_invariants(polygon) infers descending and verifies all
    // consecutive pairs.
    auto poly = run_polygon();
    Submap s = build_run_submap(poly, /*swap_p2_p3=*/false);
    s.check_invariants(poly);

    std::printf("  [PASS] check_invariants_polygon_size_3_run\n");
}

// ════════════════════════════════════════════════════════════════
//  10. check_invariants(polygon) — wrapped (LEFT→RIGHT) endpoint arc
//      ([C91 §2.4 tex 142] double-backing)
// ════════════════════════════════════════════════════════════════

static void test_check_invariants_polygon_wrapped_arc() {
    // [C91 §2.4 tex 142]: "an arc might wrap around both sides of C,
    // something we call double-backing."  Single-arc submap whose only
    // arc wraps at c_end (LEFT→RIGHT); the edge_count cache check
    // ([C91 §2.2 tex 106]) counts nonnull ∂C edges per leg.
    Polygon poly({{0,0,0}, {1,1,1}, {2,0,2}});

    Submap s;
    s.add_node();
    s.start_vertex = 0; s.end_vertex = 2;

    Arc a;
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2 * poly.count_nonnull_edges(0, 1);
    std::size_t ai0 = s.add_arc(a);
    assert(s.start_arc == ai0 && s.end_arc == ai0);

    s.check_invariants(poly);

    std::printf("  [PASS] check_invariants_polygon_wrapped_arc\n");
}

// ════════════════════════════════════════════════════════════════
//  11. Death: non-monotonic start-y in a same-first_edge run
//      ([C91 §2.4 tex 144])
// ════════════════════════════════════════════════════════════════

static void test_non_monotonic_run_fires() {
    // [C91 §2.4 tex 144]: same-first_edge run must be start-y monotonic
    // in the inferred direction.  With P3 inserted before P2 the run
    // reads 8 → 2 → 4 — descending then ascending — which the
    // run-direction inference flags as non-monotonic.
    assert(assert_fires([]{
        auto poly = run_polygon();
        Submap s = build_run_submap(poly, /*swap_p2_p3=*/true);
        s.check_invariants(poly);  // ← non-monotonic must fire
    }));
    std::printf("  [PASS] non_monotonic_run_fires\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::printf("[C91 §2.4 tests]:\n");
    test_double_identify_simple();
    test_double_identify_multi_edge();
    test_double_identify_same_edge();
    test_endpoint_pointers();
    test_arc_sequence_ordering();
    test_double_identify_extremum();
    test_double_identify_miss();
    test_check_invariants_polygon_positive();
    test_check_invariants_polygon_size_3_run();
    test_check_invariants_polygon_wrapped_arc();
    test_non_monotonic_run_fires();
    std::printf("All §2.4 tests passed.\n");
    return 0;
}
