// tests/2.2-tests.cpp — Tests for [C91 §2.2]: submaps, weight, chord removal.

#include "polygon/polygon.h"
#include "submap/submap.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <functional>
#include <sys/prctl.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace chazelle;

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
//  [C91 §2.4 tex 142]: a wrap-spanning arc is ONE arc-structure that
//  double-backs around an endpoint of C — never split.  Region 2's arc
//  passes around C's end vertex (end wrap) and region 0's around C's
//  start vertex (start wrap), so the table holds 4 arc-structures:
//    aE  L(2)→R(2)  r2 (end wrap),   a1 L[1,1] r1,
//    a4  R[1,1] r1,                  aS R(0)→L(0) r0 (start wrap).
// ════════════════════════════════════════════════════════════════

static Submap build_3region_submap() {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();

    // [C91 §2.4(iii) tex 138]: set C's endpoints before the arcs so
    // add_arc can classify the wrap arcs.
    s.start_vertex = 0;
    s.end_vertex = 4;

    // Canonical table order ([C91 §2.4(iii) tex 138]): LEFT-starting
    // arcs ascending first_edge, then RIGHT-starting descending; the
    // start-turn arc sorts last.

    Arc a1; a1.first_edge = 1; a1.last_edge = 1; a1.first_side = LEFT; a1.last_side = LEFT;
    a1.region_node = r1; a1.edge_count = 1;
    std::size_t ai1 = s.add_arc(a1);

    // Region 2's arc: from chord 1's LEFT endpoint up edge 2–3, around
    // C's end vertex, back down to chord 1's RIGHT endpoint —
    // ONE end-wrap structure ([C91 §2.4 tex 142]); ᾱ = edges [2,3].
    Arc aE; aE.first_edge = 2; aE.last_edge = 2; aE.first_side = LEFT; aE.last_side = RIGHT;
    aE.region_node = r2; aE.edge_count = 2;
    std::size_t aiE = s.add_arc(aE);

    Arc a4; a4.first_edge = 1; a4.last_edge = 1; a4.first_side = RIGHT; a4.last_side = RIGHT;
    a4.region_node = r1; a4.edge_count = 1;
    std::size_t ai4 = s.add_arc(a4);

    // Region 0's arc: from chord 0's RIGHT endpoint down edge 0, around
    // C's start vertex, back up to chord 0's LEFT endpoint — ONE
    // start-wrap structure; ᾱ = edge [0,0].
    Arc aS; aS.first_edge = 0; aS.last_edge = 0; aS.first_side = RIGHT; aS.last_side = LEFT;
    aS.region_node = r0; aS.edge_count = 1;
    std::size_t aiS = s.add_arc(aS);

    // Both chords sit at polygon vertices (y_tag 1 and 2 match vertices
    // of test_polygon), so each endpoint records ONE adj arc — the
    // before-arc in ∂C traversal order ([C91 §2.2 tex 94 + §2.4(ii)]
    // 1-slot convention; the after-arc is found by remove_chord's
    // junction scan).

    // Chord 0: between r0 and r1 at vertex 1.  LEFT before-arc aS (its
    // LEFT leg ends at v1); RIGHT before-arc a4 (RIGHT traversal
    // descends: [1,1] ends at v1).
    Chord c0;
    c0.region[0] = r0; c0.region[1] = r1;
    c0.left_adj = {{aiS}, 1}; c0.right_adj = {{ai4}, 1};
    c0.left_edge = 0; c0.right_edge = 0; c0.y = 1.0; c0.y_tag = 1;
    s.add_chord(c0);

    // Chord 1: between r1 and r2 at vertex 2.  LEFT before-arc a1;
    // RIGHT before-arc aE (its RIGHT leg [2,3] descending ends at v2).
    Chord c1;
    c1.region[0] = r1; c1.region[1] = r2;
    c1.left_adj = {{ai1}, 1}; c1.right_adj = {{aiE}, 1};
    c1.left_edge = 1; c1.right_edge = 1; c1.y = 2.0; c1.y_tag = 2;
    s.add_chord(c1);

    // [C91 §2.4(iii) tex 138]: add_arc auto-registered the endpoint
    // pointers from the wrap classes.
    assert(s.start_arc == aiS && s.end_arc == aiE);

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
    // [C91 §2.2 tex 96 + §2.4 tex 142]: arc-structure count == number of
    // arcs — one per boundary alternation slot, wrap-spanning arcs
    // included as single structures.
    assert(s.num_arcs() == 4);

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

    // r0's single start-wrap arc has ec=1 → weight = 1
    assert(s.region_weight(0) == 1);
    // r1 has arcs a1(ec=1) and a4(ec=1) → weight = 1
    assert(s.region_weight(1) == 1);
    // r2's single end-wrap arc spans edges [2,3] (ec=2) → weight = 2
    // ([C91 §2.2 tex 106]: single count over ᾱ, [C91 §2.4 tex 142]).
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

static void test_chordless_region_weight() {
    // [C91 §2.2 tex 106 + §2.4 tex 142]: the chordless submap is one
    // region bounded by the single closed arc (all of ∂C), stored cut
    // at C's start turnaround; its weight is the arc's nonnull count.
    // ([C91 §2.2]'s "weight 0 if the region is empty" is exercised by
    // the null-length chord tests below.)
    auto poly = test_polygon();
    Submap s;
    std::size_t r0 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a;
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;   // closed-arc encoding
    a.region_node = r0;
    a.edge_count = poly.count_nonnull_edges(0, 3);
    std::size_t ai = s.add_arc(a);
    assert(s.start_arc == ai && s.end_arc == ai &&
           "[C91 §2.4(iii) tex 138]: the closed arc is both endpoint arcs");

    assert(s.region_weight(r0) == 4);
    s.check_invariants(poly);

    std::printf("  [PASS] chordless_region_weight\n");
}

// ════════════════════════════════════════════════════════════════
//  8. remove_chord arc merging verification
// ════════════════════════════════════════════════════════════════

static void test_remove_chord_4_adj_arcs() {
    // Build a submap where a chord's endpoints are NOT polygon vertices
    // (y=1.5, between polygon vertices at y=1 and y=2).  Removing it
    // merges the adjacent arcs at BOTH endpoints ([C91 §2.2 tex 94]);
    // both glue mates are wrap-spanning single structures ([C91 §2.4
    // tex 142]), so the merged results are wrap arcs too.  A second
    // chord (at vertex 3) keeps the submap chordful so the mid-edge
    // glue path — not the last-chord closure — is exercised.
    //
    // Layout (test_polygon, C = vertices 0..4):
    //   r0: start-wrap arc S  R(1)→L(1)  [chord_mid's endpoints]
    //   r1: B L[1,2] and D R[2,1]        [between the two chords]
    //   r2: end-wrap arc E  L(3)→R(3)    [around C's end vertex]
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a;
    a = {}; a.first_edge = 1; a.last_edge = 2; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 2;
    std::size_t B = s.add_arc(a);
    a = {}; a.first_edge = 3; a.last_edge = 3; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = r2; a.edge_count = 1;          // end wrap: ᾱ = [3,3]
    std::size_t E = s.add_arc(a);
    a = {}; a.first_edge = 2; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r1; a.edge_count = 2;
    std::size_t D = s.add_arc(a);
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 2;          // start wrap: ᾱ = [0,1]
    std::size_t S = s.add_arc(a);

    // chord_mid at y=1.5 on edge 1 (both sides, mid-edge → 2 slots each,
    // [C91 §2.2 tex 94]): LEFT endpoint {before=S, after=B}; RIGHT
    // endpoint {before=D, after=S}.
    Chord c;
    c.region[0] = r0; c.region[1] = r1;
    c.left_adj = {{S, B}, 2}; c.right_adj = {{D, S}, 2};
    c.left_edge = 1; c.right_edge = 1;
    c.left_side = LEFT; c.right_side = RIGHT;
    c.y = 1.5; c.y_tag = 99; // tag 99 doesn't match any vertex index
    s.add_chord(c);

    // chord_v at vertex 3 (y=3, tag 3; vertex endpoints → 1 slot each,
    // the before-arc): LEFT before = B (ends at v3); RIGHT before = E
    // (its RIGHT leg [3,3] descending ends at v3).
    Chord cv;
    cv.region[0] = r1; cv.region[1] = r2;
    cv.left_adj = {{B}, 1}; cv.right_adj = {{E}, 1};
    cv.left_edge = 2; cv.right_edge = 2;
    cv.left_side = LEFT; cv.right_side = RIGHT;
    cv.y = 3.0; cv.y_tag = 3;
    s.add_chord(cv);

    auto poly = test_polygon();
    s.remove_chord(0, poly);

    // Both endpoints glue ([C91 §2.2 tex 94]): LEFT merges S+B, RIGHT
    // merges (D + the S∪B survivor) → ONE start-wrap arc R(2)→L(2) with
    // ᾱ = [0,2] (ec = 3), plus E.  [C91 §2.4 tex 142]: still single
    // structures.
    assert(s.num_live_arcs() == 2 && "two glues chain S,B,D into one arc");
    bool found_merged = false;
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        if (s.arc(i).first_side == RIGHT && s.arc(i).last_side == LEFT) {
            assert(s.arc(i).first_edge == 2 && s.arc(i).last_edge == 2 &&
                   s.arc(i).edge_count == 3 &&
                   "merged start-wrap arc must span ᾱ = [0,2]");
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

// ════════════════════════════════════════════════════════════════
//  8b. remove_chord — arcs merge at vertex endpoints too
// ════════════════════════════════════════════════════════════════

static void test_remove_chord_merge_at_vertex() {
    // Chord endpoint IS a polygon vertex → the arcs meeting there fuse:
    // [C91 §2.2 tex 96] region boundaries alternate exit chords with
    // arcs (Fig. 2.4's region II passes through removed chords' vertex
    // endpoints inside single arcs), and tex 94's vertex/non-vertex
    // distinction governs only ∂C's vertex set (the vertex persists as
    // an interior vertex of the fused arc).  Vertex endpoints carry ONE
    // adj arc (the before-arc, [C91 §2.4(ii)] 1-slot convention); the
    // after-arc is found by remove_chord's junction scan.
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0; s.end_vertex = 2;

    // C = vertices 0..2.  r1's arc double-backs around C's end vertex
    // (E: L(1)→R(1)); r0's around C's start vertex (S: R(0)→L(0)) —
    // single structures per [C91 §2.4 tex 142].
    Arc a;
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = r1; a.edge_count = 1;
    std::size_t E = s.add_arc(a);
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 1;
    std::size_t S = s.add_arc(a);

    Chord c;
    c.region[0] = r0; c.region[1] = r1;
    // Vertex endpoints: 1 adj arc per side, the before-arc — S's LEFT
    // leg ends at v1; E's RIGHT leg ends at v1.
    c.left_adj = {{S}, 1}; c.right_adj = {{E}, 1};
    // Chord at y=1.0, tag=1 — matches polygon vertex 1 at (1,1,1).
    c.left_edge = 0; c.right_edge = 0;
    c.left_side = LEFT; c.right_side = RIGHT;
    c.y = 1.0; c.y_tag = 1;
    s.add_chord(c);

    auto poly = test_polygon();
    // Pre-removal: full normal-form validation (endpoint↔count rule).
    s.check_invariants(poly);

    s.remove_chord(0, poly);

    // [C91 §2.2 tex 94/96]: the chord was the submap's LAST chord, so
    // gluing at both vertex endpoints closes ∂C into the single closed
    // arc covering all of C ([C91 §2.4 tex 142]), stored cut at C's
    // start turnaround.
    assert(s.num_live_arcs() == 1 &&
           "[C91 §2.2 tex 96]: vertex endpoints must glue their arcs");
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 &&
               "closed arc is cut at C's start turnaround "
               "([C91 §2.4(iii) tex 138])");
        assert(s.arc(i).edge_count == 2 &&
               "closed arc spans all of C through the vertex");
    }

    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] remove_chord_merge_at_vertex\n");
}


// ════════════════════════════════════════════════════════════════
//  8c. remove_chord with 2 adj_arcs (chord at C's start vertex)
// ════════════════════════════════════════════════════════════════

static void test_remove_chord_2_adj_arcs() {
    // [C91 §2.4(ii) tex 137]: "two... arcs adjacent to it."  Total
    // adj_arcs = 2 happens when both endpoints are the companion
    // duplicates of one C endpoint ([C91 §2.1 tex 72] case 3: "possibly
    // the same chord for the two companions").  The chord separates the
    // turnaround cap — bounded by the zero-length wrap arc Z passing
    // through C's start turnaround ([C91 §2.4 tex 142]) — from the main
    // region, whose single arc M double-backs around C's END vertex.
    Submap s;
    std::size_t r_main = s.add_node();
    std::size_t r_cap = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a;
    // M: from the LEFT start-companion up the whole LEFT side, around
    // C's end vertex, down the whole RIGHT side to the RIGHT
    // start-companion — one end-wrap structure, ᾱ = all of C.
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = r_main; a.edge_count = 4;
    std::size_t M = s.add_arc(a);

    // Z: the zero-length passage between the two start companions,
    // through C's start turnaround ([C91 §2.1 tex 72] case 3).
    a = {}; a.first_edge = 0; a.last_edge = 0; a.first_side = RIGHT; a.last_side = LEFT;
    a.region_node = r_cap; a.edge_count = 0;
    std::size_t Z = s.add_arc(a);

    // The companion chord at vertex 0 (y=0, tag=0): LEFT endpoint's
    // before-arc = Z (ends at the LEFT companion); RIGHT endpoint's
    // before-arc = M (ends at the RIGHT companion).
    Chord c;
    c.region[0] = r_main; c.region[1] = r_cap;
    c.left_adj = {{Z}, 1}; c.right_adj = {{M}, 1};
    c.left_edge = 0; c.right_edge = 0;
    c.left_side = LEFT; c.right_side = RIGHT;
    c.y = 0.0; c.y_tag = 0;
    s.add_chord(c);

    assert(s.start_arc == Z && s.end_arc == M &&
           "[C91 §2.4(iii) tex 138]: endpoint arcs auto-registered");

    auto poly = test_polygon();
    s.remove_chord(0, poly);

    // [C91 §2.2 tex 94]: removing the LAST chord glues at both
    // companions, closing ∂C into the single closed arc.
    assert(s.num_live_arcs() == 1 &&
           "companion-chord removal closes ∂C into one arc");
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 &&
               s.arc(i).edge_count == 4 &&
               "[C91 §2.4 tex 142]: closed arc covering all of C");
    }
    s.assert_tree_property();

    s.compact();
    s.check_invariants(poly);

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
    // [C91 §2.4(ii) tex 137]: "three... arcs adjacent to it."
    // LEFT endpoint at vertex 0 (C's start LEFT companion) → count=1.
    // RIGHT endpoint at y=0 on edge 2 (non-vertex) → count=2.
    //
    // r0's arc double-backs around C's END vertex (E: from the LEFT
    // companion of v0 all the way up and back down to the chord's
    // mid-edge endpoint on RIGHT edge 2); r1's arc double-backs around
    // C's START vertex (S: from the mid-edge endpoint down to v0 and
    // through the turnaround, ending at the LEFT companion) — single
    // structures per [C91 §2.4 tex 142].
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a;
    // E: L(0)→R(2); ᾱ = [0,3] (4 nonnull edges).
    a = {}; a.first_edge = 0; a.last_edge = 2; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = r0; a.edge_count = 4;
    std::size_t E = s.add_arc(a);

    // S: R(2)→L(0); ᾱ = [0,2] (3 nonnull edges).
    a = {}; a.first_edge = 2; a.last_edge = 0; a.first_side = RIGHT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 3;
    std::size_t S = s.add_arc(a);

    // Chord at y=0, tag=0.
    // LEFT endpoint on edge 0 → vertex 0 (match), count=1: before-arc =
    // S (its LEFT leg ends at the companion).
    // RIGHT endpoint on edge 2 → v2(y=-5,tag=2) and v3(y=5,tag=3),
    //   neither matches (0.0, 0), so non-vertex, count=2: before = E,
    //   after = S.
    Chord c;
    c.region[0] = r0; c.region[1] = r1;
    c.left_adj = {{S}, 1};
    c.right_adj = {{E, S}, 2};
    c.left_edge = 0; c.left_side = LEFT;
    c.right_edge = 2; c.right_side = RIGHT;
    c.y = 0.0; c.y_tag = 0;
    s.add_chord(c);

    assert(s.start_arc == S && s.end_arc == E);

    auto poly = non_monotone_polygon();
    s.remove_chord(0, poly);

    // The chord is the submap's LAST chord: gluing at the vertex
    // endpoint and the mid-edge endpoint chains E and S into the single
    // closed arc ([C91 §2.2 tex 94] / [C91 §2.4 tex 142]).
    assert(s.num_live_arcs() == 1 && "last-chord removal closes ∂C");
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 &&
               s.arc(i).edge_count == 4 &&
               "closed arc covering all of C");
    }
    s.assert_tree_property();

    std::printf("  [PASS] remove_chord_3_adj_arcs\n");
}

// ════════════════════════════════════════════════════════════════
//  8e. remove_chord — leaf region with a single arc adjacent at BOTH
//      (non-vertex) endpoints ([C91 §2.2 tex 94] + [C91 §2.3 tex 121])
// ════════════════════════════════════════════════════════════════

static void test_remove_chord_shared_arc_left() {
    // [C91 §2.3 tex 121]: granularity condition (ii) contracts edges
    // incident on degree-<3 nodes; a degree-1 leaf region bounded by one
    // chord and ONE arc is the canonical target.  Both chord endpoints
    // sit mid-edge on the LEFT side, so r0's single arc W runs the LONG
    // way around — through BOTH of C's endpoint turnarounds — and is
    // stored as ONE DOUBLE-WRAP structure ([C91 §2.4 tex 142]: "an arc
    // might wrap around both sides of C"): first=(3,LEFT) > last=
    // (1,LEFT) encodes the double-backing; ᾱ = all of C.
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0; s.end_vertex = 4;

    Arc a;
    // A: the leaf arc between the chord's endpoints (r1).
    a = {}; a.first_edge = 1; a.last_edge = 3; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 3;
    std::size_t A = s.add_arc(a);
    // W: double-wrap arc of r0 — LEFT [3,3] + all of RIGHT + LEFT [0,1].
    a = {}; a.first_edge = 3; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 4;
    std::size_t W = s.add_arc(a);

    Chord c;
    c.region[0] = r0; c.region[1] = r1;
    c.left_edge = 1;  c.left_side = LEFT;
    c.right_edge = 3; c.right_side = LEFT;
    c.y = 1.5; c.y_tag = 99;          // matches no vertex → both non-vertex
    c.left_adj  = {{W, A}, 2};        // W ends at chord, A starts
    c.right_adj = {{A, W}, 2};        // A ends at chord, W starts
    s.add_chord(c);

    // [C91 §2.4 tex 142/138]: the double-wrap arc passes through BOTH
    // turnarounds, so it is both endpoint arcs.
    assert(s.start_arc == W && s.end_arc == W);

    auto poly = test_polygon();

    // [C91 §2.3 tex 121]: contracting the LAST chord closes ∂C — the
    // merged weight is the nonnull count of all of C.
    assert(s.simulated_contraction_weight(0, poly) == 4);

    s.remove_chord(0, poly);
    s.assert_tree_property();
    assert(s.num_live_arcs() == 1 &&
           "shared-arc removal of the last chord closes ∂C into one arc");
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 &&
               s.arc(i).edge_count == 4 &&
               "closed arc must span the full glued range");
    }

    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] remove_chord_shared_arc_left\n");
}

static void test_remove_chord_shared_arc_right() {
    // Mirror orientation: the leaf arc A lies on the RIGHT side, so
    // r0's double-wrap arc W has the RIGHT-side encoding —
    // first=(1,RIGHT) < last=(3,RIGHT) ([C91 §2.4 tex 142]): RIGHT
    // [0,1] + all of LEFT + RIGHT [3,3]; ᾱ = all of C.
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0; s.end_vertex = 4;

    Arc a;
    // A: the leaf arc between the chord's endpoints (r1, RIGHT side).
    a = {}; a.first_edge = 3; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r1; a.edge_count = 3;
    std::size_t A = s.add_arc(a);
    // W: double-wrap arc of r0 (no LEFT-starting arcs exist).
    a = {}; a.first_edge = 1; a.last_edge = 3; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r0; a.edge_count = 4;
    std::size_t W = s.add_arc(a);

    Chord c;
    c.region[0] = r0; c.region[1] = r1;
    c.left_edge = 3;  c.left_side = RIGHT;
    c.right_edge = 1; c.right_side = RIGHT;
    c.y = 1.5; c.y_tag = 99;
    c.left_adj  = {{W, A}, 2};        // W ends at chord, A starts
    c.right_adj = {{A, W}, 2};        // A ends at chord, W starts
    s.add_chord(c);

    assert(s.start_arc == W && s.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

    auto poly = test_polygon();
    assert(s.simulated_contraction_weight(0, poly) == 4);

    s.remove_chord(0, poly);
    s.assert_tree_property();
    assert(s.num_live_arcs() == 1);
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 &&
               s.arc(i).edge_count == 4 &&
               "closed arc must span the full glued range");
    }

    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] remove_chord_shared_arc_right\n");
}

// ════════════════════════════════════════════════════════════════
//  8f. Removing a chord whose endpoints are ∂C's ONLY subdivision
//      points closes ∂C ([C91 §2.2 tex 94] + [C91 §2.4 tex 142])
// ════════════════════════════════════════════════════════════════

static void test_remove_chord_fully_wrapped_closes() {
    // Both regions are single-arc leaves (each arc wraps around one C
    // endpoint).  Removing the chord glues ∂C at both of its endpoints
    // ([C91 §2.2 tex 94]) — since they are ∂C's only subdivision
    // points, the boundary closes into the chordless submap's single
    // closed arc covering all of ∂C, stored cut at C's start
    // turnaround ([C91 §2.4 tex 142/138]).
    auto poly = test_polygon();
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0; s.end_vertex = 4;

    // A: LEFT edge 1 → wraps C's end vertex → RIGHT edge 2 (r1).
    Arc a;
    a = {}; a.first_edge = 1; a.last_edge = 2;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = r1;
    a.edge_count = poly.count_nonnull_edges(1, 3);
    std::size_t A = s.add_arc(a);
    // B: RIGHT edge 2 → wraps C's start vertex → LEFT edge 1 (r0).
    a = {}; a.first_edge = 2; a.last_edge = 1;
    a.first_side = RIGHT; a.last_side = LEFT;
    a.region_node = r0;
    a.edge_count = poly.count_nonnull_edges(0, 2);
    std::size_t B = s.add_arc(a);

    Chord c;
    c.region[0] = r0; c.region[1] = r1;
    c.left_edge = 1;  c.left_side = LEFT;
    c.right_edge = 2; c.right_side = RIGHT;
    c.y = 1.5; c.y_tag = 99;
    c.left_adj  = {{B, A}, 2};    // B ends at chord, A starts
    c.right_adj = {{A, B}, 2};    // A ends at chord, B starts
    s.add_chord(c);

    assert(s.start_arc == B && s.end_arc == A &&
           "[C91 §2.4(iii) tex 138]: wrap arcs auto-registered");

    s.remove_chord(0, poly);

    assert(s.num_live_nodes() == 1 && s.num_live_chords() == 0);
    assert(s.num_live_arcs() == 1 &&
           "[C91 §2.2 tex 94]: last-chord removal closes ∂C");
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 &&
               s.arc(i).edge_count == 4 &&
               "[C91 §2.4 tex 142]: the closed arc covers all of C");
    }
    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] remove_chord_fully_wrapped_closes\n");
}

// ════════════════════════════════════════════════════════════════
//  9. Null-length chord
// ════════════════════════════════════════════════════════════════

static void test_null_length_chord() {
    // [C91 §2.1 tex 72]: "one of these pairs... gives rise to a chord of
    // null length."  Null-length chords create empty regions.  C =
    // vertices 0..2 with the null chord at vertex 1: the outer region's
    // single arc W runs from the null position the long way around ∂C —
    // through BOTH of C's endpoint turnarounds — back to the null
    // position: ONE double-wrap structure ([C91 §2.4 tex 142]),
    // first=(1,LEFT) > last=(0,LEFT).
    Submap s;
    std::size_t r0 = s.add_node(); // main region
    std::size_t r1 = s.add_node(); // empty region inside the null-length chord
    s.start_vertex = 0;
    s.end_vertex = 2;

    // Zero-length arc inside the null-length chord's empty region.
    Arc a;
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 0;
    std::size_t N = s.add_arc(a);

    // W: outer double-wrap arc — LEFT [1,1] + all of RIGHT + LEFT [0,0].
    a = {}; a.first_edge = 1; a.last_edge = 0; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 2;
    std::size_t W = s.add_arc(a);

    // Null-length chord at vertex 1: one adj arc per slot ([C91 §2.1
    // tex 72]); the outer slot holds W (ends at the first duplicate),
    // the other the inner null arc.
    Chord c;
    c.region[0] = r0; c.region[1] = r1;
    c.left_adj = {{W}, 1}; c.right_adj = {{N}, 1};
    c.left_edge = 1; c.right_edge = 1;
    c.left_side = LEFT; c.right_side = LEFT;
    c.is_null_length = true;
    c.y = 1.0; c.y_tag = 1;
    s.add_chord(c);

    assert(s.start_arc == W && s.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

    assert(s.chord(0).is_null_length);
    assert(s.node(r1).degree() == 1); // pendant

    // [C91 §2.2 tex 106]: the null-length chord's interior is empty
    // (weight 0).
    assert(s.region_weight(r1) == 0);

    // Tree property holds.
    s.assert_tree_property();

    // Symbolic y accessor.
    SymbolicY sy = s.chord(0).symbolic_y();
    assert(sy.y == 1.0 && sy.tag == 1);

    // [C91 §2.2 tex 108]: "once removed, a chord of zero length ceases
    // to separate any arcs" — removing the LAST chord (a null one)
    // closes ∂C into the single closed arc.
    auto poly = Polygon({{0,0,0}, {1,1,1}, {2,2,2}});
    s.remove_chord(0, poly);
    assert(s.num_live_arcs() == 1 && s.num_live_nodes() == 1);
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        assert(s.arc(i).first_side == LEFT && s.arc(i).last_side == RIGHT &&
               s.arc(i).first_edge == 0 && s.arc(i).last_edge == 0 &&
               s.arc(i).edge_count == 2);
    }
    s.compact();
    s.check_invariants(poly);

    std::printf("  [PASS] null_length_chord\n");
}

// ════════════════════════════════════════════════════════════════
//  10. Null-length chord — RIGHT side
// ════════════════════════════════════════════════════════════════

static void test_null_length_chord_right_side() {
    // [C91 §2.1 tex 72]: null-length chords arise at y-extrema; the
    // "inside" pair of duplicate vertices can land on either ∂C side
    // depending on which way the curve turns.  Mirror the LEFT-side
    // test on the RIGHT side to verify the asserts don't accidentally
    // hardcode LEFT.
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 2;

    // Zero-length arc inside the null-length chord's empty region.
    Arc a;
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r1; a.edge_count = 0;
    std::size_t N = s.add_arc(a);

    // W: outer double-wrap arc — RIGHT [0,0] + all of LEFT + RIGHT
    // [1,1] ([C91 §2.4 tex 142]; the vertex-adjacent pointers record
    // the incident edge reflecting each end's true extent direction).
    a = {}; a.first_edge = 0; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r0; a.edge_count = 2;
    std::size_t W = s.add_arc(a);

    Chord c;
    c.region[0] = r0; c.region[1] = r1;
    c.left_adj = {{W}, 1}; c.right_adj = {{N}, 1};
    c.left_edge = 1; c.right_edge = 1;
    c.left_side = RIGHT; c.right_side = RIGHT;
    c.is_null_length = true;
    c.y = 1.0; c.y_tag = 1;
    s.add_chord(c);

    assert(s.start_arc == W && s.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");
    assert(s.chord(0).is_null_length);
    assert(s.node(r1).degree() == 1);
    assert(s.region_weight(r1) == 0);
    s.assert_tree_property();

    std::printf("  [PASS] null_length_chord_right_side\n");
}

// ════════════════════════════════════════════════════════════════
//  11. Death: assert_tree_property on empty submap ([C91 §2.2 tex 96])
// ════════════════════════════════════════════════════════════════

static void test_empty_submap_fires() {
    // [C91 §2.2 tex 96]: every submap is a polygonal subdivision ⟹ ≥1 region.
    assert(assert_fires([]{
        Submap s;  // 0 nodes, 0 chords.
        s.assert_tree_property();
    }));
    std::printf("  [PASS] empty_submap_fires\n");
}

// ════════════════════════════════════════════════════════════════
//  12. Deaths: add_chord null-length chord structural invariants ([C91 §2.1 tex 72])
// ════════════════════════════════════════════════════════════════

static void test_null_length_chord_mismatched_edges_fires() {
    // [C91 §2.1 tex 72]: null-length chord endpoints must coincide.
    assert(assert_fires([]{
        Submap s;
        std::size_t r0 = s.add_node();
        std::size_t r1 = s.add_node();

        Arc a;
        a = {}; a.first_edge = 0; a.last_edge = 0;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r0; a.edge_count = 1;
        std::size_t ai0 = s.add_arc(a);

        a = {}; a.first_edge = 1; a.last_edge = 1;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r1; a.edge_count = 0;
        std::size_t ai1 = s.add_arc(a);

        Chord c;
        c.region[0] = r0; c.region[1] = r1;
        c.left_edge = 0; c.right_edge = 1;            // ← mismatched edges
        c.left_side = LEFT; c.right_side = LEFT;
        c.is_null_length = true;
        c.y = 1.0; c.y_tag = 1;
        c.left_adj = {{ai0}, 1}; c.right_adj = {{ai1}, 1};
        s.add_chord(c);
    }));
    std::printf("  [PASS] null_length_chord_mismatched_edges_fires\n");
}

static void test_null_length_chord_mismatched_sides_fires() {
    // [C91 §2.1 tex 72]: both null-length chord endpoints on the same ∂C side.
    assert(assert_fires([]{
        Submap s;
        std::size_t r0 = s.add_node();
        std::size_t r1 = s.add_node();

        Arc a;
        a = {}; a.first_edge = 1; a.last_edge = 1;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r0; a.edge_count = 1;
        std::size_t ai0 = s.add_arc(a);

        a = {}; a.first_edge = 1; a.last_edge = 1;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r1; a.edge_count = 0;
        std::size_t ai1 = s.add_arc(a);

        Chord c;
        c.region[0] = r0; c.region[1] = r1;
        c.left_edge = 1; c.right_edge = 1;
        c.left_side = LEFT; c.right_side = RIGHT;     // ← mismatched sides
        c.is_null_length = true;
        c.y = 1.0; c.y_tag = 1;
        c.left_adj = {{ai0}, 1}; c.right_adj = {{ai1}, 1};
        s.add_chord(c);
    }));
    std::printf("  [PASS] null_length_chord_mismatched_sides_fires\n");
}

static void test_null_length_chord_missing_y_tag_fires() {
    // [C91 §2 tex 47] (SoS): null-length chord must carry the y-extremum vertex's tag.
    assert(assert_fires([]{
        Submap s;
        std::size_t r0 = s.add_node();
        std::size_t r1 = s.add_node();

        Arc a;
        a = {}; a.first_edge = 1; a.last_edge = 1;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r0; a.edge_count = 1;
        std::size_t ai0 = s.add_arc(a);

        a = {}; a.first_edge = 1; a.last_edge = 1;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r1; a.edge_count = 0;
        std::size_t ai1 = s.add_arc(a);

        Chord c;
        c.region[0] = r0; c.region[1] = r1;
        c.left_edge = 1; c.right_edge = 1;
        c.left_side = LEFT; c.right_side = LEFT;
        c.is_null_length = true;
        c.y = 1.0;
        // c.y_tag left at default NONE (== SOS_NONE)  ← missing tag
        c.left_adj = {{ai0}, 1}; c.right_adj = {{ai1}, 1};
        s.add_chord(c);
    }));
    std::printf("  [PASS] null_length_chord_missing_y_tag_fires\n");
}

static void test_null_length_chord_wrong_adj_count_fires() {
    // [C91 §2.1 tex 72]: null-length chord has exactly 1 adj arc per side.
    assert(assert_fires([]{
        Submap s;
        std::size_t r0 = s.add_node();
        std::size_t r1 = s.add_node();

        Arc a;
        a = {}; a.first_edge = 1; a.last_edge = 1;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r0; a.edge_count = 1;
        std::size_t ai0 = s.add_arc(a);

        a = {}; a.first_edge = 1; a.last_edge = 1;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r0; a.edge_count = 1;
        std::size_t ai0b = s.add_arc(a);

        a = {}; a.first_edge = 1; a.last_edge = 1;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r1; a.edge_count = 0;
        std::size_t ai1 = s.add_arc(a);

        Chord c;
        c.region[0] = r0; c.region[1] = r1;
        c.left_edge = 1; c.right_edge = 1;
        c.left_side = LEFT; c.right_side = LEFT;
        c.is_null_length = true;
        c.y = 1.0; c.y_tag = 1;
        c.left_adj = {{ai0, ai0b}, 2};                // ← count 2 (must be 1)
        c.right_adj = {{ai1}, 1};
        s.add_chord(c);
    }));
    std::printf("  [PASS] null_length_chord_wrong_adj_count_fires\n");
}

// ════════════════════════════════════════════════════════════════
//  13. Death: arc.edge_count cache mismatch ([C91 §2.2 tex 106])
// ════════════════════════════════════════════════════════════════

static void test_wrong_edge_count_fires() {
    // [C91 §2.2 tex 106]: arc.edge_count must match
    // polygon.count_nonnull_edges over its underlying edge range.
    assert(assert_fires([]{
        Polygon poly({{0,0,0}, {1,3,1}, {2,1,2}, {3,4,3}, {4,2,4}});
        Submap s;
        s.add_node();
        s.start_vertex = 0; s.end_vertex = 4;

        // The chordless submap's single closed arc ([C91 §2.4 tex
        // 142/138]) with a corrupted count cache.
        Arc a;
        a = {}; a.first_edge = 0; a.last_edge = 0;
        a.first_side = LEFT; a.last_side = RIGHT;
        a.region_node = 0;
        a.edge_count = 99;                            // ← wrong (truth = 4)
        s.add_arc(a);

        s.check_invariants(poly);                     // ← cache mismatch must fire
    }));
    std::printf("  [PASS] wrong_edge_count_fires\n");
}

// ════════════════════════════════════════════════════════════════

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
    std::printf("All §2.2 tests passed.\n");
    return 0;
}
