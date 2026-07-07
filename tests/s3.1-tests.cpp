// tests/s3.1-tests.cpp — Tests for [C91 §3.1]: Fusion of two submaps.

#include "merge/fusion.h"
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
// dumps disabled.  Returns true iff the child terminated by SIGABRT.
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

// ════════════════════════════════════════════════════════════════
//  Helpers
// ════════════════════════════════════════════════════════════════

// C₁ = 5 vertices (edges 0–3), with one chord at vertex 2.
static Polygon make_C1() {
    return Polygon({{0,0,0}, {1,2,1}, {2,4,2}, {3,1,3}, {4,3,4}});
}

// Build a 2-region submap of C₁ with one chord at vertex 2.
// [C91 §2.4 tex 142]: r1's arc double-backs around C₁'s end vertex
// (E: L(1)→R(1), ᾱ = [1,3]) and r0's around its start vertex
// (S: R(1)→L(1), ᾱ = [0,1]) — single structures, never split.
static Submap make_S1(const Polygon& C) {
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1
    s.start_vertex = 0; s.end_vertex = 4;

    Arc a{};
    a.first_edge = 1; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 3;
    std::size_t aiE = s.add_arc(a);

    a = {};
    a.first_edge = 1; a.last_edge = 1;
    a.first_side = RIGHT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t aiS = s.add_arc(a);

    // Chord at vertex 2 (y = 4): both endpoints are polygon-vertex
    // companions → ONE adj arc each, the before-arc ([C91 §2.2 tex 94 +
    // §2.4(ii)] 1-slot convention): S's LEFT leg ends at (1,LEFT,y₂);
    // E's RIGHT leg ends at (1,RIGHT,y₂).
    Chord c{};
    c.region[0] = 0; c.region[1] = 1;
    c.left_edge = 1; c.right_edge = 1;
    c.left_side = LEFT; c.right_side = RIGHT;
    c.y = C.vertex(2).y; c.y_tag = 2;
    c.left_adj = {{aiS}, 1};
    c.right_adj = {{aiE}, 1};
    s.add_chord(c);

    assert(s.start_arc == aiS && s.end_arc == aiE &&
           "[C91 §2.4(iii) tex 138]: endpoint arcs auto-registered");
    return s;
}

// Chain variant with a MIDDLE region carrying two plain arcs (for the
// local_shoot multi-arc selection tests): r0 — c0(v1) — r1 — c1(v3) — r2.
static Submap make_S1_chain(const Polygon& C) {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();
    s.start_vertex = 0; s.end_vertex = 4;

    Arc a{};
    a.first_edge = 1; a.last_edge = 2;
    a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 2;
    std::size_t a1 = s.add_arc(a);
    a = {};
    a.first_edge = 3; a.last_edge = 3;
    a.first_side = LEFT; a.last_side = RIGHT;   // end wrap: ᾱ = [3,3]
    a.region_node = r2; a.edge_count = 1;
    std::size_t aE = s.add_arc(a);
    a = {};
    a.first_edge = 2; a.last_edge = 1;
    a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r1; a.edge_count = 2;
    std::size_t a4 = s.add_arc(a);
    a = {};
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = RIGHT; a.last_side = LEFT;   // start wrap: ᾱ = [0,0]
    a.region_node = r0; a.edge_count = 1;
    std::size_t aS = s.add_arc(a);

    Chord c{};
    c.region[0] = r0; c.region[1] = r1;
    c.left_edge = 0; c.right_edge = 0;
    c.left_side = LEFT; c.right_side = RIGHT;
    c.y = C.vertex(1).y; c.y_tag = 1;
    c.left_adj = {{aS}, 1}; c.right_adj = {{a4}, 1};
    s.add_chord(c);
    c = {};
    c.region[0] = r1; c.region[1] = r2;
    c.left_edge = 2; c.right_edge = 2;
    c.left_side = LEFT; c.right_side = RIGHT;
    c.y = C.vertex(3).y; c.y_tag = 3;
    c.left_adj = {{a1}, 1}; c.right_adj = {{aE}, 1};
    s.add_chord(c);

    assert(s.start_arc == aS && s.end_arc == aE);
    return s;
}

// ════════════════════════════════════════════════════════════════
//  1. Fusion sequence — single chord submap
// ════════════════════════════════════════════════════════════════

static void test_fusion_sequence_basic() {
    auto C1 = make_C1();
    auto S1 = make_S1(C1);

    FusionState state;
    build_fusion_sequence(state, S1, C1);
    const auto& seq = state.sequence;

    // 1 chord × 2 endpoints + 2 companions = 4 vertices.
    assert(seq.size() == 4);

    // First = a₀ (RIGHT companion at junction vertex 4).
    assert(seq.front().is_companion);
    assert(seq.front().side == RIGHT);
    assert(seq.front().edge == 3); // junction_edge = n-2 = 3

    // Last = a_{m+1} (LEFT companion at junction vertex 4).
    assert(seq.back().is_companion);
    assert(seq.back().side == LEFT);
    assert(seq.back().edge == 3);

    // Both companions have junction vertex's y.
    assert(symbolic_y_equal(seq.front().y,
                            symbolic_y_of(C1.vertex(4))));
    assert(symbolic_y_equal(seq.back().y,
                            symbolic_y_of(C1.vertex(4))));

    // Middle two are the chord's LEFT and RIGHT endpoints.
    assert(!seq[1].is_companion);
    assert(!seq[2].is_companion);
    assert(seq[1].chord_idx == 0);
    assert(seq[2].chord_idx == 0);

    std::printf("  [PASS] fusion_sequence_basic\n");
}

// ════════════════════════════════════════════════════════════════
//  2. Fusion sequence — no chords (single-region submap)
// ════════════════════════════════════════════════════════════════

static void test_fusion_sequence_no_chords() {
    Polygon C({{0,0,0}, {1,1,1}, {2,2,2}});
    // Chordless submap: the single closed arc ([C91 §2.4 tex 142/138]).
    Submap s;
    s.add_node();
    s.start_vertex = 0; s.end_vertex = 2;
    Arc a{};
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    s.add_arc(a);

    FusionState state;
    build_fusion_sequence(state, s, C);
    const auto& seq = state.sequence;

    // No chords → only the two companions.
    assert(seq.size() == 2);
    assert(seq[0].is_companion && seq[0].side == RIGHT);
    assert(seq[1].is_companion && seq[1].side == LEFT);

    std::printf("  [PASS] fusion_sequence_no_chords\n");
}

// ════════════════════════════════════════════════════════════════
//  3. Fusion sequence — clockwise ordering
// ════════════════════════════════════════════════════════════════

static void test_fusion_sequence_ordering() {
    auto C1 = make_C1();
    auto S1 = make_S1(C1);

    FusionState state;
    build_fusion_sequence(state, S1, C1);
    const auto& seq = state.sequence;

    // The chord at vertex 2 has LEFT endpoint on edge 1 (LEFT side)
    // and RIGHT endpoint on edge 1 (RIGHT side).
    //
    // Clockwise from a₀ (RIGHT at edge 3):
    //   RIGHT side: edge 3→2→1→0
    //   LEFT side: edge 0→1→2→3
    //
    // RIGHT endpoint (edge 1, RIGHT): traversal pos = 3-1 = 2
    // LEFT endpoint (edge 1, LEFT): traversal pos = 4+1 = 5
    //
    // So RIGHT endpoint comes before LEFT endpoint in clockwise order.
    assert(seq[1].side == RIGHT);
    assert(seq[2].side == LEFT);

    std::printf("  [PASS] fusion_sequence_ordering\n");
}

// ════════════════════════════════════════════════════════════════
//  4. Companion vertices have correct SoS identity
// ════════════════════════════════════════════════════════════════

static void test_companion_identity() {
    auto C1 = make_C1();
    auto S1 = make_S1(C1);

    FusionState state;
    build_fusion_sequence(state, S1, C1);
    const auto& seq = state.sequence;

    // Both companions correspond to the junction vertex (index 4).
    SymbolicY junction_y = symbolic_y_of(C1.vertex(4));
    assert(symbolic_y_equal(seq.front().y, junction_y));
    assert(symbolic_y_equal(seq.back().y, junction_y));

    // Their edge is the last edge of C₁.
    assert(seq.front().edge == C1.num_edges() - 1);
    assert(seq.back().edge == C1.num_edges() - 1);

    std::printf("  [PASS] companion_identity\n");
}

// ════════════════════════════════════════════════════════════════
//  5. collect_region_arcs
// ════════════════════════════════════════════════════════════════

static void test_collect_region_arcs() {
    auto C1 = make_C1();
    auto S1 = make_S1(C1);

    // [C91 §2.2 tex 96 + §2.4 tex 142]: arc-structure count == region
    // degree — each region's wrap-spanning arc is ONE structure.
    auto arcs0 = collect_region_arcs(S1, 0);
    assert(arcs0.count == 1);

    auto arcs1 = collect_region_arcs(S1, 1);
    assert(arcs1.count == 1);

    std::printf("  [PASS] collect_region_arcs\n");
}

// ════════════════════════════════════════════════════════════════
//  6. local_shoot with stub oracle
// ════════════════════════════════════════════════════════════════

// Oracle that returns a fixed hit for any subarc on the RIGHT side.
struct FixedRayShooter : RayShootingOracle {
    RayHit shoot(Point p, Side /*direction*/,
                 std::size_t /*arc_idx*/,
                 const Subarc& target) const override {
        // Hit at the target's first edge, same y as p.
        RayHit h;
        h.hit = true;
        h.x = 5.0; // fixed x
        h.y = p.y;
        h.edge = target.first_edge;
        h.side = target.first_side;
        return h;
    }
};

static void test_local_shoot() {
    auto C1 = make_C1();
    auto S1 = make_S1_chain(C1);

    FixedRayShooter oracle;
    Point p{0.0, 2.0, 99};

    // Shoot from region 0 — should check region 0's arcs.
    auto hit = local_shoot(p, RIGHT, 0, S1, C1, oracle);
    assert(hit.hit);

    // Shoot from region 1 (two plain arcs).
    hit = local_shoot(p, RIGHT, 1, S1, C1, oracle);
    assert(hit.hit);

    std::printf("  [PASS] local_shoot\n");
}

// ════════════════════════════════════════════════════════════════
//  7. local_shoot — nearest hit selection
// ════════════════════════════════════════════════════════════════

// Oracle that returns different x values per subarc side.
struct DistanceRayShooter : RayShootingOracle {
    RayHit shoot(Point p, Side /*direction*/,
                 std::size_t /*arc_idx*/,
                 const Subarc& target) const override {
        RayHit h;
        h.hit = true;
        h.y = p.y;
        h.edge = target.first_edge;
        h.side = target.first_side;
        // LEFT-side arcs hit at x=10, RIGHT-side at x=3.
        h.x = (target.first_side == LEFT) ? 10.0 : 3.0;
        return h;
    }
};

static void test_local_shoot_nearest() {
    auto C1 = make_C1();
    auto S1 = make_S1_chain(C1);

    DistanceRayShooter oracle;
    Point p{1.0, 2.0, 99};

    // Region 1 has two plain arcs (LEFT a1, RIGHT a4).
    // Shooting RIGHT from x=1: hits at x=3 (RIGHT arc) and x=10 (LEFT arc).
    // Nearest in RIGHT direction = x=3 (closer to p.x=1).
    auto hit = local_shoot(p, RIGHT, 1, S1, C1, oracle);
    assert(hit.hit);
    assert(hit.x == 3.0);

    // Shooting LEFT from x=15: hits at x=3 and x=10.
    // Nearest in LEFT direction = x=10 (closer to p.x=15).
    Point p2{15.0, 2.0, 99};
    hit = local_shoot(p2, LEFT, 1, S1, C1, oracle);
    assert(hit.hit);
    assert(hit.x == 10.0);

    std::printf("  [PASS] local_shoot_nearest\n");
}

// ════════════════════════════════════════════════════════════════
//  8. fusion_startup — Case 1 (c₀ on ∂C₂)
// ════════════════════════════════════════════════════════════════

// Oracle that hits at a fixed x on the target's first_edge, with
// hit.side derived geometrically from the edge's orientation per
// [C91 §3.1] / fusion.cpp::shooting_direction's inverse:
//   LEFT-going ray strikes the EAST face (ascending→RIGHT, descending→LEFT).
//   RIGHT-going ray strikes the WEST face (ascending→LEFT, descending→RIGHT).
struct StartupOracle : RayShootingOracle {
    const Polygon* C;
    double hit_x;
    StartupOracle(const Polygon* c, double x) : C(c), hit_x(x) {}
    RayHit shoot(Point p, Side dir,
                 std::size_t /*arc_idx*/,
                 const Subarc& target) const override {
        RayHit h;
        h.hit = true;
        h.x = hit_x;
        h.y = p.y;
        h.edge = target.first_edge;
        const auto& e = C->edge(target.first_edge);
        bool ascending = symbolic_y_less(
            symbolic_y_of(C->vertex(e.start_idx)),
            symbolic_y_of(C->vertex(e.end_idx)));
        // Ascending: face struck is opposite of dir.  Descending: same.
        if (ascending)
            h.side = (dir == LEFT) ? RIGHT : LEFT;
        else
            h.side = dir;
        return h;
    }
};

static void test_startup_case1() {
    // C₁ and C₂ share vertex (4,3,4).
    // C₁ = {(0,0,0), (1,2,1), (2,4,2), (3,1,3), (4,3,4)}
    // C₂ = {(4,3,4), (5,5,5), (6,1,6)}
    auto C1 = make_C1();
    Polygon C2({{4,3,4}, {5,5,5}, {6,1,6}});
    auto S1 = make_S1(C1);

    // Single-region S₂: the chordless submap's single closed arc
    // ([C91 §2.4 tex 142/138]).
    Submap S2;
    S2.add_node();
    S2.start_vertex = 0; S2.end_vertex = 2;
    Arc a{};
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = S2.add_arc(a);
    assert(S2.start_arc == ai0 && S2.end_arc == ai0);

    // a₀ = RIGHT companion at vertex 4 = (4,3).  Edge 3 ascending, RIGHT
    // (east) wall → shoot RIGHT ([C91 §2.1 tex 72]: into the free region,
    // toward C₂ which extends east).  Hits must be to the RIGHT (x > 4).
    // Oracle for S₁: hits at x=18 (far right).
    // Oracle for S₂: hits at x=5 (close right).
    // c₀ on ∂C₂ (closer) → Case 1.
    StartupOracle oracle1(&C1, 18.0);
    StartupOracle oracle2(&C2, 5.0);

    FusionState state;
    build_fusion_sequence(state, S1, C1);

    std::size_t start = fusion_startup(state, S1, C1, S2, C2,
                                        oracle1, oracle2);

    // Case 1: p = a₀, main loop starts at k=1 (arc A₁ from a₀ to a₁).
    assert(start == 1);
    assert(state.s2_region != NONE);
    assert(state.s2_region == 0); // single-region S₂
    assert(!state.chords.empty());

    std::printf("  [PASS] startup_case1\n");
}

// ════════════════════════════════════════════════════════════════
//  9. fusion_startup — Case 2 (c₀ on ∂C₁)
// ════════════════════════════════════════════════════════════════

static void test_startup_case2() {
    auto C1 = make_C1();
    Polygon C2({{4,3,4}, {5,5,5}, {6,1,6}});
    auto S1 = make_S1(C1);

    // Chordless S₂: the single closed arc ([C91 §2.4 tex 142/138]).
    Submap S2;
    S2.add_node();
    S2.start_vertex = 0; S2.end_vertex = 2;
    Arc a{};
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = S2.add_arc(a);
    assert(S2.start_arc == ai0 && S2.end_arc == ai0);

    // a₀ = RIGHT companion at vertex 4 = (4,3).  Edge 3 ascending → shoot
    // RIGHT ([C91 §2.1 tex 72]).
    // Oracle for S₁: hits at x=5 (close right, on ∂C₁).
    // Oracle for S₂: hits at x=18 (far right, on ∂C₂).
    // c₀ on ∂C₁ (closer) → Case 2.
    StartupOracle oracle1(&C1, 5.0);
    StartupOracle oracle2(&C2, 18.0);

    FusionState state;
    build_fusion_sequence(state, S1, C1);

    fusion_startup(state, S1, C1, S2, C2, oracle1, oracle2);

    // Case 2: s2_region set, chord recorded.
    assert(state.s2_region != NONE);
    assert(state.s2_region == 0);
    assert(!state.chords.empty());

    std::printf("  [PASS] startup_case2\n");
}

// ════════════════════════════════════════════════════════════════
//  10. shooting_direction — all four (side, ascent) cases
// ════════════════════════════════════════════════════════════════

static void test_shooting_direction_all_cases() {
    // [C91 §2.1 tex 72]: chord direction = "left of an observer walking
    // clockwise around ∂C".  [C91 §2.2 tex 96]: the clockwise ∂C tour
    // walks each region's boundary counterclockwise w.r.t. the region,
    // so the free region — and the chord — is on the observer's LEFT.
    // The canonical order ([C91 §2.4(iii) tex 138]) walks Side::LEFT
    // WITH the curve (ascending) and Side::RIGHT against it:
    //   walk = ±curve_dir  ⟹  chord direction = perp_CCW(walk).
    //
    // Construct edges with controlled ascending / descending y, then
    // verify all four (side × ascent) combinations.
    Polygon up({{0,0,0}, {1,2,1}});             // edge 0 ascending (y: 0 → 2).
    Polygon down({{0,2,0}, {1,0,1}});           // edge 0 descending.

    // LEFT side, ascending edge: LEFT walks with curve = (1,2);
    // observer's left = perp_CCW(1,2) = (-2,1) → horizontal -x → LEFT
    // (the west wall shoots into its free region to the west).
    assert(shooting_direction(0, LEFT, up) == LEFT);

    // LEFT side, descending edge: walk with curve = (1,-2);
    // perp_CCW = (2,1) → horizontal +x → RIGHT.
    assert(shooting_direction(0, LEFT, down) == RIGHT);

    // RIGHT side, ascending: RIGHT walks against curve = (-1,-2);
    // perp_CCW = (2,-1) → horizontal +x → RIGHT.
    assert(shooting_direction(0, RIGHT, up) == RIGHT);

    // RIGHT side, descending: walk = (-1,2);
    // perp_CCW = (-2,-1) → horizontal -x → LEFT.
    assert(shooting_direction(0, RIGHT, down) == LEFT);

    std::printf("  [PASS] shooting_direction_all_cases\n");
}

// ════════════════════════════════════════════════════════════════
//  11. local_shoot — [C91 §2.1] double-boundary tie-break at same distance
// ════════════════════════════════════════════════════════════════

// Oracle returning two arcs at the *same* x but different sides.
struct TieBreakRayShooter : RayShootingOracle {
    Side first_arc_side;                        //< side of arc 0 hit
    Side second_arc_side;                       //< side of arc 1 hit
    TieBreakRayShooter(Side s0, Side s1)
        : first_arc_side(s0), second_arc_side(s1) {}
    RayHit shoot(Point p, Side /*dir*/,
                 std::size_t arc_idx,
                 const Subarc& target) const override {
        RayHit h;
        h.hit = true;
        h.x = 5.0;                              // SAME x for both → distance tie
        h.y = p.y;
        h.edge = target.first_edge;
        h.side = (arc_idx == 0) ? first_arc_side : second_arc_side;
        return h;
    }
};

static void test_local_shoot_tie_break() {
    // [C91 §2.1 tex 72] (double boundary's snake left/right): when two
    // arcs deliver hits at exactly the same x (distance tie),
    // the disambiguation in fusion.cpp:114-145 picks the face struck
    // first by an infinitesimally thick ray, computed from the edge's
    // geometric ascent and the shooting direction:
    //   struck_first_face (dist > 0, dir=RIGHT) = minus_x_face
    //   minus_x_face for ascending edge = LEFT
    auto C1 = make_C1();                        // edge 1: (1,2)→(2,4) ascending
    auto S1 = make_S1_chain(C1);

    // Region 1 has two plain arcs (a1 = arc 0, a4 = arc 2); oracle
    // returns both at x=5 with controlled sides via arc_idx.  For
    // shooting RIGHT toward an ascending-edge hit, struck_first_face =
    // LEFT (the -x face).
    //
    // First oracle: arc 0 → LEFT, others → RIGHT.  Best stays at the
    //   first hit (LEFT, struck_first_face=LEFT, no swap needed).
    // Second oracle: arc 0 → RIGHT, others → LEFT.  Tie-break swaps to
    //   the LEFT-side hit (matching struck_first_face=LEFT).
    // Both must return side=LEFT — the [C91 §2.1 tex 72] disambiguation rule.
    Point p{0.0, 2.0, 99};

    TieBreakRayShooter oracle(LEFT, RIGHT);
    RayHit h1 = local_shoot(p, RIGHT, 1, S1, C1, oracle);
    assert(h1.hit);
    assert(h1.x == 5.0);
    assert(h1.side == LEFT &&
           "[C91 §2.1 tex 72]: struck-first face for shooting RIGHT toward "
           "an ascending-edge hit must be LEFT (the -x face)");

    TieBreakRayShooter oracle_rev(RIGHT, LEFT);
    RayHit h2 = local_shoot(p, RIGHT, 1, S1, C1, oracle_rev);
    assert(h2.hit);
    assert(h2.x == 5.0);
    assert(h2.side == LEFT &&
           "[C91 §2.1 tex 72]: tie-break must select the same face regardless "
           "of which arc the oracle reports it on");

    std::printf("  [PASS] local_shoot_tie_break\n");
}

// ════════════════════════════════════════════════════════════════
//  12. fusion_startup — d1 == d2 tie defaults to Case 1 (∂C₂)
// ════════════════════════════════════════════════════════════════

static void test_startup_d1_eq_d2_defaults_to_case1() {
    // [C91 §3.1 tex 191]: "for simplicity we still go on saying that c₀
    // sees a point of ∂C₂ with respect to C."  When hit_c1 and hit_c2
    // are at equal distance (degenerate under SoS), the implementation
    // defaults to Case 1 (c₀ ∈ ∂C₂).
    auto C1 = make_C1();
    Polygon C2({{4,3,4}, {5,5,5}, {6,1,6}});
    auto S1 = make_S1(C1);

    // Chordless S₂: the single closed arc ([C91 §2.4 tex 142/138]).
    Submap S2;
    S2.add_node();
    S2.start_vertex = 0; S2.end_vertex = 2;
    Arc a{};
    a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = S2.add_arc(a);
    assert(S2.start_arc == ai0 && S2.end_arc == ai0);

    // Both oracles return hits at the SAME x (5.0, forward of a₀'s
    // eastward shot), so d1 == d2.
    StartupOracle oracle1(&C1, 5.0);
    StartupOracle oracle2(&C2, 5.0);

    FusionState state;
    build_fusion_sequence(state, S1, C1);

    std::size_t start = fusion_startup(state, S1, C1, S2, C2, oracle1, oracle2);

    // tex 191 default: tie → Case 1 → main loop starts at k=1.
    assert(start == 1);
    assert(!state.chords.empty());

    std::printf("  [PASS] startup_d1_eq_d2_defaults_to_case1\n");
}

// ════════════════════════════════════════════════════════════════
//  13. build_fusion_sequence skips null-length chords ([C91 §3.1 tex 179] + [C91 §2.2 tex 96])
// ════════════════════════════════════════════════════════════════

static void test_build_fusion_sequence_skips_null_length_chords() {
    // [C91 §3.1 tex 179]: "Let a₁, a₂, ..., aₘ be the canonical vertex
    // enumeration of S₁ ... exit chord endpoints in S₁."  Per
    // [C91 §2.2 tex 96], exit chords are distinguished from null-length
    // chords; per [C91 §3.1 tex 224], null-length chords are "carried
    // over automatically" to
    // the fused submap and must NOT appear in the main-loop sequence.
    auto C = make_C1();                         // 5 vertices, 4 edges.

    // Submap with 1 exit chord + 1 null-length chord, in the paper's
    // extremum shape ([C91 §2.1 tex 72]): vertex 2 = (2,4) is a y-max
    // of C₁; its inside pair (LEFT) carries the null chord cN, its
    // outside pair (RIGHT) the exit chord cO joining the two outside
    // companions.  Wrap-spanning arcs are single structures
    // ([C91 §2.4 tex 142]).
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r_null = s.add_node();          // inside the null chord
    std::size_t r_cap = s.add_node();           // beyond the exit chord
    s.start_vertex = 0; s.end_vertex = 4;

    Arc a{};
    // N: the inside pair's null arc (LEFT at vertex 2).
    a.first_edge = 2; a.last_edge = 2; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r_null; a.edge_count = 0;
    std::size_t N = s.add_arc(a);
    // arc2: from the inside pair down LEFT edges 2–3, around C's end
    // vertex, up RIGHT edges 3–2 to the first outside companion.
    a = {}; a.first_edge = 2; a.last_edge = 2; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = r0; a.edge_count = 2;
    std::size_t arc2 = s.add_arc(a);
    // Z: the outside pair's zero arc (RIGHT at vertex 2).
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r_cap; a.edge_count = 0;
    std::size_t Z = s.add_arc(a);
    // arc1: from the second outside companion down RIGHT edges 1–0,
    // around C's start vertex, up LEFT edges 0–1 to the inside pair.
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = RIGHT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 2;
    std::size_t arc1 = s.add_arc(a);

    // Exit chord cO: the outside pair's chord — vertex endpoints, one
    // adj arc each ([C91 §2.2 tex 94]); sourced at vertex 2
    // ([C91 §2 tex 47]).
    Chord ec;
    ec.region[0] = r0; ec.region[1] = r_cap;
    ec.left_edge = 2; ec.left_side = RIGHT; ec.left_adj = {{arc2}, 1};
    ec.right_edge = 1; ec.right_side = RIGHT; ec.right_adj = {{Z}, 1};
    ec.y = C.vertex(2).y; ec.y_tag = 2;
    s.add_chord(ec);

    // Null-length chord cN at vertex 2 (LEFT inside pair).
    Chord null_chord;
    null_chord.region[0] = r0; null_chord.region[1] = r_null;
    null_chord.left_edge = 2; null_chord.right_edge = 2;
    null_chord.left_side = LEFT; null_chord.right_side = LEFT;
    null_chord.y = C.vertex(2).y; null_chord.y_tag = 2;
    null_chord.is_null_length = true;
    null_chord.left_adj = {{arc1}, 1};
    null_chord.right_adj = {{N}, 1};
    s.add_chord(null_chord);

    assert(s.start_arc == arc1 && s.end_arc == arc2);

    FusionState state;
    build_fusion_sequence(state, s, C);

    // Sequence = a₀ + 2 exit-chord endpoints + a_{m+1} = 4 vertices.
    // (The null-length chord is excluded; including it would give 6.)
    assert(state.sequence.size() == 4 &&
           "[C91 §3.1 tex 179]: null-length chords must not appear in the "
           "canonical vertex enumeration (paper distinguishes 'exit chord "
           "endpoints' from null-length chords per [C91 §2.2 tex 96])");

    // None of the non-companion vertices should reference the null-length chord.
    std::size_t null_chord_idx = 1;             // null-length chord is chord #1
    for (std::size_t i = 0; i < state.sequence.size(); ++i) {
        const auto& v = state.sequence[i];
        if (v.is_companion) continue;
        assert(v.chord_idx != null_chord_idx);
    }

    std::printf("  [PASS] build_fusion_sequence_skips_null_length_chords\n");
}

// ════════════════════════════════════════════════════════════════
//  13. fusion_startup — vertex-to-vertex tie-break
// ════════════════════════════════════════════════════════════════
//
// Exercises the resolve_s2_region branch where both endpoints of the
// matching S₂ chord are polygon vertices (both curve-endpoint vertices
// in practice).  Construction:
//
//   C₂'s start endpoint (= junction) and C₂'s end endpoint share
//   raw y under SoS (distinct tags).  V(C₂) has a horizontal chord
//   between them at y = junction.y, surviving γ-granularity in S₂.
//
//   C₁'s last edge ascends → leaving_downward = true and a₀ (RIGHT
//   companion, east wall) shoots RIGHT ([C91 §2.1 tex 72]).  C₂'s
//   vertices arranged so that v_end is to the RIGHT of the junction
//   (higher x), matching the shoot direction.
//
// Expected: with leaving_downward=true (moving DOWN below the chord),
// we enter the region BELOW the chord (the "outside" of the curve).
static void test_startup_vertex_to_vertex_tie_break() {
    // C₁: last edge v3(3,1) → v4(4,3) ascends → leaving_downward=true,
    // a₀_dir = RIGHT.
    auto C1 = make_C1();  // {(0,0,0), (1,2,1), (2,4,2), (3,1,3), (4,3,4)}
    auto S1 = make_S1(C1);

    // C₂: start at junction (4,3,4), goes UP-RIGHT to (5,5,5), then
    // DOWN-RIGHT to end vertex (6,3,6).  Start and end share raw y=3.
    Polygon C2({{4,3,4}, {5,5,5}, {6,3,6}});

    // S₂ structure: the chord joins the RIGHT companions of C₂'s two
    // endpoints, cutting off the pocket between the chord and the
    // peak's underside.  [C91 §2.4 tex 142]: the outer region's arc
    // runs from the junction companion around BOTH turnarounds — one
    // DOUBLE-WRAP structure W (first=(0,RIGHT) < last=(1,RIGHT)); the
    // pocket's arc B is the plain RIGHT stretch between the endpoints.
    Submap S2;
    std::size_t r_outer  = S2.add_node(); // outer region (below the chord)
    std::size_t r_pocket = S2.add_node(); // pocket above the chord
    S2.start_vertex = 0; S2.end_vertex = 2;

    Arc a{};
    a.first_edge = 1; a.last_edge = 0;
    a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r_pocket;
    a.edge_count = 2;
    std::size_t B = S2.add_arc(a);

    a = {};
    a.first_edge = 0; a.last_edge = 1;
    a.first_side = RIGHT; a.last_side = RIGHT;   // double wrap
    a.region_node = r_outer;
    a.edge_count = 2;
    std::size_t W = S2.add_arc(a);

    // Matching chord at y=3, y_tag = junction's tag (4).  Both
    // endpoints are polygon vertices → one adj slot each, the arc
    // ENDING there ([C91 §2.2 tex 94]): at the junction companion that
    // is B; at v_end's companion, W.
    Chord c{};
    c.region[0] = r_outer;
    c.region[1] = r_pocket;
    c.left_edge = 0; c.left_side = RIGHT;
    c.right_edge = 1; c.right_side = RIGHT;
    c.y = 3.0; c.y_tag = 4; // junction's index
    c.left_adj = {{B}, 1};
    c.right_adj = {{W}, 1};
    S2.add_chord(c);

    assert(S2.start_arc == W && S2.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

    // a₀ at (4,3). Edge 3 ascending, RIGHT wall → shoot RIGHT.
    // Oracle for S₁: hit at x=20 (far RIGHT, not interesting).
    // Oracle for S₂: hit at x=6 (where v_end of C₂ is).
    // c₀ on ∂C₂ (closer) → Case 1.  a₀c₀ lies on the matching chord.
    StartupOracle oracle1(&C1, 20.0);
    StartupOracle oracle2(&C2, 6.0);

    FusionState state;
    build_fusion_sequence(state, S1, C1);

    std::size_t start = fusion_startup(state, S1, C1, S2, C2,
                                        oracle1, oracle2);

    // Case 1: main loop starts at k=1.
    assert(start == 1);
    // Tie-break: leaving_downward=true → we move DOWN below the chord.
    // The pocket's arc B approaches the chord endpoints from ABOVE (via
    // the peak vertex y=5), so r_pocket is the above-chord region and
    // the outer region — whose double-wrap arc W passes below through
    // the turnarounds — is elected.
    assert(state.s2_region == r_outer &&
           "vertex-to-vertex tie-break: leaving_downward=true should "
           "elect the region below the chord (r_outer)");
    assert(!state.chords.empty());

    std::printf("  [PASS] startup_vertex_to_vertex_tie_break\n");
}

// ════════════════════════════════════════════════════════════════
//  14. fusion_startup — mid-edge matching chord tie-break
// ════════════════════════════════════════════════════════════════
//
// Exercises the mid-edge branch of resolve_s2_region.  Under SoS this
// is the typical case: the junction's horizontal ray crosses an edge
// of ∂C₂ mid-edge (not at a vertex), so the matching chord has a
// vertex endpoint at the junction and a mid-edge endpoint at the
// crossing.  The two adj arcs at the mid-edge endpoint lie on opposite
// sides of the chord; the lambda must classify them correctly via the
// structural check + adjacent vertex SoS comparison.
//
// Setup: C₂'s curve goes UP from junction (y=3) to v_mid (y=5) and
// back DOWN past y=3 to v_end (y=1).  a₀ (east wall) shoots RIGHT
// ([C91 §2.1 tex 72]); the horizontal ray from the junction at y=3
// crosses edge 1 (v_mid → v_end) mid-edge on its west wall (= RIGHT
// side of the descending edge).  Two RIGHT-side adj arcs at the
// crossing — one below (from v_end up to the crossing, ending at the
// chord) and one above (from the crossing up to v_mid, starting at it).
//
// (The null-length chord at v_mid is conceptually present per [C91 §2.1 tex 72]
// but is not modeled explicitly — the lambda processes only the matching
// chord at junction's y_tag, so other chords don't enter the matching
// branch.  This minimal S₂ exercises the mid-edge path directly.)
static void test_startup_mid_edge_tie_break() {
    auto C1 = make_C1();
    auto S1 = make_S1(C1);

    // C₂: junction at y=3, peak at y=5, end at y=1.  Edge 1 descends
    // from y=5 to y=1 crossing y=3 at x=5.5.
    Polygon C2({{4,3,4}, {5,5,5}, {6,1,6}});

    Submap S2;
    std::size_t r_above = S2.add_node();
    std::size_t r_below = S2.add_node();
    S2.start_vertex = 0; S2.end_vertex = 2;

    // [C91 §2.4 tex 142]: r_above's arc U runs from the junction
    // companion through C₂'s START turnaround and over the peak to the
    // mid-edge crossing — one start-wrap structure R(0)→R... hm — the
    // crossing is on edge 1's RIGHT face, reached along the RIGHT
    // ascent from v_end; U instead covers RIGHT edge 0 + the crossing
    // approach from above: first=(1,RIGHT), last=(0,RIGHT) would be the
    // plain stretch.  Concretely: U = from the crossing up edge 1's
    // west wall (RIGHT), over v_mid, down edge 0's east wall (RIGHT) to
    // the junction companion — plain RIGHT arc [1,0] (no turnaround in
    // its span).  W = r_below's arc from the junction companion through
    // the START turn, all of LEFT, the END turn, and up edge 1's RIGHT
    // face to the crossing — one double-wrap structure
    // (first=(0,RIGHT) < last=(1,RIGHT), [C91 §2.4 tex 142]).
    Arc a{};
    a.first_edge = 1; a.last_edge = 0;
    a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r_above;
    a.edge_count = 2;
    std::size_t U = S2.add_arc(a);

    a = {};
    a.first_edge = 0; a.last_edge = 1;
    a.first_side = RIGHT; a.last_side = RIGHT;   // double wrap
    a.region_node = r_below;
    a.edge_count = 2;
    std::size_t W = S2.add_arc(a);

    // Matching chord at y=3, y_tag=4 (junction's tag).
    // Left slot: junction companion on edge 0 (vertex endpoint, one adj
    // arc — U ends there).  Right slot: mid-edge crossing on edge 1
    // (two adj arcs — slot 0 ends at the chord (W), slot 1 starts at it
    // (U), [C91 §2.4(ii) tex 137]).
    Chord c{};
    c.region[0] = r_above;
    c.region[1] = r_below;
    c.left_edge = 0; c.left_side = RIGHT;
    c.right_edge = 1; c.right_side = RIGHT;
    c.y = 3.0; c.y_tag = 4;
    c.left_adj = {{U}, 1};
    c.right_adj = {{W, U}, 2};
    S2.add_chord(c);

    assert(S2.start_arc == W && S2.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

    // a₀ shoots RIGHT and hits the mid-edge crossing on edge 1 at x=5.5.
    StartupOracle oracle1(&C1, 20.0);
    StartupOracle oracle2(&C2, 5.5);

    FusionState state;
    build_fusion_sequence(state, S1, C1);
    std::size_t start = fusion_startup(state, S1, C1, S2, C2,
                                        oracle1, oracle2);

    assert(start == 1);
    // Tie-break: leaving_downward=true → elect region BELOW chord.
    //
    // Trace through resolve_s2_region's mid-edge branch (adj = the
    // count-2 right slot {W, U}):
    //   - W: ends at the chord (its last position is the crossing).
    //     Adjacent vertex away from the chord (RIGHT, last_edge=1) =
    //     C₂.edge(1).end_idx = v_end (y=1 < 3) → arc0_above = FALSE.
    //   - U: starts at the chord (slot 1).  Adjacent vertex (RIGHT,
    //     first_edge=1) = C₂.edge(1).start_idx = v_mid (y=5 > 3) →
    //     arc1_above = TRUE.
    //   - above_r = r_above, below_r = r_below.
    //   - leaving_downward=true → return below_r = r_below.
    assert(state.s2_region == r_below &&
           "mid-edge tie-break: leaving_downward=true should elect "
           "the region below the chord (r_below)");

    std::printf("  [PASS] startup_mid_edge_tie_break\n");
}

// ════════════════════════════════════════════════════════════════
//  15. fuse_submaps — main loop smoke test, case (i) path.
// ════════════════════════════════════════════════════════════════
//
// Forward-distance oracle: hit_x is always at `forward_dist` from p in
// the shoot direction.  hit.side is derived geometrically.  Unlike
// StartupOracle (absolute hit_x), this works for any p — the main loop
// shifts p across iterations.
//
// [C91 §3.0 tex 166]: the oracle reports the true location of the hit
// on the target arc, so the reported edge must actually span the ray's
// (perturbed) y — a horizontal ray cannot hit an edge whose y-range
// excludes it.  hit.x is still synthetic (forward_dist), which is fine
// for the smoke tests' distance comparisons.
struct ForwardOracle : RayShootingOracle {
    const Polygon* C;
    double forward_dist;
    ForwardOracle(const Polygon* c, double d) : C(c), forward_dist(d) {}
    RayHit shoot(Point p, Side dir,
                 std::size_t /*arc_idx*/,
                 const Subarc& target) const override {
        RayHit h;
        // Walk the subarc's legs in cycle order ([C91 §2.4 tex 142]:
        // wrap-spanning targets decompose per leg), each leg's edges in
        // traversal order; hit the first edge whose symbolic y-range
        // contains the ray's y.
        SymbolicY py{p.y, p.index};
        ArcLeg legs[3];
        std::size_t n = subarc_legs(target, 0, C->num_vertices() - 1,
                                    legs);
        for (std::size_t li = 0; li < n; ++li) {
            for (std::size_t k = 0; k <= legs[li].hi - legs[li].lo;
                 ++k) {
                std::size_t e = (legs[li].side == LEFT)
                    ? legs[li].lo + k    // LEFT traversal ascends
                    : legs[li].hi - k;   // RIGHT descends
                const auto& ed = C->edge(e);
                SymbolicY y0 = symbolic_y_of(C->vertex(ed.start_idx));
                SymbolicY y1 = symbolic_y_of(C->vertex(ed.end_idx));
                SymbolicY lo = symbolic_y_less(y0, y1) ? y0 : y1;
                SymbolicY hi = symbolic_y_less(y0, y1) ? y1 : y0;
                // [C91 §3.0(i) tex 169]: α' is endpoint-exact — skip
                // candidates beyond an endpoint on a shared boundary
                // edge and keep scanning.
                if (!subarc_contains_point(target, *C, e, legs[li].side,
                                           py, 0, C->num_vertices() - 1))
                    continue;
                if (symbolic_y_leq(lo, py) && symbolic_y_leq(py, hi)) {
                    h.hit = true;
                    h.y = p.y;
                    h.x = (dir == LEFT) ? (p.x - forward_dist)
                                        : (p.x + forward_dist);
                    h.edge = e;
                    h.side = legs[li].side;
                    return h;
                }
            }
        }
        h.hit = false;
        return h;
    }
};

// Smoke test: main loop runs without asserting on the make_C1 setup.
// Most iterations will not fire (i) on this specific geometry (chord
// endpoints at vertex 2 are y-extremum companions; t_dist = 0), but the
// loop must walk j = k..m+1 cleanly and terminate via (iii) or via a
// fired case + break.  Verifies:
//   (1) build_fusion_sequence + fusion_startup + main loop scaffolding
//       chain together end-to-end.
//   (2) case_i_test predicate's local_shoot calls succeed with
//       require_hit=false (no spurious asserts).
//   (3) Invariant asserts at the top of each outer iteration hold.
static void test_fuse_main_loop_smoke() {
    auto C1 = make_C1();
    auto S1 = make_S1(C1);
    Polygon C2({{4,3,4}, {5,5,5}, {6,1,6}});

    // Chordless S₂: the single closed arc ([C91 §2.4 tex 142/138]).
    Submap S2;
    S2.add_node();
    S2.start_vertex = 0; S2.end_vertex = 2;
    Arc a{};
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = S2.add_arc(a);
    assert(S2.start_arc == ai0 && S2.end_arc == ai0);

    // S₂ closer than S₁ → Case 1 startup.  Forward dist same for both
    // so startup distance comparison favors S₂ (tie → ∂C₂ per tex 191).
    ForwardOracle oracle1(&C1, 5.0);
    ForwardOracle oracle2(&C2, 1.0);

    FusionState state;
    fuse_submaps(state, S1, C1, S2, C2, oracle1, oracle2);

    // At minimum, fusion_startup recorded a₀c₀ chord.
    assert(!state.chords.empty());
    std::printf("  [PASS] fuse_main_loop_smoke (chords=%zu)\n",
                state.chords.size());
}

// ════════════════════════════════════════════════════════════════
//  16. fuse_submaps — case (ii) path exercises the loop body.
// ════════════════════════════════════════════════════════════════

// Shared fixture for the case (ii) tests (16, 18, 18b).  C₂ goes UP
// from the junction (y=3) to v_mid (y=5) and back DOWN past y=3 to
// v_end (y=1); edge 1 (v_mid→v_end) crosses y=3 at x≈2.5.  S₂ carries
// one chord ab at y=3 (the junction's y): a mid-edge crossing on edge
// 1 (LEFT) and a vertex endpoint at the junction on edge 0 (RIGHT), so
// |ab| ≈ 1.5.  [C91 §2.4 tex 142]: r_above's arc Q double-backs around
// C₂'s START vertex (from the junction companion over the peak to the
// crossing); r_below's arc P around C₂'s END vertex (from the crossing
// under the curve back to the junction companion).
static Submap make_peak_S2() {
    Submap S2;
    std::size_t r_above = S2.add_node();
    std::size_t r_below = S2.add_node();
    S2.start_vertex = 0; S2.end_vertex = 2;

    Arc a{};
    a.first_edge = 1; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;    // end wrap
    a.region_node = r_below;
    a.edge_count = 2;
    std::size_t P = S2.add_arc(a);

    a = {};
    a.first_edge = 0; a.last_edge = 1;
    a.first_side = RIGHT; a.last_side = LEFT;    // start wrap
    a.region_node = r_above;
    a.edge_count = 2;
    std::size_t Q = S2.add_arc(a);

    Chord c{};
    c.region[0] = r_above;
    c.region[1] = r_below;
    c.left_edge = 1; c.left_side = LEFT;
    c.right_edge = 0; c.right_side = RIGHT;
    c.y = 3.0; c.y_tag = 4;
    c.left_adj = {{Q, P}, 2};                    // before, after
    c.right_adj = {{P}, 1};                      // vertex endpoint
    S2.add_chord(c);

    assert(S2.start_arc == Q && S2.end_arc == P);
    return S2;
}
static Polygon make_peak_C2() {
    return Polygon({{4,3,4}, {3,5,5}, {2,1,6}});
}

// Runs the full main loop against make_peak_S2.  Goal: exercise the
// case (ii) per-endpoint test code path — R has an incident exit
// chord, so the (ii) test iterates real candidates (none ultimately
// fire for this distance setting; test 18b pins the positive case).
static void test_fuse_main_loop_case_ii_smoke() {
    auto C1 = make_C1();
    auto S1 = make_S1(C1);
    Polygon C2 = make_peak_C2();
    Submap S2 = make_peak_S2();

    // Forward oracles — distance-from-p so the main loop iterates without
    // backward-hit asserts.  S₂ closer than S₁ to keep Case 1 startup.
    ForwardOracle oracle1(&C1, 5.0);
    ForwardOracle oracle2(&C2, 1.5);

    FusionState state;
    fuse_submaps(state, S1, C1, S2, C2, oracle1, oracle2);

    assert(!state.chords.empty());
    std::printf("  [PASS] fuse_main_loop_case_ii_smoke (chords=%zu)\n",
                state.chords.size());
}

// ════════════════════════════════════════════════════════════════
//  17. rebuild_submap — [C91 §3.1 tex 226] normal-form assembly
// ════════════════════════════════════════════════════════════════

// Minimal single-region conformal submap (no chords) over `poly`: one
// region bounded by the single closed arc ([C91 §2.4 tex 142/138]).
static Submap make_chordless(const Polygon& poly) {
    Submap s;
    s.add_node();
    s.start_vertex = 0; s.end_vertex = poly.num_vertices() - 1;
    Arc a{};
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count =
        2 * poly.count_nonnull_edges(0, poly.num_edges() - 1);
    std::size_t ai0 = s.add_arc(a);
    assert(s.start_arc == ai0 && s.end_arc == ai0);
    return s;
}

static void test_rebuild_no_chords() {
    // Non-extremum junction, empty chord inventory: the fused submap is a
    // single region with one full arc per ∂C side, in normal form.
    Polygon C1({{0,0,0}, {1,1,1}});
    Polygon C2({{1,1,1}, {2,2,2}});
    Polygon C ({{0,0,0}, {1,1,1}, {2,2,2}});
    Submap S1 = make_chordless(C1);
    Submap S2 = make_chordless(C2);

    FusionState st1, st2;
    Submap out;
    rebuild_submap(out, C, S1, C1, S2, C2, st1, st2);

    assert(out.num_nodes() == 1);
    assert(out.num_chords() == 0);
    // [C91 §2.4 tex 142/138]: ONE closed arc covering all of ∂C, cut at
    // C's start turnaround.
    assert(out.num_arcs() == 1);
    assert(out.arc(0).first_side == LEFT && out.arc(0).last_side == RIGHT &&
           out.arc(0).first_edge == 0 && out.arc(0).last_edge == 0);
    assert(out.start_arc == 0 && out.end_arc == 0);
    out.assert_tree_property();
    out.check_invariants(C);

    std::printf("  [PASS] rebuild_no_chords\n");
}

static void test_rebuild_junction_extremum_inside_right() {
    // [C91 §2.1 tex 72 case 2 + §3.1 tex 224]: the junction (0,1) is a
    // y-max of C; the previous branch descends to the LEFT (x = -1), the
    // next to the RIGHT (x = +1), so the "inside of the turn" is the
    // RIGHT face of C's edge 1 — and the synthesized null-length chord
    // carries the junction vertex's own SoS tag ([C91 §2 tex 47]).
    Polygon C1({{-1,0,0}, {0,1,1}});
    Polygon C2({{0,1,1}, {1,0,2}});
    Polygon C ({{-1,0,0}, {0,1,1}, {1,0,2}});
    Submap S1 = make_chordless(C1);
    Submap S2 = make_chordless(C2);

    FusionState st1, st2;
    Submap out;
    rebuild_submap(out, C, S1, C1, S2, C2, st1, st2);

    assert(out.num_nodes() == 2);
    assert(out.num_chords() == 1);
    const Chord& nc = out.chord(0);
    assert(nc.is_null_length);
    assert(nc.y_tag == 1 &&
           "[C91 §2 tex 47]: inside pair duplicates the junction vertex — "
           "its symbolic y IS the vertex's (y, index)");
    assert(nc.left_edge == 1 && nc.right_edge == 1);
    assert(nc.left_side == RIGHT && nc.right_side == RIGHT &&
           "[C91 §2.1 tex 72]: inside of the turn is the RIGHT face here");
    assert(out.region_weight(nc.region[1]) == 0 &&
           "[C91 §2.2 tex 106]: null-length chord's inner region is empty");
    out.assert_tree_property();
    out.check_invariants(C);

    std::printf("  [PASS] rebuild_junction_extremum_inside_right\n");
}

static void test_rebuild_junction_extremum_inside_left() {
    // Mirror of the previous test: previous branch on the RIGHT (x = +1),
    // next branch to the LEFT (x = -1) — inside is the LEFT face.
    Polygon C1({{1,0,0}, {0,1,1}});
    Polygon C2({{0,1,1}, {-1,0,2}});
    Polygon C ({{1,0,0}, {0,1,1}, {-1,0,2}});
    Submap S1 = make_chordless(C1);
    Submap S2 = make_chordless(C2);

    FusionState st1, st2;
    Submap out;
    rebuild_submap(out, C, S1, C1, S2, C2, st1, st2);

    assert(out.num_nodes() == 2);
    assert(out.num_chords() == 1);
    const Chord& nc = out.chord(0);
    assert(nc.is_null_length);
    assert(nc.y_tag == 1);
    assert(nc.left_edge == 1 && nc.right_edge == 1);
    assert(nc.left_side == LEFT && nc.right_side == LEFT &&
           "[C91 §2.1 tex 72]: inside of the turn is the LEFT face here");
    assert(out.region_weight(nc.region[1]) == 0);
    out.assert_tree_property();
    out.check_invariants(C);

    std::printf("  [PASS] rebuild_junction_extremum_inside_left\n");
}

static void test_rebuild_discovered_chord_frames() {
    // [C91 §3.1 tex 224]: DiscoveredChord flags are WALKER-frame.  The
    // same physical chord — sourced at C₁'s vertex 1 (y = 2, its own SoS
    // tag per [C91 §2 tex 47]), from C₁'s edge 0 (vertex endpoint) to a
    // mid-edge crossing on C₂'s edge 0 (= C's edge 2) — must translate
    // identically whether discovered by pass 1 (walker = C₁) or pass 2
    // (walker = C₂).
    Polygon C1({{0,0,0}, {1,2,1}, {2,1,2}});
    Polygon C2({{2,1,2}, {3,3,3}});
    Polygon C ({{0,0,0}, {1,2,1}, {2,1,2}, {3,3,3}});
    Submap S1 = make_chordless(C1);
    Submap S2 = make_chordless(C2);

    auto run = [&](bool via_state1) -> Submap {
        FusionState st1, st2;
        FusionState::DiscoveredChord dc;
        dc.y = SymbolicY{2.0, 1};                 // C₁'s vertex 1
        dc.left_edge = 0;  dc.left_side = LEFT;   // on C₁ (C-edge 0)
        dc.right_edge = 0; dc.right_side = LEFT;  // on C₂ (C-edge 2)
        if (via_state1) {
            dc.left_on_walker = true;    // walker = C₁
            dc.right_on_walker = false;
            st1.chords.push_back(dc);
        } else {
            dc.left_on_walker = false;   // walker = C₂
            dc.right_on_walker = true;
            st2.chords.push_back(dc);
        }
        Submap out;
        rebuild_submap(out, C, S1, C1, S2, C2, st1, st2);
        return out;
    };

    for (bool via_state1 : {true, false}) {
        Submap out = run(via_state1);
        // The junction (2,1) is also a y-min of C, so rebuild synthesizes
        // its inside-pair null-length chord alongside the discovered one:
        // 3 regions, 2 chords ([C91 §2.1 tex 72 case 2 + §3.1 tex 224]).
        assert(out.num_nodes() == 3);
        assert(out.num_chords() == 2);
        std::size_t discovered = NONE, null_c = NONE;
        for (std::size_t ci = 0; ci < out.num_chords(); ++ci)
            (out.chord(ci).is_null_length ? null_c : discovered) = ci;
        assert(discovered != NONE && null_c != NONE);

        const Chord& c = out.chord(discovered);
        assert(c.left_edge == 0 &&
               "walker-frame translation: C₁ endpoint → C's edge 0");
        assert(c.right_edge == 2 &&
               "walker-frame translation: C₂ endpoint → C's edge 2");
        assert(c.y_tag == 1);
        // [C91 §2.4(ii) tex 137]: vertex endpoint → 1 adj arc; mid-edge
        // endpoint → 2.
        assert(c.left_adj.count == 1 && c.right_adj.count == 2);

        const Chord& nc = out.chord(null_c);
        assert(nc.y_tag == 2 &&
               "[C91 §2 tex 47]: junction null chord carries the junction "
               "vertex's own SoS tag");
        assert(nc.left_edge == 2 && nc.left_side == LEFT &&
               "[C91 §2.1 tex 72]: inside of the min turn is the LEFT "
               "face here (previous branch to the left)");
        assert(out.region_weight(nc.region[1]) == 0);

        out.assert_tree_property();
        out.check_invariants(C);
    }

    std::printf("  [PASS] rebuild_discovered_chord_frames\n");
}

static void test_rebuild_junction_null_in_input_fires() {
    // [C91 §2.1 tex 72 case 3]: the junction is a C-endpoint of each Cᵢ
    // and C-endpoints are never y-extrema, so Sᵢ cannot contain a
    // null-length chord sourced at the junction vertex.
    assert(assert_fires([]{
        Polygon C1({{0,0,0}, {1,1,1}});
        Polygon C2({{1,1,1}, {2,2,2}});
        Polygon C ({{0,0,0}, {1,1,1}, {2,2,2}});
        Submap S1 = make_chordless(C1);

        Submap S2;
        std::size_t r0 = S2.add_node();
        std::size_t r1 = S2.add_node();
        S2.start_vertex = 0; S2.end_vertex = 1;
        Arc a{};
        // Outer arc: wraps C₂'s end vertex ([C91 §2.4 tex 142]).
        a.first_edge = 0; a.last_edge = 0;
        a.first_side = LEFT; a.last_side = RIGHT;
        a.region_node = r0; a.edge_count = 1;
        std::size_t a0 = S2.add_arc(a);
        a = {}; a.first_edge = 0; a.last_edge = 0;
        a.first_side = RIGHT; a.last_side = LEFT;
        a.region_node = r1; a.edge_count = 0;
        std::size_t a_null = S2.add_arc(a);
        Chord nc{};
        nc.region[0] = r0; nc.region[1] = r1;
        nc.left_edge = 0; nc.right_edge = 0;
        nc.left_side = LEFT; nc.right_side = LEFT;
        nc.is_null_length = true;
        nc.y = 1.0; nc.y_tag = 1;                 // ← junction's tag
        nc.left_adj = {{a0}, 1}; nc.right_adj = {{a_null}, 1};
        S2.add_chord(nc);
        (void)a0; (void)a_null;

        FusionState st1, st2;
        Submap out;
        rebuild_submap(out, C, S1, C1, S2, C2, st1, st2);  // ← must fire
    }));
    std::printf("  [PASS] rebuild_junction_null_in_input_fires\n");
}

// ════════════════════════════════════════════════════════════════
//  18. case (ii) — [C91 §3.1 tex 222] "hit does not lie on ab"
// ════════════════════════════════════════════════════════════════

static void test_case_ii_hit_beyond_ab_disqualified() {
    // Same geometry as test 16 (fuse_main_loop_case_ii_smoke), but the
    // S₁ oracle reports hits far beyond the S₂ chord's other endpoint.
    // [C91 §3.1 tex 222]: such a hit "does not lie on ab" — the sightline
    // would pass through the other endpoint (a point of ∂C₂) — so the
    // candidate must be disqualified and no case (ii) chord recorded.
    auto C1 = make_C1();
    auto S1 = make_S1(C1);
    Polygon C2 = make_peak_C2();
    Submap S2 = make_peak_S2();

    // The S₂ chord spans x ∈ [~2.5, 4] (length ~1.5).  ForwardOracle
    // with forward_dist = 100 puts every S₁ hit far beyond the other
    // endpoint → all case (ii) candidates disqualified by the on-ab test.
    ForwardOracle oracle1(&C1, 100.0);
    ForwardOracle oracle2(&C2, 1.5);

    FusionState state;
    fuse_submaps(state, S1, C1, S2, C2, oracle1, oracle2);

    // Four chords: the startup a₀c₀, a₁'s and a₂'s case (i) discoveries
    // at the S₁ chord's y (the fixture's S₁ chord is an outside
    // duplicate pair at the y-maximum (2,4) — it runs through infinity
    // per [C91 §2.1 tex 70], so t is wrapped and the direct synthetic
    // S₂ hits precede it in the lexicographic ray order), and a_{m+1}'s
    // case (i) discovery (walker endpoint at the junction edge 3).
    // Crucially, NO case (ii) chord: without the on-ab
    // disqualification, the S₂ chord endpoints' far hits (distance
    // 100 ≫ |ab| ≈ 1.5) pass the back-shot test (both synthetic
    // distances equal 100) and would record spurious case (ii) chords /
    // derail the walk.  A case (ii) product would pair an S₂ chord
    // endpoint (y = {3, 4}) with a mid-arc walker point p'; every
    // recorded chord at that y must instead have its walker endpoint on
    // the junction edge 3 (startup / a_{m+1} companions).
    assert(state.chords.size() == 4 &&
           "[C91 §3.1 tex 222]: hits beyond ab must be disqualified — "
           "startup + a₁/a₂/a_{m+1} case (i) only");
    for (const auto& dc : state.chords) {
        if (!(dc.y.y == 3.0 && dc.y.tag == 4)) continue;   // S₂ chord's y
        assert(dc.left_on_walker != dc.right_on_walker);
        std::size_t walker_edge = dc.left_on_walker ? dc.left_edge
                                                    : dc.right_edge;
        assert(walker_edge == 3 &&
               "[C91 §3.1 tex 222]: chords at the S₂ chord's y must be "
               "junction-companion records, not case (ii) products");
    }

    std::printf("  [PASS] case_ii_hit_beyond_ab_disqualified\n");
}

// ════════════════════════════════════════════════════════════════
//  18b. case (ii) — POSITIVE firing ([C91 §3.1 tex 202/206/222])
// ════════════════════════════════════════════════════════════════

// True geometric ray shooter over Cᵢ restricted to the target subarc
// (copy of tests/s3.2-tests.cpp's GeomRayShooter).  [C91 §2.1 tex 70]:
// a ray that misses everything wraps through the point at infinity.
struct GeomOracle : RayShootingOracle {
    const Polygon* Ci;
    explicit GeomOracle(const Polygon* c) : Ci(c) {}

    static bool crossing_x(const Polygon& C, std::size_t e, SymbolicY sy,
                           double* x) {
        const auto ed = C.edge(e);
        const Point& vs = C.vertex(ed.start_idx);
        const Point& ve = C.vertex(ed.end_idx);
        SymbolicY y0 = symbolic_y_of(vs);
        SymbolicY y1 = symbolic_y_of(ve);
        if (symbolic_y_equal(sy, y0)) { *x = vs.x; return true; }
        if (symbolic_y_equal(sy, y1)) { *x = ve.x; return true; }
        bool between =
            (symbolic_y_less(y0, sy) && symbolic_y_less(sy, y1)) ||
            (symbolic_y_less(y1, sy) && symbolic_y_less(sy, y0));
        if (!between) return false;
        double t = (sy.y - vs.y) / (ve.y - vs.y);
        *x = vs.x + t * (ve.x - vs.x);
        return true;
    }

    RayHit shoot(Point p, Side dir, std::size_t /*arc_idx*/,
                 const Subarc& target) const override {
        SymbolicY sy{p.y, p.index};
        ArcLeg legs[3];
        std::size_t nl = subarc_legs(target, 0, Ci->num_vertices() - 1,
                                     legs);
        RayHit best;
        best.hit = false;
        double best_d = 0.0;
        for (std::size_t g = 0; g < nl; ++g) {
            for (std::size_t e = legs[g].lo; e <= legs[g].hi; ++e) {
                double x;
                if (!crossing_x(*Ci, e, sy, &x)) continue;
                const auto ed = Ci->edge(e);
                bool asc = symbolic_y_less(
                    symbolic_y_of(Ci->vertex(ed.start_idx)),
                    symbolic_y_of(Ci->vertex(ed.end_idx)));
                Side minus_x = asc ? LEFT : RIGHT;
                Side struck = (dir == RIGHT)
                    ? minus_x : (minus_x == LEFT ? RIGHT : LEFT);
                // [C91 §3.0(i) tex 169]: α' is endpoint-exact — skip
                // candidates off α' (wrong side OR beyond an endpoint
                // on a shared boundary edge) and keep scanning.
                if (!subarc_contains_point(target, *Ci, e, struck, sy,
                                           0, Ci->num_vertices() - 1))
                    continue;
                double d = (dir == RIGHT) ? (x - p.x) : (p.x - x);
                bool wrapped = (d <= 0.0);
                bool better;
                if (!best.hit) better = true;
                else if (wrapped != best.wrapped) better = !wrapped;
                else better = d < best_d;
                if (better) {
                    best.hit = true;
                    best.x = x;
                    best.y = p.y;
                    best.edge = e;
                    best.side = struck;
                    best.wrapped = wrapped;
                    best_d = d;
                }
            }
        }
        return best;
    }
};

static void test_case_ii_fires() {
    // POSITIVE complement of test 18, with REAL geometric oracles.
    // C₁'s chord vertex v2 sits at y = 6, ABOVE C₂'s y-range [1, 5]:
    // shots from a₁/a₂ (the S₁ chord endpoints, y = {6, 2}) miss
    // every C₂ edge, so a_j ∉ R and condition (i) FAILS there
    // ([C91 §3.1 tex 220]: "If there is no hit on any arc, a_j is not
    // in R"); the junction companion a_{m+1} sees ∂C₁ (C₁'s edge 2
    // crosses y=3 at x ≈ 2.545, nearer than C₂'s crossing at 2.5),
    // failing (i) too.  The S₂ chord ab spans (2.5, 3)–(4, 3) and
    // C₁'s edge 2 pierces it: endpoint a = (2.5, 3) sees
    // p' = (2.545, 3) on A_j — on ab, strictly after p, properly
    // oriented, back-shot-confirmed — so case (ii) must FIRE
    // ([C91 §3.1 tex 202/206/222]) and record the chord q → p'
    // pairing an S₂ exit-chord endpoint with a MID-ARC walker point.
    Polygon C1({{0,0,0}, {1,2,1}, {2,6,2}, {3,0.5,3}, {4,3,4}});
    auto S1 = make_S1(C1);
    Polygon C2 = make_peak_C2();

    // S₂ carrying the geometrically consistent DIRECT junction chord
    // of V(C₂): the junction companion (e0, LEFT) at (4, 3) sees
    // (2.5, 3) on e1's east face (e1, LEFT).  (make_peak_S2's chord
    // mixes slots of the two different real junction chords — fine for
    // the synthetic-oracle tests, unusable with real geometry.)
    // arc_pocket = (0,L)→(1,L) over the peak; arc_out = (1,L)→(0,L)
    // double-backs around BOTH C₂ endpoints ([C91 §2.4 tex 142]).
    Submap S2;
    std::size_t r_pocket = S2.add_node();
    std::size_t r_out = S2.add_node();
    S2.start_vertex = 0; S2.end_vertex = 2;
    Arc a{};
    a.first_edge = 0; a.first_side = LEFT;
    a.last_edge = 1; a.last_side = LEFT;
    a.region_node = r_pocket; a.edge_count = 2;
    std::size_t A_pocket = S2.add_arc(a);
    a = {};
    a.first_edge = 1; a.first_side = LEFT;
    a.last_edge = 0; a.last_side = LEFT;     // double wrap
    a.region_node = r_out; a.edge_count = 2;
    std::size_t A_out = S2.add_arc(a);
    assert(S2.start_arc == A_out && S2.end_arc == A_out);
    Chord cc{};
    cc.region[0] = r_pocket; cc.region[1] = r_out;
    cc.left_edge = 1; cc.left_side = LEFT;     // (2.5, 3) mid-edge
    cc.right_edge = 0; cc.right_side = LEFT;   // (4, 3) junction vertex
    cc.y = 3.0; cc.y_tag = 4;
    cc.left_adj = {{A_pocket, A_out}, 2};
    cc.right_adj = {{A_out}, 1};
    S2.add_chord(cc);

    GeomOracle oracle1(&C1);
    GeomOracle oracle2(&C2);

    FusionState state;
    fuse_submaps(state, S1, C1, S2, C2, oracle1, oracle2);

    // The case (ii) product: y = ab's {3, 4}; left slot = p' = C₁'s
    // (edge 2, LEFT) at x ≈ 2.545; right slot = q = b = C₂'s
    // (edge 0, LEFT) at x 4 (ascending x, [C91 §2.4(ii)]).  The other
    // endpoint a's symmetric candidate is killed by the back-shot test
    // (its sightline to (4, 3) on ∂C₁ is blocked by p' itself).
    // Startup/companion
    // records keep their walker endpoint on the junction edge 3 — a
    // mid-arc walker endpoint is case (ii)'s signature.
    bool case_ii_product = false;
    for (const auto& dc : state.chords) {
        if (!(dc.y.y == 3.0 && dc.y.tag == 4)) continue;   // ab's y
        if (dc.left_on_walker == dc.right_on_walker) continue;
        std::size_t walker_edge = dc.left_on_walker ? dc.left_edge
                                                    : dc.right_edge;
        if (walker_edge == 3) continue;         // companion record
        assert(dc.left_on_walker && dc.left_edge == 2 &&
               dc.left_side == LEFT &&
               "[C91 §3.1 tex 206]: case (ii) chord's p' lies mid-arc "
               "on C₁ edge 2 LEFT at (≈2.545, 3)");
        assert(!dc.right_on_walker && dc.right_edge == 0 &&
               dc.right_side == LEFT &&
               "[C91 §3.1 tex 206]: case (ii) chord's q = the S₂ exit "
               "chord endpoint b at (4, 3) on C₂ edge 0 LEFT");
        case_ii_product = true;
    }
    assert(case_ii_product &&
           "[C91 §3.1 tex 202/206]: an on-ab, after-p, back-shot-"
           "confirmed candidate must fire case (ii) and record the "
           "exit-chord-endpoint → p' chord");

    std::printf("  [PASS] case_ii_fires (chords=%zu)\n",
                state.chords.size());
}

// ════════════════════════════════════════════════════════════════
//  19. [C91 §3.1 tex 179] symmetric pass — junction at the FIRST vertex
// ════════════════════════════════════════════════════════════════

static void test_fusion_sequence_junction_at_start() {
    // Same submap as test 1, but walked as the SECOND pass: the junction
    // is the walked curve's FIRST vertex, so the cw tour starts at the
    // LEFT companion at edge 0, walks the LEFT half first (ascending),
    // then the RIGHT half, ending at the RIGHT companion.
    auto C1 = make_C1();
    auto S1 = make_S1(C1);

    FusionState state;
    state.junction_at_end = false;
    build_fusion_sequence(state, S1, C1);
    const auto& seq = state.sequence;

    assert(seq.size() == 4);
    assert(seq.front().is_companion && seq.front().side == LEFT &&
           "[C91 §3.1 tex 179]: junction-at-start tour begins at the "
           "LEFT companion");
    assert(seq.front().edge == 0);
    assert(seq.back().is_companion && seq.back().side == RIGHT &&
           "[C91 §3.1 tex 179]: junction-at-start tour ends at the "
           "RIGHT companion");
    assert(seq.back().edge == 0);
    // Both companions carry the junction (= vertex 0) symbolic y.
    SymbolicY jy = symbolic_y_of(C1.vertex(0));
    assert(symbolic_y_equal(seq.front().y, jy));
    assert(symbolic_y_equal(seq.back().y, jy));

    // LEFT half comes first: the chord's LEFT endpoint (edge 1, LEFT,
    // tour pos 1) precedes its RIGHT endpoint (edge 1, RIGHT, pos 6) —
    // the mirror image of the pass-1 ordering (test 3).
    assert(seq[1].side == LEFT);
    assert(seq[2].side == RIGHT);

    std::printf("  [PASS] fusion_sequence_junction_at_start\n");
}

static void test_startup_case1_junction_at_start() {
    // Walked curve W's FIRST vertex (4,3,4) is the junction; the target
    // T's LAST vertex is the same point ([C91 §3 tex 160]).
    Polygon W({{4,3,4}, {5,5,5}, {6,1,6}});
    auto T = make_C1();                 // {(0,0,0), ..., (4,3,4)}
    auto S_T = make_S1(T);

    // Chordless S_W: the single closed arc ([C91 §2.4 tex 142/138]).
    Submap S_W;
    S_W.add_node();
    S_W.start_vertex = 0; S_W.end_vertex = 2;
    Arc a{};
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = S_W.add_arc(a);
    assert(S_W.start_arc == ai0 && S_W.end_arc == ai0);

    // a₀ = LEFT companion at W's vertex 0 = (4,3).  W's edge 0 ascends
    // (3→5): the LEFT (west) wall shoots LEFT ([C91 §2.1 tex 72] — toward
    // the target T, which extends west).  Walker hit far (x=-10, d=14);
    // target hit near (x=3, d=1) → c₀ ∈ ∂C_target → Case 1.
    StartupOracle oracleW(&W, -10.0);
    StartupOracle oracleT(&T, 3.0);

    FusionState state;
    state.junction_at_end = false;
    build_fusion_sequence(state, S_W, W);
    std::size_t start = fusion_startup(state, S_W, W, S_T, T,
                                        oracleW, oracleT);

    assert(start == 1);
    // [C91 §3.1 tex 185]: a₀'s target region comes from the target's
    // C-endpoint pointers — the junction is the TARGET's LAST vertex,
    // so end_arc's region (make_S1: ai1 → region 1).
    assert(state.s2_region == S_T.arc(S_T.end_arc).region_node);
    assert(state.s2_region == 1);
    assert(!state.chords.empty());
    // Startup chord flags are walker-frame; the a₀ endpoint is on W.
    const auto& dc = state.chords[0];
    assert(dc.left_on_walker || dc.right_on_walker);

    std::printf("  [PASS] startup_case1_junction_at_start\n");
}

static void test_fuse_main_loop_smoke_junction_at_start() {
    // Full main loop in the mirrored orientation with the honest
    // ForwardOracle.  W (walked) has the junction (0,0,10) FIRST; the
    // target T ends at it.  Indices are consecutive per [C91 §2.4 tex 133]
    // (contiguous subchains of P).
    Polygon T({{-2,5,8}, {-1,1,9}, {0,0,10}});
    Polygon W({{0,0,10}, {1,2,11}, {2,4,12}, {3,1,13}, {4,3,14}});

    // make_S1-shaped 2-region submap over W (chord at W's vertex 2,
    // whose SoS index is 12): both region arcs are single wrap
    // structures ([C91 §2.4 tex 142]).
    Submap S_W;
    S_W.add_node();
    S_W.add_node();
    S_W.start_vertex = 0; S_W.end_vertex = 4;
    Arc a{};
    a.first_edge = 1; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = RIGHT;    // end wrap (r1)
    a.region_node = 1; a.edge_count = 3;
    std::size_t aiE = S_W.add_arc(a);
    a = {};
    a.first_edge = 1; a.last_edge = 1;
    a.first_side = RIGHT; a.last_side = LEFT;    // start wrap (r0)
    a.region_node = 0; a.edge_count = 2;
    std::size_t aiS = S_W.add_arc(a);
    Chord c{};
    c.region[0] = 0; c.region[1] = 1;
    c.left_edge = 1; c.right_edge = 1;
    c.left_side = LEFT; c.right_side = RIGHT;
    c.y = W.vertex(2).y; c.y_tag = W.vertex(2).index;
    c.left_adj = {{aiS}, 1};
    c.right_adj = {{aiE}, 1};
    S_W.add_chord(c);
    assert(S_W.start_arc == aiS && S_W.end_arc == aiE);

    // Chordless S_T: the single closed arc ([C91 §2.4 tex 142/138]).
    Submap S_T;
    S_T.add_node();
    S_T.start_vertex = 0; S_T.end_vertex = 2;
    a = {};
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t t0 = S_T.add_arc(a);
    assert(S_T.start_arc == t0 && S_T.end_arc == t0);

    ForwardOracle oracleW(&W, 5.0);
    ForwardOracle oracleT(&T, 1.0);

    FusionState state;
    state.junction_at_end = false;
    fuse_submaps(state, S_W, W, S_T, T, oracleW, oracleT);

    assert(!state.chords.empty());
    std::printf("  [PASS] fuse_main_loop_smoke_junction_at_start "
                "(chords=%zu)\n", state.chords.size());
}

// ════════════════════════════════════════════════════════════════
//  20. rebuild_submap — [C91 §3.1 tex 224] inventory deduplication
// ════════════════════════════════════════════════════════════════

static void test_rebuild_dedup() {
    // The junction-companion chords are shared between the two passes
    // ([C91 §3.1 tex 179]: pass 1's a_{m+1} is pass 2's a₀ as a point of
    // ∂C), so both passes can record the SAME chord.  Lemma 2.2's tree
    // structure admits each visible pair once — rebuild must dedup.
    // Same fixture as test_rebuild_discovered_chord_frames, but the
    // chord is fed through BOTH states.
    Polygon C1({{0,0,0}, {1,2,1}, {2,1,2}});
    Polygon C2({{2,1,2}, {3,3,3}});
    Polygon C ({{0,0,0}, {1,2,1}, {2,1,2}, {3,3,3}});
    Submap S1 = make_chordless(C1);
    Submap S2 = make_chordless(C2);

    FusionState st1, st2;
    st2.junction_at_end = false;
    FusionState::DiscoveredChord dc;
    dc.y = SymbolicY{2.0, 1};
    dc.left_edge = 0;  dc.left_side = LEFT;   // on C₁ (C-edge 0)
    dc.right_edge = 0; dc.right_side = LEFT;  // on C₂ (C-edge 2)
    dc.left_on_walker = true;                 // pass 1: walker = C₁
    dc.right_on_walker = false;
    st1.chords.push_back(dc);
    dc.left_on_walker = false;                // pass 2: walker = C₂
    dc.right_on_walker = true;
    st2.chords.push_back(dc);

    Submap out;
    rebuild_submap(out, C, S1, C1, S2, C2, st1, st2);

    // Identical output to the single-record run: the duplicate is
    // dropped; only the discovered chord + the junction's synthesized
    // null-length chord remain.
    assert(out.num_chords() == 2 &&
           "[C91 §3.1 tex 224]: duplicate chords must be deduplicated");
    assert(out.num_nodes() == 3);
    out.assert_tree_property();
    out.check_invariants(C);

    std::printf("  [PASS] rebuild_dedup\n");
}

// ════════════════════════════════════════════════════════════════
//  20b. rebuild_submap — cross-pass junction labels dedup
// ════════════════════════════════════════════════════════════════

static void test_rebuild_dedup_junction_cross_edge_labels() {
    // [C91 §3.0(i) tex 169]: a report names "the edge of P that
    // contains" the point hit; at a NON-extremum junction the same ∂C
    // point is contained in both incident edges, so pass 1's a_{m+1}
    // record labels it with C₁'s last edge while pass 2's startup
    // labels it with C₂'s first ([C91 §3.1 tex 179]).  Rebuild must
    // canonicalize the labels before dedup ([C91 §3.1 tex 224] +
    // Lemma 2.2: each visible pair appears once) — without it the
    // duplicate survives and the parenthesis sweep aborts.
    //
    // Geometry: junction v1=(2,2) is monotone (non-extremum); C₂ loops
    // back west so the junction's west-wall companion sees C₂'s edge
    // v3→v4 (crosses y=2 at x=−0.75).  Full two-pass fusion with real
    // geometric oracles.
    Polygon P({{0,0,0}, {2,2,1}, {4,4,2}, {0,5,3}, {-1,1,4}});
    Polygon C1 = P.subchain(0, 2);
    Polygon C2 = P.subchain(1, 4);
    Submap S1 = make_chordless(C1);
    Submap S2 = make_chordless(C2);

    GeomOracle o1(&C1);
    GeomOracle o2(&C2);

    FusionState st1;
    st1.junction_at_end = true;
    fuse_submaps(st1, S1, C1, S2, C2, o1, o2);

    FusionState st2;
    st2.junction_at_end = false;
    fuse_submaps(st2, S2, C2, S1, C1, o2, o1);

    Polygon C(C1, C2);
    Submap out;
    rebuild_submap(out, C, S1, C1, S2, C2, st1, st2);

    // One chord per junction companion side (each same-side pass 1 /
    // pass 2 record names the same ∂C point): the east walls see each
    // other across the reflex pocket, the west walls across the far
    // C₂ edge — 2 chords, 3 regions.
    assert(out.num_chords() == 2 &&
           "[C91 §3.1 tex 224]: cross-pass junction records must dedup "
           "(same ∂C point, different incident-edge labels)");
    assert(out.num_nodes() == 3);
    out.assert_tree_property();
    out.check_invariants(C);

    std::printf("  [PASS] rebuild_dedup_junction_cross_edge_labels\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::setbuf(stdout, nullptr);
    std::printf("[C91 §3.1 tests]:\n");
    test_fusion_sequence_basic();
    test_fusion_sequence_no_chords();
    test_fusion_sequence_ordering();
    test_companion_identity();
    test_collect_region_arcs();
    test_local_shoot();
    test_local_shoot_nearest();
    test_startup_case1();
    test_startup_case2();
    test_shooting_direction_all_cases();
    test_local_shoot_tie_break();
    test_startup_d1_eq_d2_defaults_to_case1();
    test_build_fusion_sequence_skips_null_length_chords();
    test_startup_vertex_to_vertex_tie_break();
    test_startup_mid_edge_tie_break();
    test_fuse_main_loop_smoke();
    test_fuse_main_loop_case_ii_smoke();
    test_rebuild_no_chords();
    test_rebuild_junction_extremum_inside_right();
    test_rebuild_junction_extremum_inside_left();
    test_rebuild_discovered_chord_frames();
    test_rebuild_junction_null_in_input_fires();
    test_case_ii_hit_beyond_ab_disqualified();
    test_case_ii_fires();
    test_fusion_sequence_junction_at_start();
    test_startup_case1_junction_at_start();
    test_fuse_main_loop_smoke_junction_at_start();
    test_rebuild_dedup();
    test_rebuild_dedup_junction_cross_edge_labels();
    std::printf("All §3.1 tests passed.\n");
    return 0;
}
