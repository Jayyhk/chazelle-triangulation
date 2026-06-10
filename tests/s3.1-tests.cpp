// tests/s3.1-tests.cpp — Tests for [C91 §3.1]: Fusion of two submaps.

#include "merge/fusion.h"
#include "polygon/polygon.h"
#include "submap/submap.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════
//  Helpers
// ════════════════════════════════════════════════════════════════

// C₁ = 5 vertices (edges 0–3), with one chord at vertex 2.
static Polygon make_C1() {
    return Polygon({{0,0,0}, {1,2,1}, {2,4,2}, {3,1,3}, {4,3,4}});
}

// Build a 2-region submap of C₁ with one chord at vertex 2.
static Submap make_S1(const Polygon& C) {
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1

    Arc a{};
    // LEFT arcs
    a.first_edge = 0; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = s.add_arc(a);

    a.first_edge = 1; a.last_edge = 3;
    a.region_node = 1; a.edge_count = 3;
    std::size_t ai1 = s.add_arc(a);

    // RIGHT arcs
    a = {};
    a.first_edge = 3; a.last_edge = 1;
    a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 3;
    std::size_t ai2 = s.add_arc(a);

    a.first_edge = 1; a.last_edge = 0;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai3 = s.add_arc(a);

    Chord c{};
    c.region[0] = 0; c.region[1] = 1;
    c.left_edge = 1; c.right_edge = 1;
    c.left_side = LEFT; c.right_side = RIGHT;
    c.y = C.vertex(2).y; c.y_tag = 2;
    c.left_adj = {{ai0, ai1}, 2};
    c.right_adj = {{ai2, ai3}, 2};
    s.add_chord(c);

    s.start_arc = ai0; s.end_arc = ai1;
    s.start_vertex = 0; s.end_vertex = 4;

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
    Submap s;
    s.add_node();
    Arc a{};
    a.first_edge = 0; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    s.add_arc(a);
    a.first_edge = 1; a.last_edge = 0;
    a.first_side = RIGHT; a.last_side = RIGHT;
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

    // Region 0 has 2 arcs (LEFT ai0 and RIGHT ai3).
    auto arcs0 = collect_region_arcs(S1, 0);
    assert(arcs0.count == 2);

    // Region 1 has 2 arcs (LEFT ai1 and RIGHT ai2).
    auto arcs1 = collect_region_arcs(S1, 1);
    assert(arcs1.count == 2);

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
    auto S1 = make_S1(C1);

    FixedRayShooter oracle;
    Point p{0.0, 2.0, 99};

    // Shoot from region 0 — should check region 0's arcs.
    auto hit = local_shoot(p, RIGHT, 0, S1, C1, oracle);
    assert(hit.hit);

    // Shoot from region 1.
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
    auto S1 = make_S1(C1);

    DistanceRayShooter oracle;
    Point p{1.0, 2.0, 99};

    // Shooting RIGHT from x=1: hits at x=3 (RIGHT arc) and x=10 (LEFT arc).
    // Nearest in RIGHT direction = x=3 (closer to p.x=1).
    auto hit = local_shoot(p, RIGHT, 0, S1, C1, oracle);
    assert(hit.hit);
    assert(hit.x == 3.0);

    // Shooting LEFT from x=15: hits at x=3 and x=10.
    // Nearest in LEFT direction = x=10 (closer to p.x=15).
    Point p2{15.0, 2.0, 99};
    hit = local_shoot(p2, LEFT, 0, S1, C1, oracle);
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

    // Single-region S₂.
    Submap S2;
    S2.add_node();
    Arc a{};
    a.first_edge = 0; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = S2.add_arc(a);
    a.first_edge = 1; a.last_edge = 0;
    a.first_side = RIGHT; a.last_side = RIGHT;
    S2.add_arc(a);
    // [C91 §2.4] (tex 144): end_arc = last LEFT arc (ai0).
    S2.start_arc = ai0; S2.end_arc = ai0;
    S2.start_vertex = 0; S2.end_vertex = 2;

    // a₀ at vertex 4 = (4,3). Edge 3 ascending → shoot LEFT → p.x ≈ 4.
    // Hits must be to the LEFT (x < 4).
    // Oracle for S₁: hits at x=-10 (far left).
    // Oracle for S₂: hits at x=3 (close left).
    // c₀ on ∂C₂ (closer) → Case 1.
    StartupOracle oracle1(&C1, -10.0);
    StartupOracle oracle2(&C2, 3.0);

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

    Submap S2;
    S2.add_node();
    Arc a{};
    a.first_edge = 0; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = S2.add_arc(a);
    a.first_edge = 1; a.last_edge = 0;
    a.first_side = RIGHT; a.last_side = RIGHT;
    S2.add_arc(a);
    // [C91 §2.4] (tex 144): end_arc = last LEFT arc (ai0).
    S2.start_arc = ai0; S2.end_arc = ai0;
    S2.start_vertex = 0; S2.end_vertex = 2;

    // a₀ at vertex 4 = (4,3). Edge 3 ascending → shoot LEFT → p.x ≈ 4.
    // Oracle for S₁: hits at x=3 (close left, on ∂C₁).
    // Oracle for S₂: hits at x=-10 (far left, on ∂C₂).
    // c₀ on ∂C₁ (closer) → Case 2.
    StartupOracle oracle1(&C1, 3.0);
    StartupOracle oracle2(&C2, -10.0);

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
    // clockwise around ∂C".  Under [C91 §2.1 tex 66]'s induced orientation,
    // RIGHT-side walks WITH curve, LEFT-side walks AGAINST curve.
    //   walk = ±curve_dir  ⟹  observer's left = perp_CCW(walk).
    //
    // Construct edges with controlled ascending / descending y, then
    // verify all four (side × ascent) combinations.
    Polygon up({{0,0,0}, {1,2,1}});             // edge 0 ascending (y: 0 → 2).
    Polygon down({{0,2,0}, {1,0,1}});           // edge 0 descending.

    // LEFT side, ascending edge: LEFT walks against curve = (-1,-2);
    // observer's left = perp_CCW(-1,-2) = (2,-1) → horizontal +x → RIGHT.
    assert(shooting_direction(0, LEFT, up) == RIGHT);

    // LEFT side, descending edge: walk against descending curve = (-1, +2);
    // perp_CCW = (-2, -1) → horizontal -x → LEFT.
    assert(shooting_direction(0, LEFT, down) == LEFT);

    // RIGHT side, ascending: RIGHT walks with curve = (1,2);
    // perp_CCW = (-2, 1) → horizontal -x → LEFT.
    assert(shooting_direction(0, RIGHT, up) == LEFT);

    // RIGHT side, descending: walk = (1, -2);
    // perp_CCW = (2, 1) → horizontal +x → RIGHT.
    assert(shooting_direction(0, RIGHT, down) == RIGHT);

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
    auto S1 = make_S1(C1);

    // Region 0 has two arcs; oracle returns both at x=5 with controlled
    // sides via arc_idx.  For shooting RIGHT toward an ascending-edge
    // hit, struck_first_face = LEFT (the -x face).
    //
    // First oracle: arc 0 → LEFT, others → RIGHT.  Best stays at the
    //   first hit (LEFT, struck_first_face=LEFT, no swap needed).
    // Second oracle: arc 0 → RIGHT, others → LEFT.  Tie-break swaps to
    //   the LEFT-side hit (matching struck_first_face=LEFT).
    // Both must return side=LEFT — the [C91 §2.1 tex 72] disambiguation rule.
    Point p{0.0, 2.0, 99};

    TieBreakRayShooter oracle(LEFT, RIGHT);
    RayHit h1 = local_shoot(p, RIGHT, 0, S1, C1, oracle);
    assert(h1.hit);
    assert(h1.x == 5.0);
    assert(h1.side == LEFT &&
           "[C91 §2.1 tex 72]: struck-first face for shooting RIGHT toward "
           "an ascending-edge hit must be LEFT (the -x face)");

    TieBreakRayShooter oracle_rev(RIGHT, LEFT);
    RayHit h2 = local_shoot(p, RIGHT, 0, S1, C1, oracle_rev);
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

    Submap S2;
    S2.add_node();
    Arc a{};
    a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = S2.add_arc(a);
    a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    S2.add_arc(a);
    S2.start_arc = ai0; S2.end_arc = ai0;
    S2.start_vertex = 0; S2.end_vertex = 2;

    // Both oracles return hits at the SAME x (3.0), so d1 == d2.
    StartupOracle oracle1(&C1, 3.0);
    StartupOracle oracle2(&C2, 3.0);

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

    // Synthetic submap with 1 exit chord + 1 null-length chord.
    // Null-length chords are at y-extrema (vertex 2 = (2, 4) is a y-max in C1).
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();              // for the exit chord
    std::size_t r_null = s.add_node();          // empty region inside the null-length chord

    Arc a{};
    a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 2;
    std::size_t ai0 = s.add_arc(a);

    // Null-length empty arc on LEFT at vertex 2 (y-max).
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r_null; a.edge_count = 0;   // null-length
    std::size_t ai_null = s.add_arc(a);

    a = {}; a.first_edge = 1; a.last_edge = 3; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 3;
    std::size_t ai1 = s.add_arc(a);

    a = {}; a.first_edge = 3; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r1; a.edge_count = 3;
    std::size_t ai2 = s.add_arc(a);

    a = {}; a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r0; a.edge_count = 2;
    std::size_t ai3 = s.add_arc(a);

    // Exit chord: r0 ↔ r1 at vertex 2's y level.
    Chord ec;
    ec.region[0] = r0; ec.region[1] = r1;
    ec.left_edge = 1; ec.right_edge = 1;
    ec.left_side = LEFT; ec.right_side = RIGHT;
    ec.y = C.vertex(2).y; ec.y_tag = 22;
    ec.left_adj = {{ai0, ai1}, 2};
    ec.right_adj = {{ai2, ai3}, 2};
    s.add_chord(ec);

    // Null-length chord: r0 ↔ r_null at vertex 2 LEFT side.
    Chord null_chord;
    null_chord.region[0] = r0; null_chord.region[1] = r_null;
    null_chord.left_edge = 1; null_chord.right_edge = 1;
    null_chord.left_side = LEFT; null_chord.right_side = LEFT;
    null_chord.y = C.vertex(2).y; null_chord.y_tag = 2;
    null_chord.is_null_length = true;
    null_chord.left_adj = {{ai0}, 1};
    null_chord.right_adj = {{ai_null}, 1};
    s.add_chord(null_chord);

    s.start_arc = ai0; s.end_arc = ai1;
    s.start_vertex = 0; s.end_vertex = 4;

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
//   C₁'s last edge ascends → leaving_downward = true and a₀ shoots LEFT.
//   C₂'s vertices arranged so that v_end is to the LEFT of junction
//   (lower x), matching the shoot direction.
//
// Expected: with leaving_downward=true (moving DOWN below the chord),
// we enter the region BELOW the chord (the "outside" of the curve).
static void test_startup_vertex_to_vertex_tie_break() {
    // C₁: last edge v3(3,1) → v4(4,3) ascends → leaving_downward=true,
    // a₀_dir = LEFT.
    auto C1 = make_C1();  // {(0,0,0), (1,2,1), (2,4,2), (3,1,3), (4,3,4)}
    auto S1 = make_S1(C1);

    // C₂: start at junction (4,3,4), goes UP-LEFT to (3,5,5), then
    // DOWN-LEFT to end vertex (2,3,6).  Start and end share raw y=3.
    Polygon C2({{4,3,4}, {3,5,5}, {2,3,6}});

    // S₂ structure (minimal: the null-length chord at v_mid is
    // γ-removed, leaving just one big arc per ∂C side and the matching
    // chord between two regions).
    Submap S2;
    std::size_t r_loop    = S2.add_node(); // contains the curve arcs (above chord)
    std::size_t r_outside = S2.add_node(); // empty region below chord

    // Big LEFT-side arc covering edges 0–1 (v4 → v5 → v6 on LEFT side).
    Arc a{};
    a.first_edge = 0; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r_loop;
    a.edge_count = 2;
    std::size_t ai_left = S2.add_arc(a);

    // Big RIGHT-side arc covering edges 1 → 0 (descending: v6 → v5 → v4).
    a = {};
    a.first_edge = 1; a.last_edge = 0;
    a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r_loop;
    a.edge_count = 2;
    std::size_t ai_right = S2.add_arc(a);

    // Matching chord at y=3, y_tag = junction's tag (4).
    // Both endpoints are polygon vertices.
    Chord c{};
    c.region[0] = r_loop;
    c.region[1] = r_outside;
    c.left_edge = 1; c.left_side = LEFT;
    c.right_edge = 0; c.right_side = RIGHT;
    c.y = 3.0; c.y_tag = 4; // junction's index
    c.left_adj = {{ai_left}, 1};
    c.right_adj = {{ai_right}, 1};
    S2.add_chord(c);

    S2.start_arc = ai_left; S2.end_arc = ai_left;
    S2.start_vertex = 0; S2.end_vertex = 2;

    // a₀ at (4,3). Edge 3 ascending → shoot LEFT (toward lower x).
    // Oracle for S₁: hit at x=-10 (far LEFT, not interesting).
    // Oracle for S₂: hit at x=2 (where v_end of C₂ is).
    // c₀ on ∂C₂ (closer) → Case 1.  a₀c₀ lies on the matching chord.
    StartupOracle oracle1(&C1, -10.0);
    StartupOracle oracle2(&C2, 2.0);

    FusionState state;
    build_fusion_sequence(state, S1, C1);

    std::size_t start = fusion_startup(state, S1, C1, S2, C2,
                                        oracle1, oracle2);

    // Case 1: main loop starts at k=1.
    assert(start == 1);
    // Tie-break: leaving_downward=true → we move DOWN below the chord.
    // Above-chord region is r_loop (containing C₂'s arcs); below is
    // r_outside.  Algorithm elects r_outside.
    assert(state.s2_region == r_outside &&
           "vertex-to-vertex tie-break: leaving_downward=true should "
           "elect the region below the chord (r_outside)");
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
// back DOWN past y=3 to v_end (y=1).  The horizontal ray from junction
// at y=3 crosses edge 1 (v_mid → v_end) mid-edge.  Two LEFT-side adj
// arcs at the crossing — one above (from v_mid_L down to crossing)
// and one below (from crossing down to v_end_L).
//
// (The null-length chord at v_mid is conceptually present per [C91 §2.1 tex 72]
// but is not modeled explicitly — the lambda processes only the matching
// chord at junction's y_tag, so other chords don't enter the matching
// branch.  This minimal S₂ exercises the mid-edge path directly.)
static void test_startup_mid_edge_tie_break() {
    auto C1 = make_C1();
    auto S1 = make_S1(C1);

    // C₂: junction at y=3, peak at y=5, end at y=1.  Edge 1 descends
    // from y=5 to y=1 crossing y=3 at x≈2.5.
    Polygon C2({{4,3,4}, {3,5,5}, {2,1,6}});

    Submap S2;
    std::size_t r_above = S2.add_node();
    std::size_t r_below = S2.add_node();

    // LEFT-side arc above the chord: from v_mid_L down edge 1 to the
    // mid-edge crossing.  Single-edge on edge 1.  Starts at v_mid_L
    // (the top companion of v_mid's null-length chord, conceptually),
    // ends at the chord crossing.
    Arc a{};
    a.first_edge = 1; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r_above;
    a.edge_count = 1;
    std::size_t ai_above_L = S2.add_arc(a);

    // LEFT-side arc below the chord: from mid-edge crossing down edge 1
    // to v_end_L.  Single-edge on edge 1.  Starts at the chord.
    a = {};
    a.first_edge = 1; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r_below;
    a.edge_count = 1;
    std::size_t ai_below_L = S2.add_arc(a);

    // RIGHT-side arc at v_junction_R (vertex endpoint at the junction).
    // Single-edge on edge 0.  RIGHT-side traversal starts at v_mid_R
    // and descends edge 0 to v_junction_R.
    a = {};
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r_above;
    a.edge_count = 1;
    std::size_t ai_right_at_junction = S2.add_arc(a);

    // Matching chord at y=3, y_tag=4 (junction's tag).
    // Left slot: mid-edge crossing on edge 1 (two adj arcs).
    // Right slot: v_junction_R on edge 0 (one adj arc).
    Chord c{};
    c.region[0] = r_above;
    c.region[1] = r_below;
    c.left_edge = 1; c.left_side = LEFT;
    c.right_edge = 0; c.right_side = RIGHT;
    c.y = 3.0; c.y_tag = 4;
    c.left_adj = {{ai_above_L, ai_below_L}, 2};
    c.right_adj = {{ai_right_at_junction}, 1};
    S2.add_chord(c);

    S2.start_arc = ai_above_L; S2.end_arc = ai_above_L;
    S2.start_vertex = 0; S2.end_vertex = 2;

    // a₀ shoots LEFT and hits the mid-edge crossing on edge 1 at x≈2.5.
    StartupOracle oracle1(&C1, -10.0);
    StartupOracle oracle2(&C2, 2.5);

    FusionState state;
    build_fusion_sequence(state, S1, C1);
    std::size_t start = fusion_startup(state, S1, C1, S2, C2,
                                        oracle1, oracle2);

    assert(start == 1);
    // Tie-break: leaving_downward=true → elect region BELOW chord.
    //
    // Trace through resolve_s2_region's mid-edge branch:
    //   - ai_above_L: single-edge on chord's edge.  arc_start_symbolic_y
    //     hits the polygon-vertex case at vertex(first_edge=1) = v_mid
    //     (y=5,tag=5) ≠ chord (y=3,tag=4) → ENDS at chord.  Previous vertex (LEFT,
    //     last_edge=1) = C₂.edge(1).start_idx = v_mid (y=5 > 3) →
    //     arc0_above = TRUE.
    //   - ai_below_L: single-edge.  arc_start_symbolic_y finds chord via
    //     c.left_adj.arcs[1] match → returns chord's (3,4) → STARTS at
    //     chord.  Next vertex (LEFT, first_edge=1) = C₂.edge(1).end_idx
    //     = v_end (y=1 < 3) → arc1_above = FALSE.
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
struct ForwardOracle : RayShootingOracle {
    const Polygon* C;
    double forward_dist;
    ForwardOracle(const Polygon* c, double d) : C(c), forward_dist(d) {}
    RayHit shoot(Point p, Side dir,
                 std::size_t /*arc_idx*/,
                 const Subarc& target) const override {
        RayHit h;
        h.hit = true;
        h.y = p.y;
        h.x = (dir == LEFT) ? (p.x - forward_dist) : (p.x + forward_dist);
        h.edge = target.first_edge;
        const auto& e = C->edge(target.first_edge);
        bool asc = symbolic_y_less(symbolic_y_of(C->vertex(e.start_idx)),
                                    symbolic_y_of(C->vertex(e.end_idx)));
        h.side = asc ? (dir == LEFT ? RIGHT : LEFT) : dir;
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

    Submap S2;
    S2.add_node();
    Arc a{};
    a.first_edge = 0; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    std::size_t ai0 = S2.add_arc(a);
    a.first_edge = 1; a.last_edge = 0;
    a.first_side = RIGHT; a.last_side = RIGHT;
    S2.add_arc(a);
    S2.start_arc = ai0; S2.end_arc = ai0;
    S2.start_vertex = 0; S2.end_vertex = 2;

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
//
// Adapts test_startup_mid_edge_tie_break's S₂ (which has a chord and
// 2 regions) and runs the full main loop.  Goal: exercise the case (ii)
// per-endpoint test code path — R has an incident exit chord, so the
// (ii) test iterates real candidates (even if none ultimately fire for
// this specific geometry, the predicate's loop body is reached).
static void test_fuse_main_loop_case_ii_smoke() {
    auto C1 = make_C1();
    auto S1 = make_S1(C1);

    // C₂ goes UP from junction (y=3) to v_mid (y=5) and back DOWN past
    // y=3 to v_end (y=1).  Edge 1 (v_mid→v_end) crosses y=3 at x≈2.5.
    Polygon C2({{4,3,4}, {3,5,5}, {2,1,6}});

    Submap S2;
    std::size_t r_above = S2.add_node();
    std::size_t r_below = S2.add_node();

    // S₂ chord at y=3 (junction's y): mid-edge crossing on edge 1 (LEFT)
    // and vertex endpoint at junction on edge 0 (RIGHT).
    Arc a{};
    a.first_edge = 1; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r_above;
    a.edge_count = 1;
    std::size_t ai_above_L = S2.add_arc(a);

    a = {};
    a.first_edge = 1; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r_below;
    a.edge_count = 1;
    std::size_t ai_below_L = S2.add_arc(a);

    a = {};
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r_above;
    a.edge_count = 1;
    std::size_t ai_right_at_junction = S2.add_arc(a);

    Chord c{};
    c.region[0] = r_above;
    c.region[1] = r_below;
    c.left_edge = 1; c.left_side = LEFT;
    c.right_edge = 0; c.right_side = RIGHT;
    c.y = 3.0; c.y_tag = 4;
    c.left_adj = {{ai_above_L, ai_below_L}, 2};
    c.right_adj = {{ai_right_at_junction}, 1};
    S2.add_chord(c);

    S2.start_arc = ai_above_L; S2.end_arc = ai_above_L;
    S2.start_vertex = 0; S2.end_vertex = 2;

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
    std::printf("All §3.1 tests passed.\n");
    return 0;
}
