/// tests/s3.1-tests.cpp — Tests for §3.1: Fusion of two submaps.

#include "merge/fusion.h"
#include "polygon/polygon.h"
#include "submap/submap.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════
//  Helpers
// ════════════════════════════════════════════════════════════════

/// C₁ = 5 vertices (edges 0–3), with one chord at vertex 2.
static Polygon make_C1() {
    return Polygon({{0,0,0}, {1,2,1}, {2,4,2}, {3,1,3}, {4,3,4}});
}

/// Build a 2-region submap of C₁ with one chord at vertex 2.
static Submap make_S1(const Polygon& C) {
    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1

    Arc a{};
    // LEFT arcs
    a.first_edge = 0; a.last_edge = 1;
    a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = 0; a.edge_count = 2;
    a.key_y = C.vertex(0).y; a.key_y_tag = 0;
    std::size_t ai0 = s.add_arc(a);

    a.first_edge = 1; a.last_edge = 3;
    a.region_node = 1; a.edge_count = 3;
    a.key_y = C.vertex(2).y; a.key_y_tag = 2;
    std::size_t ai1 = s.add_arc(a);

    // RIGHT arcs
    a = {};
    a.first_edge = 3; a.last_edge = 1;
    a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = 1; a.edge_count = 3;
    a.key_y = C.vertex(4).y; a.key_y_tag = 4;
    std::size_t ai2 = s.add_arc(a);

    a.first_edge = 1; a.last_edge = 0;
    a.region_node = 0; a.edge_count = 2;
    a.key_y = C.vertex(2).y; a.key_y_tag = 2;
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
    a.key_y_tag = 0;
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

/// Oracle that returns a fixed hit for any subarc on the RIGHT side.
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

/// Oracle that returns different x values per subarc side.
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

/// Oracle that always hits at a fixed x on the target's side.
struct StartupOracle : RayShootingOracle {
    double hit_x;
    explicit StartupOracle(double x) : hit_x(x) {}
    RayHit shoot(Point p, Side dir,
                 std::size_t /*arc_idx*/,
                 const Subarc& target) const override {
        RayHit h;
        h.hit = true;
        h.x = hit_x;
        h.y = p.y;
        h.edge = target.first_edge;
        // A horizontal ray always hits the face it approaches first.
        // Shooting LEFT → hits a LEFT-side point; RIGHT → RIGHT.
        // This ensures the discovered chord connects opposite sides
        // (§2.2 tex 96: exit chords cross the double boundary).
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
    a.key_y = C2.vertex(0).y; a.key_y_tag = 0;
    std::size_t ai0 = S2.add_arc(a);
    a.first_edge = 1; a.last_edge = 0;
    a.first_side = RIGHT; a.last_side = RIGHT;
    a.key_y = C2.vertex(2).y; a.key_y_tag = 2;
    std::size_t ai1 = S2.add_arc(a);
    // §2.4 (tex 144): end_arc = last LEFT arc (ai0).
    S2.start_arc = ai0; S2.end_arc = ai0;
    S2.start_vertex = 0; S2.end_vertex = 2;

    // a₀ at vertex 4 = (4,3). Edge 3 ascending → shoot LEFT → p.x ≈ 4.
    // Hits must be to the LEFT (x < 4).
    // Oracle for S₁: hits at x=-10 (far left).
    // Oracle for S₂: hits at x=3 (close left).
    // c₀ on ∂C₂ (closer) → Case 1.
    StartupOracle oracle1(-10.0);
    StartupOracle oracle2(3.0);

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
    a.key_y = C2.vertex(0).y; a.key_y_tag = 0;
    std::size_t ai0 = S2.add_arc(a);
    a.first_edge = 1; a.last_edge = 0;
    a.first_side = RIGHT; a.last_side = RIGHT;
    a.key_y = C2.vertex(2).y; a.key_y_tag = 2;
    std::size_t ai1 = S2.add_arc(a);
    // §2.4 (tex 144): end_arc = last LEFT arc (ai0).
    S2.start_arc = ai0; S2.end_arc = ai0;
    S2.start_vertex = 0; S2.end_vertex = 2;

    // a₀ at vertex 4 = (4,3). Edge 3 ascending → shoot LEFT → p.x ≈ 4.
    // Oracle for S₁: hits at x=3 (close left, on ∂C₁).
    // Oracle for S₂: hits at x=-10 (far left, on ∂C₂).
    // c₀ on ∂C₁ (closer) → Case 2.
    StartupOracle oracle1(3.0);
    StartupOracle oracle2(-10.0);

    FusionState state;
    build_fusion_sequence(state, S1, C1);

    std::size_t start = fusion_startup(state, S1, C1, S2, C2,
                                        oracle1, oracle2);

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
    // clockwise around ∂C".  Under §2.1 tex 66's induced orientation,
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
//  11. local_shoot — §2.1 double-boundary tie-break at same distance
// ════════════════════════════════════════════════════════════════

/// Oracle returning two arcs at the *same* x but different sides.
struct TieBreakRayShooter : RayShootingOracle {
    Side first_arc_side;                        ///< side of arc 0 hit
    Side second_arc_side;                       ///< side of arc 1 hit
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
    // [C91 §2.1 tex 72] (snake left/right) + §2.4 tex 142:
    // when two arcs deliver hits at exactly the same x (distance tie),
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
    // Both must return side=LEFT — the §2.1 tex 72 disambiguation rule.
    Point p{0.0, 2.0, 99};

    TieBreakRayShooter oracle(LEFT, RIGHT);
    RayHit h1 = local_shoot(p, RIGHT, 0, S1, C1, oracle);
    assert(h1.hit);
    assert(h1.x == 5.0);
    assert(h1.side == LEFT &&
           "§2.1 tex 72: struck-first face for shooting RIGHT toward "
           "an ascending-edge hit must be LEFT (the -x face)");

    TieBreakRayShooter oracle_rev(RIGHT, LEFT);
    RayHit h2 = local_shoot(p, RIGHT, 0, S1, C1, oracle_rev);
    assert(h2.hit);
    assert(h2.x == 5.0);
    assert(h2.side == LEFT &&
           "§2.1 tex 72: tie-break must select the same face regardless "
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
    a.key_y = C2.vertex(0).y; a.key_y_tag = 0;
    std::size_t ai0 = S2.add_arc(a);
    a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.key_y = C2.vertex(2).y; a.key_y_tag = 2;
    S2.add_arc(a);
    S2.start_arc = ai0; S2.end_arc = ai0;
    S2.start_vertex = 0; S2.end_vertex = 2;

    // Both oracles return hits at the SAME x (3.0), so d1 == d2.
    StartupOracle oracle1(3.0);
    StartupOracle oracle2(3.0);

    FusionState state;
    build_fusion_sequence(state, S1, C1);

    std::size_t start = fusion_startup(state, S1, C1, S2, C2, oracle1, oracle2);

    // tex 191 default: tie → Case 1 → main loop starts at k=1.
    assert(start == 1);
    assert(!state.chords.empty());

    std::printf("  [PASS] startup_d1_eq_d2_defaults_to_case1\n");
}

// ════════════════════════════════════════════════════════════════
//  13. build_fusion_sequence skips NLCs (§3.1 tex 179 + §2.2 tex 96)
// ════════════════════════════════════════════════════════════════

static void test_build_fusion_sequence_skips_nlcs() {
    // [C91 §3.1 tex 179]: "Let a₁, a₂, ..., aₘ be the canonical vertex
    // enumeration of S₁ ... exit chord endpoints in S₁."  Per §2.2 tex
    // 96, exit chords are distinguished from null-length chords; per §3.1
    // tex 224, NLCs are "carried over automatically" to the fused submap
    // and must NOT appear in the main-loop sequence.
    auto C = make_C1();                         // 5 vertices, 4 edges.

    // Synthetic submap with 1 exit chord + 1 NLC.  NLCs are at y-extrema
    // (vertex 2 = (2, 4) is a y-max in C1).
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();              // for exit chord
    std::size_t r_nlc = s.add_node();           // empty region of NLC

    Arc a{};
    a.first_edge = 0; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r0; a.edge_count = 2;
    a.key_y = C.vertex(0).y; a.key_y_tag = 0;
    std::size_t ai0 = s.add_arc(a);

    // NLC empty arc on LEFT at vertex 2 (y-max).
    a = {}; a.first_edge = 1; a.last_edge = 1; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r_nlc; a.edge_count = 0;    // null-length
    a.key_y = C.vertex(2).y; a.key_y_tag = 2;
    std::size_t ai_nlc = s.add_arc(a);

    a = {}; a.first_edge = 1; a.last_edge = 3; a.first_side = LEFT; a.last_side = LEFT;
    a.region_node = r1; a.edge_count = 3;
    a.key_y = C.vertex(2).y; a.key_y_tag = 22;  // distinct tag
    std::size_t ai1 = s.add_arc(a);

    a = {}; a.first_edge = 3; a.last_edge = 1; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r1; a.edge_count = 3;
    a.key_y = C.vertex(4).y; a.key_y_tag = 4;
    std::size_t ai2 = s.add_arc(a);

    a = {}; a.first_edge = 1; a.last_edge = 0; a.first_side = RIGHT; a.last_side = RIGHT;
    a.region_node = r0; a.edge_count = 2;
    a.key_y = C.vertex(2).y; a.key_y_tag = 22;
    std::size_t ai3 = s.add_arc(a);

    // Exit chord (non-NLC): r0 ↔ r1 at vertex 2's y level.
    Chord ec;
    ec.region[0] = r0; ec.region[1] = r1;
    ec.left_edge = 1; ec.right_edge = 1;
    ec.left_side = LEFT; ec.right_side = RIGHT;
    ec.y = C.vertex(2).y; ec.y_tag = 22;
    ec.left_adj = {{ai0, ai1}, 2};
    ec.right_adj = {{ai2, ai3}, 2};
    s.add_chord(ec);

    // NLC: r0 ↔ r_nlc at vertex 2 LEFT side.
    Chord nlc;
    nlc.region[0] = r0; nlc.region[1] = r_nlc;
    nlc.left_edge = 1; nlc.right_edge = 1;
    nlc.left_side = LEFT; nlc.right_side = LEFT;
    nlc.y = C.vertex(2).y; nlc.y_tag = 2;
    nlc.is_null_length = true;
    nlc.left_adj = {{ai0}, 1};
    nlc.right_adj = {{ai_nlc}, 1};
    s.add_chord(nlc);

    s.start_arc = ai0; s.end_arc = ai1;
    s.start_vertex = 0; s.end_vertex = 4;

    FusionState state;
    build_fusion_sequence(state, s, C);

    // Sequence = a₀ + 2 endpoints (exit chord only — NLC excluded) + a_{m+1}
    //          = 4 vertices.  If the NLC were included, we'd get 6.
    assert(state.sequence.size() == 4 &&
           "§3.1 tex 179: NLCs must not appear in the canonical "
           "vertex enumeration (paper distinguishes 'exit chord "
           "endpoints' from null-length chords per §2.2 tex 96)");

    // None of the non-companion vertices should reference the NLC.
    std::size_t nlc_idx = 1;                    // NLC is chord #1
    for (std::size_t i = 0; i < state.sequence.size(); ++i) {
        const auto& v = state.sequence[i];
        if (v.is_companion) continue;
        assert(v.chord_idx != nlc_idx);
    }

    std::printf("  [PASS] build_fusion_sequence_skips_nlcs\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::setbuf(stdout, nullptr);
    std::printf("§3.1 tests:\n");
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
    test_build_fusion_sequence_skips_nlcs();
    std::printf("All §3.1 tests passed.\n");
    return 0;
}
