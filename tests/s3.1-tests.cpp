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

    auto seq = build_fusion_sequence(S1, C1);

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

    auto seq = build_fusion_sequence(s, C);

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

    auto seq = build_fusion_sequence(S1, C1);

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

    auto seq = build_fusion_sequence(S1, C1);

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
    std::vector<std::size_t> arcs;
    collect_region_arcs(S1, 0, arcs);
    assert(arcs.size() == 2);

    // Region 1 has 2 arcs (LEFT ai1 and RIGHT ai2).
    arcs.clear();
    collect_region_arcs(S1, 1, arcs);
    assert(arcs.size() == 2);

    std::printf("  [PASS] collect_region_arcs\n");
}

// ════════════════════════════════════════════════════════════════
//  6. local_shoot with stub oracle
// ════════════════════════════════════════════════════════════════

/// Oracle that returns a fixed hit for any subarc on the RIGHT side.
struct FixedRayShooter : RayShootingOracle {
    RayHit shoot(Point p, Side /*direction*/,
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
    std::printf("All §3.1 tests passed.\n");
    return 0;
}
