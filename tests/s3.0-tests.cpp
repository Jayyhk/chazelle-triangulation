/// tests/s3.0-tests.cpp — Tests for §3.0: merge setup and preconditions.

#include "merge/merge.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════
//  Stub oracles for testing (§3.4 provides real implementations)
// ════════════════════════════════════════════════════════════════

struct StubRayShooter : RayShootingOracle {
    RayHit shoot(Point, Side, std::size_t, const Subarc&) const override {
        return {}; // no hit
    }
};

struct StubArcCutter : ArcCuttingOracle {
    std::vector<ArcPiece> cut(std::size_t, const Subarc&) const override {
        return {};
    }
};

static const StubRayShooter STUB_RAY;
static const StubArcCutter STUB_ARC;

/// Build a MergeInput with stub oracles.
static MergeInput make_input(const Polygon& C1, const Polygon& C2,
                              Submap& S1, Submap& S2,
                              std::size_t g1, std::size_t g2,
                              std::size_t g) {
    MergeInput in;
    in.C1 = &C1; in.C2 = &C2;
    in.S1 = &S1; in.S2 = &S2;
    in.gamma1 = g1; in.gamma2 = g2; in.gamma = g;
    in.ray_shooter_1 = &STUB_RAY;
    in.ray_shooter_2 = &STUB_RAY;
    in.arc_cutter_1  = &STUB_ARC;
    in.arc_cutter_2  = &STUB_ARC;
    return in;
}

// ════════════════════════════════════════════════════════════════
//  Helpers
// ════════════════════════════════════════════════════════════════

/// Two curves sharing vertex (5,5,2): C₁ = v0→v1→v2, C₂ = v2→v3→v4.
static Polygon make_C1() {
    return Polygon({{0,0,0}, {3,3,1}, {5,5,2}});
}
static Polygon make_C2() {
    return Polygon({{5,5,2}, {7,2,3}, {10,8,4}});
}

/// Single-region conformal submap (no chords → granular by default).
static Submap make_single_region_submap(const Polygon& poly) {
    Submap s;
    s.add_node();
    Arc a{};
    a.first_edge = 0;
    a.last_edge = poly.num_edges() - 1;
    a.first_side = LEFT;
    a.last_side = LEFT;
    a.region_node = 0;
    a.edge_count = poly.count_nonnull_edges(0, poly.num_edges() - 1);
    a.key_y = poly.vertex(0).y;
    a.key_y_tag = 0;
    std::size_t ai0 = s.add_arc(a);

    a.first_edge = poly.num_edges() - 1;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.key_y = poly.vertex(poly.num_vertices() - 1).y;
    a.key_y_tag = poly.num_vertices() - 1;
    std::size_t ai1 = s.add_arc(a);

    s.start_arc = ai0;
    s.end_arc = ai1;
    s.start_vertex = 0;
    s.end_vertex = poly.num_vertices() - 1;
    s.build_tree_decomposition();
    return s;
}

// ════════════════════════════════════════════════════════════════
//  1. Preconditions pass for valid input
// ════════════════════════════════════════════════════════════════

static void test_valid_preconditions() {
    auto C1 = make_C1();
    auto C2 = make_C2();
    auto S1 = make_single_region_submap(C1);
    auto S2 = make_single_region_submap(C2);

    auto in = make_input(C1, C2, S1, S2, 2, 2, 2);
    assert_merge_preconditions(in);

    std::printf("  [PASS] valid_preconditions\n");
}

// ════════════════════════════════════════════════════════════════
//  2. Merged curve C = C₁ ∪ C₂
// ════════════════════════════════════════════════════════════════

static void test_merged_curve() {
    auto C1 = make_C1();
    auto C2 = make_C2();
    auto S1 = make_single_region_submap(C1);
    auto S2 = make_single_region_submap(C2);

    auto in = make_input(C1, C2, S1, S2, 2, 2, 2);
    auto result = merge(in);

    // C = C₁ ∪ C₂: 3 + 3 - 1 (shared vertex) = 5 vertices.
    assert(result.C.num_vertices() == 5);
    assert(result.C.num_edges() == 4);

    // Verify vertex sequence: v0, v1, v2, v3, v4.
    assert(result.C.vertex(0).index == 0);
    assert(result.C.vertex(1).index == 1);
    assert(result.C.vertex(2).index == 2); // shared
    assert(result.C.vertex(3).index == 3);
    assert(result.C.vertex(4).index == 4);

    // Coordinates preserved.
    assert(result.C.vertex(0).x == 0.0 && result.C.vertex(0).y == 0.0);
    assert(result.C.vertex(2).x == 5.0 && result.C.vertex(2).y == 5.0);
    assert(result.C.vertex(4).x == 10.0 && result.C.vertex(4).y == 8.0);

    std::printf("  [PASS] merged_curve\n");
}

// ════════════════════════════════════════════════════════════════
//  3. γ₁ ≤ γ₂ constraint
// ════════════════════════════════════════════════════════════════

static void test_gamma_ordering() {
    auto C1 = make_C1();
    auto C2 = make_C2();
    auto S1 = make_single_region_submap(C1);
    auto S2 = make_single_region_submap(C2);

    // γ₁ = γ₂ (equality allowed).
    assert_merge_preconditions(make_input(C1, C2, S1, S2, 5, 5, 5));

    // γ₁ < γ₂ < γ.  γ₁ must still be ≥ region weight (= 2).
    assert_merge_preconditions(make_input(C1, C2, S1, S2, 2, 3, 10));

    // γ = γ₂ (equality allowed).
    assert_merge_preconditions(make_input(C1, C2, S1, S2, 2, 4, 4));

    std::printf("  [PASS] gamma_ordering\n");
}

// ════════════════════════════════════════════════════════════════
//  4. Shared vertex deduplication
// ════════════════════════════════════════════════════════════════

static void test_shared_vertex() {
    // Minimal curves: 2 vertices each, sharing one.
    Polygon C1({{0,0,0}, {1,1,1}});
    Polygon C2({{1,1,1}, {2,2,2}});
    auto S1 = make_single_region_submap(C1);
    auto S2 = make_single_region_submap(C2);

    auto in = make_input(C1, C2, S1, S2, 1, 1, 1);
    auto result = merge(in);

    // 2 + 2 - 1 = 3 vertices.
    assert(result.C.num_vertices() == 3);
    assert(result.C.num_edges() == 2);

    // Shared vertex appears exactly once.
    assert(result.C.vertex(1).x == 1.0);
    assert(result.C.vertex(1).y == 1.0);
    assert(result.C.vertex(1).index == 1);

    std::printf("  [PASS] shared_vertex\n");
}

// ════════════════════════════════════════════════════════════════
//  5. Edge count in merged curve
// ════════════════════════════════════════════════════════════════

static void test_merged_edges() {
    // C₁ has 4 vertices (3 edges), C₂ has 3 vertices (2 edges).
    // Shared vertex at (5,5,3).
    Polygon C1({{0,0,0}, {2,3,1}, {4,1,2}, {5,5,3}});
    Polygon C2({{5,5,3}, {7,2,4}, {9,6,5}});
    auto S1 = make_single_region_submap(C1);
    auto S2 = make_single_region_submap(C2);

    auto in = make_input(C1, C2, S1, S2, 3, 3, 3);
    auto result = merge(in);

    // 4 + 3 - 1 = 6 vertices, 5 edges.
    assert(result.C.num_vertices() == 6);
    assert(result.C.num_edges() == 5);

    // count_nonnull_edges on full range.
    assert(result.C.count_nonnull_edges(0, 4) == 5);

    std::printf("  [PASS] merged_edges\n");
}

// ════════════════════════════════════════════════════════════════
//  6. Conformality precondition
// ════════════════════════════════════════════════════════════════

static void test_conformality_required() {
    // Both submaps must be conformal.  Verify the check passes
    // for single-region submaps (trivially conformal, degree 0).
    auto C1 = make_C1();
    auto C2 = make_C2();
    auto S1 = make_single_region_submap(C1);
    auto S2 = make_single_region_submap(C2);

    assert(S1.is_conformal());
    assert(S2.is_conformal());

    assert_merge_preconditions(make_input(C1, C2, S1, S2, 2, 2, 2));

    std::printf("  [PASS] conformality_required\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::setbuf(stdout, nullptr);
    std::printf("§3.0 tests:\n");
    test_valid_preconditions();
    test_merged_curve();
    test_gamma_ordering();
    test_shared_vertex();
    test_merged_edges();
    test_conformality_required();
    std::printf("All §3.0 tests passed.\n");
    return 0;
}
