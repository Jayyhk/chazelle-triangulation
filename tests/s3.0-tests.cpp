// tests/s3.0-tests.cpp — Tests for [C91 §3.0]: merge setup and preconditions.

#include "merge/merge.h"
#include "merge/ray_shooting.h"

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

// ════════════════════════════════════════════════════════════════
//  Stub oracles for testing ([C91 §3.4] provides real implementations)
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

// Real [C91 §3.4] ray-shooting oracles for the tests that drive an
// actual merge (stage 1 shoots for real since the §3.4 wiring).
struct OracleRig {
    SubmapRayShooter ray1, ray2;
    OracleRig(const Submap& S1, const Polygon& C1, std::size_t g1,
              const Submap& S2, const Polygon& C2, std::size_t g2)
        : ray1(S1, C1, g1), ray2(S2, C2, g2) {}
};

// Build a MergeInput with stub oracles (precondition-only tests).
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
    // [C91 §3.0(ii) tex 170]: nominal arc-cutter bounds for the stubs.
    in.g_gamma1 = 1; in.g_gamma2 = 1;
    in.h_gamma1 = 1; in.h_gamma2 = 1;
    return in;
}

// ════════════════════════════════════════════════════════════════
//  Helpers
// ════════════════════════════════════════════════════════════════

// Two curves sharing vertex (5,5,2): C₁ = v0→v1→v2, C₂ = v2→v3→v4 —
// subchain views of the one input table of P ([C91 §2.4 tex 133]:
// "The input table is read-only: it is never to be modified or even
// copied"; [C91 §3 tex 160]: C₁ ∩ C₂ is a vertex of P).
static const Polygon& input_P() {
    static Polygon P({{0,0,0}, {3,3,1}, {5,5,2}, {7,2,3}, {10,8,4}});
    return P;
}
static Polygon make_C1() { return input_P().subchain(0, 3); }
static Polygon make_C2() { return input_P().subchain(2, 3); }

// Single-region conformal submap (no chords → granular by default):
// one region bounded by the single closed arc — all of ∂C, one
// arc-structure stored cut at C's start turnaround
// ([C91 §2.4 tex 142/138]).
static Submap make_single_region_submap(const Polygon& poly) {
    Submap s;
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = poly.num_vertices() - 1;
    Arc a{};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count =
        2 * poly.count_nonnull_edges(0, poly.num_edges() - 1);
    std::size_t ai0 = s.add_arc(a);
    assert(s.start_arc == ai0 && s.end_arc == ai0);
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

    auto in = make_input(C1, C2, S1, S2, 4, 4, 4);
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

    OracleRig rig(S1, C1, 4, S2, C2, 4);
    auto in = make_input(C1, C2, S1, S2, 4, 4, 4);
    in.ray_shooter_1 = &rig.ray1;
    in.ray_shooter_2 = &rig.ray2;
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

    // γ₁ < γ₂ < γ.  γ₁ must still be ≥ region weight (= 4,
    // [C91 §2.2 tex 106]: both sides of the closed arc's 2 edges).
    assert_merge_preconditions(make_input(C1, C2, S1, S2, 4, 5, 10));

    // γ = γ₂ (equality allowed).
    assert_merge_preconditions(make_input(C1, C2, S1, S2, 4, 5, 5));

    std::printf("  [PASS] gamma_ordering\n");
}

// ════════════════════════════════════════════════════════════════
//  4. Shared vertex deduplication
// ════════════════════════════════════════════════════════════════

static void test_shared_vertex() {
    // Minimal curves: 2 vertices each, sharing one ([C91 §2.4 tex 133]:
    // subchain views of one input table).
    Polygon P({{0,0,0}, {1,1,1}, {2,2,2}});
    Polygon C1 = P.subchain(0, 2);
    Polygon C2 = P.subchain(1, 2);
    auto S1 = make_single_region_submap(C1);
    auto S2 = make_single_region_submap(C2);

    OracleRig rig(S1, C1, 2, S2, C2, 2);
    auto in = make_input(C1, C2, S1, S2, 2, 2, 2);
    in.ray_shooter_1 = &rig.ray1;
    in.ray_shooter_2 = &rig.ray2;
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
    // Shared vertex at (5,5,3) ([C91 §2.4 tex 133]: subchain views).
    Polygon P({{0,0,0}, {2,3,1}, {4,1,2}, {5,5,3}, {7,2,4}, {9,6,5}});
    Polygon C1 = P.subchain(0, 4);
    Polygon C2 = P.subchain(3, 3);
    auto S1 = make_single_region_submap(C1);
    auto S2 = make_single_region_submap(C2);

    OracleRig rig(S1, C1, 6, S2, C2, 6);
    auto in = make_input(C1, C2, S1, S2, 6, 6, 6);
    in.ray_shooter_1 = &rig.ray1;
    in.ray_shooter_2 = &rig.ray2;
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

    assert_merge_preconditions(make_input(C1, C2, S1, S2, 4, 4, 4));

    std::printf("  [PASS] conformality_required\n");
}

// ════════════════════════════════════════════════════════════════
//  7. assert_cut_postconditions — positive case
// ════════════════════════════════════════════════════════════════

// Minimal valid 2-vertex submap for use as a piece's V(ᾱⱼ).
static Polygon make_segment_curve() {
    return Polygon({{0,0,100}, {1,1,101}});
}
static Submap make_segment_submap(const Polygon& C) {
    Submap s;
    s.add_node();
    s.start_vertex = 0; s.end_vertex = 1;
    // The single closed arc ([C91 §2.4 tex 142/138]).
    Arc a{};
    a.first_edge = 0; a.last_edge = 0; a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2 * C.count_nonnull_edges(0, 0);
    std::size_t ai0 = s.add_arc(a);
    assert(s.start_arc == ai0 && s.end_arc == ai0);
    s.build_tree_decomposition();
    return s;
}

static void test_cut_postconditions_valid() {
    // [C91 §3.0(ii) tex 170]: a valid cut() output (single piece, single
    // edge, vertex-to-vertex with a granular submap) must pass.
    Subarc target{0, LEFT, 0, LEFT};
    Polygon C = make_segment_curve();
    Submap S = make_segment_submap(C);

    ArcPiece p;
    p.subarc = target;
    p.submap = &S;
    p.curve = &C;
    p.is_boundary_piece = false;

    assert_cut_postconditions(target, &p, 1, /*max_pieces=*/4, /*h_gamma=*/1);

    std::printf("  [PASS] cut_postconditions_valid\n");
}

// ════════════════════════════════════════════════════════════════
//  8. assert_cut_postconditions — death tests for each (i)(ii)(iii)
// ════════════════════════════════════════════════════════════════

static void test_cut_count_zero_fires() {
    // [C91 §3.0(ii) tex 170]: subdivides into AT MOST g(γᵢ) pieces; zero
    // is impossible (a subarc cannot be subdivided into nothing).
    assert(assert_fires([]{
        Subarc t{0, LEFT, 0, LEFT};
        ArcPiece p{};
        assert_cut_postconditions(t, &p, /*count=*/0, 4, 1);
    }));
    std::printf("  [PASS] cut_count_zero_fires\n");
}

static void test_cut_count_exceeds_bound_fires() {
    // [C91 §3.0(ii) tex 170]: count must not exceed g(γᵢ).
    assert(assert_fires([]{
        Subarc t{0, LEFT, 0, LEFT};
        Polygon C = make_segment_curve();
        Submap S = make_segment_submap(C);
        ArcPiece pieces[5];
        for (auto& p : pieces) {
            p.subarc = {0, LEFT, 0, LEFT};
            p.submap = &S; p.curve = &C; p.is_boundary_piece = false;
        }
        assert_cut_postconditions(t, pieces, 5, /*max_pieces=*/4, 1);
    }));
    std::printf("  [PASS] cut_count_exceeds_bound_fires\n");
}

static void test_cut_first_endpoint_mismatch_fires() {
    // [C91 §3.0(ii) tex 170]: subdivision identity α' = α₁ ∪ ... ∪ αₖ —
    // first piece must start at α'.first.
    assert(assert_fires([]{
        Subarc target{0, LEFT, 1, LEFT};
        Polygon C = make_segment_curve();
        Submap S = make_segment_submap(C);
        ArcPiece p;
        p.subarc = {1, LEFT, 1, LEFT};         // ← does not start at target.first (edge 0)
        p.submap = &S; p.curve = &C; p.is_boundary_piece = false;
        assert_cut_postconditions(target, &p, 1, 4, 1);
    }));
    std::printf("  [PASS] cut_first_endpoint_mismatch_fires\n");
}

static void test_cut_last_endpoint_mismatch_fires() {
    // [C91 §3.0(ii) tex 170]: last piece must end at α'.last.
    assert(assert_fires([]{
        Subarc target{0, LEFT, 1, LEFT};
        Polygon C = make_segment_curve();
        Submap S = make_segment_submap(C);
        ArcPiece p;
        p.subarc = {0, LEFT, 0, LEFT};         // ← does not end at target.last (edge 1)
        p.submap = &S; p.curve = &C; p.is_boundary_piece = false;
        assert_cut_postconditions(target, &p, 1, 4, 1);
    }));
    std::printf("  [PASS] cut_last_endpoint_mismatch_fires\n");
}

static void test_cut_double_back_fires() {
    // [C91 §3.0(ii)(2) tex 170]: each piece on one side of C — no double-backing.
    assert(assert_fires([]{
        Subarc target{0, LEFT, 0, RIGHT};       // target itself doubles back
        Polygon C = make_segment_curve();
        Submap S = make_segment_submap(C);
        ArcPiece p;
        p.subarc = {0, LEFT, 0, RIGHT};         // ← first_side != last_side
        p.submap = &S; p.curve = &C; p.is_boundary_piece = false;
        assert_cut_postconditions(target, &p, 1, 4, 1);
    }));
    std::printf("  [PASS] cut_double_back_fires\n");
}

static void test_cut_non_clockwise_left_fires() {
    // [C91 §3.0(ii)(1) tex 170]: clockwise traversal — LEFT side ascends.
    // To exercise the LEFT-clockwise per-piece check (and not the prior
    // subdivision-identity check), target and piece share the same
    // (bad) endpoints: piece matches target so identity passes, then
    // the LEFT-side `first_edge ≤ last_edge` assertion fires.
    assert(assert_fires([]{
        Subarc target{1, LEFT, 0, LEFT};        // matches piece below
        Polygon C = make_segment_curve();
        Submap S = make_segment_submap(C);
        ArcPiece p;
        p.subarc = {1, LEFT, 0, LEFT};          // ← LEFT but first_edge > last_edge
        p.submap = &S; p.curve = &C; p.is_boundary_piece = false;
        assert_cut_postconditions(target, &p, 1, 4, 1);
    }));
    std::printf("  [PASS] cut_non_clockwise_left_fires\n");
}

static void test_cut_interior_boundary_piece_fires() {
    // [C91 §3.0(ii)(3) tex 170]: only first / last pieces may be boundary
    // pieces (ᾱ₁ and ᾱ₂ correspond to α''s two endpoints).
    assert(assert_fires([]{
        Subarc target{0, LEFT, 0, LEFT};
        Polygon C = make_segment_curve();
        Submap S = make_segment_submap(C);
        ArcPiece pieces[3];
        for (auto& p : pieces) {
            p.subarc = {0, LEFT, 0, LEFT};
            p.submap = &S; p.curve = &C; p.is_boundary_piece = false;
        }
        pieces[1].is_boundary_piece = true;     // ← interior boundary piece (j=1, count=3)
        pieces[1].submap = nullptr; pieces[1].curve = nullptr;
        assert_cut_postconditions(target, pieces, 3, 4, 1);
    }));
    std::printf("  [PASS] cut_interior_boundary_piece_fires\n");
}

static void test_cut_non_boundary_null_submap_fires() {
    // [C91 §3.0(ii)(3) tex 170]: non-boundary pieces require an h(γᵢ)-
    // granular conformal submap of V(ᾱⱼ).
    assert(assert_fires([]{
        Subarc target{0, LEFT, 0, LEFT};
        ArcPiece p;
        p.subarc = {0, LEFT, 0, LEFT};
        p.submap = nullptr;                     // ← non-boundary requires submap
        p.curve = nullptr;
        p.is_boundary_piece = false;
        assert_cut_postconditions(target, &p, 1, 4, 1);
    }));
    std::printf("  [PASS] cut_non_boundary_null_submap_fires\n");
}

static void test_cut_boundary_multi_edge_fires() {
    // [C91 §3.0(ii)(3) tex 170]: boundary pieces are single-edge pieces.
    assert(assert_fires([]{
        Subarc target{0, LEFT, 1, LEFT};
        ArcPiece pieces[2];
        pieces[0].subarc = {0, LEFT, 1, LEFT};  // ← boundary spans 2 edges
        pieces[0].is_boundary_piece = true;
        pieces[0].submap = nullptr; pieces[0].curve = nullptr;
        pieces[1].subarc = {1, LEFT, 1, LEFT};
        pieces[1].is_boundary_piece = true;
        pieces[1].submap = nullptr; pieces[1].curve = nullptr;
        assert_cut_postconditions(target, pieces, 2, 4, 1);
    }));
    std::printf("  [PASS] cut_boundary_multi_edge_fires\n");
}

// ════════════════════════════════════════════════════════════════
//  9. assert_merge_preconditions — death tests
// ════════════════════════════════════════════════════════════════

static void test_merge_preconds_null_curve_fires() {
    assert(assert_fires([]{
        auto C2 = make_C2();
        auto S2 = make_single_region_submap(C2);
        Polygon C1stub({{0,0,0},{1,1,1}});
        auto S1 = make_single_region_submap(C1stub);
        MergeInput in;
        in.C1 = nullptr;                        // ← null C₁
        in.C2 = &C2; in.S1 = &S1; in.S2 = &S2;
        in.gamma1 = 1; in.gamma2 = 1; in.gamma = 1;
        in.ray_shooter_1 = &STUB_RAY; in.ray_shooter_2 = &STUB_RAY;
        in.arc_cutter_1 = &STUB_ARC; in.arc_cutter_2 = &STUB_ARC;
        assert_merge_preconditions(in);
    }));
    std::printf("  [PASS] merge_preconds_null_curve_fires\n");
}

static void test_merge_preconds_unshared_junction_fires() {
    // [C91 §3 tex 160]: C₁ ∩ C₂ must be a vertex (last of C₁ = first of C₂).
    assert(assert_fires([]{
        Polygon C1({{0,0,0}, {1,1,1}});
        Polygon C2({{2,2,2}, {3,3,3}});         // ← C₂'s first ≠ C₁'s last
        auto S1 = make_single_region_submap(C1);
        auto S2 = make_single_region_submap(C2);
        assert_merge_preconditions(make_input(C1, C2, S1, S2, 1, 1, 1));
    }));
    std::printf("  [PASS] merge_preconds_unshared_junction_fires\n");
}

static void test_merge_preconds_gamma_order_fires() {
    // [C91 §3 tex 160]: γ₁ ≤ γ₂.
    assert(assert_fires([]{
        auto C1 = make_C1(); auto C2 = make_C2();
        auto S1 = make_single_region_submap(C1);
        auto S2 = make_single_region_submap(C2);
        assert_merge_preconditions(make_input(C1, C2, S1, S2, /*g1=*/5, /*g2=*/3, 5));
    }));
    std::printf("  [PASS] merge_preconds_gamma_order_fires\n");
}

static void test_merge_preconds_gamma_target_low_fires() {
    // [C91 §3 tex 160]: γ ≥ γ₂.
    assert(assert_fires([]{
        auto C1 = make_C1(); auto C2 = make_C2();
        auto S1 = make_single_region_submap(C1);
        auto S2 = make_single_region_submap(C2);
        assert_merge_preconditions(make_input(C1, C2, S1, S2, 2, 5, /*g=*/3));
    }));
    std::printf("  [PASS] merge_preconds_gamma_target_low_fires\n");
}

static void test_merge_preconds_null_oracle_fires() {
    // [C91 §3.0 tex 166-170]: ray-shooter and arc-cutter required for each Sᵢ.
    assert(assert_fires([]{
        auto C1 = make_C1(); auto C2 = make_C2();
        auto S1 = make_single_region_submap(C1);
        auto S2 = make_single_region_submap(C2);
        auto in = make_input(C1, C2, S1, S2, 2, 2, 2);
        in.ray_shooter_1 = nullptr;             // ← missing oracle
        assert_merge_preconditions(in);
    }));
    std::printf("  [PASS] merge_preconds_null_oracle_fires\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::setbuf(stdout, nullptr);
    std::printf("[C91 §3.0 tests]:\n");
    test_valid_preconditions();
    test_merged_curve();
    test_gamma_ordering();
    test_shared_vertex();
    test_merged_edges();
    test_conformality_required();
    test_cut_postconditions_valid();
    test_cut_count_zero_fires();
    test_cut_count_exceeds_bound_fires();
    test_cut_first_endpoint_mismatch_fires();
    test_cut_last_endpoint_mismatch_fires();
    test_cut_double_back_fires();
    test_cut_non_clockwise_left_fires();
    test_cut_interior_boundary_piece_fires();
    test_cut_non_boundary_null_submap_fires();
    test_cut_boundary_multi_edge_fires();
    test_merge_preconds_null_curve_fires();
    test_merge_preconds_unshared_junction_fires();
    test_merge_preconds_gamma_order_fires();
    test_merge_preconds_gamma_target_low_fires();
    test_merge_preconds_null_oracle_fires();
    std::printf("All §3.0 tests passed.\n");
    return 0;
}
