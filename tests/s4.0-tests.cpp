// tests/s4.0-tests.cpp — Tests for [C91 §4.0]: The Visibility
// Algorithm preamble (tex 314–323).
//
// Covers: β = 1/5 as an exact rational (tex 323), padding the input
// curve to n = 2^p + 1 with weight-free null-length edges (tex 316),
// the chain-in-grade-λ grid and its "obviously" facts (i)–(iii)
// (tex 316–321), view semantics over the one input table ([C91 §2.4
// tex 133]), O(1)-union y-extremes up the grid ([C91 §2 tex 47]), and
// the precondition asserts.

#include "visibility/chain.h"
#include "polygon/polygon.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <functional>
#include <vector>
#include <sys/prctl.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace chazelle;

namespace {

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

// A simple nonclosed staircase curve with n vertices: strictly
// increasing x, alternating y — every coordinate distinct, indices
// consecutive from 0.
std::vector<Point> staircase(std::size_t n) {
    std::vector<Point> v;
    v.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        v.push_back(Point{static_cast<double>(i),
                          (i % 2 == 0) ? static_cast<double>(i)
                                       : -static_cast<double>(i),
                          i});
    return v;
}

// The comb of tests/s3.2-tests.cpp / s3.4-tests.cpp (13 vertices).
std::vector<Point> comb() {
    return {{0,0,0}, {2,20,1}, {4,6,2}, {6,24,3}, {8,4,4},
            {10,22,5}, {12,5,6}, {14,26,7}, {16,2,8},
            {18,23,9}, {20,4.5,10}, {22,25,11}, {24,1,12}};
}

} // namespace

// ════════════════════════════════════════════════════════════════
//  β = 1/5 ([C91 §4 tex 323])
// ════════════════════════════════════════════════════════════════

static void test_beta() {
    // tex 323: "we set β = 1/5".
    assert(BETA_NUM == 1 && BETA_DEN == 5);

    // ⌈β·k⌉ in exact integer arithmetic.
    assert(ceil_beta(0) == 0);
    assert(ceil_beta(1) == 1);
    assert(ceil_beta(4) == 1);
    assert(ceil_beta(5) == 1);
    assert(ceil_beta(6) == 2);
    assert(ceil_beta(10) == 2);
    assert(ceil_beta(11) == 3);
    assert(ceil_beta(100) == 20);

    std::printf("  [PASS] beta\n");
}

// ════════════════════════════════════════════════════════════════
//  Padding ([C91 §4 tex 316])
// ════════════════════════════════════════════════════════════════

// n = 2^p + 1 after padding, p minimal; the original prefix is
// untouched and the tail duplicates the last vertex with consecutive
// SoS indices.
static void test_padding_sizes_and_content() {
    for (std::size_t n = 2; n <= 70; ++n) {
        const std::vector<Point> orig = staircase(n);
        const std::vector<Point> padded = pad_curve(orig);

        // n' = 2^p + 1 for some p (tex 316).
        const std::size_t np = padded.size();
        assert(np >= n);
        std::size_t p = 0;
        while ((std::size_t{1} << p) + 1 < np) ++p;
        assert((std::size_t{1} << p) + 1 == np &&
               "[C91 §4 tex 316]: padded n must be 2^p + 1");
        // Minimal p: "if necessary" — a curve that already fits in
        // 2^{p-1} + 1 vertices must not be padded past it.
        assert(p == 0 || (std::size_t{1} << (p - 1)) + 1 < n);

        // Original prefix untouched.
        for (std::size_t i = 0; i < n; ++i) {
            assert(padded[i].x == orig[i].x && padded[i].y == orig[i].y &&
                   padded[i].index == orig[i].index);
        }
        // Tail: coordinates of the last original vertex, consecutive
        // SoS indices ([C91 §2.4 tex 133]: the input table's indices
        // are consecutive).
        for (std::size_t i = n; i < np; ++i) {
            assert(padded[i].x == orig[n - 1].x &&
                   padded[i].y == orig[n - 1].y &&
                   "[C91 §4 tex 316]: padding duplicates the last vertex");
            assert(padded[i].index == padded[i - 1].index + 1);
        }
    }
    std::printf("  [PASS] padding_sizes_and_content\n");
}

// A curve with n = 2^p + 1 already needs no padding ("if necessary").
static void test_padding_no_op() {
    for (std::size_t n : {std::size_t{2}, std::size_t{3}, std::size_t{5},
                          std::size_t{9}, std::size_t{17}, std::size_t{65}}) {
        const std::vector<Point> orig = staircase(n);
        const std::vector<Point> padded = pad_curve(orig);
        assert(padded.size() == n &&
               "[C91 §4 tex 316]: pad only 'if necessary'");
        for (std::size_t i = 0; i < n; ++i)
            assert(padded[i].index == orig[i].index &&
                   padded[i].x == orig[i].x && padded[i].y == orig[i].y);
    }
    std::printf("  [PASS] padding_no_op\n");
}

// Padded edges are null length: they contribute no weight
// ([C91 §2.2 tex 106] counts nonnull-length edges only).
static void test_padding_null_edges() {
    const std::vector<Point> orig = comb();     // 13 vertices → 17
    const std::size_t n = orig.size();
    Polygon P(pad_curve(orig));
    assert(P.num_vertices() == 17);

    // All original edges are nonnull; every padded edge is null.
    assert(P.count_nonnull_edges(0, P.num_edges() - 1) == n - 1 &&
           "[C91 §2.2 tex 106]: padding must add zero weight");
    assert(P.count_nonnull_edges(n - 1, P.num_edges() - 1) == 0);

    // Padded vertices carry consecutive SoS tags, so the tail is
    // symbolically strictly descending ([C91 §2 tex 47]: larger tag ⟹
    // lower perturbed y) — interior padded vertices are never local
    // y-extrema.
    for (std::size_t i = n; i + 1 < P.num_vertices(); ++i)
        assert(!P.is_y_extremum(i) &&
               "[C91 §2 tex 47]: monotone null tail has no extrema");

    std::printf("  [PASS] padding_null_edges\n");
}

// Global y-extremes of the padded curve under SoS: the max is
// unaffected (padded copies sit symbolically below the vertex they
// copy); if the last vertex was the minimum, the new minimum is the
// final padded copy — same coordinates, largest tag.
static void test_padding_y_extremes() {
    // comb(): max at v7 (y=26), min at v0 (y=0) — neither is the last
    // vertex, so padding must leave both untouched.
    {
        Polygon P(pad_curve(comb()));
        assert(P.vertex(P.max_y_vertex()).index == 7);
        assert(P.vertex(P.min_y_vertex()).index == 0);
    }
    // A curve ENDING at its global minimum: the padded copies tie its
    // raw y, and the largest tag wins the symbolic minimum.
    {
        Polygon P(pad_curve({{0,5,0}, {1,3,1}, {2,4,2}, {3,0,3}}));
        assert(P.num_vertices() == 5);          // 4 → 2^2 + 1
        assert(P.vertex(P.max_y_vertex()).index == 0);
        assert(P.min_y_vertex() == P.num_vertices() - 1 &&
               P.vertex(P.min_y_vertex()).y == 0.0 &&
               "[C91 §2 tex 47]: largest tag at tied y is the symbolic "
               "min");
    }
    std::printf("  [PASS] padding_y_extremes\n");
}

// ════════════════════════════════════════════════════════════════
//  Chain grid facts (i)–(iii) ([C91 §4 tex 316–321])
// ════════════════════════════════════════════════════════════════

static void test_grid_facts() {
    const GradedCurve G(comb());                // 13 → 17 = 2^4 + 1
    assert(G.p() == 4);
    assert(G.curve().num_vertices() == 17);

    // (iii) p + 1 grades: 0, 1, ..., p (tex 320).
    assert(G.num_grades() == 5);

    for (std::size_t grade = 0; grade <= G.p(); ++grade) {
        // (ii) 2^{p−λ} chains in grade λ (tex 319).
        assert(G.num_chains(grade) ==
               (std::size_t{1} << (G.p() - grade)));

        for (std::size_t i = 0; i < G.num_chains(grade); ++i) {
            const Polygon& c = G.chain(grade, i);
            // (i) a grade-λ chain has 2^λ + 1 vertices (tex 318).
            assert(c.num_vertices() == (std::size_t{1} << grade) + 1);
            // tex 316: v_a..v_b with a − 1 = i·2^λ (1-based), i.e.
            // 0-based start i·2^λ — checked through the SoS index,
            // which equals the table position for this input.
            assert(c.vertex(0).index == i * (std::size_t{1} << grade) &&
                   "[C91 §4 tex 316]: a − 1 must be a multiple of 2^λ");
        }
    }
    std::printf("  [PASS] grid_facts\n");
}

// [C91 §2.4 tex 133]: every chain is a view into the one input table —
// vertex ADDRESSES coincide with the padded curve's.
static void test_chains_are_views() {
    const GradedCurve G(comb());
    const Polygon& P = G.curve();
    for (std::size_t grade = 0; grade <= G.p(); ++grade) {
        const std::size_t len = std::size_t{1} << grade;
        for (std::size_t i = 0; i < G.num_chains(grade); ++i) {
            const Polygon& c = G.chain(grade, i);
            for (std::size_t k = 0; k <= len; ++k)
                assert(&c.vertex(k) == &P.vertex(i * len + k) &&
                       "[C91 §2.4 tex 133]: the input table is never "
                       "copied");
        }
    }
    std::printf("  [PASS] chains_are_views\n");
}

// Consecutive chains of a grade share exactly one vertex of P, and a
// grade-λ chain is the union of its two grade-(λ−1) halves — the merge
// pattern of the up-phase ([C91 §3 tex 160]: C₁ ∩ C₂ is a vertex of P).
static void test_chain_adjacency_and_union() {
    const GradedCurve G(staircase(100));        // 100 → 129 = 2^7 + 1
    assert(G.p() == 7);

    for (std::size_t grade = 0; grade <= G.p(); ++grade) {
        for (std::size_t i = 0; i + 1 < G.num_chains(grade); ++i) {
            const Polygon& a = G.chain(grade, i);
            const Polygon& b = G.chain(grade, i + 1);
            assert(&a.vertex(a.num_vertices() - 1) == &b.vertex(0) &&
                   "[C91 §3 tex 160]: consecutive chains share one "
                   "vertex of P");
        }
    }
    for (std::size_t grade = 1; grade <= G.p(); ++grade) {
        const std::size_t half = std::size_t{1} << (grade - 1);
        for (std::size_t i = 0; i < G.num_chains(grade); ++i) {
            const Polygon& c = G.chain(grade, i);
            const Polygon& l = G.chain(grade - 1, 2 * i);
            const Polygon& r = G.chain(grade - 1, 2 * i + 1);
            assert(&c.vertex(0) == &l.vertex(0));
            assert(&c.vertex(half) == &r.vertex(0));
            assert(&c.vertex(2 * half) == &r.vertex(half));
        }
    }
    std::printf("  [PASS] chain_adjacency_and_union\n");
}

// Grade p holds the single chain covering all of P (tex 319: 2^0 = 1).
static void test_grade_p_is_whole_curve() {
    const GradedCurve G(staircase(33));         // already 2^5 + 1
    assert(G.p() == 5);
    assert(G.num_chains(G.p()) == 1);
    const Polygon& top = G.chain(G.p(), 0);
    assert(top.num_vertices() == G.curve().num_vertices());
    assert(&top.vertex(0) == &G.curve().vertex(0));
    std::printf("  [PASS] grade_p_is_whole_curve\n");
}

// The grid's y-extremes are combined in O(1) per union ([C91 §2 tex
// 47]); validate every chain's cached extremes against a direct scan.
static void test_chain_y_extremes() {
    const GradedCurve G(comb());
    for (std::size_t grade = 0; grade <= G.p(); ++grade) {
        for (std::size_t i = 0; i < G.num_chains(grade); ++i) {
            const Polygon& c = G.chain(grade, i);
            std::size_t mx = 0, mn = 0;
            for (std::size_t k = 1; k < c.num_vertices(); ++k) {
                if (point_y_above(c.vertex(k), c.vertex(mx))) mx = k;
                if (point_y_below(c.vertex(k), c.vertex(mn))) mn = k;
            }
            assert(c.max_y_vertex() == mx &&
                   "[C91 §2 tex 47]: union-combined max must equal scan");
            assert(c.min_y_vertex() == mn &&
                   "[C91 §2 tex 47]: union-combined min must equal scan");
        }
    }
    std::printf("  [PASS] chain_y_extremes\n");
}

// A larger instance: every fact of tex 317–321 across all grades.
static void test_large_instance() {
    const std::size_t n = 1000;                 // → 1025 = 2^10 + 1
    const GradedCurve G(staircase(n));
    assert(G.p() == 10);
    assert(G.curve().num_vertices() == 1025);
    assert(G.num_grades() == 11);

    std::size_t total_chains = 0;
    for (std::size_t grade = 0; grade <= G.p(); ++grade) {
        assert(G.num_chains(grade) == (std::size_t{1} << (G.p() - grade)));
        total_chains += G.num_chains(grade);
        // Spot-check first/last chain of the grade.
        const Polygon& first = G.chain(grade, 0);
        const Polygon& last = G.chain(grade, G.num_chains(grade) - 1);
        assert(first.num_vertices() == (std::size_t{1} << grade) + 1);
        assert(&first.vertex(0) == &G.curve().vertex(0));
        assert(&last.vertex(last.num_vertices() - 1) ==
               &G.curve().vertex(1024));
    }
    // Σ_{0 ≤ λ ≤ p} 2^{p−λ} = 2^{p+1} − 1 chains in the whole grid.
    assert(total_chains == (std::size_t{1} << (G.p() + 1)) - 1);
    std::printf("  [PASS] large_instance\n");
}

// Degenerate sizes: n = 2 (p = 0) and n = 3 (p = 1).
static void test_tiny_curves() {
    {
        const GradedCurve G(staircase(2));      // 2 = 2^0 + 1
        assert(G.p() == 0 && G.num_grades() == 1 && G.num_chains(0) == 1);
        assert(G.chain(0, 0).num_vertices() == 2);
    }
    {
        const GradedCurve G(staircase(3));      // 3 = 2^1 + 1
        assert(G.p() == 1 && G.num_grades() == 2);
        assert(G.num_chains(0) == 2 && G.num_chains(1) == 1);
    }
    std::printf("  [PASS] tiny_curves\n");
}

// ════════════════════════════════════════════════════════════════
//  Precondition asserts
// ════════════════════════════════════════════════════════════════

static void test_asserts_fire() {
    // Grades are 0..p (tex 320): out-of-range grade must abort.
    assert(assert_fires([] {
        const GradedCurve G(staircase(13));
        (void)G.num_chains(G.p() + 1);
    }));
    assert(assert_fires([] {
        const GradedCurve G(staircase(13));
        (void)G.chain(G.p() + 1, 0);
    }));
    // Only 2^{p−λ} chains in grade λ (tex 319).
    assert(assert_fires([] {
        const GradedCurve G(staircase(13));
        (void)G.chain(0, G.num_chains(0));
    }));
    // A curve needs ≥ 2 vertices ([C91 §2.1]).
    assert(assert_fires([] { (void)pad_curve(staircase(1)); }));
    assert(assert_fires([] { (void)GradedCurve(staircase(1)); }));

    std::printf("  [PASS] asserts_fire\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::printf("[C91 §4.0 tests]:\n");
    test_beta();
    test_padding_sizes_and_content();
    test_padding_no_op();
    test_padding_null_edges();
    test_padding_y_extremes();
    test_grid_facts();
    test_chains_are_views();
    test_chain_adjacency_and_union();
    test_grade_p_is_whole_curve();
    test_chain_y_extremes();
    test_large_instance();
    test_tiny_curves();
    test_asserts_fire();
    std::printf("All §4.0 tests passed.\n");
    return 0;
}
