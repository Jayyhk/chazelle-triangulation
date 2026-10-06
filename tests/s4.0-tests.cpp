#include "algorithm/polygon/polygon.h"
#include "algorithm/visibility/chain.h"
#include "support/assertions.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <vector>

using namespace chazelle;

using chazelle::test::require_assertion_abort;

namespace {

std::vector<Point> staircase(std::size_t n) {
    std::vector<Point> v;
    v.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        v.push_back(Point{static_cast<Exact>(i),
                          (i % 2 == 0) ? static_cast<Exact>(i) : -static_cast<Exact>(i), i});
    return v;
}

std::vector<Point> comb() {
    return {{0, 0, 0},     {2, 20, 1},   {4, 6, 2},   {6, 24, 3}, {8, 4, 4},
            {10, 22, 5},   {12, 5, 6},   {14, 26, 7}, {16, 2, 8}, {18, 23, 9},
            {20, 4.5, 10}, {22, 25, 11}, {24, 1, 12}};
}

}

static void test_beta() {
    assert(BETA_NUM == 1 && BETA_DEN == 5);

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

static void test_padding_sizes_and_content() {
    for (std::size_t n = 2; n <= 70; ++n) {
        const std::vector<Point> orig = staircase(n);
        const std::vector<Point> padded = pad_curve(orig);
        const GradedCurve graded(orig);
        const auto original_positions = graded.original_vertex_positions();
        assert(original_positions.size() == n);

        const std::size_t np = padded.size();
        assert(np >= n);
        std::size_t p = 0;
        while ((std::size_t{1} << p) + 1 < np)
            ++p;
        assert((std::size_t{1} << p) + 1 == np && "[C91 §4 tex 316]: padded n must be 2^p + 1");

        assert(p == 0 || (std::size_t{1} << (p - 1)) + 1 < n);

        const std::size_t added = np - n;
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t j = i + std::min(i, added);
            assert(original_positions[i] == j && graded.curve().vertex(j).x == orig[i].x &&
                   graded.curve().vertex(j).y == orig[i].y &&
                   "[C91 §4 tex 316]: subdivision retains the original vertex positions");
            assert(padded[j].x == orig[i].x && padded[j].y == orig[i].y);
            if (i < added) {
                assert(padded[j + 1].x == exact_midpoint(orig[i].x, orig[i + 1].x));
                assert(padded[j + 1].y == exact_midpoint(orig[i].y, orig[i + 1].y));
            }
        }
        for (std::size_t i = 0; i < np; ++i)
            assert(padded[i].index == i);
        assert(padded.front().x == orig.front().x && padded.back().x == orig.back().x &&
               padded.back().y == orig.back().y);
    }
    std::printf("  [PASS] padding_sizes_and_content\n");
}

static void test_padding_no_op() {
    for (std::size_t n : {std::size_t{2}, std::size_t{3}, std::size_t{5}, std::size_t{9},
                          std::size_t{17}, std::size_t{65}}) {
        const std::vector<Point> orig = staircase(n);
        const std::vector<Point> padded = pad_curve(orig);
        assert(padded.size() == n && "[C91 §4 tex 316]: pad only 'if necessary'");
        for (std::size_t i = 0; i < n; ++i)
            assert(padded[i].index == orig[i].index && padded[i].x == orig[i].x &&
                   padded[i].y == orig[i].y);
    }
    std::printf("  [PASS] padding_no_op\n");
}

static void test_padding_nonnull_edges() {
    Polygon input_curve(pad_curve(comb()));
    assert(input_curve.num_vertices() == 17);
    assert(input_curve.count_nonnull_edges(0, input_curve.num_edges() - 1) ==
           input_curve.num_edges());
    for (std::size_t i = 0; i < 4; ++i)
        assert(!input_curve.is_y_extremum(2 * i + 1));
    std::printf("  [PASS] padding_nonnull_edges\n");
}

static void test_padding_y_extremes() {
    {
        Polygon input_curve(pad_curve(comb()));
        assert(input_curve.vertex(input_curve.max_y_vertex()).x == 14 &&
               input_curve.vertex(input_curve.max_y_vertex()).y == 26);
        assert(input_curve.vertex(input_curve.min_y_vertex()).index == 0);
    }

    {
        Polygon input_curve(pad_curve({{0, 5, 0}, {1, 3, 1}, {2, 4, 2}, {3, 0, 3}}));
        assert(input_curve.num_vertices() == 5);
        assert(input_curve.vertex(input_curve.max_y_vertex()).index == 0);
        assert(input_curve.min_y_vertex() == input_curve.num_vertices() - 1 &&
               input_curve.vertex(input_curve.min_y_vertex()).y == 0.0 &&
               "[C91 §4 tex 316]: padding preserves the endpoint");
    }
    std::printf("  [PASS] padding_y_extremes\n");
}

static void test_grid_facts() {
    const GradedCurve G(comb());
    assert(G.maximum_grade() == 4);
    assert(G.curve().num_vertices() == 17);

    assert(G.num_grades() == 5);

    for (std::size_t grade = 0; grade <= G.maximum_grade(); ++grade) {
        assert(G.num_chains(grade) == (std::size_t{1} << (G.maximum_grade() - grade)));

        for (std::size_t i = 0; i < G.num_chains(grade); ++i) {
            const Polygon& c = G.chain(grade, i);

            assert(c.num_vertices() == (std::size_t{1} << grade) + 1);

            assert(c.vertex(0).index == i * (std::size_t{1} << grade) &&
                   "[C91 §4 tex 316]: a − 1 must be a multiple of 2^λ");
        }
    }
    std::printf("  [PASS] grid_facts\n");
}

static void test_chains_are_views() {
    const GradedCurve G(comb());
    const Polygon& input_curve = G.curve();
    for (std::size_t grade = 0; grade <= G.maximum_grade(); ++grade) {
        const std::size_t len = std::size_t{1} << grade;
        for (std::size_t i = 0; i < G.num_chains(grade); ++i) {
            const Polygon& c = G.chain(grade, i);
            for (std::size_t k = 0; k <= len; ++k)
                assert(&c.vertex(k) == &input_curve.vertex(i * len + k) &&
                       "[C91 §2.4 tex 133]: the input table is never "
                       "copied");
        }
    }
    std::printf("  [PASS] chains_are_views\n");
}

static void test_chain_adjacency_and_union() {
    const GradedCurve G(staircase(100));
    assert(G.maximum_grade() == 7);

    for (std::size_t grade = 0; grade <= G.maximum_grade(); ++grade) {
        for (std::size_t i = 0; i + 1 < G.num_chains(grade); ++i) {
            const Polygon& a = G.chain(grade, i);
            const Polygon& b = G.chain(grade, i + 1);
            assert(&a.vertex(a.num_vertices() - 1) == &b.vertex(0) &&
                   "[C91 §3 tex 160]: consecutive chains share one "
                   "vertex of P");
        }
    }
    for (std::size_t grade = 1; grade <= G.maximum_grade(); ++grade) {
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

static void test_grade_p_is_whole_curve() {
    const GradedCurve G(staircase(33));
    assert(G.maximum_grade() == 5);
    assert(G.num_chains(G.maximum_grade()) == 1);
    const Polygon& top = G.chain(G.maximum_grade(), 0);
    assert(top.num_vertices() == G.curve().num_vertices());
    assert(&top.vertex(0) == &G.curve().vertex(0));
    std::printf("  [PASS] grade_p_is_whole_curve\n");
}

static void test_chain_y_extremes() {
    const GradedCurve G(comb());
    for (std::size_t grade = 0; grade <= G.maximum_grade(); ++grade) {
        for (std::size_t i = 0; i < G.num_chains(grade); ++i) {
            const Polygon& c = G.chain(grade, i);
            std::size_t mx = 0, mn = 0;
            for (std::size_t k = 1; k < c.num_vertices(); ++k) {
                if (point_y_above(c.vertex(k), c.vertex(mx)))
                    mx = k;
                if (point_y_below(c.vertex(k), c.vertex(mn)))
                    mn = k;
            }
            assert(c.max_y_vertex() == mx && "[C91 §2 tex 47]: union-combined max must equal scan");
            assert(c.min_y_vertex() == mn && "[C91 §2 tex 47]: union-combined min must equal scan");
        }
    }
    std::printf("  [PASS] chain_y_extremes\n");
}

static void test_large_instance() {
    const std::size_t n = 1000;
    const GradedCurve G(staircase(n));
    assert(G.maximum_grade() == 10);
    assert(G.curve().num_vertices() == 1025);
    assert(G.num_grades() == 11);

    std::size_t total_chains = 0;
    for (std::size_t grade = 0; grade <= G.maximum_grade(); ++grade) {
        assert(G.num_chains(grade) == (std::size_t{1} << (G.maximum_grade() - grade)));
        total_chains += G.num_chains(grade);

        const Polygon& first = G.chain(grade, 0);
        const Polygon& last = G.chain(grade, G.num_chains(grade) - 1);
        assert(first.num_vertices() == (std::size_t{1} << grade) + 1);
        assert(&first.vertex(0) == &G.curve().vertex(0));
        assert(&last.vertex(last.num_vertices() - 1) == &G.curve().vertex(1024));
    }

    assert(total_chains == (std::size_t{1} << (G.maximum_grade() + 1)) - 1);
    std::printf("  [PASS] large_instance\n");
}

static void test_tiny_curves() {
    {
        const GradedCurve G(staircase(2));
        assert(G.maximum_grade() == 0 && G.num_grades() == 1 && G.num_chains(0) == 1);
        assert(G.chain(0, 0).num_vertices() == 2);
    }
    {
        const GradedCurve G(staircase(3));
        assert(G.maximum_grade() == 1 && G.num_grades() == 2);
        assert(G.num_chains(0) == 2 && G.num_chains(1) == 1);
    }
    std::printf("  [PASS] tiny_curves\n");
}

static void test_asserts_fire() {
    require_assertion_abort([] {
        const GradedCurve G(staircase(13));
        (void)G.num_chains(G.maximum_grade() + 1);
    });
    require_assertion_abort([] {
        const GradedCurve G(staircase(13));
        (void)G.chain(G.maximum_grade() + 1, 0);
    });

    require_assertion_abort([] {
        const GradedCurve G(staircase(13));
        (void)G.chain(0, G.num_chains(0));
    });

    require_assertion_abort([] { (void)pad_curve(staircase(1)); });
    require_assertion_abort([] { (void)GradedCurve(staircase(1)); });

    std::printf("  [PASS] asserts_fire\n");
}

int main() {
    std::printf("[C91 §4.0 tests]:\n");
    test_beta();
    test_padding_sizes_and_content();
    test_padding_no_op();
    test_padding_nonnull_edges();
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
