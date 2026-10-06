#include "algorithm/merge/oracle.h"
#include "algorithm/polygon/polygon.h"
#include "algorithm/submap/submap.h"
#include "algorithm/visibility/chain.h"
#include "algorithm/visibility/naive_visibility.h"
#include "algorithm/visibility/up_phase.h"
#include "support/random.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <vector>

using namespace chazelle;
using chazelle::test::DeterministicRandomGenerator;

namespace {

std::vector<Point> comb() {
    return {{0, 0, 0},     {2, 20, 1},   {4, 6, 2},   {6, 24, 3}, {8, 4, 4},
            {10, 22, 5},   {12, 5, 6},   {14, 26, 7}, {16, 2, 8}, {18, 23, 9},
            {20, 4.5, 10}, {22, 25, 11}, {24, 1, 12}};
}

std::vector<Point> zigzag18() {
    return {{0.0, 14.35, 0},   {1.1, 12.28, 1},   {2.4, 11.46, 2},   {3.6, 22.60, 3},
            {5.1, 37.35, 4},   {6.3, 6.76, 5},    {7.6, 38.54, 6},   {8.5, 9.75, 7},
            {10.2, 23.86, 8},  {11.4, 29.38, 9},  {12.9, 34.53, 10}, {14.1, 3.11, 11},
            {15.5, 18.42, 12}, {16.8, 27.90, 13}, {18.0, 8.02, 14},  {19.3, 31.77, 15},
            {20.7, 16.09, 16}, {22.0, 24.71, 17}};
}

void check_all_grades(const UpPhase& up) {
    for (std::size_t lam = 0; lam <= up.graded().maximum_grade(); ++lam) {
        for (std::size_t i = 0; i < up.graded().num_chains(lam); ++i) {
            const Polygon& c = up.graded().chain(lam, i);
            const Submap& s = up.chain_submap(lam, i);
            s.check_invariants(c);
            assert(s.is_conformal() && "[C91 §4.1 tex 327]: canonical ⟹ conformal");
            assert(s.is_granular(UpPhase::grade_granularity(lam), c) &&
                   "[C91 §4.1 tex 327/330]: canonical ⟹ "
                   "2^{⌈βλ⌉}-granular");
        }
    }
}

}

static void test_grade_gamma() {
    assert(UpPhase::grade_granularity(0) == 1);
    assert(UpPhase::grade_granularity(1) == 2);
    assert(UpPhase::grade_granularity(5) == 2);
    assert(UpPhase::grade_granularity(6) == 4);
    assert(UpPhase::grade_granularity(10) == 4);
    assert(UpPhase::grade_granularity(11) == 8);

    assert(UpPhase::NAIVE_MAX_GRADE == 1);
    assert(!(ceil_beta(1) < 1) && (ceil_beta(2) < 2));

    assert(canonical_granularity(1) == 1);
    assert(canonical_granularity(2) == 2);
    assert(canonical_granularity(32) == 2);
    assert(canonical_granularity(33) == 4);

    std::printf("  [PASS] grade_gamma\n");
}

static void test_naive_single_edge_vc() {
    Polygon curve({{0, 0, 0}, {1, 5, 1}});
    Submap submap = build_full_visibility_map(curve);
    submap.check_invariants(curve);
    assert(submap.is_conformal());

    std::size_t live_chords = 0;
    for (std::size_t ci = 0; ci < submap.num_chords(); ++ci) {
        const Chord& ch = submap.chord(ci);
        if (ch.dead)
            continue;
        ++live_chords;

        assert(ch.left_edge == 0 && ch.right_edge == 0);
        assert(ch.left_side != ch.right_side);
        assert(!ch.is_null_length);
        assert(ch.y_tag == 0 || ch.y_tag == 1);
    }
    assert(live_chords == 2 && "[C91 §2.1 tex 72] case 3: one chord per endpoint duplicate "
                               "pair of a single-edge curve");
    assert(submap.num_nodes() == 3 && "two wrap chords at distinct levels split the sphere into "
                                      "3 regions");

    std::printf("  [PASS] naive_single_edge_vc\n");
}

static void test_naive_small_chains() {
    Polygon curve({{0, 1, 0}, {1, 6, 1}, {2, 2, 2}});
    Submap submap = build_full_visibility_map(curve);
    submap.check_invariants(curve);
    assert(submap.is_conformal());
    assert(submap.is_semigranular(1) && "[C91 §2.3]: the full V(C) is 1-semigranular");

    bool has_null = false;
    for (std::size_t ci = 0; ci < submap.num_chords(); ++ci)
        if (!submap.chord(ci).dead && submap.chord(ci).is_null_length) {
            assert(submap.chord(ci).y_tag == 1 &&
                   "[C91 §2.1 tex 70/72]: the null chord sits at the "
                   "interior extremum's inside pair");
            has_null = true;
        }
    assert(has_null);

    Submap Sc = build_canonical_submap_naive(curve);
    Sc.check_invariants(curve);
    assert(Sc.is_conformal());
    assert(Sc.is_granular(2, curve));

    std::printf("  [PASS] naive_small_chains\n");
}

static void test_chain_grid() {
    std::vector<Point> in;
    in.reserve(7);
    for (std::size_t i = 0; i < 7; ++i)
        in.push_back({Exact(i), Exact((i * 3) % 7), i});
    GradedCurve g(in);
    assert(g.maximum_grade() == 3);
    assert(g.curve().num_vertices() == 9);
    assert(g.curve().count_nonnull_edges(0, 7) == 8);
    assert(g.curve().vertex(8).x == in.back().x && g.curve().vertex(8).y == in.back().y);

    for (std::size_t lam = 0; lam <= 3; ++lam) {
        assert(g.num_chains(lam) == (std::size_t{1} << (3 - lam)));
        for (std::size_t i = 0; i < g.num_chains(lam); ++i) {
            const Polygon& c = g.chain(lam, i);
            assert(c.num_vertices() == (std::size_t{1} << lam) + 1);
            assert(c.table_offset() == i * (std::size_t{1} << lam));
        }
    }
    assert(&g.chain(3, 0) == &g.curve());

    std::printf("  [PASS] chain_grid\n");
}

static void test_up_phase_comb() {
    UpPhase up(comb());
    assert(up.graded().maximum_grade() == 4);
    check_all_grades(up);
    std::printf("  [PASS] up_phase_comb (p=%zu)\n", up.graded().maximum_grade());
}

static void test_up_phase_zigzag() {
    UpPhase up(zigzag18());
    assert(up.graded().maximum_grade() == 5);
    check_all_grades(up);
    std::printf("  [PASS] up_phase_zigzag (p=%zu)\n", up.graded().maximum_grade());
}

static void test_canonical_portion() {
    UpPhase up(zigzag18());
    const std::size_t np = up.graded().curve().num_vertices();

    struct Range {
        std::size_t a, b;
    };
    const Range ranges[] = {{3, 9}, {4, 8}, {0, 5}, {1, 16}, {5, 21}, {13, 29}, {2, 32}, {20, 32}};
    for (const Range& r : ranges) {
        assert(r.b < np);
        auto pr = up.compute_canonical_portion(r.a, r.b);
        pr.submap.check_invariants(pr.curve);
        assert(pr.submap.is_conformal() && "[C91 Lemma 4.1 tex 347]: the portion's submap is "
                                           "canonical ⟹ conformal");
        assert(pr.curve.num_vertices() == r.b - r.a + 1);
        assert(pr.curve.table_offset() == r.a);

        assert(pr.submap.is_granular(canonical_granularity(r.b - r.a), pr.curve) &&
               "[C91 §4.1 tex 327]: canonical granularity of the "
               "portion");
    }

    std::printf("  [PASS] canonical_portion\n");
}

static void test_up_phase_cutter() {
    UpPhase up(zigzag18());
    const std::size_t lam = 3;
    const Polygon& curve = up.graded().chain(lam, 0);
    const Submap& submap = up.chain_submap(lam, 0);

    UpPhaseArcCutter cutter(up, submap, curve, lam);
    bool saw_smaller_granularity = false;
    for (std::size_t ai = 0; ai < submap.num_arcs(); ++ai) {
        if (submap.arc(ai).dead)
            continue;
        const Arc& a = submap.arc(ai);
        Subarc full;
        full.first_edge = a.first_edge;
        full.first_side = a.first_side;
        full.last_edge = a.last_edge;
        full.last_side = a.last_side;
        full.first_y = submap.arc_start_symbolic_y(ai, curve);
        full.last_y = submap.arc_end_symbolic_y(ai, curve);

        auto pieces = cutter.cut(ai, full);

        assert_cut_postconditions(curve, full, pieces.data(), pieces.size(),
                                  UpPhaseArcCutter::piece_count_bound(lam),
                                  UpPhaseArcCutter::piece_granularity_bound(lam));
        assert(pieces.size() <= UpPhaseArcCutter::piece_count_bound(lam) &&
               "[C91 §4.1 tex 343–346]: at most g(γ) = O(λ) pieces");
        for (const ArcPiece& p : pieces) {
            if (p.submap == nullptr)
                continue;
            assert(p.submap->is_granular(p.granularity, *p.curve));
            assert(p.submap->is_semigranular(UpPhaseArcCutter::piece_granularity_bound(lam)));
            assert(p.granularity <= UpPhaseArcCutter::piece_granularity_bound(lam) &&
                   "[C91 §4.1 tex 343–346]: piece granularity ≤ "
                   "2^{⌈β⌈βλ⌉⌉}");

            if (p.granularity == 1) {
                assert(!p.is_boundary_piece && p.curve->num_edges() == 1);
                assert(p.submap->num_live_chords() == 2);
                assert(!p.submap->is_granular(UpPhaseArcCutter::piece_granularity_bound(lam),
                                              *p.curve));
                saw_smaller_granularity = true;
            }
        }
    }
    assert(saw_smaller_granularity && "fixture must witness the uniform-h contract conflict");

    std::printf("  [PASS] up_phase_cutter\n");
}

static void test_double_identify_exhaustive() {
    Polygon curve(zigzag18());
    Submap submap = build_full_visibility_map(curve);
    submap.compact();
    submap.check_invariants(curve);

    auto reference = [&](std::size_t e, const SymbolicY& qy) {
        std::vector<std::size_t> out;
        for (std::size_t ai = 0; ai < submap.num_arcs(); ++ai) {
            if (submap.arc(ai).dead)
                continue;
            const Arc& a = submap.arc(ai);
            Subarc full;
            full.first_edge = a.first_edge;
            full.first_side = a.first_side;
            full.last_edge = a.last_edge;
            full.last_side = a.last_side;
            full.first_y = submap.arc_start_symbolic_y(ai, curve);
            full.last_y = submap.arc_end_symbolic_y(ai, curve);
            if (subarc_contains_point(full, curve, e, LEFT, qy, 0, curve.num_vertices() - 1) ||
                subarc_contains_point(full, curve, e, RIGHT, qy, 0, curve.num_vertices() - 1))
                out.push_back(ai);
        }
        return out;
    };

    for (std::size_t e = 0; e < curve.num_edges(); ++e) {
        SymbolicY lo = symbolic_y_of(curve.vertex(curve.edge(e).start_idx));
        SymbolicY hi = symbolic_y_of(curve.vertex(curve.edge(e).end_idx));
        if (symbolic_y_greater(lo, hi))
            std::swap(lo, hi);
        for (std::size_t v = 0; v < curve.num_vertices(); ++v) {
            SymbolicY qy = symbolic_y_of(curve.vertex(v));
            if (symbolic_y_less(qy, lo) || symbolic_y_greater(qy, hi))
                continue;

            auto res = submap.double_identify(e, qy, curve);
            auto ref = reference(e, qy);
            assert(res.count <= 6 && "[C91 §2.4 tex 144]: at most six arcs contain an edge");
            assert(res.count == ref.size() &&
                   "[C91 §2.4 tex 144]: double_identify must return EVERY "
                   "arc containing the queried edge (exhaustiveness vs the "
                   "full V(C))");
            for (std::size_t ai : ref) {
                bool found = false;
                for (std::size_t k = 0; k < res.count; ++k)
                    found = found || res.arcs[k] == ai;
                assert(found && "[C91 §2.4 tex 144]: reference arc missing from "
                                "double_identify's result");
            }
        }
    }

    std::printf("  [PASS] double_identify_exhaustive\n");
}

static void test_duplicate_raw_ys() {
    std::vector<Point> v;
    const Exact ys[] = {4, 6, 0, 7, 0, 1, 3, 7, 1, 7, 2, 5, 4, 5, 3, 2, 2, 2, 0, 7};
    Exact x = 0.0;
    for (std::size_t i = 0; i < 20; ++i) {
        x += 1.0 + 0.1 * Exact(i % 3);
        v.push_back({x, ys[i], i});
    }
    UpPhase up(v);
    assert(up.graded().maximum_grade() == 5);
    check_all_grades(up);

    const std::size_t np = up.graded().curve().num_vertices();
    struct Range {
        std::size_t a, b;
    };
    const Range ranges[] = {{0, 8}, {3, 9}, {5, 15}, {8, 16}, {2, 18}, {17, 32}, {24, 32}};
    for (const Range& r : ranges) {
        assert(r.b < np);
        auto pr = up.compute_canonical_portion(r.a, r.b);
        pr.submap.check_invariants(pr.curve);
        assert(pr.submap.is_conformal());
        assert(pr.submap.is_granular(canonical_granularity(r.b - r.a), pr.curve));
    }

    std::printf("  [PASS] duplicate_raw_ys (p=%zu)\n", up.graded().maximum_grade());
}

static void test_deep_grades_invariant_b() {
    DeterministicRandomGenerator rng(1);
    const std::size_t n = 60 + rng.next() % 200;
    assert(n == 248);
    std::vector<Point> v;
    Exact x = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        x += rng.uniform(0.5, 2.0);
        v.push_back({x, rng.uniform(0.0, 40.0), i});
    }
    UpPhase up(v);
    assert(up.graded().maximum_grade() == 8);
    check_all_grades(up);

    auto pr = up.compute_canonical_portion(208, 224);
    pr.submap.check_invariants(pr.curve);
    assert(pr.submap.is_conformal() && pr.submap.is_granular(canonical_granularity(16), pr.curve));
    auto pr2 = up.compute_canonical_portion(3, 200);
    pr2.submap.check_invariants(pr2.curve);
    assert(pr2.submap.is_conformal() &&
           pr2.submap.is_granular(canonical_granularity(197), pr2.curve));

    while (v.size() < 257) {
        x += 1.0;
        v.push_back({x, rng.uniform(0.0, 40.0), v.size()});
    }
    UpPhase historical(v);
    check_all_grades(historical);
    auto pinned = historical.compute_canonical_portion(208, 224);
    pinned.submap.check_invariants(pinned.curve);
    assert(pinned.submap.is_conformal() &&
           pinned.submap.is_granular(canonical_granularity(16), pinned.curve));

    std::printf("  [PASS] deep_grades_invariant_b (p=%zu)\n", up.graded().maximum_grade());
}

int main() {
    std::printf("[C91 §4.1 tests]:\n");

    UpPhase vertical_tail({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {2, 3, 3}});
    check_all_grades(vertical_tail);
    test_grade_gamma();
    test_naive_single_edge_vc();
    test_naive_small_chains();
    test_chain_grid();
    test_up_phase_comb();
    test_up_phase_zigzag();
    test_canonical_portion();
    test_up_phase_cutter();
    test_double_identify_exhaustive();
    test_duplicate_raw_ys();
    test_deep_grades_invariant_b();
    std::printf("[C91 §4.1 tests]: all passed\n");
    return 0;
}
