// ════════════════════════════════════════════════════════════════
//  [C91 §4.1 tex 336–357] — The Up-Phase: canonical submaps per
//  grade, Lemma 4.1's partition/merge, and the tex 350–353 oracles.
// ════════════════════════════════════════════════════════════════

#include "../src/visibility/up_phase.h"
#include "../src/visibility/naive_visibility.h"
#include "../src/visibility/chain.h"
#include "../src/merge/oracle.h"
#include "../src/polygon/polygon.h"
#include "../src/submap/submap.h"

#include <cassert>
#include <cstdio>
#include <vector>

using namespace chazelle;

namespace {

// The comb of tests/s3.2-tests.cpp / s4.0-tests.cpp (13 vertices).
std::vector<Point> comb() {
    return {{0,0,0}, {2,20,1}, {4,6,2}, {6,24,3}, {8,4,4},
            {10,22,5}, {12,5,6}, {14,26,7}, {16,2,8},
            {18,23,9}, {20,4.5,10}, {22,25,11}, {24,1,12}};
}

// A fixed monotone-x polyline with distinct raw ys ([C91 §2 tex 45]),
// 18 vertices → padded to n = 2^5 + 1 = 33 (p = 5).
std::vector<Point> zigzag18() {
    return {{0.0,14.35,0},  {1.1,12.28,1},  {2.4,11.46,2},
            {3.6,22.60,3},  {5.1,37.35,4},  {6.3,6.76,5},
            {7.6,38.54,6},  {8.5,9.75,7},   {10.2,23.86,8},
            {11.4,29.38,9}, {12.9,34.53,10},{14.1,3.11,11},
            {15.5,18.42,12},{16.8,27.90,13},{18.0,8.02,14},
            {19.3,31.77,15},{20.7,16.09,16},{22.0,24.71,17}};
}

// Every chain of every processed grade carries a canonical submap:
// 2^{⌈βλ⌉}-granular, conformal, normal-form ([C91 §4.1 tex 338]).
void check_all_grades(const UpPhase& up) {
    for (std::size_t lam = 0; lam <= up.graded().p(); ++lam) {
        for (std::size_t i = 0; i < up.graded().num_chains(lam); ++i) {
            const Polygon& c = up.graded().chain(lam, i);
            const Submap& s = up.chain_submap(lam, i);
            s.check_invariants(c);
            assert(s.is_conformal() &&
                   "[C91 §4.1 tex 338]: canonical ⟹ conformal");
            assert(s.is_granular(UpPhase::grade_gamma(lam), c) &&
                   "[C91 §4.1 tex 338/342]: canonical ⟹ "
                   "2^{⌈βλ⌉}-granular");
        }
    }
}

} // namespace

// ════════════════════════════════════════════════════════════════
//  1. Grade granularities γ = 2^{⌈βλ⌉} ([C91 §4.1 tex 338/342])
// ════════════════════════════════════════════════════════════════

static void test_grade_gamma() {
    // β = 1/5 (tex 323): ⌈λ/5⌉ doublings.
    assert(UpPhase::grade_gamma(0) == 1);
    assert(UpPhase::grade_gamma(1) == 2);
    assert(UpPhase::grade_gamma(5) == 2);
    assert(UpPhase::grade_gamma(6) == 4);
    assert(UpPhase::grade_gamma(10) == 4);
    assert(UpPhase::grade_gamma(11) == 8);

    // [C91 §4.1 tex 345/352]: Lemma 4.1 needs ⌈βλ⌉ < λ — false exactly
    // for λ ≤ 1, the naively-computed early grades.
    assert(UpPhase::NAIVE_MAX_GRADE == 1);
    assert(!(ceil_beta(1) < 1) && (ceil_beta(2) < 2));

    // [C91 §4.1 tex 338]: canonical granularity of an m-edge curve is
    // 2^{⌈β⌈log m⌉⌉}.
    assert(canonical_granularity(1) == 1);
    assert(canonical_granularity(2) == 2);
    assert(canonical_granularity(32) == 2);   // ⌈log 32⌉ = 5, ⌈5/5⌉ = 1
    assert(canonical_granularity(33) == 4);   // ⌈log 33⌉ = 6, ⌈6/5⌉ = 2

    std::printf("  [PASS] grade_gamma\n");
}

// ════════════════════════════════════════════════════════════════
//  2. Full V(C) of a single-edge chain ([C91 §2.1 tex 68–72])
// ════════════════════════════════════════════════════════════════

static void test_naive_single_edge_vc() {
    // One edge, two curve endpoints.  tex 72 case 3: each endpoint has
    // two companions, one chord each; the four shots pair up into the
    // TWO wrap chords {(0,L),(0,R)} at the two endpoint levels,
    // splitting the plane into 3 regions.
    Polygon C({{0, 0, 0}, {1, 5, 1}});
    Submap S = build_full_visibility_map(C);
    S.check_invariants(C);
    assert(S.is_conformal());

    std::size_t live_chords = 0;
    for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
        const Chord& ch = S.chord(ci);
        if (ch.dead) continue;
        ++live_chords;
        // Both chords join the edge's two faces through infinity.
        assert(ch.left_edge == 0 && ch.right_edge == 0);
        assert(ch.left_side != ch.right_side);
        assert(!ch.is_null_length);
        assert(ch.y_tag == 0 || ch.y_tag == 1);
    }
    assert(live_chords == 2 &&
           "[C91 §2.1 tex 72] case 3: one chord per endpoint duplicate "
           "pair of a single-edge curve");
    assert(S.num_nodes() == 3 &&
           "two wrap chords at distinct levels split the sphere into "
           "3 regions");

    std::printf("  [PASS] naive_single_edge_vc\n");
}

// ════════════════════════════════════════════════════════════════
//  3. Full V(C) + canonical submap of small chains
//     ([C91 §2.1 tex 68–72] + [C91 §4.1 tex 338])
// ════════════════════════════════════════════════════════════════

static void test_naive_small_chains() {
    // 2-edge chain with an interior extremum: v1 is a local max, so
    // V(C) carries its inside-pair null chord ([C91 §2.1 tex 70/72]).
    // (Built under CHAZELLE_EXPENSIVE_ASSERTS, every chord is verified
    // mutually visible against naive_first_contact.)
    Polygon C({{0, 1, 0}, {1, 6, 1}, {2, 2, 2}});
    Submap S = build_full_visibility_map(C);
    S.check_invariants(C);
    assert(S.is_conformal());
    assert(S.is_semigranular(1) &&
           "[C91 §2.3]: the full V(C) is 1-semigranular");

    bool has_null = false;
    for (std::size_t ci = 0; ci < S.num_chords(); ++ci)
        if (!S.chord(ci).dead && S.chord(ci).is_null_length) {
            assert(S.chord(ci).y_tag == 1 &&
                   "[C91 §2.1 tex 70/72]: the null chord sits at the "
                   "interior extremum's inside pair");
            has_null = true;
        }
    assert(has_null);

    // The canonical submap enforces γ = canonical_granularity(2) = 2
    // on top ([C91 §4.1 tex 338] + §3.3).
    Submap Sc = build_canonical_submap_naive(C);
    Sc.check_invariants(C);
    assert(Sc.is_conformal());
    assert(Sc.is_granular(2, C));

    std::printf("  [PASS] naive_small_chains\n");
}

// ════════════════════════════════════════════════════════════════
//  4. Padding + chain grid ([C91 §4 tex 316–321])
// ════════════════════════════════════════════════════════════════

static void test_chain_grid() {
    // 7 input vertices → padded to n = 2^3 + 1 = 9 (tex 316), null
    // suffix cloning the last real vertex.
    std::vector<Point> in;
    for (std::size_t i = 0; i < 7; ++i)
        in.push_back({double(i), double((i * 3) % 7), i});
    GradedCurve g(in);
    assert(g.p() == 3);
    assert(g.curve().num_vertices() == 9);
    for (std::size_t v = 7; v < 9; ++v) {
        assert(g.curve().vertex(v).x == g.curve().vertex(6).x &&
               g.curve().vertex(v).y == g.curve().vertex(6).y &&
               "[C91 §4 tex 316]: null-length padding clones the last "
               "real vertex");
    }

    // tex 319 (ii): 2^{p−λ} chains in grade λ; grade-λ chain i spans
    // vertices [i·2^λ, (i+1)·2^λ] of P; consecutive chains share one
    // vertex (tex 320).
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

// ════════════════════════════════════════════════════════════════
//  5. Up-phase end-to-end on the comb ([C91 §4.1 tex 340–345])
// ════════════════════════════════════════════════════════════════

static void test_up_phase_comb() {
    UpPhase up(comb());                      // 13 → n = 17, p = 4
    assert(up.graded().p() == 4);
    check_all_grades(up);
    std::printf("  [PASS] up_phase_comb (p=%zu)\n", up.graded().p());
}

// ════════════════════════════════════════════════════════════════
//  6. Up-phase end-to-end on a general curve ([C91 §4.1 tex 340–345])
// ════════════════════════════════════════════════════════════════

static void test_up_phase_zigzag() {
    UpPhase up(zigzag18());                  // 18 → n = 33, p = 5
    assert(up.graded().p() == 5);
    check_all_grades(up);
    std::printf("  [PASS] up_phase_zigzag (p=%zu)\n", up.graded().p());
}

// ════════════════════════════════════════════════════════════════
//  7. Lemma 4.1 general portions ([C91 §4.1 tex 347–356])
// ════════════════════════════════════════════════════════════════

static void test_canonical_portion() {
    UpPhase up(zigzag18());
    const std::size_t np = up.graded().curve().num_vertices();   // 33

    // Portions v_a..v_b with 2^{λ−1} < b − a ≤ 2^λ for several λ and
    // alignments — including padding-spanning suffixes.
    struct Range { std::size_t a, b; };
    const Range ranges[] = {{3, 9},  {4, 8},   {0, 5},  {1, 16},
                            {5, 21}, {13, 29}, {2, 32}, {20, 32}};
    for (const Range& r : ranges) {
        assert(r.b < np);
        auto pr = up.compute_canonical_portion(r.a, r.b);
        pr.submap.check_invariants(pr.curve);
        assert(pr.submap.is_conformal() &&
               "[C91 Lemma 4.1 tex 347]: the portion's submap is "
               "canonical ⟹ conformal");
        assert(pr.curve.num_vertices() == r.b - r.a + 1);
        assert(pr.curve.table_offset() == r.a);
        // λ = ⌈log(b − a)⌉; canonical ⟹ 2^{⌈βλ⌉}-granular (tex 338).
        assert(pr.submap.is_granular(
                   canonical_granularity(r.b - r.a), pr.curve) &&
               "[C91 §4.1 tex 338]: canonical granularity of the "
               "portion");
    }

    std::printf("  [PASS] canonical_portion\n");
}

// ════════════════════════════════════════════════════════════════
//  8. The tex 352–353 arc-cutter bounds ([C91 §3.0(ii) tex 170])
// ════════════════════════════════════════════════════════════════

static void test_up_phase_cutter() {
    UpPhase up(zigzag18());
    const std::size_t lam = 3;               // a mid grade (> naive)
    const Polygon& C = up.graded().chain(lam, 0);
    const Submap& S = up.chain_submap(lam, 0);

    UpPhaseArcCutter cutter(up, S, C, lam);
    for (std::size_t ai = 0; ai < S.num_arcs(); ++ai) {
        if (S.arc(ai).dead) continue;
        const Arc& a = S.arc(ai);
        Subarc full;
        full.first_edge = a.first_edge;
        full.first_side = a.first_side;
        full.last_edge = a.last_edge;
        full.last_side = a.last_side;
        full.first_y = S.arc_start_symbolic_y(ai, C);
        full.last_y = S.arc_end_symbolic_y(ai, C);

        auto pieces = cutter.cut(ai, full);
        // [C91 §4.1 tex 353]: g(γ) pieces, each an h(γ)-granular
        // normal-form conformal submap with h(γ) ≤ 2^{⌈β⌈βλ⌉⌉}.
        // (assert_cut_postconditions re-validates each piece under
        // CHAZELLE_EXPENSIVE_ASSERTS inside cut().)
        assert(pieces.size() <= UpPhaseArcCutter::g_bound(lam) &&
               "[C91 §4.1 tex 353]: at most g(γ) = O(λ) pieces");
        for (const ArcPiece& p : pieces) {
            if (p.submap == nullptr) continue;   // boundary segment
            assert(p.granularity <= UpPhaseArcCutter::h_bound(lam) &&
                   "[C91 §4.1 tex 353]: piece granularity ≤ "
                   "2^{⌈β⌈βλ⌉⌉}");
        }
    }

    std::printf("  [PASS] up_phase_cutter\n");
}

// ════════════════════════════════════════════════════════════════
//  9. double_identify exhaustiveness against the full V(C)
//     ([C91 §2.4 tex 144]) — closes the deferral that waited on §4's
//     symbolic reference map.
// ════════════════════════════════════════════════════════════════

static void test_double_identify_exhaustive() {
    Polygon C(zigzag18());
    Submap S = build_full_visibility_map(C);
    S.compact();
    S.check_invariants(C);

    // Position-exact reference: the live arcs whose endpoint-precise
    // span ([C91 §3.0(i) tex 169]) contains either face of the query
    // edge at the query level.
    auto reference = [&](std::size_t e, SymbolicY qy) {
        std::vector<std::size_t> out;
        for (std::size_t ai = 0; ai < S.num_arcs(); ++ai) {
            if (S.arc(ai).dead) continue;
            const Arc& a = S.arc(ai);
            Subarc full;
            full.first_edge = a.first_edge;
            full.first_side = a.first_side;
            full.last_edge = a.last_edge;
            full.last_side = a.last_side;
            full.first_y = S.arc_start_symbolic_y(ai, C);
            full.last_y = S.arc_end_symbolic_y(ai, C);
            if (subarc_contains_point(full, C, e, LEFT, qy, 0,
                                      C.num_vertices() - 1) ||
                subarc_contains_point(full, C, e, RIGHT, qy, 0,
                                      C.num_vertices() - 1))
                out.push_back(ai);
        }
        return out;
    };

    for (std::size_t e = 0; e < C.num_edges(); ++e) {
        // Query at every vertex level within the edge's span — real
        // SoS tags ([C91 §2 tex 47]; double_identify requires one),
        // covering arc boundaries and interiors alike.
        const double y0 = C.vertex(C.edge(e).start_idx).y;
        const double y1 = C.vertex(C.edge(e).end_idx).y;
        const double lo = y0 < y1 ? y0 : y1;
        const double hi = y0 < y1 ? y1 : y0;
        for (std::size_t v = 0; v < C.num_vertices(); ++v) {
        SymbolicY qy = symbolic_y_of(C.vertex(v));
        if (qy.y < lo || qy.y > hi) continue;

        auto res = S.double_identify(e, qy, C);
        auto ref = reference(e, qy);
        assert(res.count <= 6 &&
               "[C91 §2.4 tex 144]: at most six arcs contain an edge");
        assert(res.count == ref.size() &&
               "[C91 §2.4 tex 144]: double_identify must return EVERY "
               "arc containing the queried edge (exhaustiveness vs the "
               "full V(C))");
        for (std::size_t ai : ref) {
            bool found = false;
            for (std::size_t k = 0; k < res.count; ++k)
                found = found || res.arcs[k] == ai;
            assert(found &&
                   "[C91 §2.4 tex 144]: reference arc missing from "
                   "double_identify's result");
        }
        }
    }

    std::printf("  [PASS] double_identify_exhaustive\n");
}

// ════════════════════════════════════════════════════════════════
//  10. Duplicate raw y-coordinates ([C91 §2 tex 45]: "we can easily
//      get around this assumption by applying the symbolic
//      perturbation techniques of [10] and [31]") — the y-only SoS
//      lift, exercised end-to-end.
// ════════════════════════════════════════════════════════════════

static void test_duplicate_raw_ys() {
    // A curve drawing ys from {0..7}: raw-tied vertices across the
    // curve, horizontal edges (consecutive equal ys), raw-tied wrap
    // chords, and padding clusters raw-tied with real vertices.  The
    // perturbed level order is the tag order ([C91 §2 tex 47]; larger
    // tag ⟹ lower), and every raw coincidence resolves through the
    // perturbed x-offsets (polygon.h::perturbed_x_offset).
    std::vector<Point> v;
    const double ys[] = {4, 6, 0, 7, 0, 1, 3, 7, 1, 7, 2, 5, 4,
                         5, 3, 2, 2, 2, 0, 7};
    double x = 0.0;
    for (std::size_t i = 0; i < 20; ++i) {
        x += 1.0 + 0.1 * double(i % 3);
        v.push_back({x, ys[i], i});
    }
    UpPhase up(v);                            // 20 → n = 33, p = 5
    assert(up.graded().p() == 5);
    check_all_grades(up);

    // General portions, including raw-tie-heavy and padding-adjacent
    // windows.
    const std::size_t np = up.graded().curve().num_vertices();
    struct Range { std::size_t a, b; };
    const Range ranges[] = {{0, 8}, {3, 9}, {5, 15}, {8, 16},
                            {2, 18}, {17, 32}, {24, 32}};
    for (const Range& r : ranges) {
        assert(r.b < np);
        auto pr = up.compute_canonical_portion(r.a, r.b);
        pr.submap.check_invariants(pr.curve);
        assert(pr.submap.is_conformal());
        assert(pr.submap.is_granular(
            canonical_granularity(r.b - r.a), pr.curve));
    }

    std::printf("  [PASS] duplicate_raw_ys (p=%zu)\n", up.graded().p());
}

// ════════════════════════════════════════════════════════════════
//  11. Scale pin — deep merge trees (p = 8) and the invariant-B tie
//      election ([C91 §3.1 tex 195]).
// ════════════════════════════════════════════════════════════════

static void test_deep_grades_invariant_b() {
    // Deterministic 248-vertex curve (xorshift, seed 1) that once
    // produced an invalid carried-over wrap chord: at the [208,224]
    // merge, a case (i) hit landed exactly on an S₂ chord endpoint,
    // p's sight rode that chord, and without invariant (B)'s tie
    // election ("the one that we enter as we locally leave p
    // clockwise", tex 195) the walk kept the wrong region and never
    // re-validated a stale wrap chord (tex 224).
    struct Rng {
        unsigned long long s;
        explicit Rng(unsigned long long seed)
            : s(seed ^ 0x9e3779b97f4a7c15ull) {}
        unsigned long long next() {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s;
        }
        double uniform(double lo, double hi) {
            return lo + (hi - lo) *
                            double(next() % 1000000ull) / 1000000.0;
        }
    };
    Rng rng(1);
    const std::size_t n = 60 + rng.next() % 200;    // 248
    assert(n == 248);
    std::vector<Point> v;
    double x = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        x += rng.uniform(0.5, 2.0);
        v.push_back({x, rng.uniform(0.0, 40.0), i});
    }
    UpPhase up(v);                                  // n → 257, p = 8
    assert(up.graded().p() == 8);
    check_all_grades(up);

    // The historically failing portion, plus a broad one.
    auto pr = up.compute_canonical_portion(208, 224);
    pr.submap.check_invariants(pr.curve);
    assert(pr.submap.is_conformal() &&
           pr.submap.is_granular(canonical_granularity(16), pr.curve));
    auto pr2 = up.compute_canonical_portion(3, 200);
    pr2.submap.check_invariants(pr2.curve);
    assert(pr2.submap.is_conformal() &&
           pr2.submap.is_granular(canonical_granularity(197),
                                  pr2.curve));

    std::printf("  [PASS] deep_grades_invariant_b (p=%zu)\n",
                up.graded().p());
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::printf("[C91 §4.1 tests]:\n");
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
