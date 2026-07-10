// tests/s3-e2e-tests.cpp — End-to-end tests for [C91 §3]: Merging Two
// Submaps (§3.0 oracle contracts, §3.1 fusion, §3.2 conformality,
// §3.3 granularity, §3.4 ray-shooting), driven with the REAL
// production oracles ([C91 §4.1 tex 341]) and validated against the
// naive reference shooter ([C91 §2.1 tex 68–72]).
//
// The §2 e2e suite validates submap surgery against brute force; this
// suite validates the whole MERGE pipeline the same way: every chord
// of every merged submap must be a genuine visible pair with respect
// to the merged curve ([C91 §2.2 tex 90]), the output must be a
// canonical submap ([C91 §4.1 tex 338]), and the Lemma 3.6 structure
// built over it must agree with the naive shooter ray-for-ray.

#include "../src/visibility/up_phase.h"
#include "../src/visibility/naive_visibility.h"
#include "../src/merge/merge.h"
#include "../src/merge/granularity.h"
#include "../src/merge/ray_shooting.h"
#include "../src/merge/fusion.h"
#include "../src/polygon/polygon.h"
#include "../src/submap/submap.h"

#include <cassert>
#include <cstdio>
#include <vector>

using namespace chazelle;

namespace {

// Deterministic xorshift (the stress harness's generator).
struct Rng {
    unsigned long long s;
    explicit Rng(unsigned long long seed)
        : s(seed ^ 0x9e3779b97f4a7c15ull) {}
    unsigned long long next() {
        s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s;
    }
    double uniform(double lo, double hi) {
        return lo + (hi - lo) * double(next() % 1000000ull) / 1000000.0;
    }
};

// Monotone-x polyline (always simple).  dup_ys draws raw ys from
// {0..7} — raw ties across the curve, horizontal edges, and padding
// clusters raw-tied with real vertices ([C91 §2 tex 45/47]).
std::vector<Point> random_curve(Rng& rng, std::size_t n, bool dup_ys) {
    std::vector<Point> v;
    double x = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        x += rng.uniform(0.5, 2.0);
        double y = dup_ys ? double(rng.next() % 8)
                          : rng.uniform(0.0, 40.0);
        v.push_back({x, y, i});
    }
    return v;
}

// [C91 §2.2 tex 90 / §2.1 tex 74]: every non-null chord of S joins a
// mutually visible pair with respect to C — validated with the naive
// reference shooter from BOTH endpoints (unconditionally, not gated
// behind CHAZELLE_EXPENSIVE_ASSERTS).
void assert_chords_are_visible_pairs(const Submap& S, const Polygon& C) {
    for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
        const Chord& ch = S.chord(ci);
        if (ch.dead || ch.is_null_length) continue;
        struct End { std::size_t e; Side s; };
        const End ends[2] = {{ch.left_edge, ch.left_side},
                             {ch.right_edge, ch.right_side}};
        for (int k = 0; k < 2; ++k) {
            const End& a = ends[k];
            const End& b = ends[1 - k];
            const double xa = edge_x_at_y(C, a.e, ch.symbolic_y());
            const double xb = edge_x_at_y(C, b.e, ch.symbolic_y());
            Point p{xa, ch.symbolic_y().y, ch.symbolic_y().tag};
            RayHit h = naive_first_contact(
                C, p, ch.symbolic_y(),
                shooting_direction(a.e, a.s, C), /*source_edge=*/a.e);
            // The first contact is the partner — same position; the
            // edge label may be the adjacent one when the partner sits
            // at a non-extremum vertex ([C91 §3.0(i) tex 169]).
            [[maybe_unused]] const bool label_ok =
                h.edge == b.e ||
                (symbolic_y_equal(
                     ch.symbolic_y(),
                     symbolic_y_of(C.vertex(std::max(h.edge, b.e)))) &&
                 (h.edge + 1 == b.e || b.e + 1 == h.edge));
            assert(h.hit && h.x == xb && label_ok &&
                   "[C91 §2.2 tex 90]: merged-submap chord endpoints "
                   "must see each other with respect to C");
        }
    }
}

// [C91 §4.1 tex 338] canonical-form postconditions + invariants.
void assert_canonical(const Submap& S, const Polygon& C,
                      std::size_t gamma) {
    S.check_invariants(C);
    assert(S.is_conformal() &&
           "[C91 §3.2 Lemma 3.4]: merge output is conformal");
    assert(S.is_granular(gamma, C) &&
           "[C91 §3.3 Lemma 3.5]: merge output is γ-granular");
    assert(!S.tree_decomposition().empty() &&
           "[C91 §2.4(iv)]: normal form carries the tree decomposition");
}

} // namespace

// ════════════════════════════════════════════════════════════════
//  1. Merge pipeline vs the naive reference, across random curves
//     and portions ([C91 §3.0–3.4 + §4.1 Lemma 4.1]) — both y-modes.
// ════════════════════════════════════════════════════════════════

static void test_merged_portions_vs_naive(bool dup_ys) {
    for (unsigned seed = 1; seed <= 4; ++seed) {
        Rng rng(seed);
        const std::size_t n = 20 + rng.next() % 60;
        UpPhase up(random_curve(rng, n, dup_ys));
        const Polygon& P = up.graded().curve();

        // Every chain's canonical submap (built through §3 merges for
        // grades ≥ 2) has honest chords and canonical form.
        for (std::size_t lam = 0; lam <= up.graded().p(); ++lam)
            for (std::size_t i = 0; i < up.graded().num_chains(lam);
                 ++i) {
                const Polygon& c = up.graded().chain(lam, i);
                const Submap& s = up.chain_submap(lam, i);
                assert_canonical(s, c, UpPhase::grade_gamma(lam));
                assert_chords_are_visible_pairs(s, c);
            }

        // General portions (Lemma 4.1): each drives a fresh cascade of
        // §3 merges with the production oracles.
        const std::size_t np = P.num_vertices();
        for (int t = 0; t < 3; ++t) {
            std::size_t a = rng.next() % (np - 4);
            std::size_t len = 4 + rng.next() % (np - a - 4);
            if (len < 3) continue;
            auto pr = up.compute_canonical_portion(a, a + len);
            std::size_t lam = 0;
            while ((std::size_t{1} << lam) < len) ++lam;
            assert_canonical(pr.submap, pr.curve,
                             UpPhase::grade_gamma(lam));
            assert_chords_are_visible_pairs(pr.submap, pr.curve);
        }
    }
    std::printf("  [PASS] merged_portions_vs_naive (dup_ys=%d)\n",
                (int)dup_ys);
}

// ════════════════════════════════════════════════════════════════
//  2. One explicit merge() call ([C91 §3 preamble tex 160–165]):
//     two adjacent chains, granularity reset (§3.3), real oracles
//     (§4.1 tex 341), full validation of the result.
// ════════════════════════════════════════════════════════════════

static void test_direct_merge(bool dup_ys) {
    Rng rng(7);
    UpPhase up(random_curve(rng, 40, dup_ys));   // n = 65, p = 6
    assert(up.graded().p() == 6);

    // Merge the two grade-2 chains forming grade-3 chain 1
    // (vertices [8,16]), at λ = 3 ⟹ γ = 2^{⌈3/5⌉} = 2.
    const std::size_t lambda = 3;
    const std::size_t gamma = UpPhase::grade_gamma(lambda);
    Polygon c1 = up.graded().chain(2, 2);
    Polygon c2 = up.graded().chain(2, 3);
    Submap s1 = up.chain_submap(2, 2);           // copies
    Submap s2 = up.chain_submap(2, 3);

    // [C91 §4.1 tex 349]: granularity reset via §3.3.
    enforce_granularity(s1, c1, gamma);
    s1.normalize(c1);
    enforce_granularity(s2, c2, gamma);
    s2.normalize(c2);

    UpPhaseRayShooter ray1(up, s1, c1, lambda);
    UpPhaseRayShooter ray2(up, s2, c2, lambda);
    UpPhaseArcCutter cut1(up, s1, c1, lambda);
    UpPhaseArcCutter cut2(up, s2, c2, lambda);

    MergeInput in;
    in.C1 = &c1;
    in.C2 = &c2;
    in.S1 = &s1;
    in.S2 = &s2;
    in.gamma1 = in.gamma2 = in.gamma = gamma;
    in.ray_shooter_1 = &ray1;
    in.ray_shooter_2 = &ray2;
    in.arc_cutter_1 = &cut1;
    in.arc_cutter_2 = &cut2;
    in.g_gamma1 = in.g_gamma2 = UpPhaseArcCutter::g_bound(lambda);
    in.h_gamma1 = in.h_gamma2 = UpPhaseArcCutter::h_bound(lambda);

    MergeResult res = merge(in);
    assert(res.C.num_vertices() ==
               c1.num_vertices() + c2.num_vertices() - 1 &&
           res.C.table_offset() == c1.table_offset() &&
           "[C91 §3 tex 160]: C = C₁ ∪ C₂ sharing the junction vertex");
    assert_canonical(res.S, res.C, gamma);
    assert_chords_are_visible_pairs(res.S, res.C);

    std::printf("  [PASS] direct_merge (dup_ys=%d, chords=%zu)\n",
                (int)dup_ys, res.S.num_live_chords());
}

// ════════════════════════════════════════════════════════════════
//  3. [C91 §3.4 Lemma 3.6] end-to-end: the ray-shooting structure
//     over merged submaps agrees with the naive shooter ray-for-ray
//     (faces, dual graph, LT79 separators, vertical line, double
//     identification — the whole apparatus).
// ════════════════════════════════════════════════════════════════

static void test_ray_shooting_vs_naive(bool dup_ys) {
    for (unsigned seed = 11; seed <= 13; ++seed) {
        Rng rng(seed);
        const std::size_t n = 20 + rng.next() % 40;
        UpPhase up(random_curve(rng, n, dup_ys));

        for (std::size_t lam = 2; lam <= up.graded().p(); ++lam) {
            const std::size_t i =
                rng.next() % up.graded().num_chains(lam);
            const Polygon& c = up.graded().chain(lam, i);
            const RayShootingStructure& rs = up.chain_structure(lam, i);

            // Sources: every wall's midpoint (mid-edge raw level, tag
            // SOS_NONE — symbolically just below any raw-tied vertex
            // cluster, [C91 §2 tex 47]), plus every vertex level from
            // both walls of each incident edge.
            auto compare = [&](Point p, SymbolicY sy, Side dir,
                               std::size_t src_edge) {
                const double off =
                    perturbed_x_offset(c, sy, src_edge);
                RayHit a = rs.shoot_toward_boundary(p, dir, off);
                RayHit b = naive_first_contact(c, p, sy, dir, src_edge);
                assert(a.hit == b.hit &&
                       "[C91 §3.4 Lemma 3.6]: structure and naive "
                       "shooter must agree on hit existence");
                if (a.hit) {
                    // [C91 §3.0(i) tex 169]: a contact at a
                    // NON-extremum interior vertex "is contained in
                    // both incident edges" — the two shooters may pick
                    // either admissible label for the SAME wall.
                    bool same_label = (a.edge == b.edge);
                    if (!same_label && a.x == b.x) {
                        const std::size_t v = std::max(a.edge, b.edge);
                        same_label =
                            (a.edge + 1 == b.edge ||
                             b.edge + 1 == a.edge) &&
                            v < c.num_vertices() &&
                            symbolic_y_equal(sy,
                                             symbolic_y_of(
                                                 c.vertex(v))) &&
                            !c.is_y_extremum(v);
                    }
                    const bool ok = a.x == b.x && same_label &&
                                    a.side == b.side &&
                                    a.wrapped == b.wrapped;
                    if (!ok)
                        std::fprintf(stderr,
                            "MISMATCH p=(%.6g,%.6g@%zu) dir=%c src=%zu "
                            "struct=(%zu,%c,x=%.6g,w=%d) "
                            "naive=(%zu,%c,x=%.6g,w=%d)\n",
                            p.x, sy.y, sy.tag, dir == LEFT ? 'L' : 'R',
                            src_edge, a.edge,
                            a.side == LEFT ? 'L' : 'R', a.x,
                            (int)a.wrapped, b.edge,
                            b.side == LEFT ? 'L' : 'R', b.x,
                            (int)b.wrapped);
                    assert(ok &&
                           "[C91 §3.4 Lemma 3.6]: structure and naive "
                           "shooter must report the same first contact");
                }
            };

            for (std::size_t e = 0; e < c.num_edges(); ++e) {
                if (c.edge_is_null(e)) continue;
                const auto& ed = c.edge(e);
                const double y0 = c.vertex(ed.start_idx).y;
                const double y1 = c.vertex(ed.end_idx).y;
                SymbolicY mid{(y0 + y1) / 2.0, SOS_NONE};
                double x;
                if (!edge_crossing_x(c, e, mid, &x)) continue;
                for (Side s : {LEFT, RIGHT}) {
                    Point p{x, mid.y, mid.tag};
                    compare(p, mid, shooting_direction(e, s, c), e);
                }
            }
            for (std::size_t v = 0; v < c.num_vertices(); ++v) {
                SymbolicY vy = symbolic_y_of(c.vertex(v));
                for (std::size_t e :
                     {v > 0 ? v - 1 : NONE,
                      v < c.num_edges() ? v : NONE}) {
                    if (e == NONE) continue;
                    for (Side s : {LEFT, RIGHT}) {
                        if (is_inside_companion(c, e, s, v))
                            continue;   // sees only its d=0 sibling
                        Point p{c.vertex(v).x, vy.y, vy.tag};
                        compare(p, vy, shooting_direction(e, s, c), e);
                    }
                }
            }
        }
    }
    std::printf("  [PASS] ray_shooting_vs_naive (dup_ys=%d)\n",
                (int)dup_ys);
}

// ════════════════════════════════════════════════════════════════
//  4. double_identify on MERGED submaps vs the position-exact
//     reference ([C91 §2.4 tex 144]) — merged arcs wrap, double-back,
//     and carry junction labels the full-V(C) test never produces.
// ════════════════════════════════════════════════════════════════

static void test_double_identify_on_merged(bool dup_ys) {
    Rng rng(17);
    UpPhase up(random_curve(rng, 50, dup_ys));   // n = 65, p = 6

    const Polygon& c = up.graded().chain(up.graded().p(), 0);
    Submap s = up.chain_submap(up.graded().p(), 0);   // copy
    s.compact();

    auto reference = [&](std::size_t e, SymbolicY qy) {
        std::vector<std::size_t> out;
        for (std::size_t ai = 0; ai < s.num_arcs(); ++ai) {
            if (s.arc(ai).dead) continue;
            const Arc& a = s.arc(ai);
            Subarc full;
            full.first_edge = a.first_edge;
            full.first_side = a.first_side;
            full.last_edge = a.last_edge;
            full.last_side = a.last_side;
            full.first_y = s.arc_start_symbolic_y(ai, c);
            full.last_y = s.arc_end_symbolic_y(ai, c);
            if (subarc_contains_point(full, c, e, LEFT, qy, 0,
                                      c.num_vertices() - 1) ||
                subarc_contains_point(full, c, e, RIGHT, qy, 0,
                                      c.num_vertices() - 1))
                out.push_back(ai);
        }
        return out;
    };

    for (std::size_t e = 0; e < c.num_edges(); ++e) {
        const auto& ed = c.edge(e);
        // On-edge test must be SYMBOLIC ([C91 §2 tex 45]): with
        // duplicate raw ys a raw-tied query can fall symbolically
        // outside the edge's span (larger tag ⟹ lower).
        SymbolicY sy0 = symbolic_y_of(c.vertex(ed.start_idx));
        SymbolicY sy1 = symbolic_y_of(c.vertex(ed.end_idx));
        if (symbolic_y_greater(sy0, sy1)) std::swap(sy0, sy1);
        for (std::size_t v = 0; v < c.num_vertices(); ++v) {
            SymbolicY qy = symbolic_y_of(c.vertex(v));
            if (symbolic_y_less(qy, sy0) || symbolic_y_greater(qy, sy1))
                continue;
            auto res = s.double_identify(e, qy, c);
            auto ref = reference(e, qy);
            assert(res.count == ref.size() &&
                   "[C91 §2.4 tex 144]: double_identify must return "
                   "EXACTLY the arcs containing the queried edge");
            for (std::size_t ai : ref) {
                bool found = false;
                for (std::size_t k = 0; k < res.count; ++k)
                    found = found || res.arcs[k] == ai;
                assert(found);
            }
        }
    }
    std::printf("  [PASS] double_identify_on_merged (dup_ys=%d)\n",
                (int)dup_ys);
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::printf("[C91 §3 e2e tests]:\n");
    for (bool dup : {false, true}) {
        test_merged_portions_vs_naive(dup);
        test_direct_merge(dup);
        test_ray_shooting_vs_naive(dup);
        test_double_identify_on_merged(dup);
    }
    std::printf("[C91 §3 e2e tests]: all passed\n");
    return 0;
}
