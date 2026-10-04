#include "merge/fusion.h"
#include "merge/granularity.h"
#include "merge/merge.h"
#include "merge/ray_shooting.h"
#include "polygon/polygon.h"
#include "submap/submap.h"
#include "support/random.h"
#include "visibility/naive_visibility.h"
#include "visibility/up_phase.h"

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <vector>

using namespace chazelle;
using chazelle::test::DeterministicRandomGenerator;

namespace {

enum class CurveTransform : std::uint8_t { none, swap_axes, rotate, half_turn };

std::vector<Point> random_curve(DeterministicRandomGenerator& rng, std::size_t n, bool dup_ys,
                                CurveTransform transform = CurveTransform::none) {
    std::vector<Point> v;
    Exact x = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        x += rng.uniform(0.5, 2.0);
        Exact y = dup_ys ? Exact(rng.next() % 8) : rng.uniform(0.0, 40.0);
        Point p{x, y, i};
        switch (transform) {
        case CurveTransform::none:
            break;
        case CurveTransform::swap_axes:
            p.x = y;
            p.y = x;
            break;
        case CurveTransform::rotate:
            p.x = x + y;
            p.y = y - x;
            break;
        case CurveTransform::half_turn:
            p.x = -x;
            p.y = -y;
            break;
        }
        v.push_back(p);
    }
    return v;
}

void assert_chords_are_visible_pairs(const Submap& submap, const Polygon& curve) {
    for (std::size_t ci = 0; ci < submap.num_chords(); ++ci) {
        const Chord& ch = submap.chord(ci);
        if (ch.dead || ch.is_null_length)
            continue;
        struct End {
            std::size_t e;
            Side s;
        };
        const End ends[2] = {{ch.left_edge, ch.left_side}, {ch.right_edge, ch.right_side}};
        for (int k = 0; k < 2; ++k) {
            const End& a = ends[k];
            const End& b = ends[1 - k];
            const Exact xa = edge_x_at_y(curve, a.e, ch.symbolic_y());
            const Exact xb = edge_x_at_y(curve, b.e, ch.symbolic_y());
            Point p{xa, ch.symbolic_y().y, ch.symbolic_y().tag};
            RayHit h = naive_first_contact(curve, p, ch.symbolic_y(),
                                           shooting_direction(a.e, a.s, curve), a.e);

            [[maybe_unused]] const bool label_ok =
                h.edge == b.e ||
                (symbolic_y_equal(ch.symbolic_y(),
                                  symbolic_y_of(curve.vertex(std::max(h.edge, b.e)))) &&
                 !curve.is_y_extremum(std::max(h.edge, b.e)) &&
                 (h.edge + 1 == b.e || b.e + 1 == h.edge));
            assert(h.hit && h.x == xb && h.side == b.s && label_ok &&
                   "[C91 §2.2 tex 90]: merged-submap chord endpoints "
                   "must see each other with respect to C");
        }
    }
}

void assert_canonical(const Submap& submap, const Polygon& curve, std::size_t granularity) {
    submap.check_invariants(curve);
    assert(submap.is_conformal() && "[C91 §3.2 Lemma 3.4]: merge output is conformal");
    assert(submap.is_granular(granularity, curve) &&
           "[C91 §3.3 Lemma 3.5]: merge output is γ-granular");
    assert(!submap.tree_decomposition().empty() &&
           "[C91 §2.4(iv)]: normal form carries the tree decomposition");
}

}

static void test_merged_portions_vs_naive(bool dup_ys,
                                          CurveTransform transform = CurveTransform::none) {
    for (unsigned seed = 1; seed <= 4; ++seed) {
        DeterministicRandomGenerator rng(seed);
        const std::size_t n = 20 + rng.next() % 60;
        UpPhase up(random_curve(rng, n, dup_ys, transform));
        const Polygon& input_curve = up.graded().curve();

        for (std::size_t lam = 0; lam <= up.graded().maximum_grade(); ++lam)
            for (std::size_t i = 0; i < up.graded().num_chains(lam); ++i) {
                const Polygon& c = up.graded().chain(lam, i);
                const Submap& s = up.chain_submap(lam, i);
                assert_canonical(s, c, UpPhase::grade_granularity(lam));
                assert_chords_are_visible_pairs(s, c);
            }

        const std::size_t np = input_curve.num_vertices();
        for (int t = 0; t < 3; ++t) {
            std::size_t a = rng.next() % (np - 4);
            std::size_t len = 4 + rng.next() % (np - a - 4);
            if (len < 3)
                continue;
            auto pr = up.compute_canonical_portion(a, a + len);
            std::size_t lam = 0;
            while ((std::size_t{1} << lam) < len)
                ++lam;
            assert_canonical(pr.submap, pr.curve, UpPhase::grade_granularity(lam));
            assert_chords_are_visible_pairs(pr.submap, pr.curve);
        }
    }
    std::printf("  [PASS] merged_portions_vs_naive (dup_ys=%d)\n", (int)dup_ys);
}

static void test_direct_merge(bool dup_ys) {
    DeterministicRandomGenerator rng(7);
    UpPhase up(random_curve(rng, 40, dup_ys));
    assert(up.graded().maximum_grade() == 6);

    const std::size_t lambda = 3;
    const std::size_t granularity = UpPhase::grade_granularity(lambda);
    Polygon c1 = up.graded().chain(2, 2);
    Polygon c2 = up.graded().chain(2, 3);
    Submap s1 = up.chain_submap(2, 2);
    Submap s2 = up.chain_submap(2, 3);

    enforce_granularity(s1, c1, granularity);
    s1.normalize(c1);
    enforce_granularity(s2, c2, granularity);
    s2.normalize(c2);

    UpPhaseRayShooter first_ray_shooter(up, s1, c1, lambda);
    UpPhaseRayShooter second_ray_shooter(up, s2, c2, lambda);
    UpPhaseArcCutter first_arc_cutter(up, s1, c1, lambda);
    UpPhaseArcCutter second_arc_cutter(up, s2, c2, lambda);

    MergeInput in;
    in.first_curve = &c1;
    in.second_curve = &c2;
    in.first_submap = &s1;
    in.second_submap = &s2;
    in.first_granularity = in.second_granularity = in.granularity = granularity;
    in.first_ray_shooter = &first_ray_shooter;
    in.second_ray_shooter = &second_ray_shooter;
    in.first_arc_cutter = &first_arc_cutter;
    in.second_arc_cutter = &second_arc_cutter;
    in.first_piece_count_bound = in.second_piece_count_bound =
        UpPhaseArcCutter::piece_count_bound(lambda);
    in.first_piece_granularity_bound = in.second_piece_granularity_bound =
        UpPhaseArcCutter::piece_granularity_bound(lambda);

    MergeResult res = merge(in);
    assert(res.curve.num_vertices() == c1.num_vertices() + c2.num_vertices() - 1 &&
           res.curve.table_offset() == c1.table_offset() &&
           "[C91 §3 tex 160]: C = C₁ ∪ C₂ sharing the junction vertex");
    assert_canonical(res.submap, res.curve, granularity);
    assert_chords_are_visible_pairs(res.submap, res.curve);

    std::printf("  [PASS] direct_merge (dup_ys=%d, chords=%zu)\n", (int)dup_ys,
                res.submap.num_live_chords());
}

static void test_ray_shooting_vs_naive(bool dup_ys,
                                       CurveTransform transform = CurveTransform::none) {
    for (unsigned seed = 11; seed <= 13; ++seed) {
        DeterministicRandomGenerator rng(seed);
        const std::size_t n = 20 + rng.next() % 40;
        UpPhase up(random_curve(rng, n, dup_ys, transform));

        for (std::size_t lam = 2; lam <= up.graded().maximum_grade(); ++lam) {
            const std::size_t i = rng.next() % up.graded().num_chains(lam);
            const Polygon& c = up.graded().chain(lam, i);
            const RayShootingStructure& rs = up.chain_structure(lam, i);

            auto compare = [&](const Point& p, const SymbolicY& sy, Side dir,
                               std::size_t src_edge) {
                const Exact off = perturbed_x_offset(c, sy, src_edge);
                RayHit a = rs.shoot_toward_boundary(p, dir, off);
                RayHit b = naive_first_contact(c, p, sy, dir, src_edge);
                assert(a.hit == b.hit && "[C91 §3.4 Lemma 3.6]: structure and naive "
                                         "shooter must agree on hit existence");
                if (a.hit) {
                    bool same_label = (a.edge == b.edge);
                    if (!same_label && a.x == b.x) {
                        const std::size_t v = std::max(a.edge, b.edge);
                        same_label = (a.edge + 1 == b.edge || b.edge + 1 == a.edge) &&
                                     v < c.num_vertices() &&
                                     symbolic_y_equal(sy, symbolic_y_of(c.vertex(v))) &&
                                     !c.is_y_extremum(v);
                    }
                    assert(a.x == b.x && same_label && a.side == b.side && a.wrapped == b.wrapped &&
                           "[C91 §3.4 Lemma 3.6]: structure and naive "
                           "shooter must report the same first contact");
                }
            };

            for (std::size_t e = 0; e < c.num_edges(); ++e) {
                if (c.edge_is_null(e))
                    continue;
                const auto& ed = c.edge(e);
                const Exact y0 = c.vertex(ed.start_idx).y;
                const Exact y1 = c.vertex(ed.end_idx).y;
                SymbolicY mid{(y0 + y1) / 2.0, SOS_NONE};
                Exact x;
                if (!edge_crossing_x(c, e, mid, &x))
                    continue;
                for (Side s : {LEFT, RIGHT}) {
                    Point p{x, mid.y, mid.tag};
                    compare(p, mid, shooting_direction(e, s, c), e);
                }
            }
            for (std::size_t v = 0; v < c.num_vertices(); ++v) {
                SymbolicY vy = symbolic_y_of(c.vertex(v));
                for (std::size_t e : {v > 0 ? v - 1 : NONE, v < c.num_edges() ? v : NONE}) {
                    if (e == NONE)
                        continue;
                    for (Side s : {LEFT, RIGHT}) {
                        if (is_inside_companion(c, e, s, v))
                            continue;
                        Point p{c.vertex(v).x, vy.y, vy.tag};
                        compare(p, vy, shooting_direction(e, s, c), e);
                    }
                }
            }
        }
    }
    std::printf("  [PASS] ray_shooting_vs_naive (dup_ys=%d)\n", (int)dup_ys);
}

static void test_double_identify_on_merged(bool dup_ys) {
    DeterministicRandomGenerator rng(17);
    UpPhase up(random_curve(rng, 50, dup_ys));

    const Polygon& c = up.graded().chain(up.graded().maximum_grade(), 0);
    Submap s = up.chain_submap(up.graded().maximum_grade(), 0);
    s.compact();

    auto reference = [&](std::size_t e, const SymbolicY& qy) {
        std::vector<std::size_t> out;
        for (std::size_t ai = 0; ai < s.num_arcs(); ++ai) {
            if (s.arc(ai).dead)
                continue;
            const Arc& a = s.arc(ai);
            Subarc full;
            full.first_edge = a.first_edge;
            full.first_side = a.first_side;
            full.last_edge = a.last_edge;
            full.last_side = a.last_side;
            full.first_y = s.arc_start_symbolic_y(ai, c);
            full.last_y = s.arc_end_symbolic_y(ai, c);
            if (subarc_contains_point(full, c, e, LEFT, qy, 0, c.num_vertices() - 1) ||
                subarc_contains_point(full, c, e, RIGHT, qy, 0, c.num_vertices() - 1))
                out.push_back(ai);
        }
        return out;
    };

    for (std::size_t e = 0; e < c.num_edges(); ++e) {
        const auto& ed = c.edge(e);

        SymbolicY sy0 = symbolic_y_of(c.vertex(ed.start_idx));
        SymbolicY sy1 = symbolic_y_of(c.vertex(ed.end_idx));
        if (symbolic_y_greater(sy0, sy1))
            std::swap(sy0, sy1);
        for (std::size_t v = 0; v < c.num_vertices(); ++v) {
            SymbolicY qy = symbolic_y_of(c.vertex(v));
            if (symbolic_y_less(qy, sy0) || symbolic_y_greater(qy, sy1))
                continue;
            auto res = s.double_identify(e, qy, c);
            auto ref = reference(e, qy);
            assert(res.count == ref.size() && "[C91 §2.4 tex 144]: double_identify must return "
                                              "EXACTLY the arcs containing the queried edge");
            for (std::size_t ai : ref) {
                bool found = false;
                for (std::size_t k = 0; k < res.count; ++k)
                    found = found || res.arcs[k] == ai;
                assert(found);
            }
        }
    }
    std::printf("  [PASS] double_identify_on_merged (dup_ys=%d)\n", (int)dup_ys);
}

int main() {
    std::printf("[C91 §3 e2e tests]:\n");
    for (bool dup : {false, true}) {
        for (CurveTransform transform : {CurveTransform::none, CurveTransform::swap_axes,
                                         CurveTransform::rotate, CurveTransform::half_turn}) {
            test_merged_portions_vs_naive(dup, transform);
            test_ray_shooting_vs_naive(dup, transform);
        }
        test_direct_merge(dup);
        test_double_identify_on_merged(dup);
    }
    std::printf("[C91 §3 e2e tests]: all passed\n");
    return 0;
}
