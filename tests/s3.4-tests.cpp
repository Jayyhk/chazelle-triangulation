#include "algorithm/merge/conformality.h"
#include "algorithm/merge/fusion.h"
#include "algorithm/merge/granularity.h"
#include "algorithm/merge/merge.h"
#include "algorithm/merge/ray_shooting.h"
#include "algorithm/polygon/polygon.h"
#include "algorithm/submap/chord_inventory.h"
#include "algorithm/submap/submap.h"
#include "algorithm/visibility/up_phase.h"
#include "support/arc_ray_shooter.h"
#include "support/assertions.h"
#include "support/random.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <memory>
#include <vector>

using namespace chazelle;
using chazelle::test::DeterministicRandomGenerator;
using chazelle::test::require_assertion_abort;

static const Polygon& input_P() {
    static Polygon input_curve({{0, 0, 0},
                                {2, 20, 1},
                                {4, 6, 2},
                                {6, 24, 3},
                                {8, 4, 4},
                                {10, 22, 5},
                                {12, 5, 6},
                                {14, 26, 7},
                                {16, 2, 8},
                                {18, 23, 9},
                                {20, 4.5, 10},
                                {22, 25, 11},
                                {24, 1, 12}});
    return input_curve;
}
static Polygon make_C() {
    return input_P().subchain(0, 13);
}
static Polygon make_C1() {
    return input_P().subchain(0, 8);
}
static Polygon make_C2() {
    return input_P().subchain(7, 6);
}

static Submap make_chordless_normal(const Polygon& poly) {
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
    a.edge_count = 2 * poly.count_nonnull_edges(0, poly.num_edges() - 1);
    std::size_t ai0 = s.add_arc(a);
    assert(s.start_arc == ai0 && s.end_arc == ai0);
    s.build_tree_decomposition();
    return s;
}

struct CombFixture {
    Polygon curve, first_curve, second_curve;
    Submap first_submap, second_submap;
    Submap submap;
    CombFixture()
        : curve(make_C()), first_curve(make_C1()), second_curve(make_C2()),
          first_submap(make_chordless_normal(first_curve)),
          second_submap(make_chordless_normal(second_curve)) {
        build();
    }
    void build() {
        for (int i = 0; i < 9; ++i)
            submap.add_node();
        auto add = [&](std::size_t fe, Side fs, std::size_t le, Side ls, std::size_t region,
                       std::size_t count) -> std::size_t {
            Arc a{};
            a.first_edge = fe;
            a.first_side = fs;
            a.last_edge = le;
            a.last_side = ls;
            a.region_node = region;
            a.edge_count = count;
            return submap.add_arc(a);
        };
        submap.start_vertex = 0;
        submap.end_vertex = 12;

        add(7, LEFT, 7, LEFT, 7, 0);
        add(7, LEFT, 11, RIGHT, 0, 6);
        add(11, RIGHT, 10, RIGHT, 6, 2);
        add(9, RIGHT, 9, RIGHT, 0, 0);
        add(9, RIGHT, 8, RIGHT, 5, 2);
        add(8, RIGHT, 7, RIGHT, 0, 2);
        add(7, RIGHT, 7, RIGHT, 4, 1);
        add(7, RIGHT, 7, RIGHT, 8, 0);
        add(6, RIGHT, 6, RIGHT, 4, 1);
        add(5, RIGHT, 5, RIGHT, 0, 0);
        add(5, RIGHT, 4, RIGHT, 3, 2);
        add(4, RIGHT, 3, RIGHT, 0, 2);
        add(3, RIGHT, 2, RIGHT, 2, 2);
        add(1, RIGHT, 1, RIGHT, 0, 0);
        add(1, RIGHT, 0, RIGHT, 1, 2);
        add(0, RIGHT, 6, LEFT, 0, 8);
        auto chord = [&](std::size_t le, Side lsd, Chord::AdjArcs ladj, std::size_t re, Side rsd,
                         Chord::AdjArcs radj, const Exact& y, std::size_t tag, std::size_t r0,
                         std::size_t r1, bool null_len = false) {
            Chord c{};
            c.left_edge = le;
            c.left_side = lsd;
            c.left_adj = ladj;
            c.right_edge = re;
            c.right_side = rsd;
            c.right_adj = radj;
            c.y = y;
            c.y_tag = tag;
            c.region[0] = r0;
            c.region[1] = r1;
            c.is_null_length = null_len;
            submap.add_chord(c);
        };
        chord(0, RIGHT, {{14, 15}, 2}, 1, RIGHT, {{13}, 1}, 6.0, 2, 0, 1);
        chord(2, RIGHT, {{12}, 1}, 3, RIGHT, {{11, 12}, 2}, 6.0, 2, 0, 2);
        chord(4, RIGHT, {{10, 11}, 2}, 5, RIGHT, {{9}, 1}, 5.0, 6, 0, 3);
        chord(6, RIGHT, {{8}, 1}, 7, RIGHT, {{5, 6}, 2}, 5.0, 6, 0, 4);
        chord(8, RIGHT, {{4, 5}, 2}, 9, RIGHT, {{3}, 1}, 4.5, 10, 0, 5);
        chord(10, RIGHT, {{2}, 1}, 11, RIGHT, {{1, 2}, 2}, 4.5, 10, 0, 6);
        chord(6, LEFT, {{15}, 1}, 7, LEFT, {{0}, 1}, 26.0, 7, 0, 7);
        chord(7, RIGHT, {{6}, 1}, 7, RIGHT, {{7}, 1}, 26.0, 7, 4, 8, true);

        assert(submap.start_arc == 15 && submap.end_arc == 1);
    }
};

struct GeomRayShooter : RayShootingOracle {
    const Polygon* input_curve;
    explicit GeomRayShooter(const Polygon* c) : input_curve(c) {}

    static bool crossing_x(const Polygon& curve, std::size_t e, const SymbolicY& sy, Exact* x) {
        const auto& ed = curve.edge(e);
        const Point& vs = curve.vertex(ed.start_idx);
        const Point& ve = curve.vertex(ed.end_idx);
        SymbolicY y0 = symbolic_y_of(vs);
        SymbolicY y1 = symbolic_y_of(ve);
        if (symbolic_y_equal(sy, y0)) {
            *x = vs.x;
            return true;
        }
        if (symbolic_y_equal(sy, y1)) {
            *x = ve.x;
            return true;
        }
        bool between = (symbolic_y_less(y0, sy) && symbolic_y_less(sy, y1)) ||
                       (symbolic_y_less(y1, sy) && symbolic_y_less(sy, y0));
        if (!between)
            return false;
        Exact t = (sy.y - vs.y) / (ve.y - vs.y);
        *x = vs.x + t * (ve.x - vs.x);
        return true;
    }

    RayHit shoot(Point p, Side dir, std::size_t, const Subarc& target,
                 SourceOffset = SOURCE_OFFSET_NONE) const override {
        SymbolicY sy{p.y, p.index};

        ArcSideRange ranges[3];
        std::size_t nl = subarc_side_ranges(target, 0, input_curve->num_vertices() - 1, ranges);
        RayHit best;
        best.hit = false;
        Exact best_d = 0.0;
        for (std::size_t g = 0; g < nl; ++g) {
            for (std::size_t e = ranges[g].first_edge; e <= ranges[g].last_edge; ++e) {
                Exact x;
                if (!crossing_x(*input_curve, e, sy, &x))
                    continue;
                const auto& ed = input_curve->edge(e);
                bool asc = symbolic_y_less(symbolic_y_of(input_curve->vertex(ed.start_idx)),
                                           symbolic_y_of(input_curve->vertex(ed.end_idx)));
                Side minus_x = asc ? LEFT : RIGHT;
                Side struck = (dir == RIGHT) ? minus_x : (minus_x == LEFT ? RIGHT : LEFT);

                if (!subarc_contains_point(target, *input_curve, e, struck, sy, 0,
                                           input_curve->num_vertices() - 1))
                    continue;
                Exact d = (dir == RIGHT) ? (x - p.x) : (p.x - x);
                bool wrapped = (d <= 0.0);
                bool better;
                if (!best.hit)
                    better = true;
                else if (wrapped != best.wrapped)
                    better = !wrapped;
                else
                    better = d < best_d;
                if (better) {
                    best.hit = true;
                    best.x = x;
                    best.y = p.y;
                    best.edge = e;
                    best.side = struck;
                    best.wrapped = wrapped;
                    best_d = d;
                }
            }
        }
        return best;
    }
};

struct SafeArcCutter : ArcCuttingOracle {
    const Polygon* input_curve;
    std::size_t piece_len;
    mutable std::vector<std::unique_ptr<Polygon>> curves;
    mutable std::vector<std::unique_ptr<Submap>> submaps;

    SafeArcCutter(const Polygon* ci, std::size_t pl) : input_curve(ci), piece_len(pl) {}

    std::vector<ArcPiece> cut_side_range(const Subarc& leg, bool target_first,
                                         bool target_last) const {
        const Side s = leg.first_side;
        const std::size_t lo = std::min(leg.first_edge, leg.last_edge);
        const std::size_t hi = std::max(leg.first_edge, leg.last_edge);
        const bool blo = (s == LEFT) ? target_first : target_last;
        const bool bhi = (s == LEFT) ? target_last : target_first;

        auto boundary_piece = [&](std::size_t e) {
            ArcPiece p;
            p.subarc = Subarc{e, s, e, s};
            p.is_boundary_piece = true;
            return p;
        };
        auto middle_piece = [&](std::size_t a, std::size_t b) {
            std::vector<Point> vs;
            for (std::size_t v = a; v <= b + 1; ++v)
                vs.push_back(input_curve->vertex(v));
            curves.push_back(std::make_unique<Polygon>(std::move(vs)));
            submaps.push_back(std::make_unique<Submap>(make_chordless_normal(*curves.back())));
            ArcPiece p;
            p.subarc = (s == LEFT) ? Subarc{a, s, b, s} : Subarc{b, s, a, s};
            p.curve = curves.back().get();
            p.submap = submaps.back().get();

            p.granularity = 2 * (b - a + 1);
            return p;
        };

        std::vector<ArcPiece> out;
        if (lo == hi && (blo || bhi)) {
            out.push_back(boundary_piece(lo));
            return out;
        }
        std::size_t mlo = lo + (blo ? 1 : 0);
        std::size_t mhi = hi - (bhi ? 1 : 0);
        std::vector<std::size_t> chunk_lo, chunk_hi;
        for (std::size_t a = mlo; a <= mhi; a += piece_len) {
            chunk_lo.push_back(a);
            chunk_hi.push_back(std::min(a + piece_len - 1, mhi));
        }
        if (s == LEFT) {
            if (blo)
                out.push_back(boundary_piece(lo));
            for (std::size_t k = 0; k < chunk_lo.size(); ++k)
                out.push_back(middle_piece(chunk_lo[k], chunk_hi[k]));
            if (bhi)
                out.push_back(boundary_piece(hi));
        } else {
            if (bhi)
                out.push_back(boundary_piece(hi));
            for (std::size_t k = chunk_lo.size(); k-- > 0;)
                out.push_back(middle_piece(chunk_lo[k], chunk_hi[k]));
            if (blo)
                out.push_back(boundary_piece(lo));
        }
        return out;
    }

    std::vector<ArcPiece> cut(std::size_t, const Subarc& target) const override {
        ArcSideRange ranges[3];
        std::size_t nl = subarc_side_ranges(target, 0, input_curve->num_vertices() - 1, ranges);
        std::vector<ArcPiece> out;
        for (std::size_t g = 0; g < nl; ++g) {
            Subarc subarc_on_side =
                (ranges[g].side == LEFT)
                    ? Subarc{ranges[g].first_edge, LEFT, ranges[g].last_edge, LEFT}
                    : Subarc{ranges[g].last_edge, RIGHT, ranges[g].first_edge, RIGHT};
            auto pieces = cut_side_range(subarc_on_side, g == 0, g + 1 == nl);
            out.insert(out.end(), pieces.begin(), pieces.end());
        }
        test_set_cut_endpoints(*input_curve, target, out);
        return out;
    }
};

static RayHit brute_shoot(const Polygon& curve, const Point& p, Side dir) {
    SymbolicY sy{p.y, p.index};
    RayHit best;
    Exact best_d = 0.0;
    for (std::size_t e = 0; e < curve.num_edges(); ++e) {
        Exact x;
        if (!GeomRayShooter::crossing_x(curve, e, sy, &x))
            continue;
        const auto& ed = curve.edge(e);
        bool asc = symbolic_y_less(symbolic_y_of(curve.vertex(ed.start_idx)),
                                   symbolic_y_of(curve.vertex(ed.end_idx)));
        Side minus_x = asc ? LEFT : RIGHT;
        Side struck = (dir == RIGHT) ? minus_x : (minus_x == LEFT ? RIGHT : LEFT);
        Exact d = (dir == RIGHT) ? (x - p.x) : (p.x - x);
        bool wrapped = (d <= 0.0);
        bool better;
        if (!best.hit)
            better = true;
        else if (wrapped != best.wrapped)
            better = !wrapped;
        else
            better = d < best_d;
        if (better) {
            best.hit = true;
            best.x = x;
            best.y = p.y;
            best.edge = e;
            best.side = struck;
            best.wrapped = wrapped;
            best_d = d;
        }
    }
    return best;
}

static void check_against_brute(const RayShootingStructure& rs, const Polygon& curve,
                                const Point& p, Side dir) {
    RayHit got = rs.shoot_toward_boundary(p, dir);
    RayHit want = brute_shoot(curve, p, dir);
    assert(got.hit == want.hit && "[C91 Lemma 3.6 tex 310–312]: structure and brute force must "
                                  "agree on hit existence");
    if (!got.hit)
        return;
    assert(got.wrapped == want.wrapped && got.x == want.x &&
           "[C91 Lemma 3.6 tex 310–312]: structure must report the first "
           "contact in the wrap metric");

    if (got.edge != want.edge) {
        Exact xa, xb;
        assert(GeomRayShooter::crossing_x(curve, got.edge, SymbolicY{p.y, p.index}, &xa) &&
               GeomRayShooter::crossing_x(curve, want.edge, SymbolicY{p.y, p.index}, &xb) &&
               xa == xb && "differing edges only at a coincident contact");
    } else {
        assert(got.side == want.side);
    }
}

struct GranularComb {
    CombFixture fx;
    GeomRayShooter first_ray_shooter, second_ray_shooter;
    SafeArcCutter first_arc_cutter, second_arc_cutter;
    std::size_t granularity = 0;

    GranularComb()
        : first_ray_shooter(&fx.first_curve), second_ray_shooter(&fx.second_curve),
          first_arc_cutter(&fx.first_curve, 64), second_arc_cutter(&fx.second_curve, 64) {
        ConformalityOracles o;
        o.first_submap = &fx.first_submap;
        o.second_submap = &fx.second_submap;
        o.first_curve = &fx.first_curve;
        o.second_curve = &fx.second_curve;
        o.first_ray_shooter = &first_ray_shooter;
        o.second_ray_shooter = &second_ray_shooter;
        o.first_arc_cutter = &first_arc_cutter;
        o.second_arc_cutter = &second_arc_cutter;
        o.first_piece_count_bound = 64;
        o.second_piece_count_bound = 64;
        o.first_piece_granularity_bound = 64;
        o.second_piece_granularity_bound = 64;
        restore_conformality(fx.submap, fx.curve, o);

        for (std::size_t r = 0; r < fx.submap.num_nodes(); ++r)
            if (!fx.submap.node(r).dead)
                granularity = std::max(granularity, fx.submap.region_weight(r));
        enforce_granularity(fx.submap, fx.curve, granularity);
        fx.submap.normalize(fx.curve);
    }
};

struct RecursiveFixture {
    static std::vector<Point> vertices() {
        std::vector<Point> v;
        v.reserve(129);
        for (std::size_t i = 0; i < 129; ++i)
            v.push_back({Exact(i), (i % 2 ? 20.0 : 0.0) + Exact(i) / 1024.0, i});
        return v;
    }
    UpPhase up{vertices()};
    const Polygon& curve = up.graded().curve();
    const Submap& submap = up.chain_submap(up.graded().maximum_grade(), 0);
    std::size_t granularity = UpPhase::grade_granularity(up.graded().maximum_grade());
};

static void test_structure_trivial_mu1() {
    Polygon curve = make_C();
    Submap submap = make_chordless_normal(curve);

    RayShootingStructure rs(submap, curve, 2 * curve.num_edges());
    assert(rs.num_faces() == 1 && "[C91 §3.4 tex 297]: the chordless submap has one face");

    DeterministicRandomGenerator rng;
    for (std::size_t v = 0; v < curve.num_vertices(); ++v) {
        check_against_brute(rs, curve, curve.vertex(v), LEFT);
        check_against_brute(rs, curve, curve.vertex(v), RIGHT);
    }
    for (int i = 0; i < 200; ++i) {
        Point p{rng.uniform(-3.0, 27.0), rng.uniform(-2.0, 28.0), SOS_NONE};
        check_against_brute(rs, curve, p, (i & 1) ? LEFT : RIGHT);
    }

    RayHit h = rs.shoot_toward_boundary(Point{5, 100, SOS_NONE}, LEFT);
    assert(!h.hit && "[C91 §2.1 tex 70]: no contact above C's y-range");

    std::printf("  [PASS] structure_trivial_mu1\n");
}

static void test_structure_comb() {
    GranularComb gc;
    const Submap& submap = gc.fx.submap;
    const Polygon& curve = gc.fx.curve;

    RayShootingStructure rs(submap, curve, gc.granularity);

    std::size_t live = submap.num_live_nodes();
    std::size_t nulls = 0;
    for (std::size_t ci = 0; ci < submap.num_chords(); ++ci)
        if (!submap.chord(ci).dead && submap.chord(ci).is_null_length)
            ++nulls;
    assert(rs.num_faces() == live - nulls &&
           "[C91 §3.4 tex 286]: empty regions have no associated faces");
    assert(rs.num_faces() > 1 && "the granular comb retains structure");

    assert(rs.vertical_line().size() == 1 &&
           "[C91 §3.4 tex 306]: exactly the apex outside-pair chord "
           "crosses the vertical line");
    {
        const auto& lc = rs.vertical_line()[0];
        assert(symbolic_y_equal(lc.y, SymbolicY{26.0, 7}));
        assert(submap.region_weight(lc.region_above) == 0 &&
               "the region above the topmost crossing is the polar cap");
    }

    DeterministicRandomGenerator rng;
    for (std::size_t v = 0; v < curve.num_vertices(); ++v) {
        check_against_brute(rs, curve, curve.vertex(v), LEFT);
        check_against_brute(rs, curve, curve.vertex(v), RIGHT);
    }
    for (int i = 0; i < 500; ++i) {
        Point p{rng.uniform(-3.0, 27.0), rng.uniform(-2.0, 28.0), SOS_NONE};
        check_against_brute(rs, curve, p, (i & 1) ? LEFT : RIGHT);
    }

    std::printf("  [PASS] structure_comb (mu=%zu, gamma=%zu)\n", rs.num_faces(), gc.granularity);
}

static void test_structure_separator_recursion() {
    RecursiveFixture cc;
    const Submap& submap = cc.submap;
    const Polygon& curve = cc.curve;

    submap.check_invariants(curve);
    assert(submap.is_conformal());
    assert(submap.is_granular(cc.granularity, curve));

    RayShootingStructure rs(submap, curve, cc.granularity);

    assert(rs.num_faces() > 4 && "the conformal comb keeps many regions (μ > μ^{2/3} leaf "
                                 "threshold)");
    assert(rs.decomposition().dstar_size >= 1 &&
           "[C91 §3.4 tex 304]: a non-leaf μ forces a nonempty D*");
    assert(rs.decomposition().num_subsets >= 2 &&
           "[C91 §3.4 tex 304]: the separator partitions G into ≥ 2 D_i");

    DeterministicRandomGenerator rng(0xC0FFEEu);
    for (std::size_t v = 0; v < curve.num_vertices(); ++v) {
        check_against_brute(rs, curve, curve.vertex(v), LEFT);
        check_against_brute(rs, curve, curve.vertex(v), RIGHT);
    }
    for (int i = 0; i < 3000; ++i) {
        Point p{rng.uniform(-3.0, 131.0), rng.uniform(-2.0, 23.0), SOS_NONE};
        check_against_brute(rs, curve, p, (i & 1) ? LEFT : RIGHT);
    }

    std::printf("  [PASS] structure_separator_recursion "
                "(mu=%zu, |D*|=%zu, subsets=%zu)\n",
                rs.num_faces(), rs.decomposition().dstar_size, rs.decomposition().num_subsets);
}

static std::vector<std::size_t> region_table_arcs(const Submap& submap, std::size_t r) {
    std::vector<std::size_t> out;
    for (std::size_t ai = 0; ai < submap.num_arcs(); ++ai)
        if (!submap.arc(ai).dead && submap.arc(ai).region_node == r)
            out.push_back(ai);
    return out;
}

static void test_collect_region_arcs_wrap_straddle() {
    GranularComb cc;
    const Submap& submap = cc.fx.submap;

    bool saw_wrap_region = false;
    for (std::size_t r = 0; r < submap.num_nodes(); ++r) {
        if (submap.node(r).dead)
            continue;
        auto want = region_table_arcs(submap, r);
        if (want.empty())
            continue;

        assert(want.size() <= 4 && "[C91 §2.3 tex 114]: conformal region has ≤ 4 "
                                   "arc-structures — wrap-spanning arcs are never split");
        for (std::size_t ai : want)
            if (submap.arc(ai).wraps())
                saw_wrap_region = true;

        RegionArcs got = collect_region_arcs(submap, r);
        assert(got.count == want.size() && "[C91 §3.1 tex 181]: collect_region_arcs must gather "
                                           "every arc-structure of the region");
        for (std::size_t ai : want) {
            bool found = false;
            for (std::size_t g : got)
                if (g == ai) {
                    found = true;
                    break;
                }
            assert(found && "[C91 §3.1 tex 181]: collect_region_arcs reaches "
                            "every arc-structure through chord adjacency");
        }
    }
    assert(saw_wrap_region && "[C91 §2.4 tex 142]: the conformal comb must contain a "
                              "region with a double-backing arc (else this test proves "
                              "nothing about the single-structure representation)");

    std::printf("  [PASS] collect_region_arcs_wrap_straddle\n");
}

static void test_region_weight_wrap_straddle() {
    GranularComb cc;
    const Submap& submap = cc.fx.submap;

    bool saw_wrap = false;
    for (std::size_t r = 0; r < submap.num_nodes(); ++r) {
        if (submap.node(r).dead)
            continue;
        auto arcs = region_table_arcs(submap, r);
        std::size_t truth = 0;
        for (std::size_t ai : arcs) {
            truth = std::max(truth, submap.arc(ai).edge_count);
            if (submap.arc(ai).wraps())
                saw_wrap = true;
        }
        assert(submap.region_weight(r) == truth &&
               "[C91 §2.2 tex 106]: region_weight must equal the max "
               "edge_count over ALL the region's arc-structures — "
               "including double-backing ones");
    }
    assert(saw_wrap && "the conformal comb must contain a double-backing arc");

    std::printf("  [PASS] region_weight_wrap_straddle\n");
}

static void test_structure_no_wrapped_chords() {
    Polygon input_curve({{0, 0, 0}, {2, 3, 1}, {4, 5, 2}, {6, 2, 3}, {8, 9, 4}});
    Polygon first_curve = input_curve.subchain(0, 3);
    Polygon second_curve = input_curve.subchain(2, 3);
    Submap first_submap = make_chordless_normal(first_curve);
    Submap second_submap = make_chordless_normal(second_curve);

    const std::size_t g = 2 * std::max(first_curve.num_edges(), second_curve.num_edges());
    TestArcRayShooter r1(first_submap, first_curve, g), r2(second_submap, second_curve, g);
    SafeArcCutter c1(&first_curve, 64), c2(&second_curve, 64);

    MergeInput in;
    in.first_curve = &first_curve;
    in.second_curve = &second_curve;
    in.first_submap = &first_submap;
    in.second_submap = &second_submap;
    in.first_granularity = g;
    in.second_granularity = g;
    in.granularity = g;
    in.first_ray_shooter = &r1;
    in.second_ray_shooter = &r2;
    in.first_arc_cutter = &c1;
    in.second_arc_cutter = &c2;
    in.first_piece_count_bound = 64;
    in.second_piece_count_bound = 64;
    in.first_piece_granularity_bound = 64;
    in.second_piece_granularity_bound = 64;
    MergeResult res = merge(in);

    std::size_t granularity = 0;
    for (std::size_t r = 0; r < res.submap.num_nodes(); ++r)
        if (!res.submap.node(r).dead)
            granularity = std::max(granularity, res.submap.region_weight(r));

    RayShootingStructure rs(res.submap, res.curve, granularity);
    assert(rs.num_faces() > 1 && "the merge must retain structure");
    assert(rs.vertical_line().empty() && "[C91 §2.1 tex 70]: no chord wraps ⟹ empty vertical line");
    assert(rs.region_at_infinity() != NONE && rs.region_at_infinity() < res.submap.num_nodes() &&
           !res.submap.node(rs.region_at_infinity()).dead &&
           "[C91 §3.4 tex 306]: a live polar region is elected when no "
           "chord crosses the vertical line");

    DeterministicRandomGenerator rng(0x5EEDu);
    for (std::size_t v = 0; v < res.curve.num_vertices(); ++v) {
        check_against_brute(rs, res.curve, res.curve.vertex(v), LEFT);
        check_against_brute(rs, res.curve, res.curve.vertex(v), RIGHT);
    }
    for (int i = 0; i < 2000; ++i) {
        Point p{rng.uniform(-2.0, 10.0), rng.uniform(-2.0, 11.0), SOS_NONE};
        check_against_brute(rs, res.curve, p, (i & 1) ? LEFT : RIGHT);
    }

    std::printf("  [PASS] structure_no_wrapped_chords (mu=%zu, "
                "region_infinity path)\n",
                rs.num_faces());
}

static void test_local_min_junction() {
    Polygon input_curve({{0, 1, 0}, {2, 3, 1}, {4, 2, 2}, {6, 4, 3}, {8, 7, 4}});
    Polygon first_curve = input_curve.subchain(0, 3);
    Polygon second_curve = input_curve.subchain(2, 3);
    Submap first_submap = make_chordless_normal(first_curve);
    Submap second_submap = make_chordless_normal(second_curve);

    const std::size_t g = 2 * std::max(first_curve.num_edges(), second_curve.num_edges());
    TestArcRayShooter r1(first_submap, first_curve, g), r2(second_submap, second_curve, g);
    SafeArcCutter c1(&first_curve, 64), c2(&second_curve, 64);

    MergeInput in;
    in.first_curve = &first_curve;
    in.second_curve = &second_curve;
    in.first_submap = &first_submap;
    in.second_submap = &second_submap;
    in.first_granularity = g;
    in.second_granularity = g;
    in.granularity = g;
    in.first_ray_shooter = &r1;
    in.second_ray_shooter = &r2;
    in.first_arc_cutter = &c1;
    in.second_arc_cutter = &c2;
    in.first_piece_count_bound = 64;
    in.second_piece_count_bound = 64;
    in.first_piece_granularity_bound = 64;
    in.second_piece_granularity_bound = 64;
    MergeResult res = merge(in);

    res.submap.check_invariants(res.curve);
    assert(res.submap.is_conformal());
    assert(res.submap.is_granular(g, res.curve));

    std::size_t granularity = 0;
    for (std::size_t r = 0; r < res.submap.num_nodes(); ++r)
        if (!res.submap.node(r).dead)
            granularity = std::max(granularity, res.submap.region_weight(r));
    RayShootingStructure rs(res.submap, res.curve, granularity);

    DeterministicRandomGenerator rng(0xA11CEu);
    for (std::size_t v = 0; v < res.curve.num_vertices(); ++v) {
        check_against_brute(rs, res.curve, res.curve.vertex(v), LEFT);
        check_against_brute(rs, res.curve, res.curve.vertex(v), RIGHT);
    }
    for (int i = 0; i < 2000; ++i) {
        Point p{rng.uniform(-2.0, 10.0), rng.uniform(-1.0, 9.0), SOS_NONE};
        check_against_brute(rs, res.curve, p, (i & 1) ? LEFT : RIGHT);
    }

    std::printf("  [PASS] local_min_junction\n");
}

static void test_target_only_oracle() {
    Polygon first_curve = make_C1();
    Submap first_submap = make_chordless_normal(first_curve);
    TestArcRayShooter shooter(first_submap, first_curve, 2 * first_curve.num_edges());

    const Point& p = first_curve.vertex(4);
    Subarc right_whole{6,
                       RIGHT,
                       0,
                       RIGHT,
                       symbolic_y_of(first_curve.vertex(7)),
                       symbolic_y_of(first_curve.vertex(0))};
    RayHit h = shooter.shoot(p, LEFT, 0, right_whole);
    assert(h.hit && !h.wrapped && h.edge == 0 && h.side == RIGHT && h.x == Exact{2} / 5);

    Subarc left_whole{0,
                      LEFT,
                      6,
                      LEFT,
                      symbolic_y_of(first_curve.vertex(0)),
                      symbolic_y_of(first_curve.vertex(7))};
    h = shooter.shoot(p, LEFT, 0, left_whole);
    assert(h.hit && !h.wrapped && h.edge == 0 && h.side == RIGHT && h.x == Exact{2} / 5);

    Polygon curve({{0, 0, 0}, {2, 4, 1}, {4, 1, 2}, {6, 5, 3}});
    TestArcRayShooter target_shooter(&curve);
    Subarc target{2, LEFT, 2, LEFT, symbolic_y_of(curve.vertex(2)), symbolic_y_of(curve.vertex(3))};
    h = target_shooter.shoot(Point{-1, 2, SOS_NONE}, RIGHT, 0, target);
    assert(h.hit && !h.wrapped && h.edge == 2 && h.side == LEFT && h.x == 4.5);

    const Point& apex_w = first_curve.vertex(1);
    h = shooter.shoot(apex_w, LEFT, 0, right_whole);
    assert(h.hit && h.wrapped && h.edge == 6 && h.side == RIGHT);

    std::printf("  [PASS] target_only_oracle\n");
}

static void test_fusion_wrapped_junction() {
    Polygon first_curve = make_C1(), second_curve = make_C2();
    Submap first_submap = make_chordless_normal(first_curve);
    Submap second_submap = make_chordless_normal(second_curve);

    TestArcRayShooter first_ray_shooter(first_submap, first_curve, 14),
        second_ray_shooter(second_submap, second_curve, 10);

    FusionState st2;
    st2.junction_at_end = false;
    fuse_submaps(st2, second_submap, second_curve, first_submap, first_curve, second_ray_shooter,
                 first_ray_shooter);
    bool found_apex = false;
    for (const auto& dc : st2.chords) {
        if (!symbolic_y_equal(dc.y, SymbolicY{26.0, 7}))
            continue;

        std::size_t le = dc.left_on_first_curve ? dc.left_edge + 7 : dc.left_edge;
        std::size_t re = dc.right_on_first_curve ? dc.right_edge + 7 : dc.right_edge;
        if ((le == 6 && re == 7) || (le == 7 && re == 6))
            found_apex = true;
    }
    assert(found_apex && "[C91 §2.1 tex 70] + [C91 §3.1 tex 191]: the wrapped startup "
                         "discovers the apex outside-pair chord");

    std::printf("  [PASS] fusion_wrapped_junction\n");
}

static MergeResult run_comb_merge(std::size_t granularity) {
    Polygon first_curve = make_C1(), second_curve = make_C2();
    static std::vector<std::unique_ptr<Polygon>> keep_p;
    static std::vector<std::unique_ptr<Submap>> keep_s;
    keep_p.push_back(std::make_unique<Polygon>(first_curve));
    keep_p.push_back(std::make_unique<Polygon>(second_curve));
    keep_s.push_back(std::make_unique<Submap>(make_chordless_normal(*keep_p[keep_p.size() - 2])));
    keep_s.push_back(std::make_unique<Submap>(make_chordless_normal(*keep_p[keep_p.size() - 1])));
    const Polygon& c1 = *keep_p[keep_p.size() - 2];
    const Polygon& c2 = *keep_p[keep_p.size() - 1];
    Submap& s1 = *keep_s[keep_s.size() - 2];
    Submap& s2 = *keep_s[keep_s.size() - 1];

    static std::vector<std::unique_ptr<TestArcRayShooter>> keep_r;
    static std::vector<std::unique_ptr<SafeArcCutter>> keep_c;

    keep_r.push_back(std::make_unique<TestArcRayShooter>(s1, c1, 14));
    keep_r.push_back(std::make_unique<TestArcRayShooter>(s2, c2, 14));
    keep_c.push_back(std::make_unique<SafeArcCutter>(&c1, 64));
    keep_c.push_back(std::make_unique<SafeArcCutter>(&c2, 64));

    MergeInput in;
    in.first_curve = &c1;
    in.second_curve = &c2;
    in.first_submap = &s1;
    in.second_submap = &s2;
    in.first_granularity = 14;
    in.second_granularity = 14;
    in.granularity = granularity;
    in.first_ray_shooter = keep_r[keep_r.size() - 2].get();
    in.second_ray_shooter = keep_r[keep_r.size() - 1].get();
    in.first_arc_cutter = keep_c[keep_c.size() - 2].get();
    in.second_arc_cutter = keep_c[keep_c.size() - 1].get();
    in.first_piece_count_bound = 64;
    in.second_piece_count_bound = 64;
    in.first_piece_granularity_bound = 64;
    in.second_piece_granularity_bound = 64;
    return merge(in);
}

static void test_merge_comb_e2e() {
    {
        MergeResult r = run_comb_merge(24);
        r.submap.check_invariants(r.curve);
        assert(r.submap.is_conformal());
        assert(r.submap.is_granular(24, r.curve));
        assert(r.submap.num_live_chords() == 0 && "γ ≥ total weight ⟹ γ-granular ⟺ chordless");
        assert(r.submap.num_live_nodes() == 1);
    }

    {
        MergeResult r = run_comb_merge(14);
        r.submap.check_invariants(r.curve);
        assert(r.submap.is_conformal() && "[C91 Lemma 3.5 tex 279]: merge output is conformal");
        assert(r.submap.is_granular(14, r.curve) &&
               "[C91 Lemma 3.5 tex 279]: merge output is γ-granular");
        assert(!r.submap.tree_decomposition().empty() &&
               "[C91 §2.4(iv)]: normal form includes the tree "
               "decomposition");
        assert(r.submap.num_live_chords() >= 1 &&
               "γ = 14 < total weight forces surviving exit chords");

        for (std::size_t rg = 0; rg < r.submap.num_nodes(); ++rg)
            if (!r.submap.node(rg).dead)
                assert(r.submap.region_weight(rg) <= 14);
    }

    std::printf("  [PASS] merge_comb_e2e\n");
}

static void test_structure_preconditions() {
    require_assertion_abort([] {
        CombFixture fx;
        RayShootingStructure rs(fx.submap, fx.curve, 100);
        (void)rs;
    });

    require_assertion_abort([] {
        Polygon curve({{0, 0, 0}, {0, 4, 1}, {4, 5, 2}, {4, -1, 3}});
        PendingChord ch{{2.0, 42}, 0, RIGHT, 2, RIGHT, false};
        Submap submap;
        build_submap_from_chords(submap, curve, {ch});
        RayShootingStructure rs(submap, curve, 100);
        (void)rs;
    });

    require_assertion_abort([] {
        Polygon curve = make_C();
        Submap submap = make_chordless_normal(curve);
        RayShootingStructure rs(submap, curve, 1);
        (void)rs;
    });
    std::printf("  [PASS] structure_preconditions\n");
}

int main() {
    std::setbuf(stdout, nullptr);
    std::printf("[C91 §3.4 tests]:\n");
    test_structure_trivial_mu1();
    test_structure_comb();
    test_structure_separator_recursion();
    test_collect_region_arcs_wrap_straddle();
    test_region_weight_wrap_straddle();
    test_structure_no_wrapped_chords();
    test_local_min_junction();
    test_target_only_oracle();
    test_fusion_wrapped_junction();
    test_merge_comb_e2e();
    test_structure_preconditions();
    std::printf("All §3.4 tests passed.\n");
    return 0;
}
