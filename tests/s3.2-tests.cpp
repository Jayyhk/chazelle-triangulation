#include "merge/conformality.h"
#include "merge/fusion.h"
#include "polygon/polygon.h"
#include "submap/submap.h"
#include "support/arc_ray_shooter.h"
#include "support/assertions.h"

#include <cassert>
#include <cstdio>
#include <memory>

using namespace chazelle;

using chazelle::test::require_assertion_abort;

static Polygon make_C() {
    return Polygon({{0, 0, 0},
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
}
static Polygon make_C1() {
    return Polygon({{0, 0, 0},
                    {2, 20, 1},
                    {4, 6, 2},
                    {6, 24, 3},
                    {8, 4, 4},
                    {10, 22, 5},
                    {12, 5, 6},
                    {14, 26, 7}});
}
static Polygon make_C2() {
    return Polygon(
        {{14, 26, 7}, {16, 2, 8}, {18, 23, 9}, {20, 4.5, 10}, {22, 25, 11}, {24, 1, 12}});
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

using GeomRayShooter = TestArcRayShooter;

struct GeomArcCutter : ArcCuttingOracle {
    const Polygon* input_curve;
    std::size_t piece_len;
    mutable std::vector<std::unique_ptr<Polygon>> curves;
    mutable std::vector<std::unique_ptr<Submap>> submaps;

    GeomArcCutter(const Polygon* ci, std::size_t pl) : input_curve(ci), piece_len(pl) {}

    Submap build_piece_submap(const Polygon& alpha) const {
        return make_chordless_normal(alpha);
    }

    std::vector<ArcPiece> cut(std::size_t arc_idx, const Subarc& target) const override {
        if (target.first_side != target.last_side ||
            (target.first_side == LEFT && target.first_edge > target.last_edge) ||
            (target.first_side == RIGHT && target.first_edge < target.last_edge)) {
            ArcSideRange ranges[3];
            std::size_t nl = subarc_side_ranges(target, 0, input_curve->num_vertices() - 1, ranges);
            std::vector<ArcPiece> out;
            for (std::size_t g = 0; g < nl; ++g) {
                Subarc subarc_on_side =
                    (ranges[g].side == LEFT)
                        ? Subarc{ranges[g].first_edge, LEFT, ranges[g].last_edge, LEFT}
                        : Subarc{ranges[g].last_edge, RIGHT, ranges[g].first_edge, RIGHT};
                subarc_on_side =
                    test_full_subarc(*input_curve, subarc_on_side.first_edge,
                                     subarc_on_side.first_side, subarc_on_side.last_edge);
                if (g == 0)
                    subarc_on_side.first_y = target.first_y;
                if (g + 1 == nl)
                    subarc_on_side.last_y = target.last_y;
                auto pieces = cut(arc_idx, subarc_on_side);
                out.insert(out.end(), pieces.begin(), pieces.end());
            }
            test_set_cut_endpoints(*input_curve, target, out);
            return out;
        }
        const Side s = target.first_side;
        const std::size_t lo = std::min(target.first_edge, target.last_edge);
        const std::size_t hi = std::max(target.first_edge, target.last_edge);

        const std::size_t first_e = target.first_edge;
        const std::size_t last_e = target.last_edge;

        const bool first_partial = !symbolic_y_equal(
            target.first_y, symbolic_y_of(input_curve->vertex(first_e + (s == RIGHT))));
        const bool last_partial = !symbolic_y_equal(
            target.last_y, symbolic_y_of(input_curve->vertex(last_e + (s == LEFT))));
        const bool low_partial = (s == LEFT) ? first_partial : last_partial;
        const bool high_partial = (s == LEFT) ? last_partial : first_partial;

        std::size_t mlo = lo + (low_partial ? 1 : 0);
        std::size_t mhi_p1 = hi + 1 - (high_partial ? 1 : 0);

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
            submaps.push_back(std::make_unique<Submap>(build_piece_submap(*curves.back())));
            ArcPiece p;
            p.subarc = (s == LEFT) ? Subarc{a, s, b, s} : Subarc{b, s, a, s};
            p.curve = curves.back().get();
            p.submap = submaps.back().get();

            p.granularity = 2 * (b - a + 1);
            return p;
        };

        std::vector<std::size_t> chunk_lo, chunk_hi;
        for (std::size_t a = mlo; a < mhi_p1; a += piece_len) {
            chunk_lo.push_back(a);
            chunk_hi.push_back(std::min(a + piece_len - 1, mhi_p1 - 1));
        }

        std::vector<ArcPiece> out;
        if (lo == hi && low_partial && high_partial) {
            out.push_back(boundary_piece(lo));
            test_set_cut_endpoints(*input_curve, target, out);
            return out;
        }
        if (s == LEFT) {
            if (low_partial)
                out.push_back(boundary_piece(lo));
            for (std::size_t k = 0; k < chunk_lo.size(); ++k)
                out.push_back(middle_piece(chunk_lo[k], chunk_hi[k]));
            if (high_partial)
                out.push_back(boundary_piece(hi));
        } else {
            if (high_partial)
                out.push_back(boundary_piece(hi));
            for (std::size_t k = chunk_lo.size(); k-- > 0;)
                out.push_back(middle_piece(chunk_lo[k], chunk_hi[k]));
            if (low_partial)
                out.push_back(boundary_piece(lo));
        }
        test_set_cut_endpoints(*input_curve, target, out);
        return out;
    }
};

static void test_cutter_wrapped_partial_endpoints() {
    Polygon c({{0, 0, 0}, {1, 2, 1}, {2, 4, 2}});
    GeomArcCutter cutter(&c, 2);
    const Subarc targets[] = {{0, LEFT, 1, RIGHT, {1, 42}, {3, 43}},
                              {1, RIGHT, 0, LEFT, {3, 42}, {1, 43}},
                              {1, LEFT, 0, LEFT, {3, 42}, {1, 43}}};
    for (const Subarc& target : targets) {
        const auto pieces = cutter.cut(0, target);
        assert_cut_postconditions(c, target, pieces.data(), pieces.size(), 6, 4);
        assert(pieces.size() == 3 && pieces.front().is_boundary_piece &&
               pieces.back().is_boundary_piece && !pieces[1].is_boundary_piece);
    }
    std::printf("  [PASS] cutter_wrapped_partial_endpoints\n");
}

struct CombRig {
    CombFixture fx;
    GeomRayShooter first_ray_shooter, second_ray_shooter;
    GeomArcCutter first_arc_cutter, second_arc_cutter;
    ConformalityOracles oracles;

    explicit CombRig(std::size_t piece_len = 1)
        : first_ray_shooter(&fx.first_curve), second_ray_shooter(&fx.second_curve),
          first_arc_cutter(&fx.first_curve, piece_len),
          second_arc_cutter(&fx.second_curve, piece_len) {
        oracles.first_submap = &fx.first_submap;
        oracles.second_submap = &fx.second_submap;
        oracles.first_curve = &fx.first_curve;
        oracles.second_curve = &fx.second_curve;
        oracles.first_ray_shooter = &first_ray_shooter;
        oracles.second_ray_shooter = &second_ray_shooter;
        oracles.first_arc_cutter = &first_arc_cutter;
        oracles.second_arc_cutter = &second_arc_cutter;
        oracles.first_piece_count_bound = 64;
        oracles.second_piece_count_bound = 64;
        oracles.first_piece_granularity_bound = 64;
        oracles.second_piece_granularity_bound = 64;
    }
};

static std::vector<std::size_t> region_arcs_of(const Submap& submap, std::size_t region) {
    std::vector<std::size_t> out;
    for (std::size_t ai = 0; ai < submap.num_arcs(); ++ai)
        if (!submap.arc(ai).dead && submap.arc(ai).region_node == region)
            out.push_back(ai);
    return out;
}

static void assert_chord_geometrically_valid(const Polygon& curve, const Chord& c) {
    SymbolicY y = c.symbolic_y();
    Exact x1 = edge_x_at_y(curve, c.left_edge, y);
    Exact x2 = edge_x_at_y(curve, c.right_edge, y);
    const bool wraps = chord_runs_through_infinity(curve, c);
    if (x1 == x2 && !wraps)
        return;
    Exact lo = std::min(x1, x2);
    Exact hi = std::max(x1, x2);
    for (std::size_t e = 0; e < curve.num_edges(); ++e) {
        Exact x;
        if (!GeomRayShooter::crossing_x(curve, e, y, &x))
            continue;
        if (wraps)
            assert(!(x < lo || x > hi) && "[C91 §2.1 tex 74]: a wrapping chord's through-infinity "
                                          "segment must not cross C (mutual visibility)");
        else
            assert(!(x > lo && x < hi) && "[C91 §2.1 tex 74]: a direct chord's open interior must "
                                          "not cross C (mutual visibility)");
    }
}

static void test_fixture_and_arc_end_y() {
    CombFixture fx;
    fx.submap.check_invariants(fx.curve);
    fx.submap.assert_tree_property();
    assert(!fx.submap.is_conformal() && "R_big has degree 7 — not conformal yet");

    SymbolicY e5 = fx.submap.arc_end_symbolic_y(5, fx.curve);
    assert(symbolic_y_equal(e5, SymbolicY{5.0, 6}));

    SymbolicY e15 = fx.submap.arc_end_symbolic_y(15, fx.curve);
    assert(symbolic_y_equal(e15, SymbolicY{26.0, 7}));

    SymbolicY s15 = fx.submap.arc_start_symbolic_y(15, fx.curve);
    assert(symbolic_y_equal(s15, SymbolicY{6.0, 2}));

    SymbolicY e7 = fx.submap.arc_end_symbolic_y(7, fx.curve);
    assert(symbolic_y_equal(e7, SymbolicY{26.0, 7}));

    std::printf("  [PASS] fixture_and_arc_end_y\n");
}

static void test_arc_sources() {
    CombFixture fx;
    auto prov = identify_arc_sources(fx.submap, fx.curve, fx.first_submap, fx.first_curve,
                                     fx.second_submap, fx.second_curve);
    assert(prov.size() == 16);

    auto expect = [&](std::size_t ai, bool on_first_curve, std::size_t si_arc) {
        assert(prov[ai].on_first_curve == on_first_curve);
        assert(prov[ai].input_arc == si_arc);
    };
    expect(0, false, 0);
    expect(1, false, 0);
    expect(5, false, 0);
    expect(7, false, 0);
    expect(8, true, 0);
    expect(11, true, 0);
    expect(15, true, 0);

    std::printf("  [PASS] arc_sources\n");
}

static void test_fused_region_cycle() {
    CombFixture fx;

    auto cycle = fused_region_cycle(fx.submap, fx.curve, 0, region_arcs_of(fx.submap, 0));
    assert(cycle.count == 7 && "[C91 §3.2 tex 238]: the fused big region has 7 arcs");

    assert(cycle.arcs[0].arc == 1 && !cycle.arcs[0].is_zero_length);
    assert(cycle.arcs[1].arc == 3 && cycle.arcs[1].is_zero_length);
    assert(cycle.arcs[2].arc == 5 && !cycle.arcs[2].is_zero_length);
    assert(cycle.arcs[3].arc == 9 && cycle.arcs[3].is_zero_length);
    assert(cycle.arcs[4].arc == 11);
    assert(cycle.arcs[5].arc == 13 && cycle.arcs[5].is_zero_length);
    assert(cycle.arcs[6].arc == 15 && !cycle.arcs[6].is_zero_length);

    auto p4 = fused_region_cycle(fx.submap, fx.curve, 4, region_arcs_of(fx.submap, 4));
    assert(p4.count == 2);

    auto az = fused_region_cycle(fx.submap, fx.curve, 7, region_arcs_of(fx.submap, 7));
    assert(az.count == 1 && az.arcs[0].is_zero_length);

    std::printf("  [PASS] fused_region_cycle\n");
}

static void test_local_shoot_fused() {
    CombRig rig;
    auto& fx = rig.fx;
    auto prov = identify_arc_sources(fx.submap, fx.curve, fx.first_submap, fx.first_curve,
                                     fx.second_submap, fx.second_curve);
    auto cycle = fused_region_cycle(fx.submap, fx.curve, 0, region_arcs_of(fx.submap, 0));

    FusedShootContext ctx;
    ctx.submap = &fx.submap;
    ctx.curve = &fx.curve;
    ctx.first_curve = &fx.first_curve;
    ctx.second_curve = &fx.second_curve;
    ctx.first_ray_shooter = &rig.first_ray_shooter;
    ctx.second_ray_shooter = &rig.second_ray_shooter;
    ctx.arc_sources = &prov;

    RayHit h = local_shoot_fused(Point{8, 4, 4}, SymbolicY{4.0, 4}, LEFT, cycle, ctx);
    assert(h.hit && !h.wrapped);
    assert(h.hit_arc_idx == 15);
    assert(h.edge == 0 && h.side == RIGHT);
    assert(h.x == Exact{2} / 5);

    h = local_shoot_fused(Point{16, 2, 8}, SymbolicY{2.0, 8}, RIGHT, cycle, ctx);
    assert(h.hit && !h.wrapped);
    assert(h.hit_arc_idx == 1);
    assert(h.edge == 11 && h.side == RIGHT);

    h = local_shoot_fused(Point{2, 20, 1}, SymbolicY{20.0, 1}, LEFT, cycle, ctx);
    assert(h.hit && h.wrapped);
    assert(h.hit_arc_idx == 1);
    assert(h.edge == 11 && h.side == LEFT);

    std::printf("  [PASS] local_shoot_fused\n");
}

static void test_insert_chord() {
    CombFixture fx;

    std::size_t flat[] = {1, 3, 5, 9, 11, 13, 15};
    Submap::ChordPointSpec p{11, 3, RIGHT, 8.0};
    Submap::ChordPointSpec q{15, 0, RIGHT, Exact{2} / 5};
    auto res = fx.submap.insert_chord(p, q, SymbolicY{4.0, 4}, 0, flat, 7, fx.curve);

    assert(res.chord_idx == 8);
    assert(res.new_region == 9);
    assert(fx.submap.num_nodes() == 10);
    assert(fx.submap.num_chords() == 9);
    assert(fx.submap.num_arcs() == 18);
    fx.submap.assert_tree_property();

    const Chord& nc = fx.submap.chord(8);

    assert(nc.left_edge == 0 && nc.left_side == RIGHT);
    assert(nc.left_adj.count == 2);
    assert(nc.left_adj.arcs[0] == 15 && nc.left_adj.arcs[1] == res.q_after_arc);
    assert(nc.right_edge == 3 && nc.right_side == RIGHT);
    assert(nc.right_adj.count == 1 && nc.right_adj.arcs[0] == 11);
    assert(nc.region[0] == 0 && nc.region[1] == 9);

    assert(fx.submap.arc(res.p_after_arc).region_node == 9);
    assert(fx.submap.arc(13).region_node == 9);
    assert(fx.submap.arc(15).region_node == 9);
    assert(fx.submap.arc(res.q_after_arc).region_node == 0);
    assert(fx.submap.arc(11).region_node == 0);

    assert(fx.submap.arc(11).last_edge == 4 && fx.submap.arc(11).first_edge == 4);
    assert(fx.submap.arc(11).edge_count == 1);
    assert(fx.submap.arc(res.p_after_arc).first_edge == 3 &&
           fx.submap.arc(res.p_after_arc).last_edge == 3);
    assert(fx.submap.arc(res.p_after_arc).edge_count == 1);

    assert(fx.submap.arc(15).first_side == RIGHT && fx.submap.arc(15).last_side == RIGHT &&
           fx.submap.arc(15).edge_count == 1);
    assert(fx.submap.arc(res.q_after_arc).first_side == RIGHT &&
           fx.submap.arc(res.q_after_arc).last_side == LEFT &&
           fx.submap.arc(res.q_after_arc).edge_count == 8);
    assert(fx.submap.start_arc == res.q_after_arc);

    assert(fx.submap.chord(1).right_adj.arcs[0] == res.p_after_arc);
    assert(fx.submap.chord(1).region[0] == 9);

    assert(fx.submap.chord(2).left_adj.arcs[1] == 11);

    assert(fx.submap.chord(0).region[0] == 9);

    auto c9 = fused_region_cycle(fx.submap, fx.curve, 9, region_arcs_of(fx.submap, 9));
    assert(c9.count == 3);
    auto c0 = fused_region_cycle(fx.submap, fx.curve, 0, region_arcs_of(fx.submap, 0));
    assert(c0.count == 6);

    std::printf("  [PASS] insert_chord\n");
}

static void test_insert_chord_equal_x_wrap() {
    Polygon curve({{0, 0, 0}, {4, 10, 1}, {2, 20, 2}, {6, 30, 3}});

    Submap submap;
    for (int i = 0; i < 3; ++i)
        submap.add_node();
    submap.start_vertex = 0;
    submap.end_vertex = 3;

    auto add = [&](std::size_t fe, Side fs, std::size_t le, Side ls, std::size_t region,
                   std::size_t count) {
        Arc a{};
        a.first_edge = fe;
        a.first_side = fs;
        a.last_edge = le;
        a.last_side = ls;
        a.region_node = region;
        a.edge_count = count;
        return submap.add_arc(a);
    };

    add(1, LEFT, 2, LEFT, 0, 2);
    add(2, LEFT, 2, RIGHT, 1, 1);
    add(2, RIGHT, 1, RIGHT, 0, 2);
    add(0, RIGHT, 0, LEFT, 2, 1);
    assert(submap.end_arc == 1 && submap.start_arc == 3);

    auto chord = [&](std::size_t le, Side lsd, Chord::AdjArcs ladj, std::size_t re, Side rsd,
                     Chord::AdjArcs radj, const Exact& y, std::size_t tag, std::size_t r0,
                     std::size_t r1) {
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
        c.is_null_length = false;
        submap.add_chord(c);
    };
    chord(2, LEFT, {{0, 1}, 2}, 2, RIGHT, {{1, 2}, 2}, 25.0, 99, 1, 0);
    chord(0, LEFT, {{3}, 1}, 1, RIGHT, {{2}, 1}, 10.0, 1, 2, 0);

    std::size_t cyc[] = {0, 2};
    Submap::ChordPointSpec p{0, 1, LEFT, 2.0};
    Submap::ChordPointSpec q{2, 1, RIGHT, 2.0};
    auto res = submap.insert_chord(p, q, SymbolicY{20.0, 2}, 0, cyc, 2, curve);

    assert(submap.num_chords() == 3);
    submap.assert_tree_property();

    const Chord& nc = submap.chord(res.chord_idx);

    assert(nc.left_edge == 1 && nc.left_side == LEFT);
    assert(nc.right_edge == 1 && nc.right_side == RIGHT);

    assert(nc.left_adj.count == 1 && nc.left_adj.arcs[0] == 0);
    assert(nc.right_adj.count == 1 && nc.right_adj.arcs[0] == 2);

    assert(chord_runs_through_infinity(curve, nc));

    assert(submap.arc(0).first_edge == 1 && submap.arc(0).last_edge == 1 &&
           submap.arc(0).edge_count == 1);
    assert(submap.arc(res.p_after_arc).first_edge == 2 &&
           submap.arc(res.p_after_arc).last_edge == 2 &&
           submap.arc(res.p_after_arc).edge_count == 1);
    assert(submap.arc(2).first_edge == 2 && submap.arc(2).last_edge == 2 &&
           submap.arc(2).edge_count == 1);
    assert(submap.arc(res.q_after_arc).first_edge == 1 &&
           submap.arc(res.q_after_arc).last_edge == 1 &&
           submap.arc(res.q_after_arc).edge_count == 1);

    assert(submap.arc(res.p_after_arc).region_node == res.new_region);
    assert(submap.arc(2).region_node == res.new_region);
    assert(submap.arc(0).region_node == 0);
    assert(submap.arc(res.q_after_arc).region_node == 0);

    assert(submap.chord(0).region[0] == 1 && submap.chord(0).region[1] == res.new_region);
    assert(submap.chord(0).left_adj.arcs[0] == res.p_after_arc);
    assert(submap.chord(0).right_adj.arcs[1] == 2);

    assert(submap.chord(1).right_adj.arcs[0] == res.q_after_arc);
    assert(submap.chord(1).region[0] == 2 && submap.chord(1).region[1] == 0);

    std::printf("  [PASS] insert_chord_equal_x_wrap\n");
}

static void test_find_visible_point() {
    CombRig rig;
    auto& fx = rig.fx;
    auto prov = identify_arc_sources(fx.submap, fx.curve, fx.first_submap, fx.first_curve,
                                     fx.second_submap, fx.second_curve);
    auto cycle = fused_region_cycle(fx.submap, fx.curve, 0, region_arcs_of(fx.submap, 0));

    VisiblePoint vp = find_visible_point(fx.submap, fx.curve, 0, cycle.arcs[2].arc,
                                         cycle.arcs[0].arc, cycle, prov, rig.oracles);
    assert(vp.found);
    assert(symbolic_y_equal(vp.y, SymbolicY{2.0, 8}));
    assert(vp.p_table_arc == 5 && vp.p_edge == 8 && vp.p_side == RIGHT);
    assert(vp.q_table_arc == 1 && vp.q_edge == 11 && vp.q_side == RIGHT);

    vp = find_visible_point(fx.submap, fx.curve, 0, cycle.arcs[4].arc, cycle.arcs[6].arc, cycle,
                            prov, rig.oracles);
    assert(vp.found);
    assert(symbolic_y_equal(vp.y, SymbolicY{4.0, 4}));
    assert(vp.q_table_arc == 15 && vp.q_x == Exact{2} / 5);

    vp = find_visible_point(fx.submap, fx.curve, 0, cycle.arcs[6].arc, cycle.arcs[2].arc, cycle,
                            prov, rig.oracles);
    assert(!vp.found);

    require_assertion_abort([&] {
        find_visible_point(fx.submap, fx.curve, 0, cycle.arcs[1].arc, cycle.arcs[4].arc, cycle,
                           prov, rig.oracles);
    });

    std::printf("  [PASS] find_visible_point\n");
}

struct DescentCutter : ArcCuttingOracle {
    const Polygon* first_curve;
    mutable std::vector<std::unique_ptr<Polygon>> curves;
    mutable std::vector<std::unique_ptr<Submap>> submaps;
    explicit DescentCutter(const Polygon* c1) : first_curve(c1) {}

    std::vector<ArcPiece> cut(std::size_t, const Subarc& target) const override {
        std::vector<ArcPiece> out;
        assert(target.first_edge == 0 && target.first_side == RIGHT && target.last_edge == 6 &&
               target.last_side == LEFT && "descent test cuts the wrap-spanning S*");

        {
            ArcPiece p;
            p.subarc = Subarc{0, RIGHT, 0, RIGHT};
            p.is_boundary_piece = true;
            out.push_back(p);
        }

        std::vector<Point> vs;
        for (std::size_t v = 0; v <= 7; ++v)
            vs.push_back(first_curve->vertex(v));
        curves.push_back(std::make_unique<Polygon>(std::move(vs)));

        auto sm = std::make_unique<Submap>();
        std::size_t r_out = sm->add_node();
        std::size_t r_pocket = sm->add_node();
        sm->start_vertex = 0;
        sm->end_vertex = 7;

        Arc a{};
        a.first_edge = 3;
        a.last_edge = 6;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = r_pocket;
        a.edge_count = 4;
        std::size_t l_pocket = sm->add_arc(a);
        a = {};
        a.first_edge = 6;
        a.last_edge = 2;
        a.first_side = LEFT;
        a.last_side = LEFT;

        a.region_node = r_out;
        a.edge_count = 11;
        std::size_t w_out = sm->add_arc(a);

        Chord c{};
        c.region[0] = r_out;
        c.region[1] = r_pocket;
        c.left_edge = 3;
        c.left_side = LEFT;
        c.right_edge = 6;
        c.right_side = LEFT;
        c.y = 24.0;
        c.y_tag = 3;
        c.left_adj = {{w_out}, 1};
        c.right_adj = {{l_pocket, w_out}, 2};
        sm->add_chord(c);
        assert(sm->start_arc == w_out && sm->end_arc == w_out &&
               "[C91 §2.4 tex 142]: the double-wrap arc is both "
               "endpoint arcs");
        sm->build_tree_decomposition();
        submaps.push_back(std::move(sm));

        ArcPiece p;
        p.subarc = Subarc{0, LEFT, 6, LEFT};
        p.curve = curves.back().get();
        p.submap = submaps.back().get();

        p.granularity = 11;
        out.push_back(p);
        test_set_cut_endpoints(*first_curve, target, out);
        return out;
    }
};

static void test_descent() {
    CombRig rig;
    auto& fx = rig.fx;
    DescentCutter dcut(&fx.first_curve);
    ConformalityOracles oracles = rig.oracles;
    oracles.first_arc_cutter = &dcut;
    oracles.first_piece_count_bound = 4;

    oracles.first_piece_granularity_bound = 11;

    auto prov = identify_arc_sources(fx.submap, fx.curve, fx.first_submap, fx.first_curve,
                                     fx.second_submap, fx.second_curve);
    auto cycle = fused_region_cycle(fx.submap, fx.curve, 0, region_arcs_of(fx.submap, 0));

    VisiblePoint vp = find_visible_point(fx.submap, fx.curve, 0, cycle.arcs[6].arc,
                                         cycle.arcs[2].arc, cycle, prov, oracles);
    assert(!vp.found);

    vp = find_visible_point(fx.submap, fx.curve, 0, cycle.arcs[6].arc, cycle.arcs[4].arc, cycle,
                            prov, oracles);
    assert(!vp.found);

    std::printf("  [PASS] descent\n");
}

static void run_restore_and_check(std::size_t piece_len) {
    CombRig rig(piece_len);
    auto& fx = rig.fx;
    const std::size_t chords_before = fx.submap.num_chords();
    const std::size_t nodes_before = fx.submap.num_nodes();

    restore_conformality(fx.submap, fx.curve, rig.oracles);

    assert(fx.submap.is_conformal());
    fx.submap.assert_tree_property();

    for (std::size_t r = 0; r < fx.submap.num_nodes(); ++r) {
        if (fx.submap.node(r).dead)
            continue;
        auto inv = region_arcs_of(fx.submap, r);
        if (inv.empty())
            continue;
        auto cyc = fused_region_cycle(fx.submap, fx.curve, r, inv);
        assert(cyc.count <= 4);
    }

    std::size_t added = fx.submap.num_chords() - chords_before;
    assert(added >= 2 && added <= 3);
    assert(fx.submap.num_nodes() - nodes_before == added);

    for (std::size_t ci = chords_before; ci < fx.submap.num_chords(); ++ci) {
        const Chord& c = fx.submap.chord(ci);
        assert(!c.dead && !c.is_null_length);
        assert_chord_geometrically_valid(fx.curve, c);

        assert(c.y_tag < fx.curve.num_vertices());
    }

    for (std::size_t ci = 0; ci < chords_before; ++ci)
        assert(!fx.submap.chord(ci).dead);
}

static void test_restore_conformality_e2e() {
    run_restore_and_check(1);
    std::printf("  [PASS] restore_conformality_e2e (piece_len=1)\n");
    run_restore_and_check(64);
    std::printf("  [PASS] restore_conformality_e2e (piece_len=64)\n");
}

static void test_insert_chord_asserts() {
    require_assertion_abort([] {
        CombFixture fx;
        std::size_t flat[] = {1, 3, 5, 9, 11, 13, 15};

        Submap::ChordPointSpec p{15, 0, RIGHT, Exact{3} / 5};
        Submap::ChordPointSpec q{11, 3, RIGHT, 7.0};
        fx.submap.insert_chord(p, q, SymbolicY{6.0, 2}, 0, flat, 7, fx.curve);
    });

    require_assertion_abort([] {
        CombFixture fx;
        std::size_t flat[] = {1, 3, 5, 9, 11, 13, 15};
        Submap::ChordPointSpec p{11, 4, RIGHT, 8.2};
        Submap::ChordPointSpec q{11, 3, RIGHT, 7.0};
        fx.submap.insert_chord(p, q, SymbolicY{4.7, 4}, 0, flat, 7, fx.curve);
    });

    std::printf("  [PASS] insert_chord_asserts\n");
}

int main() {
    std::setbuf(stdout, nullptr);
    std::printf("[C91 §3.2 tests]:\n");
    test_fixture_and_arc_end_y();
    test_cutter_wrapped_partial_endpoints();
    test_arc_sources();
    test_fused_region_cycle();
    test_local_shoot_fused();
    test_insert_chord();
    test_insert_chord_equal_x_wrap();
    test_find_visible_point();
    test_descent();
    test_restore_conformality_e2e();
    test_insert_chord_asserts();
    std::printf("All §3.2 tests passed.\n");
    return 0;
}
