#include "algorithm/merge/conformality.h"
#include "algorithm/merge/granularity.h"
#include "algorithm/polygon/polygon.h"
#include "algorithm/submap/submap.h"
#include "support/arc_ray_shooter.h"
#include "support/assertions.h"

#include <cassert>
#include <cstdio>
#include <memory>
#include <vector>

using namespace chazelle;

using chazelle::test::require_assertion_abort;

namespace {

std::size_t brute_region_weight(const Submap& s, std::size_t region) {
    std::size_t w = 0;
    for (std::size_t ai = 0; ai < s.num_arcs(); ++ai) {
        const Arc& a = s.arc(ai);
        if (a.dead || a.region_node != region)
            continue;
        if (a.edge_count > w)
            w = a.edge_count;
    }
    return w;
}
}

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
            submaps.push_back(std::make_unique<Submap>(make_chordless_normal(*curves.back())));
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

static Polygon chain_polygon() {
    return Polygon({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4}});
}

static Submap build_3region_submap() {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a{};
    auto add = [&](std::size_t fe, Side fs, std::size_t le, Side ls, std::size_t region,
                   std::size_t count) {
        a = {};
        a.first_edge = fe;
        a.last_edge = le;
        a.first_side = fs;
        a.last_side = ls;
        a.region_node = region;
        a.edge_count = count;
        return s.add_arc(a);
    };
    std::size_t a1 = add(1, LEFT, 1, LEFT, r1, 1);
    std::size_t aE = add(2, LEFT, 2, RIGHT, r2, 4);
    std::size_t a4 = add(1, RIGHT, 1, RIGHT, r1, 1);
    std::size_t aS = add(0, RIGHT, 0, LEFT, r0, 2);

    Chord c0{};
    c0.region[0] = r0;
    c0.region[1] = r1;
    c0.left_adj = {{aS}, 1};
    c0.right_adj = {{a4}, 1};
    c0.left_edge = 0;
    c0.right_edge = 0;
    c0.y = 1.0;
    c0.y_tag = 1;
    s.add_chord(c0);

    Chord c1{};
    c1.region[0] = r1;
    c1.region[1] = r2;
    c1.left_adj = {{a1}, 1};
    c1.right_adj = {{aE}, 1};
    c1.left_edge = 1;
    c1.right_edge = 1;
    c1.y = 2.0;
    c1.y_tag = 2;
    s.add_chord(c1);

    assert(s.start_arc == aS && s.end_arc == aE &&
           "[C91 §2.4(iii) tex 138]: endpoint arcs auto-registered");
    return s;
}

static void test_enforce_3region_partial() {
    Submap s = build_3region_submap();
    Polygon poly = chain_polygon();
    s.check_invariants(poly);

    enforce_granularity(s, poly, 4);

    assert(s.num_live_nodes() == 2 && s.num_live_chords() == 1);
    assert(!s.chord(1).dead && "c1's contraction weight 8 > γ = 4: kept");

    assert(s.num_live_arcs() == 2);
    std::size_t four_edge = 0;
    for (std::size_t ai = 0; ai < s.num_arcs(); ++ai)
        if (!s.arc(ai).dead && s.arc(ai).edge_count == 4)
            ++four_edge;
    assert(four_edge == 2 && "S+a1+a4 glued into one wrap arc; E untouched");

    s.normalize(poly);
    s.check_invariants(poly);
    assert(!s.tree_decomposition().empty());
    assert(s.is_conformal());
    assert(s.is_granular(4, poly) && "[C91 Lemma 3.5]: output must be γ-granular");

    for (std::size_t r = 0; r < s.num_nodes(); ++r) {
        if (s.node(r).dead)
            continue;
        assert(s.region_weight(r) == brute_region_weight(s, r));
    }

    auto res = s.double_identify(1, SymbolicY{2.0, 2}, poly);
    assert(res.count >= 1 && res.count <= 6);

    std::printf("  [PASS] enforce_3region_partial\n");
}

static void test_enforce_3region_full() {
    Submap s = build_3region_submap();
    Polygon poly = chain_polygon();

    enforce_granularity(s, poly, 8);

    assert(s.num_live_nodes() == 1 && s.num_live_chords() == 0);

    assert(s.num_live_arcs() == 1 && "full contraction leaves the single closed arc");
    for (std::size_t ai = 0; ai < s.num_arcs(); ++ai)
        if (!s.arc(ai).dead) {
            assert(s.arc(ai).first_side == LEFT && s.arc(ai).last_side == RIGHT &&
                   s.arc(ai).first_edge == 0 && s.arc(ai).last_edge == 0);
            assert(s.arc(ai).edge_count == 8);
        }

    s.normalize(poly);
    s.check_invariants(poly);

    assert(s.is_granular(8, poly));
    assert(s.region_weight(0) == 8 && s.region_weight(0) == brute_region_weight(s, 0));

    std::printf("  [PASS] enforce_3region_full\n");
}

static void test_comb_glue_chains() {
    CombFixture fx;
    fx.submap.check_invariants(fx.curve);

    fx.submap.remove_chord(3, fx.curve);
    fx.submap.assert_tree_property();
    assert(fx.submap.node(4).dead && "pocket P4 merged into the outer region");
    assert(!fx.submap.arc(8).dead && fx.submap.arc(8).edge_count == 1 &&
           "P4b + zero-length Z3: ec unchanged");
    assert(fx.submap.arc(9).dead && "Z3 marked as removed by the vertex glue");
    assert(!fx.submap.arc(5).dead && fx.submap.arc(5).edge_count == 2 &&
           "R2 + P4a glued at the mid-edge endpoint (ec 2)");
    assert(fx.submap.arc(6).dead && "P4a marked as removed by the mid-edge glue");

    fx.submap.remove_chord(7, fx.curve);
    fx.submap.assert_tree_property();
    assert(fx.submap.node(8).dead && "n7's empty inner region merged");
    assert(!fx.submap.arc(5).dead && fx.submap.arc(5).edge_count == 3 &&
           "[C91 §2.2 tex 108]: null chord removal fuses "
           "before + null + after into one arc");
    assert(fx.submap.arc(7).dead && fx.submap.arc(8).dead);

    assert(fx.submap.region_weight(0) == brute_region_weight(fx.submap, 0));

    std::printf("  [PASS] comb_glue_chains\n");
}

static void test_precondition_asserts() {
    require_assertion_abort([] {
        CombFixture fx;
        enforce_granularity(fx.submap, fx.curve, 100);
    });

    require_assertion_abort([] {
        Submap s = build_3region_submap();
        Polygon poly = chain_polygon();
        enforce_granularity(s, poly, 1);
    });

    require_assertion_abort([] {
        CombFixture fx;
        fx.submap.normalize(fx.curve);
    });

    std::printf("  [PASS] precondition_asserts\n");
}

static void run_pipeline_and_check(std::size_t granularity, bool expect_full_contraction) {
    CombRig rig;
    auto& fx = rig.fx;

    restore_conformality(fx.submap, fx.curve, rig.oracles);
    assert(fx.submap.is_conformal());

    enforce_granularity(fx.submap, fx.curve, granularity);
    fx.submap.normalize(fx.curve);

    fx.submap.check_invariants(fx.curve);
    assert(fx.submap.is_conformal());
    assert(fx.submap.is_granular(granularity, fx.curve));
    assert(!fx.submap.tree_decomposition().empty());

    if (expect_full_contraction) {
        assert(fx.submap.num_live_chords() == 0 && fx.submap.num_live_nodes() == 1 &&
               "γ ≥ total weight forces full contraction");
    }

    for (std::size_t r = 0; r < fx.submap.num_nodes(); ++r) {
        if (fx.submap.node(r).dead)
            continue;
        assert(fx.submap.region_weight(r) == brute_region_weight(fx.submap, r));
        assert(fx.submap.region_weight(r) <= granularity &&
               "[C91 §2.3 tex 120]: criterion (i) after §3.3");
    }

    auto res = fx.submap.double_identify(3, SymbolicY{6.0, 2}, fx.curve);
    assert(res.count >= 1 && res.count <= 6);
    for (std::size_t ai : res) {
        const Arc& a = fx.submap.arc(ai);
        auto [lo, hi] = a.underlying_edge_range(0, 12);
        assert(lo <= 3 && 3 <= hi && "double_identify result must contain the query edge");
    }
}

static void test_pipeline_gammas() {
    run_pipeline_and_check(8, false);
    std::printf("  [PASS] pipeline (gamma=8)\n");
    run_pipeline_and_check(10, false);
    std::printf("  [PASS] pipeline (gamma=10)\n");
    run_pipeline_and_check(100, true);
    std::printf("  [PASS] pipeline (gamma=100)\n");
}

static void test_normalize_repairs_table() {
    CombRig rig;
    auto& fx = rig.fx;
    std::size_t arcs_before = fx.submap.num_arcs();

    restore_conformality(fx.submap, fx.curve, rig.oracles);
    assert(fx.submap.num_arcs() > arcs_before && "§3.2 appends split arc halves");

    {
        require_assertion_abort([] {
            CombRig r2;
            restore_conformality(r2.fx.submap, r2.fx.curve, r2.oracles);
            (void)r2.fx.submap.double_identify(3, SymbolicY{6.0, 2}, r2.fx.curve);
        });
    }

    fx.submap.normalize(fx.curve);
    fx.submap.check_invariants(fx.curve);
    auto res = fx.submap.double_identify(3, SymbolicY{6.0, 2}, fx.curve);
    assert(res.count >= 1);

    fx.submap.normalize(fx.curve);
    fx.submap.check_invariants(fx.curve);

    std::printf("  [PASS] normalize_repairs_table\n");
}

static void test_degree_drop_recheck() {
    CombRig rig;
    auto& fx = rig.fx;
    restore_conformality(fx.submap, fx.curve, rig.oracles);
    fx.submap.normalize(fx.curve);

    std::size_t nch = fx.submap.num_chords();
    assert(nch >= 10 && "comb + §3.2 chords: 8 fused + ≥2 added");

    Submap clone;
    for (std::size_t r = 0; r < fx.submap.num_nodes(); ++r) {
        clone.add_node();
        clone.node(r).dead = fx.submap.node(r).dead;
    }
    clone.start_vertex = fx.submap.start_vertex;
    clone.end_vertex = fx.submap.end_vertex;
    for (std::size_t ai = 0; ai < fx.submap.num_arcs(); ++ai)
        clone.add_arc(fx.submap.arc(ai));
    for (std::size_t ci = nch; ci-- > 0;) {
        Chord c = fx.submap.chord(ci);
        clone.add_chord(c);
    }
    clone.start_arc = fx.submap.start_arc;
    clone.end_arc = fx.submap.end_arc;
    clone.start_vertex = fx.submap.start_vertex;
    clone.end_vertex = fx.submap.end_vertex;
    clone.check_invariants(fx.curve);

    enforce_granularity(clone, fx.curve, 100);
    assert(clone.num_live_chords() == 0 && clone.num_live_nodes() == 1 &&
           "[C91 Lemma 3.5]: γ ≥ total weight must fully contract");

    clone.normalize(fx.curve);
    assert(clone.is_granular(100, fx.curve));

    std::printf("  [PASS] degree_drop_recheck\n");
}

int main() {
    std::setbuf(stdout, nullptr);
    std::printf("[C91 §3.3 tests]:\n");
    test_enforce_3region_partial();
    test_enforce_3region_full();
    test_comb_glue_chains();
    test_precondition_asserts();
    test_pipeline_gammas();
    test_normalize_repairs_table();
    test_degree_drop_recheck();
    std::printf("All §3.3 tests passed.\n");
    return 0;
}
