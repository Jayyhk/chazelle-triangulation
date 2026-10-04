#include "merge/merge.h"
#include "merge/ray_shooting.h"
#include "support/arc_ray_shooter.h"
#include "support/assertions.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

using chazelle::test::require_assertion_abort;

struct StubRayShooter : RayShootingOracle {
    RayHit shoot(Point, Side, std::size_t, const Subarc&,
                 SourceOffset = SOURCE_OFFSET_NONE) const override {
        return {};
    }
};

struct StubArcCutter : ArcCuttingOracle {
    std::vector<ArcPiece> cut(std::size_t, const Subarc&) const override {
        return {};
    }
};

static const StubRayShooter STUB_RAY;
static const StubArcCutter STUB_ARC;

struct OracleRig {
    TestArcRayShooter first_ray_shooter, second_ray_shooter;
    OracleRig(const Submap& first_submap, const Polygon& first_curve, std::size_t g1,
              const Submap& second_submap, const Polygon& second_curve, std::size_t g2)
        : first_ray_shooter(first_submap, first_curve, g1),
          second_ray_shooter(second_submap, second_curve, g2) {}
};

static MergeInput make_input(const Polygon& first_curve, const Polygon& second_curve,
                             Submap& first_submap, Submap& second_submap, std::size_t g1,
                             std::size_t g2, std::size_t g) {
    MergeInput in;
    in.first_curve = &first_curve;
    in.second_curve = &second_curve;
    in.first_submap = &first_submap;
    in.second_submap = &second_submap;
    in.first_granularity = g1;
    in.second_granularity = g2;
    in.granularity = g;
    in.first_ray_shooter = &STUB_RAY;
    in.second_ray_shooter = &STUB_RAY;
    in.first_arc_cutter = &STUB_ARC;
    in.second_arc_cutter = &STUB_ARC;

    in.first_piece_count_bound = 1;
    in.second_piece_count_bound = 1;
    in.first_piece_granularity_bound = 1;
    in.second_piece_granularity_bound = 1;
    return in;
}

static const Polygon& input_P() {
    static Polygon input_curve({{0, 0, 0}, {3, 3, 1}, {5, 5, 2}, {7, 2, 3}, {10, 8, 4}});
    return input_curve;
}
static Polygon make_C1() {
    return input_P().subchain(0, 3);
}
static Polygon make_C2() {
    return input_P().subchain(2, 3);
}

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
    a.edge_count = 2 * poly.count_nonnull_edges(0, poly.num_edges() - 1);
    std::size_t ai0 = s.add_arc(a);
    assert(s.start_arc == ai0 && s.end_arc == ai0);
    s.build_tree_decomposition();
    return s;
}

static void test_valid_preconditions() {
    auto first_curve = make_C1();
    auto second_curve = make_C2();
    auto first_submap = make_single_region_submap(first_curve);
    auto second_submap = make_single_region_submap(second_curve);

    auto in = make_input(first_curve, second_curve, first_submap, second_submap, 4, 4, 4);
    assert_merge_preconditions(in);

    std::printf("  [PASS] valid_preconditions\n");
}

static void test_merged_curve() {
    auto first_curve = make_C1();
    auto second_curve = make_C2();
    auto first_submap = make_single_region_submap(first_curve);
    auto second_submap = make_single_region_submap(second_curve);

    OracleRig rig(first_submap, first_curve, 4, second_submap, second_curve, 4);
    auto in = make_input(first_curve, second_curve, first_submap, second_submap, 4, 4, 4);
    in.first_ray_shooter = &rig.first_ray_shooter;
    in.second_ray_shooter = &rig.second_ray_shooter;
    auto result = merge(in);

    assert(result.curve.num_vertices() == 5);
    assert(result.curve.num_edges() == 4);

    assert(result.curve.vertex(0).index == 0);
    assert(result.curve.vertex(1).index == 1);
    assert(result.curve.vertex(2).index == 2);
    assert(result.curve.vertex(3).index == 3);
    assert(result.curve.vertex(4).index == 4);

    assert(result.curve.vertex(0).x == 0.0 && result.curve.vertex(0).y == 0.0);
    assert(result.curve.vertex(2).x == 5.0 && result.curve.vertex(2).y == 5.0);
    assert(result.curve.vertex(4).x == 10.0 && result.curve.vertex(4).y == 8.0);

    std::printf("  [PASS] merged_curve\n");
}

static void test_gamma_ordering() {
    auto first_curve = make_C1();
    auto second_curve = make_C2();
    auto first_submap = make_single_region_submap(first_curve);
    auto second_submap = make_single_region_submap(second_curve);

    assert_merge_preconditions(
        make_input(first_curve, second_curve, first_submap, second_submap, 5, 5, 5));

    assert_merge_preconditions(
        make_input(first_curve, second_curve, first_submap, second_submap, 4, 5, 10));

    assert_merge_preconditions(
        make_input(first_curve, second_curve, first_submap, second_submap, 4, 5, 5));

    std::printf("  [PASS] gamma_ordering\n");
}

static void test_shared_vertex() {
    Polygon input_curve({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}});
    Polygon first_curve = input_curve.subchain(0, 2);
    Polygon second_curve = input_curve.subchain(1, 2);
    auto first_submap = make_single_region_submap(first_curve);
    auto second_submap = make_single_region_submap(second_curve);

    OracleRig rig(first_submap, first_curve, 2, second_submap, second_curve, 2);
    auto in = make_input(first_curve, second_curve, first_submap, second_submap, 2, 2, 2);
    in.first_ray_shooter = &rig.first_ray_shooter;
    in.second_ray_shooter = &rig.second_ray_shooter;
    auto result = merge(in);

    assert(result.curve.num_vertices() == 3);
    assert(result.curve.num_edges() == 2);

    assert(result.curve.vertex(1).x == 1.0);
    assert(result.curve.vertex(1).y == 1.0);
    assert(result.curve.vertex(1).index == 1);

    std::printf("  [PASS] shared_vertex\n");
}

static void test_merged_edges() {
    Polygon input_curve({{0, 0, 0}, {2, 3, 1}, {4, 1, 2}, {5, 5, 3}, {7, 2, 4}, {9, 6, 5}});
    Polygon first_curve = input_curve.subchain(0, 4);
    Polygon second_curve = input_curve.subchain(3, 3);
    auto first_submap = make_single_region_submap(first_curve);
    auto second_submap = make_single_region_submap(second_curve);

    OracleRig rig(first_submap, first_curve, 6, second_submap, second_curve, 6);
    auto in = make_input(first_curve, second_curve, first_submap, second_submap, 6, 6, 6);
    in.first_ray_shooter = &rig.first_ray_shooter;
    in.second_ray_shooter = &rig.second_ray_shooter;
    auto result = merge(in);

    assert(result.curve.num_vertices() == 6);
    assert(result.curve.num_edges() == 5);

    assert(result.curve.count_nonnull_edges(0, 4) == 5);

    std::printf("  [PASS] merged_edges\n");
}

static void test_conformality_required() {
    auto first_curve = make_C1();
    auto second_curve = make_C2();
    auto first_submap = make_single_region_submap(first_curve);
    auto second_submap = make_single_region_submap(second_curve);

    assert(first_submap.is_conformal());
    assert(second_submap.is_conformal());

    assert_merge_preconditions(
        make_input(first_curve, second_curve, first_submap, second_submap, 4, 4, 4));

    std::printf("  [PASS] conformality_required\n");
}

static Polygon make_segment_curve() {
    return Polygon({{0, 0, 100}, {1, 1, 101}});
}
static Submap make_segment_submap(const Polygon& curve) {
    Submap s;
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 1;

    Arc a{};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2 * curve.count_nonnull_edges(0, 0);
    std::size_t ai0 = s.add_arc(a);
    assert(s.start_arc == ai0 && s.end_arc == ai0);
    s.build_tree_decomposition();
    return s;
}

static void test_cut_postconditions_valid() {
    Polygon curve = make_segment_curve();
    Subarc target = test_full_subarc(curve, 0, LEFT, 0);
    Submap submap = make_segment_submap(curve);

    ArcPiece p;
    p.subarc = target;
    p.submap = &submap;
    p.curve = &curve;
    p.is_boundary_piece = false;

    p.granularity = 2;
    assert_cut_postconditions(curve, target, &p, 1, 4, 2);

    std::printf("  [PASS] cut_postconditions_valid\n");
}

static void test_cut_count_zero_fires() {
    require_assertion_abort([] {
        Subarc t{0, LEFT, 0, LEFT};
        Polygon curve = make_segment_curve();
        ArcPiece p{};
        assert_cut_postconditions(curve, t, &p, 0, 4, 1);
    });
    std::printf("  [PASS] cut_count_zero_fires\n");
}

static void test_cut_count_exceeds_bound_fires() {
    require_assertion_abort([] {
        Subarc t{0, LEFT, 0, LEFT};
        Polygon curve = make_segment_curve();
        Submap submap = make_segment_submap(curve);
        ArcPiece pieces[5];
        for (auto& p : pieces) {
            p.subarc = {0, LEFT, 0, LEFT};
            p.submap = &submap;
            p.curve = &curve;
            p.is_boundary_piece = false;
        }
        assert_cut_postconditions(curve, t, pieces, 5, 4, 1);
    });
    std::printf("  [PASS] cut_count_exceeds_bound_fires\n");
}

static void test_cut_first_endpoint_mismatch_fires() {
    require_assertion_abort([] {
        Subarc target{0, LEFT, 1, LEFT};
        Polygon curve = make_segment_curve();
        Submap submap = make_segment_submap(curve);
        ArcPiece p;
        p.subarc = {1, LEFT, 1, LEFT};
        p.submap = &submap;
        p.curve = &curve;
        p.is_boundary_piece = false;
        assert_cut_postconditions(curve, target, &p, 1, 4, 1);
    });
    std::printf("  [PASS] cut_first_endpoint_mismatch_fires\n");
}

static void test_cut_last_endpoint_mismatch_fires() {
    require_assertion_abort([] {
        Subarc target{0, LEFT, 1, LEFT};
        Polygon curve = make_segment_curve();
        Submap submap = make_segment_submap(curve);
        ArcPiece p;
        p.subarc = {0, LEFT, 0, LEFT};
        p.submap = &submap;
        p.curve = &curve;
        p.is_boundary_piece = false;
        assert_cut_postconditions(curve, target, &p, 1, 4, 1);
    });
    std::printf("  [PASS] cut_last_endpoint_mismatch_fires\n");
}

static void test_cut_double_back_fires() {
    require_assertion_abort([] {
        Subarc target{0, LEFT, 0, RIGHT};
        Polygon curve = make_segment_curve();
        Submap submap = make_segment_submap(curve);
        ArcPiece p;
        p.subarc = {0, LEFT, 0, RIGHT};
        p.submap = &submap;
        p.curve = &curve;
        p.is_boundary_piece = false;
        assert_cut_postconditions(curve, target, &p, 1, 4, 1);
    });
    std::printf("  [PASS] cut_double_back_fires\n");
}

static void test_cut_non_clockwise_left_fires() {
    require_assertion_abort([] {
        Subarc target{1, LEFT, 0, LEFT};
        Polygon curve = make_segment_curve();
        Submap submap = make_segment_submap(curve);
        ArcPiece p;
        p.subarc = {1, LEFT, 0, LEFT};
        p.submap = &submap;
        p.curve = &curve;
        p.is_boundary_piece = false;
        assert_cut_postconditions(curve, target, &p, 1, 4, 1);
    });
    std::printf("  [PASS] cut_non_clockwise_left_fires\n");
}

static void test_cut_interior_boundary_piece_fires() {
    require_assertion_abort([] {
        Subarc target{0, LEFT, 0, LEFT};
        Polygon curve = make_segment_curve();
        Submap submap = make_segment_submap(curve);
        ArcPiece pieces[3];
        for (auto& p : pieces) {
            p.subarc = {0, LEFT, 0, LEFT};
            p.submap = &submap;
            p.curve = &curve;
            p.is_boundary_piece = false;
        }
        pieces[1].is_boundary_piece = true;
        pieces[1].submap = nullptr;
        pieces[1].curve = nullptr;
        assert_cut_postconditions(curve, target, pieces, 3, 4, 1);
    });
    std::printf("  [PASS] cut_interior_boundary_piece_fires\n");
}

static void test_cut_non_boundary_null_submap_fires() {
    require_assertion_abort([] {
        Subarc target{0, LEFT, 0, LEFT};
        Polygon curve = make_segment_curve();
        ArcPiece p;
        p.subarc = {0, LEFT, 0, LEFT};
        p.submap = nullptr;
        p.curve = nullptr;
        p.is_boundary_piece = false;
        assert_cut_postconditions(curve, target, &p, 1, 4, 1);
    });
    std::printf("  [PASS] cut_non_boundary_null_submap_fires\n");
}

static void test_cut_boundary_multi_edge_fires() {
    require_assertion_abort([] {
        Subarc target{0, LEFT, 1, LEFT};
        Polygon curve = make_segment_curve();
        ArcPiece pieces[2];
        pieces[0].subarc = {0, LEFT, 1, LEFT};
        pieces[0].is_boundary_piece = true;
        pieces[0].submap = nullptr;
        pieces[0].curve = nullptr;
        pieces[1].subarc = {1, LEFT, 1, LEFT};
        pieces[1].is_boundary_piece = true;
        pieces[1].submap = nullptr;
        pieces[1].curve = nullptr;
        assert_cut_postconditions(curve, target, pieces, 2, 4, 1);
    });
    std::printf("  [PASS] cut_boundary_multi_edge_fires\n");
}

static void test_cut_exact_endpoint_asserts() {
    require_assertion_abort([] {
        Polygon curve({{0, 0, 0}, {4, 4, 1}});
        Subarc target = test_full_subarc(curve, 0, LEFT, 0);
        ArcPiece p;
        p.is_boundary_piece = true;
        p.subarc = target;
        p.subarc.last_y = {3.0, 42};
        assert_cut_postconditions(curve, target, &p, 1, 2, 2);
    });
    require_assertion_abort([] {
        Polygon curve({{0, 0, 0}, {4, 4, 1}});
        Subarc target = test_full_subarc(curve, 0, LEFT, 0);
        ArcPiece p[2];
        for (auto& x : p) {
            x.is_boundary_piece = true;
            x.subarc = target;
        }
        p[0].subarc.last_y = {1.0, 42};
        p[1].subarc.first_y = {2.0, 43};
        assert_cut_postconditions(curve, target, p, 2, 2, 2);
    });
    std::printf("  [PASS] cut_exact_endpoint_asserts\n");
}

static void test_merge_preconds_null_curve_fires() {
    require_assertion_abort([] {
        auto second_curve = make_C2();
        auto second_submap = make_single_region_submap(second_curve);
        Polygon C1stub({{0, 0, 0}, {1, 1, 1}});
        auto first_submap = make_single_region_submap(C1stub);
        MergeInput in;
        in.first_curve = nullptr;
        in.second_curve = &second_curve;
        in.first_submap = &first_submap;
        in.second_submap = &second_submap;
        in.first_granularity = 1;
        in.second_granularity = 1;
        in.granularity = 1;
        in.first_ray_shooter = &STUB_RAY;
        in.second_ray_shooter = &STUB_RAY;
        in.first_arc_cutter = &STUB_ARC;
        in.second_arc_cutter = &STUB_ARC;
        assert_merge_preconditions(in);
    });
    std::printf("  [PASS] merge_preconds_null_curve_fires\n");
}

static void test_merge_preconds_unshared_junction_fires() {
    require_assertion_abort([] {
        Polygon first_curve({{0, 0, 0}, {1, 1, 1}});
        Polygon second_curve({{2, 2, 2}, {3, 3, 3}});
        auto first_submap = make_single_region_submap(first_curve);
        auto second_submap = make_single_region_submap(second_curve);
        assert_merge_preconditions(
            make_input(first_curve, second_curve, first_submap, second_submap, 1, 1, 1));
    });
    std::printf("  [PASS] merge_preconds_unshared_junction_fires\n");
}

static void test_merge_preconds_gamma_order_fires() {
    require_assertion_abort([] {
        auto first_curve = make_C1();
        auto second_curve = make_C2();
        auto first_submap = make_single_region_submap(first_curve);
        auto second_submap = make_single_region_submap(second_curve);
        assert_merge_preconditions(
            make_input(first_curve, second_curve, first_submap, second_submap, 5, 3, 5));
    });
    std::printf("  [PASS] merge_preconds_gamma_order_fires\n");
}

static void test_merge_preconds_gamma_target_low_fires() {
    require_assertion_abort([] {
        auto first_curve = make_C1();
        auto second_curve = make_C2();
        auto first_submap = make_single_region_submap(first_curve);
        auto second_submap = make_single_region_submap(second_curve);
        assert_merge_preconditions(
            make_input(first_curve, second_curve, first_submap, second_submap, 2, 5, 3));
    });
    std::printf("  [PASS] merge_preconds_gamma_target_low_fires\n");
}

static void test_merge_preconds_null_oracle_fires() {
    require_assertion_abort([] {
        auto first_curve = make_C1();
        auto second_curve = make_C2();
        auto first_submap = make_single_region_submap(first_curve);
        auto second_submap = make_single_region_submap(second_curve);
        auto in = make_input(first_curve, second_curve, first_submap, second_submap, 2, 2, 2);
        in.first_ray_shooter = nullptr;
        assert_merge_preconditions(in);
    });
    std::printf("  [PASS] merge_preconds_null_oracle_fires\n");
}

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
    test_cut_exact_endpoint_asserts();
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
