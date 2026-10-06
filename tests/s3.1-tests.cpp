#include "algorithm/merge/fusion.h"
#include "algorithm/polygon/polygon.h"
#include "algorithm/submap/submap.h"
#include "support/arc_ray_shooter.h"
#include "support/assertions.h"

#include <cassert>
#include <cstdio>
#include <utility>

using namespace chazelle;

using chazelle::test::require_assertion_abort;

static Polygon make_C1() {
    return Polygon({{0, 0, 0}, {1, 2, 1}, {2, 4, 2}, {3, 1, 3}, {4, 3, 4}});
}

static Submap make_S1(const Polygon& curve) {
    Submap s;
    s.add_node();
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a{};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 1;
    a.edge_count = 3;
    std::size_t aiE = s.add_arc(a);

    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = 0;
    a.edge_count = 2;
    std::size_t aiS = s.add_arc(a);

    Chord c{};
    c.region[0] = 0;
    c.region[1] = 1;
    c.left_edge = 1;
    c.right_edge = 1;
    c.left_side = LEFT;
    c.right_side = RIGHT;
    c.y = curve.vertex(2).y;
    c.y_tag = 2;
    c.left_adj = {{aiS}, 1};
    c.right_adj = {{aiE}, 1};
    s.add_chord(c);

    assert(s.start_arc == aiS && s.end_arc == aiE &&
           "[C91 §2.4(iii) tex 138]: endpoint arcs auto-registered");
    return s;
}

static Submap make_S1_chain(const Polygon& curve) {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a{};
    a.first_edge = 1;
    a.last_edge = 2;
    a.first_side = LEFT;
    a.last_side = LEFT;
    a.region_node = r1;
    a.edge_count = 2;
    std::size_t a1 = s.add_arc(a);
    a = {};
    a.first_edge = 3;
    a.last_edge = 3;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = r2;
    a.edge_count = 1;
    std::size_t aE = s.add_arc(a);
    a = {};
    a.first_edge = 2;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r1;
    a.edge_count = 2;
    std::size_t a4 = s.add_arc(a);
    a = {};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = r0;
    a.edge_count = 1;
    std::size_t aS = s.add_arc(a);

    Chord c{};
    c.region[0] = r0;
    c.region[1] = r1;
    c.left_edge = 0;
    c.right_edge = 0;
    c.left_side = LEFT;
    c.right_side = RIGHT;
    c.y = curve.vertex(1).y;
    c.y_tag = 1;
    c.left_adj = {{aS}, 1};
    c.right_adj = {{a4}, 1};
    s.add_chord(c);
    c = {};
    c.region[0] = r1;
    c.region[1] = r2;
    c.left_edge = 2;
    c.right_edge = 2;
    c.left_side = LEFT;
    c.right_side = RIGHT;
    c.y = curve.vertex(3).y;
    c.y_tag = 3;
    c.left_adj = {{a1}, 1};
    c.right_adj = {{aE}, 1};
    s.add_chord(c);

    assert(s.start_arc == aS && s.end_arc == aE);
    return s;
}

static void test_fusion_sequence_basic() {
    auto first_curve = make_C1();
    auto first_submap = make_S1(first_curve);

    FusionState state;
    build_fusion_sequence(state, first_submap, first_curve);
    const auto& seq = state.sequence;

    assert(seq.size() == 4);

    assert(seq.front().is_companion);
    assert(seq.front().side == RIGHT);
    assert(seq.front().edge == 3);

    assert(seq.back().is_companion);
    assert(seq.back().side == LEFT);
    assert(seq.back().edge == 3);

    assert(symbolic_y_equal(seq.front().y, symbolic_y_of(first_curve.vertex(4))));
    assert(symbolic_y_equal(seq.back().y, symbolic_y_of(first_curve.vertex(4))));

    assert(!seq[1].is_companion);
    assert(!seq[2].is_companion);
    assert(seq[1].chord_idx == 0);
    assert(seq[2].chord_idx == 0);

    std::printf("  [PASS] fusion_sequence_basic\n");
}

static void test_fusion_sequence_no_chords() {
    Polygon curve({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}});

    Submap s;
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 2;
    Arc a{};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2;
    s.add_arc(a);

    FusionState state;
    build_fusion_sequence(state, s, curve);
    const auto& seq = state.sequence;

    assert(seq.size() == 2);
    assert(seq[0].is_companion && seq[0].side == RIGHT);
    assert(seq[1].is_companion && seq[1].side == LEFT);

    std::printf("  [PASS] fusion_sequence_no_chords\n");
}

static void test_fusion_sequence_ordering() {
    auto first_curve = make_C1();
    auto first_submap = make_S1(first_curve);

    FusionState state;
    build_fusion_sequence(state, first_submap, first_curve);
    const auto& seq = state.sequence;

    assert(seq[1].side == RIGHT);
    assert(seq[2].side == LEFT);

    std::printf("  [PASS] fusion_sequence_ordering\n");
}

static void test_companion_identity() {
    auto first_curve = make_C1();
    auto first_submap = make_S1(first_curve);

    FusionState state;
    build_fusion_sequence(state, first_submap, first_curve);
    const auto& seq = state.sequence;

    SymbolicY junction_y = symbolic_y_of(first_curve.vertex(4));
    assert(symbolic_y_equal(seq.front().y, junction_y));
    assert(symbolic_y_equal(seq.back().y, junction_y));

    assert(seq.front().edge == first_curve.num_edges() - 1);
    assert(seq.back().edge == first_curve.num_edges() - 1);

    std::printf("  [PASS] companion_identity\n");
}

static void test_collect_region_arcs() {
    auto first_curve = make_C1();
    auto first_submap = make_S1(first_curve);

    auto arcs0 = collect_region_arcs(first_submap, 0);
    assert(arcs0.count == 1);

    auto arcs1 = collect_region_arcs(first_submap, 1);
    assert(arcs1.count == 1);

    std::printf("  [PASS] collect_region_arcs\n");
}

struct FixedRayShooter : RayShootingOracle {
    RayHit shoot(Point p, Side, std::size_t, const Subarc& target,
                 SourceOffset = SOURCE_OFFSET_NONE) const override {
        RayHit h;
        h.hit = true;
        h.x = 5.0 + (Exact)target.first_edge + (target.first_side == RIGHT ? 0.25 : 0.0);
        h.y = p.y;
        h.edge = target.first_edge;
        h.side = target.first_side;
        return h;
    }
};

static void test_local_shoot() {
    auto first_curve = make_C1();
    auto first_submap = make_S1_chain(first_curve);

    FixedRayShooter oracle;
    Point p{0.0, 2.0, 99};

    auto hit = local_shoot(p, RIGHT, 0, first_submap, first_curve, oracle, false);
    assert(hit.hit);

    hit = local_shoot(p, RIGHT, 1, first_submap, first_curve, oracle, false);
    assert(hit.hit);

    std::printf("  [PASS] local_shoot\n");
}

struct DistanceRayShooter : RayShootingOracle {
    RayHit shoot(Point p, Side, std::size_t, const Subarc& target,
                 SourceOffset = SOURCE_OFFSET_NONE) const override {
        RayHit h;
        h.hit = true;
        h.y = p.y;
        h.edge = target.first_edge;
        h.side = target.first_side;

        h.x = (target.first_side == LEFT) ? 10.0 : 3.0;
        return h;
    }
};

static void test_local_shoot_nearest() {
    auto first_curve = make_C1();
    auto first_submap = make_S1_chain(first_curve);

    DistanceRayShooter oracle;
    Point p{1.0, 2.0, 99};

    auto hit = local_shoot(p, RIGHT, 1, first_submap, first_curve, oracle, false);
    assert(hit.hit);
    assert(hit.x == 3.0);

    Point p2{15.0, 2.0, 99};
    hit = local_shoot(p2, LEFT, 1, first_submap, first_curve, oracle, false);
    assert(hit.hit);
    assert(hit.x == 10.0);

    std::printf("  [PASS] local_shoot_nearest\n");
}

struct StartupOracle : RayShootingOracle {
    const Polygon* curve;
    Exact hit_x;
    std::size_t prefer_edge;
    StartupOracle(const Polygon* c, Exact x, std::size_t prefer = NONE)
        : curve(c), hit_x(std::move(x)), prefer_edge(prefer) {}
    RayHit shoot(Point p, Side dir, std::size_t, const Subarc& target,
                 SourceOffset = SOURCE_OFFSET_NONE) const override {
        SymbolicY sy{p.y, p.index};
        ArcSideRange ranges[3];
        std::size_t nl = subarc_side_ranges(target, 0, curve->num_vertices() - 1, ranges);
        for (int pass = 0; pass < 2; ++pass)
            for (std::size_t g = 0; g < nl; ++g) {
                for (std::size_t e = ranges[g].first_edge; e <= ranges[g].last_edge; ++e) {
                    if (pass == 0 && prefer_edge != NONE && e != prefer_edge)
                        continue;
                    if (!subarc_contains_point(target, *curve, e, ranges[g].side, sy, 0,
                                               curve->num_vertices() - 1))
                        continue;
                    RayHit h;
                    h.hit = true;
                    h.x = hit_x;
                    h.y = p.y;
                    h.edge = e;
                    h.side = ranges[g].side;
                    Exact d = (dir == RIGHT) ? (h.x - p.x) : (p.x - h.x);
                    h.wrapped = (d <= 0.0);
                    return h;
                }
            }
        return {};
    }
};

static void test_startup_case1() {
    Polygon input_curve(
        {{0, 0, 0}, {1, 2, 1}, {2, 4, 2}, {3, 1, 3}, {4, 3, 4}, {5, 5, 5}, {6, 1, 6}});
    Polygon first_curve = input_curve.subchain(0, 5);
    Polygon second_curve = input_curve.subchain(4, 3);
    auto first_submap = make_S1(first_curve);

    Submap second_submap;
    second_submap.add_node();
    second_submap.start_vertex = 0;
    second_submap.end_vertex = 2;
    Arc a{};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2;
    std::size_t ai0 = second_submap.add_arc(a);
    assert(second_submap.start_arc == ai0 && second_submap.end_arc == ai0);

    StartupOracle oracle1(&first_curve, 18.0);
    StartupOracle oracle2(&second_curve, 5.0);

    FusionState state;
    build_fusion_sequence(state, first_submap, first_curve);

    std::size_t start = fusion_startup(state, first_submap, first_curve, second_submap,
                                       second_curve, oracle1, oracle2);

    assert(start == 1);
    assert(state.current_region != NONE);
    assert(state.current_region == 0);
    assert(!state.chords.empty());

    std::printf("  [PASS] startup_case1\n");
}

static void test_startup_case2() {
    Polygon input_curve(
        {{0, 0, 0}, {1, 2, 1}, {2, 4, 2}, {3, 1, 3}, {4, 3, 4}, {5, 5, 5}, {6, 1, 6}});
    Polygon first_curve = input_curve.subchain(0, 5);
    Polygon second_curve = input_curve.subchain(4, 3);
    auto first_submap = make_S1(first_curve);

    Submap second_submap;
    second_submap.add_node();
    second_submap.start_vertex = 0;
    second_submap.end_vertex = 2;
    Arc a{};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2;
    std::size_t ai0 = second_submap.add_arc(a);
    assert(second_submap.start_arc == ai0 && second_submap.end_arc == ai0);

    StartupOracle oracle1(&first_curve, 5.0);
    StartupOracle oracle2(&second_curve, 18.0);

    FusionState state;
    build_fusion_sequence(state, first_submap, first_curve);

    fusion_startup(state, first_submap, first_curve, second_submap, second_curve, oracle1, oracle2);

    assert(state.current_region != NONE);
    assert(state.current_region == 0);
    assert(!state.chords.empty());

    std::printf("  [PASS] startup_case2\n");
}

static void test_shooting_direction_all_cases() {
    Polygon up({{0, 0, 0}, {1, 2, 1}});
    Polygon down({{0, 2, 0}, {1, 0, 1}});

    assert(shooting_direction(0, LEFT, up) == LEFT);

    assert(shooting_direction(0, LEFT, down) == RIGHT);

    assert(shooting_direction(0, RIGHT, up) == RIGHT);

    assert(shooting_direction(0, RIGHT, down) == LEFT);

    std::printf("  [PASS] shooting_direction_all_cases\n");
}

static void test_ray_contact_tie_break() {
    {
        Polygon curve({{0, 0, 0}, {2, 4, 1}, {4, 0, 2}});
        SymbolicY sy = symbolic_y_of(curve.vertex(1));
        for (Side dir : {LEFT, RIGHT}) {
            auto opposing = [&](std::size_t e) -> Side {
                const auto& ed = curve.edge(e);
                bool asc = symbolic_y_less(symbolic_y_of(curve.vertex(ed.start_idx)),
                                           symbolic_y_of(curve.vertex(ed.end_idx)));
                Side minus_x = asc ? LEFT : RIGHT;
                return (dir == RIGHT) ? minus_x : (minus_x == LEFT ? RIGHT : LEFT);
            };
            Side s0 = opposing(0), s1 = opposing(1);
            bool i0 = is_inside_companion(curve, 0, s0, 1);
            bool i1 = is_inside_companion(curve, 1, s1, 1);
            assert(i0 != i1 && "[C91 §2.1 tex 72]: exactly one apex wall is the "
                               "inside-of-turn duplicate");

            bool new0_wins = ray_contact_precedes(curve, sy, dir, 0, s0, 1, s1);
            bool new1_wins = ray_contact_precedes(curve, sy, dir, 1, s1, 0, s0);
            assert(new0_wins == (i1 && !i0));
            assert(new1_wins == (i0 && !i1));
            assert(!(new0_wins && new1_wins) && "antisymmetry");
        }
    }

    {
        Polygon curve({{0, 0, 0}, {4, -3, 1}, {8, 2, 2}, {8, 2, 3}, {8, 2, 4}});
        assert(curve.edge_is_null(2) && curve.edge_is_null(3) && !curve.edge_is_null(1));
        SymbolicY sy = symbolic_y_of(curve.vertex(3));

        assert(symbolic_y_less(symbolic_y_of(curve.vertex(1)), sy) &&
               symbolic_y_less(sy, symbolic_y_of(curve.vertex(2))));

        for (std::size_t null_e : {std::size_t{2}, std::size_t{3}}) {
            for (Side ns : {LEFT, RIGHT}) {
                for (Side es : {LEFT, RIGHT}) {
                    assert(ray_contact_precedes(curve, sy, RIGHT, 1, es, null_e, ns));
                    assert(!ray_contact_precedes(curve, sy, RIGHT, null_e, ns, 1, es));

                    assert(ray_contact_precedes(curve, sy, LEFT, null_e, ns, 1, es));
                    assert(!ray_contact_precedes(curve, sy, LEFT, 1, es, null_e, ns));
                }
            }
        }
    }

    {
        Polygon curve({{0, 0, 0}, {8, -3, 1}, {8, 2, 2}, {8, 2, 3}, {8, 2, 4}});
        SymbolicY sy = symbolic_y_of(curve.vertex(3));
        assert(symbolic_y_less(symbolic_y_of(curve.vertex(1)), sy) &&
               symbolic_y_less(sy, symbolic_y_of(curve.vertex(2))));
        for (Side dir : {LEFT, RIGHT}) {
            assert(ray_contact_precedes(curve, sy, dir, 1, LEFT, 2, LEFT));
            assert(!ray_contact_precedes(curve, sy, dir, 2, LEFT, 1, LEFT));
        }
    }

    std::printf("  [PASS] ray_contact_tie_break\n");
}

static void test_startup_d1_eq_d2_defaults_to_case1() {
    Polygon input_curve(
        {{0, 0, 0}, {1, 2, 1}, {2, 4, 2}, {3, 1, 3}, {4, 3, 4}, {5, 5, 5}, {6, 1, 6}});
    Polygon first_curve = input_curve.subchain(0, 5);
    Polygon second_curve = input_curve.subchain(4, 3);
    auto first_submap = make_S1(first_curve);

    Submap second_submap;
    second_submap.add_node();
    second_submap.start_vertex = 0;
    second_submap.end_vertex = 2;
    Arc a{};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2;
    std::size_t ai0 = second_submap.add_arc(a);
    assert(second_submap.start_arc == ai0 && second_submap.end_arc == ai0);

    StartupOracle oracle1(&first_curve, 4.0, 3);
    StartupOracle oracle2(&second_curve, 4.0, 0);

    FusionState state;
    build_fusion_sequence(state, first_submap, first_curve);

    std::size_t start = fusion_startup(state, first_submap, first_curve, second_submap,
                                       second_curve, oracle1, oracle2);

    assert(start == 1);
    assert(!state.chords.empty());

    std::printf("  [PASS] startup_d1_eq_d2_defaults_to_case1\n");
}

static void test_build_fusion_sequence_skips_null_length_chords() {
    auto curve = make_C1();

    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r_null = s.add_node();
    std::size_t r_cap = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a{};

    a.first_edge = 2;
    a.last_edge = 2;
    a.first_side = LEFT;
    a.last_side = LEFT;
    a.region_node = r_null;
    a.edge_count = 0;
    std::size_t N = s.add_arc(a);

    a = {};
    a.first_edge = 2;
    a.last_edge = 2;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = r0;
    a.edge_count = 2;
    std::size_t arc2 = s.add_arc(a);

    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r_cap;
    a.edge_count = 0;
    std::size_t Z = s.add_arc(a);

    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = r0;
    a.edge_count = 2;
    std::size_t arc1 = s.add_arc(a);

    Chord ec;
    ec.region[0] = r0;
    ec.region[1] = r_cap;
    ec.left_edge = 2;
    ec.left_side = RIGHT;
    ec.left_adj = {{arc2}, 1};
    ec.right_edge = 1;
    ec.right_side = RIGHT;
    ec.right_adj = {{Z}, 1};
    ec.y = curve.vertex(2).y;
    ec.y_tag = 2;
    s.add_chord(ec);

    Chord null_chord;
    null_chord.region[0] = r0;
    null_chord.region[1] = r_null;
    null_chord.left_edge = 2;
    null_chord.right_edge = 2;
    null_chord.left_side = LEFT;
    null_chord.right_side = LEFT;
    null_chord.y = curve.vertex(2).y;
    null_chord.y_tag = 2;
    null_chord.is_null_length = true;
    null_chord.left_adj = {{arc1}, 1};
    null_chord.right_adj = {{N}, 1};
    s.add_chord(null_chord);

    assert(s.start_arc == arc1 && s.end_arc == arc2);

    FusionState state;
    build_fusion_sequence(state, s, curve);

    assert(state.sequence.size() == 4 &&
           "[C91 §3.1 tex 179]: null-length chords must not appear in the "
           "canonical vertex enumeration (paper distinguishes 'exit chord "
           "endpoints' from null-length chords per [C91 §2.2 tex 96])");

    std::size_t null_chord_idx = 1;
    for (std::size_t i = 0; i < state.sequence.size(); ++i) {
        const auto& v = state.sequence[i];
        if (v.is_companion)
            continue;
        assert(v.chord_idx != null_chord_idx);
    }

    std::printf("  [PASS] build_fusion_sequence_skips_null_length_chords\n");
}

static void test_startup_vertex_to_vertex_tie_break() {
    Polygon Pfx({{0, 0, 0}, {1, 2, 1}, {2, 4, 2}, {3, 1, 3}, {4, 3, 4}, {5, 5, 5}, {6, 3, 6}});
    Polygon first_curve = Pfx.subchain(0, 5);
    auto first_submap = make_S1(first_curve);

    Polygon second_curve = Pfx.subchain(4, 3);

    Submap second_submap;
    std::size_t r_outer = second_submap.add_node();
    std::size_t r_pocket = second_submap.add_node();
    second_submap.start_vertex = 0;
    second_submap.end_vertex = 2;

    Arc a{};
    a.first_edge = 1;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r_pocket;
    a.edge_count = 2;
    std::size_t B = second_submap.add_arc(a);

    a = {};
    a.first_edge = 0;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r_outer;
    a.edge_count = 2;
    std::size_t W = second_submap.add_arc(a);

    Chord c{};
    c.region[0] = r_outer;
    c.region[1] = r_pocket;
    c.left_edge = 0;
    c.left_side = RIGHT;
    c.right_edge = 1;
    c.right_side = RIGHT;
    c.y = 3.0;
    c.y_tag = 4;
    c.left_adj = {{B}, 1};
    c.right_adj = {{W}, 1};
    second_submap.add_chord(c);

    assert(second_submap.start_arc == W && second_submap.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

    StartupOracle oracle1(&first_curve, 20.0);
    StartupOracle oracle2(&second_curve, 6.0);

    FusionState state;
    build_fusion_sequence(state, first_submap, first_curve);

    std::size_t start = fusion_startup(state, first_submap, first_curve, second_submap,
                                       second_curve, oracle1, oracle2);

    assert(start == 1);

    assert(state.current_region == r_outer &&
           "vertex-to-vertex tie-break: leaving_downward=true should "
           "elect the region below the chord (r_outer)");
    assert(!state.chords.empty());

    std::printf("  [PASS] startup_vertex_to_vertex_tie_break\n");
}

static void test_startup_mid_edge_tie_break() {
    auto first_curve = make_C1();
    auto first_submap = make_S1(first_curve);

    Polygon second_curve({{4, 3, 4}, {5, 5, 5}, {6, 1, 6}});

    Submap second_submap;
    std::size_t r_above = second_submap.add_node();
    std::size_t r_below = second_submap.add_node();
    second_submap.start_vertex = 0;
    second_submap.end_vertex = 2;

    Arc a{};
    a.first_edge = 1;
    a.last_edge = 0;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r_above;
    a.edge_count = 2;
    std::size_t U = second_submap.add_arc(a);

    a = {};
    a.first_edge = 0;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = RIGHT;
    a.region_node = r_below;
    a.edge_count = 2;
    std::size_t W = second_submap.add_arc(a);

    Chord c{};
    c.region[0] = r_above;
    c.region[1] = r_below;
    c.left_edge = 0;
    c.left_side = RIGHT;
    c.right_edge = 1;
    c.right_side = RIGHT;
    c.y = 3.0;
    c.y_tag = 4;
    c.left_adj = {{U}, 1};
    c.right_adj = {{W, U}, 2};
    second_submap.add_chord(c);

    assert(second_submap.start_arc == W && second_submap.end_arc == W &&
           "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

    StartupOracle oracle1(&first_curve, 20.0);
    StartupOracle oracle2(&second_curve, 5.5);

    FusionState state;
    build_fusion_sequence(state, first_submap, first_curve);
    std::size_t start = fusion_startup(state, first_submap, first_curve, second_submap,
                                       second_curve, oracle1, oracle2);

    assert(start == 1);

    assert(state.current_region == r_below &&
           "mid-edge tie-break: leaving_downward=true should elect "
           "the region below the chord (r_below)");

    std::printf("  [PASS] startup_mid_edge_tie_break\n");
}

struct ForwardOracle : RayShootingOracle {
    const Polygon* curve;
    Exact forward_dist;
    ForwardOracle(const Polygon* c, Exact d) : curve(c), forward_dist(std::move(d)) {}
    RayHit shoot(Point p, Side dir, std::size_t, const Subarc& target,
                 SourceOffset = SOURCE_OFFSET_NONE) const override {
        RayHit h;

        SymbolicY py{p.y, p.index};
        ArcSideRange ranges[3];
        std::size_t n = subarc_side_ranges(target, 0, curve->num_vertices() - 1, ranges);
        for (std::size_t li = 0; li < n; ++li) {
            for (std::size_t k = 0; k <= ranges[li].last_edge - ranges[li].first_edge; ++k) {
                std::size_t e = (ranges[li].side == LEFT) ? ranges[li].first_edge + k
                                                          : ranges[li].last_edge - k;
                const auto& ed = curve->edge(e);
                SymbolicY y0 = symbolic_y_of(curve->vertex(ed.start_idx));
                SymbolicY y1 = symbolic_y_of(curve->vertex(ed.end_idx));
                SymbolicY lo = symbolic_y_less(y0, y1) ? y0 : y1;
                SymbolicY hi = symbolic_y_less(y0, y1) ? y1 : y0;

                if (!subarc_contains_point(target, *curve, e, ranges[li].side, py, 0,
                                           curve->num_vertices() - 1))
                    continue;
                if (symbolic_y_leq(lo, py) && symbolic_y_leq(py, hi)) {
                    h.hit = true;
                    h.y = p.y;
                    h.x = (dir == LEFT) ? (p.x - forward_dist) : (p.x + forward_dist);
                    h.edge = e;
                    h.side = ranges[li].side;
                    return h;
                }
            }
        }
        h.hit = false;
        return h;
    }
};

static void test_fuse_main_loop_smoke() {
    Polygon Pfx({{0, 0, 0}, {1, 2, 1}, {2, 4, 2}, {3, 1, 3}, {4, 3, 4}, {5, 5, 5}, {6, 1, 6}});
    Polygon first_curve = Pfx.subchain(0, 5);
    auto first_submap = make_S1(first_curve);
    Polygon second_curve = Pfx.subchain(4, 3);

    Submap second_submap;
    second_submap.add_node();
    second_submap.start_vertex = 0;
    second_submap.end_vertex = 2;
    Arc a{};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2;
    std::size_t ai0 = second_submap.add_arc(a);
    assert(second_submap.start_arc == ai0 && second_submap.end_arc == ai0);

    ForwardOracle oracle1(&first_curve, 5.0);
    ForwardOracle oracle2(&second_curve, 1.0);

    FusionState state;
    fuse_submaps(state, first_submap, first_curve, second_submap, second_curve, oracle1, oracle2);

    assert(!state.chords.empty());
    std::printf("  [PASS] fuse_main_loop_smoke (chords=%zu)\n", state.chords.size());
}

static Submap make_peak_S2() {
    Submap second_submap;
    std::size_t r_above = second_submap.add_node();
    std::size_t r_below = second_submap.add_node();
    second_submap.start_vertex = 0;
    second_submap.end_vertex = 2;

    Arc a{};
    a.first_edge = 1;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = r_below;
    a.edge_count = 2;
    std::size_t input_curve = second_submap.add_arc(a);

    a = {};
    a.first_edge = 0;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = r_above;
    a.edge_count = 2;
    std::size_t Q = second_submap.add_arc(a);

    Chord c{};
    c.region[0] = r_above;
    c.region[1] = r_below;
    c.left_edge = 1;
    c.left_side = LEFT;
    c.right_edge = 0;
    c.right_side = RIGHT;
    c.y = 3.0;
    c.y_tag = 4;
    c.left_adj = {{Q, input_curve}, 2};
    c.right_adj = {{input_curve}, 1};
    second_submap.add_chord(c);

    assert(second_submap.start_arc == Q && second_submap.end_arc == input_curve);
    return second_submap;
}

static void test_fuse_main_loop_case_ii_smoke() {
    Polygon Pfx({{0, 0, 0}, {1, 2, 1}, {2, 4, 2}, {3, 1, 3}, {4, 3, 4}, {3, 5, 5}, {2, 1, 6}});
    Polygon first_curve = Pfx.subchain(0, 5);
    auto first_submap = make_S1(first_curve);
    Polygon second_curve = Pfx.subchain(4, 3);
    Submap second_submap = make_peak_S2();

    ForwardOracle oracle1(&first_curve, 5.0);
    ForwardOracle oracle2(&second_curve, 1.5);

    FusionState state;
    fuse_submaps(state, first_submap, first_curve, second_submap, second_curve, oracle1, oracle2);

    assert(!state.chords.empty());
    std::printf("  [PASS] fuse_main_loop_case_ii_smoke (chords=%zu)\n", state.chords.size());
}

static Submap make_chordless(const Polygon& poly) {
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
    return s;
}

static void test_rebuild_no_chords() {
    Polygon Pfx({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}});
    Polygon first_curve = Pfx.subchain(0, 2);
    Polygon second_curve = Pfx.subchain(1, 2);
    Polygon curve({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}});
    Submap first_submap = make_chordless(first_curve);
    Submap second_submap = make_chordless(second_curve);

    FusionState st1, st2;
    st1.invalidated_first_chords.assign(first_submap.num_chords(), false);
    st1.invalidated_second_chords.assign(second_submap.num_chords(), false);
    st2.invalidated_first_chords.assign(second_submap.num_chords(), false);
    st2.invalidated_second_chords.assign(first_submap.num_chords(), false);
    Submap out;
    rebuild_submap(out, curve, first_submap, first_curve, second_submap, second_curve, st1, st2);

    assert(out.num_nodes() == 1);
    assert(out.num_chords() == 0);

    assert(out.num_arcs() == 1);
    assert(out.arc(0).first_side == LEFT && out.arc(0).last_side == RIGHT &&
           out.arc(0).first_edge == 0 && out.arc(0).last_edge == 0);
    assert(out.start_arc == 0 && out.end_arc == 0);
    out.assert_tree_property();
    out.check_invariants(curve);

    std::printf("  [PASS] rebuild_no_chords\n");
}

static void test_rebuild_junction_extremum_inside_right() {
    Polygon Pfx({{-1, 0, 0}, {0, 1, 1}, {1, 0, 2}});
    Polygon first_curve = Pfx.subchain(0, 2);
    Polygon second_curve = Pfx.subchain(1, 2);
    Polygon curve({{-1, 0, 0}, {0, 1, 1}, {1, 0, 2}});
    Submap first_submap = make_chordless(first_curve);
    Submap second_submap = make_chordless(second_curve);

    FusionState st1, st2;
    st1.invalidated_first_chords.assign(first_submap.num_chords(), false);
    st1.invalidated_second_chords.assign(second_submap.num_chords(), false);
    st2.invalidated_first_chords.assign(second_submap.num_chords(), false);
    st2.invalidated_second_chords.assign(first_submap.num_chords(), false);
    Submap out;
    rebuild_submap(out, curve, first_submap, first_curve, second_submap, second_curve, st1, st2);

    assert(out.num_nodes() == 2);
    assert(out.num_chords() == 1);
    const Chord& nc = out.chord(0);
    assert(nc.is_null_length);
    assert(nc.y_tag == 1 && "[C91 §2 tex 47]: inside pair duplicates the junction vertex — "
                            "its symbolic y IS the vertex's (y, index)");
    assert(nc.left_edge == 1 && nc.right_edge == 1);
    assert(nc.left_side == RIGHT && nc.right_side == RIGHT &&
           "[C91 §2.1 tex 72]: inside of the turn is the RIGHT face here");
    assert(out.region_weight(nc.region[1]) == 0 &&
           "[C91 §2.2 tex 106]: null-length chord's inner region is empty");
    out.assert_tree_property();
    out.check_invariants(curve);

    std::printf("  [PASS] rebuild_junction_extremum_inside_right\n");
}

static void test_rebuild_junction_extremum_inside_left() {
    Polygon Pfx({{1, 0, 0}, {0, 1, 1}, {-1, 0, 2}});
    Polygon first_curve = Pfx.subchain(0, 2);
    Polygon second_curve = Pfx.subchain(1, 2);
    Polygon curve({{1, 0, 0}, {0, 1, 1}, {-1, 0, 2}});
    Submap first_submap = make_chordless(first_curve);
    Submap second_submap = make_chordless(second_curve);

    FusionState st1, st2;
    st1.invalidated_first_chords.assign(first_submap.num_chords(), false);
    st1.invalidated_second_chords.assign(second_submap.num_chords(), false);
    st2.invalidated_first_chords.assign(second_submap.num_chords(), false);
    st2.invalidated_second_chords.assign(first_submap.num_chords(), false);
    Submap out;
    rebuild_submap(out, curve, first_submap, first_curve, second_submap, second_curve, st1, st2);

    assert(out.num_nodes() == 2);
    assert(out.num_chords() == 1);
    const Chord& nc = out.chord(0);
    assert(nc.is_null_length);
    assert(nc.y_tag == 1);
    assert(nc.left_edge == 1 && nc.right_edge == 1);
    assert(nc.left_side == LEFT && nc.right_side == LEFT &&
           "[C91 §2.1 tex 72]: inside of the turn is the LEFT face here");
    assert(out.region_weight(nc.region[1]) == 0);
    out.assert_tree_property();
    out.check_invariants(curve);

    std::printf("  [PASS] rebuild_junction_extremum_inside_left\n");
}

static void test_rebuild_discovered_chord_frames() {
    Polygon Pfx({{0, 0, 0}, {1, 2, 1}, {2, 1, 2}, {3, 3, 3}});
    Polygon first_curve = Pfx.subchain(0, 3);
    Polygon second_curve = Pfx.subchain(2, 2);
    Polygon curve({{0, 0, 0}, {1, 2, 1}, {2, 1, 2}, {3, 3, 3}});
    Submap first_submap = make_chordless(first_curve);
    Submap second_submap = make_chordless(second_curve);

    auto run = [&](bool via_state1, bool wrong_wall = false) -> Submap {
        FusionState st1, st2;
        st1.invalidated_first_chords.assign(first_submap.num_chords(), false);
        st1.invalidated_second_chords.assign(second_submap.num_chords(), false);
        st2.invalidated_first_chords.assign(second_submap.num_chords(), false);
        st2.invalidated_second_chords.assign(first_submap.num_chords(), false);
        FusionState::DiscoveredChord dc;
        dc.y = SymbolicY{2.0, 1};
        dc.left_edge = 0;
        dc.left_side = LEFT;
        dc.right_edge = 0;
        dc.right_side = wrong_wall ? LEFT : RIGHT;
        if (via_state1) {
            dc.left_on_first_curve = true;
            dc.right_on_first_curve = false;
            st1.chords.push_back(dc);
        } else {
            dc.left_on_first_curve = false;
            dc.right_on_first_curve = true;
            st2.chords.push_back(dc);
        }
        Submap out;
        rebuild_submap(out, curve, first_submap, first_curve, second_submap, second_curve, st1,
                       st2);
        return out;
    };

    for (bool via_state1 : {true, false}) {
        Submap out = run(via_state1);

        assert(out.num_nodes() == 3);
        assert(out.num_chords() == 2);
        std::size_t discovered = NONE, null_c = NONE;
        for (std::size_t ci = 0; ci < out.num_chords(); ++ci)
            (out.chord(ci).is_null_length ? null_c : discovered) = ci;
        assert(discovered != NONE && null_c != NONE);

        const Chord& c = out.chord(discovered);
        assert(c.left_edge == 0 && "input-curve edge translation: C₁ endpoint → C's edge 0");
        assert(c.right_edge == 2 && "input-curve edge translation: C₂ endpoint → C's edge 2");
        assert(c.y_tag == 1);

        assert(c.left_adj.count == 1 && c.right_adj.count == 2);

        const Chord& nc = out.chord(null_c);
        assert(nc.y_tag == 2 && "[C91 §2 tex 47]: junction null chord carries the junction "
                                "vertex's own SoS tag");
        assert(nc.left_edge == 2 && nc.left_side == LEFT &&
               "[C91 §2.1 tex 72]: inside of the min turn is the LEFT "
               "face here (previous branch to the left)");
        assert(out.region_weight(nc.region[1]) == 0);

        out.assert_tree_property();
        out.check_invariants(curve);
    }

    require_assertion_abort([&] { (void)run(true, true); });
    std::printf("  [PASS] rebuild_discovered_chord_frames\n");
}

static void test_rebuild_junction_null_in_input_fires() {
    require_assertion_abort([] {
        Polygon Pfx({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}});
        Polygon first_curve = Pfx.subchain(0, 2);
        Polygon second_curve = Pfx.subchain(1, 2);
        Polygon curve({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}});
        Submap first_submap = make_chordless(first_curve);

        Submap second_submap;
        std::size_t r0 = second_submap.add_node();
        std::size_t r1 = second_submap.add_node();
        second_submap.start_vertex = 0;
        second_submap.end_vertex = 1;
        Arc a{};

        a.first_edge = 0;
        a.last_edge = 0;
        a.first_side = LEFT;
        a.last_side = RIGHT;
        a.region_node = r0;
        a.edge_count = 1;
        std::size_t a0 = second_submap.add_arc(a);
        a = {};
        a.first_edge = 0;
        a.last_edge = 0;
        a.first_side = RIGHT;
        a.last_side = LEFT;
        a.region_node = r1;
        a.edge_count = 0;
        std::size_t a_null = second_submap.add_arc(a);
        Chord nc{};
        nc.region[0] = r0;
        nc.region[1] = r1;
        nc.left_edge = 0;
        nc.right_edge = 0;
        nc.left_side = LEFT;
        nc.right_side = LEFT;
        nc.is_null_length = true;
        nc.y = 1.0;
        nc.y_tag = 1;
        nc.left_adj = {{a0}, 1};
        nc.right_adj = {{a_null}, 1};
        second_submap.add_chord(nc);
        (void)a0;
        (void)a_null;

        FusionState st1, st2;
        st1.invalidated_first_chords.assign(first_submap.num_chords(), false);
        st1.invalidated_second_chords.assign(second_submap.num_chords(), false);
        st2.invalidated_first_chords.assign(second_submap.num_chords(), false);
        st2.invalidated_second_chords.assign(first_submap.num_chords(), false);
        Submap out;
        rebuild_submap(out, curve, first_submap, first_curve, second_submap, second_curve, st1,
                       st2);
    });
    std::printf("  [PASS] rebuild_junction_null_in_input_fires\n");
}

static void test_case_ii_hit_beyond_ab_disqualified() {
    Polygon Pfx({{0, 0, 0}, {1, 2, 1}, {2, 4, 2}, {3, 1, 3}, {4, 3, 4}, {3, 5, 5}, {2, 1, 6}});
    Polygon first_curve = Pfx.subchain(0, 5);
    auto first_submap = make_S1(first_curve);
    Polygon second_curve = Pfx.subchain(4, 3);
    Submap second_submap = make_peak_S2();

    ForwardOracle oracle1(&first_curve, 100.0);
    ForwardOracle oracle2(&second_curve, 1.5);

    FusionState state;
    fuse_submaps(state, first_submap, first_curve, second_submap, second_curve, oracle1, oracle2);

    assert(state.chords.size() == 4 && "[C91 §3.1 tex 222]: hits beyond ab must be disqualified — "
                                       "startup + a₁/a₂/a_{m+1} case (i) only");
    for (const auto& dc : state.chords) {
        if (!(dc.y.y == 3.0 && dc.y.tag == 4))
            continue;
        assert(dc.left_on_first_curve != dc.right_on_first_curve);
        std::size_t first_curve_edge = dc.left_on_first_curve ? dc.left_edge : dc.right_edge;
        assert(first_curve_edge == 3 && "[C91 §3.1 tex 222]: chords at the S₂ chord's y must be "
                                        "junction-companion records, not case (ii) products");
    }

    std::printf("  [PASS] case_ii_hit_beyond_ab_disqualified\n");
}

using GeomOracle = TestArcRayShooter;

static void test_case_ii_fires() {
    Polygon Pfx({{0, 0, 0}, {1, 2, 1}, {2, 6, 2}, {3, 0.5, 3}, {4, 3, 4}, {3, 5, 5}, {2, 1, 6}});
    Polygon first_curve = Pfx.subchain(0, 5);
    auto first_submap = make_S1(first_curve);
    Polygon second_curve = Pfx.subchain(4, 3);

    Submap second_submap;
    std::size_t r_pocket = second_submap.add_node();
    std::size_t r_out = second_submap.add_node();
    second_submap.start_vertex = 0;
    second_submap.end_vertex = 2;
    Arc a{};
    a.first_edge = 0;
    a.first_side = LEFT;
    a.last_edge = 1;
    a.last_side = LEFT;
    a.region_node = r_pocket;
    a.edge_count = 2;
    std::size_t A_pocket = second_submap.add_arc(a);
    a = {};
    a.first_edge = 1;
    a.first_side = LEFT;
    a.last_edge = 0;
    a.last_side = LEFT;
    a.region_node = r_out;
    a.edge_count = 2;
    std::size_t A_out = second_submap.add_arc(a);
    assert(second_submap.start_arc == A_out && second_submap.end_arc == A_out);
    Chord cc{};
    cc.region[0] = r_pocket;
    cc.region[1] = r_out;
    cc.left_edge = 1;
    cc.left_side = LEFT;
    cc.right_edge = 0;
    cc.right_side = LEFT;
    cc.y = 3.0;
    cc.y_tag = 4;
    cc.left_adj = {{A_pocket, A_out}, 2};
    cc.right_adj = {{A_out}, 1};
    second_submap.add_chord(cc);

    GeomOracle oracle1(&first_curve);
    GeomOracle oracle2(&second_curve);

    FusionState state;
    fuse_submaps(state, first_submap, first_curve, second_submap, second_curve, oracle1, oracle2);

    bool case_ii_product = false;
    bool case_ii_product_a = false;
    for (const auto& dc : state.chords) {
        if (!(dc.y.y == 3.0 && dc.y.tag == 4))
            continue;
        if (dc.left_on_first_curve == dc.right_on_first_curve)
            continue;
        std::size_t first_curve_edge = dc.left_on_first_curve ? dc.left_edge : dc.right_edge;
        if (first_curve_edge == 3)
            continue;
        if (dc.left_on_first_curve && dc.left_edge == 2 && dc.left_side == LEFT &&
            !dc.right_on_first_curve && dc.right_edge == 0 && dc.right_side == LEFT) {
            case_ii_product = true;
        } else if (!dc.left_on_first_curve && dc.left_edge == 1 && dc.left_side == LEFT &&
                   dc.right_on_first_curve && dc.right_edge == 2 && dc.right_side == RIGHT) {
            case_ii_product_a = true;
        } else {
            assert(false && "[C91 §3.1 tex 206]: only the two true ab-endpoint "
                            "chords may be recorded at ab's level");
        }
    }
    assert(case_ii_product_a && "[C91 §3.1 tex 206/222]: a's eastward chord to C₁ edge 2 "
                                "RIGHT must be recorded");
    assert(case_ii_product && "[C91 §3.1 tex 202/206]: an on-ab, after-p, back-shot-"
                              "confirmed candidate must fire case (ii) and record the "
                              "exit-chord-endpoint → p' chord");

    std::printf("  [PASS] case_ii_fires (chords=%zu)\n", state.chords.size());
}

static void test_fusion_sequence_junction_at_start() {
    auto first_curve = make_C1();
    auto first_submap = make_S1(first_curve);

    FusionState state;
    state.junction_at_end = false;
    build_fusion_sequence(state, first_submap, first_curve);
    const auto& seq = state.sequence;

    assert(seq.size() == 4);
    assert(seq.front().is_companion && seq.front().side == LEFT &&
           "[C91 §3.1 tex 179]: junction-at-start tour begins at the "
           "LEFT companion");
    assert(seq.front().edge == 0);
    assert(seq.back().is_companion && seq.back().side == RIGHT &&
           "[C91 §3.1 tex 179]: junction-at-start tour ends at the "
           "RIGHT companion");
    assert(seq.back().edge == 0);

    SymbolicY jy = symbolic_y_of(first_curve.vertex(0));
    assert(symbolic_y_equal(seq.front().y, jy));
    assert(symbolic_y_equal(seq.back().y, jy));

    assert(seq[1].side == LEFT);
    assert(seq[2].side == RIGHT);

    std::printf("  [PASS] fusion_sequence_junction_at_start\n");
}

static void test_startup_case1_junction_at_start() {
    Polygon W({{4, 3, 4}, {5, 5, 5}, {6, 1, 6}});
    auto T = make_C1();
    auto S_T = make_S1(T);

    Submap S_W;
    S_W.add_node();
    S_W.start_vertex = 0;
    S_W.end_vertex = 2;
    Arc a{};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2;
    std::size_t ai0 = S_W.add_arc(a);
    assert(S_W.start_arc == ai0 && S_W.end_arc == ai0);

    StartupOracle oracleW(&W, -10.0);
    StartupOracle oracleT(&T, 3.0);

    FusionState state;
    state.junction_at_end = false;
    build_fusion_sequence(state, S_W, W);
    std::size_t start = fusion_startup(state, S_W, W, S_T, T, oracleW, oracleT);

    assert(start == 1);

    assert(state.current_region == S_T.arc(S_T.end_arc).region_node);
    assert(state.current_region == 1);
    assert(!state.chords.empty());

    const auto& dc = state.chords[0];
    assert(dc.left_on_first_curve || dc.right_on_first_curve);

    std::printf("  [PASS] startup_case1_junction_at_start\n");
}

static void test_fuse_main_loop_smoke_junction_at_start() {
    Polygon Pfx(
        {{-2, 5, 8}, {-1, 1, 9}, {0, 0, 10}, {1, 2, 11}, {2, 4, 12}, {3, 1, 13}, {4, 3, 14}});
    Polygon T = Pfx.subchain(0, 3);
    Polygon W = Pfx.subchain(2, 5);

    Submap S_W;
    S_W.add_node();
    S_W.add_node();
    S_W.start_vertex = 0;
    S_W.end_vertex = 4;
    Arc a{};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 1;
    a.edge_count = 3;
    std::size_t aiE = S_W.add_arc(a);
    a = {};
    a.first_edge = 1;
    a.last_edge = 1;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = 0;
    a.edge_count = 2;
    std::size_t aiS = S_W.add_arc(a);
    Chord c{};
    c.region[0] = 0;
    c.region[1] = 1;
    c.left_edge = 1;
    c.right_edge = 1;
    c.left_side = LEFT;
    c.right_side = RIGHT;
    c.y = W.vertex(2).y;
    c.y_tag = W.vertex(2).index;
    c.left_adj = {{aiS}, 1};
    c.right_adj = {{aiE}, 1};
    S_W.add_chord(c);
    assert(S_W.start_arc == aiS && S_W.end_arc == aiE);

    Submap S_T;
    S_T.add_node();
    S_T.start_vertex = 0;
    S_T.end_vertex = 2;
    a = {};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 2;
    std::size_t t0 = S_T.add_arc(a);
    assert(S_T.start_arc == t0 && S_T.end_arc == t0);

    ForwardOracle oracleW(&W, 5.0);
    ForwardOracle oracleT(&T, 1.0);

    FusionState state;
    state.junction_at_end = false;
    fuse_submaps(state, S_W, W, S_T, T, oracleW, oracleT);

    assert(!state.chords.empty());
    std::printf("  [PASS] fuse_main_loop_smoke_junction_at_start "
                "(chords=%zu)\n",
                state.chords.size());
}

static void test_rebuild_dedup() {
    Polygon Pfx({{0, 0, 0}, {1, 2, 1}, {2, 1, 2}, {3, 3, 3}});
    Polygon first_curve = Pfx.subchain(0, 3);
    Polygon second_curve = Pfx.subchain(2, 2);
    Polygon curve({{0, 0, 0}, {1, 2, 1}, {2, 1, 2}, {3, 3, 3}});
    Submap first_submap = make_chordless(first_curve);
    Submap second_submap = make_chordless(second_curve);

    FusionState st1, st2;
    st1.invalidated_first_chords.assign(first_submap.num_chords(), false);
    st1.invalidated_second_chords.assign(second_submap.num_chords(), false);
    st2.invalidated_first_chords.assign(second_submap.num_chords(), false);
    st2.invalidated_second_chords.assign(first_submap.num_chords(), false);
    st2.junction_at_end = false;
    FusionState::DiscoveredChord dc;
    dc.y = SymbolicY{2.0, 1};
    dc.left_edge = 0;
    dc.left_side = LEFT;
    dc.right_edge = 0;
    dc.right_side = RIGHT;
    dc.left_on_first_curve = true;
    dc.right_on_first_curve = false;
    st1.chords.push_back(dc);
    dc.left_on_first_curve = false;
    dc.right_on_first_curve = true;
    st2.chords.push_back(dc);

    Submap out;
    rebuild_submap(out, curve, first_submap, first_curve, second_submap, second_curve, st1, st2);

    assert(out.num_chords() == 2 && "[C91 §3.1 tex 224]: duplicate chords must be deduplicated");
    assert(out.num_nodes() == 3);
    out.assert_tree_property();
    out.check_invariants(curve);

    std::printf("  [PASS] rebuild_dedup\n");
}

static void test_rebuild_dedup_junction_cross_edge_labels() {
    Polygon input_curve({{0, 0, 0}, {2, 2, 1}, {4, 4, 2}, {0, 5, 3}, {-1, 1, 4}});
    Polygon first_curve = input_curve.subchain(0, 2);
    Polygon second_curve = input_curve.subchain(1, 4);
    Submap first_submap = make_chordless(first_curve);
    Submap second_submap = make_chordless(second_curve);

    GeomOracle o1(&first_curve);
    GeomOracle o2(&second_curve);

    FusionState st1;
    st1.junction_at_end = true;
    fuse_submaps(st1, first_submap, first_curve, second_submap, second_curve, o1, o2);

    FusionState st2;
    st2.junction_at_end = false;
    fuse_submaps(st2, second_submap, second_curve, first_submap, first_curve, o2, o1);

    Polygon curve(first_curve, second_curve);
    Submap out;
    rebuild_submap(out, curve, first_submap, first_curve, second_submap, second_curve, st1, st2);

    assert(out.num_chords() == 2 && "[C91 §3.1 tex 224]: cross-pass junction records must dedup "
                                    "(same ∂C point, different incident-edge labels)");
    assert(out.num_nodes() == 3);
    out.assert_tree_property();
    out.check_invariants(curve);

    std::printf("  [PASS] rebuild_dedup_junction_cross_edge_labels\n");
}

int main() {
    std::setbuf(stdout, nullptr);
    std::printf("[C91 §3.1 tests]:\n");
    test_fusion_sequence_basic();
    test_fusion_sequence_no_chords();
    test_fusion_sequence_ordering();
    test_companion_identity();
    test_collect_region_arcs();
    test_local_shoot();
    test_local_shoot_nearest();
    test_startup_case1();
    test_startup_case2();
    test_shooting_direction_all_cases();
    test_ray_contact_tie_break();
    test_startup_d1_eq_d2_defaults_to_case1();
    test_build_fusion_sequence_skips_null_length_chords();
    test_startup_vertex_to_vertex_tie_break();
    test_startup_mid_edge_tie_break();
    test_fuse_main_loop_smoke();
    test_fuse_main_loop_case_ii_smoke();
    test_rebuild_no_chords();
    test_rebuild_junction_extremum_inside_right();
    test_rebuild_junction_extremum_inside_left();
    test_rebuild_discovered_chord_frames();
    test_rebuild_junction_null_in_input_fires();
    test_case_ii_hit_beyond_ab_disqualified();
    test_case_ii_fires();
    test_fusion_sequence_junction_at_start();
    test_startup_case1_junction_at_start();
    test_fuse_main_loop_smoke_junction_at_start();
    test_rebuild_dedup();
    test_rebuild_dedup_junction_cross_edge_labels();
    std::printf("All §3.1 tests passed.\n");
    return 0;
}
