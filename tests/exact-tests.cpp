#include "algorithm/polygon/polygon.h"
#include "algorithm/visibility/chain.h"
#include "algorithm/visibility/naive_visibility.h"
#include "algorithm/visibility/up_phase.h"
#include "support/arc_ray_shooter.h"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <type_traits>
#include <utility>

using namespace chazelle;

static_assert(!std::is_convertible_v<Exact, double>);
static_assert(!std::is_constructible_v<Exact, long double>);

static void test_rational_arithmetic() {
    const Exact third = Exact{1} / 3;
    assert(third * 3 == 1 && third + third + third == 1);
    const Exact largest_integer{std::numeric_limits<std::uint64_t>::max()};
    assert(largest_integer + 1 == Exact{std::uint64_t{1} << 63} * 2);

    Exact a = third;
    a += a;
    assert(a == third * 2);
    a *= a;
    assert(a == Exact{4} / 9);
    Exact b = std::move(a);
    assert(b == Exact{4} / 9);
    a = b;
    a /= a;
    assert(a == 1 && b == Exact{4} / 9);
}

static void test_exact_binary64_import() {
    const double adjacent = std::nextafter(1.0, 2.0);
    const Exact ulp = Exact{adjacent} - 1;
    assert(ulp == Exact{1} / (Exact{std::uint64_t{1} << 52}));
    const Exact halfway = exact_midpoint(1, adjacent);
    assert(halfway > 1 && halfway < adjacent && halfway - 1 == ulp / 2);

    const Exact tiny{std::numeric_limits<double>::denorm_min()};
    assert(tiny > 0 && tiny / 2 > 0 && (tiny / 2) * 2 == tiny);
    const Exact huge{std::numeric_limits<double>::max()};
    assert(huge * 2 > huge && (huge * 2) / 2 == huge);
}

static void test_exact_subdivision_and_crossing() {
    const double adjacent = std::nextafter(1.0, 2.0);
    const std::vector<Point> input{{1, 0, 0}, {adjacent, 2, 1}, {2, 1, 2}, {3, 3, 3}};
    const auto padded = pad_curve(input);
    const Point& mid = padded[1];
    const Exact determinant = (input[1].x - input[0].x) * (mid.y - input[0].y) -
                              (input[1].y - input[0].y) * (mid.x - input[0].x);
    assert(determinant == 0 && mid.x > input[0].x && mid.x < input[1].x);

    Polygon c(input);
    Exact crossing;
    assert(edge_crossing_x(c, 0, {1, 2}, &crossing));
    assert(crossing == mid.x);
    assert(edge_x_at_y(c, 0, {1, 2}) == crossing);
    UpPhase up(input);
    assert(up.chain_submap(up.graded().maximum_grade(), 0)
               .is_granular(UpPhase::grade_granularity(up.graded().maximum_grade()),
                            up.graded().curve()));

    std::vector<Point> diagonal;
    double q = 1.0;
    for (std::size_t i = 0; i < 4; ++i) {
        diagonal.push_back({q, q, i});
        q = std::nextafter(q, 2.0);
    }
    UpPhase diagonal_up(diagonal);
    const Polygon& d = diagonal_up.graded().curve();
    for (std::size_t i = 0; i < d.num_vertices(); ++i)
        assert(d.vertex(i).x == d.vertex(i).y);
    assert(d.vertex(0).x < d.vertex(1).x && d.vertex(1).x < d.vertex(2).x);
}

static void test_exact_branch_order() {
    const double huge = std::numeric_limits<double>::max();
    const Point u{huge, huge, 0}, v{0, 0, 1};
    const Point w{std::nextafter(huge, 0.0), huge, 2};

    assert(!extremum_prev_branch_left(u, v, w));
    const double tiny = std::numeric_limits<double>::denorm_min();
    const Point a{tiny, tiny, 0}, b{0, 0, 1}, c{2 * tiny, tiny, 2};
    assert(extremum_prev_branch_left(a, b, c));
}

static void test_free_source_at_tied_level() {
    Polygon c({{0, 0, 0}, {1, 1, 1}, {2, 0, 2}});
    Submap s = build_canonical_submap_naive(c);
    RayShootingStructure structure(s, c, 2);
    TestArcRayShooter arc_shooter(&c);
    const Subarc target{0, LEFT, 1, LEFT, {0, 0}, {0, 2}};

    const Point p{1, 1, SOS_NONE};
    for (Side direction : {LEFT, RIGHT}) {
        const std::size_t expected_edge = direction == LEFT ? 0 : 1;
        const RayHit hits[] = {naive_first_contact(c, p, {1, SOS_NONE}, direction),
                               structure.shoot_toward_boundary(p, direction),
                               arc_shooter.shoot(p, direction, 0, target)};
        for (const RayHit& hit : hits)
            assert(hit.hit && !hit.wrapped && hit.x == 1 && hit.edge == expected_edge &&
                   hit.side == RIGHT);
    }
}

int main() {
    test_rational_arithmetic();
    test_exact_binary64_import();
    test_exact_subdivision_and_crossing();
    test_exact_branch_order();
    test_free_source_at_tied_level();
    std::printf("All exact-geometry tests passed.\n");
}
