#include "support/triangulation_production_checks.h"
#include "support/unimonotone_fixture.h"

#include <algorithm>
#include <cstdio>
#include <numeric>

using namespace chazelle;

namespace {

void check_chain(std::vector<Point> points, bool clockwise, bool expect_backtracking) {
    UnimonotoneDecomposition partition;
    UnimonotonePolygon polygon;
    polygon.vertices.resize(points.size());
    std::iota(polygon.vertices.begin(), polygon.vertices.end(), 0);
    if (!clockwise)
        std::reverse(polygon.vertices.begin(), polygon.vertices.end());
    polygon.top_vertex = 0;
    polygon.bottom_vertex = 0;
    for (std::size_t vertex = 0; vertex < points.size(); ++vertex) {
        if (point_y_above(points[vertex], points[polygon.top_vertex]))
            polygon.top_vertex = vertex;
        if (point_y_below(points[vertex], points[polygon.bottom_vertex]))
            polygon.bottom_vertex = vertex;
    }
    partition.polygons.push_back(std::move(polygon));
    const auto result = triangulate_unimonotone(points, partition);
    test::check_production_triangulation(points, result);
    if (expect_backtracking)
        test::require_triangulation(result.work.backward_steps > points.size() / 8,
                                    "[FM84 Algorithm 3]: exercise a linear number of backtracks");
}

void check_long_chains() {
    constexpr std::size_t count = 200003;
    std::vector<Point> points;
    points.reserve(count);
    for (std::size_t vertex = 0; vertex < count; ++vertex)
        points.push_back({Exact(vertex), Exact(vertex) * Exact(vertex), vertex});
    check_chain(std::move(points), false, false);
    constexpr std::size_t reflex_count = 100003;
    points.clear();
    points.reserve(reflex_count);
    points.push_back({0, 0, 0});
    for (std::size_t vertex = 1; vertex + 1 < reflex_count; ++vertex)
        points.push_back({vertex % 2 == 0 ? 1 : 3, Exact(vertex), vertex});
    points.push_back({0, Exact(reflex_count), reflex_count - 1});
    check_chain(std::move(points), false, true);
    points.clear();
    points.reserve(reflex_count);
    points.push_back({0, Exact(reflex_count), 0});
    for (std::size_t vertex = 1; vertex + 1 < reflex_count; ++vertex)
        points.push_back({vertex % 2 == 0 ? 1 : 3, Exact(reflex_count - vertex / 2), vertex});
    points.push_back({0, 0, reflex_count - 1});
    check_chain(std::move(points), true, true);
}

}

int main() {
    const auto input = test::alternating_chains(50000);
    test::check_production_triangulation(input.vertices,
                                         triangulate_trapezoids(input.vertices, input.trapezoids));
    check_long_chains();
    std::puts("[FM84 Algorithm 3 production tests]: all passed");
}
