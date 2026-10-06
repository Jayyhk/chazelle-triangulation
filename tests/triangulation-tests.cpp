#include "algorithm/visibility/naive_visibility.h"
#include "support/assertions.h"
#include "support/triangulation_checks.h"
#include "support/unimonotone_fixture.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <numeric>
#include <random>

using namespace chazelle;

namespace {

UnimonotoneDecomposition one_polygon(std::span<const Point> points) {
    UnimonotoneDecomposition result;
    UnimonotonePolygon polygon;
    polygon.vertices.resize(points.size());
    std::iota(polygon.vertices.begin(), polygon.vertices.end(), 0);
    Exact area = 0;
    polygon.top_vertex = 0;
    polygon.bottom_vertex = 0;
    for (std::size_t vertex = 0; vertex < points.size(); ++vertex) {
        const Point& a = points[vertex];
        const Point& b = points[(vertex + 1) % points.size()];
        area += a.x * b.y - a.y * b.x;
        if (test::triangle_height(a) > test::triangle_height(points[polygon.top_vertex]))
            polygon.top_vertex = vertex;
        if (test::triangle_height(a) < test::triangle_height(points[polygon.bottom_vertex]))
            polygon.bottom_vertex = vertex;
    }
    if (area > 0)
        std::reverse(polygon.vertices.begin(), polygon.vertices.end());
    result.polygons.push_back(std::move(polygon));
    return result;
}

void check_endpoint_backtracking() {
    const std::vector<Point> points{{0, 8, 0}, {6, 6, 1}, {4, 4, 2}, {1, 2, 3}, {0, 0, 4}};
    const auto result = triangulate_unimonotone(points, one_polygon(points));
    test::check_triangulation(points, result);
    assert((result.triangles[0].vertices == std::array<std::size_t, 3>{0, 1, 2}));
    assert((result.triangles[1].vertices == std::array<std::size_t, 3>{0, 2, 3}));
    assert((result.triangles[2].vertices == std::array<std::size_t, 3>{0, 3, 4}));
    assert(test::triangle_orientation(points[4], points[0], points[3]) < 0 &&
           test::triangle_orientation(points[0], points[2], points[3]) < 0 &&
           test::triangle_orientation(points[2], points[4], points[3]) < 0 &&
           "[FM84 Algorithm 3 printed backtrack]: removing vertex 0 would enclose vertex 3");
}

void check_chain_orders() {
    const std::vector<std::vector<Point>> fixtures{
        {{0, 0, 0}, {4, 2, 1}, {0, 4, 2}},
        {{0, 0, 0}, {3, 1, 1}, {1, 2, 2}, {5, 3, 3}, {1, 4, 4}, {6, 5, 5}, {0, 6, 6}},
        {{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {0, 5, 4}},
        {{0, 0, 0}, {3, 1, 1}, {3, 2, 2}, {3, 3, 3}, {2, 4, 4}, {0, 5, 5}}};
    for (const auto& fixture : fixtures) {
        for (bool reverse : {false, true}) {
            for (bool reflect : {false, true}) {
                for (std::size_t first = 0; first < fixture.size(); ++first) {
                    auto points = fixture;
                    if (reverse)
                        std::reverse(points.begin(), points.end());
                    std::rotate(points.begin(), points.begin() + static_cast<std::ptrdiff_t>(first),
                                points.end());
                    for (std::size_t vertex = 0; vertex < points.size(); ++vertex) {
                        points[vertex].index = SOS_NONE - 97 + vertex;
                        if (reflect)
                            points[vertex].x = -points[vertex].x;
                    }
                    const auto partition = one_polygon(points);
                    const auto result = triangulate_unimonotone(points, partition);
                    test::check_triangulation(points, result);
                    auto rotated = partition;
                    auto& boundary = rotated.polygons[0].vertices;
                    std::rotate(boundary.begin(), boundary.begin() + 1, boundary.end());
                    const auto other = triangulate_unimonotone(points, rotated);
                    assert(other.triangles.size() == result.triangles.size());
                    for (std::size_t index = 0; index < result.triangles.size(); ++index)
                        assert(other.triangles[index].vertices == result.triangles[index].vertices);
                }
            }
        }
    }
    std::mt19937 random(1601984);
    std::size_t forward = 0;
    for (std::size_t sample = 0; sample < 250; ++sample) {
        const std::size_t count = 3 + random() % 25;
        std::vector<Point> points{{0, 0, 0}};
        for (std::size_t vertex = 1; vertex + 1 < count; ++vertex)
            points.push_back({Exact(1 + random() % 40) / 7, Exact(vertex), vertex});
        points.push_back({0, Exact(count), count - 1});
        if (sample % 2 == 0)
            std::reverse(points.begin(), points.end());
        for (std::size_t vertex = 0; vertex < count; ++vertex)
            points[vertex].index = vertex + 71;
        const auto result = triangulate_unimonotone(points, one_polygon(points));
        forward += result.work.forward_steps;
        test::check_triangulation(points, result);
    }
    assert(forward > 0);
}

void check_trapezoid_pipeline() {
    const std::vector<std::vector<Point>> fixtures{
        {{0, 0, 0}, {4, 0, 1}, {4, 4, 2}, {0, 4, 3}},
        {{0, 0, 0}, {2, 0, 1}, {4, 0, 2}, {4, 2, 3}, {4, 4, 4}, {0, 4, 5}},
        {{0, 0, 0}, {6, 0, 1}, {6, 6, 2}, {4, 6, 3}, {4, 2, 4}, {2, 2, 5}, {2, 6, 6}, {0, 6, 7}},
        {{0, -3, 0},
         {0, 2, 1},
         {1, 8, 2},
         {9, 9, 3},
         {8, 3, 4},
         {8, -2, 5},
         {6, 0, 6},
         {4, 4, 7},
         {2, 1, 8}}};
    for (const auto& fixture : fixtures) {
        for (bool reverse : {false, true}) {
            for (std::size_t first = 0; first < fixture.size(); ++first) {
                auto points = fixture;
                if (reverse)
                    std::reverse(points.begin(), points.end());
                std::rotate(points.begin(), points.begin() + static_cast<std::ptrdiff_t>(first),
                            points.end());
                for (std::size_t vertex = 0; vertex < points.size(); ++vertex)
                    points[vertex].index = vertex + 113;
                const Polygon curve(points);
                std::vector<std::size_t> identity(points.size());
                std::iota(identity.begin(), identity.end(), 0);
                const auto trapezoids =
                    extract_trapezoids(curve, build_full_visibility_map(curve), identity);
                test::check_triangulation(points, triangulate_trapezoids(points, trapezoids));
                if (first == 0)
                    test::check_triangulation(points, triangulate_polygon(points));
            }
        }
    }
    auto points = fixtures.back();
    for (auto& point : points) {
        point.x *= Exact(0x1p-400);
        point.y *= Exact(0x1p400);
    }
    test::check_triangulation(points, triangulate_polygon(points));
    for (std::size_t size : {1U, 2U, 3U, 5U, 17U, 31U}) {
        const auto input = test::alternating_chains(size);
        test::check_triangulation(input.vertices,
                                  triangulate_trapezoids(input.vertices, input.trapezoids), false);
    }
}

void check_premises() {
    const std::vector<Point> points{{0, 8, 0}, {6, 6, 1}, {4, 4, 2}, {1, 2, 3}, {0, 0, 4}};
    const auto valid = one_polygon(points);
    test::require_assertion_abort([&] {
        auto invalid = valid;
        invalid.polygons.clear();
        triangulate_unimonotone(points, invalid);
    });
    test::require_assertion_abort([&] {
        auto invalid = valid;
        invalid.polygons[0].vertices[1] = points.size();
        triangulate_unimonotone(points, invalid);
    });
    test::require_assertion_abort([&] {
        auto invalid = valid;
        invalid.polygons[0].vertices[1] = 0;
        triangulate_unimonotone(points, invalid);
    });
    test::require_assertion_abort([&] {
        auto invalid = valid;
        invalid.polygons[0].top_vertex = 1;
        triangulate_unimonotone(points, invalid);
    });
    test::require_assertion_abort([&] {
        auto invalid = valid;
        std::reverse(invalid.polygons[0].vertices.begin(), invalid.polygons[0].vertices.end());
        triangulate_unimonotone(points, invalid);
    });
    test::require_assertion_abort([&] {
        auto invalid = valid;
        std::swap(invalid.polygons[0].vertices[1], invalid.polygons[0].vertices[2]);
        triangulate_unimonotone(points, invalid);
    });
}

}

int main() {
    check_endpoint_backtracking();
    check_chain_orders();
    check_trapezoid_pipeline();
    check_premises();
    std::puts("[FM84 Algorithm 3 tests]: all passed");
}
