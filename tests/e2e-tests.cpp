#include "support/pipeline_fixtures.h"
#include "support/random.h"
#include "support/triangulation_geometry.h"
#include "support/triangulation_production_checks.h"

#include <algorithm>
#include <cstdio>

using namespace chazelle;

namespace {

void check_geometry(std::span<const Point> points, const Triangulation& result) {
    Exact original_area = 0;
    Exact symbolic_area = 0;
    for (std::size_t vertex = 0; vertex < points.size(); ++vertex) {
        const Point& a = points[vertex];
        const Point& b = points[(vertex + 1) % points.size()];
        original_area += a.x * b.y - a.y * b.x;
        symbolic_area += a.x * test::triangle_height(b) - test::triangle_height(a) * b.x;
    }
    Exact triangle_area = 0;
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    edges.reserve(3 * result.triangles.size());
    for (const Triangle& triangle : result.triangles) {
        const auto& vertices = triangle.vertices;
        const Point& a = points[vertices[0]];
        const Point& b = points[vertices[1]];
        const Point& c = points[vertices[2]];
        const Exact orientation = test::triangle_orientation(a, b, c);
        test::require_triangulation(orientation <= 0,
                                    "[FM84 Algorithm 3]: exact symbolic triangle orientation");
        triangle_area += orientation;
        for (std::size_t side = 0; side < 3; ++side)
            edges.emplace_back(std::min(vertices[side], vertices[(side + 1) % 3]),
                               std::max(vertices[side], vertices[(side + 1) % 3]));
        for (const Point& point : points)
            test::require_triangulation(
                !(test::triangle_orientation(a, b, point) < 0 &&
                  test::triangle_orientation(b, c, point) < 0 &&
                  test::triangle_orientation(c, a, point) < 0),
                "[FM84 tex 448-452]: no original vertex lies strictly inside a triangle");
        if (orientation < 0) {
            const Exact x = (a.x + b.x + c.x) / 3;
            const Exact y =
                (test::triangle_height(a) + test::triangle_height(b) + test::triangle_height(c)) /
                3;
            bool inside = false;
            for (std::size_t vertex = 0; vertex < points.size(); ++vertex) {
                const Point& p = points[vertex];
                const Point& q = points[(vertex + 1) % points.size()];
                const Exact py = test::triangle_height(p);
                const Exact qy = test::triangle_height(q);
                if ((py > y) != (qy > y) && x < p.x + (y - py) * (q.x - p.x) / (qy - py))
                    inside = !inside;
            }
            test::require_triangulation(
                inside, "[C91 Theorem 4.3]: every nondegenerate triangle lies inside the polygon");
        }
    }
    test::require_triangulation(triangle_area ==
                                    (original_area < 0 ? symbolic_area : -symbolic_area),
                                "[C91 Theorem 4.3]: exact symbolic area preservation");
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    for (std::size_t first = 0; first < edges.size(); ++first)
        for (std::size_t second = first + 1; second < edges.size(); ++second)
            test::require_triangulation(!test::proper_triangle_crossing(
                                            points[edges[first].first], points[edges[first].second],
                                            points[edges[second].first],
                                            points[edges[second].second]),
                                        "[C91 tex 26]: output triangle edges do not cross");
}

Triangulation check_polygon(const std::vector<Point>& points, bool geometry) {
    const auto result = triangulate_polygon(points);
    test::check_production_triangulation(points, result);
    if (geometry)
        check_geometry(points, result);
    return result;
}

void check_orders_and_transforms() {
    for (const auto& fixture : test::polygon_fixtures()) {
        for (bool reverse : {false, true}) {
            for (std::size_t first = 0; first < fixture.size(); ++first) {
                auto points = test::boundary_order(fixture, reverse, first, SOS_NONE - 257);
                check_polygon(points, true);
            }
        }
        for (std::size_t transform = 0; transform < 3; ++transform) {
            auto points = fixture;
            for (Point& point : points) {
                if (transform == 0) {
                    std::swap(point.x, point.y);
                } else if (transform == 1) {
                    point.x = Exact(0x1p-400) * point.x + 3;
                    point.y = Exact(0x1p400) * point.y - 7;
                } else {
                    const Exact x = point.x;
                    point.x = x + point.y;
                    point.y -= x;
                }
            }
            check_polygon(test::boundary_order(points, false, 0, 71), true);
        }
    }
}

void check_random_polygons() {
    test::DeterministicRandomGenerator random(431984);
    for (std::size_t sample = 0; sample < 24; ++sample) {
        const std::size_t count = 5 + random.next() % 14;
        std::vector<Point> points{{-1, -3, 0}};
        points.reserve(count);
        for (std::size_t vertex = 1; vertex + 1 < count; ++vertex)
            points.push_back({Exact(vertex) / 3, Exact(random.next() % 11) / 5, vertex});
        points.push_back({Exact(count), -3, count - 1});
        check_polygon(test::boundary_order(points, sample % 2 == 0, random.next() % count), true);
    }
}

void check_grade_boundaries() {
    for (std::size_t count : {5U, 17U, 33U, 65U, 257U, 1026U}) {
        std::vector<Point> points;
        points.reserve(count);
        for (std::size_t vertex = 0; vertex < count; ++vertex)
            points.push_back({Exact(vertex), Exact(vertex) * Exact(vertex), vertex});
        points = test::boundary_order(points, true, count / 3, 101);
        const auto result = check_polygon(points, count <= 33);
        const auto top = std::max_element(points.begin(), points.end(),
                                          [](const Point& a, const Point& b) { return a.y < b.y; });
        const auto index = static_cast<std::size_t>(top - points.begin());
        test::require_triangulation(
            result.vertex_triangles[index].size() == count - 2,
            "[FM84 Algorithm 3]: this convex unimonotone chain produces the expected fan");
    }
}

void check_subdivided_concave_polygon() {
    const auto boundary = test::polygon_fixtures()[3];
    for (std::size_t subdivisions : {3U, 129U}) {
        std::vector<Point> points;
        points.reserve(boundary.size() * subdivisions);
        for (std::size_t edge = 0; edge < boundary.size(); ++edge) {
            const Point& a = boundary[edge];
            const Point& b = boundary[(edge + 1) % boundary.size()];
            for (std::size_t position = 0; position < subdivisions; ++position) {
                const Exact parameter = Exact(position) / Exact(subdivisions);
                points.push_back(
                    {a.x + parameter * (b.x - a.x), a.y + parameter * (b.y - a.y), points.size()});
            }
        }
        points = test::boundary_order(points, true, points.size() / 3, SOS_NONE - 2049);
        const auto result = check_polygon(points, subdivisions == 3);
        test::require_triangulation(
            result.work.polygon_vertices > points.size(),
            "[FM84 §3]: this concave polygon requires class-B diagonals before triangulation");
    }
}

}

int main() {
    check_orders_and_transforms();
    check_random_polygons();
    check_grade_boundaries();
    check_subdivided_concave_polygon();
    std::puts("Full polygon triangulation end-to-end tests passed");
}
