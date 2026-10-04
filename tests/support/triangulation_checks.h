#pragma once

#include "triangulation/triangulation.h"
#include "triangulation_geometry.h"

#include <algorithm>
#include <cassert>
#include <map>
#include <utility>

namespace chazelle::test {

inline void check_triangulation(std::span<const Point> points, const Triangulation& result,
                                bool check_geometry = true) {
    const std::size_t count = points.size();
    assert(count >= 3);
    assert(result.triangles.size() == count - 2 && result.vertex_triangles.size() == count);
    Exact area = 0;
    Exact symbolic_area = 0;
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        const Point& a = points[vertex];
        const Point& b = points[(vertex + 1) % count];
        area += a.x * b.y - a.y * b.x;
        symbolic_area += a.x * triangle_height(b) - triangle_height(a) * b.x;
    }
    assert(area != 0);
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edges;
    std::vector<std::size_t> incidents(result.triangles.size(), 0);
    Exact triangle_area = 0;
    Exact triangle_symbolic_area = 0;
    for (const Triangle& triangle : result.triangles) {
        const auto& vertices = triangle.vertices;
        assert(vertices[0] < count && vertices[1] < count && vertices[2] < count);
        assert(vertices[0] != vertices[1] && vertices[1] != vertices[2] &&
               vertices[2] != vertices[0]);
        const Point& a = points[vertices[0]];
        const Point& b = points[vertices[1]];
        const Point& c = points[vertices[2]];
        const Exact raw = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        const Exact symbolic = triangle_orientation(a, b, c);
        assert(raw <= 0 && symbolic <= 0);
        triangle_area += raw;
        triangle_symbolic_area += symbolic;
        for (std::size_t side = 0; side < 3; ++side)
            ++edges[{vertices[side], vertices[(side + 1) % 3]}];
        if (check_geometry) {
            for (std::size_t vertex = 0; vertex < count; ++vertex) {
                if (std::find(vertices.begin(), vertices.end(), vertex) != vertices.end())
                    continue;
                const Point& point = points[vertex];
                assert(!(triangle_orientation(a, b, point) < 0 &&
                         triangle_orientation(b, c, point) < 0 &&
                         triangle_orientation(c, a, point) < 0));
            }
            if (symbolic < 0) {
                const Exact x = (a.x + b.x + c.x) / 3;
                const Exact y = (triangle_height(a) + triangle_height(b) + triangle_height(c)) / 3;
                bool inside = false;
                for (std::size_t vertex = 0; vertex < count; ++vertex) {
                    const Point& p = points[vertex];
                    const Point& q = points[(vertex + 1) % count];
                    const Exact py = triangle_height(p);
                    const Exact qy = triangle_height(q);
                    if ((py > y) != (qy > y) && x < p.x + (y - py) * (q.x - p.x) / (qy - py))
                        inside = !inside;
                }
                assert(inside);
            }
        }
    }
    assert(triangle_area == (area < 0 ? area : -area));
    assert(triangle_symbolic_area == (area < 0 ? symbolic_area : -symbolic_area));
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        std::size_t previous = NONE;
        assert(!result.vertex_triangles[vertex].empty());
        for (const std::size_t triangle : result.vertex_triangles[vertex]) {
            assert(triangle < result.triangles.size() && (previous == NONE || previous < triangle));
            const auto& corners = result.triangles[triangle].vertices;
            assert(std::count(corners.begin(), corners.end(), vertex) == 1);
            ++incidents[triangle];
            previous = triangle;
        }
        const std::size_t next = (vertex + 1) % count;
        const auto edge = area < 0 ? std::pair{vertex, next} : std::pair{next, vertex};
        assert((edges[edge] == 1 && edges[std::pair{edge.second, edge.first}] == 0));
        edges.erase(edge);
        edges.erase({edge.second, edge.first});
    }
    assert(std::all_of(incidents.begin(), incidents.end(),
                       [](std::size_t value) { return value == 3; }));
    assert(edges.size() == 2 * (count - 3));
    for (const auto& [edge, occurrences] : edges) {
        assert(occurrences == 1 && edges.at({edge.second, edge.first}) == 1);
        if (check_geometry) {
            for (const auto& [other, unused] : edges) {
                (void)unused;
                assert(!proper_triangle_crossing(points[edge.first], points[edge.second],
                                                 points[other.first], points[other.second]));
            }
            for (std::size_t vertex = 0; vertex < count; ++vertex)
                assert(!proper_triangle_crossing(points[edge.first], points[edge.second],
                                                 points[vertex], points[(vertex + 1) % count]));
        }
    }
    assert(result.work.polygon_vertices <= 3 * count - 6);
    assert(result.work.removed_vertices == result.triangles.size());
    assert(result.work.forward_steps == result.work.backward_steps);
    assert(result.work.backward_steps <= count - 3);
    assert(result.work.convexity_tests == count - 2 + result.work.forward_steps);
    assert(result.work.convexity_tests <= 2 * (count - 2));
    assert(result.work.triangle_incidents == 3 * (count - 2));
}

}
