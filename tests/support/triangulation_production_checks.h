#pragma once

#include "triangulation/triangulation.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <tuple>

namespace chazelle::test {

inline void require_triangulation(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "%s\n", message);
        std::abort();
    }
}

inline void check_production_triangulation(std::span<const Point> points,
                                           const Triangulation& result) {
    const std::size_t count = points.size();
    require_triangulation(count >= 3 && result.triangles.size() == count - 2 &&
                              result.vertex_triangles.size() == count,
                          "[FM84 Theorem 3]: n-2 triangles, with original vertex adjacency");
    Exact area = 0;
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        const Point& a = points[vertex];
        const Point& b = points[(vertex + 1) % count];
        area += a.x * b.y - a.y * b.x;
    }
    require_triangulation(area != 0, "[FM84 input]: the simple polygon encloses positive area");
    std::vector<std::tuple<std::size_t, std::size_t, bool>> edges;
    edges.reserve(3 * result.triangles.size());
    Exact triangle_area = 0;
    for (const Triangle& triangle : result.triangles) {
        const auto& vertices = triangle.vertices;
        require_triangulation(vertices[0] < count && vertices[1] < count && vertices[2] < count &&
                                  vertices[0] != vertices[1] && vertices[1] != vertices[2] &&
                                  vertices[2] != vertices[0],
                              "[FM84 Algorithm 3]: distinct original triangle corners");
        const Point& a = points[vertices[0]];
        const Point& b = points[vertices[1]];
        const Point& c = points[vertices[2]];
        const Exact turn = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        require_triangulation(turn <= 0, "[FM84 Algorithm 3]: clockwise triangles");
        if (turn == 0) {
            const Exact coefficient = (b.x - a.x) * (Exact(a.index) - Exact(c.index)) -
                                      (Exact(a.index) - Exact(b.index)) * (c.x - a.x);
            require_triangulation(coefficient <= 0,
                                  "[FM84 Algorithm 3]: exact convexity at tied coordinates");
        }
        triangle_area += turn;
        for (std::size_t side = 0; side < 3; ++side) {
            const std::size_t first = vertices[side];
            const std::size_t last = vertices[(side + 1) % 3];
            edges.emplace_back(std::min(first, last), std::max(first, last), first < last);
        }
    }
    require_triangulation(triangle_area == (area < 0 ? area : -area),
                          "[FM84 Algorithm 3]: exact original area preservation");
    std::sort(edges.begin(), edges.end());
    std::size_t boundary = 0;
    std::size_t internal = 0;
    for (std::size_t first = 0; first < edges.size();) {
        const auto [a, b, direction] = edges[first];
        std::size_t last = first + 1;
        while (last < edges.size() && std::get<0>(edges[last]) == a &&
               std::get<1>(edges[last]) == b)
            ++last;
        if (b == a + 1 || (a == 0 && b == count - 1)) {
            require_triangulation(last == first + 1 && direction == ((area < 0) == (b == a + 1)),
                                  "[FM84 Algorithm 3]: each original boundary edge occurs once");
            ++boundary;
        } else {
            require_triangulation(
                last == first + 2 && direction != std::get<2>(edges[first + 1]),
                "[FM84 Algorithm 3]: each diagonal occurs in opposite directions");
            ++internal;
        }
        first = last;
    }
    require_triangulation(boundary == count && internal == count - 3,
                          "[FM84 Algorithm 3]: complete boundary and n-3 internal diagonals");
    std::vector<std::size_t> incidents(result.triangles.size(), 0);
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        require_triangulation(
            !result.vertex_triangles[vertex].empty(),
            "[FM84 Algorithm 3 output]: every original vertex has triangle adjacency");
        std::size_t previous = NONE;
        for (const std::size_t triangle : result.vertex_triangles[vertex]) {
            require_triangulation(triangle < result.triangles.size() &&
                                      (previous == NONE || previous < triangle),
                                  "[FM84 Algorithm 3 output]: valid, distinct incident triangles");
            const auto& corners = result.triangles[triangle].vertices;
            require_triangulation(std::count(corners.begin(), corners.end(), vertex) == 1,
                                  "[FM84 Algorithm 3 output]: adjacency refers to this vertex");
            ++incidents[triangle];
            previous = triangle;
        }
    }
    require_triangulation(
        std::all_of(incidents.begin(), incidents.end(),
                    [](std::size_t value) { return value == 3; }),
        "[FM84 Algorithm 3 output]: exactly three adjacency entries per triangle");
    require_triangulation(result.work.polygon_vertices <= 3 * count - 6 &&
                              result.work.removed_vertices == count - 2 &&
                              result.work.forward_steps == result.work.backward_steps &&
                              result.work.backward_steps <= count - 3 &&
                              result.work.convexity_tests ==
                                  count - 2 + result.work.forward_steps &&
                              result.work.convexity_tests <= 2 * (count - 2) &&
                              result.work.triangle_incidents == 3 * (count - 2),
                          "[FM84 Theorems 2-3]: linear work, including triangle adjacency");
}

}
