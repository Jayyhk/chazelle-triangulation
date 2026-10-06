#include "algorithm/triangulation/unimonotone.h"
#include "algorithm/visibility/naive_visibility.h"
#include "support/assertions.h"
#include "support/unimonotone_fixture.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <utility>
#include <vector>

using namespace chazelle;

namespace {

Exact height(const Point& point) {
    return point.y - Exact(point.index) * Exact::infinitesimal(3);
}

Exact orientation(const Point& a, const Point& b, const Point& c) {
    return (b.x - a.x) * (height(c) - height(a)) - (height(b) - height(a)) * (c.x - a.x);
}

bool on_segment(const Point& a, const Point& b, const Point& point) {
    return orientation(a, b, point) == 0 && std::min(a.x, b.x) <= point.x &&
           point.x <= std::max(a.x, b.x) && std::min(height(a), height(b)) <= height(point) &&
           height(point) <= std::max(height(a), height(b));
}

bool intersects(const Point& a, const Point& b, const Point& c, const Point& d) {
    const Exact abc = orientation(a, b, c);
    const Exact abd = orientation(a, b, d);
    const Exact cda = orientation(c, d, a);
    const Exact cdb = orientation(c, d, b);
    return (((abc > 0 && abd < 0) || (abc < 0 && abd > 0)) &&
            ((cda > 0 && cdb < 0) || (cda < 0 && cdb > 0))) ||
           on_segment(a, b, c) || on_segment(a, b, d) || on_segment(c, d, a) || on_segment(c, d, b);
}

bool midpoint_inside(std::span<const Point> vertices, std::size_t a, std::size_t b) {
    const Exact x = (vertices[a].x + vertices[b].x) / 2;
    const Exact y = (height(vertices[a]) + height(vertices[b])) / 2;
    bool inside = false;
    for (std::size_t edge = 0; edge < vertices.size(); ++edge) {
        const Point& first = vertices[edge];
        const Point& last = vertices[(edge + 1) % vertices.size()];
        const Exact first_y = height(first);
        const Exact last_y = height(last);
        if ((first_y > y) != (last_y > y)) {
            const Exact crossing =
                first.x + (y - first_y) * (last.x - first.x) / (last_y - first_y);
            if (x < crossing)
                inside = !inside;
        }
    }
    return inside;
}

TrapezoidDecomposition extract(std::span<const Point> vertices) {
    const Polygon curve(std::vector<Point>(vertices.begin(), vertices.end()));
    std::vector<std::size_t> identity(vertices.size());
    std::iota(identity.begin(), identity.end(), 0);
    return extract_trapezoids(curve, build_full_visibility_map(curve), identity);
}

void check_partition(std::span<const Point> points, const TrapezoidDecomposition& trapezoids,
                     const UnimonotoneDecomposition& result) {
    const std::size_t count = points.size();
    assert(count >= 3);
    Exact original_area = 0;
    Exact original_symbolic_area = 0;
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        const Point& a = points[vertex];
        const Point& b = points[(vertex + 1) % count];
        original_area += a.x * b.y - b.x * a.y;
        original_symbolic_area += a.x * height(b) - b.x * height(a);
    }
    std::vector<bool> chosen(trapezoids.trapezoids.size(), false);
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> expected_edges;
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        const std::size_t next = (vertex + 1) % count;
        ++expected_edges[original_area < 0 ? std::pair{vertex, next} : std::pair{next, vertex}];
    }
    for (const TrapezoidDiagonal& diagonal : result.diagonals) {
        assert(diagonal.trapezoid < trapezoids.trapezoids.size());
        assert(!chosen[diagonal.trapezoid]);
        chosen[diagonal.trapezoid] = true;
        const Trapezoid& trapezoid = trapezoids.trapezoids[diagonal.trapezoid];
        const std::size_t top = diagonal.top_vertex;
        const std::size_t bottom = diagonal.bottom_vertex;
        assert(top == trapezoid.top_vertex && bottom == trapezoid.bottom_vertex);
        assert((top + 1) % count != bottom && (bottom + 1) % count != top);
        ++expected_edges[{top, bottom}];
        ++expected_edges[{bottom, top}];
        assert(midpoint_inside(points, top, bottom) &&
               "[FM84 tex 328-332]: a class-B diagonal lies inside the polygon");
        for (std::size_t edge = 0; edge < count; ++edge) {
            const std::size_t next = (edge + 1) % count;
            if (edge == top || edge == bottom || next == top || next == bottom)
                continue;
            assert(!intersects(points[top], points[bottom], points[edge], points[next]) &&
                   "[FM84 tex 328-332]: a diagonal does not cross the boundary");
        }
    }
    for (std::size_t index = 0; index < trapezoids.trapezoids.size(); ++index) {
        const Trapezoid& trapezoid = trapezoids.trapezoids[index];
        const bool class_b = (trapezoid.top_vertex + 1) % count != trapezoid.bottom_vertex &&
                             (trapezoid.bottom_vertex + 1) % count != trapezoid.top_vertex;
        assert(chosen[index] == class_b &&
               "[FM84 tex 328-339]: all and only class-B trapezoids supply diagonals");
    }
    for (std::size_t i = 0; i < result.diagonals.size(); ++i) {
        const auto& a = result.diagonals[i];
        for (std::size_t j = 0; j < i; ++j) {
            const auto& b = result.diagonals[j];
            if (a.top_vertex == b.top_vertex || a.top_vertex == b.bottom_vertex ||
                a.bottom_vertex == b.top_vertex || a.bottom_vertex == b.bottom_vertex)
                continue;
            assert(!intersects(points[a.top_vertex], points[a.bottom_vertex], points[b.top_vertex],
                               points[b.bottom_vertex]) &&
                   "[FM84 tex 328-339]: class-B diagonals do not cross one another");
        }
    }
    Exact area = 0;
    Exact symbolic_area = 0;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> actual_edges;
    std::size_t emitted = 0;
    for (const auto& polygon : result.polygons) {
        const auto& boundary = polygon.vertices;
        assert(boundary.size() >= 3 && boundary.size() <= count);
        assert(std::set<std::size_t>(boundary.begin(), boundary.end()).size() == boundary.size());
        emitted += boundary.size();
        Exact piece_area = 0;
        Exact piece_symbolic_area = 0;
        std::size_t top = 0;
        std::size_t bottom = 0;
        for (std::size_t position = 0; position < boundary.size(); ++position) {
            const std::size_t next = (position + 1) % boundary.size();
            assert(boundary[position] < count);
            const Point& a = points[boundary[position]];
            const Point& b = points[boundary[next]];
            piece_area += a.x * b.y - b.x * a.y;
            piece_symbolic_area += a.x * height(b) - b.x * height(a);
            ++actual_edges[{boundary[position], boundary[next]}];
            if (height(a) > height(points[boundary[top]]))
                top = position;
            if (height(a) < height(points[boundary[bottom]]))
                bottom = position;
        }
        assert(piece_area <= 0 && piece_symbolic_area < 0);
        assert(polygon.top_vertex == boundary[top] && polygon.bottom_vertex == boundary[bottom]);
        assert((top + 1) % boundary.size() == bottom || (bottom + 1) % boundary.size() == top);
        const bool increasing = (top + 1) % boundary.size() == bottom;
        const std::size_t end = increasing ? top : bottom;
        std::size_t position = increasing ? bottom : top;
        while (position != end) {
            const std::size_t next = (position + 1) % boundary.size();
            assert(increasing
                       ? height(points[boundary[position]]) < height(points[boundary[next]])
                       : height(points[boundary[position]]) > height(points[boundary[next]]));
            position = next;
        }
        area += piece_area;
        symbolic_area += piece_symbolic_area;
    }
    assert(area == (original_area < 0 ? original_area : -original_area));
    assert(symbolic_area == (original_area < 0 ? original_symbolic_area : -original_symbolic_area));
    assert(actual_edges == expected_edges);
    const std::size_t diagonals = result.diagonals.size();
    assert(diagonals <= count - 3 && result.polygons.size() == diagonals + 1);
    assert(emitted == count + 2 * diagonals && result.work.emitted_vertices == emitted);
    assert(result.work.vertices_initialized == count);
    assert(result.work.vertex_visits <= count + 2 * diagonals);
    assert(result.work.trapezoids_examined == trapezoids.trapezoids.size());
    assert(result.work.maximum_stack_size <= diagonals);
}

std::vector<Point> two_trapezoids_at_one_vertex() {
    return {{0, -3, 0}, {0, 2, 1}, {1, 8, 2}, {9, 9, 3}, {8, 3, 4},
            {8, -2, 5}, {6, 0, 6}, {4, 4, 7}, {2, 1, 8}};
}

void check_both_trapezoids() {
    const auto points = two_trapezoids_at_one_vertex();
    const auto trapezoids = extract(points);
    const auto& entry = trapezoids.vertex_trapezoids[7];
    assert(entry.count == 2);
    assert(trapezoids.trapezoids[entry.trapezoids[0]].bottom_vertex == 1);
    assert(trapezoids.trapezoids[entry.trapezoids[1]].bottom_vertex == 4);
    const auto result = decompose_unimonotone(points, trapezoids);
    check_partition(points, trapezoids, result);
    std::vector<std::size_t> bottoms;
    for (const auto& diagonal : result.diagonals)
        if (diagonal.top_vertex == 7)
            bottoms.push_back(diagonal.bottom_vertex);
    assert((bottoms == std::vector<std::size_t>{4, 1}) &&
           "[FM84 Algorithm 2]: clockwise splitting consumes the rightmost trapezoid first");
}

void check_orders_and_degeneracies() {
    const std::vector<std::vector<Point>> fixtures{
        {{0, 0, 0}, {4, 0, 1}, {2, 3, 2}},
        {{0, 0, 0}, {2, 0, 1}, {4, 0, 2}, {4, 2, 3}, {4, 4, 4}, {0, 4, 5}},
        {{0, 0, 0}, {6, 0, 1}, {6, 6, 2}, {4, 6, 3}, {4, 2, 4}, {2, 2, 5}, {2, 6, 6}, {0, 6, 7}},
        two_trapezoids_at_one_vertex()};
    for (const auto& fixture : fixtures) {
        for (bool reverse : {false, true}) {
            for (std::size_t first = 0; first < fixture.size(); ++first) {
                auto points = fixture;
                if (reverse)
                    std::reverse(points.begin(), points.end());
                std::rotate(points.begin(), points.begin() + static_cast<std::ptrdiff_t>(first),
                            points.end());
                for (std::size_t vertex = 0; vertex < points.size(); ++vertex)
                    points[vertex].index = SOS_NONE - 127 + vertex;
                const auto trapezoids = extract(points);
                check_partition(points, trapezoids, decompose_unimonotone(points, trapezoids));
            }
        }
    }
    for (auto points : fixtures) {
        for (auto& point : points) {
            point.x *= Exact(0x1p-400);
            point.y *= Exact(0x1p400);
        }
        const auto trapezoids = extract(points);
        check_partition(points, trapezoids, decompose_unimonotone(points, trapezoids));
    }
}

void check_random_polygons() {
    std::mt19937 random(9837);
    for (std::size_t test = 0; test < 200; ++test) {
        const std::size_t count = 3 + random() % 30;
        std::vector<Point> points{{-1, -3, 0}};
        for (std::size_t vertex = 1; vertex + 1 < count; ++vertex)
            points.push_back({Exact(vertex) / 7, Exact(random() % 10) / 3, vertex});
        points.push_back({Exact(count), -3, count - 1});
        if (test % 2 == 0)
            std::reverse(points.begin(), points.end());
        std::rotate(points.begin(), points.begin() + static_cast<std::ptrdiff_t>(random() % count),
                    points.end());
        for (std::size_t vertex = 0; vertex < count; ++vertex) {
            if (test % 3 == 0) {
                const Exact x = points[vertex].x;
                points[vertex].x = -points[vertex].y;
                points[vertex].y = x;
            }
            points[vertex].index = vertex + 19;
        }
        const auto trapezoids = extract(points);
        check_partition(points, trapezoids, decompose_unimonotone(points, trapezoids));
    }
}

void check_analytic_fixture() {
    for (std::size_t size : {1U, 2U, 3U, 5U, 17U, 31U}) {
        const auto polygon = test::alternating_chains(size);
        const auto expected = extract(polygon.vertices);
        assert(expected.trapezoids.size() == polygon.trapezoids.trapezoids.size());
        std::set<std::vector<std::size_t>> actual_keys;
        std::set<std::vector<std::size_t>> expected_keys;
        for (const auto& t : polygon.trapezoids.trapezoids)
            actual_keys.insert({t.top_vertex, t.bottom_vertex, t.left_edge, t.right_edge});
        for (const auto& t : expected.trapezoids)
            expected_keys.insert({t.top_vertex, t.bottom_vertex, t.left_edge, t.right_edge});
        assert(actual_keys == expected_keys);
        const auto result = decompose_unimonotone(polygon.vertices, polygon.trapezoids);
        check_partition(polygon.vertices, polygon.trapezoids, result);
        assert(result.diagonals.size() == polygon.vertices.size() - 3);
        assert(std::all_of(result.polygons.begin(), result.polygons.end(),
                           [](const auto& piece) { return piece.vertices.size() == 3; }));
    }
}

void check_pipeline() {
    const auto points = two_trapezoids_at_one_vertex();
    const auto trapezoids = compute_trapezoid_decomposition(points);
    check_partition(points, trapezoids, compute_unimonotone_decomposition(points));
}

void check_premises() {
    const auto points = two_trapezoids_at_one_vertex();
    const auto valid = extract(points);
    test::require_assertion_abort([&] {
        auto invalid = valid;
        invalid.vertex_trapezoids.pop_back();
        decompose_unimonotone(points, invalid);
    });
    test::require_assertion_abort([&] {
        auto invalid = valid;
        invalid.vertex_trapezoids[7].count = 3;
        decompose_unimonotone(points, invalid);
    });
    test::require_assertion_abort([&] {
        auto invalid = valid;
        invalid.vertex_trapezoids[7].trapezoids[1] = invalid.vertex_trapezoids[7].trapezoids[0];
        decompose_unimonotone(points, invalid);
    });
    test::require_assertion_abort([&] {
        auto invalid = valid;
        invalid.trapezoids[0].bottom_vertex = points.size();
        decompose_unimonotone(points, invalid);
    });
    test::require_assertion_abort([&] {
        auto invalid = valid;
        invalid.vertex_trapezoids[7].count = 0;
        decompose_unimonotone(points, invalid);
    });
}

}

int main() {
    check_both_trapezoids();
    check_orders_and_degeneracies();
    check_random_polygons();
    check_analytic_fixture();
    check_pipeline();
    check_premises();
    std::puts("[FM84 Algorithm 2 tests]: all passed");
}
