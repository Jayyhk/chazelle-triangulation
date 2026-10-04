#include "support/unimonotone_fixture.h"
#include "triangulation/unimonotone.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace chazelle;

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "%s\n", message);
        std::abort();
    }
}

void check_partition(const test::TrapezoidizedPolygon& input, std::size_t expected_diagonals,
                     bool triangles) {
    const auto& points = input.vertices;
    const std::size_t count = points.size();
    const auto result = decompose_unimonotone(points, input.trapezoids);
    require(result.diagonals.size() == expected_diagonals &&
                result.polygons.size() == expected_diagonals + 1,
            "[FM84 Algorithm 2]: exact number of diagonals and unimonotone polygons");
    require(result.work.vertices_initialized == count &&
                result.work.trapezoids_examined == input.trapezoids.trapezoids.size() &&
                result.work.vertex_visits <= count + 2 * expected_diagonals &&
                result.work.emitted_vertices == count + 2 * expected_diagonals &&
                result.work.maximum_stack_size <= expected_diagonals,
            "[FM84 tex 466-474]: linear work independent of polygon size or split depth");
    std::vector<std::array<std::size_t, 3>> incident(count);
    std::vector<std::size_t> degrees(count, 0);
    std::vector<bool> selected(input.trapezoids.trapezoids.size(), false);
    for (std::size_t index = 0; index < result.diagonals.size(); ++index) {
        const auto& diagonal = result.diagonals[index];
        require(diagonal.trapezoid < selected.size() && !selected[diagonal.trapezoid],
                "[FM84 Algorithm 2]: each class-B trapezoid supplies exactly one diagonal");
        selected[diagonal.trapezoid] = true;
        const auto& trapezoid = input.trapezoids.trapezoids[diagonal.trapezoid];
        require(diagonal.top_vertex == trapezoid.top_vertex &&
                    diagonal.bottom_vertex == trapezoid.bottom_vertex,
                "[FM84 tex 328-332]: a diagonal joins its trapezoid's defining vertices");
        for (const std::size_t vertex : {diagonal.top_vertex, diagonal.bottom_vertex}) {
            require(degrees[vertex] < 3, "[FM84 tex 165-174]: bounded trapezoid incidence");
            incident[vertex][degrees[vertex]++] = index;
        }
    }
    for (std::size_t index = 0; index < selected.size(); ++index) {
        const auto& trapezoid = input.trapezoids.trapezoids[index];
        const bool class_b = (trapezoid.top_vertex + 1) % count != trapezoid.bottom_vertex &&
                             (trapezoid.bottom_vertex + 1) % count != trapezoid.top_vertex;
        require(selected[index] == class_b,
                "[FM84 Algorithm 2]: all class-B diagonals are present");
    }
    Exact input_area = 0;
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        const Point& a = points[vertex];
        const Point& b = points[(vertex + 1) % count];
        input_area += a.x * b.y - a.y * b.x;
    }
    Exact output_area = 0;
    std::vector<std::size_t> boundary_edges(count, 0);
    std::vector<std::array<std::size_t, 2>> diagonal_edges(result.diagonals.size());
    std::vector<bool> seen(count, false);
    for (const auto& polygon : result.polygons) {
        const auto& vertices = polygon.vertices;
        require(vertices.size() >= 3 && (!triangles || vertices.size() == 3),
                "[FM84 Algorithm 2]: valid piece size for this polygon");
        Exact area = 0;
        std::size_t top = 0;
        std::size_t bottom = 0;
        for (std::size_t position = 0; position < vertices.size(); ++position) {
            const std::size_t a = vertices[position];
            const std::size_t b = vertices[(position + 1) % vertices.size()];
            require(a < count && b < count, "[FM84 Algorithm 2]: original vertex indices");
            seen[a] = true;
            area += points[a].x * points[b].y - points[a].y * points[b].x;
            if (point_y_above(points[a], points[vertices[top]]))
                top = position;
            if (point_y_below(points[a], points[vertices[bottom]]))
                bottom = position;
            if ((a + 1) % count == b || (b + 1) % count == a) {
                const std::size_t edge = (a + 1) % count == b ? a : b;
                require(a == (input_area < 0 ? edge : (edge + 1) % count),
                        "[FM84 Algorithm 1 input]: every piece is clockwise");
                ++boundary_edges[edge];
            } else {
                bool found = false;
                for (std::size_t slot = 0; slot < degrees[a]; ++slot) {
                    const std::size_t index = incident[a][slot];
                    const auto& diagonal = result.diagonals[index];
                    const std::size_t target =
                        a == diagonal.top_vertex ? diagonal.bottom_vertex : diagonal.top_vertex;
                    if (target == b) {
                        require(!found, "[FM84 Algorithm 2]: unique diagonal endpoints");
                        ++diagonal_edges[index][a == diagonal.top_vertex ? 0 : 1];
                        found = true;
                    }
                }
                require(found, "[FM84 Algorithm 2]: no edge outside the prescribed partition");
            }
        }
        require(area <= 0, "[FM84 Algorithm 1 input]: clockwise or zero-area symbolic piece");
        output_area += area;
        require(
            polygon.top_vertex == vertices[top] && polygon.bottom_vertex == vertices[bottom] &&
                ((top + 1) % vertices.size() == bottom || (bottom + 1) % vertices.size() == top),
            "[FM84 tex 334-339]: adjacent y extrema");
        const bool increasing = (top + 1) % vertices.size() == bottom;
        const std::size_t end = increasing ? top : bottom;
        std::size_t position = increasing ? bottom : top;
        while (position != end) {
            const std::size_t next = (position + 1) % vertices.size();
            require(increasing ? point_y_above(points[vertices[next]], points[vertices[position]])
                               : point_y_below(points[vertices[next]], points[vertices[position]]),
                    "[FM84 tex 334-339]: a monotone remaining chain");
            position = next;
        }
    }
    require(output_area == (input_area < 0 ? input_area : -input_area),
            "[FM84 Algorithm 2]: exact area preservation");
    require(std::all_of(seen.begin(), seen.end(), [](bool present) { return present; }) &&
                std::all_of(boundary_edges.begin(), boundary_edges.end(),
                            [](std::size_t occurrences) { return occurrences == 1; }) &&
                std::all_of(diagonal_edges.begin(), diagonal_edges.end(),
                            [](const auto& occurrences) {
                                return occurrences[0] == 1 && occurrences[1] == 1;
                            }),
            "[FM84 Algorithm 2]: original boundary once, each diagonal in both directions");
}

void change_boundary_order(test::TrapezoidizedPolygon& polygon, bool reverse, std::size_t first) {
    const std::size_t count = polygon.vertices.size();
    if (reverse)
        std::reverse(polygon.vertices.begin(), polygon.vertices.end());
    std::rotate(polygon.vertices.begin(),
                polygon.vertices.begin() + static_cast<std::ptrdiff_t>(first),
                polygon.vertices.end());
    std::vector<std::size_t> positions(count);
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        positions[polygon.vertices[vertex].index] = vertex;
        polygon.vertices[vertex].index = vertex + 71;
    }
    polygon.trapezoids.vertex_trapezoids.assign(count, {});
    for (std::size_t index = 0; index < polygon.trapezoids.trapezoids.size(); ++index) {
        auto& trapezoid = polygon.trapezoids.trapezoids[index];
        trapezoid.top_vertex = positions[trapezoid.top_vertex];
        trapezoid.bottom_vertex = positions[trapezoid.bottom_vertex];
        trapezoid.left_edge = positions[(trapezoid.left_edge + (reverse ? 1 : 0)) % count];
        trapezoid.right_edge = positions[(trapezoid.right_edge + (reverse ? 1 : 0)) % count];
        auto& entry = polygon.trapezoids.vertex_trapezoids[trapezoid.top_vertex];
        require(entry.count == 0, "analytic fixture has at most one trapezoid per vertex");
        entry = {{{index, NONE}}, 1};
    }
}

void check_large_alternating_chains() {
    for (bool reverse : {false, true}) {
        auto polygon = test::alternating_chains(50000);
        change_boundary_order(polygon, reverse, reverse ? 431 : 0);
        check_partition(polygon, polygon.vertices.size() - 3, true);
    }
}

void check_large_unimonotone_polygon() {
    test::TrapezoidizedPolygon polygon;
    constexpr std::size_t count = 200003;
    polygon.vertices.reserve(count);
    polygon.trapezoids.trapezoids.reserve(count - 1);
    polygon.trapezoids.vertex_trapezoids.resize(count);
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        polygon.vertices.push_back({Exact(vertex), Exact(vertex) * Exact(vertex), vertex});
        if (vertex > 0) {
            polygon.trapezoids.trapezoids.push_back({vertex, vertex - 1, count - 1, vertex - 1});
            polygon.trapezoids.vertex_trapezoids[vertex] = {{{vertex - 1, NONE}}, 1};
        }
    }
    check_partition(polygon, 0, false);
}

}

int main() {
    check_large_alternating_chains();
    check_large_unimonotone_polygon();
    std::puts("[FM84 Algorithm 2 production tests]: all passed");
}
