#include "algorithm/triangulation/trapezoids.h"
#include "support/triangulation_production_checks.h"

#include <algorithm>
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

void check_convex_polygon(std::size_t count, bool reverse, std::size_t first) {
    std::vector<Point> vertices;
    vertices.reserve(count);
    for (std::size_t vertex = 0; vertex < count; ++vertex)
        vertices.push_back({Exact(vertex), Exact(vertex) * Exact(vertex), vertex});
    if (reverse)
        std::reverse(vertices.begin(), vertices.end());
    std::rotate(vertices.begin(), vertices.begin() + static_cast<std::ptrdiff_t>(first),
                vertices.end());
    std::vector<std::size_t> original_positions(count);
    std::vector<std::size_t> original_vertices(count);
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        original_positions[vertices[vertex].index] = vertex;
        original_vertices[vertex] = vertices[vertex].index;
        vertices[vertex].index = vertex + 17;
    }
    const auto result = compute_trapezoid_decomposition(vertices);
    require(result.trapezoids.size() == count - 1,
            "[FM84 Theorem 5]: this convex polygon has one trapezoid between consecutive heights");
    require(result.vertex_trapezoids.size() == count,
            "[FM84 Algorithm 1 output]: preserve the original vertex table");
    std::vector<bool> seen(count, false);
    for (std::size_t index = 0; index < result.trapezoids.size(); ++index) {
        const Trapezoid& trapezoid = result.trapezoids[index];
        require(trapezoid.top_vertex < count && trapezoid.bottom_vertex < count &&
                    trapezoid.left_edge < count && trapezoid.right_edge < count,
                "[FM84 Algorithm 1 output]: no padding index escapes");
        const std::size_t original_top = original_vertices[trapezoid.top_vertex];
        require(original_top > 0 && !seen[original_top],
                "[FM84 Algorithm 1 output]: each trapezoid appears once at its top vertex");
        seen[original_top] = true;
        require(trapezoid.bottom_vertex == original_positions[original_top - 1] &&
                    trapezoid.left_edge == original_positions[reverse ? 0 : count - 1] &&
                    trapezoid.right_edge ==
                        original_positions[reverse ? original_top : original_top - 1],
                "[FM84 §2]: exact convex-polygon trapezoid boundaries");
        const auto& entry = result.vertex_trapezoids[trapezoid.top_vertex];
        require(entry.count == 1 && entry.trapezoids[0] == index,
                "[FM84 Algorithm 1 output]: correct vertex-to-trapezoid references");
    }
    require(result.vertex_trapezoids[original_positions[0]].count == 0,
            "[FM84 Algorithm 1 output]: no trapezoid begins at the bottommost vertex");
    require(result.work.boundary_edges < 2 * count && result.work.regions <= 4 * count &&
                result.work.arcs <= 16 * count && result.work.arc_edges <= 48 * count &&
                result.work.chords <= 4 * count && result.work.joined_regions <= 4 * count,
            "[C91 §2.1, FM84 Theorem 5]: linear trapezoid extraction work");
    test::check_production_triangulation(vertices, triangulate_trapezoids(vertices, result));
}

}

int main() {
    check_convex_polygon(257, false, 0);
    check_convex_polygon(1024, true, 413);
    check_convex_polygon(2049, false, 0);
    check_convex_polygon(2050, true, 837);
    std::puts("Production trapezoid extraction tests passed");
}
