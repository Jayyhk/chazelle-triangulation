#include "triangulation.h"

#include <cassert>

namespace chazelle {

namespace {

struct VertexLinks {
    std::size_t previous;
    std::size_t next;
};

bool convex(const Point& a, const Point& b, const Point& c) {
    const Exact ab_x = b.x - a.x;
    const Exact ac_x = c.x - a.x;
    const Exact turn = ab_x * (c.y - a.y) - (b.y - a.y) * ac_x;
    if (turn != 0)
        return turn < 0;
    const Exact height_term =
        ab_x * (Exact(a.index) - Exact(c.index)) - (Exact(a.index) - Exact(b.index)) * ac_x;
    return height_term <= 0;
}

void append_triangle(const Triangle& triangle, Triangulation& result) {
    const std::size_t index = result.triangles.size();
    result.triangles.push_back(triangle);
    for (const std::size_t vertex : triangle.vertices) {
        result.vertex_triangles[vertex].push_back(index);
        ++result.work.triangle_incidents;
    }
}

}

Triangulation triangulate_unimonotone(std::span<const Point> points,
                                      const UnimonotoneDecomposition& decomposition) {
    assert(points.size() >= 3 && !decomposition.polygons.empty() &&
           "[FM84 Algorithm 3 input]: a simple polygon's unimonotone partition is required");
    Triangulation result;
    result.triangles.reserve(points.size() - 2);
    result.vertex_triangles.resize(points.size());
    std::vector<std::size_t> seen(points.size(), NONE);
    std::vector<VertexLinks> links;
    for (std::size_t piece = 0; piece < decomposition.polygons.size(); ++piece) {
        const UnimonotonePolygon& polygon = decomposition.polygons[piece];
        const auto& vertices = polygon.vertices;
        const std::size_t count = vertices.size();
        assert(count >= 3 && count <= points.size() &&
               "[FM84 Algorithm 3 input]: each piece is a simple polygon");
        links.resize(count);
        std::size_t top = 0;
        std::size_t bottom = 0;
        Exact twice_area = 0;
        Exact height_term = 0;
        for (std::size_t position = 0; position < count; ++position) {
            const std::size_t vertex = vertices[position];
            const std::size_t next = (position + 1) % count;
            assert(vertex < points.size() && vertices[next] < points.size() &&
                   seen[vertex] != piece &&
                   "[FM84 Algorithm 3 input]: each piece has distinct original vertex positions");
            seen[vertex] = piece;
            const Point& a = points[vertex];
            const Point& b = points[vertices[next]];
            assert(a.index != SOS_NONE && "[C91 §2 tex 47]: symbolic height requires a vertex tag");
            twice_area += a.x * b.y - a.y * b.x;
            height_term += Exact(a.index) * b.x - a.x * Exact(b.index);
            links[position] = {position == 0 ? count - 1 : position - 1, next};
            if (point_y_above(a, points[vertices[top]]))
                top = position;
            if (point_y_below(a, points[vertices[bottom]]))
                bottom = position;
        }
        assert((twice_area < 0 || (twice_area == 0 && height_term < 0)) &&
               "[FM84 tex 182-186]: each piece is clockwise in symbolic coordinates");
        assert(polygon.top_vertex == vertices[top] && polygon.bottom_vertex == vertices[bottom] &&
               (links[top].next == bottom || links[bottom].next == top) &&
               "[FM84 tex 334-337]: the height extrema are adjacent");
        const std::size_t start = links[top].next == bottom ? bottom : top;
        [[maybe_unused]] const std::size_t end = links[start].previous;
        for (std::size_t position = start; position != end; position = links[position].next) {
            [[maybe_unused]] const Point& a = points[vertices[position]];
            [[maybe_unused]] const Point& b = points[vertices[links[position].next]];
            assert((start == top ? point_y_below(b, a) : point_y_above(b, a)) &&
                   "[FM84 tex 334-337]: the other chain is monotone in symbolic height");
        }
        result.work.polygon_vertices += count;
        [[maybe_unused]] const std::size_t first_triangle = result.triangles.size();
        std::size_t remaining = count;
        std::size_t current = links[start].next;
        while (remaining >= 3) {
            assert(current != start && current != end &&
                   "[FM84 Algorithm 3 analysis tex 453-459]: preserve both extrema");
            const std::size_t previous = links[current].previous;
            const std::size_t next = links[current].next;
            ++result.work.convexity_tests;
            if (convex(points[vertices[previous]], points[vertices[current]],
                       points[vertices[next]])) {
                append_triangle({{vertices[previous], vertices[current], vertices[next]}}, result);
                links[previous].next = next;
                links[next].previous = previous;
                --remaining;
                ++result.work.removed_vertices;
                if (previous == start) {
                    current = next;
                } else {
                    current = previous;
                    ++result.work.backward_steps;
                }
            } else {
                current = next;
                ++result.work.forward_steps;
            }
        }
        assert(links[start].next == end && links[end].next == start &&
               result.triangles.size() - first_triangle == count - 2 &&
               "[FM84 Algorithm 3]: m-2 triangles leave the two extrema");
    }
    assert(result.triangles.size() == points.size() - 2 &&
           result.work.polygon_vertices ==
               points.size() + 2 * (decomposition.polygons.size() - 1) &&
           "[FM84 Theorem 3 tex 467-472]: the complete partition yields n-2 triangles");
    assert(result.work.forward_steps == result.work.backward_steps &&
           result.work.backward_steps <= result.work.removed_vertices &&
           result.work.convexity_tests ==
               result.work.removed_vertices + result.work.forward_steps &&
           result.work.triangle_incidents == 3 * result.triangles.size() &&
           "[FM84 Algorithm 3 analysis tex 457-459]: each deletion pays for at most one backtrack");
    return result;
}

Triangulation triangulate_trapezoids(std::span<const Point> vertices,
                                     const TrapezoidDecomposition& trapezoids) {
    return triangulate_unimonotone(vertices, decompose_unimonotone(vertices, trapezoids));
}

Triangulation triangulate_polygon(const std::vector<Point>& vertices) {
    return triangulate_trapezoids(vertices, compute_trapezoid_decomposition(vertices));
}

}
