#include "unimonotone.h"
#include "../trace.h"

#include <algorithm>
#include <cassert>
#include <utility>

namespace chazelle::animation {

namespace {

struct Vertex {
    std::size_t previous = NONE;
    std::size_t next = NONE;
    std::size_t remaining_trapezoids = 0;
    bool done = false;
};

struct Split {
    std::size_t top;
    std::size_t bottom;
    std::size_t saved_next;
    std::size_t saved_previous;
};

void append_polygon(std::size_t first, [[maybe_unused]] std::size_t last,
                    std::span<const Point> points, std::span<const Vertex> vertices,
                    UnimonotoneDecomposition& result) {
    assert(vertices[last].next == first && vertices[first].previous == last &&
           "[FM84 Algorithm 2 tex 353-380]: first and last close the current polygon");
    UnimonotonePolygon polygon;
    polygon.top_vertex = first;
    polygon.bottom_vertex = first;
    std::size_t current = first;
    do {
        assert(polygon.vertices.size() < points.size() &&
               "[FM84 Algorithm 2]: each subpolygon contains a vertex at most once");
        polygon.vertices.push_back(current);
        ++result.work.emitted_vertices;
        if (point_y_above(points[current], points[polygon.top_vertex]))
            polygon.top_vertex = current;
        if (point_y_below(points[current], points[polygon.bottom_vertex]))
            polygon.bottom_vertex = current;
        current = vertices[current].next;
    } while (current != first);
    assert(polygon.vertices.size() >= 3 &&
           "[FM84 Algorithm 2 tex 328-339]: a nonadjacent diagonal produces two polygons");
    const std::size_t top = polygon.top_vertex;
    const std::size_t bottom = polygon.bottom_vertex;
    assert((vertices[top].next == bottom || vertices[bottom].next == top) &&
           "[FM84 tex 334-339, 386-396]: a unimonotone polygon has adjacent y extrema");
    const bool increasing = vertices[top].next == bottom;
    const std::size_t end = increasing ? top : bottom;
    current = increasing ? bottom : top;
    while (current != end) {
        assert((increasing ? point_y_above(points[vertices[current].next], points[current])
                           : point_y_below(points[vertices[current].next], points[current])) &&
               "[FM84 tex 334-339]: the remaining chain is monotone in symbolic y");
        current = vertices[current].next;
    }
    if (auto* trace = AnimationTrace::current())
        trace->piece(polygon.vertices);
    result.polygons.push_back(std::move(polygon));
}

}

UnimonotoneDecomposition decompose_unimonotone(std::span<const Point> points,
                                               const TrapezoidDecomposition& trapezoids) {
    animation_checkpoint("unimonotone");
    const std::size_t count = points.size();
    assert(count >= 3 && trapezoids.vertex_trapezoids.size() == count &&
           "[FM84 Algorithm 2 input]: a simple polygon and its trapezoid table are required");
    Exact twice_area = 0;
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        const Point& a = points[vertex];
        const Point& b = points[(vertex + 1) % count];
        assert((a.x != b.x || a.y != b.y) &&
               "[FM84 Algorithm 1 input]: consecutive polygon vertices are distinct");
        twice_area += a.x * b.y - a.y * b.x;
    }
    assert(twice_area != 0 && "[FM84 Algorithm 1 input]: a simple polygon encloses positive area");
    std::vector<Vertex> vertices(count);
    std::vector<bool> seen(trapezoids.trapezoids.size(), false);
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        const std::size_t before = vertex == 0 ? count - 1 : vertex - 1;
        const std::size_t after = (vertex + 1) % count;
        vertices[vertex].previous = twice_area < 0 ? before : after;
        vertices[vertex].next = twice_area < 0 ? after : before;
        const VertexTrapezoids& entry = trapezoids.vertex_trapezoids[vertex];
        assert(entry.count <= 2 && "[FM84 tex 192-195]: a vertex points to at most two trapezoids");
        vertices[vertex].remaining_trapezoids = entry.count;
        for (std::size_t slot = 0; slot < entry.count; ++slot) {
            const std::size_t index = entry.trapezoids[slot];
            assert(index < trapezoids.trapezoids.size() && !seen[index] &&
                   "[FM84 tex 192-195]: each trapezoid appears once in the vertex table");
            seen[index] = true;
            [[maybe_unused]] const Trapezoid& trapezoid = trapezoids.trapezoids[index];
            assert(trapezoid.top_vertex == vertex && trapezoid.bottom_vertex < count &&
                   trapezoid.left_edge < count && trapezoid.right_edge < count &&
                   "[FM84 Algorithm 1 output]: trapezoid indices refer to the input polygon");
            assert(point_y_above(points[vertex], points[trapezoid.bottom_vertex]) &&
                   "[FM84 tex 192-195]: each trapezoid is stored under its top vertex");
        }
    }
    assert(std::all_of(seen.begin(), seen.end(), [](bool present) { return present; }) &&
           "[FM84 tex 192-195]: the vertex table contains every trapezoid");

    UnimonotoneDecomposition result;
    result.work.vertices_initialized = count;
    std::vector<Split> splits;
    std::size_t first = 0;
    std::size_t last = vertices[first].previous;
    std::size_t current = first;
    for (;;) {
        while (!vertices[current].done) {
            vertices[current].done = true;
            ++result.work.vertex_visits;
            if (auto* trace = AnimationTrace::current())
                trace->record("partition_vertex", {{"vertex", current}});
            const VertexTrapezoids& entry = trapezoids.vertex_trapezoids[current];
            std::size_t diagonal = NONE;
            while (vertices[current].remaining_trapezoids > 0) {
                const std::size_t index =
                    entry.trapezoids[--vertices[current].remaining_trapezoids];
                ++result.work.trapezoids_examined;
                const std::size_t bottom = trapezoids.trapezoids[index].bottom_vertex;
                const std::size_t original_previous = current == 0 ? count - 1 : current - 1;
                const std::size_t original_next = (current + 1) % count;
                const bool needs_diagonal = bottom != original_previous && bottom != original_next;
                if (auto* trace = AnimationTrace::current())
                    trace->record("trapezoid_test", {{"vertex", current},
                                                     {"trapezoid", index},
                                                     {"bottom", bottom},
                                                     {"split", needs_diagonal}});
                if (needs_diagonal) {
                    assert(bottom != vertices[current].next &&
                           bottom != vertices[current].previous &&
                           "[FM84 tex 328-339]: a class-B diagonal has not already been inserted");
                    diagonal = index;
                    break;
                }
            }
            if (diagonal == NONE) {
                current = vertices[current].next;
                continue;
            }
            const std::size_t bottom = trapezoids.trapezoids[diagonal].bottom_vertex;
            result.diagonals.push_back({current, bottom, diagonal});
            if (auto* trace = AnimationTrace::current())
                trace->diagonal(current, bottom);
            splits.push_back({current, bottom, vertices[current].next, vertices[bottom].previous});
            result.work.maximum_stack_size =
                std::max(result.work.maximum_stack_size, splits.size());
            vertices[current].next = bottom;
            vertices[bottom].previous = current;
            first = bottom;
            last = current;
            current = first;
        }
        append_polygon(first, last, points, vertices, result);
        if (splits.empty())
            break;
        const Split split = splits.back();
        splits.pop_back();
        if (auto* trace = AnimationTrace::current())
            trace->record("partition_return",
                          {{"top", split.top}, {"bottom", split.bottom}, {"depth", splits.size()}});
        vertices[split.top].done = false;
        vertices[split.bottom].done = false;
        vertices[split.top].next = split.saved_next;
        vertices[split.bottom].previous = split.saved_previous;
        vertices[split.bottom].next = split.top;
        vertices[split.top].previous = split.bottom;
        first = split.top;
        last = split.bottom;
        current = first;
    }
    assert(result.diagonals.size() <= count - 3 &&
           "[FM84 Algorithm 2]: noncrossing diagonals number at most n-3");
    assert(result.polygons.size() == result.diagonals.size() + 1 &&
           result.work.emitted_vertices == count + 2 * result.diagonals.size() &&
           "[FM84 Algorithm 2]: each split adds one polygon and two boundary occurrences");
    assert(result.work.vertex_visits <= count + 2 * result.diagonals.size() &&
           result.work.trapezoids_examined == trapezoids.trapezoids.size() &&
           "[FM84 tex 466-474]: splitting examines each trapezoid once and takes O(n)");
    return result;
}

UnimonotoneDecomposition compute_unimonotone_decomposition(const std::vector<Point>& vertices) {
    return decompose_unimonotone(vertices, compute_trapezoid_decomposition(vertices));
}

}
