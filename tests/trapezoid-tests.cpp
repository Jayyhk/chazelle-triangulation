#include "algorithm/triangulation/trapezoids.h"
#include "algorithm/visibility/naive_visibility.h"
#include "algorithm/visibility/up_phase.h"
#include "support/triangulation_checks.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <map>
#include <numeric>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

using namespace chazelle;

namespace {

using TrapezoidKey = std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>;

std::vector<TrapezoidKey> keys(const std::vector<Trapezoid>& trapezoids) {
    std::vector<TrapezoidKey> result;
    result.reserve(trapezoids.size());
    for (const Trapezoid& trapezoid : trapezoids)
        result.emplace_back(trapezoid.top_vertex, trapezoid.bottom_vertex, trapezoid.left_edge,
                            trapezoid.right_edge);
    std::sort(result.begin(), result.end());
    return result;
}

Exact height(const Point& point) {
    return point.y - Exact(point.index) * Exact::infinitesimal(3);
}

Exact side_x(const Polygon& curve, std::span<const std::size_t> original_positions,
             std::size_t edge, const Exact& y, bool perturbed) {
    const Point& a = curve.vertex(original_positions[edge]);
    const Point& b = curve.vertex(original_positions[(edge + 1) % original_positions.size()]);
    const Exact ay = perturbed ? height(a) : a.y;
    const Exact by = perturbed ? height(b) : b.y;
    assert(ay != by);
    return a.x + (y - ay) * (b.x - a.x) / (by - ay);
}

std::vector<Trapezoid> reference_trapezoids(const Polygon& curve,
                                            std::span<const std::size_t> original_positions) {
    const std::size_t count = original_positions.size();
    std::vector<std::size_t> order(count);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return height(curve.vertex(original_positions[a])) >
               height(curve.vertex(original_positions[b]));
    });
    std::vector<Trapezoid> result;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> previous;
    for (std::size_t band = 0; band + 1 < count; ++band) {
        const std::size_t top = order[band];
        const std::size_t bottom = order[band + 1];
        const Exact y = (height(curve.vertex(original_positions[top])) +
                         height(curve.vertex(original_positions[bottom]))) /
                        2;
        std::vector<std::pair<Exact, std::size_t>> crossings;
        crossings.reserve(count);
        for (std::size_t edge = 0; edge < count; ++edge) {
            const Point& a = curve.vertex(original_positions[edge]);
            const Point& b = curve.vertex(original_positions[(edge + 1) % count]);
            const Exact ay = height(a);
            const Exact by = height(b);
            if ((ay < y && y < by) || (by < y && y < ay))
                crossings.emplace_back(a.x + (y - ay) * (b.x - a.x) / (by - ay), edge);
        }
        std::sort(crossings.begin(), crossings.end());
        assert(crossings.size() % 2 == 0);
        std::map<std::pair<std::size_t, std::size_t>, std::size_t> current;
        for (std::size_t crossing = 0; crossing < crossings.size(); crossing += 2) {
            const auto sides =
                std::pair{crossings[crossing].second, crossings[crossing + 1].second};
            assert(crossings[crossing].first < crossings[crossing + 1].first);
            const auto old = previous.find(sides);
            if (old == previous.end()) {
                current[sides] = result.size();
                result.push_back({top, bottom, sides.first, sides.second});
            } else {
                current[sides] = old->second;
                result[old->second].bottom_vertex = bottom;
            }
        }
        previous = std::move(current);
    }
    return result;
}

void check_result(const Polygon& curve, std::span<const std::size_t> original_positions,
                  const Submap& map, const TrapezoidDecomposition& result) {
    const auto expected = reference_trapezoids(curve, original_positions);
    if (keys(result.trapezoids) != keys(expected)) {
        std::fprintf(stderr, "Mismatch for %zu original vertices (%zu padded):\n",
                     original_positions.size(), curve.num_vertices());
        for (const auto& [top, bottom, left, right] : keys(result.trapezoids))
            std::fprintf(stderr, "actual %zu %zu %zu %zu\n", top, bottom, left, right);
        for (const auto& [top, bottom, left, right] : keys(expected))
            std::fprintf(stderr, "expect %zu %zu %zu %zu\n", top, bottom, left, right);
        assert(false &&
               "[FM84 Theorem 5]: visibility extraction equals exact polygon decomposition");
    }
    assert(result.vertex_trapezoids.size() == original_positions.size());
    std::vector<bool> seen(result.trapezoids.size(), false);
    Exact trapezoid_area = 0;
    Exact polygon_twice_area = 0;
    for (std::size_t vertex = 0; vertex < original_positions.size(); ++vertex) {
        const Point& a = curve.vertex(original_positions[vertex]);
        const Point& b = curve.vertex(original_positions[(vertex + 1) % original_positions.size()]);
        polygon_twice_area += a.x * b.y - a.y * b.x;
        const auto& entry = result.vertex_trapezoids[vertex];
        assert(entry.count <= 2);
        for (std::size_t slot = 0; slot < entry.count; ++slot) {
            const std::size_t index = entry.trapezoids[slot];
            assert(index < result.trapezoids.size() && !seen[index]);
            seen[index] = true;
            const auto& trapezoid = result.trapezoids[index];
            assert(trapezoid.top_vertex == vertex &&
                   trapezoid.bottom_vertex < original_positions.size() &&
                   trapezoid.left_edge < original_positions.size() &&
                   trapezoid.right_edge < original_positions.size());
            assert(point_y_above(curve.vertex(original_positions[trapezoid.top_vertex]),
                                 curve.vertex(original_positions[trapezoid.bottom_vertex])));
            const Exact top_y = curve.vertex(original_positions[trapezoid.top_vertex]).y;
            const Exact bottom_y = curve.vertex(original_positions[trapezoid.bottom_vertex]).y;
            if (top_y > bottom_y) {
                const Exact top_width =
                    side_x(curve, original_positions, trapezoid.right_edge, top_y, false) -
                    side_x(curve, original_positions, trapezoid.left_edge, top_y, false);
                const Exact bottom_width =
                    side_x(curve, original_positions, trapezoid.right_edge, bottom_y, false) -
                    side_x(curve, original_positions, trapezoid.left_edge, bottom_y, false);
                assert(top_width >= 0 && bottom_width >= 0);
                trapezoid_area += (top_y - bottom_y) * (top_width + bottom_width) / 2;
            }
        }
        if (entry.count == 2) {
            const Exact y = height(a) - Exact::infinitesimal(2);
            const Trapezoid& left = result.trapezoids[entry.trapezoids[0]];
            const Trapezoid& right = result.trapezoids[entry.trapezoids[1]];
            assert(side_x(curve, original_positions, left.left_edge, y, true) <
                   side_x(curve, original_positions, right.left_edge, y, true));
        }
    }
    assert(2 * trapezoid_area ==
           (polygon_twice_area < 0 ? -polygon_twice_area : polygon_twice_area));
    assert(std::all_of(seen.begin(), seen.end(), [](bool present) { return present; }));
    assert(result.work.boundary_edges == curve.num_edges());
    assert(result.work.regions == map.num_live_nodes());
    assert(result.work.arcs <= map.num_live_arcs());
    assert(result.work.arc_edges <= 4 * map.num_live_arcs());
    assert(result.work.chords <= map.num_live_chords());
    assert(result.work.joined_regions <= map.num_live_nodes());
}

void check_polygon(const std::vector<Point>& vertices, bool check_pipeline) {
    const Polygon curve(pad_curve(vertices));
    const std::size_t subdivisions = curve.num_vertices() - vertices.size();
    std::vector<std::size_t> original_positions;
    original_positions.reserve(vertices.size());
    for (std::size_t vertex = 0; vertex < vertices.size(); ++vertex)
        original_positions.push_back(vertex + std::min(vertex, subdivisions));
    const Submap map = build_full_visibility_map(curve);
    const auto result = extract_trapezoids(curve, map, original_positions);
    check_result(curve, original_positions, map, result);
    test::check_triangulation(vertices, triangulate_trapezoids(vertices, result), false);
    const Polygon unpadded(vertices);
    std::vector<std::size_t> identity(vertices.size());
    std::iota(identity.begin(), identity.end(), 0);
    const Submap unpadded_map = build_full_visibility_map(unpadded);
    check_result(unpadded, identity, unpadded_map,
                 extract_trapezoids(unpadded, unpadded_map, identity));
    if (check_pipeline) {
        const auto pipeline = compute_trapezoid_decomposition(vertices);
        assert(keys(pipeline.trapezoids) == keys(result.trapezoids));
        const GradedCurve graded(vertices);
        assert(std::equal(original_positions.begin(), original_positions.end(),
                          graded.original_vertex_positions().begin(),
                          graded.original_vertex_positions().end()));
        test::check_triangulation(vertices, triangulate_trapezoids(vertices, pipeline), false);
    }
}

std::vector<std::vector<Point>> fixtures() {
    return {
        {{0, 0, 0}, {4, 0, 1}, {2, 3, 2}},
        {{0, 0, 0}, {4, 0, 1}, {4, 4, 2}, {0, 4, 3}},
        {{0, 0, 0}, {2, 0, 1}, {4, 0, 2}, {4, 2, 3}, {4, 4, 4}, {0, 4, 5}},
        {{0, 0, 0}, {6, 0, 1}, {6, 6, 2}, {4, 6, 3}, {4, 2, 4}, {2, 2, 5}, {2, 6, 6}, {0, 6, 7}},
        {{0, 0, 0},
         {8, 0, 1},
         {8, 8, 2},
         {6, 8, 3},
         {6, 2, 4},
         {2, 2, 5},
         {2, 6, 6},
         {4, 6, 7},
         {4, 8, 8},
         {0, 8, 9}},
        {{0, 0, 0}, {6, 1, 1}, {4, 5, 2}, {3, 2, 3}, {1, 6, 4}, {-1, 3, 5}},
    };
}

void check_fixture_orders() {
    for (const auto& fixture : fixtures()) {
        for (bool reverse : {false, true}) {
            for (std::size_t first = 0; first < fixture.size(); ++first) {
                std::vector<Point> vertices = fixture;
                if (reverse)
                    std::reverse(vertices.begin(), vertices.end());
                std::rotate(vertices.begin(), vertices.begin() + static_cast<std::ptrdiff_t>(first),
                            vertices.end());
                for (std::size_t vertex = 0; vertex < vertices.size(); ++vertex)
                    vertices[vertex].index = vertex + 71;
                check_polygon(vertices, first == 0);
            }
        }
    }
}

void check_random_polygons() {
    std::mt19937 random(337);
    for (std::size_t test = 0; test < 120; ++test) {
        const std::size_t count = 3 + random() % 25;
        std::vector<Point> vertices;
        vertices.reserve(count);
        vertices.push_back({-1, -3, 0});
        for (std::size_t vertex = 1; vertex + 1 < count; ++vertex)
            vertices.push_back({Exact(vertex) / 3, Exact(random() % 17) / 5, vertex});
        vertices.push_back({Exact(count), -3, count - 1});
        if (test % 2 == 0)
            std::reverse(vertices.begin(), vertices.end());
        std::rotate(vertices.begin(),
                    vertices.begin() + static_cast<std::ptrdiff_t>(random() % count),
                    vertices.end());
        for (std::size_t vertex = 0; vertex < count; ++vertex) {
            if (test % 3 == 0) {
                const Exact x = vertices[vertex].x;
                vertices[vertex].x = -vertices[vertex].y;
                vertices[vertex].y = x;
            }
            vertices[vertex].index = vertex + 103;
        }
        check_polygon(vertices, test < 8);
    }
}

void check_radial_polygons() {
    std::mt19937 random(7349);
    for (std::size_t test = 0; test < 120; ++test) {
        const std::size_t side_count = 2 + random() % 7;
        std::vector<Point> vertices;
        vertices.reserve(4 * side_count);
        for (std::size_t quadrant = 0; quadrant < 4; ++quadrant) {
            for (std::size_t step = 0; step < side_count; ++step) {
                const Exact scale = Exact(1 + random() % 20) / 7;
                Exact x = Exact(side_count - step) * scale;
                Exact y = Exact(step) * scale;
                for (std::size_t turn = 0; turn < quadrant; ++turn) {
                    const Exact old_x = x;
                    x = -y;
                    y = old_x;
                }
                vertices.push_back({x, y, vertices.size()});
            }
        }
        if (test % 2 == 0)
            std::reverse(vertices.begin(), vertices.end());
        std::rotate(vertices.begin(),
                    vertices.begin() + static_cast<std::ptrdiff_t>(random() % vertices.size()),
                    vertices.end());
        for (std::size_t vertex = 0; vertex < vertices.size(); ++vertex)
            vertices[vertex].index = vertex + 233;
        check_polygon(vertices, test < 8);
    }
}

void check_exhaustive_heights() {
    for (std::size_t sequence = 0; sequence < 81; ++sequence) {
        std::size_t heights = sequence;
        std::vector<Point> vertices{{-1, -3, 0}};
        vertices.reserve(6);
        for (std::size_t vertex = 0; vertex < 4; ++vertex) {
            vertices.push_back({Exact(vertex), Exact(heights % 3), vertex + 1});
            heights /= 3;
        }
        vertices.push_back({4, -3, 5});
        check_polygon(vertices, false);
    }
}

void check_extreme_coordinates() {
    for (const auto& fixture : fixtures()) {
        for (bool transpose : {false, true}) {
            auto vertices = fixture;
            for (Point& point : vertices) {
                point.x *= Exact(0x1p-400);
                point.y *= Exact(0x1p400);
                if (transpose)
                    std::swap(point.x, point.y);
            }
            check_polygon(vertices, true);
        }
    }
}

void check_multiple_subdivisions() {
    for (const auto& fixture : fixtures()) {
        std::vector<Point> vertices;
        vertices.reserve(5 * fixture.size());
        std::vector<std::size_t> original_positions;
        original_positions.reserve(fixture.size());
        for (std::size_t edge = 0; edge + 1 < fixture.size(); ++edge) {
            original_positions.push_back(vertices.size());
            const std::size_t pieces = 3 + edge % 3;
            const Point& a = fixture[edge];
            const Point& b = fixture[edge + 1];
            for (std::size_t piece = 0; piece < pieces; ++piece) {
                const Exact t = Exact(piece) / Exact(pieces);
                vertices.push_back({a.x + t * (b.x - a.x), a.y + t * (b.y - a.y),
                                    SOS_NONE - 257 + vertices.size()});
            }
        }
        original_positions.push_back(vertices.size());
        const Point& last = fixture.back();
        vertices.push_back({last.x, last.y, SOS_NONE - 257 + vertices.size()});
        const Polygon curve(std::move(vertices));
        const Submap map = build_full_visibility_map(curve);
        check_result(curve, original_positions, map,
                     extract_trapezoids(curve, map, original_positions));
    }
}

}

int main() {
    check_fixture_orders();
    check_random_polygons();
    check_radial_polygons();
    check_exhaustive_heights();
    check_extreme_coordinates();
    check_multiple_subdivisions();
    std::puts("[FM84 trapezoid extraction tests]: all passed");
}
