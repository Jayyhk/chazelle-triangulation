#include "support/triangulation_production_checks.h"
#include "visualizer/triangulation_svg.h"

#include <cstdio>
#include <sstream>
#include <string_view>

using namespace chazelle;

namespace {

std::size_t occurrences(std::string_view text, std::string_view pattern) {
    std::size_t count = 0;
    std::size_t position = 0;
    while ((position = text.find(pattern, position)) != std::string_view::npos) {
        ++count;
        position += pattern.size();
    }
    return count;
}

std::string render(std::span<const Point> points, std::span<const Triangle> triangles) {
    std::ostringstream output;
    write_triangulation_svg(output, points, triangles);
    test::require_triangulation(static_cast<bool>(output), "SVG output stream succeeded");
    return output.str();
}

void check_exact_normalization() {
    const std::vector<Point> points{{0, 0, 0}, {4, 0, 1}, {4, 4, 2}, {0, 4, 3}};
    const std::array<Triangle, 2> triangles{{{{0, 3, 2}}, {{0, 2, 1}}}};
    const auto expected = render(points, triangles);
    test::require_triangulation(
        expected.find("id=\"triangle-0\" points=\"40,960 40,40 960,40 \"") != std::string::npos &&
            expected.find("id=\"triangle-1\" points=\"40,960 960,40 960,960 \"") !=
                std::string::npos &&
            expected.find("points=\"40,960 960,960 960,40 40,40 \"") != std::string::npos,
        "SVG preserves original triangle corners, boundary order, aspect ratio, and upward y");
    Exact large = Exact(0x1p400) * Exact(0x1p400) * Exact(0x1p400);
    for (const Exact& scale : {Exact(1), large, Exact(1) / large}) {
        auto transformed = points;
        for (Point& point : transformed) {
            point.x = large + scale * point.x;
            point.y = -large + scale * point.y;
        }
        test::require_triangulation(
            render(transformed, triangles) == expected,
            "Exact normalization precedes display conversion, preserving extreme scales and offsets");
    }
    auto decimal = points;
    for (Point& point : decimal) {
        point.x /= 10;
        point.y /= 10;
    }
    test::require_triangulation(render(decimal, triangles) == expected,
                                "Rational coordinates retain the same normalized image");
    test::require_triangulation(occurrences(expected, "<circle ") == points.size() &&
                                    occurrences(expected, "<text ") == points.size() &&
                                    expected.find(">0</text>") != std::string::npos &&
                                    expected.find(">3</text>") != std::string::npos,
                                "The image labels every original vertex");
}

void check_large_image() {
    const std::size_t count = 8193;
    std::vector<Point> points;
    std::vector<Triangle> triangles;
    points.reserve(count);
    triangles.reserve(count - 2);
    for (std::size_t vertex = 0; vertex < count; ++vertex)
        points.push_back({Exact(vertex), Exact(vertex) * Exact(vertex), vertex});
    for (std::size_t vertex = 0; vertex + 2 < count; ++vertex)
        triangles.push_back({{count - 1, vertex + 1, vertex}});
    const auto image = render(points, triangles);
    test::require_triangulation(
        occurrences(image, "<polygon id=\"triangle-") == count - 2 &&
            occurrences(image, "<circle ") == count && occurrences(image, "<text ") == count &&
            occurrences(image, "id=\"boundary\"") == 1 && image.size() < 1024 * count,
        "SVG output contains the complete triangulation with linear size");
}

}

int main() {
    check_exact_normalization();
    check_large_image();
    std::puts("Triangulation visualizer tests passed");
}
