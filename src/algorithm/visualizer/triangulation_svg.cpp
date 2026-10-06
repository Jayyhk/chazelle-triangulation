#include "triangulation_svg.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <ostream>
#include <vector>

namespace chazelle {

namespace {

struct ImagePoint {
    double x;
    double y;
};

std::vector<ImagePoint> image_points(std::span<const Point> vertices) {
    assert(vertices.size() >= 3);
    Exact left = vertices.front().x;
    Exact right = left;
    Exact bottom = vertices.front().y;
    Exact top = bottom;
    for (const Point& point : vertices) {
        left = std::min(left, point.x);
        right = std::max(right, point.x);
        bottom = std::min(bottom, point.y);
        top = std::max(top, point.y);
    }
    const Exact width = right - left;
    const Exact height = top - bottom;
    const Exact extent = std::max(width, height);
    assert(extent > 0 && "[FM84 Algorithm 1 input]: a simple polygon encloses nonzero area");
    std::vector<ImagePoint> result;
    result.reserve(vertices.size());
    for (const Point& point : vertices) {
        const Exact x = (point.x - left - width / 2) / extent;
        const Exact y = (point.y - bottom - height / 2) / extent;
        assert(x >= Exact(-1) / 2 && x <= Exact(1) / 2 && y >= Exact(-1) / 2 && y <= Exact(1) / 2);
        const ImagePoint pixel{500 + 920 * x.rational_approximation(),
                               500 - 920 * y.rational_approximation()};
        assert(std::isfinite(pixel.x) && std::isfinite(pixel.y));
        result.push_back(pixel);
    }
    return result;
}

}

void write_triangulation_svg(std::ostream& output, std::span<const Point> vertices,
                             std::span<const Triangle> triangles) {
    assert(vertices.size() >= 3 && triangles.size() == vertices.size() - 2 &&
           "[FM84 Theorem 3]: the triangulation contains n-2 original-vertex triangles");
    const auto points = image_points(vertices);
    output.imbue(std::locale::classic());
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    output << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\" "
              "viewBox=\"0 0 1000 1000\">\n"
              "<title>Polygon triangulation</title>\n"
              "<rect width=\"1000\" height=\"1000\" fill=\"white\"/>\n"
              "<g id=\"triangles\" fill=\"#f3f6fa\" stroke=\"#6b809a\" stroke-width=\"1.5\" "
              "stroke-linejoin=\"round\">\n";
    for (std::size_t index = 0; index < triangles.size(); ++index) {
        output << "<polygon id=\"triangle-" << index << "\" points=\"";
        for (const std::size_t vertex : triangles[index].vertices) {
            assert(vertex < points.size() && "[FM84 Algorithm 3]: original vertex indices");
            output << points[vertex].x << ',' << points[vertex].y << ' ';
        }
        output << "\"/>\n";
    }
    output << "</g>\n<polygon id=\"boundary\" fill=\"none\" stroke=\"#17212f\" "
              "stroke-width=\"3\" stroke-linejoin=\"round\" points=\"";
    for (const ImagePoint& point : points)
        output << point.x << ',' << point.y << ' ';
    output << "\"/>\n<g id=\"vertices\" fill=\"#17212f\" font-family=\"sans-serif\" "
              "font-size=\"14\">\n";
    for (std::size_t vertex = 0; vertex < points.size(); ++vertex) {
        const ImagePoint& point = points[vertex];
        output << "<circle id=\"vertex-" << vertex << "\" cx=\"" << point.x << "\" cy=\"" << point.y
               << "\" r=\"3\"/>\n<text x=\"" << point.x + 7 << "\" y=\"" << point.y - 7 << "\">"
               << vertex << "</text>\n";
    }
    output << "</g>\n</svg>\n";
}

}
