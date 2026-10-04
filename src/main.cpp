#include "triangulation/triangulation.h"
#include "visualizer/triangulation_svg.h"

#include <charconv>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

std::filesystem::path dated_image_path(const std::filesystem::path& filename) {
    const std::time_t now = std::time(nullptr);
    const std::tm* calendar = now == -1 ? nullptr : std::localtime(&now);
    if (!calendar)
        throw std::runtime_error("Failed to determine the local date for the image filename.");
    std::ostringstream date;
    date << std::put_time(calendar, "%Y-%m-%d");
    return std::filesystem::path("images") / (filename.stem().string() + '-' + date.str() + ".svg");
}

chazelle::Exact power_of_ten(std::size_t exponent) {
    chazelle::Exact result = 1;
    chazelle::Exact factor = 10;
    while (exponent != 0) {
        if (exponent % 2 != 0)
            result *= factor;
        exponent /= 2;
        if (exponent != 0)
            factor *= factor;
    }
    return result;
}

std::optional<chazelle::Exact> decimal_coordinate(std::string_view text) {
    std::size_t position = 0;
    const bool negative = !text.empty() && text.front() == '-';
    if (!text.empty() && (text.front() == '+' || text.front() == '-'))
        ++position;
    chazelle::Exact value = 0;
    std::size_t fractional_digits = 0;
    bool decimal_point = false;
    bool has_digits = false;
    while (position < text.size()) {
        const char character = text[position];
        if (character >= '0' && character <= '9') {
            value = value * 10 + (character - '0');
            fractional_digits += decimal_point;
            has_digits = true;
        } else if (character == '.' && !decimal_point) {
            decimal_point = true;
        } else {
            break;
        }
        ++position;
    }
    if (!has_digits)
        return std::nullopt;
    std::size_t exponent = 0;
    bool negative_exponent = false;
    if (position < text.size()) {
        if (text[position] != 'e' && text[position] != 'E')
            return std::nullopt;
        ++position;
        if (position < text.size() && (text[position] == '+' || text[position] == '-')) {
            negative_exponent = text[position] == '-';
            ++position;
        }
        const char* end = text.data() + text.size();
        const auto conversion = std::from_chars(text.data() + position, end, exponent);
        if (conversion.ec != std::errc{} || conversion.ptr != end)
            return std::nullopt;
    }
    if (negative_exponent && exponent > std::numeric_limits<std::size_t>::max() - fractional_digits)
        return std::nullopt;
    if (value != 0) {
        if (negative_exponent)
            value /= power_of_ten(fractional_digits + exponent);
        else if (exponent < fractional_digits)
            value /= power_of_ten(fractional_digits - exponent);
        else
            value *= power_of_ten(exponent - fractional_digits);
    }
    return negative ? -value : value;
}

int triangulate_input(const std::optional<std::filesystem::path>& image_path) {
    std::string token;
    std::size_t count = 0;
    std::vector<chazelle::Point> vertices;
    if (!(std::cin >> token)) {
        std::cerr << "Expected a vertex count.\n";
        return 1;
    }
    const auto conversion = std::from_chars(token.data(), token.data() + token.size(), count);
    if (conversion.ec != std::errc{} || conversion.ptr != token.data() + token.size() ||
        count < 3 || count > vertices.max_size()) {
        std::cerr << "Expected an integer vertex count of at least three.\n";
        return 1;
    }
    vertices.reserve(count);
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        std::string x;
        std::string y;
        if (!(std::cin >> x >> y)) {
            std::cerr << "Expected two coordinates for vertex " << vertex << ".\n";
            return 1;
        }
        auto exact_x = decimal_coordinate(x);
        auto exact_y = decimal_coordinate(y);
        if (!exact_x || !exact_y) {
            std::cerr << "Expected decimal coordinates for vertex " << vertex << ".\n";
            return 1;
        }
        vertices.push_back({std::move(*exact_x), std::move(*exact_y), vertex});
    }
    if (std::cin >> token) {
        std::cerr << "Unexpected data after the last vertex.\n";
        return 1;
    }
    if (std::cin.bad()) {
        std::cerr << "Failed to read the polygon.\n";
        return 1;
    }
    const auto result = chazelle::triangulate_polygon(vertices);
    if (image_path) {
        std::filesystem::create_directories(image_path->parent_path());
        std::ofstream image(*image_path);
        if (!image) {
            std::cerr << "Failed to open image file: " << *image_path << '\n';
            return 1;
        }
        chazelle::write_triangulation_svg(image, vertices, result.triangles);
        image.close();
        if (!image) {
            std::cerr << "Failed to write image file: " << *image_path << '\n';
            return 1;
        }
    }
    std::cout << result.triangles.size() << '\n';
    for (const auto& triangle : result.triangles)
        std::cout << triangle.vertices[0] << ' ' << triangle.vertices[1] << ' '
                  << triangle.vertices[2] << '\n';
    if (!std::cout) {
        std::cerr << "Failed to write the triangles.\n";
        return 1;
    }
    return 0;
}

}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    try {
        std::optional<std::filesystem::path> image_path;
        if (argc == 2 && std::string_view(argv[1]) == "--help") {
            std::cout
                << "Usage: chazelle [--visualize [name.svg]]\n"
                   "Read n, then n pairs of exact decimal coordinates, from standard input.\n"
                   "Write n-2, then triples of zero-based triangle indices.\n"
                   "Save images as images/<name>-YYYY-MM-DD.svg; the default name is triangulation.\n";
            return 0;
        }
        if (argc != 1) {
            if (argc > 3 || std::string_view(argv[1]) != "--visualize") {
                std::cerr << "Usage: chazelle [--visualize [name.svg]]\n";
                return 1;
            }
            const std::filesystem::path filename = argc == 3 ? argv[2] : "triangulation.svg";
            if (filename.extension() != ".svg") {
                std::cerr << "The visualization output must use the .svg extension.\n";
                return 1;
            }
            if (filename.has_parent_path()) {
                std::cerr
                    << "Provide an image filename without a directory; images are saved in images/.\n";
                return 1;
            }
            image_path = dated_image_path(filename);
        }
        return triangulate_input(image_path);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
