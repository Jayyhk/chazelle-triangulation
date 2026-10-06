#include "algorithm/triangulation/triangulation.h"
#include "algorithm/visualizer/triangulation_svg.h"
#include "animation/render.h"
#include "animation/triangulate.h"

#include <cassert>
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

std::filesystem::path dated_path(const std::filesystem::path& directory,
                                 const std::filesystem::path& filename) {
    const std::time_t now = std::time(nullptr);
    const std::tm* calendar = now == -1 ? nullptr : std::localtime(&now);
    if (!calendar)
        throw std::runtime_error("Failed to determine the local date for the output filename.");
    std::ostringstream date;
    date << std::put_time(calendar, "%Y-%m-%d");
    return directory /
           (filename.stem().string() + '-' + date.str() + filename.extension().string());
}

template <class Coordinate> Coordinate power_of_ten(std::size_t exponent) {
    Coordinate result = 1;
    Coordinate factor = 10;
    while (exponent != 0) {
        if (exponent % 2 != 0)
            result *= factor;
        exponent /= 2;
        if (exponent != 0)
            factor *= factor;
    }
    return result;
}

template <class Coordinate> std::optional<Coordinate> decimal_coordinate(std::string_view text) {
    std::size_t position = 0;
    const bool negative = !text.empty() && text.front() == '-';
    if (!text.empty() && (text.front() == '+' || text.front() == '-'))
        ++position;
    Coordinate value = 0;
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
            value /= power_of_ten<Coordinate>(fractional_digits + exponent);
        else if (exponent < fractional_digits)
            value /= power_of_ten<Coordinate>(fractional_digits - exponent);
        else
            value *= power_of_ten<Coordinate>(exponent - fractional_digits);
    }
    return negative ? -value : value;
}

struct Options {
    std::optional<std::filesystem::path> image;
    std::optional<std::filesystem::path> animation;
    std::optional<std::filesystem::path> trace;
};

int triangulate_input(const Options& options) {
    std::string token;
    std::size_t count = 0;
    std::vector<chazelle::Point> vertices;
    std::vector<chazelle::animation::Point> animation_vertices;
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
    if (options.trace)
        animation_vertices.reserve(count);
    for (std::size_t vertex = 0; vertex < count; ++vertex) {
        std::string x;
        std::string y;
        if (!(std::cin >> x >> y)) {
            std::cerr << "Expected two coordinates for vertex " << vertex << ".\n";
            return 1;
        }
        auto exact_x = decimal_coordinate<chazelle::Exact>(x);
        auto exact_y = decimal_coordinate<chazelle::Exact>(y);
        if (!exact_x || !exact_y) {
            std::cerr << "Expected decimal coordinates for vertex " << vertex << ".\n";
            return 1;
        }
        vertices.push_back({std::move(*exact_x), std::move(*exact_y), vertex});
        if (options.trace) {
            auto animation_x = decimal_coordinate<chazelle::animation::Exact>(x);
            auto animation_y = decimal_coordinate<chazelle::animation::Exact>(y);
            assert(animation_x && animation_y);
            animation_vertices.push_back(
                {std::move(*animation_x), std::move(*animation_y), vertex});
        }
    }
    if (std::cin >> token) {
        std::cerr << "Unexpected data after the last vertex.\n";
        return 1;
    }
    if (std::cin.bad()) {
        std::cerr << "Failed to read the polygon.\n";
        return 1;
    }
    std::ofstream trace_output;
    if (options.trace) {
        if (options.trace->has_parent_path())
            std::filesystem::create_directories(options.trace->parent_path());
        trace_output.open(*options.trace);
        if (!trace_output)
            throw std::runtime_error("Failed to open animation trace: " + options.trace->string());
    }
    const auto result =
        options.trace
            ? chazelle::animation::triangulate_with_trace(animation_vertices, trace_output)
            : chazelle::triangulate_polygon(vertices);
    if (options.trace) {
        trace_output.close();
        if (!trace_output)
            throw std::runtime_error("Failed to close the animation trace.");
    }
    const auto& image_path = options.image;
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
    if (options.animation) {
        std::filesystem::create_directories(options.animation->parent_path());
        chazelle::render_animation(*options.trace, *options.animation);
        std::cerr << "Animation saved to " << *options.animation << '\n';
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
        constexpr std::string_view usage =
            "Usage: chazelle [--visualize [name.svg]] [--animate [name.mp4]]\n"
            "                [--trace path.json]\n";
        if (argc == 2 && std::string_view(argv[1]) == "--help") {
            std::cout
                << usage
                << "Read n, then n pairs of exact decimal coordinates, from standard input.\n"
                   "Write n-2, then triples of zero-based triangle indices.\n"
                   "Save SVGs in media/images/ and MP4s in media/animations/, "
                   "with -YYYY-MM-DD filenames.\n"
                   "--animate replays the algorithm.\n"
                   "Exact traces are retained in .cache/chazelle/ unless --trace is provided.\n";
            return 0;
        }
        Options options;
        for (int i = 1; i < argc; ++i) {
            const std::string_view argument = argv[i];
            if (argument == "--trace") {
                if (options.trace || i + 1 == argc)
                    throw std::runtime_error("Provide one --trace path.json.");
                options.trace = argv[++i];
                if (options.trace->extension() != ".json")
                    throw std::runtime_error("The animation trace must use the .json extension.");
            } else if (argument == "--visualize" || argument == "--animate") {
                const bool animate = argument == "--animate";
                auto& output = animate ? options.animation : options.image;
                if (output)
                    throw std::runtime_error("Each output flag may be used only once.");
                std::filesystem::path filename = animate ? "algorithm.mp4" : "triangulation.svg";
                if (i + 1 < argc && !std::string_view(argv[i + 1]).starts_with("--"))
                    filename = argv[++i];
                if (filename.extension() != (animate ? ".mp4" : ".svg"))
                    throw std::runtime_error(
                        animate ? "The animation output must use the .mp4 extension."
                                : "The visualization output must use the .svg extension.");
                if (filename.has_parent_path())
                    throw std::runtime_error(animate
                                                 ? "Provide a video filename without a directory; "
                                                   "videos are saved in media/animations/."
                                                 : "Provide an image filename without a directory; "
                                                   "images are saved in media/images/.");
                output = dated_path(
                    std::filesystem::path("media") / (animate ? "animations" : "images"), filename);
            } else {
                std::cerr << usage;
                return 1;
            }
        }
        if (options.animation && !options.trace) {
            auto filename = options.animation->filename();
            filename.replace_extension(".json");
            options.trace = std::filesystem::path(".cache/chazelle") / filename;
        }
        return triangulate_input(options);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
