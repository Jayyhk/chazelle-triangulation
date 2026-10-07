#include "animation/trace.h"
#include "animation/triangulate.h"
#include "support/pipeline_fixtures.h"
#include "support/triangulation_production_checks.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

using namespace chazelle;
using chazelle::animation::AnimationTrace;

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

void check_trace(const std::vector<Point>& points,
                 const std::vector<animation::Point>& animation_points,
                 const std::filesystem::path& path) {
    test::require_triangulation(AnimationTrace::current() == nullptr, "Tracing starts disabled");
    const auto ordinary = triangulate_polygon(points);
    std::ostringstream output;
    const auto recorded = animation::triangulate_with_trace(animation_points, output);
    test::check_production_triangulation(points, recorded);
    test::require_triangulation(ordinary.triangles.size() == recorded.triangles.size(),
                                "Tracing preserves the triangle count");
    for (std::size_t i = 0; i < ordinary.triangles.size(); ++i)
        test::require_triangulation(ordinary.triangles[i].vertices ==
                                        recorded.triangles[i].vertices,
                                    "Tracing preserves every triangle and its output order");
    test::require_triangulation(ordinary.vertex_triangles == recorded.vertex_triangles,
                                "The animation copy preserves every vertex's triangle adjacency");
    test::require_triangulation(
        ordinary.work.polygon_vertices == recorded.work.polygon_vertices &&
            ordinary.work.convexity_tests == recorded.work.convexity_tests &&
            ordinary.work.forward_steps == recorded.work.forward_steps &&
            ordinary.work.backward_steps == recorded.work.backward_steps &&
            ordinary.work.removed_vertices == recorded.work.removed_vertices &&
            ordinary.work.triangle_incidents == recorded.work.triangle_incidents,
        "The animation copy performs the same triangulation operations as the pure algorithm");
    test::require_triangulation(AnimationTrace::current() == nullptr,
                                "The trace session restores disabled tracing");
    const auto json = output.str();
    const auto events = occurrences(json, "{\"seq\":");
    test::require_triangulation(occurrences(json, "\"kind\":\"triangle\"") == points.size() - 2,
                                "Every emitted triangle is recorded exactly once");
    test::require_triangulation(
        occurrences(json, "\"kind\":\"convexity_test\"") == recorded.work.convexity_tests &&
            occurrences(json, "\"kind\":\"vertex_remove\"") == recorded.work.removed_vertices &&
            occurrences(json, "\"kind\":\"triangle_cursor\"") == recorded.work.convexity_tests,
        "The replay records every real convexity test, deletion, and cursor step");
    test::require_triangulation(occurrences(json, "\"kind\":\"search_begin\"") ==
                                    occurrences(json, "\"kind\":\"search_end\""),
                                "Every nested ray search has its own completed execution");
    test::require_triangulation(occurrences(json, "\"kind\":\"boundary\"") == 1,
                                "The padded boundary is recorded once");
    std::ofstream file(path);
    file << json;
    test::require_triangulation(static_cast<bool>(file), "The test trace was saved");
    test::require_triangulation(events <= 8192 * points.size(),
                                "The fixture's trace fits a linear event envelope");
}

void check_exact_terms() {
    using animation::Exact;
    using animation::Rational;
    const Exact epsilon = Exact::infinitesimal(0);
    const Exact value = Exact(7) / 3 + epsilon / 5;
    std::size_t numerator = 0;
    std::size_t denominator = 0;
    bool constant = false;
    bool perturbation = false;
    value.visit_terms([&](bool is_denominator, const auto& powers, const Rational& coefficient) {
        if (is_denominator) {
            ++denominator;
            test::require_triangulation(powers == std::array<int, 4>{} &&
                                            coefficient.to_string() == "1",
                                        "The normalized denominator is exact");
        } else {
            ++numerator;
            constant |= powers == std::array<int, 4>{} && coefficient.to_string() == "7/3";
            perturbation |=
                powers == std::array<int, 4>{1, 0, 0, 0} && coefficient.to_string() == "1/5";
        }
    });
    test::require_triangulation(numerator == 2 && denominator == 1 && constant && perturbation,
                                "The visitor preserves rational and infinitesimal coefficients");
    test::require_triangulation(Rational(-13).to_string() == "-13",
                                "Negative integers serialize exactly");
}

void check_sessions() {
    const auto points = test::polygon_fixtures<animation::Point>().front();
    std::ostringstream first;
    std::ostringstream second;
    AnimationTrace outer(first, points);
    const auto events = outer.event_count();
    const auto pure = triangulate_polygon(test::polygon_fixtures().front());
    test::require_triangulation(
        pure.triangles.size() == points.size() - 2 && outer.event_count() == events,
        "The pure algorithm cannot emit events into an active animation trace");
    {
        AnimationTrace inner(second, points);
        test::require_triangulation(AnimationTrace::current() == &inner,
                                    "The nested trace is active");
        inner.finish();
    }
    test::require_triangulation(AnimationTrace::current() == &outer,
                                "The enclosing trace is restored");
    {
        AnimationTrace::QueryRecording disabled(false);
        outer.record("search_scan");
        {
            AnimationTrace::QueryRecording nested(true);
            outer.record("boundary_search");
        }
        test::require_triangulation(outer.event_count() == events,
                                    "Verification searches remain suppressed through nested calls");
        outer.record("search_structure");
        test::require_triangulation(outer.event_count() == events + 1,
                                    "Preprocessed structures remain available to later searches");
    }
    test::require_triangulation(outer.record("search_scan") == events + 1,
                                "The algorithm's query recording resumes after verification");
    outer.finish();
    std::ostringstream failed;
    AnimationTrace trace(failed, points);
    failed.setstate(std::ios::badbit);
    bool rejected = false;
    try {
        trace.finish();
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    test::require_triangulation(rejected, "A failed trace stream cannot report success");
}

}

int main(int argc, char* argv[]) {
    test::require_triangulation(argc == 2, "Provide a test trace directory");
    const std::filesystem::path directory(argv[1]);
    std::filesystem::create_directories(directory);
    check_exact_terms();
    check_sessions();
    std::size_t fixture = 0;
    const auto fixtures = test::polygon_fixtures();
    const auto animation_fixtures = test::polygon_fixtures<animation::Point>();
    for (std::size_t index = 0; index < fixtures.size(); ++index) {
        const auto& points = fixtures[index];
        for (const bool reverse : {false, true}) {
            const auto ordered = test::boundary_order(points, reverse, points.size() / 2, 71);
            const auto animation_ordered =
                test::boundary_order(animation_fixtures[index], reverse, points.size() / 2, 71);
            check_trace(ordered, animation_ordered,
                        directory / ("fixture-" + std::to_string(fixture++) + ".json"));
        }
    }
    for (const std::size_t count : {17U, 65U, 257U, 1026U}) {
        std::vector<Point> points;
        std::vector<animation::Point> animation_points;
        for (std::size_t vertex = 0; vertex < count; ++vertex) {
            points.push_back({Exact(vertex), Exact(vertex) * Exact(vertex), vertex});
            animation_points.push_back({animation::Exact(vertex),
                                        animation::Exact(vertex) * animation::Exact(vertex),
                                        vertex});
        }
        check_trace(points, animation_points,
                    directory / ("grade-" + std::to_string(count) + ".json"));
    }
    auto points = test::polygon_fixtures()[1];
    const Exact large = Exact(0x1p400) * Exact(0x1p400) * Exact(0x1p400);
    const animation::Exact animation_large =
        animation::Exact(0x1p400) * animation::Exact(0x1p400) * animation::Exact(0x1p400);
    for (const bool reciprocal : {false, true}) {
        const Exact scale = reciprocal ? Exact(1) / large : large;
        const animation::Exact animation_scale =
            reciprocal ? animation::Exact(1) / animation_large : animation_large;
        auto transformed = points;
        for (Point& point : transformed) {
            point.x = large + point.x * scale;
            point.y = -large + point.y * scale;
        }
        auto animation_transformed = animation_fixtures[1];
        for (animation::Point& point : animation_transformed) {
            point.x = animation_large + point.x * animation_scale;
            point.y = -animation_large + point.y * animation_scale;
        }
        check_trace(transformed, animation_transformed,
                    directory / ("exact-" + std::to_string(fixture++) + ".json"));
    }
}
