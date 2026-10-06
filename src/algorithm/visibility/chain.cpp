#include "chain.h"

#include <utility>

namespace chazelle {

namespace {

std::size_t smallest_covering_grade(std::size_t n) {
    assert(n >= 2 && "[C91 §2.1]: curve needs ≥ 2 vertices");
    std::size_t maximum_grade = 0;
    while ((std::size_t{1} << maximum_grade) + 1 < n)
        ++maximum_grade;
    return maximum_grade;
}

struct PaddedCurve {
    std::vector<Point> vertices;
    std::vector<std::size_t> original_vertex_positions;
};

PaddedCurve pad_with_positions(std::vector<Point> vertices) {
    const std::size_t maximum_grade = smallest_covering_grade(vertices.size());
    const std::size_t padded_vertex_count = (std::size_t{1} << maximum_grade) + 1;

    assert((maximum_grade == 0 || (std::size_t{1} << (maximum_grade - 1)) + 1 < vertices.size()) &&
           "[C91 §4 tex 316]: pad only 'if necessary' — p is minimal");

    assert(padded_vertex_count < 2 * vertices.size() &&
           "[C91 §4 tex 316]: minimal padding less than doubles the curve");

    const std::size_t subdivision_count = padded_vertex_count - vertices.size();
    std::vector<std::size_t> positions;
    positions.reserve(vertices.size());
    const std::size_t first_vertex_tag = vertices.front().index;
    assert(first_vertex_tag != SOS_NONE && padded_vertex_count - 1 < SOS_NONE - first_vertex_tag);
    std::vector<Point> padded;
    padded.reserve(padded_vertex_count);
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        Point vertex = vertices[i];
        assert(i < SOS_NONE - first_vertex_tag && vertex.index == first_vertex_tag + i &&
               "[C91 §2.4 tex 133]: original input tags follow table order");
        vertex.index = first_vertex_tag + padded.size();
        positions.push_back(padded.size());
        padded.push_back(vertex);
        if (i < subdivision_count) {
            const Point& next = vertices[i + 1];
            Point middle{exact_midpoint(vertex.x, next.x), exact_midpoint(vertex.y, next.y),
                         first_vertex_tag + padded.size()};
            assert((middle.x != vertex.x || middle.y != vertex.y) &&
                   (middle.x != next.x || middle.y != next.y) &&
                   "[C91 tex 68/316]: the exact midpoint of a nonnull "
                   "edge is distinct from both endpoints");
            assert((middle.x - vertex.x) * (next.y - vertex.y) ==
                       (middle.y - vertex.y) * (next.x - vertex.x) &&
                   "[C91 tex 316]: subdivision stays on its original edge");
            padded.push_back(middle);
        }
    }
    assert(padded.size() == padded_vertex_count);
    return {std::move(padded), std::move(positions)};
}

}

std::vector<Point> pad_curve(std::vector<Point> vertices) {
    return pad_with_positions(std::move(vertices)).vertices;
}

GradedCurve::GradedCurve(std::vector<Point> vertices) {
    PaddedCurve padded = pad_with_positions(std::move(vertices));
    original_vertex_positions_ = std::move(padded.original_vertex_positions);
    Polygon input_curve(std::move(padded.vertices));

    const std::size_t n = input_curve.num_vertices();
    maximum_grade_ = smallest_covering_grade(n);
    assert((std::size_t{1} << maximum_grade_) + 1 == n &&
           "[C91 §4 tex 316]: padded curve has n = 2^p + 1 vertices");

    chains_.resize(maximum_grade_ + 1);

    chains_[0].reserve(std::size_t{1} << maximum_grade_);
    for (std::size_t i = 0; i < (std::size_t{1} << maximum_grade_); ++i)
        chains_[0].push_back(input_curve.subchain(i, 2));

    for (std::size_t grade = 1; grade <= maximum_grade_; ++grade) {
        const std::size_t count = std::size_t{1} << (maximum_grade_ - grade);
        chains_[grade].reserve(count);
        for (std::size_t i = 0; i < count; ++i)
            chains_[grade].push_back(
                Polygon(chains_[grade - 1][2 * i], chains_[grade - 1][2 * i + 1]));
    }

    for (std::size_t grade = 0; grade <= maximum_grade_; ++grade) {
        assert(chains_[grade].size() == (std::size_t{1} << (maximum_grade_ - grade)) &&
               "[C91 §4 tex 319 (ii)]: 2^{p−λ} chains in grade λ");
        for ([[maybe_unused]] const Polygon& c : chains_[grade])
            assert(c.num_vertices() == (std::size_t{1} << grade) + 1 &&
                   "[C91 §4 tex 318 (i)]: a grade-λ chain has 2^λ + 1 "
                   "vertices");
    }

    assert(chains_[maximum_grade_].size() == 1 &&
           "[C91 §4 tex 319 (ii)]: grade p has a single chain");
    assert(chains_[maximum_grade_][0].num_vertices() == n &&
           &chains_[maximum_grade_][0].vertex(0) == &input_curve.vertex(0) &&
           "[C91 §4 tex 316]: the grade-p chain is the whole curve P");
    assert(chains_[maximum_grade_][0].max_y_vertex() == input_curve.max_y_vertex() &&
           chains_[maximum_grade_][0].min_y_vertex() == input_curve.min_y_vertex() &&
           "[C91 §2 tex 47]: combined y-extremes match the direct scan");
}

}
