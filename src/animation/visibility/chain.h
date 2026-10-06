#pragma once

#include "../polygon/polygon.h"

#include <cassert>
#include <cstddef>
#include <span>
#include <vector>

namespace chazelle::animation {

inline constexpr std::size_t BETA_NUM = 1;
inline constexpr std::size_t BETA_DEN = 5;

inline constexpr std::size_t ceil_beta(std::size_t grade) noexcept {
    return (grade * BETA_NUM + BETA_DEN - 1) / BETA_DEN;
}

std::vector<Point> pad_curve(std::vector<Point> vertices);

class GradedCurve {
public:
    explicit GradedCurve(std::vector<Point> vertices);

    const Polygon& curve() const noexcept {
        return chains_[maximum_grade_][0];
    }

    std::span<const std::size_t> original_vertex_positions() const noexcept {
        return original_vertex_positions_;
    }

    std::size_t maximum_grade() const noexcept {
        return maximum_grade_;
    }

    std::size_t num_grades() const noexcept {
        return maximum_grade_ + 1;
    }

    std::size_t num_chains(std::size_t grade) const noexcept {
        assert(grade <= maximum_grade_ && "[C91 §4 tex 320 (iii)]: grades are 0, 1, ..., p");
        return std::size_t{1} << (maximum_grade_ - grade);
    }

    const Polygon& chain(std::size_t grade, std::size_t i) const noexcept {
        assert(grade <= maximum_grade_ && "[C91 §4 tex 320 (iii)]: grades are 0, 1, ..., p");
        assert(i < num_chains(grade) && "[C91 §4 tex 319 (ii)]: only 2^{p−λ} chains in grade λ");
        return chains_[grade][i];
    }

private:
    std::size_t maximum_grade_ = 0;

    std::vector<std::vector<Polygon>> chains_;
    std::vector<std::size_t> original_vertex_positions_;
};

}
