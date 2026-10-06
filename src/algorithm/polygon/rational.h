#pragma once

#include <cassert>
#include <charconv>
#include <cmath>
#include <concepts>
#include <gmp.h>
#include <limits>
#include <utility>

namespace chazelle {

class Rational {
public:
    Rational() {
        mpq_init(value_);
    }

    Rational(double value) : Rational() {
        assert(std::isfinite(value) && "input coordinates must be finite");
        mpq_set_d(value_, value);
        mpq_canonicalize(value_);
    }
    Rational(long double) = delete;

    template <std::integral Integer> Rational(Integer value) : Rational() {
        char decimal_digits[std::numeric_limits<Integer>::digits10 + 4];
        auto conversion =
            std::to_chars(decimal_digits, decimal_digits + sizeof(decimal_digits) - 1, +value);
        assert(conversion.ec == std::errc{});
        *conversion.ptr = '\0';
        [[maybe_unused]] int parse_status = mpq_set_str(value_, decimal_digits, 10);
        assert(parse_status == 0);
    }

    Rational(const Rational& other) : Rational() {
        mpq_set(value_, other.value_);
    }
    Rational(Rational&& other) noexcept : Rational() {
        mpq_swap(value_, other.value_);
    }
    ~Rational() {
        mpq_clear(value_);
    }

    Rational& operator=(const Rational& other) {
        if (this != &other)
            mpq_set(value_, other.value_);
        return *this;
    }
    Rational& operator=(Rational&& other) noexcept {
        mpq_swap(value_, other.value_);
        return *this;
    }

    double to_double() const noexcept {
        return mpq_get_d(value_);
    }

    Rational& operator+=(const Rational& rhs) {
        mpq_add(value_, value_, rhs.value_);
        return *this;
    }
    Rational& operator-=(const Rational& rhs) {
        mpq_sub(value_, value_, rhs.value_);
        return *this;
    }
    Rational& operator*=(const Rational& rhs) {
        mpq_mul(value_, value_, rhs.value_);
        return *this;
    }
    Rational& operator/=(const Rational& rhs) {
        assert(mpq_sgn(rhs.value_) != 0 && "exact division requires nonzero divisor");
        mpq_div(value_, value_, rhs.value_);
        return *this;
    }

    friend Rational operator+(Rational lhs, const Rational& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend Rational operator-(Rational lhs, const Rational& rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend Rational operator*(Rational lhs, const Rational& rhs) {
        lhs *= rhs;
        return lhs;
    }
    friend Rational operator/(Rational lhs, const Rational& rhs) {
        lhs /= rhs;
        return lhs;
    }
    friend Rational operator-(Rational value) {
        mpq_neg(value.value_, value.value_);
        return value;
    }

    friend bool operator==(const Rational& lhs, const Rational& rhs) noexcept {
        return mpq_equal(lhs.value_, rhs.value_) != 0;
    }
    friend bool operator<(const Rational& lhs, const Rational& rhs) noexcept {
        return mpq_cmp(lhs.value_, rhs.value_) < 0;
    }
    friend bool operator>(const Rational& lhs, const Rational& rhs) noexcept {
        return rhs < lhs;
    }
    friend bool operator<=(const Rational& lhs, const Rational& rhs) noexcept {
        return !(rhs < lhs);
    }
    friend bool operator>=(const Rational& lhs, const Rational& rhs) noexcept {
        return !(lhs < rhs);
    }

private:
    mpq_t value_;
};

}
