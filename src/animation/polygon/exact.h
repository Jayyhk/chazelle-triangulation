#pragma once

#include "rational.h"

#include <array>
#include <map>
#include <memory>

namespace chazelle::animation {

class Exact {
    using Powers = std::array<int, 4>;
    using Polynomial = std::map<Powers, Rational>;

    struct Fraction {
        Polynomial numerator;
        Polynomial denominator;
    };

public:
    Exact() = default;
    Exact(double value) : rational_(value) {}
    Exact(long double) = delete;
    template <std::integral Integer> Exact(Integer value) : rational_(value) {}

    double rational_approximation() const {
        assert(!fraction_ && "display coordinates must be finite rational numbers");
        return rational_.to_double();
    }

    template <class Visitor> void visit_terms(Visitor&& visitor) const {
        if (fraction_) {
            for (const auto& [powers, coefficient] : fraction_->numerator)
                visitor(false, powers, coefficient);
            for (const auto& [powers, coefficient] : fraction_->denominator)
                visitor(true, powers, coefficient);
        } else {
            visitor(false, Powers{}, rational_);
            visitor(true, Powers{}, Rational{1});
        }
    }

    static Exact infinitesimal(std::size_t level) {
        assert(level < 4);
        Powers powers{};
        powers[level] = 1;
        return from_fraction({{{powers, Rational{1}}}, constant(Rational{1})});
    }

    Exact& operator+=(const Exact& rhs) {
        if (!fraction_ && !rhs.fraction_) {
            rational_ += rhs.rational_;
            return *this;
        }
        const Fraction left = as_fraction();
        const Fraction right = rhs.as_fraction();
        if (left.denominator == right.denominator)
            *this = from_fraction({add(left.numerator, right.numerator), left.denominator});
        else
            *this = from_fraction({add(multiply(left.numerator, right.denominator),
                                       multiply(right.numerator, left.denominator)),
                                   multiply(left.denominator, right.denominator)});
        return *this;
    }

    Exact& operator-=(const Exact& rhs) {
        return *this += -rhs;
    }

    Exact& operator*=(const Exact& rhs) {
        if (!fraction_ && !rhs.fraction_) {
            rational_ *= rhs.rational_;
            return *this;
        }
        const Fraction left = as_fraction();
        const Fraction right = rhs.as_fraction();
        *this = from_fraction({multiply(left.numerator, right.numerator),
                               multiply(left.denominator, right.denominator)});
        return *this;
    }

    Exact& operator/=(const Exact& rhs) {
        assert(rhs != 0 && "exact division requires a nonzero divisor");
        if (!fraction_ && !rhs.fraction_) {
            rational_ /= rhs.rational_;
            return *this;
        }
        const Fraction left = as_fraction();
        const Fraction right = rhs.as_fraction();
        *this = from_fraction({multiply(left.numerator, right.denominator),
                               multiply(left.denominator, right.numerator)});
        return *this;
    }

    friend Exact operator+(Exact lhs, const Exact& rhs) {
        return lhs += rhs;
    }
    friend Exact operator-(Exact lhs, const Exact& rhs) {
        return lhs -= rhs;
    }
    friend Exact operator*(Exact lhs, const Exact& rhs) {
        return lhs *= rhs;
    }
    friend Exact operator/(Exact lhs, const Exact& rhs) {
        return lhs /= rhs;
    }
    friend Exact operator-(Exact value) {
        if (value.fraction_) {
            Fraction fraction = *value.fraction_;
            for (auto& [powers, coefficient] : fraction.numerator)
                coefficient = -coefficient;
            value = from_fraction(std::move(fraction));
        } else {
            value.rational_ = -value.rational_;
        }
        return value;
    }

    friend bool operator==(const Exact& lhs, const Exact& rhs) {
        return compare(lhs, rhs) == 0;
    }
    friend bool operator<(const Exact& lhs, const Exact& rhs) {
        return compare(lhs, rhs) < 0;
    }
    friend bool operator>(const Exact& lhs, const Exact& rhs) {
        return rhs < lhs;
    }
    friend bool operator<=(const Exact& lhs, const Exact& rhs) {
        return !(rhs < lhs);
    }
    friend bool operator>=(const Exact& lhs, const Exact& rhs) {
        return !(lhs < rhs);
    }

private:
    Rational rational_;
    std::shared_ptr<const Fraction> fraction_;

    static Polynomial constant(const Rational& value) {
        return value == 0 ? Polynomial{} : Polynomial{{Powers{}, value}};
    }

    Fraction as_fraction() const {
        return fraction_ ? *fraction_ : Fraction{constant(rational_), constant(Rational{1})};
    }

    static int sign(const Polynomial& polynomial) {
        if (polynomial.empty())
            return 0;
        return polynomial.begin()->second > 0 ? 1 : -1;
    }

    static Polynomial add(Polynomial left, const Polynomial& right) {
        for (const auto& [powers, coefficient] : right) {
            auto [entry, inserted] = left.try_emplace(powers, coefficient);
            if (!inserted) {
                entry->second += coefficient;
                if (entry->second == 0)
                    left.erase(entry);
            }
        }
        return left;
    }

    static Polynomial multiply(const Polynomial& left, const Polynomial& right) {
        Polynomial result;
        for (const auto& [left_powers, left_coefficient] : left) {
            for (const auto& [right_powers, right_coefficient] : right) {
                Powers powers{};
                for (std::size_t i = 0; i < powers.size(); ++i)
                    powers[i] = left_powers[i] + right_powers[i];
                auto [entry, inserted] =
                    result.try_emplace(powers, left_coefficient * right_coefficient);
                if (!inserted) {
                    entry->second += left_coefficient * right_coefficient;
                    if (entry->second == 0)
                        result.erase(entry);
                }
            }
        }
        return result;
    }

    static Exact from_fraction(Fraction fraction) {
        assert(!fraction.denominator.empty());
        Exact result;
        if (fraction.numerator.empty())
            return result;
        if (fraction.numerator == fraction.denominator) {
            result.rational_ = 1;
            return result;
        }
        if (fraction.numerator.size() == 1 && fraction.denominator.size() == 1 &&
            fraction.numerator.begin()->first == fraction.denominator.begin()->first) {
            result.rational_ =
                fraction.numerator.begin()->second / fraction.denominator.begin()->second;
            return result;
        }
        const Powers denominator_shift = fraction.denominator.begin()->first;
        const Rational denominator_scale = fraction.denominator.begin()->second;
        auto normalize = [&](const Polynomial& polynomial) {
            Polynomial normalized;
            for (const auto& [powers, coefficient] : polynomial) {
                Powers shifted = powers;
                for (std::size_t i = 0; i < shifted.size(); ++i)
                    shifted[i] -= denominator_shift[i];
                normalized.emplace(shifted, coefficient / denominator_scale);
            }
            return normalized;
        };
        fraction.numerator = normalize(fraction.numerator);
        fraction.denominator = normalize(fraction.denominator);
        assert(fraction.denominator.begin()->first == Powers{} &&
               fraction.denominator.begin()->second == 1 &&
               "normalized exact denominators have leading term one");
        result.fraction_ = std::make_shared<Fraction>(std::move(fraction));
        return result;
    }

    static int compare(const Exact& left, const Exact& right) {
        if (!left.fraction_ && !right.fraction_)
            return left.rational_ == right.rational_ ? 0
                                                     : (left.rational_ < right.rational_ ? -1 : 1);
        if (left.fraction_ == right.fraction_ && left.fraction_)
            return 0;
        const Fraction left_rational = left.fraction_ ? Fraction{} : left.as_fraction();
        const Fraction right_rational = right.fraction_ ? Fraction{} : right.as_fraction();
        const Fraction& lhs = left.fraction_ ? *left.fraction_ : left_rational;
        const Fraction& rhs = right.fraction_ ? *right.fraction_ : right_rational;
        if (lhs.numerator.empty())
            return -sign(rhs.numerator);
        if (rhs.numerator.empty())
            return sign(lhs.numerator);
        const auto& [left_power, left_coefficient] = *lhs.numerator.begin();
        const auto& [right_power, right_coefficient] = *rhs.numerator.begin();
        if (left_power != right_power)
            return left_power < right_power ? sign(lhs.numerator) : -sign(rhs.numerator);
        if (left_coefficient != right_coefficient)
            return left_coefficient < right_coefficient ? -1 : 1;
        Polynomial right_numerator = rhs.numerator;
        for (auto& [powers, coefficient] : right_numerator)
            coefficient = -coefficient;
        if (lhs.denominator == rhs.denominator)
            return sign(add(lhs.numerator, right_numerator)) * sign(lhs.denominator);
        const Polynomial difference = add(multiply(lhs.numerator, rhs.denominator),
                                          multiply(right_numerator, lhs.denominator));
        return sign(difference) * sign(lhs.denominator) * sign(rhs.denominator);
    }
};

inline Exact exact_midpoint(const Exact& a, const Exact& b) {
    return (a + b) / 2;
}

}
