#pragma once

/// Common constants and sentinel values used throughout the codebase.

#include <cmath>
#include <cstddef>
#include <limits>

namespace chazelle {

/// Sentinel value meaning "no valid index."
inline constexpr std::size_t NONE = std::numeric_limits<std::size_t>::max();

/// Granularity parameter β = 1/5 (Chazelle §2.3).
inline constexpr double BETA = 1.0 / 5.0;

/// Compute granularity γ = 2^⌈βλ⌉ = 2^⌈λ/5⌉ for a given grade λ.
/// Uses integer arithmetic to avoid floating-point rounding issues.
inline std::size_t compute_granularity(std::size_t grade) {
    std::size_t exponent = (grade + 4) / 5;  // ⌈grade/5⌉
    return std::size_t(1) << exponent;
}

} // namespace chazelle
