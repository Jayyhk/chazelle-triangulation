#pragma once

// Common constants used throughout the codebase.

#include <cstddef>
#include <limits>

namespace chazelle {

// Sentinel value meaning "no valid index."
inline constexpr std::size_t NONE = std::numeric_limits<std::size_t>::max();

} // namespace chazelle
