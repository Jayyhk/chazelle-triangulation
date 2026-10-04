#pragma once

#include "point.h"

#include <cstddef>
#include <limits>

namespace chazelle {

static constexpr std::size_t SOS_NONE = (std::numeric_limits<std::size_t>::max)();

struct SymbolicY {
    Exact y = 0.0;
    std::size_t tag = SOS_NONE;
};

inline bool symbolic_y_less(const SymbolicY& a, const SymbolicY& b) noexcept {
    if (a.y < b.y)
        return true;
    if (b.y < a.y)
        return false;
    return a.tag > b.tag;
}

inline bool symbolic_y_leq(const SymbolicY& a, const SymbolicY& b) noexcept {
    return !symbolic_y_less(b, a);
}

inline bool symbolic_y_greater(const SymbolicY& a, const SymbolicY& b) noexcept {
    return symbolic_y_less(b, a);
}

inline bool symbolic_y_geq(const SymbolicY& a, const SymbolicY& b) noexcept {
    return !symbolic_y_less(a, b);
}

inline bool symbolic_y_equal(const SymbolicY& a, const SymbolicY& b) noexcept {
    return a.y == b.y && a.tag == b.tag;
}

inline int symbolic_y_compare(const SymbolicY& a, const SymbolicY& b) noexcept {
    if (symbolic_y_less(a, b))
        return -1;
    if (symbolic_y_less(b, a))
        return 1;
    return 0;
}

inline SymbolicY symbolic_y_of(const Point& p) noexcept {
    return SymbolicY{p.y, p.index};
}

inline bool point_y_below(const Point& a, const Point& b) noexcept {
    return symbolic_y_less(symbolic_y_of(a), symbolic_y_of(b));
}

inline bool point_y_above(const Point& a, const Point& b) noexcept {
    return symbolic_y_less(symbolic_y_of(b), symbolic_y_of(a));
}

inline int point_y_order(const Point& a, const Point& b) noexcept {
    return symbolic_y_compare(symbolic_y_of(a), symbolic_y_of(b));
}

inline bool is_local_y_minimum(const Point& prev, const Point& curr, const Point& next) noexcept {
    return point_y_below(curr, prev) && point_y_below(curr, next);
}

inline bool is_local_y_maximum(const Point& prev, const Point& curr, const Point& next) noexcept {
    return point_y_above(curr, prev) && point_y_above(curr, next);
}

inline bool is_local_y_extremum(const Point& prev, const Point& curr, const Point& next) noexcept {
    return is_local_y_minimum(prev, curr, next) || is_local_y_maximum(prev, curr, next);
}

}
