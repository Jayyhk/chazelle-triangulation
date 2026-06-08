#pragma once

// [C91 §2 tex 47]: Y-only Simulation of Simplicity.
//
// Paper restricts to polygons with distinct y-coordinates, then waves at
// Edelsbrunner [10] / Yap [31] for the perturbation lift.  [10] perturbs
// coord j of point i by ε^(2^i·δ+j); with ε infinitesimal this means
// lower index → larger perturbation → higher perturbed y.  Equivalently
// for y-only comparisons: larger tag → lower perturbed y.
//
// We perform no actual perturbation — just use the tag as the tie-breaker.

#include "point.h"

#include <cstddef>
#include <limits>

namespace chazelle {

static constexpr std::size_t SOS_NONE =
    (std::numeric_limits<std::size_t>::max)();

// Numeric y + SoS tie-break tag.
struct SymbolicY {
    double y = 0.0;
    std::size_t tag = SOS_NONE;
};

// ── SymbolicY comparators ───────────────────────────────────────

// Y-tied: larger tag ⟹ lower perturbed y (per [10]).
inline bool symbolic_y_less(SymbolicY a, SymbolicY b) noexcept {
    if (a.y < b.y) return true;
    if (b.y < a.y) return false;
    return a.tag > b.tag;
}

inline bool symbolic_y_leq(SymbolicY a, SymbolicY b) noexcept {
    return !symbolic_y_less(b, a);
}

inline bool symbolic_y_greater(SymbolicY a, SymbolicY b) noexcept {
    return symbolic_y_less(b, a);
}

inline bool symbolic_y_geq(SymbolicY a, SymbolicY b) noexcept {
    return !symbolic_y_less(a, b);
}

inline bool symbolic_y_equal(SymbolicY a, SymbolicY b) noexcept {
    return a.y == b.y && a.tag == b.tag;
}

inline int symbolic_y_compare(SymbolicY a, SymbolicY b) noexcept {
    if (symbolic_y_less(a, b)) return -1;
    if (symbolic_y_less(b, a)) return  1;
    return 0;
}

// ── Point-level helpers ─────────────────────────────────────────

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

// ── Local y-extremum classification ─────────────────────────────
//
// [C91 §2.1] Fig 2.2: a vertex falls into one of three cases —
//   1. non-extremum     → 2 ∂C companions
//   2. extremum non-endpoint → 4 ∂C companions + NLC
//   3. endpoint         → 2 ∂C companions (duplicates)

inline bool is_local_y_minimum(const Point& prev, const Point& curr,
                               const Point& next) noexcept {
    return point_y_below(curr, prev) && point_y_below(curr, next);
}

inline bool is_local_y_maximum(const Point& prev, const Point& curr,
                               const Point& next) noexcept {
    return point_y_above(curr, prev) && point_y_above(curr, next);
}

inline bool is_local_y_extremum(const Point& prev, const Point& curr,
                                const Point& next) noexcept {
    return is_local_y_minimum(prev, curr, next) ||
           is_local_y_maximum(prev, curr, next);
}

} // namespace chazelle
