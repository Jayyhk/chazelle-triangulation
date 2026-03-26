#pragma once

/// [C91 §2]: Y-only symbolic perturbation (SoS).
///
/// "Following the tradition of visibility algorithms, we begin by
/// restricting their applicability to polygons where no two distinct
/// vertices have the same y-coordinate. Of course, the standard
/// excuse still works: we can easily get around this assumption by
/// applying the symbolic perturbation techniques of [10] and [31]."
///
/// Edelsbrunner [10] perturbs coordinate j of point i by ε^(2^i·δ+j).
/// Since ε is infinitesimally small, ε^(small exponent) > ε^(large exponent).
/// Lower index i → smaller exponent → larger perturbation → higher perturbed y.
///
/// For y-only comparison (j fixed), this simplifies to:
///   lower index → higher perturbed y
/// Equivalently:
///   larger tag → lower perturbed y
///
/// No actual perturbation is performed — we just use the tag as a
/// tie-breaker in comparisons.

#include "point.h"

#include <cstddef>
#include <limits>

namespace chazelle {

static constexpr std::size_t SOS_NONE =
    (std::numeric_limits<std::size_t>::max)();

/// A symbolic y-coordinate: numeric y plus an SoS tie-break tag.
struct SymbolicY {
    double y = 0.0;
    std::size_t tag = SOS_NONE;
};

// ── SymbolicY comparators ───────────────────────────────────────

/// True if a's perturbed y is strictly less than b's.
/// When y ties, larger tag → lower perturbed y (per [10]).
inline bool symbolic_y_less(SymbolicY a, SymbolicY b) noexcept {
    if (a.y < b.y) return true;
    if (b.y < a.y) return false;
    return a.tag > b.tag;
}

/// True if a's perturbed y is less than or equal to b's.
inline bool symbolic_y_leq(SymbolicY a, SymbolicY b) noexcept {
    return !symbolic_y_less(b, a);
}

/// True if a's perturbed y is strictly greater than b's.
inline bool symbolic_y_greater(SymbolicY a, SymbolicY b) noexcept {
    return symbolic_y_less(b, a);
}

/// True if a's perturbed y is greater than or equal to b's.
inline bool symbolic_y_geq(SymbolicY a, SymbolicY b) noexcept {
    return !symbolic_y_less(a, b);
}

/// True if a and b represent the identical symbolic y-coordinate.
inline bool symbolic_y_equal(SymbolicY a, SymbolicY b) noexcept {
    return a.y == b.y && a.tag == b.tag;
}

/// -1 if a < b, +1 if a > b, 0 if identical.
inline int symbolic_y_compare(SymbolicY a, SymbolicY b) noexcept {
    if (symbolic_y_less(a, b)) return -1;
    if (symbolic_y_less(b, a)) return  1;
    return 0;
}

// ── Point-level helpers ─────────────────────────────────────────

/// Build a SymbolicY from a Point (tag = vertex index).
inline SymbolicY symbolic_y_of(const Point& p) noexcept {
    return SymbolicY{p.y, p.index};
}

/// True if point a's perturbed y is strictly below point b's.
inline bool point_y_below(const Point& a, const Point& b) noexcept {
    return symbolic_y_less(symbolic_y_of(a), symbolic_y_of(b));
}

/// True if point a's perturbed y is strictly above point b's.
inline bool point_y_above(const Point& a, const Point& b) noexcept {
    return symbolic_y_less(symbolic_y_of(b), symbolic_y_of(a));
}

/// -1 if a below b, +1 if a above b, 0 if same vertex.
inline int point_y_order(const Point& a, const Point& b) noexcept {
    return symbolic_y_compare(symbolic_y_of(a), symbolic_y_of(b));
}

// ── Local y-extremum classification ([C91 §2.1]) ───────────────
//
// [C91 §2.1]: "Each vertex of C that is not a local extremum in the
// y-direction gives rise to two companion vertices in ∂C...
// But what about local extrema? For each such vertex of C, if it is
// not one of the two endpoints, we create a total of four vertices
// in ∂C."
//
// These tests determine which case applies (Fig. 2.2):
//   Case 1 (non-extremum): 2 ∂C vertices
//   Case 2 (y-extremum, non-endpoint): 4 ∂C vertices + NLC
//   Case 3 (endpoint): 2 ∂C vertices (duplicates)

/// True if curr is a local y-minimum (strictly below both neighbors).
inline bool is_local_y_minimum(const Point& prev, const Point& curr,
                               const Point& next) noexcept {
    return point_y_below(curr, prev) && point_y_below(curr, next);
}

/// True if curr is a local y-maximum (strictly above both neighbors).
inline bool is_local_y_maximum(const Point& prev, const Point& curr,
                               const Point& next) noexcept {
    return point_y_above(curr, prev) && point_y_above(curr, next);
}

/// True if curr is a local y-extremum (minimum or maximum).
inline bool is_local_y_extremum(const Point& prev, const Point& curr,
                                const Point& next) noexcept {
    return is_local_y_minimum(prev, curr, next) ||
           is_local_y_maximum(prev, curr, next);
}

} // namespace chazelle
