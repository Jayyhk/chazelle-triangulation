#pragma once

/// Symbolic y-perturbation and core geometric predicates.
///
/// Chazelle §2 assumes no two vertices share the same y-coordinate.
/// Rather than perturbing actual coordinates, we use lexicographic
/// comparison (y, x, vertex-index) so that ties are broken consistently
/// and deterministically.  All geometric predicates that depend on
/// vertical ordering must route through these functions.

#include "point.h"

namespace chazelle {

// ── Perturbed y-comparison ──────────────────────────────────────────

/// Lexicographic order on (y, x, index).
/// Returns  < 0  if a is "below" b,
///          = 0  if a and b are the same vertex,
///          > 0  if a is "above" b.
inline int perturbed_y_compare(const Point& a, const Point& b) noexcept {
    if (a.y < b.y) return -1;
    if (a.y > b.y) return  1;
    if (a.x < b.x) return -1;
    if (a.x > b.x) return  1;
    if (a.index < b.index) return -1;
    if (a.index > b.index) return  1;
    return 0;
}

/// Strict "below" test under symbolic perturbation.
inline bool perturbed_y_less(const Point& a, const Point& b) noexcept {
    return perturbed_y_compare(a, b) < 0;
}

/// Functor for use with STL ordered containers / algorithms.
struct PerturbedYOrder {
    constexpr bool operator()(const Point& a, const Point& b) const noexcept {
        if (a.y != b.y) return a.y < b.y;
        if (a.x != b.x) return a.x < b.x;
        return a.index < b.index;
    }
};

// ── Local y-extremum tests ──────────────────────────────────────────

/// True if `curr` is a local y-minimum (both neighbors are above it).
inline bool is_local_y_minimum(const Point& prev, const Point& curr,
                               const Point& next) noexcept {
    return perturbed_y_compare(curr, prev) < 0
        && perturbed_y_compare(curr, next) < 0;
}

/// True if `curr` is a local y-maximum (both neighbors are below it).
inline bool is_local_y_maximum(const Point& prev, const Point& curr,
                               const Point& next) noexcept {
    return perturbed_y_compare(curr, prev) > 0
        && perturbed_y_compare(curr, next) > 0;
}

/// True if `curr` is either a local y-minimum or y-maximum.
inline bool is_local_y_extremum(const Point& prev, const Point& curr,
                                const Point& next) noexcept {
    return is_local_y_minimum(prev, curr, next)
        || is_local_y_maximum(prev, curr, next);
}

// ── Orientation predicate ───────────────────────────────────────────

/// Returns the (unscaled) signed area of triangle (a, b, c):
///   > 0  ⟹  counter-clockwise (left turn)
///   = 0  ⟹  collinear
///   < 0  ⟹  clockwise (right turn)
inline double orient2d(const Point& a, const Point& b,
                       const Point& c) noexcept {
    return (b.x - a.x) * (c.y - a.y)
         - (b.y - a.y) * (c.x - a.x);
}

/// True if the turn a → b → c is strictly convex (counter-clockwise).
inline bool is_convex(const Point& a, const Point& b,
                      const Point& c) noexcept {
    return orient2d(a, b, c) > 0.0;
}

/// True if the turn a → b → c is strictly reflex (clockwise).
inline bool is_reflex(const Point& a, const Point& b,
                      const Point& c) noexcept {
    return orient2d(a, b, c) < 0.0;
}

// ── Horizontal ray intersection ─────────────────────────────────────

/// Given an edge from `start` to `end` and a horizontal ray shooting
/// rightward from y-coordinate `ray_y`, compute the x-coordinate where
/// the ray intersects the edge's supporting line.
///
/// Precondition: the edge is not horizontal (start.y ≠ end.y).
inline double horizontal_ray_x_intercept(const Point& start,
                                         const Point& end,
                                         double ray_y) noexcept {
    double dy = end.y - start.y;
    double t  = (ray_y - start.y) / dy;
    return start.x + t * (end.x - start.x);
}

} // namespace chazelle
