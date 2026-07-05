// src/polygon/polygon.cpp

#include "polygon.h"

#include <cassert>
#include <cmath>
#include <unordered_set>

namespace chazelle {

Polygon::Polygon(std::vector<Point> vertices) {
    assert(vertices.size() >= 2 && "[C91 §2.1]: curve needs ≥ 2 vertices");

    // [C91 §2 SoS / Edelsbrunner [10]]: vertex indices must be distinct.
    // This also implies "nonclosed" ([C91 §2.1]): front.index != back.index
    // ⟹ symbolically distinct, regardless of geometric coords.
    // [C91 §2 tex 66]: spherical-plane reference point at (0, +∞) — the
    // curve must avoid that point (Jordan-curve orientation argument).
    {
        std::unordered_set<std::size_t> seen;
        for (const auto& v : vertices) {
            assert(seen.insert(v.index).second &&
                   "[C91 §2 (SoS)]: vertex indices must be distinct");
            assert(!(v.x == 0.0 && std::isinf(v.y) && v.y > 0.0) &&
                   "[C91 §2 tex 66]: curve must avoid the reference point (0, +∞)");
        }
        // [C91 §2.4 tex 133]: every curve is a contiguous subchain of the
        // input P (whose vertices carry consecutive SoS indices v₁..vₙ),
        // so a curve's indices are consecutive from vertices[0].index.
        // Callers rely on this to translate SoS tags → table positions.
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            assert(vertices[i].index == vertices[0].index + i &&
                   "[C91 §2.4 tex 133]: curve must be a contiguous "
                   "subchain of P (consecutive SoS indices)");
        }
    }

    vertices_ = std::move(vertices);
    build_edges();
    build_nonnull_prefix();
    find_y_extremes();
}

void Polygon::find_y_extremes() {
    // [C91 §3.4 tex 306] + [C91 §2 tex 47]: global y-extreme vertices
    // under the perturbed (SoS) order; one O(n) pass at input setup.
    max_y_vertex_ = 0;
    min_y_vertex_ = 0;
    for (std::size_t i = 1; i < vertices_.size(); ++i) {
        if (point_y_above(vertices_[i], vertices_[max_y_vertex_]))
            max_y_vertex_ = i;
        if (point_y_below(vertices_[i], vertices_[min_y_vertex_]))
            min_y_vertex_ = i;
    }
}

void Polygon::build_edges() {
    const std::size_t n = vertices_.size();
    edges_.clear();
    edges_.reserve(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        edges_.push_back(Edge{
            .start_idx = i,
            .end_idx   = i + 1,
        });
    }
}

void Polygon::build_nonnull_prefix() {
    // [C91 §2.2]: an edge has nonzero length iff its endpoints differ in
    // at least one coordinate.
    const std::size_t m = edges_.size();
    nonnull_prefix_.resize(m + 1, 0);
    for (std::size_t i = 0; i < m; ++i) {
        const auto& s = vertices_[edges_[i].start_idx];
        const auto& e = vertices_[edges_[i].end_idx];
        nonnull_prefix_[i + 1] = nonnull_prefix_[i]
            + ((s.x != e.x || s.y != e.y) ? 1 : 0);
    }
}

std::size_t
Polygon::count_nonnull_edges(std::size_t lo,
                              std::size_t hi) const noexcept {
    assert(lo <= hi && "[C91 §2.2]: lo must not exceed hi");
    assert(hi < edges_.size() && "[C91 §2.2]: hi must be within bounds");
    return nonnull_prefix_[hi + 1] - nonnull_prefix_[lo];
}

} // namespace chazelle
