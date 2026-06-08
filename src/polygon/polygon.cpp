// src/polygon/polygon.cpp

#include "polygon.h"

#include <cassert>
#include <unordered_set>

namespace chazelle {

Polygon::Polygon(std::vector<Point> vertices) {
    assert(vertices.size() >= 2 && "§2.1: curve needs ≥ 2 vertices");

    // [§2 SoS / Edelsbrunner [10]]: vertex indices must be distinct.
    // This also implies "nonclosed" (§2.1): front.index != back.index
    // ⟹ symbolically distinct, regardless of geometric coords.
    {
        std::unordered_set<std::size_t> seen;
        for (const auto& v : vertices) {
            assert(seen.insert(v.index).second &&
                   "§2 (SoS): vertex indices must be distinct");
        }
    }

    vertices_ = std::move(vertices);
    build_edges();
    build_nonnull_prefix();
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
    // [§2.2]: an edge has nonzero length iff its endpoints differ in
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
    assert(lo <= hi && "§2.2: lo must not exceed hi");
    assert(hi < edges_.size() && "§2.2: hi must be within bounds");
    return nonnull_prefix_[hi + 1] - nonnull_prefix_[lo];
}

} // namespace chazelle
