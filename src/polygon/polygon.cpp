/// src/polygon/polygon.cpp

#include "polygon.h"

#include <cassert>
#include <unordered_set>

namespace chazelle {

Polygon::Polygon(std::vector<Point> vertices) {
    assert(vertices.size() >= 2 &&
           "§2.1: a polygonal curve must have at least 2 vertices");

    // [C91 §2.1]: "simple (nonclosed) polygonal curve."  Under SoS,
    // "nonclosed" means front and back are distinct vertices, which
    // is implied by the distinct-index check below (front.index !=
    // back.index ⟹ symbolically distinct vertices, regardless of
    // whether they share geometric coordinates).
    //
    // [C91 §2] / [10]: SoS requires distinct vertex indices.
    {
        std::unordered_set<std::size_t> seen;
        for (const auto& v : vertices) {
            assert(seen.insert(v.index).second &&
                   "§2 (SoS): all vertex indices must be distinct");
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
    // [C91 §2.2]: "the maximum number of nonnull length edges."
    // An edge has nonzero geometric length iff its two endpoints
    // differ in at least one coordinate.
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
