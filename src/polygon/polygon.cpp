/// src/polygon/polygon.cpp

#include "polygon.h"

#include <cassert>

namespace chazelle {

Polygon::Polygon(std::vector<Point> vertices) {
    assert(vertices.size() >= 2 &&
           "§2.1: a polygonal curve must have at least 2 vertices");
    vertices_ = std::move(vertices);
    build_edges();
}

void Polygon::build_edges() {
    const std::size_t n = vertices_.size();
    edges_.clear();
    edges_.reserve(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        edges_.push_back(Edge{
            .index     = i,
            .start_idx = i,
            .end_idx   = i + 1,
        });
    }
}

} // namespace chazelle
