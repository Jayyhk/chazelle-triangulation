#pragma once

// [C91 §2.1]: a simple (nonclosed) polygonal curve C with vertices
// v1, ..., vn and its directed edge table.

#include "point.h"
#include "edge.h"
#include "perturbation.h"

#include <cassert>
#include <cstddef>
#include <vector>

namespace chazelle {

class Polygon {
public:
    // Edge table is built from the vertex sequence.
    explicit Polygon(std::vector<Point> vertices);

    std::size_t num_vertices() const noexcept { return vertices_.size(); }
    std::size_t num_edges()    const noexcept { return edges_.size(); }

    const Point& vertex(std::size_t i) const noexcept {
        assert(i < vertices_.size());
        return vertices_[i];
    }

    const Edge& edge(std::size_t i) const noexcept {
        assert(i < edges_.size());
        return edges_[i];
    }

    // [§2.1 Fig 2.2.3]: one of C's two endpoints.
    bool is_endpoint(std::size_t vertex_index) const noexcept {
        assert(vertex_index < vertices_.size() && "§2.1: invalid vertex");
        return vertex_index == 0
            || vertex_index == vertices_.size() - 1;
    }

    // [§2.1 Fig 2.2.2]: local y-extremum under SoS.  Endpoints are
    // case 3, not case 2 — never extrema.
    bool is_y_extremum(std::size_t vertex_index) const noexcept {
        assert(vertex_index < vertices_.size() && "§2.1: invalid vertex");
        if (is_endpoint(vertex_index)) return false;
        return is_local_y_extremum(
            vertices_[vertex_index - 1],
            vertices_[vertex_index],
            vertices_[vertex_index + 1]);
    }

    // [§2.2]: count of nonnull-length edges in [lo, hi].  O(1) via prefix
    // sums; used by region_weight ("max nonnull edges in any of its arcs").
    std::size_t count_nonnull_edges(std::size_t lo,
                                     std::size_t hi) const noexcept;

private:
    std::vector<Point> vertices_;
    std::vector<Edge>  edges_;
    // nonnull_prefix_[i] = number of nonnull edges in [0, i).
    std::vector<std::size_t> nonnull_prefix_;

    void build_edges();
    void build_nonnull_prefix();
};

} // namespace chazelle
