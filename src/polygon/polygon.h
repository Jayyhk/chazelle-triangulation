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

    // [C91 §2.1 Fig 2.2.3]: one of C's two endpoints.
    bool is_endpoint(std::size_t vertex_index) const noexcept {
        assert(vertex_index < vertices_.size() && "[C91 §2.1]: invalid vertex");
        return vertex_index == 0
            || vertex_index == vertices_.size() - 1;
    }

    // [C91 §2.1 Fig 2.2.2]: local y-extremum under SoS.  Endpoints are
    // case 3, not case 2 — never extrema.
    bool is_y_extremum(std::size_t vertex_index) const noexcept {
        assert(vertex_index < vertices_.size() && "[C91 §2.1]: invalid vertex");
        if (is_endpoint(vertex_index)) return false;
        return is_local_y_extremum(
            vertices_[vertex_index - 1],
            vertices_[vertex_index],
            vertices_[vertex_index + 1]);
    }

    // [C91 §2.2]: count of nonnull-length edges in [lo, hi].  O(1) via prefix
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

// [C91 §2 tex 47]: SoS-aware interpolation parameter on an edge.
//
// Returns t ∈ [0, 1] such that the (perturbed) horizontal line at
// `target_y` crosses edge `edge_idx` at the convex combination
// (1-t)·vs + t·ve.  When `target_y`'s SymbolicY matches an endpoint's
// (same raw y and tag), the perturbed crossing IS that endpoint —
// returns 0.0 or 1.0 exactly.  Otherwise raw y's strictly bracket
// target_y (else SoS would have forced a tag-match), and standard
// linear interpolation is exact.
//
// Required by the paper's "distinct y" relaxation under SoS: with
// repeated raw y-coords a horizontal edge has `vs.y == ve.y` and
// only the symbolic short-circuit produces a defined answer.
inline double edge_t_at_y(const Polygon& C, std::size_t edge_idx,
                          SymbolicY target_y) {
    const auto& e = C.edge(edge_idx);
    const Point& vs = C.vertex(e.start_idx);
    const Point& ve = C.vertex(e.end_idx);
    if (symbolic_y_equal(target_y, symbolic_y_of(vs))) return 0.0;
    if (symbolic_y_equal(target_y, symbolic_y_of(ve))) return 1.0;
    // [C91 §2 tex 47]: with no tag-match, the perturbed horizontal line
    // crosses the OPEN edge, i.e. target_y lies symbolically strictly
    // between the endpoint ys.  (A raw-y tie with an endpoint is fine —
    // the crossing is infinitesimally beside that vertex and the linear
    // interpolation below returns the exact limit point.)
    assert(((symbolic_y_less(symbolic_y_of(vs), target_y) &&
             symbolic_y_less(target_y, symbolic_y_of(ve))) ||
            (symbolic_y_less(symbolic_y_of(ve), target_y) &&
             symbolic_y_less(target_y, symbolic_y_of(vs)))) &&
           "[C91 §2 tex 47]: query must lie strictly inside the edge's "
           "perturbed y-range (no crossing otherwise)");
    // Raw-horizontal edge: its endpoints are consecutive vertices of C
    // carrying consecutive SoS indices, so no third tag can fall
    // symbolically strictly between them — unreachable given the assert
    // above.
    assert(vs.y != ve.y &&
           "[C91 §2 tex 47]: horizontal edge requires SoS tag-match at "
           "one endpoint; strictly-between target is unreachable");
    return (target_y.y - vs.y) / (ve.y - vs.y);
}

// [C91 §2 tex 47]: SoS-aware x of where the perturbed horizontal line
// at `target_y` crosses edge `edge_idx`.  Companion to `edge_t_at_y`.
inline double edge_x_at_y(const Polygon& C, std::size_t edge_idx,
                          SymbolicY target_y) {
    const auto& e = C.edge(edge_idx);
    const Point& vs = C.vertex(e.start_idx);
    const Point& ve = C.vertex(e.end_idx);
    double t = edge_t_at_y(C, edge_idx, target_y);
    return vs.x + t * (ve.x - vs.x);
}

} // namespace chazelle
