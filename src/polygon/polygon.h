#pragma once

// [C91 §2.1]: a simple (nonclosed) polygonal curve C with vertices
// v1, ..., vn and its directed edge table.
//
// [C91 §2.4 tex 133]: "Let P be the input polygonal curve (the one
// whose visibility map is sought) and let C be the subchain of P whose
// visibility map (or submap) we wish to represent. ... We assume that
// the edges of P are stored in a table (the input table) in the order
// in which they occur along the boundary of P. ... The input table is
// read-only: it is never to be modified or even copied."  A Polygon is
// therefore a VIEW (offset + length) into one shared immutable input
// table; constructing the merged curve C = C₁ ∪ C₂ ([C91 §3 tex 160])
// is O(1) index arithmetic, never a copy — Lemma 4.1's sublinear
// per-merge budget ([C91 §4.1 tex 347–350]) forbids anything else.

#include "point.h"
#include "edge.h"
#include "perturbation.h"

#include <cassert>
#include <cstddef>
#include <memory>
#include <vector>

namespace chazelle {

class Polygon {
public:
    // Root curve: allocates the input table (done once per input
    // setup; [C91 §2.4 tex 133]).  O(n).
    explicit Polygon(std::vector<Point> vertices);

    // [C91 §3 tex 160]: the merged curve C = C₁ ∪ C₂, where "C₁ ∩ C₂
    // is a vertex of P" — C₁ and C₂ must be contiguous views into the
    // SAME input table sharing exactly that junction vertex.  O(1):
    // the table is never copied ([C91 §2.4 tex 133]) and the y-extreme
    // vertices combine from the two halves.
    Polygon(const Polygon& c1, const Polygon& c2);

    // Contiguous subchain view [first, first + count) of this curve
    // ([C91 §2.4 tex 133]: every curve is a subchain of P).  O(count)
    // for the y-extreme scan — used at input setup, not inside merges.
    Polygon subchain(std::size_t first, std::size_t count) const;

    std::size_t num_vertices() const noexcept { return len_; }
    std::size_t num_edges()    const noexcept { return len_ - 1; }

    const Point& vertex(std::size_t i) const noexcept {
        assert(i < len_);
        return table_->vertices[offset_ + i];
    }

    // Directed edge i runs from local vertex i to i + 1 ([C91 §2.1]).
    Edge edge(std::size_t i) const noexcept {
        assert(i + 1 < len_);
        return Edge{i, i + 1};
    }

    // [C91 §2.1 Fig 2.2.3]: one of C's two endpoints.
    bool is_endpoint(std::size_t vertex_index) const noexcept {
        assert(vertex_index < len_ && "[C91 §2.1]: invalid vertex");
        return vertex_index == 0 || vertex_index == len_ - 1;
    }

    // [C91 §2.1 Fig 2.2.2]: local y-extremum under SoS.  Endpoints are
    // case 3, not case 2 — never extrema.
    bool is_y_extremum(std::size_t vertex_index) const noexcept {
        assert(vertex_index < len_ && "[C91 §2.1]: invalid vertex");
        if (is_endpoint(vertex_index)) return false;
        return is_local_y_extremum(vertex(vertex_index - 1),
                                   vertex(vertex_index),
                                   vertex(vertex_index + 1));
    }

    // [C91 §2.2]: count of nonnull-length edges in [lo, hi].  O(1) via
    // the table-wide prefix sums; used by region_weight ("max nonnull
    // edges in any of its arcs").
    std::size_t count_nonnull_edges(std::size_t lo,
                                    std::size_t hi) const noexcept {
        assert(lo <= hi && "[C91 §2.2]: lo must not exceed hi");
        assert(hi < num_edges() && "[C91 §2.2]: hi must be within bounds");
        return table_->nonnull_prefix[offset_ + hi + 1] -
               table_->nonnull_prefix[offset_ + lo];
    }

    // [C91 §3.4 tex 306]: the ray-shooting structure's vertical line is
    // anchored "to the right of all the vertices of P", and its polar
    // segments are delimited by the curve's global y-extremes; the
    // extreme vertices (under the SoS order of [C91 §2 tex 47]) are
    // scanned once at input setup and combined in O(1) per merge.
    std::size_t max_y_vertex() const noexcept { return max_y_vertex_; }
    std::size_t min_y_vertex() const noexcept { return min_y_vertex_; }

private:
    // The one input table of P: immutable, shared by every subchain
    // view, never copied ([C91 §2.4 tex 133]).
    struct Table {
        std::vector<Point> vertices;
        // nonnull_prefix[i] = number of nonnull edges among the table's
        // edges [0, i).
        std::vector<std::size_t> nonnull_prefix;
    };

    Polygon() = default;                       // internal (subchain)

    void find_y_extremes();                    // O(len_) scan

    std::shared_ptr<const Table> table_;
    std::size_t offset_ = 0;                   // first vertex in table
    std::size_t len_ = 0;                      // number of vertices
    std::size_t max_y_vertex_ = 0;             // local indices
    std::size_t min_y_vertex_ = 0;
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
    const auto e = C.edge(edge_idx);
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
    const auto e = C.edge(edge_idx);
    const Point& vs = C.vertex(e.start_idx);
    const Point& ve = C.vertex(e.end_idx);
    double t = edge_t_at_y(C, edge_idx, target_y);
    return vs.x + t * (ve.x - vs.x);
}

} // namespace chazelle
