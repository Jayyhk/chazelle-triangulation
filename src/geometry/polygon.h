#pragma once

/// The read-only "input table" for Chazelle's algorithm.
///
/// A simple polygon is stored as an open polygonal curve
/// C = (v₀, v₁, …, v_{n−1}) with n−1 directed edges.
/// If the input is a closed polygon, one edge is conceptually
/// removed ("punctured") to form the open curve.
///
/// The vertex count is padded to n = 2^p + 1 (by inserting
/// collinear vertices on the last edge) so that grades
/// decompose cleanly into dyadic chains.

#include "point.h"
#include "edge.h"

#include <cstddef>
#include <span>
#include <vector>

namespace chazelle {

class Polygon {
public:
    /// Construct from a vertex list in boundary order.
    /// For a closed simple polygon the closing edge (v_{n-1} → v₀)
    /// is conceptually removed; the curve runs v₀ … v_{n-1}.
    explicit Polygon(std::vector<Point> vertices);

    // ── Sizes ───────────────────────────────────────────────────

    /// Number of vertices after padding: n = 2^p + 1.
    std::size_t num_vertices() const noexcept { return vertices_.size(); }

    /// Number of edges after padding: n − 1 = 2^p.
    std::size_t num_edges() const noexcept { return edges_.size(); }

    /// Number of grades: p = ⌈log₂(n−1)⌉ (before padding this is
    /// rounded up; after padding n−1 is exactly 2^p).
    std::size_t num_grades() const noexcept { return num_grades_; }

    /// Original vertex count before padding.
    std::size_t original_size() const noexcept { return original_size_; }

    // ── Element access ──────────────────────────────────────────

    const Point& vertex(std::size_t i) const { return vertices_[i]; }
    const Edge&  edge(std::size_t i)   const { return edges_[i]; }

    std::span<const Point> vertices() const noexcept { return vertices_; }
    std::span<const Edge>  edges()    const noexcept { return edges_; }

    // ── Chain / grade queries ───────────────────────────────────

    /// A half-open vertex range [start_vertex, end_vertex] (both inclusive)
    /// defining a chain at a given grade.
    struct ChainRange {
        std::size_t start_vertex; ///< First vertex index (inclusive).
        std::size_t end_vertex;   ///< Last  vertex index (inclusive).

        std::size_t num_vertices() const noexcept {
            return end_vertex - start_vertex + 1;
        }
        std::size_t num_edges() const noexcept {
            return end_vertex - start_vertex;
        }
    };

    /// Number of chains at grade λ: 2^(p − λ).
    std::size_t num_chains(std::size_t grade) const noexcept;

    /// Vertex range for chain `chain_index` at grade `grade`.
    /// Chain j at grade λ covers vertices [j·2^λ, (j+1)·2^λ].
    ChainRange chain_range(std::size_t grade,
                           std::size_t chain_index) const noexcept;

    // ── Vertex classification ───────────────────────────────────

    /// True if `vertex_index` is the first or last vertex of C.
    bool is_endpoint(std::size_t vertex_index) const noexcept {
        return vertex_index == 0
            || vertex_index == vertices_.size() - 1;
    }

    /// True if the vertex at `vertex_index` is a local y-extremum
    /// (under symbolic perturbation).  Endpoints are never extrema.
    bool is_y_extremum(std::size_t vertex_index) const noexcept;

    /// Compute the x-coordinate of the point on edge `edge_idx` at
    /// height `y`.  Used for chord endpoint positions: a chord stores
    /// (edge, y) and this gives the exact x.  If the edge is
    /// horizontal, returns the midpoint x.
    double edge_x_at_y(std::size_t edge_idx, double y) const noexcept;

private:
    std::vector<Point> vertices_;
    std::vector<Edge>  edges_;
    std::size_t num_grades_    = 0; ///< p
    std::size_t original_size_ = 0; ///< vertex count before padding

    void pad_to_power_of_two_plus_one();
    void build_edges();
};

} // namespace chazelle
