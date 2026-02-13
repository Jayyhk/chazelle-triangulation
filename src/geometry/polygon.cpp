#include "geometry/polygon.h"
#include "geometry/perturbation.h"

#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace chazelle {

// ── Construction ────────────────────────────────────────────────────

Polygon::Polygon(std::vector<Point> vertices)
    : original_size_{vertices.size()}
{
    if (vertices.size() < 3) {
        throw std::invalid_argument(
            "Polygon must have at least 3 vertices");
    }

    // Assign vertex indices (identity mapping before padding).
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].index = i;
    }

    vertices_ = std::move(vertices);

    pad_to_power_of_two_plus_one();
    build_edges();
}

// ── Padding ─────────────────────────────────────────────────────────

void Polygon::pad_to_power_of_two_plus_one() {
    // We need n = 2^p + 1 vertices  ⟹  2^p edges.
    // Find smallest p such that 2^p ≥ (current vertex count − 1).
    const std::size_t n = vertices_.size();
    const std::size_t num_edges_needed = n - 1;

    // std::bit_ceil gives the smallest power-of-two ≥ its argument.
    const auto power =
        std::bit_ceil(static_cast<std::uint64_t>(num_edges_needed));

    num_grades_ = static_cast<std::size_t>(
        std::countr_zero(static_cast<std::uint64_t>(power)));

    const std::size_t target_n = static_cast<std::size_t>(power) + 1;

    if (target_n <= n) {
        // Already the right size (or exact power-of-two + 1).
        return;
    }

    // Insert collinear vertices on the last edge
    // (the edge from vertices_[n-2] to vertices_[n-1]).
    const std::size_t to_add = target_n - n;
    const Point  second_last = vertices_[n - 2]; // copy — reserve may invalidate refs
    const Point  last        = vertices_[n - 1]; // copy — we'll modify the vector

    // Remove the old last vertex; we'll re-add it after the new ones.
    vertices_.pop_back();
    vertices_.reserve(target_n);

    for (std::size_t i = 1; i <= to_add; ++i) {
        const double t =
            static_cast<double>(i) / static_cast<double>(to_add + 1);
        Point p;
        p.x     = second_last.x + t * (last.x - second_last.x);
        p.y     = second_last.y + t * (last.y - second_last.y);
        p.index = vertices_.size();
        vertices_.push_back(p);
    }

    // Re-add the original last vertex with its new index.
    Point last_copy  = last;
    last_copy.index  = vertices_.size();
    vertices_.push_back(last_copy);
}

// ── Edge table ──────────────────────────────────────────────────────

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

// ── Chain / grade queries ───────────────────────────────────────────

std::size_t Polygon::num_chains(std::size_t grade) const noexcept {
    assert(grade <= num_grades_);
    // 2^(p − λ)
    return std::size_t{1} << (num_grades_ - grade);
}

Polygon::ChainRange
Polygon::chain_range(std::size_t grade,
                     std::size_t chain_index) const noexcept {
    assert(grade <= num_grades_);
    assert(chain_index < num_chains(grade));

    const std::size_t chain_edges = std::size_t{1} << grade; // 2^λ
    return ChainRange{
        .start_vertex = chain_index * chain_edges,
        .end_vertex   = (chain_index + 1) * chain_edges,
    };
}

// ── Vertex classification ───────────────────────────────────────────

bool Polygon::is_y_extremum(std::size_t vertex_index) const noexcept {
    if (is_endpoint(vertex_index)) {
        return false; // Endpoints of C are not considered extrema.
    }

    return is_local_y_extremum(
        vertices_[vertex_index - 1],
        vertices_[vertex_index],
        vertices_[vertex_index + 1]);
}

double Polygon::edge_x_at_y(std::size_t edge_idx, double y) const noexcept {
    if (edge_idx >= edges_.size()) return 0.0;
    const auto& e = edges_[edge_idx];
    const auto& p1 = vertices_[e.start_idx];
    const auto& p2 = vertices_[e.end_idx];
    double dy = p2.y - p1.y;
    if (std::abs(dy) < 1e-15) {
        // Horizontal edge — return midpoint x.
        return (p1.x + p2.x) * 0.5;
    }
    double t = (y - p1.y) / dy;
    return p1.x + t * (p2.x - p1.x);
}

} // namespace chazelle
