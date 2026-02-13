#pragma once

/// Trapezoid record for the Fournier-Montuno triangulation step.
///
/// Each trapezoid is defined by a top and bottom vertex (the two
/// horizontal scan lines) and left/right bounding polygon edges.
/// Degenerate trapezoids (where top or bottom is a single point)
/// are treated the same way — they're just triangles.

#include <cstddef>
#include <limits>

namespace chazelle {

inline constexpr std::size_t FM_NONE = std::numeric_limits<std::size_t>::max();

struct Trapezoid {
    std::size_t top_vertex    = FM_NONE; ///< Vertex index with higher y.
    std::size_t bottom_vertex = FM_NONE; ///< Vertex index with lower y.
    std::size_t left_edge     = FM_NONE; ///< Left bounding edge (by its top vertex index).
    std::size_t right_edge    = FM_NONE; ///< Right bounding edge (by its top vertex index).
};

} // namespace chazelle
