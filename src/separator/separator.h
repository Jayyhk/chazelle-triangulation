#pragma once

/// Lipton-Tarjan O(n) planar separator (1979).
///
/// find_separator():  the 10-step algorithm producing (A, B, C) with
///   |C| ≤ 2√(2n),  cost(A) ≤ 2/3,  cost(B) ≤ 2/3.
///
/// iterated_separator():  recursively partition into pieces of size
///   ≤ target_piece_size, returning the separator nodes and pieces.

#include "planar_graph.h"
#include <cstddef>
#include <vector>

namespace chazelle {

/// Result of find_separator: three disjoint vertex sets.
struct SeparatorResult {
    std::vector<std::size_t> A; ///< Vertices on one side.
    std::vector<std::size_t> B; ///< Vertices on the other side.
    std::vector<std::size_t> C; ///< Separator vertices.
};

/// Find a 2/3-balanced separator of the planar graph.
/// Costs are taken from vertex.cost fields (should sum to 1.0 over
/// the vertices to be separated, or the algorithm normalises).
///
/// Implements Lipton & Tarjan 1979, §3, Steps 1–10.
SeparatorResult find_separator(PlanarGraph& graph);

/// Result of iterated_separator: pieces + global separator set.
struct IteratedSeparatorResult {
    /// separator_nodes[i] is true if vertex i is in the separator.
    std::vector<bool> separator_nodes;

    /// Each piece is a list of (non-separator) vertex indices.
    std::vector<std::vector<std::size_t>> pieces;

    /// Which piece each vertex belongs to (NONE if separator).
    std::vector<std::size_t> vertex_piece;
};

/// Recursively apply find_separator until every piece has
/// ≤ `max_piece_size` vertices.
IteratedSeparatorResult iterated_separator(PlanarGraph& graph,
                                           std::size_t max_piece_size);

} // namespace chazelle
