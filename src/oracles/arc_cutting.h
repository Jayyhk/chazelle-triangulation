#pragma once

/// Arc-Cutting Oracle — Chazelle §3.4.
///
/// Given a subarc of ∂C, decompose it into O(log γ) chains from prior
/// grades whose canonical submaps have already been computed in the up-phase.
///
/// This is a thin utility function (~10 lines of logic): binary decomposition
/// of the arc's vertex range into dyadic intervals, each matching a chain
/// from some grade λ' < λ.

#include "geometry/polygon.h"
#include "visibility/arc_structure.h"

#include <cstddef>
#include <vector>

namespace chazelle {

/// One piece of an arc decomposition.
struct ArcPiece {
    std::size_t grade;       ///< Grade of the chain this piece belongs to.
    std::size_t chain_index; ///< Chain index within that grade.
    std::size_t start_vertex; ///< Start vertex of this piece.
    std::size_t end_vertex;   ///< End vertex of this piece.
};

/// Decompose an arc's vertex range [start, end) into O(log γ) dyadic
/// intervals.  Each interval corresponds to a chain at some grade λ' < λ
/// whose canonical submap is available in grade_storage.
///
/// O(log γ) total — each piece is computed in O(1) via bit manipulation,
/// and there are at most O(log γ) pieces in the greedy decomposition.
///
/// @param start_vertex  First vertex index of the arc (inclusive).
/// @param end_vertex    Past-the-end vertex index of the arc (exclusive).
/// @param polygon       The polygon (for chain decomposition info).
/// @return              List of ArcPieces covering [start, end).
inline std::vector<ArcPiece>
cut_arc(std::size_t start_vertex, std::size_t end_vertex,
        const Polygon& polygon) {
    std::vector<ArcPiece> pieces;

    if (start_vertex >= end_vertex) return pieces;

    // Binary decomposition into dyadic intervals.
    // A dyadic interval [a, a + 2^k] is a chain at grade k with index a / 2^k.
    // We greedily find the largest dyadic interval starting at `pos`.
    std::size_t pos = start_vertex;
    while (pos < end_vertex) {
        std::size_t remaining = end_vertex - pos;

        // Find the largest k such that:
        //   1. 2^k divides pos (alignment constraint).
        //   2. 2^k ≤ remaining (doesn't overshoot end_vertex).
        //   3. k ≤ num_grades (valid grade).
        // Use bit manipulation for O(1):
        //   max_k_alignment = number of trailing zeros of pos (if pos > 0),
        //                     or a large number (if pos == 0, all k work).
        //   max_k_remaining = floor(log2(remaining)).
        std::size_t max_k_align = (pos == 0)
            ? static_cast<std::size_t>(64)
            : static_cast<std::size_t>(__builtin_ctzll(
                  static_cast<unsigned long long>(pos)));
        std::size_t max_k_fit = (remaining == 0)
            ? static_cast<std::size_t>(0)
            : static_cast<std::size_t>(63 - __builtin_clzll(
                  static_cast<unsigned long long>(remaining)));
        std::size_t k = std::min(max_k_align,
                        std::min(max_k_fit, polygon.num_grades()));
        std::size_t len = std::size_t(1) << k;

        ArcPiece piece;
        piece.grade = k;
        piece.chain_index = pos / len;
        piece.start_vertex = pos;
        piece.end_vertex = pos + len;
        pieces.push_back(piece);

        pos += len;
    }

    return pieces;
}

} // namespace chazelle
