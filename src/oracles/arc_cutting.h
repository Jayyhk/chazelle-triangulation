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
        // Find the largest power of 2 that:
        //   1. Divides pos (so [pos, pos+2^k] is a valid chain boundary).
        //   2. Doesn't exceed end_vertex (pos + 2^k ≤ end_vertex).
        //   3. Doesn't exceed num_grades (valid grade).
        std::size_t remaining = end_vertex - pos;
        std::size_t k = 0;
        std::size_t len = 1;

        while (len * 2 <= remaining &&
               (pos % (len * 2)) == 0 &&
               k + 1 <= polygon.num_grades()) {
            len *= 2;
            ++k;
        }

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
