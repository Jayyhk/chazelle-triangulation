#pragma once

/// [C91 §3.1]: Fusion of two submaps.
///
/// "By symmetry, we may limit our discussion to the problem of
/// fusing S₁ into S₂, that is, determining the points of ∂C that
/// are seen by the endpoints of the exit chords of S₁ and by the
/// companion vertices resulting from the duplication of C₁ ∩ C₂."

#include "../polygon/polygon.h"
#include "../polygon/perturbation.h"
#include "../submap/submap.h"
#include "oracle.h"

#include <cassert>
#include <cstddef>
#include <vector>

namespace chazelle {

/// [C91 §3.1]: A vertex in the fusion sequence.
///
/// The sequence is a₀, a₁, ..., aₘ, a_{m+1} where:
///   - a₀ = companion vertex at junction (start of clockwise tour)
///   - a₁, ..., aₘ = exit chord endpoints of S₁ in clockwise ∂C₁ order
///   - a_{m+1} = companion vertex at junction (end of clockwise tour)
struct FusionVertex {
    SymbolicY y;               ///< The vertex's symbolic y-coordinate.
    std::size_t edge;          ///< Edge of P containing this vertex.
    Side side;                 ///< Which side of ∂C₁.
    std::size_t chord_idx;     ///< Chord index in S₁ (NONE for companions).
    bool is_left_endpoint;     ///< True if this is the chord's LEFT endpoint.
    bool is_companion;         ///< True for a₀ or a_{m+1}.
};

/// [C91 §3.1]: Build the fusion vertex sequence for fusing S₁ into S₂.
///
/// "Let a_{m+1} and a₀ be the companion vertices, as they appear next
/// to each other clockwise around ∂C₁, resulting from the duplication
/// of the vertex C₁ ∩ C₂ in ∂C₁.  Let a₁, a₂, ..., aₘ be the
/// canonical vertex enumeration of S₁."
///
/// "A clockwise tour around ∂C₁ that begins at a₀ thus ends at a_{m+1}."
///
/// The clockwise ∂C₁ order starting from the junction (c_end of C₁) is:
///   RIGHT side (c_end → c_start, descending edges), then
///   LEFT side (c_start → c_end, ascending edges).
///
/// @param S       The submap S₁.
/// @param C       The curve C₁.
/// @return The fusion vertex sequence a₀, a₁, ..., aₘ, a_{m+1}.
std::vector<FusionVertex> build_fusion_sequence(const Submap& S,
                                                 const Polygon& C);

} // namespace chazelle
