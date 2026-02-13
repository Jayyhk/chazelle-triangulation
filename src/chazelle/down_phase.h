#pragma once

/// Down-Phase — §4.2 of Chazelle 1991.
///
/// Input:
///   - The canonical submap of V(P) from the top grade (λ = p) of the up-phase.
///   - The entire GradeStorage array (all canonical submaps and ray-shooting
///     structures from every chain at every grade).
///
/// Algorithm:
///   Process grades λ = p, p−1, …, 1.  At each grade, refine the current
///   submap by filling in details within each region:
///     1. For each region R, enumerate its arcs (≤ 4 by conformality).
///     2. Symbolic tilting: exit chords bounding R are temporarily treated
///        as curve edges.
///     3. Use arc-cutting oracle to decompose arcs into O(log γ) prior-grade
///        chains.
///     4. Retrieve canonical submaps for these chains from grade_storage.
///     5. Merge into a coarser submap of V(R*).
///     6. Extract exit chords within R; refine the global submap.
///     7. Recurse with reduced granularity from λ to ⌈βλ⌉.
///
/// Output: the complete visibility map V(P) as a normal-form Submap with
/// all chords retained (every region is a single trapezoid).

#include "chazelle/grade_storage.h"
#include "geometry/polygon.h"
#include "visibility/submap.h"

namespace chazelle {

/// Run the down-phase.
///
/// @param polygon  The input polygon.
/// @param storage  Grade storage from the up-phase (mutable — oracles may
///                 be rebound during re-merge).
/// @return         The complete visibility map V(P).
Submap down_phase(const Polygon& polygon, GradeStorage& storage);

} // namespace chazelle
