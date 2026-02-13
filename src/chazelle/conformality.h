#pragma once

/// Conformality Restoration — Stage 2 of the three-stage merge (§3.2).
///
/// After fusion, the submap may have regions with > 4 arcs (tree nodes
/// with degree > 4).  For each such region R, Lemma 3.3 guarantees
/// two non-consecutive arcs Aᵢ, Aⱼ with a vertex of ∂C on Aᵢ visible
/// to Aⱼ.  Find such a visible pair via the arc-cutting oracle and
/// binary search through the tree decomposition of a prior-grade
/// canonical submap (Lemma 2.4).  Insert the discovered chord,
/// splitting R into two regions each with fewer arcs.
///
/// This procedure is called from BOTH the up-phase merge AND the
/// down-phase refinement (§4.2).

#include "visibility/submap.h"
#include "chazelle/grade_storage.h"
#include "geometry/polygon.h"

#include <cstddef>

namespace chazelle {

/// Restore conformality: ensure every tree node has degree ≤ 4.
///
/// @param submap       The submap to make conformal (modified in place).
/// @param storage      Grade storage for arc-cutting oracle lookups.
/// @param polygon      The input polygon.
/// @param granularity  Current granularity level.
void restore_conformality(Submap& submap,
                          const GradeStorage& storage,
                          const Polygon& polygon,
                          std::size_t granularity);

} // namespace chazelle
