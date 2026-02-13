#pragma once

/// Granularity Enforcement — Stage 3 of the three-stage merge (§3.3).
///
/// After conformality restoration, the submap is conformal and γ-semigranular.
/// This step enforces the maximality condition of γ-granularity:
///
/// For each exit chord, check whether contracting the corresponding tree edge
/// (merging the two adjacent regions) would produce a node with weight ≤ γ.
/// If so, the chord is redundant: remove it.
///
/// This is a single pass over the submap tree. Only removes chords, never
/// inserts any.

#include "visibility/submap.h"

#include <cstddef>

namespace chazelle {

/// Enforce granularity: remove redundant chords so that contracting any
/// edge incident on a degree-<3 node would exceed γ.
///
/// @param submap              The submap to make granular (modified in place).
/// @param gamma               The target granularity.
/// @param protect_null_length When true, null-length chords at y-extrema are
///                            never removed.  Set this during the down-phase
///                            where the goal is the complete V(P); leave false
///                            during the up-phase merge.
/// @param max_chord_idx       Only consider chords with index < this value.
///                            Use NONE (default) to consider all chords.
void enforce_granularity(Submap& submap, std::size_t gamma,
                         bool protect_null_length = false,
                         std::size_t max_chord_idx = NONE);

} // namespace chazelle
