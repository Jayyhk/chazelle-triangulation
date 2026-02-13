#pragma once

/// Three-stage merge — combines fusion, conformality, and granularity.
///
/// Given two canonical submaps S₁, S₂ of adjacent chains C₁, C₂:
///   1. Fuse (§3.1): discover cross-chain chords.
///   2. Restore conformality (§3.2): ensure degree ≤ 4.
///   3. Enforce granularity (§3.3): remove redundant chords.
///
/// Result: a canonical (conformal, γ-granular, normal-form) submap of V(C₁ ∪ C₂).

#include "chazelle/grade_storage.h"
#include "geometry/polygon.h"

#include <cstddef>

namespace chazelle {

/// Three-stage merge: fuse → conformality → granularity.
///
/// @param grade        The grade being computed.
/// @param chain_index  The chain index within the grade.
/// @param storage      Grade storage with prior-grade canonical submaps.
/// @param polygon      The input polygon.
/// @return             The merged canonical submap.
CanonicalSubmap merge_submaps(std::size_t grade,
                               std::size_t chain_index,
                               GradeStorage& storage,
                               const Polygon& polygon);

} // namespace chazelle
