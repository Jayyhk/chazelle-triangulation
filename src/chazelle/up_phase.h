#pragma once

/// Up-phase — §4.1 of Chazelle 1991.
///
/// Process grades λ = 0, 1, …, p.
///   - Grade 0: trivial single-edge canonical submaps.
///   - Grade λ > 0: merge two grade-(λ−1) canonical submaps using the
///     three-stage merge.
///
/// All canonical submaps and ray-shooting structures are stored in
/// GradeStorage for use by the down-phase and arc-cutting oracle.

#include "chazelle/grade_storage.h"
#include "geometry/polygon.h"

namespace chazelle {

/// Run the up-phase.
///
/// @param polygon  The input polygon (padded to 2^p + 1 vertices).
/// @param storage  Output: filled with canonical submaps for all grades.
void up_phase(const Polygon& polygon, GradeStorage& storage);

} // namespace chazelle
