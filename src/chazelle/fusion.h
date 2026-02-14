#pragma once

/// Fusion — Stage 1 of the three-stage merge (§3.1).
///
/// Given two canonical submaps S₁ (of chain C₁) and S₂ (of chain C₂),
/// discover all new chords of V(C₁ ∪ C₂) that cross between C₁ and C₂.
///
/// Two symmetric passes:
///   Pass 1: Walk S₁'s exit-chord endpoints, query S₂ via ray-shooting.
///   Pass 2: Walk S₂'s exit-chord endpoints, query S₁ via ray-shooting.
///
/// Output: the fusion submap S — a valid submap of V(C) but not yet
/// conformal or appropriately granular.

#include "visibility/submap.h"
#include "oracles/ray_shooting.h"
#include "chazelle/grade_storage.h"
#include "geometry/polygon.h"

#include <cstddef>

namespace chazelle {

/// Fuse two canonical submaps into one.
///
/// @param s1         Canonical submap of chain C₁.
/// @param oracle1    Ray-shooting oracle for S₁.
/// @param s2         Canonical submap of chain C₂.
/// @param oracle2    Ray-shooting oracle for S₂.
/// @param polygon    The input polygon.
/// @return           The fused submap (not yet conformal or granular).
/// @param junction_vertex  Index of the vertex C₁ ∩ C₂ (the last vertex
///                        of C₁ = first vertex of C₂).  Per §3.1, the
///                        companion vertices from duplicating this vertex
///                        are always included in the fusion enumeration.
Submap fuse(const Submap& s1, const RayShootingOracle& oracle1,
            const Submap& s2, const RayShootingOracle& oracle2,
            const Polygon& polygon,
            std::size_t junction_vertex,
            const GradeStorage& storage);

} // namespace chazelle
