#include "chazelle/merge.h"
#include "chazelle/fusion.h"
#include "chazelle/conformality.h"
#include "chazelle/granularity.h"
#include "common.h"

#include <cassert>
#include <cmath>

namespace chazelle {

CanonicalSubmap merge_submaps(std::size_t grade,
                               std::size_t chain_index,
                               GradeStorage& storage,
                               const Polygon& polygon) {
    // Chain at grade λ with index j is the union of two grade-(λ−1) chains:
    //   C₁ = chain (2j)   at grade λ−1
    //   C₂ = chain (2j+1) at grade λ−1
    assert(grade > 0);
    std::size_t prev_grade = grade - 1;
    std::size_t child1_idx = 2 * chain_index;
    std::size_t child2_idx = 2 * chain_index + 1;

    assert(child1_idx < storage.num_chains(prev_grade));
    assert(child2_idx < storage.num_chains(prev_grade));

    auto& s1 = storage.get(prev_grade, child1_idx);
    auto& s2 = storage.get(prev_grade, child2_idx);

    std::size_t gamma = compute_granularity(grade);

    // Rebuild the oracles for the two children so they point to the
    // submaps' current addresses.  (Moving a CanonicalSubmap into
    // storage invalidates the oracle's internal pointer.)
    std::size_t gamma1 = compute_granularity(prev_grade);
    if (!s1.oracle.is_built()) s1.oracle.build(s1.submap, polygon, gamma1);
    else s1.oracle.rebind_submap(s1.submap);
    if (!s2.oracle.is_built()) s2.oracle.build(s2.submap, polygon, gamma1);
    else s2.oracle.rebind_submap(s2.submap);

    // The junction vertex is the shared vertex between C₁ and C₂.
    // C₁ = chain (2j) at grade λ-1, C₂ = chain (2j+1) at grade λ-1.
    // Junction = last vertex of C₁ = first vertex of C₂.
    std::size_t junction_vertex = s1.end_vertex;

    // Stage 1: Fusion.
    Submap fused = fuse(s1.submap, s1.oracle, s2.submap, s2.oracle,
                        polygon, junction_vertex);

    // Stage 2: Conformality restoration.
    restore_conformality(fused, storage, polygon, gamma);

    // Stage 3: Granularity enforcement.
    // Protect null-length chords: they mark y-extrema and must be
    // preserved for the complete visibility map V(P).
    enforce_granularity(fused, gamma, /*protect_null_length=*/true);

    // Compute chain vertex range.
    auto range = polygon.chain_range(grade, chain_index);

    // Package result.
    CanonicalSubmap result;
    result.submap = std::move(fused);
    result.grade = grade;
    result.chain_index = chain_index;
    result.start_vertex = range.start_vertex;
    result.end_vertex = range.end_vertex;

    // Build ray-shooting oracle for the new canonical submap.
    result.oracle.build(result.submap, polygon, gamma);

    return result;
}

} // namespace chazelle
