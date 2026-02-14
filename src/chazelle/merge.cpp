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

    // §4.1 / Lemma 4.1: "we can trivially reset the granularity of
    // each S_i to γ (Section 3.3)."
    // The children have granularity 2^⌈β(λ-1)⌉; the target is
    // γ = 2^⌈βλ⌉ ≥ that.  Resetting to the coarser γ may remove
    // chords that are now redundant at the new granularity.
    // We work on copies so the stored originals remain intact for
    // the down-phase's arc-cutting lookups.
    Submap sub1 = s1.submap;
    Submap sub2 = s2.submap;
    // Propagate chain vertex range so merge_arcs_at_vertex respects
    // the §2.2 C-vertex guard during granularity enforcement.
    sub1.set_chain_info(s1.start_vertex, s1.end_vertex, &polygon);
    sub2.set_chain_info(s2.start_vertex, s2.end_vertex, &polygon);
    enforce_granularity(sub1, gamma);
    enforce_granularity(sub2, gamma);

    // Build oracles on the (possibly reduced) copies.
    RayShootingOracle ora1, ora2;
    ora1.build(sub1, polygon, gamma);
    ora2.build(sub2, polygon, gamma);

    // The junction vertex is the shared vertex between C₁ and C₂.
    // C₁ = chain (2j) at grade λ-1, C₂ = chain (2j+1) at grade λ-1.
    // Junction = last vertex of C₁ = first vertex of C₂.
    std::size_t junction_vertex = s1.end_vertex;

    // Stage 1: Fusion.
    Submap fused = fuse(sub1, ora1, sub2, ora2,
                        polygon, junction_vertex);

    // Stage 2: Conformality restoration.
    restore_conformality(fused, storage, polygon, gamma);

    // Stage 3: Granularity enforcement (§3.3).
    enforce_granularity(fused, gamma);

    // §2.3: Put the submap in normal form.  The arc-sequence table may
    // have been left unsorted by split_arc_at_vertex (conformality) and
    // merge_arcs_at_vertex (granularity).  Sorting restores condition
    // (iii) so that double_identify's O(log m) binary search is sound.
    fused.normalize();

    // Compute chain vertex range.
    auto range = polygon.chain_range(grade, chain_index);

    // Package result.
    CanonicalSubmap result;
    result.submap = std::move(fused);
    result.grade = grade;
    result.chain_index = chain_index;
    result.start_vertex = range.start_vertex;
    result.end_vertex = range.end_vertex;

    // Store chain vertex range and polygon pointer for §2.2 / §2.4.
    result.submap.set_chain_info(range.start_vertex, range.end_vertex,
                                 &polygon);

    // Build ray-shooting oracle for the new canonical submap.
    result.oracle.build(result.submap, polygon, gamma);

    return result;
}

} // namespace chazelle
