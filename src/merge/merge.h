#pragma once

/// [C91 §3]: Merging two conformal submaps.
///
/// "Let C₁ and C₂ be two polygonal curves of n₁ and n₂ vertices
/// respectively, whose union C forms a connected vertex-to-vertex
/// piece of the input (simple and nonclosed) polygonal curve P;
/// we assume that C₁ ∩ C₂ is a vertex of P.  Let Sᵢ (i=1,2) be
/// a γᵢ-granular conformal submap of V(Cᵢ), with γ₁ ≤ γ₂.
/// Given any integer γ ≥ γ₂, to merge S₁ and S₂ (where γ is
/// understood) means to compute a normal-form γ-granular conformal
/// submap of V(C)."

#include "../polygon/polygon.h"
#include "../submap/submap.h"
#include "oracle.h"

#include <cassert>
#include <cstddef>

namespace chazelle {

/// [C91 §3]: Input to the merge operation.
struct MergeInput {
    /// The two polygonal curves.  C₁ ∩ C₂ must be a single vertex
    /// of P (the last vertex of C₁ = the first vertex of C₂).
    const Polygon* C1 = nullptr;
    const Polygon* C2 = nullptr;

    /// Normal-form γᵢ-granular conformal submaps of V(Cᵢ).
    Submap* S1 = nullptr;
    Submap* S2 = nullptr;

    /// Granularities.  γ₁ ≤ γ₂.
    std::size_t gamma1 = 0;
    std::size_t gamma2 = 0;

    /// Target granularity.  γ ≥ γ₂.
    std::size_t gamma = 0;

    /// [C91 §3]: Oracle primitives (provided by §3.4 / §4).
    const RayShootingOracle* ray_shooter = nullptr;
    const ArcCuttingOracle* arc_cutter = nullptr;
};

/// [C91 §3]: Result of the merge operation.
struct MergeResult {
    /// The merged curve C = C₁ ∪ C₂.
    Polygon C;

    /// Normal-form γ-granular conformal submap of V(C).
    Submap S;

    explicit MergeResult(Polygon c) : C(std::move(c)) {}
};

/// [C91 §3]: Validate merge preconditions.
///
/// Asserts all invariants from the paper's §3 preamble:
///   - C₁, C₂ are valid polygonal curves
///   - C₁ ∩ C₂ is a vertex (last of C₁ = first of C₂)
///   - S₁ is conformal
///   - S₂ is conformal
///   - γ₁ ≤ γ₂
///   - γ ≥ γ₂
inline void assert_merge_preconditions(const MergeInput& in) {
    assert(in.C1 != nullptr && in.C2 != nullptr &&
           "§3: merge requires two polygonal curves");
    assert(in.S1 != nullptr && in.S2 != nullptr &&
           "§3: merge requires two submaps");

    // [C91 §3]: "whose union C forms a connected vertex-to-vertex
    // piece of the input polygonal curve P; we assume that
    // C₁ ∩ C₂ is a vertex of P."
    // The last vertex of C₁ must equal the first vertex of C₂.
    assert(in.C1->num_vertices() >= 2 &&
           "§3: C₁ must have at least 2 vertices");
    assert(in.C2->num_vertices() >= 2 &&
           "§3: C₂ must have at least 2 vertices");
    {
        const auto& c1_last = in.C1->vertex(in.C1->num_vertices() - 1);
        const auto& c2_first = in.C2->vertex(0);
        assert(c1_last.x == c2_first.x &&
               c1_last.y == c2_first.y &&
               c1_last.index == c2_first.index &&
               "§3: C₁ ∩ C₂ must be a vertex of P "
               "(last vertex of C₁ = first vertex of C₂)");
    }

    // [C91 §3]: "Sᵢ be a γᵢ-granular conformal submap of V(Cᵢ)."
    assert(in.S1->is_conformal() &&
           "§3: S₁ must be conformal");
    assert(in.S2->is_conformal() &&
           "§3: S₂ must be conformal");
    assert(in.S1->is_granular(in.gamma1, *in.C1) &&
           "§3: S₁ must be γ₁-granular");
    assert(in.S2->is_granular(in.gamma2, *in.C2) &&
           "§3: S₂ must be γ₂-granular");

    // [C91 §3]: "we assume that we have at our disposal two primitives."
    assert(in.ray_shooter != nullptr &&
           "§3: merge requires a ray-shooting oracle");
    assert(in.arc_cutter != nullptr &&
           "§3: merge requires an arc-cutting oracle");

    // [C91 §3]: "with γ₁ ≤ γ₂."
    assert(in.gamma1 <= in.gamma2 &&
           "§3: γ₁ must be ≤ γ₂");

    // [C91 §3]: "Given any integer γ ≥ γ₂."
    assert(in.gamma >= in.gamma2 &&
           "§3: target γ must be ≥ γ₂");
}

/// [C91 §3]: Merge S₁ and S₂ to produce a normal-form γ-granular
/// conformal submap of V(C).
///
/// The merge proceeds in three stages (§3.1–§3.3):
///   1. Fusion: discover new chords, create submap S of V(C).
///   2. Conformality: add chords to ensure degree ≤ 4.
///   3. Granularity: remove chords to enforce γ-granularity.
///
/// TODO: (§3.1–§3.3) implement the three stages.
MergeResult merge(const MergeInput& in);

} // namespace chazelle
