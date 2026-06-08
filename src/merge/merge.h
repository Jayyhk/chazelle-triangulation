#pragma once

// [C91 §3]: Merge two conformal γᵢ-granular submaps S₁, S₂ of V(C₁), V(C₂)
// (γ₁ ≤ γ₂; C₁ ∩ C₂ a single vertex of P) into a normal-form γ-granular
// conformal submap of V(C), with γ ≥ γ₂.

#include "../polygon/polygon.h"
#include "../submap/submap.h"
#include "oracle.h"

#include <cassert>
#include <cstddef>

namespace chazelle {

struct MergeInput {
    // C₁ ∩ C₂ must be a single vertex of P (last of C₁ = first of C₂).
    const Polygon* C1 = nullptr;
    const Polygon* C2 = nullptr;

    // [§3]: Sᵢ are merge INPUTS — read by all stages, never written.
    const Submap* S1 = nullptr;
    const Submap* S2 = nullptr;

    std::size_t gamma1 = 0;     // γ₁ ≤ γ₂.
    std::size_t gamma2 = 0;
    std::size_t gamma  = 0;     // Target γ ≥ γ₂.

    // [§3.0 tex 166–170]: Per-submap oracles (§3.4 Lemma 3.6 builds
    // ray-shooting per S).  tex 220: f(γ₁) for S₁ arcs, f(γ₂) for S₂.
    const RayShootingOracle* ray_shooter_1 = nullptr;
    const RayShootingOracle* ray_shooter_2 = nullptr;
    const ArcCuttingOracle*  arc_cutter_1  = nullptr;
    const ArcCuttingOracle*  arc_cutter_2  = nullptr;
};

struct MergeResult {
    Polygon C;                  // Merged curve C = C₁ ∪ C₂.
    Submap S;                   // Normal-form γ-granular conformal V(C).
    explicit MergeResult(Polygon c) : C(std::move(c)) {}
};

// [C91 §3]: Validate every §3-preamble precondition.
inline void assert_merge_preconditions(const MergeInput& in) {
    assert(in.C1 != nullptr && in.C2 != nullptr &&
           "§3: merge requires two polygonal curves");
    assert(in.S1 != nullptr && in.S2 != nullptr &&
           "§3: merge requires two submaps");

    assert(in.C1->num_vertices() >= 2 && "§3: C₁ must have ≥ 2 vertices");
    assert(in.C2->num_vertices() >= 2 && "§3: C₂ must have ≥ 2 vertices");

    // [§3 tex 160]: "C₁ ∩ C₂ is a vertex of P."  Vertex identity under SoS
    // is the .index field alone (perturbation.h:34–47; polygon.cpp:21–25
    // asserts uniqueness).  Adding .x/.y equality would over-specify.
    {
        const auto& c1_last = in.C1->vertex(in.C1->num_vertices() - 1);
        const auto& c2_first = in.C2->vertex(0);
        assert(c1_last.index == c2_first.index &&
               "§3 tex 160: C₁ ∩ C₂ must be a vertex of P "
               "(last of C₁ = first of C₂ by SoS .index)");
    }

    // [§3 tex 166]: "each Sᵢ is given in normal form."  check_invariants
    // verifies §2.4(i)–(iii) + key_y monotonicity (tex 144) + edge_count
    // cache (tex 106).
    in.S1->check_invariants(*in.C1);
    in.S2->check_invariants(*in.C2);

    // [§3]: "Sᵢ a γᵢ-granular conformal submap of V(Cᵢ)."
    assert(in.S1->is_conformal() && "§3: S₁ must be conformal");
    assert(!in.S1->tree_decomposition().empty() &&
           "§2.4(iv): S₁ TD must be available");
    assert(in.S2->is_conformal() && "§3: S₂ must be conformal");
    assert(!in.S2->tree_decomposition().empty() &&
           "§2.4(iv): S₂ TD must be available");
    assert(in.S1->is_granular(in.gamma1, *in.C1) &&
           "§3: S₁ must be γ₁-granular");
    assert(in.S2->is_granular(in.gamma2, *in.C2) &&
           "§3: S₂ must be γ₂-granular");

    assert(in.ray_shooter_1 && in.ray_shooter_2 &&
           in.arc_cutter_1 && in.arc_cutter_2 &&
           "§3.0 tex 166–170: all four per-submap oracles required");

    assert(in.gamma1 <= in.gamma2 && "§3: γ₁ ≤ γ₂");
    assert(in.gamma  >= in.gamma2 && "§3: target γ ≥ γ₂");
}

// [C91 §3]: Three-stage merge: (§3.1) fuse via discovered chords,
// (§3.2) restore conformality, (§3.3) enforce γ-granularity.
// TODO: implement.
MergeResult merge(const MergeInput& in);

} // namespace chazelle
