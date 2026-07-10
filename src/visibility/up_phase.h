#pragma once

// [C91 §4.1 tex 336–357]: The Up-Phase.
//
// "Given a curve C consisting of m contiguous edges of P, we say that a
// submap of V(C) is canonical if it is 2^{⌈β⌈log m⌉⌉}-granular,
// conformal, and represented in normal form.  Note that a canonical
// submap for a chain in grade λ is 2^{⌈βλ⌉}-granular.  For
// λ = 0, 1, ..., p in that order, we process grade λ, which means:
//   (i)  We compute a canonical submap of V(C) for each chain C in that
//        grade.
//   (ii) We preprocess each canonical submap for ray-shooting along the
//        lines of Lemma 3.6, setting γ to the value 2^{⌈βλ⌉}."
// (tex 338–343.)  "This work can be done naively for the early grades"
// (tex 345); past them, each chain's canonical submap comes from
// Lemma 4.1 (tex 347–356): partition into ≤ 2λ chains of lower grades,
// reset their granularities to γ = 2^{⌈βλ⌉} (§3.3), and merge two-by-two
// along a perfectly balanced binary tree, with the ray-shooting and
// arc-cutting oracles realized from the lower grades' structures
// (tex 350–353) — which is also where [C91 §3.4 tex 284]'s postponed
// arc-cutter and [C91 §3.0(i) tex 169]'s per-piece ray-shooter live.
//
// Total up-phase time: "processing grade λ takes O(n·2^{−λ/76}) time,
// therefore processing all p + 1 grades, and thereby completing the
// up-phase, takes linear time" (tex 356–357).

#include "chain.h"
#include "naive_visibility.h"
#include "../merge/merge.h"
#include "../merge/oracle.h"
#include "../merge/ray_shooting.h"
#include "../polygon/polygon.h"
#include "../submap/submap.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace chazelle {

class UpPhase {
public:
    // Runs the whole up-phase (tex 340–345) on the input polygonal
    // curve: pads to n = 2^p + 1 ([C91 §4 tex 316]), builds the chain
    // grid, and processes grades 0, 1, ..., p in that order.  O(n).
    explicit UpPhase(std::vector<Point> vertices);

    const GradedCurve& graded() const noexcept { return graded_; }

    // [C91 §4.1 tex 338/342]: the canonical granularity of grade λ,
    // γ = 2^{⌈βλ⌉}, in exact integer arithmetic.
    static std::size_t grade_gamma(std::size_t lambda) noexcept {
        return std::size_t{1} << ceil_beta(lambda);
    }

    // [C91 §4.1 tex 345 / tex 352]: Lemma 4.1 assumes ⌈βλ⌉ < λ ("which
    // is true for λ large enough") — false exactly for λ ≤ 1 with
    // β = 1/5, so grades 0 and 1 are the naively-computed early grades.
    static constexpr std::size_t NAIVE_MAX_GRADE = 1;

    // (i) of tex 341: the canonical submap of V(chain(grade, i)).
    const Submap& chain_submap(std::size_t grade, std::size_t i) const {
        assert(grade < submaps_.size() && i < submaps_[grade].size() &&
               "[C91 §4.1 tex 345]: grade must be processed already");
        return submaps_[grade][i];
    }

    // (ii) of tex 342: the Lemma 3.6 ray-shooting structure of
    // chain(grade, i), preprocessed with γ = 2^{⌈βλ⌉}.
    const RayShootingStructure& chain_structure(std::size_t grade,
                                                std::size_t i) const {
        assert(grade < structures_.size() &&
               i < structures_[grade].size() &&
               "[C91 §4.1 tex 345]: grade must be processed already");
        return *structures_[grade][i];
    }

    // [C91 Lemma 4.1 tex 347]: "given any portion D of P of the form
    // v_a, ..., v_b, where 2^{λ−1} < b − a ≤ 2^λ, we can compute a
    // canonical submap of V(D) in time proportional to
    // λ²(log λ)2^{λ(1−β/3+4β²/3)}."  `a` and `b` are 0-based table
    // vertex indices of the padded curve P (the paper's v_{a} .. v_{b}
    // are 1-based); requires ⌈βλ⌉ < λ and all grades < λ processed.
    struct PortionResult {
        Polygon curve;
        Submap submap;
    };
    PortionResult compute_canonical_portion(std::size_t a,
                                            std::size_t b) const;

private:
    GradedCurve graded_;
    // Indexed [grade][chain]: the canonical submaps (i) and their
    // Lemma 3.6 shooting structures (ii) — "the data structures left
    // behind during the top-phase (the submaps for the chains and their
    // ray-shooting structures)" that the down-phase reuses
    // ([C91 §4 tex 323]).
    std::vector<std::vector<Submap>> submaps_;
    std::vector<std::vector<std::unique_ptr<RayShootingStructure>>>
        structures_;
};

// ── [C91 §4.1 tex 350–353]: the up-phase oracles ─────────────────
//
// "Any subarc α' ⊆ α can be subdivided into a constant number of
// contiguous pieces (with no double-backing) whose corresponding
// portions of P consist of single line segments (at most two of them)
// and vertex-to-vertex pieces of P, each with at most 2^{⌈βλ⌉} edges.
// Each of these pieces can be partitioned into a collection of O(λ)
// chains in grades at most ⌈βλ⌉.  ...  Thus, to shoot a ray toward α',
// we shoot toward each of the O(λ) subarcs of its decomposition and
// determine the closest hit (if any)."  (tex 350–351.)
//
// Both oracles share that decomposition.  `Si`/`Ci` is one INPUT of a
// Lemma 4.1 merge at level λ (an intermediate union of partition
// chains); subarcs arrive in Cᵢ-local edge indices and pieces are
// reported the same way; the chain grid is P-global, aligned via
// Polygon::table_offset.
//
// Padding note ([C91 §4 tex 316] + [C91 §2.2 tex 106]): the paper's
// edge counts are NONNULL counts, so a piece with ≤ 2^{⌈βλ⌉} (nonnull)
// edges may span arbitrarily many null padding edges.  The cover
// therefore admits ALL-NULL chains of any processed grade — their
// canonical submaps are chordless (total weight 0) and h-granular for
// every h, so tex 353's h(γ) ≤ 2^{⌈β⌈βλ⌉⌉} is preserved — while every
// chain containing a nonnull edge is kept in grade ≤ ⌈βλ⌉ exactly as
// tex 352 states.

class UpPhaseRayShooter final : public RayShootingOracle {
public:
    UpPhaseRayShooter(const UpPhase& up, const Submap& Si,
                      const Polygon& Ci, std::size_t lambda);

    // [C91 §3.0(i) tex 169] realized per tex 351: the first contact of
    // the ray with ᾱ' (the lexicographic (wrapped, distance) minimum
    // over the pieces' first contacts) is reported iff it lies ON α'.
    // Time O(f(γ)) with f(γ) = λ·2^{β²λ/3 + 2βλ/3} (tex 352).
    RayHit shoot(Point p, Side direction, std::size_t arc_idx,
                 const Subarc& target,
                 double source_x_offset = SOURCE_OFFSET_NONE) const override;

private:
    const UpPhase* up_;
    const Submap* Si_;
    const Polygon* Ci_;
    std::size_t lambda_;
};

class UpPhaseArcCutter final : public ArcCuttingOracle {
public:
    UpPhaseArcCutter(const UpPhase& up, const Submap& Si,
                     const Polygon& Ci, std::size_t lambda);

    // [C91 §3.0(ii) tex 170] realized per tex 352–353 — the arc-cutter
    // whose discussion [C91 §3.4 tex 284] postponed to the up-phase:
    // ≤ 2 single-edge boundary pieces at α''s endpoints plus O(λ)
    // chains whose canonical submaps and curves are the up-phase's own
    // stores, so g(γ) = O(λ) and h(γ) ≤ 2^{⌈β⌈βλ⌉⌉} (tex 353).
    std::vector<ArcPiece> cut(std::size_t arc_idx,
                              const Subarc& target) const override;

    // The concrete piece-count bound g(γ) = O(λ) advertised for level
    // λ (tex 353); MergeInput.g_gamma for Lemma 4.1's merges.
    static std::size_t g_bound(std::size_t lambda) noexcept {
        // ≤ 2 boundary pieces + per leg (≤ 3): the ≤ 2(λ+1)-chain
        // dyadic cover with the one padding-boundary chain split into
        // ≤ λ extra all-null chains + 2 final parts.
        return 2 + 3 * (2 * (lambda + 1) + lambda + 2);
    }

    // The piece-granularity bound h(γ) ≤ 2^{⌈β⌈βλ⌉⌉} of tex 353;
    // MergeInput.h_gamma for Lemma 4.1's merges.
    static std::size_t h_bound(std::size_t lambda) noexcept {
        return UpPhase::grade_gamma(ceil_beta(lambda));
    }

private:
    const UpPhase* up_;
    const Submap* Si_;
    const Polygon* Ci_;
    std::size_t lambda_;
};

} // namespace chazelle
