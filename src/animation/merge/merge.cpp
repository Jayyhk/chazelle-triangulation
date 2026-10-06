#include "merge.h"
#include "../trace.h"
#include "conformality.h"
#include "fusion.h"
#include "granularity.h"

namespace chazelle::animation {

void assert_merge_preconditions([[maybe_unused]] const MergeInput& input) {
#ifndef NDEBUG
    assert(input.first_curve != nullptr && input.second_curve != nullptr &&
           "[C91 §3]: merge requires two polygonal curves");
    assert(input.first_submap != nullptr && input.second_submap != nullptr &&
           "[C91 §3]: merge requires two submaps");

    assert(input.first_curve->num_vertices() >= 2 && "[C91 §3]: C₁ must have ≥ 2 vertices");
    assert(input.second_curve->num_vertices() >= 2 && "[C91 §3]: C₂ must have ≥ 2 vertices");

    {
        const auto& first_curve_end =
            input.first_curve->vertex(input.first_curve->num_vertices() - 1);
        const auto& second_curve_start = input.second_curve->vertex(0);
        assert(first_curve_end.index == second_curve_start.index &&
               "[C91 §3 tex 160]: C₁ ∩ C₂ must be a vertex of P "
               "(last of C₁ = first of C₂ by SoS .index)");
    }

    input.first_submap->check_invariants(*input.first_curve);
    input.second_submap->check_invariants(*input.second_curve);

    assert(input.first_submap->is_conformal() && "[C91 §3]: S₁ must be conformal");
    assert(!input.first_submap->tree_decomposition().empty() &&
           "[C91 §2.4(iv)]: S₁ tree decomposition must be available");
    assert(input.second_submap->is_conformal() && "[C91 §3]: S₂ must be conformal");
    assert(!input.second_submap->tree_decomposition().empty() &&
           "[C91 §2.4(iv)]: S₂ tree decomposition must be available");
    assert(input.first_submap->is_granular(input.first_granularity, *input.first_curve) &&
           "[C91 §3]: S₁ must be γ₁-granular");
    assert(input.second_submap->is_granular(input.second_granularity, *input.second_curve) &&
           "[C91 §3]: S₂ must be γ₂-granular");

    assert(input.first_ray_shooter && input.second_ray_shooter && input.first_arc_cutter &&
           input.second_arc_cutter &&
           "[C91 §3.0 tex 166–170]: all four per-submap oracles required");

    assert(input.first_piece_count_bound >= 1 && input.second_piece_count_bound >= 1 &&
           input.first_piece_granularity_bound >= 1 && input.second_piece_granularity_bound >= 1 &&
           "[C91 §3.0(ii) tex 170]: arc-cutter bounds g(γᵢ), h(γᵢ) required");

    assert(input.first_granularity <= input.second_granularity && "[C91 §3]: γ₁ ≤ γ₂");
    assert(input.granularity >= input.second_granularity && "[C91 §3]: target γ ≥ γ₂");
#endif
}

static void fuse_both_directions(Submap& submap, const MergeInput& input, const Polygon& curve) {
    FusionState first_pass;
    first_pass.junction_at_end = true;
    fuse_submaps(first_pass, *input.first_submap, *input.first_curve, *input.second_submap,
                 *input.second_curve, *input.first_ray_shooter, *input.second_ray_shooter);

    FusionState second_pass;
    second_pass.junction_at_end = false;
    fuse_submaps(second_pass, *input.second_submap, *input.second_curve, *input.first_submap,
                 *input.first_curve, *input.second_ray_shooter, *input.first_ray_shooter);

    rebuild_submap(submap, curve, *input.first_submap, *input.first_curve, *input.second_submap,
                   *input.second_curve, first_pass, second_pass);
}

static void restore_conformality(Submap& submap, const MergeInput& input, const Polygon& curve) {
    const ConformalityOracles oracles{
        .first_submap = input.first_submap,
        .second_submap = input.second_submap,
        .first_curve = input.first_curve,
        .second_curve = input.second_curve,
        .first_ray_shooter = input.first_ray_shooter,
        .second_ray_shooter = input.second_ray_shooter,
        .first_arc_cutter = input.first_arc_cutter,
        .second_arc_cutter = input.second_arc_cutter,
        .first_piece_count_bound = input.first_piece_count_bound,
        .second_piece_count_bound = input.second_piece_count_bound,
        .first_piece_granularity_bound = input.first_piece_granularity_bound,
        .second_piece_granularity_bound = input.second_piece_granularity_bound,
    };
    restore_conformality(submap, curve, oracles);
}

static void maintain_granularity(Submap& submap, const MergeInput& input, const Polygon& curve) {
    enforce_granularity(submap, curve, input.granularity);

    submap.normalize(curve);

    assert(submap.is_conformal() && "[C91 Lemma 3.5 tex 279]: merge output must be conformal");
    assert(submap.is_granular(input.granularity, curve) &&
           "[C91 Lemma 3.5 tex 279]: merge output must be γ-granular");
    assert(!submap.tree_decomposition().empty() &&
           "[C91 §2.4(iv) tex 139]: normal form includes the tree "
           "decomposition");
}

MergeResult merge(const MergeInput& input) {
    assert_merge_preconditions(input);
    if (auto* trace = AnimationTrace::current()) {
        trace->merge(*input.first_curve, *input.second_curve, input.granularity);
    }

    MergeResult result(Polygon(*input.first_curve, *input.second_curve));

    fuse_both_directions(result.submap, input, result.curve);
    restore_conformality(result.submap, input, result.curve);
    maintain_granularity(result.submap, input, result.curve);
    animation_checkpoint("merge_end");

    return result;
}

}
