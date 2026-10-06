#pragma once

#include "bounded_regions.h"
#include "up_phase.h"

namespace chazelle::animation {

struct DownPhaseResult {
    Submap visibility_map;
    RegionCompletionWork completion_work;
    std::size_t refinement_rounds = 0;
    std::size_t region_boundaries = 0;
};

Submap refine_visibility_submap(const UpPhase& up_phase, const Polygon& curve, const Submap& submap,
                                std::size_t grade);

DownPhaseResult compute_visibility_map(const UpPhase& up_phase, std::size_t grade,
                                       std::size_t chain_index);

DownPhaseResult compute_visibility_map(const UpPhase& up_phase);

}
