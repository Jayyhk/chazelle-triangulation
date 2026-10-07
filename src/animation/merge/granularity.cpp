#include "granularity.h"
#include "../trace.h"

#include <algorithm>
#include <cassert>
#include <vector>

namespace chazelle::animation {

namespace {

std::vector<std::size_t> collect_region_weights(const Submap& submap) {
    std::vector<std::size_t> weights(submap.num_nodes(), 0);
    for (std::size_t arc_index = 0; arc_index < submap.num_arcs(); ++arc_index) {
        const Arc& arc = submap.arc(arc_index);
        if (arc.dead)
            continue;
        assert(arc.region_node < weights.size());
        weights[arc.region_node] = std::max(weights[arc.region_node], arc.edge_count);
    }
    return weights;
}

bool has_endpoint_below_degree_three(const Submap& submap, std::size_t chord_index) {
    const Chord& chord = submap.chord(chord_index);
    return submap.node(chord.region[0]).degree() < 3 || submap.node(chord.region[1]).degree() < 3;
}

}

void enforce_granularity(Submap& submap, const Polygon& curve, std::size_t granularity) {
    assert(granularity >= 1 && "[C91 §2.3]: granularity parameter must be positive");

    assert(submap.is_conformal() && "[C91 §3.3 tex 276]: S must be conformal after §3.2");
    submap.assert_tree_property();

    std::vector<std::size_t> region_weights = collect_region_weights(submap);

#ifndef NDEBUG
    for (std::size_t r = 0; r < submap.num_nodes(); ++r) {
        if (submap.node(r).dead)
            continue;
        assert(region_weights[r] <= granularity &&
               "[C91 §3.3 tex 276]: S must be γ-semigranular on entry");
    }
#endif

    auto contraction_weight = [&](std::size_t chord_index) -> std::size_t {
        const Chord& chord = submap.chord(chord_index);
        return submap.simulated_contraction_weight(
            chord_index, curve, region_weights[chord.region[0]], region_weights[chord.region[1]]);
    };

    std::vector<std::size_t> pending_chords;
    pending_chords.reserve(submap.num_chords());
    for (std::size_t chord_index = submap.num_chords(); chord_index-- > 0;)
        pending_chords.push_back(chord_index);

    std::size_t removed_chord_count = 0;
    while (!pending_chords.empty()) {
        const std::size_t chord_index = pending_chords.back();
        pending_chords.pop_back();
        if (submap.chord(chord_index).dead)
            continue;
        const bool eligible = has_endpoint_below_degree_three(submap, chord_index);
        const std::size_t merged_weight = eligible ? contraction_weight(chord_index) : NONE;
        if (auto* trace = AnimationTrace::current()) {
            const auto& chord = submap.chord(chord_index);
            trace->record("granularity_test",
                          {{"owner", trace->map_id(submap)},
                           {"chord", chord_index},
                           {"first", chord.region[0]},
                           {"second", chord.region[1]},
                           {"first_degree", submap.node(chord.region[0]).degree()},
                           {"second_degree", submap.node(chord.region[1]).degree()},
                           {"eligible", eligible},
                           {"weight", merged_weight},
                           {"limit", granularity},
                           {"accepted", eligible && merged_weight <= granularity}});
        }
        if (!eligible)
            continue;
        if (merged_weight > granularity)
            continue;

        const std::size_t surviving_region = submap.remove_chord(chord_index, curve);
        ++removed_chord_count;
        assert(removed_chord_count <= submap.num_chords() && "each chord is removed once");

        assert(surviving_region < region_weights.size());
        region_weights[surviving_region] = merged_weight;

        assert(submap.node(surviving_region).degree() <= 4 &&
               "[C91 §3.3 tex 276]: contraction must preserve conformality");

        for (std::size_t incident_chord : submap.node(surviving_region).incident_chords) {
            assert(!submap.chord(incident_chord).dead && "incident_chords holds live chords only");
            pending_chords.push_back(incident_chord);
        }
    }

#ifndef NDEBUG
    assert(submap.is_conformal() && "[C91 §3.3 tex 276]: removals keep the submap conformal");
    submap.assert_tree_property();

    {
        const std::vector<std::size_t> recomputed_weights = collect_region_weights(submap);
        for (std::size_t r = 0; r < submap.num_nodes(); ++r) {
            if (submap.node(r).dead)
                continue;
            assert(recomputed_weights[r] == region_weights[r] &&
                   "[C91 §3.3]: maintained weights must match the table");

            assert(recomputed_weights[r] <= granularity &&
                   "[C91 §3.3 tex 276]: S stays γ-semigranular");
        }
    }

    for (std::size_t ci = 0; ci < submap.num_chords(); ++ci) {
        if (submap.chord(ci).dead)
            continue;
        if (!has_endpoint_below_degree_three(submap, ci))
            continue;
        assert(contraction_weight(ci) > granularity &&
               "[C91 §2.3 tex 121]: contracting any edge incident upon a "
               "degree-<3 node must exceed γ (Lemma 3.5 postcondition)");
    }
#endif
}

}
