#include "chazelle/granularity.h"
#include "geometry/polygon.h"

#include <algorithm>
#include <cassert>
#include <vector>

namespace chazelle {

void enforce_granularity(Submap& submap, std::size_t gamma,
                         std::size_t max_chord_idx) {
    // §3.3: Single-pass chord removal.  Chords need be processed only
    // once since removals cannot make any chord removable if it was not
    // already so before (see paper).
    //
    // A chord is removable if:
    //   1. At least one of its endpoint regions has degree < 3.
    //   2. Merging the two regions would produce weight ≤ γ.

    std::size_t chord_limit = (max_chord_idx != NONE)
                              ? std::min(max_chord_idx, submap.num_chords())
                              : submap.num_chords();

    // §3.3: "Chords need be processed only once since the removals
    // cannot make any chord removable if it was not already so
    // before."  A single pass suffices.
    for (std::size_t ci = 0; ci < chord_limit; ++ci) {
        auto& c = submap.chord(ci);
        if (c.region[0] == NONE || c.region[1] == NONE) continue;



        auto& n0 = submap.node(c.region[0]);
        auto& n1 = submap.node(c.region[1]);
        if (n0.deleted || n1.deleted) continue;

        // Condition 1: at least one endpoint has degree < 3.
        if (n0.degree() >= 3 && n1.degree() >= 3) continue;

        // Condition 2: merged weight ≤ γ.
        // When a chord is removed, arcs that were separated by the chord
        // may merge into a single larger arc.  Per §2.3, "one or both
        // endpoints of the chord might not be vertices of ∂C and might
        // thus disappear," causing adjacent arcs to merge.  The true
        // merged weight must account for these potential merges.
        //
        // We iterate over the *current* arcs of both regions (not the
        // chord's adj_arcs, which may be stale after prior removals in
        // this pass).  For each chord endpoint vertex v, we check if v
        // is still used by another chord in the survivor's incident
        // list.  If not, v would disappear and arcs meeting at v would
        // merge — their edge_counts sum.
        //
        // merged_weight = max over all resulting (merged or individual)
        // arc edge_counts.

        // Determine which chord endpoint edges would disappear.
        // An edge endpoint disappears if (1) it is not at a C-vertex
        // AND (2) it is not referenced by any other chord in either region.
        // §2.2: "those endpoints that are not vertices of C" disappear;
        // vertices of C never disappear.
        auto is_at_c_vertex = [&](std::size_t e, double y) -> bool {
            if (e == NONE) return false;
            if (submap.polygon_ == nullptr) return false;
            if (e >= submap.polygon_->num_edges()) return false;
            const auto& edge = submap.polygon_->edge(e);
            const auto& p_start = submap.polygon_->vertex(edge.start_idx);
            const auto& p_end   = submap.polygon_->vertex(edge.end_idx);
            if (p_start.y == y) return true;
            if (p_end.y   == y) return true;
            return false;
        };

        auto edge_survives = [&](std::size_t e, double y) -> bool {
            if (e == NONE) return false;
            // A C-vertex endpoint always survives.
            if (is_at_c_vertex(e, y)) return true;
            for (std::size_t oci : n0.incident_chords) {
                if (oci == ci) continue;
                const auto& oc = submap.chord(oci);
                if (oc.region[0] == NONE) continue;
                if (oc.left_edge == e || oc.right_edge == e) return true;
            }
            for (std::size_t oci : n1.incident_chords) {
                if (oci == ci) continue;
                const auto& oc = submap.chord(oci);
                if (oc.region[0] == NONE) continue;
                if (oc.left_edge == e || oc.right_edge == e) return true;
            }
            return false;
        };

        bool le_disappears = (c.left_edge != NONE && !edge_survives(c.left_edge, c.y));
        bool re_disappears = (c.right_edge != NONE && !edge_survives(c.right_edge, c.y));

        // For each disappearing vertex, find the two arcs (across both
        // regions) that meet at it and sum their edge_counts.  For
        // non-disappearing endpoints, arcs remain separate.
        //
        // We accumulate groups: arcs touching a disappearing vertex get
        // their edge_counts summed; all others contribute individually.
        std::size_t left_group = 0;
        std::size_t right_group = 0;
        std::size_t bridge_group = 0;  // arcs touching BOTH endpoints
        std::size_t merged_weight = 0;

        auto scan_arcs = [&](const SubmapNode& nd) {
            for (std::size_t ai : nd.arcs) {
                const auto& a = submap.arc(ai);
                if (a.first_edge == NONE) continue;
                std::size_t ec = a.edge_count;

                std::size_t alo = std::min(a.first_edge, a.last_edge);
                std::size_t ahi = std::max(a.first_edge, a.last_edge);

                bool touches_left = le_disappears && (
                    c.left_edge >= alo && c.left_edge <= ahi + 1);
                bool touches_right = re_disappears && (
                    c.right_edge >= alo && c.right_edge <= ahi + 1);

                if (touches_left && touches_right) {
                    // §3.3: arc spans both disappearing endpoints.
                    // Track separately to avoid double-counting when
                    // computing merged weight.
                    bridge_group += ec;
                    left_group   += ec;
                    right_group  += ec;
                } else if (touches_left) {
                    left_group  += ec;
                } else if (touches_right) {
                    right_group += ec;
                } else {
                    merged_weight = std::max(merged_weight, ec);
                }
            }
        };
        scan_arcs(n0);
        scan_arcs(n1);

        // §3.3: When both endpoints disappear and an arc bridges them,
        // all arcs touching either endpoint merge into a single arc.
        // Weight = left_group + right_group − bridge_group (de-dup).
        if (bridge_group > 0) {
            std::size_t combined = left_group + right_group - bridge_group;
            merged_weight = std::max(merged_weight, combined);
        } else {
            merged_weight = std::max(merged_weight,
                                     std::max(left_group, right_group));
        }

        if (merged_weight > gamma) continue;

        submap.remove_chord(ci);
    }
}

} // namespace chazelle
