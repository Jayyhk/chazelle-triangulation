#include "visibility/submap.h"

#include <algorithm>
#include <cassert>

namespace chazelle {

// --- Construction ---

std::size_t Submap::add_node() {
    std::size_t idx = nodes_.size();
    nodes_.emplace_back();
    return idx;
}

std::size_t Submap::add_arc(ArcStructure arc) {
    std::size_t idx = arc_sequence_.size();
    arc_sequence_.push_back(std::move(arc));
    return idx;
}

std::size_t Submap::add_chord(Chord chord) {
    std::size_t idx = chords_.size();

    // Register this chord with its two endpoint regions.
    for (int s = 0; s < 2; ++s) {
        std::size_t r = chord.region[s];
        if (r != NONE && r < nodes_.size()) {
            nodes_[r].incident_chords.push_back(idx);
        }
    }

    // Populate adj_arcs: find arcs in the endpoint regions whose edge
    // ranges are adjacent to the chord's vertex endpoints.  An arc is
    // adjacent to a chord endpoint vertex v if v is at the boundary of
    // the arc's edge range: v == lo (first vertex) or v == hi + 1 (last
    // vertex), where [lo, hi] is the arc's edge range.
    chord.num_adj_arcs = 0;
    for (int s = 0; s < 2; ++s) {
        std::size_t r = chord.region[s];
        if (r == NONE || r >= nodes_.size()) continue;
        for (std::size_t ai : nodes_[r].arcs) {
            if (chord.num_adj_arcs >= 4) break;
            const auto& a = arc_sequence_[ai];
            if (a.first_edge == NONE) continue;
            std::size_t alo = std::min(a.first_edge, a.last_edge);
            std::size_t ahi = std::max(a.first_edge, a.last_edge);
            // Check if either chord endpoint is at a boundary of this arc.
            bool adj = false;
            if (chord.left_vertex != NONE) {
                if (chord.left_vertex == alo || chord.left_vertex == ahi + 1)
                    adj = true;
            }
            if (chord.right_vertex != NONE) {
                if (chord.right_vertex == alo || chord.right_vertex == ahi + 1)
                    adj = true;
            }
            if (adj) {
                chord.adj_arcs[chord.num_adj_arcs++] = ai;
            }
        }
    }

    chords_.push_back(std::move(chord));
    return idx;
}

bool Submap::split_arc_at_vertex(std::size_t arc_idx, std::size_t vertex_idx) {
    assert(arc_idx < arc_sequence_.size());
    auto& arc = arc_sequence_[arc_idx];
    if (arc.first_edge == NONE) return false;

    std::size_t lo = std::min(arc.first_edge, arc.last_edge);
    std::size_t hi = std::max(arc.first_edge, arc.last_edge);

    // vertex_idx corresponds to the start of edge vertex_idx (if
    // vertex_idx > 0) or the boundary vertex.  An arc covering edges
    // [lo, hi] spans vertices [lo, hi+1].  Splitting at vertex_idx
    // means the left arc covers edges [lo, vertex_idx-1] (vertices
    // [lo, vertex_idx]) and the right arc covers edges [vertex_idx, hi]
    // (vertices [vertex_idx, hi+1]).
    //
    // We only split if vertex_idx is strictly interior: lo < vertex_idx <= hi.
    if (vertex_idx <= lo || vertex_idx > hi) return false;

    // The split edge boundary: left arc gets [lo, vertex_idx - 1],
    // right arc gets [vertex_idx, hi].
    std::size_t region = arc.region_node;

    // Create the right half arc.
    // The right half is near the end of the original arc, so its
    // first_side should reflect the side at the split point (which is
    // on the same side as the original arc's last_side for the portion
    // closest to last_edge), and last_side stays as the original's.
    ArcStructure right_arc;
    right_arc.first_edge = vertex_idx;
    right_arc.last_edge = hi;
    right_arc.first_side = arc.last_side;
    right_arc.last_side = arc.last_side;
    right_arc.region_node = region;
    right_arc.edge_count = hi - vertex_idx + 1;

    // Truncate the original arc to be the left half.
    // The left half is near the start, so first_side stays and
    // last_side should match first_side (interior to the arc's
    // starting portion).
    arc.first_edge = lo;
    arc.last_edge = vertex_idx - 1;
    arc.last_side = arc.first_side;
    arc.edge_count = vertex_idx - lo;

    // Add the right arc and register it with the region.
    std::size_t new_ai = arc_sequence_.size();
    arc_sequence_.push_back(right_arc);
    if (region != NONE && region < nodes_.size()) {
        nodes_[region].arcs.push_back(new_ai);
    }

    return true;
}

std::size_t Submap::remove_chord(std::size_t chord_idx) {
    assert(chord_idx < chords_.size());
    Chord& c = chords_[chord_idx];

    std::size_t r0 = c.region[0];
    std::size_t r1 = c.region[1];
    assert(r0 != NONE && r1 != NONE);

    // Determine survivor (lower index).
    std::size_t survivor = std::min(r0, r1);
    std::size_t absorbed = std::max(r0, r1);

    // Guard: if both regions are the same, just remove the chord
    // from the incident list and return.
    if (survivor == absorbed) {
        auto& sinc = nodes_[survivor].incident_chords;
        for (std::size_t i = 0; i < sinc.size(); ++i) {
            if (sinc[i] == chord_idx) {
                sinc[i] = sinc.back();
                sinc.pop_back();
                break;
            }
        }
        c.region[0] = NONE;
        c.region[1] = NONE;
        c.left_vertex = NONE;
        c.right_vertex = NONE;
        return survivor;
    }

    // Transfer all chords from absorbed to survivor.
    for (std::size_t ci : nodes_[absorbed].incident_chords) {
        if (ci == chord_idx) continue; // skip the chord being removed
        // Update the chord's region pointer.
        for (int s = 0; s < 2; ++s) {
            if (chords_[ci].region[s] == absorbed) {
                chords_[ci].region[s] = survivor;
            }
        }
        nodes_[survivor].incident_chords.push_back(ci);
    }

    // Transfer arcs from absorbed to survivor.
    for (std::size_t ai : nodes_[absorbed].arcs) {
        arc_sequence_[ai].region_node = survivor;
        nodes_[survivor].arcs.push_back(ai);
    }

    // Remove the chord from survivor's incident list.
    // Use swap-and-pop for O(1) removal instead of O(k) erase.
    auto& sinc = nodes_[survivor].incident_chords;
    for (std::size_t i = 0; i < sinc.size(); ++i) {
        if (sinc[i] == chord_idx) {
            sinc[i] = sinc.back();
            sinc.pop_back();
            break;
        }
    }

    // Mark absorbed node as deleted.
    nodes_[absorbed].deleted = true;
    nodes_[absorbed].incident_chords.clear();
    nodes_[absorbed].arcs.clear();

    // Fully invalidate the removed chord so that it is not counted
    // by any subsequent deduplication check (e.g. parent_chords in
    // the down-phase).
    std::size_t removed_lv = c.left_vertex;
    std::size_t removed_rv = c.right_vertex;
    c.region[0] = NONE;
    c.region[1] = NONE;
    c.left_vertex = NONE;
    c.right_vertex = NONE;

    // §2.3: "once removed, a chord of zero length ceases to separate
    // any arcs."  More generally, if a chord endpoint is not a vertex
    // of ∂C shared with another chord, it disappears and the arcs that
    // met at it merge into a single arc.
    //
    // For each chord endpoint vertex v, check if v is still referenced
    // by any remaining chord.  If not, find the two arcs in the
    // survivor's arc list that meet at v and merge them.
    auto merge_arcs_at_vertex = [&](std::size_t v) {
        if (v == NONE) return;
        // Check if v is still used by another chord.
        for (std::size_t ci : nodes_[survivor].incident_chords) {
            const auto& ch = chords_[ci];
            if (ch.region[0] == NONE) continue; // already removed
            if (ch.left_vertex == v || ch.right_vertex == v) return;
        }
        // v has disappeared — find the two arcs meeting at v and merge.
        std::size_t a_idx = NONE, b_idx = NONE;
        auto& arcs = nodes_[survivor].arcs;
        for (std::size_t i = 0; i < arcs.size(); ++i) {
            std::size_t ai = arcs[i];
            auto& a = arc_sequence_[ai];
            if (a.first_edge == NONE) continue;
            std::size_t alo = std::min(a.first_edge, a.last_edge);
            std::size_t ahi = std::max(a.first_edge, a.last_edge);
            if (v == alo || v == ahi + 1) {
                if (a_idx == NONE) a_idx = i;
                else { b_idx = i; break; }
            }
        }
        if (a_idx == NONE || b_idx == NONE) return;

        std::size_t ai_a = arcs[a_idx];
        std::size_t ai_b = arcs[b_idx];
        auto& arcA = arc_sequence_[ai_a];
        auto& arcB = arc_sequence_[ai_b];

        // Determine which end of each arc touches v to establish the
        // correct first/last edge/side of the merged arc.
        std::size_t ahi_a = std::max(arcA.first_edge, arcA.last_edge);
        std::size_t ahi_b = std::max(arcB.first_edge, arcB.last_edge);

        // A touches v at its "high" end if v == ahi+1, at its "low" end
        // if v == alo.  B likewise.  The merged arc's first edge/side
        // comes from A's non-touching end and last edge/side from B's
        // non-touching end (or vice versa).
        bool a_high = (v == ahi_a + 1);
        bool b_high = (v == ahi_b + 1);

        std::size_t new_first, new_last;
        Side new_first_side, new_last_side;
        if (a_high && !b_high) {
            // A's low end → B's high end
            new_first = arcA.first_edge;
            new_first_side = arcA.first_side;
            new_last = arcB.last_edge;
            new_last_side = arcB.last_side;
        } else if (!a_high && b_high) {
            // B's low end → A's high end
            new_first = arcB.first_edge;
            new_first_side = arcB.first_side;
            new_last = arcA.last_edge;
            new_last_side = arcA.last_side;
        } else if (a_high && b_high) {
            // Both touch at high end — take A's low end as first, B's
            // low end as last (order doesn't matter for weight).
            new_first = arcA.first_edge;
            new_first_side = arcA.first_side;
            new_last = arcB.first_edge;
            new_last_side = arcB.first_side;
        } else {
            // Both touch at low end.
            new_first = arcA.last_edge;
            new_first_side = arcA.last_side;
            new_last = arcB.last_edge;
            new_last_side = arcB.last_side;
        }

        // Merge into arcA, invalidate arcB.
        arcA.first_edge = new_first;
        arcA.first_side = new_first_side;
        arcA.last_edge = new_last;
        arcA.last_side = new_last_side;
        arcA.edge_count = arcA.edge_count + arcB.edge_count;

        // Invalidate arcB.
        arcB.first_edge = NONE;
        arcB.last_edge = NONE;
        arcB.edge_count = 0;
        arcB.region_node = NONE;

        // Remove b_idx from the survivor's arcs list.
        arcs[b_idx] = arcs.back();
        arcs.pop_back();
    };

    merge_arcs_at_vertex(removed_lv);
    merge_arcs_at_vertex(removed_rv);

    // Recompute survivor's weight.
    recompute_weight(survivor);

    return survivor;
}

// --- Weight ---

void Submap::recompute_weight(std::size_t node_idx) {
    assert(node_idx < nodes_.size());
    SubmapNode& nd = nodes_[node_idx];
    nd.weight = 0;
    for (std::size_t ai : nd.arcs) {
        nd.weight = std::max(nd.weight, arc_sequence_[ai].edge_count);
    }
}

void Submap::recompute_all_weights() {
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (!nodes_[i].deleted) {
            recompute_weight(i);
        }
    }
}

// --- Queries ---

bool Submap::is_conformal() const {
    for (auto& nd : nodes_) {
        if (!nd.deleted && nd.degree() > 4) return false;
    }
    return true;
}

bool Submap::is_semigranular(std::size_t gamma) const {
    for (auto& nd : nodes_) {
        if (!nd.deleted && nd.weight > gamma) return false;
    }
    return true;
}

std::size_t Submap::max_degree() const {
    std::size_t md = 0;
    for (auto& nd : nodes_) {
        if (!nd.deleted) md = std::max(md, nd.degree());
    }
    return md;
}

// --- Double Identification (§2.4) ---

std::vector<std::size_t> Submap::double_identify(std::size_t edge_idx,
                                                  double /*y*/) const {
    // Per §2.4: since we know the location of the two endpoints of C in
    // the arc-sequence table, we break the circular arc sequence into two
    // linear subsequences and binary search each using the edge name.
    // If the search lands on a contiguous interval (multiple arcs sharing
    // the same edge), the paper says to disambiguate by y-coordinate via
    // a second binary search.  Since the Submap doesn't store the polygon
    // geometry, full y-disambiguation is performed by the caller (which
    // has access to the polygon and can compare y against arc edge
    // y-ranges).  Here we do structural disambiguation: prefer arcs
    // where edge_idx is interior, and use the Side flag to distinguish
    // copies of the double boundary ∂C.

    std::vector<std::size_t> result;
    if (arc_sequence_.empty()) return result;

    // Helper: check whether arc at index i contains edge_idx.
    auto arc_contains = [&](std::size_t i) -> bool {
        const auto& a = arc_sequence_[i];
        if (a.first_edge == NONE) return false;
        std::size_t lo = std::min(a.first_edge, a.last_edge);
        std::size_t hi = std::max(a.first_edge, a.last_edge);
        return edge_idx >= lo && edge_idx <= hi;
    };

    // Helper: representative edge of arc i for ordering along ∂C.
    auto arc_rep = [&](std::size_t i) -> std::size_t {
        const auto& a = arc_sequence_[i];
        return (a.first_edge != NONE) ? a.first_edge : 0;
    };

    // Binary search in range [lo_idx, hi_idx) for an arc containing
    // edge_idx.  Arcs are ordered by their edge range along ∂C, so
    // we search by comparing edge_idx against each arc's representative.
    auto binary_search_arcs = [&](std::size_t lo_idx, std::size_t hi_idx) {
        if (lo_idx >= hi_idx) return;

        // Standard binary search for the first arc whose representative
        // edge is >= edge_idx.
        std::size_t left = lo_idx, right = hi_idx;
        while (left < right) {
            std::size_t mid = left + (right - left) / 2;
            if (arc_rep(mid) < edge_idx) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        // The matching arcs form a contiguous interval around `left`.
        // Check left and its neighbors for containment.
        // Expand leftward.
        std::size_t scan_lo = left;
        while (scan_lo > lo_idx && arc_contains(scan_lo - 1)) {
            --scan_lo;
        }
        // Expand rightward.
        std::size_t scan_hi = left;
        while (scan_hi < hi_idx && arc_contains(scan_hi)) {
            ++scan_hi;
        }
        // Also check one position before scan_lo (binary search may
        // land just past a multi-edge arc that contains edge_idx).
        if (scan_lo > lo_idx && scan_lo == left && arc_contains(scan_lo - 1)) {
            --scan_lo;
        }

        for (std::size_t i = scan_lo; i < scan_hi; ++i) {
            if (arc_contains(i)) {
                result.push_back(i);
            }
        }
    };

    // Split the arc-sequence at the C-endpoint arcs (start_arc, end_arc)
    // into two monotone subsequences along ∂C.
    if (start_arc != NONE && end_arc != NONE &&
        start_arc < arc_sequence_.size() && end_arc < arc_sequence_.size()) {
        std::size_t split1 = std::min(start_arc, end_arc);
        std::size_t split2 = std::max(start_arc, end_arc);
        // Search each half.
        binary_search_arcs(0, split1 + 1);
        binary_search_arcs(split1 + 1, split2 + 1);
        if (split2 + 1 < arc_sequence_.size()) {
            binary_search_arcs(split2 + 1, arc_sequence_.size());
        }
    } else {
        // No endpoint info — search the full table.
        binary_search_arcs(0, arc_sequence_.size());
    }

    // Structural disambiguation: when multiple arcs contain edge_idx,
    // prefer those where edge_idx is strictly interior (not at a
    // boundary), since boundary arcs only partially overlap the edge.
    // Then deduplicate by region_node.
    if (result.size() > 1) {
        std::vector<std::size_t> interior;
        for (std::size_t i : result) {
            const auto& a = arc_sequence_[i];
            if (a.first_edge == NONE) continue;
            std::size_t lo = std::min(a.first_edge, a.last_edge);
            std::size_t hi = std::max(a.first_edge, a.last_edge);
            if (edge_idx > lo && edge_idx < hi) {
                interior.push_back(i);
            }
        }
        if (!interior.empty()) result = std::move(interior);
    }
    if (result.size() > 1) {
        // Remove exact duplicates (same region_node) to avoid redundancy.
        std::vector<std::size_t> deduped;
        for (std::size_t i : result) {
            bool dup = false;
            for (std::size_t j : deduped) {
                if (arc_sequence_[i].region_node == arc_sequence_[j].region_node) {
                    dup = true;
                    break;
                }
            }
            if (!dup) deduped.push_back(i);
        }
        if (!deduped.empty()) result = std::move(deduped);
    }

    return result;
}

} // namespace chazelle
