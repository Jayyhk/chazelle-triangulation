#include "visibility/submap.h"
#include "geometry/polygon.h"

#include <algorithm>
#include <cassert>
#include <limits>

namespace chazelle {

void Submap::set_chain_info(std::size_t c_start, std::size_t c_end,
                           const Polygon* poly) {
    c_start_vertex = c_start;
    c_end_vertex   = c_end;
    polygon_       = poly;
}

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
    // ranges contain (or border) the chord's endpoint edges.  Per §3.1
    // Remark 1, chord endpoints are (edge, y) pairs.  An arc covering
    // edges [alo, ahi] is adjacent to a chord endpoint on edge e if e
    // is in [alo, ahi] (the endpoint lies on one of the arc's edges)
    // or e == ahi + 1 (the endpoint is at the arc's upper boundary).
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
            // Check if either chord endpoint edge is at/within the arc's range.
            bool adj = false;
            if (chord.left_edge != NONE) {
                if (chord.left_edge >= alo && chord.left_edge <= ahi + 1)
                    adj = true;
            }
            if (chord.right_edge != NONE) {
                if (chord.right_edge >= alo && chord.right_edge <= ahi + 1)
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

    // Determine the correct side flags for each half.
    // For arcs that do NOT double-back (first_side == last_side), both
    // halves inherit the same side.  For double-backing arcs
    // (first_side != last_side), the arc wraps around an endpoint of C.
    // The transition point is at c_start_vertex or c_end_vertex,
    // whichever falls within [lo, hi].  Edges before the transition
    // have first_side; edges at/after the transition have last_side.
    Side left_last_side  = arc.first_side;
    Side right_first_side = arc.last_side;
    if (arc.first_side != arc.last_side) {
        // Find the transition edge: the C-endpoint that falls within
        // the arc's edge range [lo, hi].  The transition vertex is the
        // C-endpoint; edges before it are on first_side, at/after on
        // last_side.
        std::size_t transition = NONE;
        if (c_start_vertex != NONE && c_start_vertex > lo && c_start_vertex <= hi)
            transition = c_start_vertex;
        if (c_end_vertex != NONE && c_end_vertex > lo && c_end_vertex <= hi)
            transition = c_end_vertex;

        if (transition != NONE) {
            // vertex_idx relative to the transition determines which
            // side the split point is on.
            if (vertex_idx < transition) {
                // Split point is before the transition: left half
                // is entirely on first_side.  Right half still crosses.
                left_last_side   = arc.first_side;
                right_first_side = arc.first_side;
            } else if (vertex_idx > transition) {
                // Split point is after the transition: left half
                // crosses.  Right half is entirely on last_side.
                left_last_side   = arc.last_side;
                right_first_side = arc.last_side;
            } else {
                // Split is exactly at the transition.
                left_last_side   = arc.first_side;
                right_first_side = arc.last_side;
            }
        }
    }

    // Create the right half arc.
    ArcStructure right_arc;
    right_arc.first_edge = vertex_idx;
    right_arc.last_edge = hi;
    right_arc.first_side = right_first_side;
    right_arc.last_side = arc.last_side;
    right_arc.region_node = region;
    right_arc.edge_count = hi - vertex_idx + 1;

    // Truncate the original arc to be the left half.
    arc.first_edge = lo;
    arc.last_edge = vertex_idx - 1;
    arc.last_side = left_last_side;
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
        c.left_edge = NONE;
        c.right_edge = NONE;
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
    std::size_t removed_le = c.left_edge;
    std::size_t removed_re = c.right_edge;
    c.region[0] = NONE;
    c.region[1] = NONE;
    c.left_edge = NONE;
    c.right_edge = NONE;

    // §2.2: "the removal of a chord entails removing not only the
    // chord itself but also those endpoints that are not vertices
    // of C, and glueing back ∂C at those points."  Vertices of C
    // never disappear; only non-C chord endpoints (introduced by
    // augmented visibility maps) can disappear.
    //
    // For each chord endpoint vertex v:
    //   - If v is a vertex of C (v ∈ [c_start_vertex, c_end_vertex]),
    //     it stays; arcs remain separate.
    //   - Else, if v is no longer referenced by any remaining chord,
    //     find the two arcs meeting at v and merge them.
    // §2.2: "the removal of a chord entails removing not only the
    // chord itself but also those endpoints that are not vertices of
    // C, and glueing back ∂C at those points."
    // Per §3.1 Remark 1, chord endpoints are (edge, y) pairs.
    // When an endpoint disappears, arcs sharing that edge merge.
    auto merge_arcs_at_edge = [&](std::size_t edge_e) {
        if (edge_e == NONE) return;
        // §2.2: Vertices of C are never removed.  For edge-based
        // endpoints, an endpoint at a C-vertex boundary (edge_e ==
        // c_start_vertex or c_end_vertex) is a C-vertex and stays.
        if (c_start_vertex != NONE && c_end_vertex != NONE &&
            edge_e >= c_start_vertex && edge_e <= c_end_vertex) {
            return;
        }
        // Check if this edge is still used by another chord endpoint.
        for (std::size_t ci : nodes_[survivor].incident_chords) {
            const auto& ch = chords_[ci];
            if (ch.region[0] == NONE) continue; // already removed
            if (ch.left_edge == edge_e || ch.right_edge == edge_e) return;
        }
        // edge_e endpoint has disappeared — find two arcs meeting
        // at this edge and merge.  With edge-based endpoints, both
        // arcs include edge_e at their boundary.
        std::size_t a_idx = NONE, b_idx = NONE;
        auto& arcs = nodes_[survivor].arcs;
        for (std::size_t i = 0; i < arcs.size(); ++i) {
            std::size_t ai = arcs[i];
            auto& a = arc_sequence_[ai];
            if (a.first_edge == NONE) continue;
            std::size_t alo = std::min(a.first_edge, a.last_edge);
            std::size_t ahi = std::max(a.first_edge, a.last_edge);
            if (edge_e == alo || edge_e == ahi || edge_e == ahi + 1) {
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

        // A touches edge_e at its "high" end if edge_e == ahi+1 or
        // edge_e == ahi, at its "low" end if edge_e == alo.  B likewise.
        // The merged arc's first edge/side comes from A's non-touching
        // end and last edge/side from B's non-touching end.
        bool a_high = (edge_e == ahi_a + 1 || edge_e == ahi_a);
        bool b_high = (edge_e == ahi_b + 1 || edge_e == ahi_b);

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

    merge_arcs_at_edge(removed_le);
    merge_arcs_at_edge(removed_re);

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
        if (!nd.deleted && nd.arcs.size() > 4) return false;
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

// --- Normal Form (§2.3) ---

void Submap::normalize() {
    if (arc_sequence_.empty()) return;

    const std::size_t n = arc_sequence_.size();

    // Build sort permutation: sorted indices by position along ∂C.
    // Sort key = min(first_edge, last_edge).  Invalidated arcs
    // (first_edge == NONE) are pushed to the end.
    std::vector<std::size_t> perm(n);
    for (std::size_t i = 0; i < n; ++i) perm[i] = i;
    std::sort(perm.begin(), perm.end(), [&](std::size_t a, std::size_t b) {
        const auto& aa = arc_sequence_[a];
        const auto& bb = arc_sequence_[b];
        std::size_t ka = (aa.first_edge != NONE)
            ? std::min(aa.first_edge, aa.last_edge)
            : std::numeric_limits<std::size_t>::max();
        std::size_t kb = (bb.first_edge != NONE)
            ? std::min(bb.first_edge, bb.last_edge)
            : std::numeric_limits<std::size_t>::max();
        return ka < kb;
    });

    // Build the inverse permutation: inv[old_index] = new_index.
    std::vector<std::size_t> inv(n);
    for (std::size_t i = 0; i < n; ++i) {
        inv[perm[i]] = i;
    }

    // Permute the arc_sequence_ table.
    std::vector<ArcStructure> sorted(n);
    for (std::size_t i = 0; i < n; ++i) {
        sorted[i] = arc_sequence_[perm[i]];
    }
    arc_sequence_ = std::move(sorted);

    // Update arc indices in all nodes.
    for (auto& nd : nodes_) {
        if (nd.deleted) continue;
        for (auto& ai : nd.arcs) {
            ai = inv[ai];
        }
    }

    // Update adj_arcs in all chords.
    for (auto& c : chords_) {
        for (std::size_t k = 0; k < c.num_adj_arcs; ++k) {
            if (c.adj_arcs[k] < n) {
                c.adj_arcs[k] = inv[c.adj_arcs[k]];
            }
        }
    }

    // Update start_arc / end_arc.
    if (start_arc != NONE && start_arc < n) {
        start_arc = inv[start_arc];
    }
    if (end_arc != NONE && end_arc < n) {
        end_arc = inv[end_arc];
    }

    // Strip invalidated arcs (those with first_edge == NONE) from the
    // end of the now-sorted table.
    while (!arc_sequence_.empty() &&
           arc_sequence_.back().first_edge == NONE) {
        arc_sequence_.pop_back();
    }
}

// --- Double Identification (§2.4) ---

std::vector<std::size_t> Submap::double_identify(std::size_t edge_idx,
                                                  double y) const {
    // Per §2.4: since we know the location of the two endpoints of C in
    // the arc-sequence table, we break the circular arc sequence into two
    // linear subsequences and binary search each using the edge name.
    // If the search lands on a contiguous interval (multiple arcs sharing
    // the same edge), disambiguate by y-coordinate via a second binary
    // search.

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

    // §2.4: "we can conceptually break up the circular arc sequence into
    // two linear sequences and perform in each of them a binary search,
    // using the name of the containing edge as a query."
    //
    // Each of the two C-sides is monotone in edge index — one ascending,
    // the other descending.  We binary-search a virtual sequence of
    // `count` elements starting at real index `offset`, wrapping around
    // the table when necessary.  `ascending` gives the direction.
    auto binary_search_arcs_seq = [&](std::size_t offset, std::size_t count,
                                      bool ascending) {
        if (count == 0) return;
        std::size_t n_arcs = arc_sequence_.size();
        auto real = [&](std::size_t vi) -> std::size_t {
            return (offset + vi) % n_arcs;
        };

        std::size_t left = 0, right = count;
        while (left < right) {
            std::size_t mid = left + (right - left) / 2;
            bool go_right = ascending
                ? (arc_rep(real(mid)) < edge_idx)
                : (arc_rep(real(mid)) > edge_idx);
            if (go_right) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        // Scan around the landing point for containment.
        std::size_t scan_lo = left;
        while (scan_lo > 0 && arc_contains(real(scan_lo - 1))) {
            --scan_lo;
        }
        std::size_t scan_hi = left;
        while (scan_hi < count && arc_contains(real(scan_hi))) {
            ++scan_hi;
        }

        for (std::size_t i = scan_lo; i < scan_hi; ++i) {
            if (arc_contains(real(i))) {
                result.push_back(real(i));
            }
        }
    };

    // Break the circular arc-sequence at the two C-endpoint arcs into
    // exactly two linear sequences (one per side of C) and binary-
    // search each.
    if (start_arc != NONE && end_arc != NONE &&
        start_arc < arc_sequence_.size() && end_arc < arc_sequence_.size()) {
        std::size_t n_arcs = arc_sequence_.size();

        // Sequence 1: from start_arc to end_arc going forward in the
        // table (wrapping if start_arc > end_arc).
        std::size_t count1 = (end_arc >= start_arc)
            ? (end_arc - start_arc + 1)
            : (n_arcs - start_arc + end_arc + 1);
        bool seq1_asc = true;
        if (count1 >= 2) {
            std::size_t first_r = start_arc;
            std::size_t last_r  = (start_arc + count1 - 1) % n_arcs;
            seq1_asc = (arc_rep(first_r) <= arc_rep(last_r));
        }
        binary_search_arcs_seq(start_arc, count1, seq1_asc);

        // Sequence 2: the complement — from end_arc+1 back to
        // start_arc−1 (the other side of C, also wrapping).
        std::size_t seq2_start = (end_arc + 1) % n_arcs;
        std::size_t count2 = n_arcs - count1;
        if (count2 > 0) {
            bool seq2_asc = true;
            if (count2 >= 2) {
                std::size_t first_r = seq2_start;
                std::size_t last_r  = (seq2_start + count2 - 1) % n_arcs;
                seq2_asc = (arc_rep(first_r) <= arc_rep(last_r));
            }
            binary_search_arcs_seq(seq2_start, count2, seq2_asc);
        }
    } else {
        // No endpoint info — linear scan (degenerate case).
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            if (arc_contains(i)) result.push_back(i);
        }
    }

    // §2.4 y-coordinate disambiguation: "We can disambiguate by
    // pursuing the binary search, now using, say, the y-coordinate
    // of q as a query."
    //
    // When multiple arcs contain edge_idx, use the polygon's vertex
    // y-coordinates to determine which arc(s) the point (edge_idx, y)
    // actually lies on.  An arc covering edges [lo, hi] spans the
    // y-range of the polygon vertices from lo to hi+1.  A point at
    // y is on this arc if y falls within that range on edge_idx.
    if (result.size() > 1 && polygon_ != nullptr) {
        std::vector<std::size_t> y_filtered;
        for (std::size_t i : result) {
            const auto& a = arc_sequence_[i];
            if (a.first_edge == NONE) continue;
            // The arc covers edge_idx.  Check if y falls within the
            // y-range of edge_idx's endpoints.
            if (edge_idx < polygon_->num_edges()) {
                const auto& e = polygon_->edge(edge_idx);
                double y0 = polygon_->vertex(e.start_idx).y;
                double y1 = polygon_->vertex(e.end_idx).y;
                double ylo = std::min(y0, y1);
                double yhi = std::max(y0, y1);
                // The point y is on this arc's edge if it's within
                // the y-range (with tolerance for endpoints).
                if (y >= ylo - 1e-12 && y <= yhi + 1e-12) {
                    y_filtered.push_back(i);
                }
            } else {
                y_filtered.push_back(i);
            }
        }
        if (!y_filtered.empty()) result = std::move(y_filtered);
    }

    // Structural disambiguation: when multiple arcs still remain,
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
