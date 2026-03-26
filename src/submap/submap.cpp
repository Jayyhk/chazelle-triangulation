/// src/submap/submap.cpp

#include "submap.h"
#include "../polygon/polygon.h"

#include <algorithm>

namespace chazelle {

std::size_t Submap::add_node() {
    std::size_t idx = nodes_.size();
    nodes_.push_back(SubmapNode{});
    return idx;
}

std::size_t Submap::add_arc(Arc arc) {
    std::size_t idx = arc_sequence_.size();
    // [C91 §2.4]: Maintain cached LEFT/RIGHT boundary.
    if (arc.first_side == LEFT) {
        left_right_boundary_ = idx + 1;
    }
    arc_sequence_.push_back(arc);
    return idx;
}

std::size_t Submap::add_chord(Chord chord) {
    std::size_t idx = chords_.size();
    chords_.push_back(chord);

    // [C91 §2.4(i)]: Update adjacency — each region knows its
    // incident chords.
    for (std::size_t r : chord.region) {
        if (r != NONE && r < nodes_.size()) {
            nodes_[r].incident_chords.push_back(idx);
        }
    }

    return idx;
}

std::size_t Submap::remove_chord(std::size_t chord_idx,
                                  const Polygon& polygon) {
    assert(chord_idx < chords_.size());
    auto& c = chords_[chord_idx];
    assert(c.region[0] != NONE && c.region[1] != NONE);

    std::size_t r0 = c.region[0];
    std::size_t r1 = c.region[1];

    // [C91 §2.2]: "the removal of a chord entails removing not only
    // the chord itself but also those endpoints that are not vertices
    // of C, and glueing back ∂C at those points."
    //
    // "A chord has two endpoints; none, one, or two of which are
    // vertices of C."  Only NON-vertex endpoints get glued back
    // (arcs merged).  Vertex endpoints remain as legitimate ∂C
    // boundaries.
    //
    // adj_arcs are in pairs: [0]+[1] meet at one chord endpoint,
    // [2]+[3] meet at the other.  We determine which endpoint each
    // pair is at by checking the chord's left_edge/right_edge
    // against polygon vertices.

    // Check if each chord endpoint is a polygon vertex.
    // A chord endpoint at edge e with y matching vertex v's y
    // is at vertex v if the symbolic y equals the vertex's symbolic y.
    auto endpoint_is_polygon_vertex = [&](std::size_t edge,
                                          double ey,
                                          std::size_t ey_tag) -> bool {
        // Check if the chord's y matches either endpoint of the edge.
        if (edge >= polygon.num_edges()) return false;
        const auto& e = polygon.edge(edge);
        const auto& v_start = polygon.vertex(e.start_idx);
        const auto& v_end   = polygon.vertex(e.end_idx);
        SymbolicY chord_y{ey, ey_tag};
        return symbolic_y_equal(chord_y, symbolic_y_of(v_start)) ||
               symbolic_y_equal(chord_y, symbolic_y_of(v_end));
    };

    bool left_is_vertex  = endpoint_is_polygon_vertex(c.left_edge, c.y, c.y_tag);
    bool right_is_vertex = endpoint_is_polygon_vertex(c.right_edge, c.y, c.y_tag);

    // [C91 §2.2]: "glueing back ∂C at those points."
    // [C91 §2.4(ii)]: "two, three, or four arcs adjacent to it."
    //
    // adj_arcs are stored in pairs per chord endpoint:
    //   4 adj_arcs (2+2): [L0, L1, R0, R1] — merge [L0+L1] and [R0+R1]
    //   3 adj_arcs (2+1): [P0, P1, S0, -]  — merge [P0+P1] only
    //   2 adj_arcs (1+1): [S0, S1, -, -]   — no merging
    //
    // With 2: each endpoint has 1 arc → nothing to merge.
    // With 4: pairs [0]+[1] at left endpoint, [2]+[3] at right.
    // With 3: [0]+[1] at the 2-arc endpoint, [2] at the 1-arc endpoint.
    //
    // Only merge at non-vertex endpoints.
    // Collect arc indices to erase (from merging).
    std::vector<std::size_t> arcs_to_erase;

    {
        auto do_merge = [&](std::size_t ai, std::size_t aj) {
            if (ai == NONE || aj == NONE) return;
            if (ai >= arc_sequence_.size() || aj >= arc_sequence_.size()) return;

            auto& a_keep = arc_sequence_[ai];
            auto& a_dead = arc_sequence_[aj];

            a_keep.last_edge = a_dead.last_edge;
            a_keep.last_side = a_dead.last_side;
            a_keep.edge_count += a_dead.edge_count;

            arcs_to_erase.push_back(aj);

            // Update adj_arcs of neighboring chords.  O(1): degree ≤ 4.
            for (std::size_t ri : {r0, r1}) {
                if (ri >= nodes_.size()) continue;
                for (std::size_t ci : nodes_[ri].incident_chords) {
                    if (ci == chord_idx) continue;
                    auto& other = chords_[ci];
                    if (other.region[0] == NONE) continue;
                    for (std::size_t k = 0; k < other.num_adj_arcs; ++k) {
                        if (other.adj_arcs[k] == aj)
                            other.adj_arcs[k] = ai;
                    }
                }
            }
        };

        if (c.num_adj_arcs == 4) {
            if (!left_is_vertex)
                do_merge(c.adj_arcs[0], c.adj_arcs[1]);
            if (!right_is_vertex)
                do_merge(c.adj_arcs[2], c.adj_arcs[3]);
        } else if (c.num_adj_arcs == 3) {
            if (!left_is_vertex)
                do_merge(c.adj_arcs[0], c.adj_arcs[1]);
        }
    }

    // Erase dead arcs from arc_sequence_ and update all indices.
    // Sort erase list descending so erasing doesn't shift earlier indices.
    std::sort(arcs_to_erase.begin(), arcs_to_erase.end(),
              std::greater<std::size_t>());
    for (std::size_t dead_idx : arcs_to_erase) {
        arc_sequence_.erase(arc_sequence_.begin()
                            + static_cast<long>(dead_idx));

        // Adjust all arc indices > dead_idx by -1.
        for (std::size_t ci = 0; ci < chords_.size(); ++ci) {
            auto& ch = chords_[ci];
            for (std::size_t k = 0; k < ch.num_adj_arcs; ++k) {
                if (ch.adj_arcs[k] != NONE && ch.adj_arcs[k] > dead_idx)
                    --ch.adj_arcs[k];
            }
        }
        if (start_arc != NONE && start_arc > dead_idx) --start_arc;
        if (end_arc != NONE && end_arc > dead_idx) --end_arc;
        if (left_right_boundary_ > dead_idx) --left_right_boundary_;
    }

    // Reassign arcs from r1 to r0.
    for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
        if (arc_sequence_[i].region_node == r1)
            arc_sequence_[i].region_node = r0;
    }

    // Move r1's incident chords (other than this one) to r0.
    for (std::size_t ci : nodes_[r1].incident_chords) {
        if (ci == chord_idx) continue;
        nodes_[r0].incident_chords.push_back(ci);
    }

    // Remove this chord from r0's incident list.
    {
        auto& ic = nodes_[r0].incident_chords;
        ic.erase(std::remove(ic.begin(), ic.end(), chord_idx), ic.end());
    }

    // Erase the dead node (r1) from nodes_.
    // Adjust all node indices > r1 by -1.
    nodes_.erase(nodes_.begin() + static_cast<long>(r1));

    auto adjust_node = [&](std::size_t& idx) {
        if (idx == NONE) return;
        if (idx == r1) idx = r0;       // shouldn't happen, but safe
        else if (idx > r1) --idx;
    };

    // Fix r0 if it was after r1.
    if (r0 > r1) --r0;

    for (auto& arc : arc_sequence_)
        adjust_node(arc.region_node);
    for (auto& ch : chords_) {
        adjust_node(ch.region[0]);
        adjust_node(ch.region[1]);
    }
    for (auto& nd : nodes_) {
        for (auto& ci : nd.incident_chords) {
            // chord indices don't change (we haven't erased the chord yet)
        }
    }

    // Erase the dead chord from chords_.
    // Adjust all chord indices > chord_idx by -1.
    chords_.erase(chords_.begin() + static_cast<long>(chord_idx));

    for (auto& nd : nodes_) {
        for (auto& ci : nd.incident_chords) {
            if (ci > chord_idx) --ci;
        }
    }
    for (auto& ch : chords_) {
        for (std::size_t k = 0; k < ch.num_adj_arcs; ++k) {
            // adj_arcs reference arc indices, not chord indices — no change
        }
    }

    return r0;
}

void Submap::assert_tree_property() const {
    // [C91 §2.2]: "the dual graph of a submap is itself a tree."
    // For a tree: num_nodes == num_edges + 1.
    // No dead entries — nodes and chords are actually erased.
    assert(nodes_.size() == chords_.size() + 1 &&
           "§2.2: submap tree property violated "
           "(num_regions ≠ num_chords + 1)");
}

void Submap::check_invariants() const {
    assert_tree_property();

    // Every chord's regions must be valid nodes.
    for (std::size_t i = 0; i < chords_.size(); ++i) {
        const auto& c = chords_[i];
        assert(c.region[0] < nodes_.size() &&
               "§2.2: chord region[0] out of range");
        assert(c.region[1] < nodes_.size() &&
               "§2.2: chord region[1] out of range");
    }

    // Every arc's region_node must be a valid node.
    for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
        assert(arc_sequence_[i].region_node < nodes_.size() &&
               "§2.2: arc region_node out of range");
    }

    // [C91 §2.4(ii)]: Each chord's adj_arcs must point to valid
    // arcs belonging to one of the chord's two regions.
    for (std::size_t i = 0; i < chords_.size(); ++i) {
        const auto& c = chords_[i];
        for (std::size_t j = 0; j < c.num_adj_arcs; ++j) {
            std::size_t ai = c.adj_arcs[j];
            assert(ai < arc_sequence_.size() &&
                   "§2.4(ii): adj_arc index out of range");
            const auto& a = arc_sequence_[ai];
            assert((a.region_node == c.region[0] ||
                    a.region_node == c.region[1]) &&
                   "§2.4(ii): adj_arc must belong to one of the "
                   "chord's endpoint regions");
        }
    }

    // [C91 §2.4(iii)]: Arc-sequence table must be in ∂C order:
    // LEFT arcs first, then RIGHT arcs.
    {
        bool seen_right = false;
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            if (arc_sequence_[i].first_side == RIGHT) seen_right = true;
            if (seen_right) {
                assert(arc_sequence_[i].first_side == RIGHT &&
                       "§2.4(iii): LEFT arc after RIGHT arc in "
                       "arc-sequence table violates ∂C order");
            }
        }
    }

    // [C91 §2.4(iii)]: Endpoint pointers, if set, must be valid.
    if (start_arc != NONE) {
        assert(start_arc < arc_sequence_.size() &&
               "§2.4(iii): start_arc out of range");
    }
    if (end_arc != NONE) {
        assert(end_arc < arc_sequence_.size() &&
               "§2.4(iii): end_arc out of range");
    }
}

Submap::DoubleIdentifyResult
Submap::double_identify(std::size_t edge_idx, SymbolicY y) const {
    DoubleIdentifyResult result;
    if (arc_sequence_.empty()) return result;

    // [C91 §2.4]: "we can conceptually break up the circular arc
    // sequence into two linear sequences and perform in each of
    // them a binary search, using the name of the containing edge
    // as a query."
    //
    // The arc-sequence table is in ∂C order: LEFT arcs first
    // (ascending first_edge), then RIGHT arcs (descending first_edge).

    // [C91 §2.4] (tex 144): Use cached LEFT/RIGHT boundary — O(1).
    std::size_t left_begin = 0;
    std::size_t left_end   = left_right_boundary_;
    std::size_t right_begin = left_end;
    std::size_t right_end   = arc_sequence_.size();

    // Binary search helper: find an arc whose edge range contains
    // edge_idx, then expand to collect all such arcs.
    auto search_half = [&](std::size_t lo, std::size_t hi, bool ascending) {
        if (lo >= hi) return;

        // Binary search for an arc whose first_edge is close to edge_idx.
        std::size_t blo = lo, bhi = hi;
        while (blo < bhi) {
            std::size_t mid = blo + (bhi - blo) / 2;
            if (ascending) {
                if (arc_sequence_[mid].first_edge < edge_idx)
                    blo = mid + 1;
                else
                    bhi = mid;
            } else {
                if (arc_sequence_[mid].first_edge > edge_idx)
                    blo = mid + 1;
                else
                    bhi = mid;
            }
        }

        // Expand around blo to collect all arcs containing edge_idx.
        auto arc_contains_edge = [&](std::size_t ai) -> bool {
            const auto& a = arc_sequence_[ai];
            std::size_t elo = std::min(a.first_edge, a.last_edge);
            std::size_t ehi = std::max(a.first_edge, a.last_edge);
            return edge_idx >= elo && edge_idx <= ehi;
        };

        // Collect candidate arcs (same edge) into a local buffer,
        // then disambiguate by y if needed.
        std::size_t candidates[DoubleIdentifyResult::MAX];
        std::size_t num_cand = 0;

        if (blo > lo) {
            for (std::size_t i = blo; i-- > lo; ) {
                if (arc_contains_edge(i)) {
                    if (num_cand < DoubleIdentifyResult::MAX)
                        candidates[num_cand++] = i;
                } else break;
            }
        }
        for (std::size_t i = blo; i < hi; ++i) {
            if (arc_contains_edge(i)) {
                if (num_cand < DoubleIdentifyResult::MAX)
                    candidates[num_cand++] = i;
            } else break;
        }

        // [C91 §2.4]: "We can disambiguate by pursuing the binary
        // search, now using, say, the y-coordinate of q as a query."
        // If only 1-2 candidates, no disambiguation needed (at most
        // two arcs per side at a non-chord point).  Otherwise,
        // filter by y: keep arcs whose key_y bracket contains y.
        if (num_cand <= 2 || y.tag == SOS_NONE) {
            for (std::size_t ci = 0; ci < num_cand; ++ci) {
                if (result.count < DoubleIdentifyResult::MAX)
                    result.push(candidates[ci]);
            }
        } else {
            // Multiple arcs at the same edge: find arcs whose
            // y-range covers the query y.  Adjacent arcs in the
            // arc-sequence are ordered by key_y, so the query y
            // falls between two consecutive key_y values.  Keep
            // the arc just before and just after the query y.
            for (std::size_t ci = 0; ci < num_cand; ++ci) {
                std::size_t ai = candidates[ci];
                const auto& a = arc_sequence_[ai];
                // An arc spans from its key_y to the next arc's
                // key_y.  Accept if the query y is close to this
                // arc's key_y, or if this is the last arc on the
                // edge (it extends to the edge boundary).
                // Conservative: accept all candidates within ±1
                // position of the y-match.
                if (symbolic_y_equal(a.key_symbolic_y(), y)) {
                    if (result.count < DoubleIdentifyResult::MAX)
                        result.push(ai);
                } else if (ci + 1 < num_cand) {
                    // Check if y falls between this arc and the next.
                    const auto& next_a = arc_sequence_[candidates[ci + 1]];
                    SymbolicY lo_y = a.key_symbolic_y();
                    SymbolicY hi_y = next_a.key_symbolic_y();
                    if (symbolic_y_less(lo_y, hi_y)) {
                        if (symbolic_y_leq(lo_y, y) && symbolic_y_leq(y, hi_y)) {
                            if (result.count < DoubleIdentifyResult::MAX)
                                result.push(ai);
                        }
                    } else {
                        if (symbolic_y_leq(hi_y, y) && symbolic_y_leq(y, lo_y)) {
                            if (result.count < DoubleIdentifyResult::MAX)
                                result.push(ai);
                        }
                    }
                } else {
                    // Last candidate — accept (boundary arc).
                    if (result.count < DoubleIdentifyResult::MAX)
                        result.push(ai);
                }
            }
        }
    };

    // [C91 §2.4]: LEFT half — ascending first_edge.
    search_half(left_begin, left_end, /*ascending=*/true);

    // [C91 §2.4]: RIGHT half — descending first_edge.
    search_half(right_begin, right_end, /*ascending=*/false);

    return result;
}

void Submap::build_tree_decomposition() {
    tree_decomp_.build(*this);
}

std::size_t Submap::region_weight(std::size_t node_idx) const noexcept {
    assert(node_idx < nodes_.size());
    // [C91 §2.2]: "the maximum number of nonnull length edges in
    // any of its arcs."
    std::size_t max_count = 0;
    for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
        if (arc_sequence_[i].region_node == node_idx) {
            if (arc_sequence_[i].edge_count > max_count)
                max_count = arc_sequence_[i].edge_count;
        }
    }
    return max_count;
}

bool Submap::is_conformal() const noexcept {
    // [C91 §2.3]: "conformal submaps [are] those with node-degree
    // at most 4."
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i].degree() > 4)
            return false;
    }
    return true;
}

bool Submap::is_semigranular(std::size_t gamma) const noexcept {
    // [C91 §2.3]: "every node of its tree has weight at most γ."
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (region_weight(i) > gamma)
            return false;
    }
    return true;
}

std::size_t Submap::simulated_contraction_weight(
        std::size_t chord_idx) const noexcept {
    assert(chord_idx < chords_.size());
    const auto& c = chords_[chord_idx];
    if (c.region[0] == NONE || c.region[1] == NONE) return 0;

    // [C91 §2.3]: "contracting any edge... produces a new node whose
    // weight."  When a chord is removed, adjacent arc pairs merge
    // ("glueing back ∂C").  The merged arc's edge_count can exceed
    // either original's.
    //
    // Start with the max of all non-adjacent arcs.
    std::size_t max_count = 0;
    for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
        const auto& a = arc_sequence_[i];
        if (a.region_node == c.region[0] ||
            a.region_node == c.region[1]) {
            if (a.edge_count > max_count)
                max_count = a.edge_count;
        }
    }

    // Now account for adj_arc merging: pairs of adjacent arcs that
    // would be glued together when the chord is removed.
    // adj_arcs are stored in pairs: [0]+[1] merge, [2]+[3] merge.
    for (std::size_t i = 0; i + 1 < c.num_adj_arcs; i += 2) {
        std::size_t ai = c.adj_arcs[i];
        std::size_t aj = c.adj_arcs[i + 1];
        if (ai == NONE || aj == NONE) continue;
        if (ai >= arc_sequence_.size() || aj >= arc_sequence_.size()) continue;
        std::size_t merged = arc_sequence_[ai].edge_count
                           + arc_sequence_[aj].edge_count;
        if (merged > max_count)
            max_count = merged;
    }

    return max_count;
}

bool Submap::is_granular(std::size_t gamma) const noexcept {
    // [C91 §2.3] condition (i): all weights ≤ γ.
    if (!is_semigranular(gamma)) return false;

    // [C91 §2.3]: "by default, if (i) holds but the submap has no
    // exit chord, it is still said to be γ-granular."
    if (chords_.empty()) return true;

    // [C91 §2.3] condition (ii): "contracting any edge incident upon
    // at least one node of degree less than 3 produces a new node
    // whose weight exceeds γ."
    for (std::size_t ci = 0; ci < chords_.size(); ++ci) {
        const auto& c = chords_[ci];

        std::size_t d0 = nodes_[c.region[0]].degree();
        std::size_t d1 = nodes_[c.region[1]].degree();

        // Only check edges incident upon at least one degree < 3 node.
        if (d0 >= 3 && d1 >= 3) continue;

        // Contraction must produce weight > γ.
        if (simulated_contraction_weight(ci) <= gamma)
            return false;
    }
    return true;
}

} // namespace chazelle
