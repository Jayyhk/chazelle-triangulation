/// src/submap/submap.cpp

#include "submap.h"
#include "../polygon/polygon.h"

#include <algorithm>

namespace chazelle {

// ── SubmapNode ──────────────────────────────────────────────────

std::size_t SubmapNode::degree() const noexcept {
    // Invariant: remove_chord actively maintains incident_chords —
    // dead chord indices are erased at removal time and only live
    // chords from merged regions are copied.  So incident_chords
    // contains only live entries for live nodes.
    return incident_chords.size();
}

// ── Construction ────────────────────────────────────────────────

std::size_t Submap::add_node() {
    std::size_t idx = nodes_.size();
    nodes_.push_back(SubmapNode{});
    return idx;
}

std::size_t Submap::add_arc(Arc arc) {
    std::size_t idx = arc_sequence_.size();

    // [C91 §2.4(iii)] (tex 138): "The arc-structures are stored in a
    // table in the order corresponding to a canonical traversal of ∂C."
    // LEFT arcs (ascending first_edge) must precede all RIGHT arcs
    // (descending first_edge).  left_right_boundary_ is computed as a
    // side effect of insertion order, so out-of-order insertion silently
    // corrupts double_identify's binary search ranges.
    if (!arc_sequence_.empty()) {
        Side prev_side = arc_sequence_.back().first_side;
        if (arc.first_side == LEFT) {
            assert(prev_side == LEFT &&
                   "§2.4(iii): LEFT arc added after RIGHT arc "
                   "violates ∂C order");
            assert(arc.first_edge >= arc_sequence_.back().first_edge &&
                   "§2.4(iii): LEFT arcs must have ascending "
                   "first_edge");
        } else {
            if (prev_side == RIGHT) {
                assert(arc.first_edge <= arc_sequence_.back().first_edge &&
                       "§2.4(iii): RIGHT arcs must have descending "
                       "first_edge");
            }
        }
    }

    // [C91 §2.4]: Maintain cached LEFT/RIGHT boundary.
    if (arc.first_side == LEFT) {
        left_right_boundary_ = idx + 1;
    }
    arc_sequence_.push_back(arc);
    return idx;
}

std::size_t Submap::add_chord(Chord chord) {
    std::size_t idx = chords_.size();

    // [C91 §2.4(ii)]: Validate per-endpoint adj arcs.
    assert(chord.left_adj.count >= 1 && chord.left_adj.count <= 2 &&
           "§2.4(ii): LEFT endpoint must have 1 or 2 adjacent arcs");
    assert(chord.right_adj.count >= 1 && chord.right_adj.count <= 2 &&
           "§2.4(ii): RIGHT endpoint must have 1 or 2 adjacent arcs");
    auto check_adj = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            assert(adj.arcs[k] != NONE &&
                   adj.arcs[k] < arc_sequence_.size() &&
                   "§2.4(ii): adj_arc must be a valid arc index");
            assert(!arc_sequence_[adj.arcs[k]].dead &&
                   "§2.4(ii): adj_arc must be live");
            assert((arc_sequence_[adj.arcs[k]].region_node ==
                        chord.region[0] ||
                    arc_sequence_[adj.arcs[k]].region_node ==
                        chord.region[1]) &&
                   "§2.4(ii): adj_arc must belong to one of the "
                   "chord's two regions");
        }
    };
    check_adj(chord.left_adj);
    check_adj(chord.right_adj);

    // [C91 §2.2] (tex 102): dual graph is a tree — no self-loops.
    assert(chord.region[0] != chord.region[1] &&
           "§2.2: chord must connect two distinct regions "
           "(tree, no self-loops)");

    chords_.push_back(chord);

    // [C91 §2.4(i)]: Update adjacency — each region knows its
    // incident chords.  Both regions must be valid.
    for (std::size_t r : chord.region) {
        assert(r != NONE && r < nodes_.size() &&
               "§2.4(i): chord must connect two valid regions");
        nodes_[r].incident_chords.push_back(idx);
    }

    return idx;
}

// ── Chord removal (O(1) via tombstones) ─────────────────────────

// [C91 §2.2] (tex 94): "the removal of a chord entails removing not
// only the chord itself but also those endpoints that are not vertices
// of C, and glueing back ∂C at those points."
//
// [C91 §3.3]: "can be enforced in time linear in the size of the
// submap tree, that is, O(n₁/γ₁ + n₂/γ₂ + 1)."  This requires O(1)
// per removal.  We achieve this by tombstoning dead entries instead
// of erasing from vectors.  Indices remain stable.  Dead entries are
// stripped by compact() before putting S in normal form.

std::size_t Submap::remove_chord(std::size_t chord_idx,
                                  const Polygon& polygon) {
    assert(chord_idx < chords_.size());
    auto& c = chords_[chord_idx];
    assert(!c.dead && "removing an already-dead chord");
    assert(c.region[0] != NONE && c.region[1] != NONE);

    std::size_t r0 = c.region[0];
    std::size_t r1 = c.region[1];
    assert(!nodes_[r0].dead && !nodes_[r1].dead);

    // Check if each chord endpoint is a polygon vertex.
    auto endpoint_is_polygon_vertex = [&](std::size_t edge,
                                          double ey,
                                          std::size_t ey_tag) -> bool {
        assert(edge < polygon.num_edges() &&
               "§2.2: chord endpoint edge must be a valid edge index");
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
    // Merge arc pairs at endpoints.  O(1): at most 2 pairs.
    auto do_merge = [&](std::size_t ai, std::size_t aj) {
        assert(ai != NONE && ai < arc_sequence_.size() &&
               "§2.4(ii): adj_arc must be valid");
        assert(aj != NONE && aj < arc_sequence_.size() &&
               "§2.4(ii): adj_arc must be valid");
        assert(!arc_sequence_[ai].dead && !arc_sequence_[aj].dead);

        auto& a_keep = arc_sequence_[ai];
        auto& a_dead = arc_sequence_[aj];

        assert(a_keep.last_edge == a_dead.first_edge &&
               "§2.4(iii): adjacent arcs at chord endpoint "
               "must share the junction edge");

        // [C91 §2.2] (tex 106): merged edge count.
        std::size_t shared_edge = a_dead.first_edge;
        std::size_t shared_nonnull =
            polygon.count_nonnull_edges(shared_edge, shared_edge);

        a_keep.last_edge = a_dead.last_edge;
        a_keep.last_side = a_dead.last_side;
        a_keep.edge_count = a_keep.edge_count + a_dead.edge_count
                          - shared_nonnull;

        // Tombstone the dead arc.
        a_dead.dead = true;
        compacted_ = false;  // arc-sequence has dead entries now

        // [C91 §2.4(iii)]: If the dead arc was start_arc or end_arc,
        // redirect to the surviving arc (which now covers the merged
        // range including the endpoint).
        if (start_arc == aj) start_arc = ai;
        if (end_arc == aj)   end_arc = ai;

        // Update adj arcs of neighboring chords to point to the
        // surviving arc.  O(degree) = O(1) for conformal submaps.
        auto replace_arc = [&](Chord::AdjArcs& adj) {
            for (std::size_t k = 0; k < adj.count; ++k)
                if (adj.arcs[k] == aj) adj.arcs[k] = ai;
        };
        for (std::size_t ri : {r0, r1}) {
            for (std::size_t ci : nodes_[ri].incident_chords) {
                if (ci == chord_idx) continue;
                auto& other = chords_[ci];
                if (other.dead) continue;
                replace_arc(other.left_adj);
                replace_arc(other.right_adj);
            }
        }
    };

    // [C91 §2.2] (tex 108): "once removed, a chord of zero length
    // ceases to separate any arcs."
    //
    // Use the stored flag — do NOT recompute from edge/side fields.
    // The default Chord has left_side=LEFT, right_side=RIGHT, so
    // any exit chord with left_edge==right_edge would be misclassified
    // as an NLC if we recomputed here.
    bool is_null_length = c.is_null_length;

    if (is_null_length) {
        assert(c.left_adj.count == 1 && c.right_adj.count == 1 &&
               "§2.2: null-length chord at vertex must have exactly "
               "1 adjacent arc per side");
        // [C91 §2.2] (tex 108): NLC endpoints are polygon vertices.
        // Vertex endpoints do NOT trigger arc merging (§2.2 tex 94).
        // The null-length arc (right_adj) is kept alive and absorbed
        // into the main region by the reassignment step below —
        // it correctly represents the y-extremum vertex on r0's boundary.
        // do_merge is NOT called here.
    } else {
        // Standard logic for chords with nonzero length: merge if NOT a vertex.
        if (!left_is_vertex) {
            // [C91 §2.2] (tex 94): a non-vertex chord endpoint always has
            // exactly 2 adjacent arcs (one ending, one starting) to be glued.
            assert(c.left_adj.count == 2 &&
                   "§2.2 (tex 94): non-vertex chord endpoint must have "
                   "exactly 2 adjacent arcs for glueing");
            do_merge(c.left_adj.arcs[0], c.left_adj.arcs[1]);
        }
        if (!right_is_vertex) {
            assert(c.right_adj.count == 2 &&
                   "§2.2 (tex 94): non-vertex chord endpoint must have "
                   "exactly 2 adjacent arcs for glueing");
            do_merge(c.right_adj.arcs[0], c.right_adj.arcs[1]);
        }
    }

    // Reassign r1's arcs to r0.  Use chord→arc adjacency for O(1).
    // [C91 §2.4(ii)] (tex 137): traverse r1's incident chords' adj arcs.
    //
    // Two call sites with different invariants:
    //   (a) Other chords' adj arcs: replace_arc already updated them
    //       to point to surviving arcs, so all must be live.
    //   (b) The removed chord's own adj arcs: do_merge may have killed
    //       one arc at each non-vertex endpoint.
    auto reassign_live = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            std::size_t ai = adj.arcs[k];
            assert(ai != NONE && ai < arc_sequence_.size() &&
                   "§2.4(ii): adj_arc index must be valid");
            assert(!arc_sequence_[ai].dead &&
                   "§2.4(ii): other chords' adj arcs must be live "
                   "(replace_arc maintains this)");
            if (arc_sequence_[ai].region_node == r1) {
                arc_sequence_[ai].region_node = r0;
            }
        }
    };
    auto reassign_maybe_dead = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            std::size_t ai = adj.arcs[k];
            assert(ai != NONE && ai < arc_sequence_.size() &&
                   "§2.4(ii): adj_arc index must be valid");
            if (arc_sequence_[ai].dead) continue;
            if (arc_sequence_[ai].region_node == r1) {
                arc_sequence_[ai].region_node = r0;
            }
        }
    };
    for (std::size_t ci : nodes_[r1].incident_chords) {
        if (ci == chord_idx) continue;
        const auto& ch = chords_[ci];
        if (ch.dead) continue;
        reassign_live(ch.left_adj);
        reassign_live(ch.right_adj);
    }
    // The removed chord's own adj arcs: do_merge may have killed one
    // arc at each non-vertex endpoint, so skip dead arcs here.
    reassign_maybe_dead(c.left_adj);
    reassign_maybe_dead(c.right_adj);

    // Move r1's incident chords (other than this one) to r0.
    // Also update the chord's region[] to point to r0.
    for (std::size_t ci : nodes_[r1].incident_chords) {
        if (ci == chord_idx) continue;
        auto& ch = chords_[ci];
        if (ch.dead) continue;
        nodes_[r0].incident_chords.push_back(ci);
        // Fix the chord's region pointer from r1 → r0.
        if (ch.region[0] == r1) ch.region[0] = r0;
        if (ch.region[1] == r1) ch.region[1] = r0;
    }

    // Remove this chord from r0's incident list.  O(degree) = O(1).
    {
        auto& ic = nodes_[r0].incident_chords;
        ic.erase(std::remove(ic.begin(), ic.end(), chord_idx), ic.end());
    }

    // Tombstone the chord and the dead region.
    c.dead = true;
    nodes_[r1].dead = true;

    return r0;
}

// ── Tombstone compaction ────────────────────────────────────────

// [C91 §3.3]: "We can now put S in normal form."  compact() strips
// dead entries and rebuilds index mappings in O(m).

void Submap::compact() {
    // [C91 §2.2] (tex 94): After cascaded chord removals, some live
    // arcs may still point to dead regions.  This happens when region
    // R2 is absorbed into R1 (chord A removed), then R1 is later
    // absorbed into R0 (chord B removed) — arcs that were adjacent to
    // chord A and belonged to R2 were reassigned to R1 by A's removal,
    // but B's removal only reassigns arcs reachable via chord adjacency,
    // missing the orphaned arcs.
    //
    // Fix: build a forwarding table from dead chords (each dead chord
    // records region[1] → region[0]), resolve chains, then fixup arcs.
    // O(m) total — does not change the paper's complexity.
    {
        std::vector<std::size_t> forward(nodes_.size(), NONE);
        for (const auto& ch : chords_) {
            if (!ch.dead) continue;
            // remove_chord always kills region[1] and keeps region[0].
            std::size_t dead_r = ch.region[1];
            std::size_t live_r = ch.region[0];
            if (dead_r < nodes_.size() && nodes_[dead_r].dead &&
                forward[dead_r] == NONE) {
                forward[dead_r] = live_r;
            }
        }
        // Resolve chains with path compression.
        auto resolve = [&](std::size_t r) -> std::size_t {
            // Phase 1: find root
            std::size_t root = r;
            while (root < nodes_.size() && nodes_[root].dead &&
                   forward[root] != NONE) {
                root = forward[root];
            }
            // Phase 2: compress path
            while (r != root && r < nodes_.size() && nodes_[r].dead) {
                std::size_t next = forward[r];
                forward[r] = root;
                r = next;
            }
            return root;
        };
        for (auto& a : arc_sequence_) {
            if (a.dead) continue;
            if (a.region_node < nodes_.size() &&
                nodes_[a.region_node].dead) {
                a.region_node = resolve(a.region_node);
            }
        }
    }

    // Build old→new index maps for each table.

    // Nodes.
    std::vector<std::size_t> node_map(nodes_.size(), NONE);
    {
        std::size_t j = 0;
        for (std::size_t i = 0; i < nodes_.size(); ++i) {
            if (!nodes_[i].dead) {
                node_map[i] = j;
                if (j != i) nodes_[j] = std::move(nodes_[i]);
                ++j;
            }
        }
        nodes_.resize(j);
    }

    // Chords.
    std::vector<std::size_t> chord_map(chords_.size(), NONE);
    {
        std::size_t j = 0;
        for (std::size_t i = 0; i < chords_.size(); ++i) {
            if (!chords_[i].dead) {
                chord_map[i] = j;
                if (j != i) chords_[j] = std::move(chords_[i]);
                ++j;
            }
        }
        chords_.resize(j);
    }

    // Arcs.
    std::vector<std::size_t> arc_map(arc_sequence_.size(), NONE);
    {
        std::size_t j = 0;
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            if (!arc_sequence_[i].dead) {
                arc_map[i] = j;
                if (j != i) arc_sequence_[j] = std::move(arc_sequence_[i]);
                ++j;
            }
        }
        arc_sequence_.resize(j);
    }

    // Rewrite all cross-references using the maps.

    // Arcs: region_node.
    for (auto& a : arc_sequence_) {
        assert(a.region_node != NONE && node_map[a.region_node] != NONE);
        a.region_node = node_map[a.region_node];
    }

    // Chords: region[], adj arcs.
    auto remap_adj = [&](Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            if (adj.arcs[k] != NONE) {
                assert(arc_map[adj.arcs[k]] != NONE);
                adj.arcs[k] = arc_map[adj.arcs[k]];
            }
        }
    };
    for (auto& ch : chords_) {
        for (auto& r : ch.region) {
            assert(r != NONE && node_map[r] != NONE);
            r = node_map[r];
        }
        remap_adj(ch.left_adj);
        remap_adj(ch.right_adj);
    }

    // Nodes: incident_chords.
    for (auto& nd : nodes_) {
        std::size_t w = 0;
        for (std::size_t ci : nd.incident_chords) {
            if (chord_map[ci] != NONE)
                nd.incident_chords[w++] = chord_map[ci];
        }
        nd.incident_chords.resize(w);
    }

    // Endpoint arc pointers.
    if (start_arc != NONE) {
        assert(arc_map[start_arc] != NONE);
        start_arc = arc_map[start_arc];
    }
    if (end_arc != NONE) {
        assert(arc_map[end_arc] != NONE);
        end_arc = arc_map[end_arc];
    }

    // Recompute left_right_boundary_.
    left_right_boundary_ = 0;
    for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
        if (arc_sequence_[i].first_side == LEFT)
            left_right_boundary_ = i + 1;
    }

    // [C91 §2.4] (tex 144): arc-sequence is now free of dead entries.
    compacted_ = true;

    // [C91 §2.4(iv)]: tree decomposition indices reference the
    // pre-compaction chord/node tables and are now stale.  Invalidate
    // to prevent silent use of wrong indices before rebuild.
    tree_decomp_ = TreeDecomposition{};
}

// ── Live counts ─────────────────────────────────────────────────

std::size_t Submap::num_live_nodes() const noexcept {
    std::size_t n = 0;
    for (const auto& nd : nodes_) if (!nd.dead) ++n;
    return n;
}

std::size_t Submap::num_live_chords() const noexcept {
    std::size_t n = 0;
    for (const auto& ch : chords_) if (!ch.dead) ++n;
    return n;
}

std::size_t Submap::num_live_arcs() const noexcept {
    std::size_t n = 0;
    for (const auto& a : arc_sequence_) if (!a.dead) ++n;
    return n;
}

// ── Invariant checks ────────────────────────────────────────────

void Submap::assert_tree_property() const {
    // [C91 §2.2]: "the dual graph of a submap is itself a tree."
    // For a tree: num_live_regions == num_live_chords + 1.
    assert(num_live_nodes() == num_live_chords() + 1 &&
           "§2.2: submap tree property violated "
           "(num_regions ≠ num_chords + 1)");
}

void Submap::check_invariants() const {
    assert_tree_property();

    // Every live chord's regions must be live nodes.
    for (std::size_t i = 0; i < chords_.size(); ++i) {
        const auto& c = chords_[i];
        if (c.dead) continue;
        assert(c.region[0] < nodes_.size() && !nodes_[c.region[0]].dead &&
               "§2.2: chord region[0] invalid or dead");
        assert(c.region[1] < nodes_.size() && !nodes_[c.region[1]].dead &&
               "§2.2: chord region[1] invalid or dead");
    }

    // Every live arc's region_node must be a live node.
    for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
        if (arc_sequence_[i].dead) continue;
        assert(arc_sequence_[i].region_node < nodes_.size() &&
               !nodes_[arc_sequence_[i].region_node].dead &&
               "§2.2: arc region_node invalid or dead");
    }

    // [C91 §2.4(ii)]: Each live chord's adj arcs must point to live
    // arcs belonging to one of the chord's two regions.
    for (std::size_t i = 0; i < chords_.size(); ++i) {
        const auto& c = chords_[i];
        if (c.dead) continue;
        auto check_adj = [&](const Chord::AdjArcs& adj) {
            for (std::size_t j = 0; j < adj.count; ++j) {
                std::size_t ai = adj.arcs[j];
                assert(ai < arc_sequence_.size() &&
                       !arc_sequence_[ai].dead &&
                       "§2.4(ii): adj_arc index invalid or dead");
                const auto& a = arc_sequence_[ai];
                assert((a.region_node == c.region[0] ||
                        a.region_node == c.region[1]) &&
                       "§2.4(ii): adj_arc must belong to one of "
                       "the chord's endpoint regions");
            }
        };
        check_adj(c.left_adj);
        check_adj(c.right_adj);
    }

    // [C91 §2.4(iii)] (tex 138): ∂C ordering of live arcs.
    // Dead arcs are gaps in the table; skip them for ordering checks.
    {
        bool seen_right = false;
        std::size_t prev_live = NONE;
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            if (arc_sequence_[i].dead) continue;
            if (arc_sequence_[i].first_side == RIGHT) seen_right = true;
            if (seen_right) {
                assert(arc_sequence_[i].first_side == RIGHT &&
                       "§2.4(iii): LEFT arc after RIGHT arc in "
                       "arc-sequence table violates ∂C order");
            }
            if (prev_live != NONE) {
                if (arc_sequence_[prev_live].first_side == LEFT &&
                    arc_sequence_[i].first_side == LEFT) {
                    assert(arc_sequence_[i].first_edge >=
                           arc_sequence_[prev_live].first_edge &&
                           "§2.4(iii): LEFT arcs must have ascending "
                           "first_edge for double_identify binary search");
                }
                if (arc_sequence_[prev_live].first_side == RIGHT &&
                    arc_sequence_[i].first_side == RIGHT) {
                    assert(arc_sequence_[i].first_edge <=
                           arc_sequence_[prev_live].first_edge &&
                           "§2.4(iii): RIGHT arcs must have descending "
                           "first_edge for double_identify binary search");
                }
            }
            prev_live = i;
        }
    }

    // [C91 §2.4(iii)] (tex 138): endpoint arc pointers.
    if (num_live_arcs() > 0) {
        assert(start_arc != NONE &&
               "§2.4(iii): start_arc must be set when arcs exist");
        assert(end_arc != NONE &&
               "§2.4(iii): end_arc must be set when arcs exist");
        assert(start_arc < arc_sequence_.size() &&
               !arc_sequence_[start_arc].dead &&
               "§2.4(iii): start_arc out of range or dead");
        assert(end_arc < arc_sequence_.size() &&
               !arc_sequence_[end_arc].dead &&
               "§2.4(iii): end_arc out of range or dead");

        assert(start_vertex != NONE &&
               "§2.4(iii): start_vertex must be set when arcs exist");
        {
            // [C91 §2.4(iii)]: start_arc must pass through C's
            // start vertex.  Vertex v is on arc [first_edge, last_edge]
            // iff edge v-1 or edge v is in the arc's range (v is the
            // start of edge v and the end of edge v-1).
            // For edge-index convention: vertex v ∈ arc iff
            //   first_edge ≤ v ≤ last_edge  (v starts edge v in range)
            //   OR v == first_edge  (v is the start of the first edge)
            // Simplified: first_edge ≤ v ≤ last_edge + 1 for non-wrapped.
            // For wrapped arcs, extend to C endpoint.
            const auto& sa = arc_sequence_[start_arc];
            std::size_t elo = std::min(sa.first_edge, sa.last_edge);
            std::size_t ehi = std::max(sa.first_edge, sa.last_edge) + 1;
            if (sa.first_side != sa.last_side)
                elo = std::min(elo, start_vertex);
            assert(start_vertex >= elo && start_vertex <= ehi &&
                   "§2.4(iii): start_arc must pass through "
                   "start_vertex");
        }

        assert(end_vertex != NONE && end_vertex > 0 &&
               "§2.4(iii): end_vertex must be set when arcs exist");
        {
            // [C91 §2.4(iii)]: end_arc must pass through C's end
            // vertex.  end_vertex is the last vertex of C; the last
            // edge of C is end_vertex - 1.
            const auto& ea = arc_sequence_[end_arc];
            std::size_t c_end_edge = end_vertex - 1;
            std::size_t elo = std::min(ea.first_edge, ea.last_edge);
            std::size_t ehi = std::max(ea.first_edge, ea.last_edge);
            if (ea.first_side != ea.last_side)
                ehi = std::max(ehi, c_end_edge);
            assert(c_end_edge >= elo && c_end_edge <= ehi &&
                   "§2.4(iii): end_arc must pass through "
                   "end_vertex");
        }
    }
}

// ── Double identification ───────────────────────────────────────

Submap::DoubleIdentifyResult
Submap::double_identify(std::size_t edge_idx, SymbolicY y,
                         const Polygon& polygon) const {
    DoubleIdentifyResult result;
    if (arc_sequence_.empty()) return result;

    // Precondition: submap must be compacted (no dead arcs).
    // The binary search assumes a contiguous sorted table; dead
    // gaps would make the expansion O(k) instead of O(1), violating
    // the paper's O(log m) bound (§2.4 tex 144).
    // O(1) flag check — maintained by remove_chord (sets false) and
    // compact() (sets true).
    assert(compacted_ &&
           "§2.4: double_identify requires compacted arc-sequence");

    // [C91 §2] (tex 47): SoS perturbation [10] ensures every point
    // on ∂C has a well-defined symbolic y.  double_identify must
    // disambiguate by y when multiple arcs share an edge (tex 144),
    // which requires a valid tag.
    assert(y.tag != SOS_NONE &&
           "§2.4: double_identify requires a valid SoS y-tag "
           "for y-disambiguation");

    // [C91 §2.4]: "we can conceptually break up the circular arc
    // sequence into two linear sequences and perform in each of
    // them a binary search, using the name of the containing edge
    // as a query."
    //
    // The arc-sequence table is in ∂C order: LEFT arcs first
    // (ascending first_edge), then RIGHT arcs (descending first_edge).

    // [C91 §2.4] (tex 144): "Since we know the location of the two 
    // endpoints of C in the arc-sequence table (i.e., which arcs pass 
    // through them) we can conceptually break up the circular arc 
    // sequence into two linear sequences..."
    //
    // The paper dictates anchoring the split at these two endpoints.
    // In our normal-form storage:
    //   - Sequence 1 (Forward/LEFT): [start_arc, ..., end_arc]
    //   - Sequence 2 (Return/RIGHT): (end_arc, ..., circular wrap to start_arc)
    //
    // For a simple chain C, this corresponds to:
    assert(start_arc == 0 &&
           "§2.4: start_arc MUST be index 0 in a normal-form canonical traversal");

    // [C91 §2.4] (tex 144): "Since we know the location of the two
    // endpoints of C in the arc-sequence table (i.e., which arcs pass
    // through them) we can conceptually break up the circular arc
    // sequence into two linear sequences."
    //
    // The split is at C's end endpoint: end_arc is the LAST LEFT arc
    // (index = left_right_boundary_ - 1).  end_arc + 1 equals
    // left_right_boundary_, giving [LEFT arcs | RIGHT arcs].
    assert(end_arc == left_right_boundary_ - 1 &&
           "§2.4 (tex 144): end_arc must be the last LEFT arc "
           "(left_right_boundary_ - 1) so the arc-sequence split "
           "produces [LEFT arcs | RIGHT arcs]");

    // Frame the two linear search sequences as prescribed by the paper:
    std::size_t left_begin  = start_arc;
    std::size_t left_end    = end_arc + 1;
    std::size_t right_begin = left_end;
    std::size_t right_end   = arc_sequence_.size();

    // [C91 §2.4] (tex 144): Two-phase binary search.
    //
    // Phase 1: binary search by first_edge → find the contiguous
    //   interval of arcs starting at edge_idx.  O(log m).
    // Phase 2: binary search by key_y within that interval → find
    //   the arc(s) containing query y.  O(log k) ≤ O(log m).
    // Plus O(1) check for a boundary arc (first_edge < edge_idx
    //   but last_edge ≥ edge_idx) at the interval's left neighbor.
    //
    // Total: O(log m).
    auto search_half = [&](std::size_t lo, std::size_t hi, bool ascending) {
        if (lo >= hi) return;

        // ── Phase 1: binary search by first_edge ─────────────────
        // Find blo = first arc with first_edge >= edge_idx (ascending)
        //   or first arc with first_edge <= edge_idx (descending).
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

        // Find the end of the same-first_edge interval: bend = first
        // arc with first_edge != edge_idx after blo.  O(log m).
        std::size_t bend = blo;
        if (bend < hi && arc_sequence_[bend].first_edge == edge_idx) {
            std::size_t slo = blo, shi = hi;
            while (slo < shi) {
                std::size_t mid = slo + (shi - slo) / 2;
                if (arc_sequence_[mid].first_edge == edge_idx)
                    slo = mid + 1;
                else
                    shi = mid;
            }
            bend = slo;
        }

        // [C91 §2.4] (tex 142): "caution must be used since an arc
        // might wrap around both sides of C, something we call
        // double-backing."  Check for at most one boundary arc at
        // blo-1 whose edge range contains edge_idx despite having
        // a different first_edge.  O(1).
        auto arc_contains_edge = [&](std::size_t ai) -> bool {
            assert(ai < arc_sequence_.size() &&
                   "§2.4: arc index must be valid");
            const auto& a = arc_sequence_[ai];
            assert(!a.dead &&
                   "§2.4: arc_contains_edge called on dead arc "
                   "(double_identify requires compacted arc-sequence)");
            std::size_t elo = std::min(a.first_edge, a.last_edge);
            std::size_t ehi = std::max(a.first_edge, a.last_edge);
            if (a.first_side != a.last_side) {
                if (ai == start_arc)
                    elo = std::min(elo, start_vertex);
                if (ai == end_arc) {
                    assert(end_vertex > 0 &&
                           "§2.1: end_vertex must be ≥ 1 "
                           "(curve has at least 2 vertices)");
                    std::size_t c_end_edge = end_vertex - 1;
                    ehi = std::max(ehi, c_end_edge);
                }
            }
            return edge_idx >= elo && edge_idx <= ehi;
        };

        std::size_t boundary_arc = NONE;
        if (blo > lo && arc_contains_edge(blo - 1))
            boundary_arc = blo - 1;

        std::size_t interval_len = bend - blo;

        // No arcs on this edge in this half.
        if (interval_len == 0 && boundary_arc == NONE)
            return;

        // Single arc (either the interval or boundary): no y-disambiguation.
        if (interval_len == 0) {
            result.push(boundary_arc);
            return;
        }
        if (interval_len == 1 && boundary_arc == NONE) {
            result.push(blo);
            return;
        }
        // Two arcs: one interval + one boundary.  Can't detect key_y
        // direction from a single interval arc, so use the edge's
        // geometric y-direction from the Polygon.
        // The junction is at the interval arc's key_y.  Together the
        // two arcs cover the entire edge.
        if (interval_len == 1 && boundary_arc != NONE) {
            SymbolicY junction_y = arc_sequence_[blo].key_symbolic_y();
            if (symbolic_y_equal(junction_y, y)) {
                // At the junction: both arcs pass through.
                result.push(blo);
                result.push(boundary_arc);
            } else {
                // y is strictly on one side of the junction.
                // The boundary arc comes BEFORE the interval arc in
                // ∂C traversal (it's at blo-1).  On the LEFT half,
                // ∂C traversal follows the edge from start→end vertex.
                // The boundary arc covers [edge_start_y, junction_y]
                // and the interval arc covers [junction_y, edge_end_y]
                // (in traversal direction).
                //
                // Determine the edge's y-direction from the Polygon.
                assert(edge_idx < polygon.num_edges());
                const auto& e = polygon.edge(edge_idx);
                SymbolicY start_y = symbolic_y_of(polygon.vertex(e.start_idx));
                bool edge_ascending = symbolic_y_less(start_y,
                    symbolic_y_of(polygon.vertex(e.end_idx)));
                // On the LEFT half, ∂C goes start→end (same as edge).
                // On the RIGHT half, ∂C goes end→start (reversed).
                bool traversal_ascending = ascending ? edge_ascending
                                                     : !edge_ascending;
                // Boundary arc is on the "start" side of junction_y.
                bool y_in_boundary = traversal_ascending
                    ? symbolic_y_less(y, junction_y)
                    : symbolic_y_greater(y, junction_y);
                if (y_in_boundary)
                    result.push(boundary_arc);
                else
                    result.push(blo);
            }
            return;
        }

        // ── Phase 2: binary search by key_y ──────────────────────
        // [C91 §2.4] (tex 144): "We can disambiguate by pursuing the
        // binary search, now using, say, the y-coordinate of q as a
        // query."
        //
        // The key_y direction within a single edge follows the ∂C
        // traversal, which depends on the edge's geometric y-direction
        // (tex 138), NOT on the LEFT/RIGHT half.  Determine it from
        // the interval endpoints.
        bool keys_ascending = (interval_len >= 2) &&
            symbolic_y_leq(arc_sequence_[blo].key_symbolic_y(),
                           arc_sequence_[bend - 1].key_symbolic_y());

        // Binary search within [blo, bend) for the arc whose key_y
        // range contains y.  O(log k).
        //
        // For ascending key_y: find the last arc with key_y <= y.
        // For descending key_y: find the last arc with key_y >= y.
        // The arc at that position (and possibly its neighbor at a
        // key_y boundary) contains the query point.
        std::size_t ylo = blo, yhi = bend;
        while (ylo < yhi) {
            std::size_t mid = ylo + (yhi - ylo) / 2;
            SymbolicY mid_y = arc_sequence_[mid].key_symbolic_y();
            if (keys_ascending) {
                // Find last arc with key_y <= y: go right if mid <= y.
                if (symbolic_y_leq(mid_y, y))
                    ylo = mid + 1;
                else
                    yhi = mid;
            } else {
                // Find last arc with key_y >= y: go right if mid >= y.
                if (symbolic_y_geq(mid_y, y))
                    ylo = mid + 1;
                else
                    yhi = mid;
            }
        }
        // ylo = first arc PAST the query position.
        // The arc containing y is at p = ylo-1 (if it exists).
        // If y exactly matches key_y[p], arc p-1 also passes through
        // (its range ends at that boundary).  Scan backward for
        // consecutive equal key_y (NLC case).  O(1) bounded by the
        // paper's max 6 arcs at any point (§2.4 tex 144).

        if (ylo > blo) {
            std::size_t p = ylo - 1;
            result.push(p);
            // At a chord boundary (y == key_y[p]), the predecessor
            // arc's range ends at this y, so it also passes through.
            // Continue backward through any consecutive equal key_y
            // (NLC duplicates at the same y).
            if (symbolic_y_equal(arc_sequence_[p].key_symbolic_y(), y)
                && p > blo) {
                result.push(p - 1);
                // Scan further only if the predecessor also has the
                // same key_y (NLC case, §2.4 tex 144: at most 6).
                for (std::size_t i = p - 1; i > blo; --i) {
                    if (!symbolic_y_equal(
                            arc_sequence_[i].key_symbolic_y(), y))
                        break;
                    result.push(i - 1);
                }
            }
        } else if (boundary_arc == NONE) {
            // y is before the first arc's key_y and no boundary arc
            // exists → first interval arc's range extends to the edge
            // boundary, so it covers y.
            result.push(blo);
            // (If boundary_arc exists, it covers this region instead
            // and is handled below.)
        }

        // Boundary arc: if it exists, it always contains edge_idx.
        // Its y-range spans from its key_y to the first interval arc's
        // key_y.  Check whether y falls in that range.  O(1).
        if (boundary_arc != NONE) {
            SymbolicY b_y = arc_sequence_[boundary_arc].key_symbolic_y();
            if (symbolic_y_equal(b_y, y)) {
                result.push(boundary_arc);
            } else {
                // The boundary arc's range is from the edge start to
                // the first interval arc's key_y (inclusive — at the
                // boundary, both the boundary arc and the first
                // interval arc pass through the point).
                SymbolicY first_y = arc_sequence_[blo].key_symbolic_y();
                bool in_boundary = keys_ascending
                    ? symbolic_y_leq(y, first_y)
                    : symbolic_y_geq(y, first_y);
                if (in_boundary)
                    result.push(boundary_arc);
            }
        }
    };

    // [C91 §2.4]: LEFT half — ascending first_edge.
    search_half(left_begin, left_end, /*ascending=*/true);
    assert(result.count <= 3 &&
           "§2.4: at most 3 arcs per ∂C half at any point "
           "(tex 144: worst case = arc + NLC + arc at y-extremum)");

    // [C91 §2.4]: RIGHT half — descending first_edge.
    std::size_t left_count = result.count;
    search_half(right_begin, right_end, /*ascending=*/false);
    assert(result.count - left_count <= 3 &&
           "§2.4: at most 3 arcs per ∂C half at any point");
    assert(result.count <= DoubleIdentifyResult::MAX &&
           "§2.4: at most 6 arcs at any point (tex 144)");

    return result;
}

// ── Tree decomposition ──────────────────────────────────────────

void Submap::build_tree_decomposition() {
    tree_decomp_.build(*this);
}

// ── Region weight (O(1) via chord→arc adjacency) ────────────────

// [C91 §2.2] (tex 106): "the maximum number of nonnull length edges
// in any of its arcs."
// [C91 §2.4(ii)] (tex 137): chord→arc adjacency enables O(1) per
// region for conformal submaps (degree ≤ 4).
// [C91 §3.3]: granularity enforcement requires O(m) total weight
// checks, so each must be O(1).

std::size_t Submap::region_weight(std::size_t node_idx) const noexcept {
    assert(node_idx < nodes_.size() && !nodes_[node_idx].dead);

    const auto& nd = nodes_[node_idx];
    std::size_t max_count = 0;

    // Enumerate this region's arcs via chord→arc adjacency:
    // region → incident_chords → each chord's adj arcs → filter by
    // region_node.  For conformal submaps: ≤ 4 chords × 2×2 adj arcs
    // = O(1) candidates.
    auto check_adj = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            std::size_t ai = adj.arcs[k];
            assert(ai != NONE && ai < arc_sequence_.size() &&
                   "§2.4(ii): chord adj_arc must be a valid arc index");
            const auto& a = arc_sequence_[ai];
            assert(!a.dead &&
                   "§2.4(ii): chord adj_arc must be live");
            if (a.region_node == node_idx) {
                if (a.edge_count > max_count)
                    max_count = a.edge_count;
            }
        }
    };
    for (std::size_t ci : nd.incident_chords) {
        assert(ci < chords_.size());
        const auto& ch = chords_[ci];
        if (ch.dead) continue;
        check_adj(ch.left_adj);
        check_adj(ch.right_adj);
    }

    // Also check arcs at C's endpoints — these might be adjacent
    // to only one chord (or none, for single-region submaps).
    // [C91 §2.4(iii)]: if set, endpoint arc pointers must be valid.
    if (start_arc != NONE) {
        assert(start_arc < arc_sequence_.size() &&
               !arc_sequence_[start_arc].dead &&
               "§2.4(iii): start_arc must be valid and live");
        if (arc_sequence_[start_arc].region_node == node_idx) {
            if (arc_sequence_[start_arc].edge_count > max_count)
                max_count = arc_sequence_[start_arc].edge_count;
        }
    }
    if (end_arc != NONE) {
        assert(end_arc < arc_sequence_.size() &&
               !arc_sequence_[end_arc].dead &&
               "§2.4(iii): end_arc must be valid and live");
        if (arc_sequence_[end_arc].region_node == node_idx) {
            if (arc_sequence_[end_arc].edge_count > max_count)
                max_count = arc_sequence_[end_arc].edge_count;
        }
    }

    return max_count;
}

// ── Conformal / semigranular / granular ─────────────────────────

bool Submap::is_conformal() const noexcept {
    // [C91 §2.3]: "conformal submaps [are] those with node-degree
    // at most 4."
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i].dead) continue;
        if (nodes_[i].degree() > 4)
            return false;
    }
    return true;
}

bool Submap::is_semigranular(std::size_t gamma) const noexcept {
    // [C91 §2.3]: "every node of its tree has weight at most γ."
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i].dead) continue;
        if (region_weight(i) > gamma)
            return false;
    }
    return true;
}

std::size_t Submap::simulated_contraction_weight(
        std::size_t chord_idx,
        const Polygon& polygon) const noexcept {
    assert(chord_idx < chords_.size() && !chords_[chord_idx].dead);
    const auto& c = chords_[chord_idx];
    assert(c.region[0] != NONE && c.region[1] != NONE &&
           "§2.4(i): chord must have valid regions");

    // [C91 §2.3] (tex 121): "contracting any edge... produces a new
    // node whose weight."
    // [C91 §2.2] (tex 94): removal only merges at non-vertex endpoints.
    //
    // Max of all individual arcs in both regions.  Use chord→arc
    // adjacency for O(1).
    std::size_t max_count = std::max(region_weight(c.region[0]),
                                      region_weight(c.region[1]));

    // Check which endpoints are polygon vertices.
    auto endpoint_is_vertex = [&](std::size_t edge,
                                   double ey,
                                   std::size_t ey_tag) -> bool {
        assert(edge < polygon.num_edges() &&
               "§2.2: chord endpoint edge must be a valid edge index");
        const auto& e = polygon.edge(edge);
        const auto& v_start = polygon.vertex(e.start_idx);
        const auto& v_end   = polygon.vertex(e.end_idx);
        SymbolicY chord_y{ey, ey_tag};
        return symbolic_y_equal(chord_y, symbolic_y_of(v_start)) ||
               symbolic_y_equal(chord_y, symbolic_y_of(v_end));
    };

    bool left_is_vertex  = endpoint_is_vertex(c.left_edge, c.y, c.y_tag);
    bool right_is_vertex = endpoint_is_vertex(c.right_edge, c.y, c.y_tag);

    // Only compute merged weights at NON-vertex endpoints.
    auto try_merge_pair = [&](std::size_t ai, std::size_t aj) {
        assert(ai < arc_sequence_.size() && !arc_sequence_[ai].dead &&
               aj < arc_sequence_.size() && !arc_sequence_[aj].dead &&
               "§2.4(ii): adj_arcs must be valid and live");
        // [C91 §2.2] (tex 94): adjacent arcs at a chord endpoint
        // share the junction edge — same invariant as remove_chord's
        // do_merge (line 158).
        assert(arc_sequence_[ai].last_edge == arc_sequence_[aj].first_edge &&
               "§2.2: adjacent arcs at chord endpoint must share "
               "the junction edge");
        std::size_t shared_edge = arc_sequence_[aj].first_edge;
        std::size_t shared_nonnull =
            polygon.count_nonnull_edges(shared_edge, shared_edge);
        std::size_t merged = arc_sequence_[ai].edge_count
                           + arc_sequence_[aj].edge_count
                           - shared_nonnull;
        if (merged > max_count)
            max_count = merged;
    };

    if (!left_is_vertex && c.left_adj.count == 2)
        try_merge_pair(c.left_adj.arcs[0], c.left_adj.arcs[1]);
    if (!right_is_vertex && c.right_adj.count == 2)
        try_merge_pair(c.right_adj.arcs[0], c.right_adj.arcs[1]);

    return max_count;
}

bool Submap::is_granular(std::size_t gamma,
                          const Polygon& polygon) const noexcept {
    // [C91 §2.3] condition (i): all weights ≤ γ.
    if (!is_semigranular(gamma)) return false;

    // [C91 §2.3]: "by default, if (i) holds but the submap has no
    // exit chord, it is still said to be γ-granular."
    if (num_live_chords() == 0) return true;

    // [C91 §2.3] condition (ii): "contracting any edge incident upon
    // at least one node of degree less than 3 produces a new node
    // whose weight exceeds γ."
    for (std::size_t ci = 0; ci < chords_.size(); ++ci) {
        const auto& c = chords_[ci];
        if (c.dead) continue;

        std::size_t d0 = nodes_[c.region[0]].degree();
        std::size_t d1 = nodes_[c.region[1]].degree();

        if (d0 >= 3 && d1 >= 3) continue;

        if (simulated_contraction_weight(ci, polygon) <= gamma)
            return false;
    }

    // [C91 §2.3 Lemma 2.3] (tex 126): "If C is a polygonal curve with n vertices, any γ-granular conformal submap
    // of the visibility map of C has O(n/γ + 1) regions and each region is bounded by O(γ) edges."
    // This bound is conditional on γ-granularity — assert only after both conditions (i) and (ii) are confirmed.
    //
    // Mathematically bounded constant derivation:
    // 1. Edges E spanning nodes of deg < 3 is > N/2 of total tree edges (N is regions).
    // 2. Converged capacity: sum(W_u + W_v) > |E|*gamma. Max vertex duplicated overlap <= 4. Maximum graph boundary <= 4n.
    // 3. Solving: 16n >= 4*sum(W_v) >= sum(W_u+W_v) > |E|*gamma >= (N/2)*gamma => N < 32*(n/gamma).
    // Note: Mathematically this is strictly N < 32(n/gamma), but because C++ performs integer
    // division truncation on (n/gamma), we use <= 32*(n/gamma) + 32 to safely preserve the ceiling bounds.
    assert(num_live_nodes() <= 32 * (polygon.num_vertices() / (gamma > 0 ? gamma : 1)) + 32 &&
           "§2.3 Lemma 2.3: Violates strict mathematically proven N < 32(n/gamma) bound");

    return true;
}

} // namespace chazelle
