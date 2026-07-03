// src/submap/submap.cpp

#include "submap.h"
#include "../polygon/polygon.h"

#include <algorithm>

namespace chazelle {

// ── SubmapNode ──────────────────────────────────────────────────

std::size_t SubmapNode::degree() const noexcept {
    // remove_chord erases dead chord indices at removal time, so
    // incident_chords contains only live entries for live nodes.
    return incident_chords.size();
}

// ── Construction ────────────────────────────────────────────────

std::size_t Submap::add_node() {
    std::size_t idx = nodes_.size();
    nodes_.push_back(SubmapNode{});
    tree_decomp_dirty_ = true;     // [C91 §2.4(iv)]: tree decomposition indexes nodes_.
    return idx;
}

std::size_t Submap::add_arc(Arc arc) {
    // [C91 §2.4(ii) tex 137]: arc must point to a live region.
    assert(arc.region_node != NONE &&
           arc.region_node < nodes_.size() &&
           !nodes_[arc.region_node].dead &&
           "[C91 §2.4(ii) tex 137]: arc.region_node must be live");

    std::size_t idx = arc_sequence_.size();

    // [C91 §2.4 tex 133]: arc edges "in clockwise order along the double
    // boundary" — LEFT ascends within the arc, RIGHT descends.  Wrapping
    // arcs (first_side != last_side) have their legs validated inline by
    // Arc::underlying_edge_range's asserts at every call site.
    if (arc.first_side == arc.last_side) {
        if (arc.first_side == LEFT) {
            assert(arc.first_edge <= arc.last_edge &&
                   "[C91 §2.4 tex 133]: LEFT arc edges must ascend (cw on ∂C)");
        } else {
            assert(arc.first_edge >= arc.last_edge &&
                   "[C91 §2.4 tex 133]: RIGHT arc edges must descend (cw on ∂C)");
        }
    }

    // [C91 §2.4(iii) tex 138]: arc-sequence is in canonical ∂C order — LEFT
    // ascending first_edge, then RIGHT descending.  left_right_boundary_
    // tracks this implicitly, so out-of-order insertion would silently
    // break double_identify's binary search.
    if (!arc_sequence_.empty()) {
        Side prev_side = arc_sequence_.back().first_side;
        if (arc.first_side == LEFT) {
            assert(prev_side == LEFT &&
                   "[C91 §2.4(iii)]: LEFT arc after RIGHT violates ∂C order");
            assert(arc.first_edge >= arc_sequence_.back().first_edge &&
                   "[C91 §2.4(iii)]: LEFT arcs must ascend in first_edge");
        } else if (prev_side == RIGHT) {
            assert(arc.first_edge <= arc_sequence_.back().first_edge &&
                   "[C91 §2.4(iii)]: RIGHT arcs must descend in first_edge");
        }
    }

    if (arc.first_side == LEFT)
        left_right_boundary_ = idx + 1;
    arc_sequence_.push_back(arc);
    tree_decomp_dirty_ = true;
    return idx;
}

std::size_t Submap::add_chord(Chord chord) {
    std::size_t idx = chords_.size();

    // [C91 §2.4(ii) tex 137]: per-endpoint adj-arc count ∈ {1, 2}, total
    // across both endpoints ∈ {2, 3, 4}.  Polygon-vertex endpoints have
    // one adj arc; non-vertex (mid-edge) endpoints have two ([C91 §2.2 tex 94]).
    assert(chord.left_adj.count >= 1 && chord.left_adj.count <= 2 &&
           "[C91 §2.4(ii)]: LEFT endpoint adj count ∈ [1,2]");
    assert(chord.right_adj.count >= 1 && chord.right_adj.count <= 2 &&
           "[C91 §2.4(ii)]: RIGHT endpoint adj count ∈ [1,2]");
    assert(chord.left_adj.count + chord.right_adj.count >= 2 &&
           chord.left_adj.count + chord.right_adj.count <= 4 &&
           "[C91 §2.4(ii) tex 137]: total adj count ∈ [2,4]");
    auto check_adj = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            assert(adj.arcs[k] != NONE &&
                   adj.arcs[k] < arc_sequence_.size() &&
                   "[C91 §2.4(ii)]: adj_arc index must be valid");
            assert(!arc_sequence_[adj.arcs[k]].dead &&
                   "[C91 §2.4(ii)]: adj_arc must be live");
            assert((arc_sequence_[adj.arcs[k]].region_node ==
                        chord.region[0] ||
                    arc_sequence_[adj.arcs[k]].region_node ==
                        chord.region[1]) &&
                   "[C91 §2.4(ii)]: adj_arc must belong to one of the "
                   "chord's two regions");
        }
        // [C91 §2.2 tex 96]: at a mid-edge endpoint the chord splits the
        // edge's ∂C arcs into two halves — one per side of the chord; the
        // two adj arcs must therefore live in DIFFERENT regions of the
        // chord.  (Vertex endpoints aren't split per [C91 §2.2 tex 94].)
        if (adj.count == 2) {
            assert(arc_sequence_[adj.arcs[0]].region_node !=
                       arc_sequence_[adj.arcs[1]].region_node &&
                   "[C91 §2.2 tex 96]: mid-edge endpoint's two adj arcs "
                   "must lie in the chord's two distinct regions");
        }
    };
    check_adj(chord.left_adj);
    check_adj(chord.right_adj);

    // [C91 §2.2 tex 102]: dual graph is a tree — no self-loops.
    assert(chord.region[0] != chord.region[1] &&
           "[C91 §2.2 tex 102]: chord must connect two distinct regions");

    // [C91 §2 tex 47 + §2.1 tex 70]: every chord is horizontal at a polygon
    // vertex's y, so it must carry that vertex's SoS tag — otherwise
    // endpoint_is_polygon_vertex() and double_identify's y-disambiguation
    // silently misclassify.
    assert(chord.y_tag != SOS_NONE &&
           "[C91 §2 tex 47 (SoS)]: chord must carry the source vertex's SoS tag");

    // [C91 §2.1 tex 72 + §2.2 tex 108]: a null-length chord arises at a
    // y-extremum from the "inside" pair of duplicate vertices — both
    // endpoints at the same ∂C point (same edge/side/symbolic y), with
    // one adj arc per ∂C side.
    if (chord.is_null_length) {
        assert(chord.left_edge == chord.right_edge &&
               chord.left_side == chord.right_side &&
               "[C91 §2.1 tex 72]: null-length chord endpoints must coincide");
        assert(chord.left_adj.count == 1 && chord.right_adj.count == 1 &&
               "[C91 §2.1 tex 72]: null-length chord has 1 adj arc per ∂C side");
    }

    chords_.push_back(chord);

    // [C91 §2.4(i)]: update region adjacency; both regions must be live.
    for (std::size_t r : chord.region) {
        assert(r != NONE && r < nodes_.size() && !nodes_[r].dead &&
               "[C91 §2.4(i)]: chord must connect two LIVE regions");
        nodes_[r].incident_chords.push_back(idx);
    }

    tree_decomp_dirty_ = true;
    return idx;
}

// ── Chord removal (O(1) via tombstones) ─────────────────────────

// [C91 §2.2 tex 94]: "remove the chord and those endpoints that are not
// vertices of C, gluing back ∂C at those points."
// [C91 §3.3]: must be O(1) per removal to keep the submap-tree-linear bound.
// We tombstone instead of erasing; indices stay stable.  compact() strips
// dead entries before normal form.

std::size_t Submap::remove_chord(std::size_t chord_idx,
                                  const Polygon& polygon) {
    assert(chord_idx < chords_.size());
    auto& c = chords_[chord_idx];
    assert(!c.dead && "removing an already-dead chord");
    assert(c.region[0] != NONE && c.region[1] != NONE);

    std::size_t r0 = c.region[0];
    std::size_t r1 = c.region[1];
    assert(!nodes_[r0].dead && !nodes_[r1].dead);

    // Is an endpoint a polygon vertex? (Chord's symbolic y matches
    // either of the edge's endpoints.)
    auto endpoint_is_polygon_vertex = [&](std::size_t edge,
                                          double ey,
                                          std::size_t ey_tag) -> bool {
        assert(edge < polygon.num_edges() && "[C91 §2.2]: invalid edge index");
        const auto& e = polygon.edge(edge);
        SymbolicY chord_y{ey, ey_tag};
        return symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.start_idx))) ||
               symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.end_idx)));
    };

    bool left_is_vertex  = endpoint_is_polygon_vertex(c.left_edge, c.y, c.y_tag);
    bool right_is_vertex = endpoint_is_polygon_vertex(c.right_edge, c.y, c.y_tag);

    // [C91 §2.2 tex 94]: "glue back ∂C" — merge arc pairs at non-vertex
    // endpoints.  O(1): at most 2 pairs.
    auto do_merge = [&](std::size_t ai, std::size_t aj) {
        assert(ai != NONE && ai < arc_sequence_.size() &&
               aj != NONE && aj < arc_sequence_.size() &&
               !arc_sequence_[ai].dead && !arc_sequence_[aj].dead &&
               "[C91 §2.4(ii)]: adj_arc must be valid + live");

        // [C91 §2.2 tex 94 + §2.1 tex 72]: the two arcs glued at an
        // endpoint are always distinct.  ai == aj would require this
        // chord's two endpoints to be the ONLY subdivision points of ∂C
        // (both regions leaf-with-single-arc), which cannot happen in a
        // submap of V(C): C's endpoints give rise to ∂C vertices (tex 72
        // case 3) that are vertices of C, and removals never glue ∂C at
        // vertices of C (tex 94) — so the arc boundaries at C's endpoints
        // persist under every chord removal.
        assert(ai != aj &&
               "[C91 §2.2 tex 94]: glued arc pair must be distinct");

        auto& a_keep = arc_sequence_[ai];
        auto& a_dead = arc_sequence_[aj];

        assert(a_keep.last_edge == a_dead.first_edge &&
               "[C91 §2.4(iii)]: adj arcs at chord endpoint share the junction edge");

        a_keep.last_edge = a_dead.last_edge;
        a_keep.last_side = a_dead.last_side;
        // [C91 §2.2 tex 106]: edge_count = nonnull P-edges in the arc's
        // underlying range (single-count per "distinct vertices of C" in
        // the Lemma 2.3 proof, tex 129).  Recompute rather than use
        // `a+b−shared`: when either input wraps ([C91 §2.4 tex 142]) its
        // range extends to the C-endpoint, so the inputs overlap on more
        // than just the boundary edge and the additive formula over-counts.
        auto [lo, hi] = a_keep.underlying_edge_range(start_vertex, end_vertex);
        a_keep.edge_count = polygon.count_nonnull_edges(lo, hi);

        // Tombstone the dead arc.
        a_dead.dead = true;
        compacted_ = false;  // arc-sequence has dead entries now

        // [C91 §2.4(iii)]: If the dead arc was start_arc or end_arc,
        // redirect to the surviving arc (which now covers the merged
        // range including the endpoint).
        if (start_arc == aj) start_arc = ai;
        if (end_arc == aj)   end_arc = ai;

        // Update adj arcs of chords incident on either region to point to
        // the surviving arc.  O(degree) = O(1) for conformal submaps.
        //
        // The removed chord itself must NOT be skipped here: when a leaf
        // region is bounded by this chord and a single arc (the degree-1
        // node targeted by [C91 §2.3 tex 121]'s contraction), that arc is
        // adjacent to the chord at BOTH endpoints.  If it is killed by
        // the first glue, the second endpoint's pair must be rewritten to
        // the survivor so that "glueing back ∂C at those points"
        // ([C91 §2.2 tex 94]) chains all three arcs into one.
        auto replace_arc = [&](Chord::AdjArcs& adj) {
            for (std::size_t k = 0; k < adj.count; ++k)
                if (adj.arcs[k] == aj) adj.arcs[k] = ai;
        };
        for (std::size_t ri : {r0, r1}) {
            for (std::size_t ci : nodes_[ri].incident_chords) {
                auto& other = chords_[ci];
                if (other.dead) continue;
                replace_arc(other.left_adj);
                replace_arc(other.right_adj);
            }
        }
    };

    // [C91 §2.2 tex 108]: "once removed, a chord of zero length
    // ceases to separate any arcs."
    //
    // Use the stored flag — recomputing from edge/side would misclassify
    // exit chords whose left_edge == right_edge but sides differ.
    bool is_null_length = c.is_null_length;

    if (is_null_length) {
        assert(c.left_adj.count == 1 && c.right_adj.count == 1 &&
               "[C91 §2.2]: null-length chord must have 1 adj arc per side");
        // [C91 §2.1 tex 72]: null-length chord endpoints are y-extremum
        // duplicate vertices — both polygon vertices by construction.
        assert(left_is_vertex && right_is_vertex &&
               "[C91 §2.1 tex 72]: null-length chord endpoints are polygon vertices");
        // Vertex endpoints don't trigger merging ([C91 §2.2 tex 94]).  The
        // null-length arc is absorbed into r0 by the reassignment below.
    } else {
        // [C91 §2.2 tex 94]: non-vertex endpoint ⟹ 2 adj arcs to glue.
        if (!left_is_vertex) {
            assert(c.left_adj.count == 2 &&
                   "[C91 §2.2 tex 94]: non-vertex endpoint needs 2 adj arcs");
            do_merge(c.left_adj.arcs[0], c.left_adj.arcs[1]);
        }
        if (!right_is_vertex) {
            assert(c.right_adj.count == 2 &&
                   "[C91 §2.2 tex 94]: non-vertex endpoint needs 2 adj arcs");
            do_merge(c.right_adj.arcs[0], c.right_adj.arcs[1]);
        }
    }

    // Reassign r1's arcs to r0 by walking r1's incident chords' adj arcs
    // ([C91 §2.4(ii) tex 137]).  All entries are live: do_merge's
    // replace_arc rewrote every reference to a killed arc — including
    // the removed chord's own slots — to the surviving arc.
    auto reassign_live = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            std::size_t ai = adj.arcs[k];
            assert(ai != NONE && ai < arc_sequence_.size() &&
                   !arc_sequence_[ai].dead &&
                   "[C91 §2.4(ii)]: adj arcs must be live after glueing");
            if (arc_sequence_[ai].region_node == r1)
                arc_sequence_[ai].region_node = r0;
        }
    };
    for (std::size_t ci : nodes_[r1].incident_chords) {
        if (ci == chord_idx) continue;
        const auto& ch = chords_[ci];
        if (ch.dead) continue;
        reassign_live(ch.left_adj);
        reassign_live(ch.right_adj);
    }
    reassign_live(c.left_adj);
    reassign_live(c.right_adj);

    // Move r1's other incident chords to r0, and rewrite their region
    // pointers from r1 → r0.
    for (std::size_t ci : nodes_[r1].incident_chords) {
        if (ci == chord_idx) continue;
        auto& ch = chords_[ci];
        if (ch.dead) continue;
        nodes_[r0].incident_chords.push_back(ci);
        if (ch.region[0] == r1) ch.region[0] = r0;
        if (ch.region[1] == r1) ch.region[1] = r0;
    }

    // Drop this chord from r0's incident list.  O(degree) = O(1).
    {
        auto& ic = nodes_[r0].incident_chords;
        ic.erase(std::remove(ic.begin(), ic.end(), chord_idx), ic.end());
    }

    c.dead = true;
    nodes_[r1].dead = true;
    tree_decomp_dirty_ = true;
    return r0;
}

// ── Tombstone compaction ────────────────────────────────────────

// [C91 §3.3]: "We can now put S in normal form."  compact() strips
// dead entries and rebuilds index mappings in O(m).

void Submap::compact() {
    // Cascaded removals can leave live arcs pointing to dead regions:
    // when r2 is absorbed into r1 (chord A removed), then r1 into r0
    // (chord B), arcs that lived in r2 were reassigned to r1 by A, but B
    // only walks chord-adjacency and misses those orphans.
    //
    // Fix: build a forwarding table (each dead chord records region[1]
    // → region[0]), resolve chains, fixup arcs.  O(m), no asymptotic hit.
    {
        std::vector<std::size_t> forward(nodes_.size(), NONE);
        for (const auto& ch : chords_) {
            if (!ch.dead) continue;
            // remove_chord always kills region[1], keeps region[0].
            std::size_t dead_r = ch.region[1];
            std::size_t live_r = ch.region[0];
            if (dead_r < nodes_.size() && nodes_[dead_r].dead &&
                forward[dead_r] == NONE) {
                forward[dead_r] = live_r;
            }
        }
        auto resolve = [&](std::size_t r) -> std::size_t {
            std::size_t root = r;
            while (root < nodes_.size() && nodes_[root].dead &&
                   forward[root] != NONE)
                root = forward[root];
            // Path-compress.
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

    compacted_ = true;          // [C91 §2.4 tex 144]: no dead arcs remain.
    tree_decomp_dirty_ = true;  // [C91 §2.4(iv)]: tree decomposition's indices no longer match.
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
    // [C91 §2.2 tex 96]: a submap is a polygonal subdivision of the
    // spherical plane — always ≥1 region by definition.
    // Asserted first so an empty submap fails with this diagnostic rather
    // than the misleading "0 ≠ 1" tree-property message below.
    assert(num_live_nodes() >= 1 &&
           "[C91 §2.2]: submap must have at least one region");

    // [C91 §2.2]: "the dual graph of a submap is itself a tree."
    assert(num_live_nodes() == num_live_chords() + 1 &&
           "[C91 §2.2]: submap tree property: num_regions = num_chords + 1");
}

void Submap::check_invariants() const {
    assert_tree_property();

    // Every live chord's regions must be live nodes.
    for (std::size_t i = 0; i < chords_.size(); ++i) {
        const auto& c = chords_[i];
        if (c.dead) continue;
        assert(c.region[0] < nodes_.size() && !nodes_[c.region[0]].dead &&
               "[C91 §2.2]: chord region[0] invalid or dead");
        assert(c.region[1] < nodes_.size() && !nodes_[c.region[1]].dead &&
               "[C91 §2.2]: chord region[1] invalid or dead");
    }

    // Every live arc's region_node must be a live node.
    for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
        if (arc_sequence_[i].dead) continue;
        assert(arc_sequence_[i].region_node < nodes_.size() &&
               !nodes_[arc_sequence_[i].region_node].dead &&
               "[C91 §2.2]: arc region_node invalid or dead");
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
                       "[C91 §2.4(ii)]: adj_arc index invalid or dead");
                const auto& a = arc_sequence_[ai];
                assert((a.region_node == c.region[0] ||
                        a.region_node == c.region[1]) &&
                       "[C91 §2.4(ii)]: adj_arc must belong to one of "
                       "the chord's endpoint regions");
            }
        };
        check_adj(c.left_adj);
        check_adj(c.right_adj);
    }

    // [C91 §2.4(iii) tex 138]: ∂C ordering of live arcs (skip dead gaps).
    {
        bool seen_right = false;
        std::size_t prev_live = NONE;
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            if (arc_sequence_[i].dead) continue;
            if (arc_sequence_[i].first_side == RIGHT) seen_right = true;
            if (seen_right) {
                assert(arc_sequence_[i].first_side == RIGHT &&
                       "[C91 §2.4(iii)]: LEFT arc after RIGHT arc in "
                       "arc-sequence table violates ∂C order");
            }
            if (prev_live != NONE) {
                if (arc_sequence_[prev_live].first_side == LEFT &&
                    arc_sequence_[i].first_side == LEFT) {
                    assert(arc_sequence_[i].first_edge >=
                           arc_sequence_[prev_live].first_edge &&
                           "[C91 §2.4(iii)]: LEFT arcs must have ascending "
                           "first_edge for double_identify binary search");
                }
                if (arc_sequence_[prev_live].first_side == RIGHT &&
                    arc_sequence_[i].first_side == RIGHT) {
                    assert(arc_sequence_[i].first_edge <=
                           arc_sequence_[prev_live].first_edge &&
                           "[C91 §2.4(iii)]: RIGHT arcs must have descending "
                           "first_edge for double_identify binary search");
                }
            }
            prev_live = i;
        }
    }

    // [C91 §2.4(iii) tex 138]: endpoint arc pointers.
    if (num_live_arcs() > 0) {
        assert(start_arc != NONE &&
               "[C91 §2.4(iii)]: start_arc must be set when arcs exist");
        assert(end_arc != NONE &&
               "[C91 §2.4(iii)]: end_arc must be set when arcs exist");
        assert(start_arc < arc_sequence_.size() &&
               !arc_sequence_[start_arc].dead &&
               "[C91 §2.4(iii)]: start_arc out of range or dead");
        assert(end_arc < arc_sequence_.size() &&
               !arc_sequence_[end_arc].dead &&
               "[C91 §2.4(iii)]: end_arc out of range or dead");

        // [C91 §2.4(iii) tex 138 + tex 144]: canonical layout has LEFT arcs
        // first (ascending), then RIGHT (descending).  The endpoints sit
        // at the LEFT-half boundaries — required by double_identify's
        // two-linear-sequence binary search.
        assert(arc_sequence_[start_arc].first_side == LEFT &&
               "[C91 §2.4(iii)]: start_arc must be LEFT-side");
        assert(arc_sequence_[end_arc].first_side == LEFT &&
               "[C91 §2.4(iii)]: end_arc must be LEFT-side");
        assert(start_arc == 0 &&
               "[C91 §2.4(iii)]: start_arc must be index 0");
        assert(end_arc + 1 == left_right_boundary_ &&
               "[C91 §2.4 tex 144]: end_arc must be the last LEFT arc");

        assert(start_vertex != NONE &&
               "[C91 §2.4(iii)]: start_vertex required when arcs exist");
        {
            // [C91 §2.4(iii)]: start_arc must pass through C's start vertex.
            // Vertex v lies on arc [first_edge, last_edge] iff
            //   first_edge ≤ v ≤ last_edge + 1 (non-wrapped).
            // For wrapped arcs, the range extends through the C endpoint.
            const auto& sa = arc_sequence_[start_arc];
            std::size_t elo = std::min(sa.first_edge, sa.last_edge);
            std::size_t ehi = std::max(sa.first_edge, sa.last_edge) + 1;
            if (sa.first_side != sa.last_side)
                elo = std::min(elo, start_vertex);
            assert(start_vertex >= elo && start_vertex <= ehi &&
                   "[C91 §2.4(iii)]: start_arc must pass through start_vertex");
        }

        assert(end_vertex != NONE && end_vertex > 0 &&
               "[C91 §2.4(iii)]: end_vertex required when arcs exist");
        {
            // [C91 §2.4(iii)]: end_arc must pass through C's end vertex.
            // C's last edge is end_vertex - 1.
            const auto& ea = arc_sequence_[end_arc];
            std::size_t c_end_edge = end_vertex - 1;
            std::size_t elo = std::min(ea.first_edge, ea.last_edge);
            std::size_t ehi = std::max(ea.first_edge, ea.last_edge);
            if (ea.first_side != ea.last_side)
                ehi = std::max(ehi, c_end_edge);
            assert(c_end_edge >= elo && c_end_edge <= ehi &&
                   "[C91 §2.4(iii)]: end_arc must pass through end_vertex");
        }
    }
}

void Submap::check_invariants(const Polygon& polygon) const {
    check_invariants();

    // [C91 §2.4 tex 144]: within a same-(first_side, first_edge) run, start-y
    // must be monotonic — double_identify Phase 2 infers the direction
    // from the run's endpoints and binary-searches.  We check
    // monotonicity in the inferred direction (more permissive than strict
    // canonical ascent/descent: at a null-length chord's duplicate-vertex pair the SoS
    // tag direction may oppose geometric edge ascent, and Phase 2 only
    // needs SOME monotonic direction).
    {
        std::size_t i = 0;
        while (i < arc_sequence_.size()) {
            if (arc_sequence_[i].dead) { ++i; continue; }
            // Find end of the maximal live (first_side, first_edge) run.
            std::size_t j = i;
            while (true) {
                std::size_t next = j + 1;
                while (next < arc_sequence_.size() &&
                       arc_sequence_[next].dead) ++next;
                if (next >= arc_sequence_.size()) break;
                if (arc_sequence_[next].first_side !=
                        arc_sequence_[i].first_side ||
                    arc_sequence_[next].first_edge !=
                        arc_sequence_[i].first_edge) break;
                j = next;
            }
            // Run is arc-sequence indices in [i, j] (live, possibly
            // with dead gaps inside).  Direction is inferred from
            // first vs last live element.
            if (j > i) {
                bool asc = symbolic_y_leq(
                    arc_start_symbolic_y(i, polygon),
                    arc_start_symbolic_y(j, polygon));
                std::size_t prev_live = i;
                for (std::size_t k = i + 1; k <= j; ++k) {
                    if (arc_sequence_[k].dead) continue;
                    if (asc) {
                        assert(symbolic_y_leq(
                                arc_start_symbolic_y(prev_live, polygon),
                                arc_start_symbolic_y(k, polygon)) &&
                               "[C91 §2.4 tex 144]: same-first_edge run must "
                               "be start-y-monotonic (ascending)");
                    } else {
                        assert(symbolic_y_geq(
                                arc_start_symbolic_y(prev_live, polygon),
                                arc_start_symbolic_y(k, polygon)) &&
                               "[C91 §2.4 tex 144]: same-first_edge run must "
                               "be start-y-monotonic (descending)");
                    }
                    prev_live = k;
                }
            }
            i = j + 1;
        }
    }

    // [C91 §2.2 tex 106]: arc.edge_count caches the max-nonnull-edges count
    // used by region_weight and simulated_contraction_weight.  Validate
    // it against polygon.count_nonnull_edges — a wrong cached value silently
    // miscomputes weights, breaking γ-granularity decisions.
    //
    // [C91 §2.4 tex 133]: null-length arcs encode a single ∂C point (the
    // null-length chord's duplicate-vertex pair), not an edge span.
    // Under [C91 §2.1] every
    // polygon edge has nonzero length, so edge_count == 0 uniquely
    // identifies null-length arcs — skip the cache check for them.
    //
    // Wrapped arcs: underlying_edge_range (arc.h) covers both legs as a
    // contiguous range.  Wrapped arcs don't appear in canonical normal
    // form ([C91 §3.0(ii)](2) "no double-backing"); only start_arc / end_arc
    // may wrap, counted via the same union range.
    if (num_live_arcs() > 0) {
        assert(start_vertex != NONE && end_vertex != NONE &&
               end_vertex > 0 &&
               "[C91 §2.4(iii)]: start/end_vertex must be set when arcs exist");
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            const auto& a = arc_sequence_[i];
            if (a.dead) continue;
            // [C91 §2 tex 47]: arc start position derived from the
            // bounding chord (or C-endpoint) — both carry SoS tags.
            assert(arc_start_symbolic_y(i, polygon).tag != SOS_NONE &&
                   "[C91 §2 tex 47]: live arc's start position needs SoS tag");
            if (a.edge_count == 0) continue;     // null-length arc

            assert(a.first_edge < polygon.num_edges() &&
                   a.last_edge < polygon.num_edges() &&
                   "[C91 §2.4(iii)]: arc edges must be valid input-table indices");
            auto [lo, hi] = a.underlying_edge_range(start_vertex, end_vertex);
            std::size_t actual = polygon.count_nonnull_edges(lo, hi);
            assert(a.edge_count == actual &&
                   "[C91 §2.2 tex 106]: arc.edge_count cache must match "
                   "polygon.count_nonnull_edges over the arc's "
                   "underlying edge range");
        }
    }

    // [C91 §2.2 tex 94] consistency: remove_chord and simulated_contraction_
    // weight classify a chord endpoint as a polygon vertex iff its
    // symbolic y matches one of the underlying edge's vertices.  Vertex
    // endpoints have one adj arc (no glueing); non-vertex endpoints have
    // two (glued).  A tag/count mismatch silently breaks both directions:
    //   one adj arc with no vertex match → contraction skips a needed merge.
    //   two adj arcs with a vertex match → contraction merges arcs at a vertex.
    auto matches_an_endpoint = [&](std::size_t edge_idx,
                                    SymbolicY chord_y) -> bool {
        assert(edge_idx < polygon.num_edges());
        const auto& e = polygon.edge(edge_idx);
        return symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.start_idx))) ||
               symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.end_idx)));
    };
    for (std::size_t ci = 0; ci < chords_.size(); ++ci) {
        const auto& c = chords_[ci];
        if (c.dead) continue;
        // Null-length chords use auxiliary SoS tags (beyond polygon-vertex
        // indices) to disambiguate multiple null-length chords at the
        // same y-extremum; add_chord enforces their structural invariants
        // separately.
        if (c.is_null_length) continue;

        // [C91 §2.1 tex 70]: a non-null-length chord is horizontal at its
        // source vertex's y, with y_tag = that vertex's index.
        assert(c.y_tag < polygon.num_vertices() &&
               "[C91 §2.1 tex 70]: exit chord y_tag must be a polygon vertex index");
        SymbolicY chord_y{c.y, c.y_tag};
        assert(symbolic_y_equal(chord_y,
                                symbolic_y_of(polygon.vertex(c.y_tag))) &&
               "[C91 §2.1 tex 70]: exit chord must be horizontal at its source vertex");

        assert((c.left_adj.count == 1) == matches_an_endpoint(c.left_edge, chord_y) &&
               "[C91 §2.2 tex 94]: LEFT endpoint count == 1 ⟺ endpoint is a polygon vertex");
        assert((c.right_adj.count == 1) == matches_an_endpoint(c.right_edge, chord_y) &&
               "[C91 §2.2 tex 94]: RIGHT endpoint count == 1 ⟺ endpoint is a polygon vertex");

        // [C91 §2.2 tex 94]: removing the chord doesn't glue ∂C at a vertex
        // endpoint — meaning the single adj arc on that side already spans
        // the polygon vertex c.y_tag (no split there).  Verify the vertex
        // lies in the adj arc's edge range.
        auto arc_spans_vertex = [&](std::size_t arc_idx) -> bool {
            const auto& a = arc_sequence_[arc_idx];
            std::size_t elo = std::min(a.first_edge, a.last_edge);
            std::size_t ehi = std::max(a.first_edge, a.last_edge) + 1;
            if (a.first_side != a.last_side) {
                if (a.first_side == LEFT) ehi = std::max(ehi, end_vertex);
                else                       elo = std::min(elo, start_vertex);
            }
            return c.y_tag >= elo && c.y_tag <= ehi;
        };
        if (c.left_adj.count == 1) {
            assert(arc_spans_vertex(c.left_adj.arcs[0]) &&
                   "[C91 §2.2 tex 94]: LEFT vertex endpoint's adj arc must span the polygon vertex");
        }
        if (c.right_adj.count == 1) {
            assert(arc_spans_vertex(c.right_adj.arcs[0]) &&
                   "[C91 §2.2 tex 94]: RIGHT vertex endpoint's adj arc must span the polygon vertex");
        }
    }
}

// ── Double identification ───────────────────────────────────────

Submap::DoubleIdentifyResult
Submap::double_identify(std::size_t edge_idx, SymbolicY y,
                         const Polygon& polygon) const {
    DoubleIdentifyResult result;
    if (arc_sequence_.empty()) return result;

    // [C91 §2.4 tex 144]: O(log m) requires a compacted arc-sequence (no
    // dead gaps).  O(1) flag — remove_chord clears it, compact() sets it.
    assert(compacted_ &&
           "[C91 §2.4]: double_identify requires compacted arc-sequence");

    assert(y.tag != SOS_NONE &&
           "[C91 §2.4]: double_identify requires a valid SoS y-tag");

    // [C91 §2.4 tex 144]: "break up the circular arc sequence into two linear
    // sequences" anchored at C's endpoints.  Normal-form layout:
    //   LEFT  half: [start_arc=0, ..., end_arc] = [0, lrb)
    //   RIGHT half: [end_arc+1, ..., end)       = [lrb, end)
    assert(start_arc == 0 && "[C91 §2.4]: start_arc must be index 0");
    assert(end_arc == left_right_boundary_ - 1 &&
           "[C91 §2.4 tex 144]: end_arc must be the last LEFT arc");

    std::size_t left_begin  = start_arc;
    std::size_t left_end    = end_arc + 1;
    std::size_t right_begin = left_end;
    std::size_t right_end   = arc_sequence_.size();

    // [C91 §2.4 tex 144]: two-phase binary search.
    //   Phase 1: bsearch by first_edge → contiguous run of arcs on edge_idx.
    //   Phase 2: bsearch by start-y within that run → arc(s) at y.
    //   Plus O(1) check at run's left neighbor for a boundary arc with
    //   first_edge < edge_idx but last_edge ≥ edge_idx.
    // Total: O(log m).
    auto search_half = [&](std::size_t lo, std::size_t hi, bool ascending) {
        if (lo >= hi) return;

        // Phase 1: bsearch for blo = first arc with first_edge ≥ edge_idx
        // (LEFT half ascending; RIGHT half descending uses the opposite).
        std::size_t blo = lo, bhi = hi;
        while (blo < bhi) {
            std::size_t mid = blo + (bhi - blo) / 2;
            bool advance = ascending
                ? (arc_sequence_[mid].first_edge < edge_idx)
                : (arc_sequence_[mid].first_edge > edge_idx);
            if (advance) blo = mid + 1; else bhi = mid;
        }

        // Find run end: bend = first arc with first_edge != edge_idx.
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

        // [C91 §2.4 tex 142]: boundary arc at blo-1 may cover edge_idx via
        // double-backing.  Only start_arc and end_arc may wrap in a
        // normal-form submap ([C91 §3.0(ii)](2) tex 170 + [C91 §2.4(iii) tex 138]).
        auto arc_contains_edge = [&](std::size_t ai) -> bool {
            assert(ai < arc_sequence_.size() && "[C91 §2.4]: invalid arc index");
            const auto& a = arc_sequence_[ai];
            assert(!a.dead && "[C91 §2.4]: arc_contains_edge on dead arc");
            std::size_t elo = std::min(a.first_edge, a.last_edge);
            std::size_t ehi = std::max(a.first_edge, a.last_edge);
            if (a.first_side != a.last_side) {
                assert((ai == start_arc || ai == end_arc) &&
                       "[C91 §3.0(ii)(2)/§2.4(iii)]: only start/end_arc may "
                       "double-back in a normal-form submap");
                if (ai == start_arc)
                    elo = std::min(elo, start_vertex);
                if (ai == end_arc) {
                    assert(end_vertex > 0 && "[C91 §2.1]: end_vertex ≥ 1");
                    ehi = std::max(ehi, end_vertex - 1);
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
        // Two arcs (one interval + one boundary): can't infer start-y
        // direction from a single interval arc, so derive it from the
        // edge's geometric y-direction.  The interval arc's start-y is
        // the shared junction; together they cover the whole edge.
        if (interval_len == 1 && boundary_arc != NONE) {
            SymbolicY junction_y = arc_start_symbolic_y(blo, polygon);
            if (symbolic_y_equal(junction_y, y)) {
                result.push(blo);
                result.push(boundary_arc);
            } else {
                // boundary_arc is BEFORE blo in ∂C traversal.
                // LEFT half: traversal = edge direction; RIGHT half: reversed.
                assert(edge_idx < polygon.num_edges());
                const auto& e = polygon.edge(edge_idx);
                bool edge_ascending = symbolic_y_less(
                    symbolic_y_of(polygon.vertex(e.start_idx)),
                    symbolic_y_of(polygon.vertex(e.end_idx)));
                bool traversal_ascending = ascending ? edge_ascending
                                                     : !edge_ascending;
                bool y_in_boundary = traversal_ascending
                    ? symbolic_y_less(y, junction_y)
                    : symbolic_y_greater(y, junction_y);
                result.push(y_in_boundary ? boundary_arc : blo);
            }
            return;
        }

        // [C91 §2.4 tex 144] Phase 2: bsearch by start-y.  Direction is
        // inferred from the run's endpoints (depends on the edge's
        // geometric y-direction per tex 138, NOT on the LEFT/RIGHT half).
        bool keys_ascending = (interval_len >= 2) &&
            symbolic_y_leq(arc_start_symbolic_y(blo, polygon),
                           arc_start_symbolic_y(bend - 1, polygon));

        std::size_t ylo = blo, yhi = bend;
        while (ylo < yhi) {
            std::size_t mid = ylo + (yhi - ylo) / 2;
            SymbolicY mid_y = arc_start_symbolic_y(mid, polygon);
            if (keys_ascending) {
                // Find last arc with start-y <= y: go right if mid <= y.
                if (symbolic_y_leq(mid_y, y))
                    ylo = mid + 1;
                else
                    yhi = mid;
            } else {
                // Find last arc with start-y >= y: go right if mid >= y.
                if (symbolic_y_geq(mid_y, y))
                    ylo = mid + 1;
                else
                    yhi = mid;
            }
        }
        // ylo = first arc past y; the arc containing y is at p = ylo-1.
        // If y exactly matches arc_start_y[p], arc p-1 also passes through
        // (its range ends at that boundary); walk back through equal-y
        // runs to pick up null-length chord duplicates.  Bounded by 6
        // ([C91 §2.4 tex 144]).
        if (ylo > blo) {
            std::size_t p = ylo - 1;
            result.push(p);
            if (symbolic_y_equal(arc_start_symbolic_y(p, polygon), y)
                && p > blo) {
                result.push(p - 1);
                for (std::size_t i = p - 1; i > blo; --i) {
                    // [C91 §2.4 tex 144]: "at most six of them" — early
                    // break keeps the fixed-size push safe in release builds
                    // (DoubleIdentifyResult::push asserts < MAX only in debug).
                    if (result.count >= DoubleIdentifyResult::MAX) break;
                    if (!symbolic_y_equal(
                            arc_start_symbolic_y(i, polygon), y))
                        break;
                    result.push(i - 1);
                }
            }
        } else if (boundary_arc == NONE) {
            // y is before the first arc's start-y; with no boundary arc,
            // the first interval arc's range extends here.
            result.push(blo);
        }

        // Boundary arc (if present) always covers edge_idx; its y-range
        // runs from its start-y to the first interval arc's start-y.
        if (boundary_arc != NONE) {
            SymbolicY b_y = arc_start_symbolic_y(boundary_arc, polygon);
            if (symbolic_y_equal(b_y, y)) {
                result.push(boundary_arc);
            } else {
                // The boundary arc's range extends from the edge start
                // to the first interval arc's start-y (inclusive — both
                // arcs pass through the boundary point).
                SymbolicY first_y = arc_start_symbolic_y(blo, polygon);
                bool in_boundary = keys_ascending
                    ? symbolic_y_leq(y, first_y)
                    : symbolic_y_geq(y, first_y);
                if (in_boundary)
                    result.push(boundary_arc);
            }
        }
    };

    // [C91 §2.4 tex 144]: at most 3 arcs per ∂C half (= arc + null-length
    // chord + arc at a y-extremum); total ≤ 6 across both halves.
    search_half(left_begin, left_end, /*ascending=*/true);
    assert(result.count <= 3 && "[C91 §2.4 tex 144]: ≤ 3 arcs per ∂C half");

    std::size_t left_count = result.count;
    search_half(right_begin, right_end, /*ascending=*/false);
    assert(result.count - left_count <= 3 && "[C91 §2.4 tex 144]: ≤ 3 arcs per ∂C half");
    assert(result.count <= DoubleIdentifyResult::MAX &&
           "[C91 §2.4 tex 144]: ≤ 6 arcs at any point");

    return result;
}

// ── Arc start symbolic y (derived from bounding chord per tex 133) ──

SymbolicY Submap::arc_start_symbolic_y(std::size_t arc_idx,
                                       const Polygon& polygon) const {
    assert(arc_idx < arc_sequence_.size() &&
           !arc_sequence_[arc_idx].dead &&
           "[C91 §2.4]: arc index must be valid + live");
    const Arc& a = arc_sequence_[arc_idx];

    // [C91 §2.4(ii) tex 137]: walk the arc's region's incident chords.
    for (std::size_t ci : nodes_[a.region_node].incident_chords) {
        const Chord& c = chords_[ci];
        if (c.dead) continue;
        if (c.is_null_length) {
            // [C91 §2.1 tex 72]: null-length chord has count==1 per side;
            // R_inner contains a single null-arc whose start is the
            // chord's y.  Either slot may hold the inner null-arc per
            // rebuild_submap's endpoint-sweep order.
            if (c.right_adj.count == 1 &&
                c.right_adj.arcs[0] == arc_idx &&
                a.region_node == c.region[1])
                return c.symbolic_y();
            if (c.left_adj.count == 1 &&
                c.left_adj.arcs[0] == arc_idx &&
                a.region_node == c.region[1])
                return c.symbolic_y();
            continue;
        }
        // [C91 §2.4 tex 137 + tex 133]: at a count==2 adj pair (mid-edge
        // split per §2.2 tex 96, or vertex endpoint with both pre/post
        // arcs listed), arcs[1] is the arc starting at the chord.
        if (c.left_adj.count == 2 && c.left_adj.arcs[1] == arc_idx)
            return c.symbolic_y();
        if (c.right_adj.count == 2 && c.right_adj.arcs[1] == arc_idx)
            return c.symbolic_y();
    }

    // [C91 §2.1 tex 72]: arc-endpoints not at chords sit at polygon
    // vertices — the outside-pair companions of an interior y-extremum
    // (tex 72 case 2), the C-endpoints' companions (case 3), or any
    // junction vertex between two same-region arcs after a non-merging
    // chord removal ([C91 §2.2 tex 94]).  The polygon's input table
    // ([C91 §2.4(iii) tex 138]) carries the y directly: edge first_edge
    // starts at vertex(first_edge) under LEFT traversal (ascending) and
    // at vertex(first_edge+1) under RIGHT traversal (descending).
    assert(a.first_edge < polygon.num_edges() &&
           "[C91 §2.4(iii)]: first_edge must be a valid edge index");
    std::size_t vidx = (a.first_side == LEFT)
        ? a.first_edge
        : a.first_edge + 1;
    assert(vidx < polygon.num_vertices() &&
           "[C91 §2.4(iii)]: arc-start polygon vertex must be valid");
    return symbolic_y_of(polygon.vertex(vidx));
}

// ── Tree decomposition ──────────────────────────────────────────

void Submap::build_tree_decomposition() {
    // If a previous tree decomposition was flagged dirty by a mutator,
    // build() throws away the old contents and rebuilds from scratch.
    tree_decomp_.build(*this);
    tree_decomp_dirty_ = false;
}

// ── Region weight ───────────────────────────────────────────────

// [C91 §2.2 tex 106]: weight = max nonnull-edge count over the region's arcs.
// O(1) per region via chord→arc adjacency ([C91 §2.4(ii) tex 137]) for conformal
// submaps (degree ≤ 4) — required by [C91 §3.3] granularity enforcement.

std::size_t Submap::region_weight(std::size_t node_idx) const noexcept {
    assert(node_idx < nodes_.size() && !nodes_[node_idx].dead);

    const auto& nd = nodes_[node_idx];
    std::size_t max_count = 0;

    // ≤ 4 incident chords × 4 adj arcs each = O(1) candidates.
    auto check_adj = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            std::size_t ai = adj.arcs[k];
            assert(ai != NONE && ai < arc_sequence_.size() &&
                   !arc_sequence_[ai].dead &&
                   "[C91 §2.4(ii)]: adj_arc must be valid + live");
            const auto& a = arc_sequence_[ai];
            if (a.region_node == node_idx && a.edge_count > max_count)
                max_count = a.edge_count;
        }
    };
    for (std::size_t ci : nd.incident_chords) {
        assert(ci < chords_.size());
        const auto& ch = chords_[ci];
        if (ch.dead) continue;
        check_adj(ch.left_adj);
        check_adj(ch.right_adj);
    }

    // [C91 §2.4(iii)]: single-region (no-chord) submap — the chord-adj
    // loop above is empty.  In multi-region submaps every arc is adj to
    // some chord, so start_arc/end_arc are already covered there.
    if (num_live_chords() == 0) {
        if (start_arc != NONE) {
            assert(start_arc < arc_sequence_.size() &&
                   !arc_sequence_[start_arc].dead &&
                   "[C91 §2.4(iii)]: start_arc must be valid and live");
            if (arc_sequence_[start_arc].region_node == node_idx &&
                arc_sequence_[start_arc].edge_count > max_count)
                max_count = arc_sequence_[start_arc].edge_count;
        }
        if (end_arc != NONE) {
            assert(end_arc < arc_sequence_.size() &&
                   !arc_sequence_[end_arc].dead &&
                   "[C91 §2.4(iii)]: end_arc must be valid and live");
            if (arc_sequence_[end_arc].region_node == node_idx &&
                arc_sequence_[end_arc].edge_count > max_count)
                max_count = arc_sequence_[end_arc].edge_count;
        }
    }

    return max_count;
}

// ── Conformal / semigranular / granular ─────────────────────────

bool Submap::is_conformal() const noexcept {
    // [C91 §2.3]: conformal = node-degree ≤ 4.
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i].dead) continue;
        if (nodes_[i].degree() > 4) return false;
    }
    return true;
}

bool Submap::is_semigranular(std::size_t gamma) const noexcept {
    // [C91 §2.3]: every region weight ≤ γ.
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i].dead) continue;
        if (region_weight(i) > gamma) return false;
    }
    return true;
}

std::size_t Submap::simulated_contraction_weight(
        std::size_t chord_idx,
        const Polygon& polygon) const noexcept {
    assert(chord_idx < chords_.size() && !chords_[chord_idx].dead);
    const auto& c = chords_[chord_idx];
    assert(c.region[0] != NONE && c.region[1] != NONE &&
           c.region[0] != c.region[1] &&
           "[C91 §2.4(i)/§2.2 tex 102]: chord regions valid + distinct (tree)");

    // [C91 §2.3 tex 121 + §2.2 tex 94]: merged region's weight = max over
    // both regions' arcs, plus merged pairs at non-vertex endpoints only.
    std::size_t max_count = std::max(region_weight(c.region[0]),
                                      region_weight(c.region[1]));

    auto endpoint_is_vertex = [&](std::size_t edge,
                                   double ey,
                                   std::size_t ey_tag) -> bool {
        assert(edge < polygon.num_edges() && "[C91 §2.2]: invalid edge");
        const auto& e = polygon.edge(edge);
        SymbolicY chord_y{ey, ey_tag};
        return symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.start_idx))) ||
               symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.end_idx)));
    };

    bool left_is_vertex  = endpoint_is_vertex(c.left_edge, c.y, c.y_tag);
    bool right_is_vertex = endpoint_is_vertex(c.right_edge, c.y, c.y_tag);

    auto try_merge_chain = [&](std::size_t a_first, std::size_t a_last) {
        assert(a_first < arc_sequence_.size() && !arc_sequence_[a_first].dead &&
               a_last < arc_sequence_.size() && !arc_sequence_[a_last].dead &&
               "[C91 §2.4(ii)]: adj_arcs must be valid + live");
        // [C91 §2.2 tex 106]: simulate the merge and count over the
        // post-merge underlying P-edge range — the additive `a+b−shared`
        // would over-count when either input wraps ([C91 §2.4 tex 142]).
        Arc merged_arc;
        merged_arc.first_edge = arc_sequence_[a_first].first_edge;
        merged_arc.first_side = arc_sequence_[a_first].first_side;
        merged_arc.last_edge  = arc_sequence_[a_last].last_edge;
        merged_arc.last_side  = arc_sequence_[a_last].last_side;
        auto [lo, hi] = merged_arc.underlying_edge_range(start_vertex, end_vertex);
        std::size_t merged = polygon.count_nonnull_edges(lo, hi);
        if (merged > max_count)
            max_count = merged;
    };
    auto assert_pair = [&](const Chord::AdjArcs& adj) {
        // [C91 §2.2 tex 94]: non-vertex endpoint ⟹ exactly 2 adj arcs,
        // sharing the junction edge (arcs[0] ends at the chord, arcs[1]
        // starts — ∂C traversal order).
        assert(adj.count == 2 &&
               "[C91 §2.2 tex 94]: non-vertex endpoint needs 2 adj arcs");
        assert(arc_sequence_[adj.arcs[0]].last_edge ==
                   arc_sequence_[adj.arcs[1]].first_edge &&
               "[C91 §2.2 tex 94]: adj arcs share junction edge");
    };

    if (!left_is_vertex)  assert_pair(c.left_adj);
    if (!right_is_vertex) assert_pair(c.right_adj);

    if (!left_is_vertex && !right_is_vertex) {
        // [C91 §2.3 tex 121 + §2.2 tex 94]: if one region is a leaf bounded
        // by this chord and a single arc, that arc is adjacent at BOTH
        // endpoints; gluing at both chains three arcs into one, and the
        // contracted weight must reflect the full chain — pairwise merges
        // would under-count it.
        bool shared_lr = (c.left_adj.arcs[1] == c.right_adj.arcs[0]);
        bool shared_rl = (c.right_adj.arcs[1] == c.left_adj.arcs[0]);
        // An arc has one end per endpoint: it cannot END (slot 0) or
        // START (slot 1) at both endpoints of the chord.
        assert(c.left_adj.arcs[0] != c.right_adj.arcs[0] &&
               c.left_adj.arcs[1] != c.right_adj.arcs[1] &&
               "[C91 §2.2 tex 94]: an arc cannot occupy the same slot at "
               "both chord endpoints");
        // Both sharings at once ⟺ ∂C has only these two subdivision
        // points — impossible; see the matching assert in remove_chord.
        assert(!(shared_lr && shared_rl) &&
               "[C91 §2.2 tex 94]: chord endpoints cannot be the only "
               "subdivision points of ∂C");
        if (shared_lr) {
            try_merge_chain(c.left_adj.arcs[0], c.right_adj.arcs[1]);
            return max_count;
        }
        if (shared_rl) {
            try_merge_chain(c.right_adj.arcs[0], c.left_adj.arcs[1]);
            return max_count;
        }
    }

    if (!left_is_vertex)
        try_merge_chain(c.left_adj.arcs[0], c.left_adj.arcs[1]);
    if (!right_is_vertex)
        try_merge_chain(c.right_adj.arcs[0], c.right_adj.arcs[1]);

    return max_count;
}

bool Submap::is_granular(std::size_t gamma,
                          const Polygon& polygon) const noexcept {
    // [C91 §2.3] (i): all weights ≤ γ.
    if (!is_semigranular(gamma)) return false;

    // [C91 §2.3]: "if (i) holds but the submap has no exit chord, it is still
    // said to be γ-granular."
    if (num_live_chords() == 0) return true;

    // [C91 §2.3] (ii): contracting any edge incident on a < 3-degree node
    // yields weight > γ.
    for (std::size_t ci = 0; ci < chords_.size(); ++ci) {
        const auto& c = chords_[ci];
        if (c.dead) continue;
        std::size_t d0 = nodes_[c.region[0]].degree();
        std::size_t d1 = nodes_[c.region[1]].degree();
        if (d0 >= 3 && d1 >= 3) continue;
        if (simulated_contraction_weight(ci, polygon) <= gamma)
            return false;
    }

    // [C91 §2.3 Lemma 2.3 tex 126,129]: γ-granular conformal V(C) submap
    // has V ≤ 2·⌊8(|C|−1)/(γ+1)⌋ regions.  Proof sketch: paper's E (tree
    // edges incident on a deg-<3 node, where (ii) merged-weight>γ applies)
    // satisfies |E|·(γ+1) ≤ 4·Σ_u W_u ≤ 8(n−1); E is a fixed fraction of
    // all tree edges (tex 129), so V = O(n/γ + 1) ≤ 2·⌊8(n−1)/(γ+1)⌋.
    // Requires conformality (is_granular alone doesn't).
    if (is_conformal()) {
        // We're past the "no chords" return, so live arcs exist and
        // [C91 §2.4(iii) tex 138] requires start_vertex/end_vertex set.
        assert(start_vertex != NONE && end_vertex != NONE &&
               "[C91 §2.4(iii) tex 138]: submap with chords needs start/end_vertex");
        std::size_t n_c = end_vertex - start_vertex + 1;
        std::size_t bound = 2 * (8 * (n_c - 1) / (gamma + 1));
        assert(num_live_nodes() <= bound &&
               "[C91 §2.3 Lemma 2.3]: V ≤ 2·⌊8(|C|−1)/(γ+1)⌋");
    }

    return true;
}

} // namespace chazelle
