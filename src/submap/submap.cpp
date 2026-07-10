// src/submap/submap.cpp

#include "submap.h"
#include <cstdio>
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
    // boundary" — LEFT ascends within the arc, RIGHT descends.  A
    // same-side arc violating the order is a DOUBLE-WRAP arc ([C91 §2.4
    // tex 142]: wraps around both endpoints of C), whose legs are
    // validated by Arc::legs' asserts at every call site.

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

    // [C91 §2.4(iii) tex 138]: "pointers to the arc-structures whose
    // corresponding arcs pass through the endpoints" — the wrap classes
    // identify them directly ([C91 §2.4 tex 142]).  Canonical insertion
    // order puts the start-turn arc last, so a later genuine end-wrap
    // arc never overwrites an earlier start claim, and the closed-arc
    // encoding (end wrap with first == last == C's first edge) claims
    // the start turn unless a later arc does.
    if (arc.wraps_end()) end_arc = idx;
    bool closed_encoding =
        arc.first_side == LEFT && arc.last_side == RIGHT &&
        arc.first_edge == arc.last_edge &&
        start_vertex != NONE && arc.first_edge == start_vertex;
    if (arc.wraps_start() || closed_encoding) start_arc = idx;

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
    ++live_chords_;

    // [C91 §2.4(i)]: update region adjacency; both regions must be live.
    for (std::size_t r : chord.region) {
        assert(r != NONE && r < nodes_.size() && !nodes_[r].dead &&
               "[C91 §2.4(i)]: chord must connect two LIVE regions");
        nodes_[r].incident_chords.push_back(idx);
    }

    tree_decomp_dirty_ = true;
    return idx;
}

// ════════════════════════════════════════════════════════════════
//  chord_regions_below_above ([C91 §2.2 tex 96])
// ════════════════════════════════════════════════════════════════

// [C91 §2.2 tex 96]: does the region of `arc_idx` lie strictly ABOVE
// the chord, judged locally at the chord endpoint (edge_at, side_at)?
// (The polygon vertex adjacent to the chord along the arc's traversal
// decides; SoS gives it a definite side.)
namespace {

bool arc_region_above_chord(const Submap& S, const Polygon& C,
                            const Chord& ch, std::size_t edge_at,
                            Side side_at, std::size_t arc_idx) {
    const Arc& a = S.arc(arc_idx);
    bool first_matches =
        (a.first_edge == edge_at && a.first_side == side_at);
    bool last_matches =
        (a.last_edge == edge_at && a.last_side == side_at);
    bool starts_at_chord;
    if (first_matches && !last_matches) {
        starts_at_chord = true;
    } else if (!first_matches && last_matches) {
        starts_at_chord = false;
    } else {
        assert(first_matches && last_matches &&
               "adj arc must touch the chord endpoint via first or last");
        starts_at_chord = symbolic_y_equal(
            S.arc_start_symbolic_y(arc_idx, C), ch.symbolic_y());
    }
    std::size_t adj_v;
    if (starts_at_chord) {
        adj_v = (a.first_side == LEFT) ? C.edge(a.first_edge).end_idx
                                       : C.edge(a.first_edge).start_idx;
    } else {
        adj_v = (a.last_side == LEFT) ? C.edge(a.last_edge).start_idx
                                      : C.edge(a.last_edge).end_idx;
    }
    const SymbolicY vy = symbolic_y_of(C.vertex(adj_v));
    if (symbolic_y_equal(vy, ch.symbolic_y())) {
        // Degenerate: the deciding arc has zero extent at the chord's
        // own level — it is the POLAR CAP's boundary arc (a region
        // bounded by this chord alone, plus ∂C points at its level).
        // Such a chord wraps at a level meeting C only at its source
        // vertex, which is therefore C's global y-extremum ([C91 §2.1
        // tex 70]: a horizontal circle strictly inside C's y-range
        // crosses the connected curve elsewhere).  The cap lies above
        // ⟺ the chord sits at the global maximum.
        const SymbolicY maxy = symbolic_y_of(C.vertex(C.max_y_vertex()));
        if (symbolic_y_equal(ch.symbolic_y(), maxy)) return true;
        assert(symbolic_y_equal(
                   ch.symbolic_y(),
                   symbolic_y_of(C.vertex(C.min_y_vertex()))) &&
               "[C91 §2.1 tex 70]: a polar-cap chord sits at C's "
               "global y-extremum");
        return false;
    }
    return symbolic_y_greater(vy, ch.symbolic_y());
}

} // namespace

// The two regions flanking a chord vertically (constant along the whole
// chord: a chord of a submap separates exactly its two regions).
void Submap::chord_regions_below_above(std::size_t ci, const Polygon& C,
                                    std::size_t* below,
                                    std::size_t* above) const {
    const Submap& S = *this;
    const Chord& c = S.chord(ci);
    const Chord::AdjArcs* adj = nullptr;
    std::size_t edge_at = NONE;
    Side side_at = LEFT;
    if (c.left_adj.count == 2) {
        adj = &c.left_adj; edge_at = c.left_edge; side_at = c.left_side;
    } else if (c.right_adj.count == 2) {
        adj = &c.right_adj; edge_at = c.right_edge; side_at = c.right_side;
    }
    if (adj) {
        bool a0 = arc_region_above_chord(S, C, c, edge_at, side_at,
                                         adj->arcs[0]);
        bool a1 = arc_region_above_chord(S, C, c, edge_at, side_at,
                                         adj->arcs[1]);
        assert(a0 != a1 &&
               "[C91 §2.2 tex 96]: the two adj arcs at a mid-edge chord "
               "endpoint lie on opposite sides of the chord");
        std::size_t r0 = S.arc(adj->arcs[0]).region_node;
        std::size_t r1 = S.arc(adj->arcs[1]).region_node;
        *above = a0 ? r0 : r1;
        *below = a0 ? r1 : r0;
        return;
    }
    // Both endpoints are polygon vertices: one adj arc per endpoint.
    bool left_above = arc_region_above_chord(S, C, c, c.left_edge,
                                             c.left_side,
                                             c.left_adj.arcs[0]);
    std::size_t r_arc = S.arc(c.left_adj.arcs[0]).region_node;
    std::size_t r_other = (c.region[0] == r_arc) ? c.region[1]
                                                 : c.region[0];
    *above = left_above ? r_arc : r_other;
    *below = left_above ? r_other : r_arc;
}


// ── Point-arc test ([C91 §2.2 tex 96/106]) ──────────────────────

bool Submap::arc_is_point(std::size_t arc_idx,
                          const Polygon& polygon) const {
    assert(arc_idx < arc_sequence_.size() && !arc_sequence_[arc_idx].dead);
    const Arc& a = arc_sequence_[arc_idx];
    // A wrap encoding always spans a turnaround crossing; a same-face
    // pair is a point iff its two endpoint positions coincide (same
    // symbolic y on the same face — [C91 §2 tex 47]).
    return a.edge_count == 0 && a.first_side == a.last_side &&
           a.first_edge == a.last_edge &&
           symbolic_y_equal(arc_start_symbolic_y(arc_idx, polygon),
                            arc_end_symbolic_y(arc_idx, polygon));
}

// ── Junction-mate lookup ([C91 §2.2 tex 96]) ────────────────────

std::size_t Submap::find_junction_arc(const Chord& c,
                                      bool query_left,
                                      std::size_t edge, Side side,
                                      std::size_t vertex_idx,
                                      bool want_after,
                                      std::size_t exclude,
                                      std::size_t exclude2,
                                      const Polygon& polygon) const {
    // [C91 §2.1 tex 72 + §2.2 tex 92]: every vertex of ∂C is incident
    // upon exactly one chord, and augmented chords land only at
    // non-vertex points — so once `c` is removed, its vertex-endpoint
    // junction is chord-free and exactly one arc starts (resp. ends)
    // there.
    assert(vertex_idx != NONE && vertex_idx < polygon.num_vertices() &&
           "[C91 §2.2 tex 96]: junction lookup requires a polygon vertex");
    assert(edge < polygon.num_edges() &&
           (vertex_idx == polygon.edge(edge).start_idx ||
            vertex_idx == polygon.edge(edge).end_idx) &&
           "[C91 §2.2 tex 96]: junction vertex must bound the endpoint edge");

    SymbolicY jy = c.symbolic_y();

    // [C91 §2.2 tex 96] alternation: pin the mate's region (see
    // submap.h).  Null-length chords keep the exclusion-based
    // disambiguation (their two slots share one position).
    std::size_t expected_region = NONE;
    if (!c.is_null_length) {
        const Chord::AdjArcs& other = query_left ? c.right_adj
                                                 : c.left_adj;
        std::size_t other_before;
        if (other.count == 1) {
            other_before = other.arcs[0];
        } else {
            const std::size_t oe = query_left ? c.right_edge
                                              : c.left_edge;
            const Side os = query_left ? c.right_side : c.left_side;
            auto starts_there = [&](std::size_t ai) {
                const Arc& a = arc_sequence_[ai];
                return a.first_edge == oe && a.first_side == os &&
                       symbolic_y_equal(
                           arc_start_symbolic_y(ai, polygon), jy);
            };
            other_before = starts_there(other.arcs[0]) ? other.arcs[1]
                                                       : other.arcs[0];
        }
        const std::size_t rb = arc_sequence_[other_before].region_node;
        assert((rb == c.region[0] || rb == c.region[1]) &&
               "[C91 §2.4(ii)]: slot arcs bound the chord's regions");
        expected_region = want_after
            ? rb
            : (c.region[0] == rb ? c.region[1] : c.region[0]);
    }

    // A mate at vertex v carries one of v's two incident edges,
    // v−1 or v (edge i runs vertex i → i+1): arcs bounded at a vertex
    // may record either adjacent edge, and zero-length arcs sit AT the
    // vertex ([C91 §2.1 tex 72] duplicate pair).  Within those two
    // edges, the symbolic y of v occurs ONLY at v itself ([C91 §2 tex
    // 47] SoS: a horizontal line meets an edge at most once, and v is
    // an endpoint of both), so edge ∈ {v−1, v} + side + symbolic-y
    // pins the ∂C point exactly.  A zero-length mate takes precedence —
    // the junction point itself precedes (resp. follows) any nonzero
    // mate along ∂C.
    std::size_t found_zero = NONE;
    std::size_t found_nonzero = NONE;

    auto consider = [&](std::size_t ai) {
        if (ai == NONE || ai == exclude || ai == exclude2) return;
        assert(ai < arc_sequence_.size());
        const Arc& a = arc_sequence_[ai];
        if (a.dead) return;
        if (a.region_node != c.region[0] && a.region_node != c.region[1])
            return;
        if (expected_region != NONE && a.region_node != expected_region)
            return;
        Side a_side = want_after ? a.first_side : a.last_side;
        std::size_t a_edge = want_after ? a.first_edge : a.last_edge;
        if (a_edge + 1 < vertex_idx || a_edge > vertex_idx) return;
        if (a_side != side) {
            // [C91 §2.1 tex 72]: at an interior-EXTREMUM junction the
            // tour crosses the vertex between corner-SIBLING walls on
            // the two adjacent edges, and the sibling's side label
            // flips whenever the edges' orientations differ (inside
            // companions face each other across the wedge; outside
            // companions round the outer turn).  Accept the mate under
            // the sibling label: other edge, same companion class.  At
            // a non-extremum vertex the geometric wall keeps its side
            // label across the corner, so a flipped label there names
            // the curve's OTHER side — never the junction's mate.
            if (!polygon.is_y_extremum(vertex_idx)) return;
            const std::size_t other_e =
                (edge == vertex_idx) ? vertex_idx - 1 : vertex_idx;
            if (a_edge != other_e) return;
            if (is_inside_companion(polygon, a_edge, a_side,
                                    vertex_idx) !=
                is_inside_companion(polygon, edge, side, vertex_idx))
                return;
        }
        SymbolicY ay = want_after ? arc_start_symbolic_y(ai, polygon)
                                  : arc_end_symbolic_y(ai, polygon);
        if (!symbolic_y_equal(ay, jy)) return;
        bool zero = (a.edge_count == 0);
        std::size_t& slot = zero ? found_zero : found_nonzero;
        assert((slot == NONE || slot == ai) &&
               "[C91 §2.2 tex 96]: at most one mate per junction "
               "(distinct arcs cannot share a ∂C start/end point)");
        slot = ai;
    };

    // Candidates: adjacency slots of both regions' incident chords
    // ([C91 §2.4(ii)]).  Exhaustive: under full junction gluing every
    // live arc ends at a live chord endpoint of its region and is that
    // chord's recorded before-arc ([C91 §2.2 tex 96] alternation), and
    // wrap-spanning arcs are single structures with chord-bounded ends
    // ([C91 §2.4 tex 142]) — no arc hides at a chord-free wrap.
    for (std::size_t r : c.region) {
        for (std::size_t ci : nodes_[r].incident_chords) {
            const Chord& ch = chords_[ci];
            if (ch.dead) continue;
            for (std::size_t k = 0; k < ch.left_adj.count; ++k)
                consider(ch.left_adj.arcs[k]);
            for (std::size_t k = 0; k < ch.right_adj.count; ++k)
                consider(ch.right_adj.arcs[k]);
        }
    }

    std::size_t result = (found_zero != NONE) ? found_zero : found_nonzero;
    assert(result != NONE &&
           "[C91 §2.2 tex 96]: an interior junction always has both a "
           "before- and an after-arc");
    return result;
}

// ── Chord removal (O(1) via tombstones) ─────────────────────────

// [C91 §2.2 tex 94]: "remove the chord and those endpoints that are not
// vertices of C, gluing back ∂C at those points."
// [C91 §2.2 tex 96]: region boundaries ALTERNATE exit chords with arcs,
// so the two arcs meeting at EVERY endpoint of the removed chord fuse
// into one (tex 108: a removed chord "ceases to separate any arcs") —
// including at C's endpoint companions, where the merged arc double-backs
// around the endpoint as ONE arc-structure ([C91 §2.4 tex 142]).
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

    // ── [C91 §2.2 tex 94/96]: removing the LAST chord closes ∂C into a
    // single arc.  The chord's ≤ 2 endpoint junctions glue its regions'
    // arcs into a chain whose two ends meet each other, so the result is
    // the chordless submap's single closed arc covering all of ∂C — one
    // arc-structure, stored cut at the canonical traversal origin (C's
    // start turnaround, [C91 §2.4(iii) tex 138] / [C91 §2.4 tex 142]).
    if (num_live_chords() == 1) {
        assert(start_vertex != NONE && end_vertex != NONE &&
               end_vertex > start_vertex &&
               "[C91 §2.4(iii) tex 138]: C endpoints must be identified");
        // The ≤ 2 live arcs (one per region, [C91 §2.2 tex 96]
        // alternation with a single chord) all appear in the chord's own
        // adjacency slots.
        std::size_t keep = NONE;
        auto sweep_slots = [&](const Chord::AdjArcs& adj) {
            for (std::size_t k = 0; k < adj.count; ++k) {
                std::size_t ai = adj.arcs[k];
                assert(ai != NONE && ai < arc_sequence_.size() &&
                       "[C91 §2.4(ii)]: adj_arc must be valid");
                if (keep == NONE) {
                    assert(!arc_sequence_[ai].dead &&
                           "[C91 §2.4(ii)]: adj_arc must be live");
                    keep = ai;
                } else if (ai != keep && !arc_sequence_[ai].dead) {
                    arc_sequence_[ai].dead = true;
                    compacted_ = false;
                }
            }
        };
        sweep_slots(c.left_adj);
        sweep_slots(c.right_adj);
        assert(keep != NONE &&
               "[C91 §2.4(ii)]: a chord references its regions' arcs");

        Arc& a = arc_sequence_[keep];
        a.first_edge = start_vertex;  a.first_side = LEFT;
        a.last_edge  = start_vertex;  a.last_side  = RIGHT;
        // [C91 §2.2 tex 106]: the closed arc is ALL of ∂C — both full
        // sides of every C edge (arc.h::arc_boundary_edge_count; no
        // corner corrections apply: both endpoints sit at the start
        // turnaround, covering their legs fully).
        a.edge_count =
            2 * polygon.count_nonnull_edges(start_vertex, end_vertex - 1);
        a.region_node = r0;
        start_arc = keep;
        end_arc = keep;

        c.dead = true;
        --live_chords_;
        nodes_[r1].dead = true;
        nodes_[r0].incident_chords.clear();
        tree_decomp_dirty_ = true;
        return r0;
    }

    // The polygon vertex at an endpoint, or NONE if the endpoint strands
    // mid-edge.  (Chord's symbolic y matches one of the edge's endpoints.)
    auto endpoint_vertex = [&](std::size_t edge,
                               double ey,
                               std::size_t ey_tag) -> std::size_t {
        assert(edge < polygon.num_edges() && "[C91 §2.2]: invalid edge index");
        const auto& e = polygon.edge(edge);
        SymbolicY chord_y{ey, ey_tag};
        if (symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.start_idx))))
            return e.start_idx;
        if (symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.end_idx))))
            return e.end_idx;
        return NONE;
    };

    std::size_t left_vertex  = endpoint_vertex(c.left_edge, c.y, c.y_tag);
    std::size_t right_vertex = endpoint_vertex(c.right_edge, c.y, c.y_tag);
    bool left_is_vertex  = left_vertex != NONE;
    bool right_is_vertex = right_vertex != NONE;

    // ── Endpoint-y derivation cache ──────────────────────────────
    // A glue extends its survivor past the removed chord's endpoint,
    // but the chord's adjacency slots still name the survivor as the
    // arc bounded there ([C91 §2.4(ii)]) until the removal completes —
    // so a LATER glue in this same removal must read the survivor's
    // true endpoints from the values the earlier glue established, not
    // from the (now stale) slot derivations.  ≤ 2 glues per removal
    // ([C91 §2.2 tex 94]).
    struct GluedYs { std::size_t arc; SymbolicY start, end; };
    GluedYs glued[2];
    std::size_t n_glued = 0;
    auto start_of = [&](std::size_t ai) -> SymbolicY {
        for (std::size_t k = 0; k < n_glued; ++k)
            if (glued[k].arc == ai) return glued[k].start;
        return arc_start_symbolic_y(ai, polygon);
    };
    auto end_of = [&](std::size_t ai) -> SymbolicY {
        for (std::size_t k = 0; k < n_glued; ++k)
            if (glued[k].arc == ai) return glued[k].end;
        return arc_end_symbolic_y(ai, polygon);
    };
    auto record_glued = [&](std::size_t ai, SymbolicY sy, SymbolicY ey) {
        for (std::size_t k = 0; k < n_glued; ++k)
            if (glued[k].arc == ai) {
                glued[k].start = sy;
                glued[k].end = ey;
                return;
            }
        assert(n_glued < 2 && "[C91 §2.2 tex 94]: at most 2 glues");
        glued[n_glued++] = {ai, sy, ey};
    };

    // [C91 §2.2 tex 94]: "glue back ∂C" — merge arc pairs at non-vertex
    // endpoints.  O(1): at most 2 pairs.
    auto glue_arcs = [&](std::size_t ai, std::size_t aj) {
        assert(ai != NONE && ai < arc_sequence_.size() &&
               aj != NONE && aj < arc_sequence_.size() &&
               !arc_sequence_[ai].dead && !arc_sequence_[aj].dead &&
               "[C91 §2.4(ii)]: adj_arc must be valid + live");

        // [C91 §2.2 tex 94 + §2.2 tex 102]: the two arcs glued at an
        // endpoint are always distinct here.  ai == aj would require
        // this chord's two endpoints to be the ONLY subdivision points
        // of ∂C — i.e. this is the submap's last chord, which is
        // handled by the closed-arc path before any gluing runs.
        assert(ai != aj &&
               "[C91 §2.2 tex 94]: glued arc pair must be distinct "
               "(the last-chord closure path handles ai == aj)");

        auto& a_keep = arc_sequence_[ai];
        auto& a_dead = arc_sequence_[aj];

        // [C91 §2.2 tex 96]: glue mates meet at one ∂C point, so they
        // share the junction's ∂C side, and their edges either coincide
        // (mid-edge split, or a zero-length mate carrying the vertex's
        // other incident edge) or are consecutive in traversal direction
        // (vertex junction between edge v−1 and edge v).
        assert(a_keep.last_side == a_dead.first_side &&
               "[C91 §2.2 tex 96]: glue mates share the junction's ∂C side");
        assert((a_keep.last_edge == a_dead.first_edge ||
                (a_keep.last_side == LEFT &&
                 a_keep.last_edge + 1 == a_dead.first_edge) ||
                (a_keep.last_side == RIGHT &&
                 a_keep.last_edge == a_dead.first_edge + 1)) &&
               "[C91 §2.2 tex 96]: glue mates' edges must coincide or be "
               "traversal-consecutive at the junction");

        // [C91 §2.2 tex 106]: edge_count = nonnull ∂C edges of the arc
        // (per leg — both sides of a doubled-over C edge count; see
        // arc.h::arc_boundary_edge_count).
        //
        // POINT pieces ([C91 §2.2 tex 96]: a single ∂C point —
        // arc_is_point, which distinguishes them from zero-WEIGHT arcs
        // spanning null-edge faces of [C91 §4 tex 316]'s padding)
        // contribute NO extent, and their recorded edge index is
        // boundary-ambiguous (either edge incident on the vertex,
        // [C91 §2.1 tex 72]), so folding it into the merged range could
        // pick up an edge the arc does not span.  Every non-point mate —
        // including a zero-weight null-run span and the zero-length WRAP
        // passage around an endpoint of C ([C91 §2.4 tex 142]: the
        // crossing must be donated or the merged structure loses its
        // double-backing) — folds its span into the encoding.
        // Point tests through the cache (a prior survivor's slot
        // derivations are stale, see above).
        auto is_point_here = [&](std::size_t idx) {
            const Arc& a = arc_sequence_[idx];
            return a.edge_count == 0 && a.first_side == a.last_side &&
                   a.first_edge == a.last_edge &&
                   symbolic_y_equal(start_of(idx), end_of(idx));
        };
        const bool dead_is_point = is_point_here(aj);
        const bool keep_is_point = is_point_here(ai);
        SymbolicY merged_start_y{}, merged_end_y{};
        if (dead_is_point) {
            // Merged arc = a_keep plus a point at its end: unchanged.
            merged_start_y = start_of(ai);
            merged_end_y   = end_of(ai);
        } else if (keep_is_point) {
            // Merged arc = a point plus a_dead's span: adopt it.
            merged_start_y = start_of(aj);
            merged_end_y   = end_of(aj);
            a_keep.first_edge = a_dead.first_edge;
            a_keep.first_side = a_dead.first_side;
            a_keep.last_edge  = a_dead.last_edge;
            a_keep.last_side  = a_dead.last_side;
            a_keep.edge_count = a_dead.edge_count;
        } else {
            // Recompute from the merged encoding rather than adding
            // `a+b−shared`: a mid-edge junction's shared edge would be
            // double-counted, and per-leg counting needs the corner
            // corrections.  Endpoint ys are read BEFORE the mutation,
            // while both arcs' bounding chords still reference them
            // (or from the cache for a prior glue's survivor).
            merged_start_y = start_of(ai);
            merged_end_y   = end_of(aj);
            a_keep.last_edge = a_dead.last_edge;
            a_keep.last_side = a_dead.last_side;
            a_keep.edge_count = arc_boundary_edge_count(
                a_keep, polygon, start_vertex, end_vertex,
                merged_start_y, merged_end_y);
        }
        record_glued(ai, merged_start_y, merged_end_y);

        // Tombstone the dead arc.
        a_dead.dead = true;
        compacted_ = false;  // arc-sequence has dead entries now

        // [C91 §2.4(iii)]: If the dead arc was a C-endpoint arc pointer,
        // redirect to the surviving arc (which now covers the merged
        // range including the endpoint's turnaround).
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

    // ── [C91 §2.2 tex 96/108]: glue the boundary back into alternating
    // form at each endpoint.
    //
    // Use the stored null-length flag — recomputing from edge/side would
    // misclassify exit chords whose left_edge == right_edge but sides
    // differ.
    if (c.is_null_length) {
        assert(c.left_adj.count == 1 && c.right_adj.count == 1 &&
               "[C91 §2.2]: null-length chord must have 1 adj arc per side");
        // [C91 §2.1 tex 72]: null-length chord endpoints are y-extremum
        // duplicate vertices — both polygon vertices by construction,
        // never C's endpoints (endpoint duplicates carry no null chord).
        assert(left_is_vertex && right_is_vertex &&
               "[C91 §2.1 tex 72]: null-length chord endpoints are polygon vertices");
        assert(left_vertex == right_vertex &&
               "[C91 §2.1 tex 72]: null-length chord endpoints share the vertex");
        assert(left_vertex != start_vertex && left_vertex != end_vertex &&
               "[C91 §2.1 tex 72]: null-length chords arise only at "
               "non-endpoint local extrema");

        // [C91 §2.2 tex 108]: "once removed, a chord of zero length
        // ceases to separate any arcs" — the outer before-arc, the inner
        // null arc, and the outer after-arc fuse into one arc containing
        // the null edge (tex 96's "three-edge arc (not a four-edge
        // arc!)").  One slot holds the inner null arc (region[1]); the
        // other an outer arc; the third piece is found by scanning.
        std::size_t sl = c.left_adj.arcs[0];
        std::size_t sr = c.right_adj.arcs[0];
        bool sl_inner = arc_sequence_[sl].region_node == c.region[1];
        bool sr_inner = arc_sequence_[sr].region_node == c.region[1];
        assert(sl_inner != sr_inner &&
               "[C91 §2.1 tex 72]: exactly one slot holds the inner null arc");
        std::size_t inner = sl_inner ? sl : sr;
        std::size_t outer = sl_inner ? sr : sl;
        assert(arc_sequence_[inner].edge_count == 0 &&
               "[C91 §2.1 tex 72]: the inner region of a null-length chord "
               "is bounded by a null arc");

        // Classify the outer slot arc: does it END at the junction
        // (before-arc) or START there (after-arc)?  A nonzero arc does
        // exactly one; a zero-length outer arc (both) is glued the same
        // way in either order.
        const Arc& oa = arc_sequence_[outer];
        bool outer_is_before =
            (oa.last_side == c.left_side) &&
            (oa.last_edge == c.left_edge ||
             oa.last_edge + 1 == c.left_edge || oa.last_edge == c.left_edge + 1) &&
            symbolic_y_equal(arc_end_symbolic_y(outer, polygon), c.symbolic_y());
        if (outer_is_before) {
            glue_arcs(outer, inner);
            std::size_t after = find_junction_arc(
                c, /*query_left=*/true, c.left_edge, c.left_side,
                left_vertex, /*want_after=*/true, outer, inner, polygon);
            glue_arcs(outer, after);
        } else {
            std::size_t before = find_junction_arc(
                c, /*query_left=*/true, c.left_edge, c.left_side,
                left_vertex, /*want_after=*/false, outer, inner, polygon);
            glue_arcs(before, inner);
            glue_arcs(before, outer);
        }
    } else {
        // [C91 §2.2 tex 94/96]: mid-edge endpoint ⟹ both halves recorded
        // (2 adj arcs) and the stranded point vanishes; vertex endpoint ⟹
        // the before-arc is recorded and the after-arc is found by
        // scanning ([C91 §2.4(ii)] 1-slot convention).  Both glue —
        // including at C's endpoint companions, where the merged arc
        // double-backs around the endpoint of C as one arc-structure
        // ([C91 §2.4 tex 142]).
        auto glue_endpoint = [&](const Chord::AdjArcs& adj, bool ql,
                                 std::size_t edge, Side side,
                                 std::size_t vtx) {
            if (vtx == NONE) {
                assert(adj.count == 2 &&
                       "[C91 §2.2 tex 94]: non-vertex endpoint needs 2 adj arcs");
                glue_arcs(adj.arcs[0], adj.arcs[1]);
                return;
            }
            assert(adj.count == 1 &&
                   "[C91 §2.2 tex 94]: vertex endpoint records one adj arc");
            std::size_t before = adj.arcs[0];
            std::size_t after = find_junction_arc(
                c, ql, edge, side, vtx, /*want_after=*/true, before,
                NONE, polygon);
            glue_arcs(before, after);
        };
        glue_endpoint(c.left_adj, true, c.left_edge, c.left_side,
                      left_vertex);
        glue_endpoint(c.right_adj, false, c.right_edge, c.right_side,
                      right_vertex);
    }

    // Reassign r1's arcs to r0 by walking r1's incident chords' adj arcs
    // ([C91 §2.4(ii) tex 137]).  All entries are live: glue_arcs's
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

    // ([C91 §2.2 tex 96]: under full junction gluing every live arc of
    // r1 ends at a live chord endpoint of r1 and is that chord's
    // recorded before-arc — the slot walk above is exhaustive, wrap
    // arcs included, [C91 §2.4 tex 142].)

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
    --live_chords_;
    nodes_[r1].dead = true;
    tree_decomp_dirty_ = true;
    return r0;
}

// ── Tombstone compaction ────────────────────────────────────────

// [C91 §3.3]: "We can now put S in normal form."  compact() strips
// dead entries and rebuilds index mappings in O(m).

void Submap::compact() {
    // [C91 §2.2 tex 96]: under full junction gluing, remove_chord's
    // slot walk over the dying region's incident chords reassigns
    // EVERY live arc of that region (each live arc is the recorded
    // before-arc of the live chord at its clockwise end, and cascaded
    // removals move inherited chords into the dying region's incident
    // list), so no live arc can point to a dead region here.
    for (const auto& a : arc_sequence_) {
        if (a.dead) continue;
        assert(a.region_node < nodes_.size() &&
               !nodes_[a.region_node].dead &&
               "[C91 §2.2 tex 96]: live arc must point to a live region "
               "(remove_chord's slot walk is exhaustive)");
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

// ── [C91 §3.3 tex 276]: Normal form ─────────────────────────────

// "We can now put S in normal form, which includes computing its tree
// decomposition.  As we discussed earlier, this can be done very simply
// in time O((n₁/γ₁ + n₂/γ₂ + 1)log(n₁ + n₂))."

void Submap::normalize(const Polygon& polygon) {
    // [C91 §2.4(iv) tex 139]: normal form includes the tree
    // decomposition, which is defined for conformal submaps ([C91 §2.3
    // tex 114]) — [C91 §3.2] has run by the time §3.3 normalizes.
    assert(is_conformal() &&
           "[C91 §2.4(iv)]: normal form requires a conformal submap");

    compact();

    const std::size_t M = arc_sequence_.size();
    if (M > 0) {
        // ── Re-sort the arc-sequence table into canonical ∂C traversal
        // order ([C91 §2.4(iii) tex 138]) — [C91 §3.2]'s insert_chord
        // appends split halves at the table's end, breaking it.
        //
        // Key: clockwise traversal position of the arc's START —
        // LEFT ascending first_edge then RIGHT descending; within one
        // (side, edge), start-y follows the traversal direction (the
        // edge's geometric y-direction on LEFT, reversed on RIGHT); a
        // zero-length arc precedes the nonzero arc sharing its start
        // point (it occupies exactly that point, [C91 §2.2 tex 96]).
        // O(M log M) comparison sort — within tex 276's normal-form
        // budget (log M ≤ log(n₁ + n₂)).
        const std::size_t n_edges = polygon.num_edges();
        auto edge_ascends = [&](std::size_t e) -> bool {
            const auto& ed = polygon.edge(e);
            return symbolic_y_less(symbolic_y_of(polygon.vertex(ed.start_idx)),
                                   symbolic_y_of(polygon.vertex(ed.end_idx)));
        };

        struct Key {
            std::size_t trav;   // cw ∂C position of the start edge
            SymbolicY   y;      // start position's symbolic y
            bool        asc;    // traversal direction along the edge
            bool        zero;   // zero-length arc
            std::size_t idx;    // original table index
        };
        std::vector<Key> keys(M);
        for (std::size_t i = 0; i < M; ++i) {
            const Arc& a = arc_sequence_[i];
            assert(!a.dead && "compact() leaves no dead arcs");
            std::size_t trav = (a.first_side == LEFT)
                ? a.first_edge
                : 2 * n_edges - 1 - a.first_edge;
            bool asc = (a.first_side == LEFT) == edge_ascends(a.first_edge);
            keys[i] = Key{trav, arc_start_symbolic_y(i, polygon), asc,
                          a.edge_count == 0, i};
        }
        std::sort(keys.begin(), keys.end(),
            [](const Key& a, const Key& b) {
                if (a.trav != b.trav) return a.trav < b.trav;
                if (!symbolic_y_equal(a.y, b.y))
                    return a.asc ? symbolic_y_less(a.y, b.y)
                                 : symbolic_y_greater(a.y, b.y);
                return a.zero && !b.zero;
            });
        // [C91 §2.4(iii)]: the traversal position of every arc start is
        // unique up to the zero/nonzero pair sharing a point — two
        // distinct nonzero arcs cannot start at the same ∂C point.
        for (std::size_t i = 0; i + 1 < M; ++i) {
            assert((keys[i].trav != keys[i + 1].trav ||
                    !symbolic_y_equal(keys[i].y, keys[i + 1].y) ||
                    (keys[i].zero && !keys[i + 1].zero)) &&
                   "[C91 §2.4(iii)]: distinct arcs cannot share a ∂C "
                   "start point and length class");
        }

        // Apply the permutation and remap all arc references.
        std::vector<std::size_t> arc_map(M, NONE);
        {
            std::vector<Arc> sorted;
            sorted.reserve(M);
            for (std::size_t p = 0; p < M; ++p) {
                arc_map[keys[p].idx] = p;
                sorted.push_back(arc_sequence_[keys[p].idx]);
            }
            arc_sequence_ = std::move(sorted);
        }
        auto remap_adj = [&](Chord::AdjArcs& adj) {
            for (std::size_t k = 0; k < adj.count; ++k) {
                assert(adj.arcs[k] != NONE && arc_map[adj.arcs[k]] != NONE);
                adj.arcs[k] = arc_map[adj.arcs[k]];
            }
        };
        for (auto& ch : chords_) {
            remap_adj(ch.left_adj);
            remap_adj(ch.right_adj);
        }
        if (start_arc != NONE) start_arc = arc_map[start_arc];
        if (end_arc != NONE)   end_arc = arc_map[end_arc];

        // [C91 §2.4 tex 144]: recompute the LEFT/RIGHT boundary.
        left_right_boundary_ = 0;
        for (std::size_t i = 0; i < M; ++i)
            if (arc_sequence_[i].first_side == LEFT)
                left_right_boundary_ = i + 1;
    }

    // [C91 §2.4(iii) tex 138 + tex 144]: full normal-form validation —
    // canonical order, endpoint pointers, start-y monotonic runs,
    // edge_count caches.
    check_invariants(polygon);

    // [C91 §2.4(iv) tex 139]: "If the submap is conformal, then its tree
    // decomposition should be available."  [C91 §2.3 tex 116]: O(M log M).
    build_tree_decomposition();
}

// ── Live counts ─────────────────────────────────────────────────

std::size_t Submap::num_live_nodes() const noexcept {
    std::size_t n = 0;
    for (const auto& nd : nodes_) if (!nd.dead) ++n;
    return n;
}

std::size_t Submap::num_live_chords() const noexcept {
#ifndef NDEBUG
    std::size_t n = 0;
    for (const auto& ch : chords_) if (!ch.dead) ++n;
    assert(n == live_chords_ && "live chord counter must match the table");
#endif
    return live_chords_;
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

    // [C91 §2.4(iii) tex 138 + §2.4 tex 142]: endpoint arc pointers.
    // C's two turnaround points are never chord endpoints ([C91 §2.1
    // tex 72] case 3: chords end at the companion vertices flanking the
    // turn), so exactly one live arc passes through each turnaround and
    // it double-backs there — a wrap-spanning arc is ONE arc-structure,
    // never split.
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
        assert(start_vertex != NONE &&
               "[C91 §2.4(iii)]: start_vertex required when arcs exist");
        assert(end_vertex != NONE && end_vertex > start_vertex &&
               "[C91 §2.4(iii)]: end_vertex required when arcs exist");

        std::size_t last_live = NONE;
        std::size_t last_live_left = NONE;
        std::size_t wrap_count = 0;
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            const Arc& a = arc_sequence_[i];
            if (a.dead) continue;
            last_live = i;
            if (a.first_side == LEFT) last_live_left = i;
            if (a.wraps()) ++wrap_count;
        }

        // [C91 §2.4(iii) tex 138]: start_arc passes through C's start
        // turnaround — the latest clockwise start position, hence the
        // LAST table entry.
        assert(start_arc == last_live &&
               "[C91 §2.4(iii) tex 138]: start_arc must be the last "
               "table entry (its cw start position is maximal)");

        // [C91 §2.4 tex 142/144]: end_arc passes through C's end
        // turnaround — the last LEFT-starting arc when one exists
        // (first_side LEFT, last position past the turn); otherwise the
        // submap's arcs all start on the RIGHT half and the double-wrap
        // arc at the table's end covers both turns.
        if (last_live_left != NONE) {
            assert(end_arc == last_live_left &&
                   "[C91 §2.4 tex 144]: end_arc must be the last "
                   "LEFT-starting arc");
        } else {
            assert(end_arc == last_live &&
                   "[C91 §2.4 tex 142]: with no LEFT-starting arc the "
                   "double-wrap arc covers both turnarounds");
        }

        const Arc& sa = arc_sequence_[start_arc];
        const Arc& ea = arc_sequence_[end_arc];

        // Wrap classes ([C91 §2.4 tex 142]): end_arc double-backs at
        // C's end; start_arc at C's start.  The closed arc of a
        // chordless submap uses the end-wrap encoding cut at the start
        // turnaround (first == last == C's first edge).
        bool closed = (num_live_chords() == 0);
        if (closed) {
            assert(num_live_arcs() == 1 && start_arc == end_arc &&
                   "[C91 §2.2 tex 96]: a chordless submap is one region "
                   "bounded by the single closed arc (all of ∂C)");
            assert(sa.first_side == LEFT && sa.last_side == RIGHT &&
                   sa.first_edge == start_vertex &&
                   sa.last_edge == start_vertex &&
                   "[C91 §2.4 tex 142/138]: the closed arc is stored cut "
                   "at C's start turnaround");
        } else {
            assert(ea.wraps_end() &&
                   "[C91 §2.4 tex 142]: end_arc must double-back around "
                   "C's end vertex (one arc-structure, never split)");
            assert(sa.wraps_start() &&
                   "[C91 §2.4 tex 142]: start_arc must double-back around "
                   "C's start vertex (one arc-structure, never split)");
        }

        // [C91 §2.4 tex 142]: the wrap-class arcs are exactly the ≤ 2
        // endpoint arcs (they coincide for a double-wrap/closed arc).
        std::size_t expected_wraps = (start_arc == end_arc) ? 1 : 2;
        assert(wrap_count == expected_wraps &&
               "[C91 §2.4 tex 142]: only the arcs passing through C's "
               "endpoints double-back");
        for (std::size_t i = 0; i < arc_sequence_.size(); ++i) {
            const Arc& a = arc_sequence_[i];
            if (a.dead || !a.wraps()) continue;
            assert((i == start_arc || i == end_arc) &&
                   "[C91 §2.4 tex 142]: a wrapping arc must be an "
                   "endpoint arc");
        }
    }
}

void Submap::check_invariants(const Polygon& polygon) const {
    check_invariants();

    // [C91 §2.4 tex 144]: within a same-(first_side, first_edge) run,
    // start-y must be monotonic IN THE EDGE'S GEOMETRIC TRAVERSAL
    // DIRECTION — the edge's own y-direction on the LEFT half, reversed
    // on the RIGHT half ([C91 §2.4(iii) tex 138]: canonical clockwise
    // order).  double_identify Phase 2 binary-searches with exactly
    // this direction (keys_ascending), so the check must use it too: a
    // run monotone only in the opposite direction would pass a weaker
    // "some direction" check yet silently break Phase 2's search.
    // Ties (leq/geq) occur at a null-length chord's duplicate-vertex
    // pair and are compatible with both directions.
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
            // with dead gaps inside).  Direction = the edge's geometric
            // y-direction on the LEFT half, reversed on the RIGHT half
            // (double_identify Phase 2's keys_ascending).
            if (j > i) {
                const auto& e = polygon.edge(arc_sequence_[i].first_edge);
                bool edge_ascending = symbolic_y_less(
                    symbolic_y_of(polygon.vertex(e.start_idx)),
                    symbolic_y_of(polygon.vertex(e.end_idx)));
                bool asc = (arc_sequence_[i].first_side == LEFT)
                    ? edge_ascending : !edge_ascending;
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

    // [C91 §2.2 tex 106]: arc.edge_count caches the arc's nonnull ∂C
    // edge count (per leg, arc.h::arc_boundary_edge_count) used by
    // region_weight and simulated_contraction_weight.  Validate the
    // cache — a wrong value silently miscomputes weights, breaking
    // γ-granularity decisions.
    //
    // [C91 §2.4 tex 133]: null-length arcs encode a single ∂C point (the
    // null-length chord's duplicate-vertex pair), not an edge span.
    // Under [C91 §2.1] every
    // polygon edge has nonzero length, so edge_count == 0 uniquely
    // identifies null-length arcs — skip the cache check for them.
    //
    // Wrapped arcs ([C91 §2.4 tex 142]): underlying_edge_range (arc.h)
    // covers all legs as one contiguous C range — exactly the arcs
    // passing through C's endpoints (start_arc / end_arc) wrap.
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
            std::size_t actual = arc_boundary_edge_count(
                a, polygon, start_vertex, end_vertex,
                arc_start_symbolic_y(i, polygon),
                arc_end_symbolic_y(i, polygon));
            assert(a.edge_count == actual &&
                   "[C91 §2.2 tex 106]: arc.edge_count cache must match "
                   "the arc's per-leg nonnull ∂C edge count "
                   "(arc_boundary_edge_count)");
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

        // [C91 §2.1 tex 70]: every chord is horizontal at its source
        // vertex's y, with y_tag = that vertex's SoS index — null-length
        // chords included ([C91 §2.1 tex 72]: the inside pair duplicates
        // the extremum vertex itself, so its symbolic y is the vertex's
        // own (y, SoS index); no other perturbation source exists).
        // [C91 §2.4 tex 133]: C is a contiguous subchain of P, so its
        // vertex indices are consecutive from vertex(0).index (asserted
        // by the Polygon constructor) — translate tag → table position.
        const std::size_t first_tag = polygon.vertex(0).index;
        assert(c.y_tag >= first_tag &&
               c.y_tag - first_tag < polygon.num_vertices() &&
               "[C91 §2.1 tex 70]: chord y_tag must be a vertex of C");
        SymbolicY chord_y{c.y, c.y_tag};
        assert(symbolic_y_equal(chord_y,
                   symbolic_y_of(polygon.vertex(c.y_tag - first_tag))) &&
               "[C91 §2.1 tex 70]: chord must be horizontal at its source vertex");

        if (c.is_null_length) {
            // [C91 §2.1 tex 72]: a null-length chord arises only from the
            // inside duplicate pair of a NON-ENDPOINT local y-extremum
            // (endpoint duplicates, case 3, carry no null chord), and its
            // coincident endpoints sit AT that vertex — on one of its two
            // incident edges.  add_chord asserts the structural half
            // (coincident endpoints, one adj arc per ∂C side).
            const std::size_t v = c.y_tag - first_tag;
            assert(polygon.is_y_extremum(v) &&
                   "[C91 §2.1 tex 72]: null-length chord source must be an "
                   "interior local y-extremum of C");
            assert((c.left_edge == v || (v >= 1 && c.left_edge == v - 1)) &&
                   "[C91 §2.1 tex 72]: null-length chord endpoints sit at "
                   "the extremum vertex (edge v-1 or v)");
            continue;
        }

        assert((c.left_adj.count == 1) == matches_an_endpoint(c.left_edge, chord_y) &&
               "[C91 §2.2 tex 94]: LEFT endpoint count == 1 ⟺ endpoint is a polygon vertex");
        assert((c.right_adj.count == 1) == matches_an_endpoint(c.right_edge, chord_y) &&
               "[C91 §2.2 tex 94]: RIGHT endpoint count == 1 ⟺ endpoint is a polygon vertex");

        // [C91 §2.2 tex 94]: removing the chord doesn't glue ∂C at a vertex
        // endpoint — meaning the single adj arc on that side already spans
        // the polygon vertex c.y_tag (no split there).  Verify the vertex
        // lies in the adj arc's edge range.
        auto arc_spans_vertex = [&](std::size_t arc_idx) -> bool {
            // [C91 §2.4 tex 142]: the underlying range ᾱ covers every C
            // edge the arc touches on either side; vertex v lies on
            // edge v−1 or v.
            const auto& a = arc_sequence_[arc_idx];
            auto [elo, ehi] =
                a.underlying_edge_range(start_vertex, end_vertex);
            // [C91 §2.4 tex 133]: the range is in C-local positions;
            // translate the global SoS tag before comparing.
            const std::size_t v = c.y_tag - first_tag;
            return v >= elo && v <= ehi + 1;
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
    //   LEFT  half: [0, lrb)   — LEFT-starting arcs, ascending first_edge
    //   RIGHT half: [lrb, end) — RIGHT-starting arcs, descending
    // The ≤ 2 arcs passing through C's endpoints double-back across the
    // seams ([C91 §2.4 tex 142]): the start-turn arc (last table entry)
    // re-enters the LEFT sequence's front, and the end-turn arc (end_arc)
    // re-enters the RIGHT sequence's front.
    assert(start_arc == arc_sequence_.size() - 1 &&
           "[C91 §2.4(iii) tex 138]: start_arc must be the last table entry");
    assert((left_right_boundary_ == 0 ||
            end_arc == left_right_boundary_ - 1) &&
           "[C91 §2.4 tex 144]: end_arc must be the last LEFT-starting arc");

    std::size_t left_begin  = 0;
    std::size_t left_end    = left_right_boundary_;
    std::size_t right_begin = left_right_boundary_;
    std::size_t right_end   = arc_sequence_.size();

    // [C91 §2.4 tex 144]: two-phase binary search.
    //   Phase 1: bsearch by first_edge → contiguous run of arcs on edge_idx.
    //   Phase 2: bsearch by start-y within that run → arc(s) at y.
    //   Plus O(1) check for a boundary arc covering edge_idx from before
    //   the run: the run's predecessor within the half, or — when the
    //   run has no predecessor — the half's SEAM arc, whose double-backed
    //   leg(s) cover the sequence's front ([C91 §2.4 tex 142/144]).
    // Total: O(log m).
    auto search_half = [&](std::size_t lo, std::size_t hi, bool ascending,
                           std::size_t seam) {
        const Side half_side = ascending ? LEFT : RIGHT;

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

        // [C91 §2.4 tex 142]: side-aware edge coverage via the arc's legs.
        auto arc_contains_edge = [&](std::size_t ai) -> bool {
            assert(ai < arc_sequence_.size() && "[C91 §2.4]: invalid arc index");
            const auto& a = arc_sequence_[ai];
            assert(!a.dead && "[C91 §2.4]: arc_contains_edge on dead arc");
            assert(!a.wraps() || ai == start_arc || ai == end_arc);
            return a.covers(edge_idx, half_side, start_vertex, end_vertex);
        };

        // The pre-run portion of edge_idx belongs to exactly one arc
        // (arcs tile ∂C): the run's in-half predecessor, or the seam arc
        // when the run sits at the sequence's front — never both (an
        // in-half predecessor's start would lie inside the seam arc's
        // leg otherwise).
        std::size_t boundary_arc = NONE;
        if (blo > lo) {
            if (arc_contains_edge(blo - 1))
                boundary_arc = blo - 1;
        } else if (seam != NONE && !arc_sequence_[seam].dead &&
                   !(seam >= blo && seam < bend) &&
                   arc_contains_edge(seam)) {
            boundary_arc = seam;
        }

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

        // [C91 §2.4 tex 144] Phase 2: bsearch by start-y.  [C91
        // §2.4(iii) tex 138]: within an edge, start-y follows the
        // clockwise traversal direction — the edge's own y-direction on
        // the LEFT half, reversed on the RIGHT half.  (Inferring the
        // direction from the run's endpoint keys instead is wrong when
        // those keys tie, e.g. when a zero-length arc at the run's
        // front shares its start-y with its successor.)
        bool keys_ascending;
        {
            assert(edge_idx < polygon.num_edges());
            const auto& e = polygon.edge(edge_idx);
            bool edge_ascending = symbolic_y_less(
                symbolic_y_of(polygon.vertex(e.start_idx)),
                symbolic_y_of(polygon.vertex(e.end_idx)));
            keys_ascending = ascending ? edge_ascending : !edge_ascending;
        }

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
        } else {
            // ylo == blo: y symbolically precedes the first interval
            // arc's start-y.  [C91 §2.4 tex 144] queries lie ON the
            // edge and arcs tile ∂C, so the piece of this edge-side
            // before blo's start belongs to an arc that started earlier
            // — the boundary arc, which must therefore exist (blo
            // starting mid-edge with no predecessor would leave ∂C
            // uncovered).  Pushing blo here would report an arc that
            // does not pass through q.
            assert(boundary_arc != NONE &&
                   "[C91 §2.4 tex 144]: query precedes the run's first "
                   "start-y with no boundary arc — off-edge query or "
                   "corrupt arc-sequence table");
        }

        // Boundary arc (if present) always covers edge_idx; its y-range
        // runs from its start-y to the first interval arc's start-y.
        if (boundary_arc != NONE) {
            const Arc& ba = arc_sequence_[boundary_arc];
            // "Query is exactly the boundary arc's start point" only
            // makes sense when that start lies on THIS edge and side —
            // a seam arc's start sits on another leg entirely
            // ([C91 §2.4 tex 142]), where an equal symbolic y names a
            // different ∂C point.
            bool starts_here = ba.first_edge == edge_idx &&
                               ba.first_side == half_side;
            if (starts_here &&
                symbolic_y_equal(arc_start_symbolic_y(boundary_arc,
                                                      polygon), y)) {
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
    // chord + arc at a y-extremum); total ≤ 6 across both halves.  Seam
    // arcs: the start-turn arc (last table entry) doubles back into the
    // LEFT sequence's front; the end-turn arc (end_arc) into the RIGHT's.
    search_half(left_begin, left_end, /*ascending=*/true,
                /*seam=*/arc_sequence_.size() - 1);
    assert(result.count <= 3 && "[C91 §2.4 tex 144]: ≤ 3 arcs per ∂C half");

    std::size_t left_count = result.count;
    search_half(right_begin, right_end, /*ascending=*/false,
                /*seam=*/end_arc);
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

    // [C91 §2.1 tex 72] case 3: a chord at a companion of C's OWN
    // endpoint.  The count-1 vertex-slot convention records only the
    // before-arc there, and the after-arc passes THROUGH the turnaround
    // (rebuild's sweep deliberately keeps the turnaround's single
    // incident edge as its first pointer, [C91 §2.4 tex 142]), so its
    // start is the C-endpoint vertex itself — NOT the face's
    // traversal-start vertex derived below.  Exactly two positions
    // behave this way: the start vertex's RIGHT companion (the tour's
    // last position — its after-arc is the start_arc) and the end
    // vertex's LEFT companion (just before the end turnaround — its
    // after-arc is the end_arc).  A chord there pins that arc's start
    // to the companion exactly: no ∂C position lies between the
    // companion and the turnaround.
    if (start_vertex != NONE && end_vertex != NONE) {
        auto companion_chord_at = [&](std::size_t edge, Side side,
                                      SymbolicY vy) -> bool {
            for (std::size_t ci : nodes_[a.region_node].incident_chords) {
                const Chord& c = chords_[ci];
                if (c.dead) continue;
                if (!symbolic_y_equal(c.symbolic_y(), vy)) continue;
                if ((c.left_edge == edge && c.left_side == side) ||
                    (c.right_edge == edge && c.right_side == side))
                    return true;
            }
            return false;
        };
        if (arc_idx == start_arc && a.first_side == RIGHT &&
            a.first_edge == start_vertex) {
            SymbolicY vy = symbolic_y_of(polygon.vertex(start_vertex));
            if (companion_chord_at(a.first_edge, RIGHT, vy)) return vy;
        }
        if (arc_idx == end_arc && a.first_side == LEFT &&
            a.first_edge + 1 == end_vertex) {
            SymbolicY vy = symbolic_y_of(polygon.vertex(end_vertex));
            if (companion_chord_at(a.first_edge, LEFT, vy)) return vy;
        }
    }

    // [C91 §2.2 tex 96]: an arc may START at a chord endpoint that
    // sits AT a polygon vertex — the count-1 vertex-slot convention
    // records only the BEFORE-arc there ([C91 §2.4(ii) tex 137]), so
    // the junction's after-arc reaches this derivation unrecorded.  It
    // starts at the face's traversal-END vertex (where the chord
    // attaches), not the traversal-start vertex of the input-table
    // derivation below.
    // The chord's endpoint label and the arc's first label both name
    // that corner point, possibly under the two adjacent edges' labels
    // ([C91 §2.1 tex 72] corner siblings; at an extremum the sibling's
    // side flips with edge orientation).
    {
        const std::size_t vend = (a.first_side == LEFT)
                                     ? a.first_edge + 1
                                     : a.first_edge;
        const SymbolicY vy = symbolic_y_of(polygon.vertex(vend));
        const std::size_t other_e = (a.first_edge == vend)
                                        ? vend - 1
                                        : vend;
        auto corner_label_matches = [&](std::size_t ce, Side cs) {
            if (ce == a.first_edge && cs == a.first_side) return true;
            if (ce != other_e || other_e >= polygon.num_edges())
                return false;
            if (polygon.is_y_extremum(vend))
                return is_inside_companion(polygon, ce, cs, vend) ==
                       is_inside_companion(polygon, a.first_edge,
                                           a.first_side, vend);
            return cs == a.first_side;
        };
        for (std::size_t ci : nodes_[a.region_node].incident_chords) {
            const Chord& c = chords_[ci];
            if (c.dead || c.is_null_length) continue;
            if (!symbolic_y_equal(c.symbolic_y(), vy)) continue;
            const Chord::AdjArcs* adj = nullptr;
            if (corner_label_matches(c.left_edge, c.left_side))
                adj = &c.left_adj;
            else if (corner_label_matches(c.right_edge, c.right_side))
                adj = &c.right_adj;
            if (!adj || adj->count != 1) continue;
            if (adj->arcs[0] == arc_idx) continue;   // the before-arc
            return c.symbolic_y();
        }
    }

    // [C91 §2.1 tex 72]: remaining arc-starts not at chords sit at
    // polygon vertices — the outside-pair companions of an interior
    // y-extremum (tex 72 case 2), the outer after-arc of a null-length
    // chord's duplicate pair, or the closed arc's canonical cut at C's
    // start turnaround ([C91 §2.4 tex 142/138]).  The polygon's input
    // table ([C91 §2.4(iii) tex 138]) carries the y directly: edge
    // first_edge starts at vertex(first_edge) under LEFT traversal
    // (ascending) and at vertex(first_edge+1) under RIGHT traversal
    // (descending).
    assert(a.first_edge < polygon.num_edges() &&
           "[C91 §2.4(iii)]: first_edge must be a valid edge index");
    std::size_t vidx = (a.first_side == LEFT)
        ? a.first_edge
        : a.first_edge + 1;
    assert(vidx < polygon.num_vertices() &&
           "[C91 §2.4(iii)]: arc-start polygon vertex must be valid");
    return symbolic_y_of(polygon.vertex(vidx));
}

SymbolicY Submap::arc_end_symbolic_y(std::size_t arc_idx,
                                     const Polygon& polygon) const {
    assert(arc_idx < arc_sequence_.size() &&
           !arc_sequence_[arc_idx].dead &&
           "[C91 §2.4]: arc index must be valid + live");
    const Arc& a = arc_sequence_[arc_idx];

    // [C91 §2.4(ii) tex 137]: walk the arc's region's incident chords.
    // Mirror of arc_start_symbolic_y — the arc ENDING at a chord endpoint
    // is slot 0 of a count-2 pair, or the single count-1 slot (the
    // "before-arc" of a polygon-vertex endpoint).
    for (std::size_t ci : nodes_[a.region_node].incident_chords) {
        const Chord& c = chords_[ci];
        if (c.dead) continue;
        if (c.is_null_length) {
            // [C91 §2.1 tex 72]: the inner null-arc both starts AND ends
            // at the null-length chord's position.
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
        // [C91 §2.4 tex 137 + tex 133]: slot 0 is the arc ending at the
        // chord; a count-1 (polygon-vertex) slot holds the before-arc,
        // which likewise ends there ([C91 §2.2 tex 94]).
        if (c.left_adj.arcs[0] == arc_idx)
            return c.symbolic_y();
        if (c.right_adj.arcs[0] == arc_idx)
            return c.symbolic_y();
    }

    // [C91 §2.1 tex 72 / §2.4(iii) tex 138]: no bounding chord at the end —
    // the arc ends at a polygon vertex (an extremum companion, or the
    // closed arc's canonical cut at C's start turnaround,
    // [C91 §2.4 tex 142]); the input table carries the y directly.  Under LEFT
    // traversal (ascending) edge last_edge ends at vertex(last_edge + 1);
    // under RIGHT traversal (descending) at vertex(last_edge).
    assert(a.last_edge < polygon.num_edges() &&
           "[C91 §2.4(iii)]: last_edge must be a valid edge index");
    std::size_t vidx = (a.last_side == LEFT)
        ? a.last_edge + 1
        : a.last_edge;
    assert(vidx < polygon.num_vertices() &&
           "[C91 §2.4(iii)]: arc-end polygon vertex must be valid");
    return symbolic_y_of(polygon.vertex(vidx));
}

// ── [C91 §2.2 tex 106]: edge_count refresh ──────────────────────

void Submap::refresh_arc_edge_counts(const Polygon& polygon) {
    assert(start_vertex != NONE && end_vertex != NONE &&
           "[C91 §2.4(iii)]: C endpoints must be identified");
    for (std::size_t ai = 0; ai < arc_sequence_.size(); ++ai) {
        Arc& a = arc_sequence_[ai];
        if (a.dead || a.edge_count == 0) continue;   // null arcs stay 0
        a.edge_count = arc_boundary_edge_count(
            a, polygon, start_vertex, end_vertex,
            arc_start_symbolic_y(ai, polygon),
            arc_end_symbolic_y(ai, polygon));
    }
}

// ── [C91 §3.2]: Chord insertion ─────────────────────────────────

// [C91 §3.2 tex 264]: "For any region with more than four arcs, let us
// apply Lemma 3.2 to every pair of nonconsecutive arcs until we find a
// chord which we can add to S."  This is that chord addition: split the
// region at the two mutually visible points p and q.

Submap::InsertChordResult Submap::insert_chord(
        const ChordPointSpec& p, const ChordPointSpec& q,
        SymbolicY y, std::size_t region,
        const std::size_t* cycle, std::size_t cycle_len,
        const Polygon& polygon) {
    assert(region < nodes_.size() && !nodes_[region].dead &&
           "[C91 §3.2]: insert_chord requires a live region");
    assert(y.tag != SOS_NONE &&
           "[C91 §2 tex 47 (SoS)]: chord must carry its source SoS tag");
    assert(p.arc != q.arc &&
           "[C91 §3.2 Lemma 3.3]: the chord connects two distinct "
           "(nonconsecutive) arcs of the region");
    assert((p.edge != q.edge || p.side != q.side) &&
           "[C91 §2.1 tex 70]: a chord connects two distinct ∂C points. "
           "Equal x IS legal: the outside duplicate pair at a y-extremum "
           "([C91 §2.1 tex 72]) sees itself through infinity ([C91 §2.1 "
           "tex 70]), giving a zero-geometric-length wrap chord between "
           "distinct ∂C points (chord_runs_through_infinity)");
    assert(cycle_len >= 2 && "[C91 §3.2]: region cycle must hold both arcs");

    // Validate the caller-supplied cycle and locate the two split arcs.
    std::size_t ip = NONE, iq = NONE;
    for (std::size_t i = 0; i < cycle_len; ++i) {
        std::size_t ai = cycle[i];
        assert(ai < arc_sequence_.size() && !arc_sequence_[ai].dead &&
               arc_sequence_[ai].region_node == region &&
               "[C91 §2.2 tex 96]: cycle entries must be live arcs of the region");
        if (ai == p.arc) ip = i;
        if (ai == q.arc) iq = i;
    }
    assert(ip != NONE && iq != NONE &&
           "[C91 §3.2]: p.arc and q.arc must appear in the region cycle");

    // [C91 §2.1 tex 70]: a ∂C point's horizontal visibility is unique, so
    // a point coinciding with an existing chord endpoint cannot source a
    // NEW chord — its visibility is already realized by that chord.
    // Position identity: (edge, side, symbolic y) — a horizontal line at
    // a fixed symbolic y crosses an edge at most once.
    auto assert_not_chord_endpoint = [&](std::size_t edge, Side side) {
        for (std::size_t ci : nodes_[region].incident_chords) {
            const Chord& c = chords_[ci];
            if (c.dead) continue;
            if (!symbolic_y_equal(c.symbolic_y(), y)) continue;
            assert(!((c.left_edge == edge && c.left_side == side) ||
                     (c.right_edge == edge && c.right_side == side)) &&
                   "[C91 §2.1 tex 70]: new chord endpoint coincides with an "
                   "existing chord endpoint (visibility already realized)");
        }
    };
    assert_not_chord_endpoint(p.edge, p.side);
    assert_not_chord_endpoint(q.edge, q.side);

    // ── Split the two arcs.  The first half keeps the table slot; the
    // second half is appended (canonical order restored by [C91 §3.3
    // tex 276] "We can now put S in normal form").
    auto split_arc = [&](const ChordPointSpec& sp) -> std::size_t {
        // Read the derived endpoint ys BEFORE mutating anything.
        SymbolicY start_y = arc_start_symbolic_y(sp.arc, polygon);
        SymbolicY end_y   = arc_end_symbolic_y(sp.arc, polygon);

        Arc before = arc_sequence_[sp.arc];
        assert(before.covers(sp.edge, sp.side, start_vertex, end_vertex) &&
               "[C91 §3.2]: split point must lie on the arc "
               "([C91 §2.4 tex 142]: wrap arcs cover per-leg)");

        bool at_start = (sp.edge == before.first_edge &&
                         sp.side == before.first_side &&
                         symbolic_y_equal(y, start_y));
        bool at_end   = (sp.edge == before.last_edge &&
                         sp.side == before.last_side &&
                         symbolic_y_equal(y, end_y));
        // Both at once would make the arc zero-length; zero-length arcs
        // ([C91 §2.2 tex 96]) carry no visibility candidates and are never
        // split ([C91 §2.1 tex 72]: their point sees only its companion).
        assert(!(at_start && at_end) &&
               "[C91 §3.2]: cannot split a zero-length arc");

        // [C91 §2.4 tex 142]: cutting an arc at a ∂C point yields two
        // arcs whose (first, last) encodings — and hence wrap classes —
        // emerge directly: the half containing a turnaround keeps its
        // double-backing.
        Arc after = before;
        after.first_edge = sp.edge;
        after.first_side = sp.side;
        before.last_edge = sp.edge;
        before.last_side = sp.side;

        // [C91 §2.4 tex 133]: the arc-structure's edge pointers are the
        // edges the arc actually spans.  A MID-EDGE split point divides
        // sp.edge between the halves (both record it), but a POLYGON-
        // VERTEX split point puts sp.edge entirely on one side:
        //   at sp.edge's traversal-END vertex → the edge belongs to the
        //     before-half; the after-half starts with the next edge.
        //   at sp.edge's traversal-START vertex → the edge belongs to
        //     the after-half; the before-half ends with the previous
        //     edge.
        // (LEFT traverses edge e from vertex e to e+1; RIGHT reversed.)
        // Without this, the vertex half's range over-counts an edge and
        // its derived start/end position (the [C91 §2.4(iii) tex 138]
        // input-table vertex derivation) reads the wrong vertex.
        //
        // Exception: when the vertex is an ENDPOINT OF C (the split
        // point is a companion vertex, [C91 §2.1 tex 72] case 3), the
        // adjacent half continues THROUGH the turnaround onto the other
        // ∂C side ([C91 §2.4 tex 142]) — its pointer stays on sp.edge
        // (the turnaround's single incident edge) and the half's wrap
        // class encodes the crossing.
        if (!at_start && !at_end) {
            const auto& pe = polygon.edge(sp.edge);
            SymbolicY y_lo = symbolic_y_of(polygon.vertex(pe.start_idx));
            SymbolicY y_hi = symbolic_y_of(polygon.vertex(pe.end_idx));
            bool at_trav_end = (sp.side == LEFT)
                ? symbolic_y_equal(y, y_hi)
                : symbolic_y_equal(y, y_lo);
            bool at_trav_start = (sp.side == LEFT)
                ? symbolic_y_equal(y, y_lo)
                : symbolic_y_equal(y, y_hi);
            assert(!(at_trav_start && at_trav_end) &&
                   "[C91 §2 tex 47 (SoS)]: an edge's endpoints have "
                   "distinct symbolic ys");
            if (at_trav_end) {
                std::size_t trav_end_v = (sp.side == LEFT)
                    ? pe.end_idx : pe.start_idx;
                if (trav_end_v != start_vertex && trav_end_v != end_vertex) {
                    after.first_edge = (sp.side == LEFT) ? sp.edge + 1
                                                         : sp.edge - 1;
                    assert(after.first_edge < polygon.num_edges() &&
                           "[C91 §3.2]: vertex split's after-half must have "
                           "a following edge (its span is nonempty)");
                }
            } else if (at_trav_start) {
                std::size_t trav_start_v = (sp.side == LEFT)
                    ? pe.start_idx : pe.end_idx;
                if (trav_start_v != start_vertex &&
                    trav_start_v != end_vertex) {
                    before.last_edge = (sp.side == LEFT) ? sp.edge - 1
                                                         : sp.edge + 1;
                    assert(before.last_edge < polygon.num_edges() &&
                           "[C91 §3.2]: vertex split's before-half must have "
                           "a preceding edge (its span is nonempty)");
                }
            }
        }

        // [C91 §2.2 tex 106]: edge_count = nonnull ∂C edges per leg
        // (arc.h::arc_boundary_edge_count) — the split point's y and
        // the original endpoint ys (read above, before mutation) give
        // the halves' exact extents; zero-length halves get 0
        // ([C91 §2.4 tex 133] null arcs).
        if (at_start) {
            before.edge_count = 0;
        } else {
            before.edge_count = arc_boundary_edge_count(
                before, polygon, start_vertex, end_vertex, start_y, y);
        }
        if (at_end) {
            after.edge_count = 0;
        } else {
            after.edge_count = arc_boundary_edge_count(
                after, polygon, start_vertex, end_vertex, y, end_y);
        }

        arc_sequence_[sp.arc] = before;
        std::size_t after_idx = arc_sequence_.size();
        arc_sequence_.push_back(after);
        return after_idx;
    };
    std::size_t p_after = split_arc(p);
    std::size_t q_after = split_arc(q);

    // ── New region: the boundary chain running clockwise (along ∂C)
    // from p to q — p's after-half, the cycle arcs strictly between,
    // and q's before-half ([C91 §2.2 tex 96] Lemma 2.2: ∂C order = region
    // boundary order).  Everything else stays in `region`.
    std::size_t r_new = add_node();
    arc_sequence_[p_after].region_node = r_new;
    for (std::size_t i = (ip + 1) % cycle_len; ; i = (i + 1) % cycle_len) {
        assert(i != ip && "[C91 §3.2]: chain walk must reach q before p");
        arc_sequence_[cycle[i]].region_node = r_new;
        if (i == iq) break;    // q's before-half (kept slot) included
    }

    // ── Re-point adjacency slots of the region's chords that referenced
    // a split arc.  [C91 §2.4(ii) tex 137] slot convention: slot 1 of a
    // count-2 pair is the arc STARTING at the chord (start unchanged —
    // the before-half keeps the slot); every other slot references the
    // arc ENDING there — now the appended after-half.
    auto repoint = [&](Chord::AdjArcs& adj, std::size_t old_arc,
                       std::size_t after_idx) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            if (adj.arcs[k] != old_arc) continue;
            if (!(adj.count == 2 && k == 1))
                adj.arcs[k] = after_idx;
        }
    };
    for (std::size_t ci : nodes_[region].incident_chords) {
        Chord& c = chords_[ci];
        assert(!c.dead && "[C91 §2.4]: incident_chords holds live chords only");
        repoint(c.left_adj,  p.arc, p_after);
        repoint(c.left_adj,  q.arc, q_after);
        repoint(c.right_adj, p.arc, p_after);
        repoint(c.right_adj, q.arc, q_after);
    }

    // ── Move chords whose boundary footprint went to the new region.
    // The new chord pq crosses no existing chord (chords are horizontal;
    // distinct symbolic ys under SoS [C91 §2 tex 47]), so each chord's
    // region-side arcs land in ONE chain — read the side off any slot.
    {
        std::vector<std::size_t> stay;
        stay.reserve(nodes_[region].incident_chords.size());
        for (std::size_t ci : nodes_[region].incident_chords) {
            Chord& c = chords_[ci];
            std::size_t side_region = NONE;
            auto scan = [&](const Chord::AdjArcs& adj) {
                for (std::size_t k = 0; k < adj.count; ++k) {
                    std::size_t rn = arc_sequence_[adj.arcs[k]].region_node;
                    if (rn != region && rn != r_new) continue;
                    assert((side_region == NONE || side_region == rn) &&
                           "[C91 §3.2]: a chord's boundary footprint must "
                           "lie entirely in one chain (pq crosses no chord)");
                    side_region = rn;
                }
            };
            scan(c.left_adj);
            scan(c.right_adj);
            assert(side_region != NONE &&
                   "[C91 §2.4(ii)]: chord must reference an arc of its region");
            if (side_region == r_new) {
                if (c.region[0] == region) c.region[0] = r_new;
                else {
                    assert(c.region[1] == region);
                    c.region[1] = r_new;
                }
                nodes_[r_new].incident_chords.push_back(ci);
            } else {
                stay.push_back(ci);
            }
        }
        nodes_[region].incident_chords = std::move(stay);
    }

    // ── The new chord itself.  [C91 §2.2 tex 94]: a polygon-vertex
    // endpoint records ONE adj arc (the before-arc); a mid-edge endpoint
    // records both halves.
    auto endpoint_is_vertex = [&](std::size_t edge) -> bool {
        assert(edge < polygon.num_edges());
        const auto& e = polygon.edge(edge);
        return symbolic_y_equal(y, symbolic_y_of(polygon.vertex(e.start_idx))) ||
               symbolic_y_equal(y, symbolic_y_of(polygon.vertex(e.end_idx)));
    };
    auto make_adj = [&](const ChordPointSpec& sp,
                        std::size_t after_idx) -> Chord::AdjArcs {
        Chord::AdjArcs adj;
        adj.arcs[0] = sp.arc;          // before-half (ends at the chord)
        if (endpoint_is_vertex(sp.edge)) {
            adj.count = 1;
        } else {
            adj.arcs[1] = after_idx;   // after-half (starts at the chord)
            adj.count = 2;
        }
        return adj;
    };

    Chord nc;
    nc.region[0] = region;
    nc.region[1] = r_new;
    nc.y = y.y;
    nc.y_tag = y.tag;
    nc.is_null_length = false;
    // Slot order along the chord: ascending x ([C91 §2.4(ii)]).  A
    // zero-geometric-length wrap chord (equal x, [C91 §2.1 tex 70/72])
    // leaves the x-order vacuous; ties are refined by clockwise ∂C
    // position — the same canonicalization as the fusion rebuild's
    // inventory dedup (fusion.cpp).
    auto tie_pos = [&](std::size_t edge, Side side) {
        return (side == LEFT) ? edge
                              : (2 * polygon.num_edges() - 1 - edge);
    };
    const bool p_left =
        p.x < q.x ||
        (p.x == q.x && tie_pos(p.edge, p.side) < tie_pos(q.edge, q.side));
    const ChordPointSpec& lp = p_left ? p : q;
    const ChordPointSpec& rp = p_left ? q : p;
    std::size_t l_after = p_left ? p_after : q_after;
    std::size_t r_after = p_left ? q_after : p_after;
    nc.left_edge  = lp.edge;  nc.left_side  = lp.side;
    nc.right_edge = rp.edge;  nc.right_side = rp.side;
    nc.left_adj  = make_adj(lp, l_after);
    nc.right_adj = make_adj(rp, r_after);
    std::size_t chord_idx = add_chord(nc);

    // ── [C91 §2.4(iii) tex 138]: keep the C-endpoint arc pointers on
    // the arcs passing through the turnarounds — splitting a wrap arc
    // leaves the double-backing in exactly one half ([C91 §2.4 tex 142]),
    // read off its wrap class.
    auto repoint_wrap = [&](const ChordPointSpec& sp,
                            std::size_t after_idx) {
        if (end_arc == sp.arc &&
            !arc_sequence_[sp.arc].wraps_end() &&
            arc_sequence_[after_idx].wraps_end())
            end_arc = after_idx;
        if (start_arc == sp.arc &&
            !arc_sequence_[sp.arc].wraps_start() &&
            arc_sequence_[after_idx].wraps_start())
            start_arc = after_idx;
    };
    repoint_wrap(p, p_after);
    repoint_wrap(q, q_after);

    // The appended after-halves break the canonical arc-sequence order
    // ([C91 §2.4(iii) tex 138]); [C91 §3.3 tex 276] re-normalizes.  Clear
    // `compacted_` so double_identify fails fast rather than silently
    // binary-searching an unordered table.
    compacted_ = false;
    tree_decomp_dirty_ = true;

    return InsertChordResult{chord_idx, r_new, p_after, q_after};
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

    // [C91 §2.2 tex 96]: under full junction gluing every live arc ends
    // at a live chord endpoint of its region and is that chord's recorded
    // before-arc — the slot walk above is exhaustive ([C91 §2.4 tex 142]:
    // wrap-spanning arcs are single chord-bounded structures).  The one
    // exception is the chordless submap's single closed arc, reachable
    // through the C-endpoint pointers ([C91 §2.4(iii) tex 138]).
    if (nd.incident_chords.empty()) {
        assert(num_live_chords() == 0 &&
               "[C91 §2.2 tex 102]: a chord-free region exists only in "
               "the chordless (single-region) submap");
        assert(start_arc != NONE && start_arc == end_arc &&
               start_arc < arc_sequence_.size() &&
               !arc_sequence_[start_arc].dead &&
               "[C91 §2.4(iii) tex 138]: the chordless submap's closed "
               "arc is the endpoint arc");
        const auto& a = arc_sequence_[start_arc];
        if (a.region_node == node_idx && a.edge_count > max_count)
            max_count = a.edge_count;
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
    return simulated_contraction_weight(chord_idx, polygon,
                                        region_weight(c.region[0]),
                                        region_weight(c.region[1]));
}

std::size_t Submap::simulated_contraction_weight(
        std::size_t chord_idx,
        const Polygon& polygon,
        std::size_t w0, std::size_t w1) const noexcept {
    assert(chord_idx < chords_.size() && !chords_[chord_idx].dead);
    const auto& c = chords_[chord_idx];
    assert(c.region[0] != NONE && c.region[1] != NONE &&
           c.region[0] != c.region[1] &&
           "[C91 §2.4(i)/§2.2 tex 102]: chord regions valid + distinct (tree)");

    // [C91 §2.3 tex 121/123]: merged region's weight = max over both
    // regions' arcs, plus the arcs glued at the removed chord's
    // endpoints — simulated with exactly remove_chord's semantics
    // ([C91 §2.2 tex 96/108]: glue at every endpoint, C-endpoint
    // companions included).
    std::size_t max_count = std::max(w0, w1);

    // [C91 §2.2 tex 94 / §2.4 tex 142]: contracting the LAST chord
    // closes ∂C into the single closed arc spanning all of C — mirror
    // of remove_chord's closure path.
    if (num_live_chords() == 1) {
        assert(start_vertex != NONE && end_vertex != NONE &&
               end_vertex > start_vertex &&
               "[C91 §2.4(iii) tex 138]: C endpoints must be identified");
        // [C91 §2.2 tex 106]: the closed arc is ALL of ∂C — both full
        // sides of every C edge (mirror of remove_chord's closure).
        std::size_t merged =
            2 * polygon.count_nonnull_edges(start_vertex, end_vertex - 1);
        return std::max(max_count, merged);
    }

    auto endpoint_vertex = [&](std::size_t edge) -> std::size_t {
        assert(edge < polygon.num_edges() && "[C91 §2.2]: invalid edge");
        const auto& e = polygon.edge(edge);
        SymbolicY chord_y{c.y, c.y_tag};
        if (symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.start_idx))))
            return e.start_idx;
        if (symbolic_y_equal(chord_y, symbolic_y_of(polygon.vertex(e.end_idx))))
            return e.end_idx;
        return NONE;
    };

    // Glue events remove_chord would perform: (before, after) per
    // glued junction.
    struct GluePair { std::size_t before, after; };
    GluePair pairs[2];
    std::size_t np = 0;

    if (c.is_null_length) {
        // [C91 §2.2 tex 108]: outer-before + inner null + outer-after
        // fuse into one arc.
        std::size_t sl = c.left_adj.arcs[0];
        std::size_t sr = c.right_adj.arcs[0];
        bool sl_inner = arc_sequence_[sl].region_node == c.region[1];
        bool sr_inner = arc_sequence_[sr].region_node == c.region[1];
        assert(sl_inner != sr_inner &&
               "[C91 §2.1 tex 72]: exactly one slot holds the inner null arc");
        std::size_t inner = sl_inner ? sl : sr;
        std::size_t outer = sl_inner ? sr : sl;
        std::size_t v = endpoint_vertex(c.left_edge);
        assert(v != NONE &&
               "[C91 §2.1 tex 72]: null-length chord endpoints are polygon vertices");
        const Arc& oa = arc_sequence_[outer];
        bool outer_is_before =
            (oa.last_side == c.left_side) &&
            (oa.last_edge == c.left_edge ||
             oa.last_edge + 1 == c.left_edge || oa.last_edge == c.left_edge + 1) &&
            symbolic_y_equal(arc_end_symbolic_y(outer, polygon), c.symbolic_y());
        std::size_t missing = find_junction_arc(
            c, /*query_left=*/true, c.left_edge, c.left_side, v,
            /*want_after=*/outer_is_before, outer, inner, polygon);
        std::size_t before = outer_is_before ? outer : missing;
        std::size_t after  = outer_is_before ? missing : outer;
        pairs[np++] = {before, inner};
        pairs[np++] = {inner, after};
    } else {
        auto collect = [&](const Chord::AdjArcs& adj, bool ql,
                           std::size_t edge, Side side) {
            std::size_t v = endpoint_vertex(edge);
            if (v == NONE) {
                assert(adj.count == 2 &&
                       "[C91 §2.2 tex 94]: non-vertex endpoint needs 2 adj arcs");
                pairs[np++] = {adj.arcs[0], adj.arcs[1]};
                return;
            }
            assert(adj.count == 1 &&
                   "[C91 §2.2 tex 94]: vertex endpoint records one adj arc");
            std::size_t after = find_junction_arc(
                c, ql, edge, side, v, /*want_after=*/true, adj.arcs[0],
                NONE, polygon);
            pairs[np++] = {adj.arcs[0], after};
        };
        collect(c.left_adj, true, c.left_edge, c.left_side);
        collect(c.right_adj, false, c.right_edge, c.right_side);
    }

    // [C91 §2.2 tex 106]: the simulated merged arc's count must equal
    // what remove_chord's glue_arcs would produce — fold the chain
    // members left to right with glue_arcs's exact zero-length algebra
    // to obtain the merged ENCODING, then count nonnull ∂C edges per
    // leg (arc.h::arc_boundary_edge_count).  The merged arc's endpoint
    // ys come from the chain's OUTER bounds, whose chords are all
    // still live during simulation.
    auto chain_count = [&](const std::size_t* members, std::size_t n) {
        for (std::size_t i = 0; i < n; ++i)
            assert(members[i] < arc_sequence_.size() &&
                   !arc_sequence_[members[i]].dead &&
                   "[C91 §2.4(ii)]: glue mates must be valid + live");

        // Fold with glue_arcs's exact algebra: POINT members
        // (arc_is_point) add nothing; every other member — including
        // zero-weight null-run spans ([C91 §4 tex 316] padding) and
        // zero-length wrap passages ([C91 §2.4 tex 142]) — folds its
        // span into the encoding.
        Arc acc = arc_sequence_[members[0]];
        bool acc_point = arc_is_point(members[0], polygon);
        for (std::size_t i = 1; i < n; ++i) {
            const Arc& d = arc_sequence_[members[i]];
            if (arc_is_point(members[i], polygon)) {
                // adds a point only
            } else if (acc_point) {
                acc = d;    // point + span: adopt wholesale
                acc_point = false;
            } else {
                acc.last_edge = d.last_edge;
                acc.last_side = d.last_side;
            }
        }
        if (acc_point) return;   // all-point chain: one ∂C point
        std::size_t merged = arc_boundary_edge_count(
            acc, polygon, start_vertex, end_vertex,
            arc_start_symbolic_y(members[0], polygon),
            arc_end_symbolic_y(members[n - 1], polygon));
        if (merged > max_count)
            max_count = merged;
    };

    if (np == 2) {
        // [C91 §2.3 tex 121 + §2.2 tex 94]: an arc adjacent to the chord
        // at BOTH endpoints (leaf region bounded by this chord and one
        // arc, or the null-chord chain) links the two glue events into a
        // single three-arc chain — pairwise counting would under-count.
        bool shared_lr = (pairs[0].after == pairs[1].before);
        bool shared_rl = (pairs[1].after == pairs[0].before);
        // An arc has one start and one end: it cannot end (resp. start)
        // at both endpoints of the chord.
        assert(pairs[0].before != pairs[1].before &&
               pairs[0].after != pairs[1].after &&
               "[C91 §2.2 tex 94]: an arc cannot take the same junction "
               "role at both chord endpoints");
        // Both sharings at once ⟺ this chord's endpoints are the ONLY
        // subdivision points of ∂C — the last-chord case, handled by the
        // closure early-return above.
        assert(!(shared_lr && shared_rl) &&
               "[C91 §2.2 tex 102]: doubly-shared glue pairs occur only "
               "for the last chord (closure early-return)");
        if (shared_lr) {
            std::size_t chain[3] = {pairs[0].before, pairs[0].after,
                                    pairs[1].after};
            chain_count(chain, 3);
            return max_count;
        }
        if (shared_rl) {
            std::size_t chain[3] = {pairs[1].before, pairs[1].after,
                                    pairs[0].after};
            chain_count(chain, 3);
            return max_count;
        }
    }
    for (std::size_t i = 0; i < np; ++i) {
        std::size_t chain[2] = {pairs[i].before, pairs[i].after};
        chain_count(chain, 2);
    }

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
