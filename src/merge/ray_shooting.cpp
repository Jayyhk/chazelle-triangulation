// src/merge/ray_shooting.cpp
//
// [C91 §3.4 tex 282–312]: the ray-shooting structure of Lemma 3.6.
// All geometry is taken over the SoS-perturbed curve ([C91 §2 tex 47]):
// the paper assumes distinct y-coordinates, and the perturbation
// restores full general position, so "zero length" below always means
// symbolic coincidence of positions on ∂C.

#include "ray_shooting.h"
#include "fusion.h"          // shooting_direction ([C91 §2.1 tex 72])

#include <algorithm>
#include <array>

namespace chazelle {

namespace {

// ════════════════════════════════════════════════════════════════
//  Positions on one side of ∂C
// ════════════════════════════════════════════════════════════════
//
// A position is (edge, symbolic y); within an edge, the traversal
// order follows the edge's y-direction on the LEFT side and reverses
// on the RIGHT ([C91 §2.4(iii) tex 138]).  All interval bookkeeping is
// done in LEFT-side coordinates (RIGHT intervals are mirrored), so the
// comparator ascends with the LEFT traversal.

struct Pos {
    std::size_t edge = NONE;
    SymbolicY y{};
};

bool edge_ascends(const Polygon& C, std::size_t e) {
    const auto& ed = C.edge(e);
    return symbolic_y_less(symbolic_y_of(C.vertex(ed.start_idx)),
                           symbolic_y_of(C.vertex(ed.end_idx)));
}

// -1 / 0 / +1 in LEFT-canonical ∂C order.
int pos_order(const Polygon& C, const Pos& a, const Pos& b) {
    if (a.edge != b.edge) return a.edge < b.edge ? -1 : 1;
    if (symbolic_y_equal(a.y, b.y)) return 0;
    bool asc = edge_ascends(C, a.edge);
    bool less = asc ? symbolic_y_less(a.y, b.y)
                    : symbolic_y_greater(a.y, b.y);
    return less ? -1 : 1;
}

// A boundary leg: a maximal one-sided piece of an arc, in LEFT-mirrored
// coordinates (lo precedes hi in LEFT-canonical order).
struct Leg {
    Side side = LEFT;
    std::size_t lo_edge = NONE, hi_edge = NONE;  // inclusive edge range
    Pos lo, hi;
};

// Decompose a live arc into its 1–3 legs ([C91 §2.4 tex 142]: wrapped
// arcs double-back around one — or both — endpoints of C).  Interior
// leg boundaries sit at the turnarounds; the arc's own endpoints carry
// the derived ys ([C91 §2.4 tex 133]).
std::size_t arc_to_legs(const Submap& S, const Polygon& C,
                        std::size_t ai, Leg out[3]) {
    const Arc& a = S.arc(ai);
    assert(S.start_vertex != NONE && S.end_vertex != NONE &&
           "[C91 §2.4(iii)]: C endpoints must be set");
    ArcLeg lg[3];
    std::size_t n = a.legs(S.start_vertex, S.end_vertex, lg);

    SymbolicY sy_start = S.arc_start_symbolic_y(ai, C);
    SymbolicY sy_end = S.arc_end_symbolic_y(ai, C);
    SymbolicY start_wrap_y = symbolic_y_of(C.vertex(S.start_vertex));
    SymbolicY end_wrap_y = symbolic_y_of(C.vertex(S.end_vertex));

    for (std::size_t i = 0; i < n; ++i) {
        Leg l;
        l.side = lg[i].side;
        l.lo_edge = lg[i].lo;
        l.hi_edge = lg[i].hi;
        // Leg boundary positions in LEFT-mirrored order: the position
        // at lo_edge's low end and hi_edge's high end.  The arc's own
        // FIRST position is the traversal start of the first leg (its
        // LEFT-order low end on LEFT, high end on RIGHT); mirrored
        // symmetrically for the last leg; every other leg boundary is a
        // turnaround.
        bool first_leg = (i == 0);
        bool last_leg = (i + 1 == n);
        if (l.side == LEFT) {
            l.lo = {l.lo_edge, first_leg ? sy_start
                                         : start_wrap_y};
            l.hi = {l.hi_edge, last_leg ? sy_end : end_wrap_y};
        } else {
            // RIGHT traversal descends; the leg's LEFT-mirrored low end
            // is its traversal END.
            l.lo = {l.lo_edge, last_leg ? sy_end : start_wrap_y};
            l.hi = {l.hi_edge, first_leg ? sy_start
                                         : end_wrap_y};
        }
        out[i] = l;
    }
    return n;
}

// Does arc ai cover ∂C position (edge, side)?
bool arc_covers(const Submap& S, const Polygon& C, std::size_t ai,
                std::size_t edge, Side side) {
    (void)C;
    assert(S.start_vertex != NONE && S.end_vertex != NONE);
    return S.arc(ai).covers(edge, side, S.start_vertex, S.end_vertex);
}

// ════════════════════════════════════════════════════════════════
//  Ray/edge contacts ([C91 §2 tex 47] symbolic crossings)
// ════════════════════════════════════════════════════════════════

// x of the crossing of the perturbed horizontal line at sy with edge e,
// if any.
bool crossing_x(const Polygon& C, std::size_t e, SymbolicY sy, double* x) {
    const auto& ed = C.edge(e);
    const Point& vs = C.vertex(ed.start_idx);
    const Point& ve = C.vertex(ed.end_idx);
    SymbolicY y0 = symbolic_y_of(vs);
    SymbolicY y1 = symbolic_y_of(ve);
    if (symbolic_y_equal(sy, y0)) { *x = vs.x; return true; }
    if (symbolic_y_equal(sy, y1)) { *x = ve.x; return true; }
    bool between = (symbolic_y_less(y0, sy) && symbolic_y_less(sy, y1)) ||
                   (symbolic_y_less(y1, sy) && symbolic_y_less(sy, y0));
    if (!between) return false;
    double t = (sy.y - vs.y) / (ve.y - vs.y);
    *x = vs.x + t * (ve.x - vs.x);
    return true;
}

// [C91 §2.1 tex 72]: the ∂C side struck by a ray traveling in `dir` —
// traveling RIGHT it arrives from −x and strikes the −x-facing wall.
Side struck_side(const Polygon& C, std::size_t e, Side dir) {
    Side minus_x = edge_ascends(C, e) ? LEFT : RIGHT;
    return (dir == RIGHT) ? minus_x : (minus_x == LEFT ? RIGHT : LEFT);
}

// A candidate contact in the wrap metric of [C91 §2.1 tex 70]:
// direct contacts (d > 0) precede all through-infinity ones (d ≤ 0);
// within a class the ray order is ascending d.
struct Contact {
    bool hit = false;
    double x = 0.0;
    std::size_t edge = NONE;
    Side side = LEFT;
    bool wrapped = false;
    double d = 0.0;

    void offer(double cx, std::size_t ce, Side cs, double cd) {
        bool cw = (cd <= 0.0);
        bool better;
        if (!hit) better = true;
        else if (cw != wrapped) better = !cw;
        else better = cd < d;
        if (better) {
            hit = true;
            x = cx; edge = ce; side = cs; wrapped = cw; d = cd;
        }
    }
};

void scan_edge_range(const Polygon& C, std::size_t lo, std::size_t hi,
                     Point p, SymbolicY sy, Side dir, Contact& best) {
    for (std::size_t e = lo; e <= hi; ++e) {
        double x;
        if (!crossing_x(C, e, sy, &x)) continue;
        double d = (dir == RIGHT) ? (x - p.x) : (p.x - x);
        best.offer(x, e, struck_side(C, e, dir), d);
    }
}

// [C91 §3.4 tex 306]: naive shooting within one region — "checking all
// the O(γ) edges of the region (and not the edges of the face, which
// might be much more numerous)".  Every geometric crossing of an edge
// bounding the region is a contact, whichever ∂C side it strikes.
void scan_region(const Submap& S, const Polygon& C,
                 const std::vector<std::size_t>& arcs,
                 Point p, SymbolicY sy, Side dir, Contact& best) {
    for (std::size_t ai : arcs) {
        Leg legs[3];
        std::size_t n = arc_to_legs(S, C, ai, legs);
        for (std::size_t i = 0; i < n; ++i)
            scan_edge_range(C, legs[i].lo_edge, legs[i].hi_edge,
                            p, sy, dir, best);
    }
}

// ════════════════════════════════════════════════════════════════
//  Chord-side helpers
// ════════════════════════════════════════════════════════════════

// [C91 §2.2 tex 96]: does the region of `arc_idx` lie strictly ABOVE
// the chord, judged locally at the chord endpoint (edge_at, side_at)?
// (The polygon vertex adjacent to the chord along the arc's traversal
// decides; SoS gives it a definite side.)
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
    return symbolic_y_greater(symbolic_y_of(C.vertex(adj_v)),
                              ch.symbolic_y());
}

// The two regions flanking a chord vertically (constant along the whole
// chord: a chord of a submap separates exactly its two regions).
void chord_vertical_regions(const Submap& S, const Polygon& C,
                            std::size_t ci, std::size_t* below,
                            std::size_t* above) {
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

// [s2-adjacency-convention]: the arc STARTING (after) / ENDING (before)
// at a chord endpoint.  Mid-edge endpoints record both; polygon-vertex
// endpoints record only the before-arc, whose ∂C successor (the next
// live table entry — the table tiles ∂C in traversal order,
// [C91 §2.4(iii) tex 138]) starts there.
//
// The starting arc is identified by its FULL start position (edge +
// side + derived start-y, [C91 §2.4 tex 133]) — the y alone is
// ambiguous when the other adj arc begins at the chord's source vertex
// (same symbolic y, different ∂C position).
bool arc_starts_at_slot(const Submap& S, const Polygon& C,
                        const Chord& c, bool left_slot,
                        std::size_t arc_idx) {
    const Arc& a = S.arc(arc_idx);
    std::size_t se = left_slot ? c.left_edge : c.right_edge;
    Side ss = left_slot ? c.left_side : c.right_side;
    return a.first_edge == se && a.first_side == ss &&
           symbolic_y_equal(S.arc_start_symbolic_y(arc_idx, C),
                            c.symbolic_y());
}
std::size_t chord_before_arc(const Submap& S, const Polygon& C,
                             const Chord& c, bool left_slot) {
    const Chord::AdjArcs& adj = left_slot ? c.left_adj : c.right_adj;
    if (adj.count == 1) return adj.arcs[0];
    assert(adj.count == 2);
    bool s0 = arc_starts_at_slot(S, C, c, left_slot, adj.arcs[0]);
    assert(s0 != arc_starts_at_slot(S, C, c, left_slot, adj.arcs[1]) &&
           "[s2-adjacency-convention]: exactly one adj arc starts at a "
           "mid-edge chord endpoint");
    return s0 ? adj.arcs[1] : adj.arcs[0];
}
std::size_t chord_after_arc(const Submap& S, const Polygon& C,
                            const Chord& c, bool left_slot) {
    const Chord::AdjArcs& adj = left_slot ? c.left_adj : c.right_adj;
    if (adj.count == 1)
        return (adj.arcs[0] + 1) % S.num_arcs();
    assert(adj.count == 2);
    bool s0 = arc_starts_at_slot(S, C, c, left_slot, adj.arcs[0]);
    assert(s0 != arc_starts_at_slot(S, C, c, left_slot, adj.arcs[1]) &&
           "[s2-adjacency-convention]: exactly one adj arc starts at a "
           "mid-edge chord endpoint");
    return s0 ? adj.arcs[0] : adj.arcs[1];
}

} // namespace

// ════════════════════════════════════════════════════════════════
//  Construction
// ════════════════════════════════════════════════════════════════

RayShootingStructure::RayShootingStructure(const Submap& S,
                                           const Polygon& C,
                                           std::size_t gamma)
    : S_(&S), C_(&C), gamma_(gamma) {
    // [C91 §3 tex 166]: "each S_i is given in normal form"; [C91 §3.4
    // tex 286]: S is a γ-granular conformal submap of V(C).
#ifdef CHAZELLE_EXPENSIVE_ASSERTS
    // Gated: O(m) full-structure validation inside preprocessing whose
    // budget is O(μ log m) ([C91 §3.4 tex 297]) — same convention as
    // assert_cut_postconditions (oracle.h).
    S.check_invariants(C);
#endif
    assert(S.is_conformal() &&
           "[C91 §3.4 tex 286]: S must be conformal");
    assert(S.is_semigranular(gamma) &&
           "[C91 §3.4 tex 286]: S must be γ-(semi)granular — the naive "
           "region scans rely on O(γ) edges per region (tex 306)");

    build_faces();
    if (mu_ > 1) {
        // [C91 §3.4 tex 297]: "If μ = 1, then ray-shooting can be done
        // trivially in O(m) time, so let us assume that μ > 1."
        build_dual_graph_and_decomposition();
        build_vertical_line();
    }
}

void RayShootingStructure::build_faces() {
    const Submap& S = *S_;
    const Polygon& C = *C_;

    arcs_of_region_.assign(S.num_nodes(), {});
    for (std::size_t ai = 0; ai < S.num_arcs(); ++ai) {
        assert(!S.arc(ai).dead &&
               "[C91 §2.4]: normal form has no dead arcs");
        arcs_of_region_[S.arc(ai).region_node].push_back(ai);
    }

    // [C91 §3.4 tex 286]: "each face of S* coincides exactly with a
    // distinct region of S ... the correspondence face/region is not
    // surjective because empty regions have no associated faces."  In
    // the perturbed model a region is empty iff its whole boundary
    // collapses to one point: every arc has symbolically-zero extent
    // and every incident chord is null-length.
    face_of_region_.assign(S.num_nodes(), NONE);
    region_of_face_.clear();
    for (std::size_t r = 0; r < S.num_nodes(); ++r) {
        if (S.node(r).dead) continue;
        bool nonempty = false;
        for (std::size_t ai : arcs_of_region_[r]) {
            Leg legs[3];
            std::size_t n = arc_to_legs(S, C, ai, legs);
            for (std::size_t i = 0; i < n && !nonempty; ++i)
                if (pos_order(C, legs[i].lo, legs[i].hi) != 0)
                    nonempty = true;
        }
        for (std::size_t ci : S.node(r).incident_chords)
            if (!S.chord(ci).dead && !S.chord(ci).is_null_length)
                nonempty = true;
        if (nonempty) {
            face_of_region_[r] = region_of_face_.size();
            region_of_face_.push_back(r);
        }
        // [C91 §2.3 tex 114]: conformal ⟹ ≤ 4 arcs (degree ≤ 4) — a
        // wrap-spanning arc is ONE arc-structure ([C91 §2.4 tex 142]).
        assert(arcs_of_region_[r].size() <= 4 &&
               "[C91 §2.3 tex 114]: conformal region has ≤ 4 arc-structures");
    }
    mu_ = region_of_face_.size();
    assert(mu_ >= 1 && "V(C) has at least one nonempty region");
}

void RayShootingStructure::build_dual_graph_and_decomposition() {
    const Submap& S = *S_;
    const Polygon& C = *C_;

    // ── Dual-graph edges ([C91 §3.4 tex 295]) ────────────────────
    // "Two faces are adjacent if and only if either they share a chord
    // or one of them has a chord endpoint that abuts on a nonnull
    // length arc of the region associated with the other face. ... It
    // can also be done faster by merging chord endpoints along both
    // sides of C."  We implement the merging variant: each live arc
    // contributes its legs as intervals on the two sides of C; a
    // symbolically-positive overlap between a LEFT and a RIGHT interval
    // is a shared S*-edge.
    struct DualEdge {
        std::size_t fa, fb;      // face ids
        bool dead = false;
    };
    std::vector<DualEdge> edges;

    struct Interval {
        Pos lo, hi;
        std::size_t arc;
        std::size_t leg;         // cycle-ordinal of the leg within the arc
        std::size_t region;
    };
    std::vector<Interval> ivl, ivr;
    for (std::size_t ai = 0; ai < S.num_arcs(); ++ai) {
        Leg legs[3];
        std::size_t n = arc_to_legs(S, C, ai, legs);
        for (std::size_t i = 0; i < n; ++i) {
            if (pos_order(C, legs[i].lo, legs[i].hi) == 0)
                continue;                       // zero length: contracted
            Interval iv{legs[i].lo, legs[i].hi, ai, i,
                        S.arc(ai).region_node};
            (legs[i].side == LEFT ? ivl : ivr).push_back(iv);
        }
    }
    auto by_lo = [&](const Interval& a, const Interval& b) {
        int o = pos_order(C, a.lo, b.lo);
        if (o != 0) return o < 0;
        return pos_order(C, a.hi, b.hi) < 0;
    };
    std::sort(ivl.begin(), ivl.end(), by_lo);
    std::sort(ivr.begin(), ivr.end(), by_lo);

    // Retain the interval lists for the query's double identification
    // ([C91 §3.4 tex 306–308]).
    auto retain = [&](const std::vector<Interval>& src,
                      std::vector<BoundaryInterval>& dst) {
        dst.clear();
        dst.reserve(src.size());
        for (const Interval& iv : src)
            dst.push_back(BoundaryInterval{iv.lo.edge, iv.hi.edge,
                                           iv.lo.y, iv.hi.y, iv.region});
    };
    retain(ivl, left_ivals_);
    retain(ivr, right_ivals_);

    // Overlap sweep; every positive overlap is one abutment edge.
    // Record it per (arc, leg) for the rotation build — a wrap-spanning
    // arc's legs occupy separate stretches of its boundary cycle
    // ([C91 §2.4 tex 142]), so their events must interleave in cycle
    // order.  LEFT lists come out in LEFT-traversal order, RIGHT lists
    // in reverse.
    std::vector<std::array<std::vector<std::size_t>, 3>> leg_events(
        S.num_arcs());
    {
        std::size_t i = 0, j = 0;
        while (i < ivl.size() && j < ivr.size()) {
            const Interval& a = ivl[i];
            const Interval& b = ivr[j];
            const Pos& lo = (pos_order(C, a.lo, b.lo) >= 0) ? a.lo : b.lo;
            const Pos& hi = (pos_order(C, a.hi, b.hi) <= 0) ? a.hi : b.hi;
            if (pos_order(C, lo, hi) < 0) {
                std::size_t fa = face_of_region_[a.region];
                std::size_t fb = face_of_region_[b.region];
                assert(fa != NONE && fb != NONE &&
                       "positive-length boundary implies nonempty regions");
                if (fa != fb) {
                    std::size_t id = edges.size();
                    edges.push_back({fa, fb, false});
                    leg_events[a.arc][a.leg].push_back(id);
                    leg_events[b.arc][b.leg].push_back(id);
                }
            }
            if (pos_order(C, a.hi, b.hi) <= 0) ++i; else ++j;
        }
    }

    // Chord adjacencies ("they share a chord", tex 295); null-length
    // chords are zero-length S*-edges, contracted away (tex 286).
    std::vector<std::size_t> chord_edge(S.num_chords(), NONE);
    for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
        const Chord& c = S.chord(ci);
        assert(!c.dead && "[C91 §2.4]: normal form has no dead chords");
        if (c.is_null_length) continue;
        std::size_t fa = face_of_region_[c.region[0]];
        std::size_t fb = face_of_region_[c.region[1]];
        assert(fa != NONE && fb != NONE &&
               "a positive-length chord bounds two nonempty regions");
        assert(fa != fb && "[C91 §2.2]: a chord separates two regions");
        chord_edge[ci] = edges.size();
        edges.push_back({fa, fb, false});
    }

    // ── Rotations: walk each face's boundary cycle ───────────────
    // Lemma 2.2 ([C91 §2.2 tex 96]): a region's arcs in canonical table
    // order ARE its boundary cycle, with its chords interleaved where
    // consecutive arcs meet at chord endpoints ([C91 §2.4 tex 142]:
    // wrap-spanning arcs are single structures, so every junction is a
    // chord endpoint).  The dual embedding inherits its rotations from
    // these cycles.
    std::vector<std::vector<std::size_t>> rot(mu_);
    for (std::size_t f = 0; f < mu_; ++f) {
        std::size_t r = region_of_face_[f];
        const auto& arcs = arcs_of_region_[r];
        assert(!arcs.empty() && "a face's region has boundary arcs");
        std::vector<std::size_t> sorted_arcs(arcs);
        std::sort(sorted_arcs.begin(), sorted_arcs.end());

        std::size_t start = sorted_arcs[0];
        std::size_t cur = start;
        std::size_t emitted = 0;
        do {
            ++emitted;
            assert(emitted <= 2 * arcs.size() + 2 &&
                   "region boundary cycle must close");
            // The arc's abutment events, in traversal order per leg —
            // legs in cycle order ([C91 §2.4 tex 142]); LEFT legs run
            // with the stored order, RIGHT legs reversed.
            {
                ArcLeg lg[3];
                std::size_t n = S.arc(cur).legs(S.start_vertex,
                                                S.end_vertex, lg);
                for (std::size_t li = 0; li < n; ++li) {
                    const auto& evs = leg_events[cur][li];
                    if (lg[li].side == LEFT) {
                        for (std::size_t id : evs)
                            rot[f].push_back(id);
                    } else {
                        for (std::size_t k = evs.size(); k-- > 0; )
                            rot[f].push_back(evs[k]);
                    }
                }
            }

            // Find the chord of r ending this arc, if any.
            std::size_t via_chord = NONE;
            bool via_left_slot = false;
            for (std::size_t ci : S.node(r).incident_chords) {
                const Chord& c = S.chord(ci);
                for (bool left_slot : {true, false}) {
                    if (chord_before_arc(S, C, c, left_slot) != cur)
                        continue;
                    // The before-arc must border r on this side.
                    assert(via_chord == NONE &&
                           "an arc ends at exactly one chord endpoint");
                    via_chord = ci;
                    via_left_slot = left_slot;
                }
            }
            // [C91 §2.2 tex 96]: the boundary alternates arcs and exit
            // chords, so every arc ends at a chord endpoint of its
            // region — wrap-spanning arcs are single structures with
            // chord-bounded ends ([C91 §2.4 tex 142]).  (A chordless
            // submap has μ == 1 and never builds the dual graph.)
            assert(via_chord != NONE &&
                   "[C91 §2.2 tex 96]: every arc ends at a chord "
                   "endpoint of its region");
            if (chord_edge[via_chord] != NONE)
                rot[f].push_back(chord_edge[via_chord]);
            std::size_t next = chord_after_arc(S, C, S.chord(via_chord),
                                               !via_left_slot);
            assert(S.arc(next).region_node == r &&
                   "[C91 §2.2 tex 96]: the boundary cycle stays in the "
                   "region");
            cur = next;
        } while (cur != start);
        assert(emitted == arcs.size() &&
               "the boundary cycle visits every arc of the region once");
    }

    // ── Dedup: [C91 §3.4 tex 295] G is simple ────────────────────
    {
        std::vector<std::size_t> order(edges.size());
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
        auto key = [&](std::size_t e) {
            auto [x, y] = std::minmax(edges[e].fa, edges[e].fb);
            return std::pair<std::size_t, std::size_t>(x, y);
        };
        std::sort(order.begin(), order.end(),
                  [&](std::size_t a, std::size_t b) {
                      return key(a) < key(b) || (key(a) == key(b) && a < b);
                  });
        for (std::size_t i = 1; i < order.size(); ++i)
            if (key(order[i]) == key(order[i - 1]))
                edges[order[i]].dead = true;
    }
    std::vector<std::size_t> edge_map(edges.size(), NONE);
    std::vector<std::pair<std::size_t, std::size_t>> live_edges;
    for (std::size_t e = 0; e < edges.size(); ++e) {
        if (edges[e].dead) continue;
        edge_map[e] = live_edges.size();
        live_edges.emplace_back(edges[e].fa, edges[e].fb);
        dual_edges_.emplace_back(edges[e].fa, edges[e].fb);
    }
    std::vector<std::vector<std::size_t>> live_rot(mu_);
    for (std::size_t f = 0; f < mu_; ++f)
        for (std::size_t id : rot[f])
            if (edge_map[id] != NONE)
                live_rot[f].push_back(edge_map[id]);

    // [C91 §3.4 tex 297]: "The first payoff is that the number of edges
    // is at most 3μ − 6."
    assert((mu_ < 3 || live_edges.size() <= 3 * mu_ - 6) &&
           "[C91 §3.4 tex 297]: |E(G)| ≤ 3μ − 6");

    EmbeddedPlanarGraph G(mu_, live_edges, live_rot);

    // G is connected (S* is a connected subdivision, tex 295–297).
    {
        std::vector<bool> vis(mu_, false);
        std::vector<std::size_t> q{0};
        vis[0] = true;
        for (std::size_t qi = 0; qi < q.size(); ++qi)
            for (std::size_t h : G.incident_halves(q[qi])) {
                std::size_t w = G.half_to(h);
                if (!vis[w]) { vis[w] = true; q.push_back(w); }
            }
        assert(q.size() == mu_ &&
               "[C91 §3.4 tex 295]: the dual graph G is connected");
    }

    // [C91 §3.4 tex 297–304]: D* and the D_i partition, δ = 2/3.
    decomp_ = build_separator_decomposition(G);
    // [C91 §3.4 tex 308]: the query "naively check[s] all the regions
    // dual to nodes in D_i, which takes O(γ μ^{2/3}) time" — enumerate
    // each D_i's members once here (O(μ), within tex 297's O(μ log m)
    // preprocessing) so the query touches only |D_i| ≤ μ^{2/3} faces.
    subset_faces_.assign(decomp_.num_subsets, {});
    for (std::size_t f = 0; f < mu_; ++f) {
        if (decomp_.subset[f] == NONE)
            dstar_faces_.push_back(f);
        else
            subset_faces_[decomp_.subset[f]].push_back(f);
    }
}

void RayShootingStructure::build_vertical_line() {
    const Submap& S = *S_;
    const Polygon& C = *C_;

    // [C91 §3.4 tex 306]: "Take a vertical line passing to the right of
    // all the vertices of P, and intersect it with the chords of the
    // regions in S."  A chord crosses that line exactly when it runs
    // through infinity ([C91 §2.1 tex 70]); its crossing sits at the
    // chord's y, and the crossing chords are found by traversing S.
    std::vector<std::size_t> wrapped;
    for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
        const Chord& c = S.chord(ci);
        if (c.is_null_length) continue;
        if (chord_runs_through_infinity(C, c)) wrapped.push_back(ci);
    }

    if (wrapped.empty()) {
        // The whole line lies in one region: the one containing the
        // polar caps.  Locate it just above the global y-maximum of C
        // ([C91 §2 tex 47] input-table aggregate).
        std::size_t vt = C.max_y_vertex();
        if (C.is_endpoint(vt)) {
            // [C91 §2.1 tex 72] case 3: the boundary doubles back
            // around the endpoint; both sides of its incident edge
            // bound the polar region there.
            region_infinity_ = (vt == 0)
                ? S.arc(S.start_arc).region_node
                : S.arc(S.end_arc).region_node;
        } else {
            assert(C.is_y_extremum(vt) &&
                   "the interior global y-max is a local extremum");
            // Outside (away-from-wedge) side of the incoming edge
            // vt−1 → vt, via the branch x-slopes at the vertex
            // ([C91 §2.1 tex 72]; same computation as the junction
            // synthesis in rebuild_submap).
            const Point& u = C.vertex(vt - 1);
            const Point& v = C.vertex(vt);
            const Point& w = C.vertex(vt + 1);
            const double t_prev = (u.x - v.x) * (v.y - w.y);
            const double t_next = (w.x - v.x) * (v.y - u.y);
            assert(t_prev != t_next &&
                   "[C91 §2 tex 47]: distinct branch slopes at a vertex "
                   "of a simple curve");
            const bool prev_left = t_prev < t_next;
            // e_in = vt−1 → vt ascends into the max; its −x face is
            // LEFT.  The wedge (inside) lies toward the next branch:
            // inside = LEFT iff the next branch is left of the previous
            // one... equivalently outside = the side of the PREV branch.
            const Side outside_in = prev_left ? LEFT : RIGHT;
            auto ids = S.double_identify(vt - 1, symbolic_y_of(v), C);
            std::size_t found = NONE;
            for (std::size_t ai : ids)
                if (arc_covers(S, C, ai, vt - 1, outside_in)) {
                    found = ai;
                    break;
                }
            assert(found != NONE &&
                   "an arc covers the outside face at the global maximum");
            region_infinity_ = S.arc(found).region_node;
        }
        assert(face_of_region_[region_infinity_] != NONE &&
               "the polar region is nonempty");
        return;
    }

    // [C91 §3.4 tex 306]: "to split up the line and identify the
    // regions cut by each segment can be done by traversing G ...
    // Since the regions cut correspond to nodes of G lying on a path,
    // sorting the intersections comes for free, and all the work can be
    // done in O(μ) time."  The crossings' regions form a simple path in
    // the submap tree (a horizontal circle cannot revisit a region:
    // that would close a walk with distinct edges in a tree), so we
    // orient each crossing (regions below/above at the far field, read
    // off the endpoint adjacencies) and chain them bottom-up.
    struct Cross {
        std::size_t chord, below, above;
    };
    std::vector<Cross> cx;
    std::vector<std::size_t> at_region(S_->num_nodes(), 0);
    for (std::size_t ci : wrapped) {
        Cross c{ci, NONE, NONE};
        chord_vertical_regions(S, C, ci, &c.below, &c.above);
        cx.push_back(c);
        ++at_region[c.below];
        ++at_region[c.above];
        assert(at_region[c.below] <= 2 && at_region[c.above] <= 2 &&
               "[C91 §3.4 tex 306]: the cut regions lie on a path");
    }
    // Chain: start at the crossing whose below-region meets no other
    // crossing (the south cap), then follow above-regions upward.
    std::vector<std::size_t> next_by_below(S_->num_nodes(), NONE);
    for (std::size_t k = 0; k < cx.size(); ++k) {
        assert(next_by_below[cx[k].below] == NONE);
        next_by_below[cx[k].below] = k;
    }
    std::size_t first = NONE;
    for (std::size_t k = 0; k < cx.size(); ++k)
        if (at_region[cx[k].below] == 1) {
            assert(first == NONE &&
                   "exactly one south-cap crossing exists");
            first = k;
        }
    assert(first != NONE);
    std::size_t curk = first;
    while (true) {
        const Cross& c = cx[curk];
        LineCrossing lc;
        lc.chord = c.chord;
        lc.y = S.chord(c.chord).symbolic_y();
        lc.region_below = c.below;
        lc.region_above = c.above;
        if (!line_.empty()) {
            assert(line_.back().region_above == c.below &&
                   "[C91 §3.4 tex 306]: consecutive segments share their "
                   "crossing chord's region");
            assert(symbolic_y_less(line_.back().y, lc.y) &&
                   "[C91 §3.4 tex 306]: sorting the intersections comes "
                   "for free");
        }
        line_.push_back(lc);
        std::size_t nk = next_by_below[c.above];
        if (nk == NONE) break;
        curk = nk;
    }
    assert(line_.size() == cx.size() &&
           "every crossing lies on the single bottom-to-top path");
}

// ════════════════════════════════════════════════════════════════
//  Query — [C91 §3.4 tex 306–308]
// ════════════════════════════════════════════════════════════════

void RayShootingStructure::regions_at_boundary(
        std::size_t edge, Side side, SymbolicY y,
        std::vector<std::size_t>& out) const {
    const Polygon& C = *C_;
    const auto& list = (side == LEFT) ? left_ivals_ : right_ivals_;
    assert(!list.empty());
    Pos pos{edge, y};

    // Last interval with lo ≤ pos ([C91 §3.4 tex 308]: "might require a
    // binary search among a large collection of collinear edges").
    std::size_t lo = 0, hi = list.size();
    while (lo < hi) {
        std::size_t mid = (lo + hi) / 2;
        Pos mlo{list[mid].lo_edge, list[mid].lo_y};
        if (pos_order(C, mlo, pos) <= 0) lo = mid + 1;
        else hi = mid;
    }
    auto contains = [&](std::size_t i) {
        Pos ilo{list[i].lo_edge, list[i].lo_y};
        Pos ihi{list[i].hi_edge, list[i].hi_y};
        return pos_order(C, ilo, pos) <= 0 && pos_order(C, pos, ihi) <= 0;
    };
    auto push = [&](std::size_t r) {
        for (std::size_t x : out)
            if (x == r) return;
        out.push_back(r);
    };
    // Every interval containing pos lies in a contiguous block ending at
    // lo − 1 (list[lo].lo > pos by the search, so contains(lo) is
    // impossible; intervals are interior-disjoint per side, so the walk
    // stops at the first non-containing predecessor).  Walk the whole
    // block: at an interval junction BOTH flanking intervals contain pos
    // (closed ends), and zero-length intervals may stack there.  [C91
    // §2.4 tex 144]: ≤ 6 arcs pass through one ∂C point, so this is O(1).
    bool any = false;
    for (std::size_t i = lo; i-- > 0 && contains(i);) {
        push(list[i].region);
        any = true;
    }
    if (!any) {
        // The contact falls in a symbolic gap between consecutive
        // intervals (a contracted point — zero-length arcs stacked at
        // one location, [C91 §3.4 tex 286]); both flanking faces are
        // candidates.
        if (lo > 0) push(list[lo - 1].region);
        if (lo < list.size()) push(list[lo].region);
    }
    assert(!out.empty() &&
           "[C91 §2.4 tex 144]: every ∂C contact identifies a region");
}

RayHit RayShootingStructure::shoot_toward_boundary(Point p,
                                                   Side dir) const {
    const Submap& S = *S_;
    const Polygon& C = *C_;
    SymbolicY sy{p.y, p.index};

    auto to_rayhit = [&](const Contact& c) {
        RayHit h;
        if (!c.hit) return h;
        h.hit = true;
        h.x = c.x;
        h.y = p.y;
        h.edge = c.edge;
        h.side = c.side;
        h.wrapped = c.wrapped;
        return h;
    };

    if (mu_ <= 1) {
        // [C91 §3.4 tex 297]: "If μ = 1, then ray-shooting can be done
        // trivially in O(m) time."
        Contact best;
        scan_edge_range(C, 0, C.num_edges() - 1, p, sy, dir, best);
        return to_rayhit(best);
    }

    // Candidate D_i regions of the query, dedup'd; SoS degeneracies
    // (the ray collinear with a chord) can yield two flanking regions —
    // scanning both subsets stays within the O(γ μ^{2/3}) budget.
    std::vector<std::size_t> subsets;
    auto add_subset = [&](std::size_t region) {
        std::size_t f = face_of_region_[region];
        assert(f != NONE);
        std::size_t sub = decomp_.subset[f];
        // A candidate region in D* contributes nothing new: its edges
        // were already scanned in step 1 (this covers the source lying
        // inside a D* region that the ray never leaves).
        if (sub == NONE) return;
        for (std::size_t s : subsets)
            if (s == sub) return;
        subsets.push_back(sub);
    };

    // 1. [C91 §3.4 tex 306]: "Our first task is to shoot within each
    // region that is dual to a node of D*."
    Contact dstar_best;
    for (std::size_t f : dstar_faces_)
        scan_region(S, C, arcs_of_region_[region_of_face_[f]], p, sy,
                    dir, dstar_best);

    if (dstar_best.hit) {
        // [C91 §3.4 tex 306–308]: "Let R be the last region of S
        // traversed before the first hit.  To identify R can be done by
        // double identification, followed by checking the local
        // orientation of the hit" — R is the region whose boundary
        // covers the struck ∂C side at the contact.
        std::vector<std::size_t> near;
        regions_at_boundary(dstar_best.edge, dstar_best.side, sy, near);
        bool all_dstar = true;
        for (std::size_t r : near) {
            std::size_t f = face_of_region_[r];
            assert(f != NONE && "a struck arc bounds a nonempty region");
            if (decomp_.subset[f] == NONE) continue;   // R ∈ D*
            all_dstar = false;
            add_subset(r);
        }
        if (all_dstar) {
            // "If R is a region dual to a node v of D*, then the
            // starting point of the ray lies in R ... and we are
            // trivially done" (tex 306).
            return to_rayhit(dstar_best);
        }
    } else {
        // 2. [C91 §3.4 tex 308]: "assume now that the ray of light hits
        // no region dual to a node in D*.  Then the ray-shooting takes
        // place entirely within the regions dual to the nodes of a
        // single D_i.  To find out which one, we shoot toward the
        // vertical line and find which segment of the line is hit.
        // This takes O(log μ) time by binary search."
        if (line_.empty()) {
            add_subset(region_infinity_);
        } else {
            // First crossing strictly above sy.
            std::size_t lo = 0, hi = line_.size();
            while (lo < hi) {
                std::size_t mid = (lo + hi) / 2;
                if (symbolic_y_less(sy, line_[mid].y)) hi = mid;
                else lo = mid + 1;
            }
            // lo = index of the first crossing STRICTLY above sy (the
            // search advances past equal keys, so an sy equal to a
            // crossing's y always lands at lo − 1).
            if (lo > 0 && symbolic_y_equal(sy, line_[lo - 1].y)) {
                add_subset(line_[lo - 1].region_below);
                add_subset(line_[lo - 1].region_above);
            } else if (lo == 0) {
                add_subset(line_[0].region_below);
            } else if (lo == line_.size()) {
                add_subset(line_.back().region_above);
            } else {
                assert(line_[lo - 1].region_above ==
                       line_[lo].region_below);
                add_subset(line_[lo].region_below);
            }
        }
    }

    // 3. [C91 §3.4 tex 308]: "We can find w, and, from there, answer
    // the ray-shooting query, by first finding D_i, which takes
    // constant time since we know R, and then naively checking all the
    // regions dual to nodes in D_i, which takes O(γ μ^{2/3}) time."
    Contact best = dstar_best;
    for (std::size_t sub : subsets) {
        const auto& members = subset_faces_[sub];
#ifndef NDEBUG
        // [C91 §3.4 tex 304]: each |D_i| ≤ μ^{2/3} (128-bit exact:
        // the cube overflows 64 bits at realistic node counts).
        __extension__ typedef unsigned __int128 u128;
        assert((u128)members.size() * members.size() * members.size() <=
               (u128)mu_ * mu_);
#endif
        for (std::size_t f : members)
            scan_region(S, C, arcs_of_region_[region_of_face_[f]], p, sy,
                        dir, best);
    }

    if (!best.hit) {
        // A horizontal circle within C's perturbed y-range crosses the
        // connected curve C; no contact means sy is above the global
        // maximum or below the global minimum.
        assert((symbolic_y_greater(sy, symbolic_y_of(
                    C.vertex(C.max_y_vertex()))) ||
                symbolic_y_less(sy, symbolic_y_of(
                    C.vertex(C.min_y_vertex())))) &&
               "[C91 §2.1 tex 70]: a wrapping ray inside C's y-range "
               "must hit C");
    }
    return to_rayhit(best);
}

// ════════════════════════════════════════════════════════════════
//  SubmapRayShooter — the [C91 §3.0(i) tex 169] interface
// ════════════════════════════════════════════════════════════════

namespace {

// Is the ∂C position (edge, side, y) on the target subarc α'?
// [C91 §3.0(i) tex 169]: α' is endpoint-exact — "specified by its two
// endpoints (along with two pointers to the input table to indicate
// the names of the edges of P that contain these two endpoints as well
// as two flags indicating which side of ∂P is to be understood)" — so
// a candidate on a shared boundary edge is clipped at α''s exact
// endpoint ys; wrapped targets decompose into their legs per
// [C91 §2.4 tex 142] (oracle.h::subarc_contains_point).
bool subarc_covers(const Submap& S, const Polygon& C, const Subarc& t,
                   std::size_t edge, Side side, SymbolicY y) {
    assert(S.start_vertex != NONE && S.end_vertex != NONE &&
           "[C91 §2.4(iii)]: C endpoints must be set");
    return subarc_contains_point(t, C, edge, side, y,
                                 S.start_vertex, S.end_vertex);
}

} // namespace

RayHit SubmapRayShooter::shoot(Point p, Side direction,
                               std::size_t arc_idx,
                               const Subarc& target) const {
    // [C91 §3.0(i) tex 166]: the arc pointer identifies α; the report
    // concerns α' ⊆ α only.
    assert(arc_idx < S_->num_arcs() && !S_->arc(arc_idx).dead);
    assert_subarc_clockwise(target);

    // Single first-contact of C, post-filtered to α' — the obstacle-C
    // semantics documented in the class comment (ray_shooting.h) and
    // TODO.md item 4; the [C91 §4.1 tex 341] per-piece structures
    // realize the tex-169 "no obstacle except α'" contract.  The hit
    // rides the ray's perturbed level ([C91 §2 tex 47]).
    RayHit h = impl_.shoot_toward_boundary(p, direction);
    if (!h.hit) return RayHit{};
    if (!subarc_covers(*S_, *C_, target, h.edge, h.side,
                       SymbolicY{p.y, p.index}))
        return RayHit{};
    return h;
}

} // namespace chazelle
