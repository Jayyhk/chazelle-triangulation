// src/merge/conformality.cpp — [C91 §3.2 tex 236–272]: Restoring Conformality.

#include "conformality.h"
#include "fusion.h"
#include "../submap/shielding.h"

#include <algorithm>

namespace chazelle {

namespace {

// ════════════════════════════════════════════════════════════════
//  Clockwise ∂C tour positions
// ════════════════════════════════════════════════════════════════

// [C91 §2.4(iii) tex 138]: canonical clockwise ∂C order — LEFT faces on
// ascending edges, then RIGHT faces on descending edges; within an edge
// the traversal follows the edge (LEFT) or reverses it (RIGHT).  Same
// keying as rebuild_submap's endpoint sort ([C91 §3.1 tex 226]).
struct TourPos {
    std::size_t tp = 0;       // LEFT: edge; RIGHT: 2n−1−edge
    SymbolicY y{};
    bool y_ascending = true;  // traversal y-direction on this edge face
};

TourPos tour_pos(const Polygon& C, std::size_t edge, Side side,
                 SymbolicY y) {
    assert(edge < C.num_edges());
    TourPos t;
    t.tp = (side == LEFT) ? edge : (2 * C.num_edges() - 1 - edge);
    const auto& e = C.edge(edge);
    bool asc = symbolic_y_less(symbolic_y_of(C.vertex(e.start_idx)),
                               symbolic_y_of(C.vertex(e.end_idx)));
    t.y_ascending = (side == LEFT) ? asc : !asc;
    t.y = y;
    return t;
}

int tour_cmp(const TourPos& a, const TourPos& b) {
    if (a.tp != b.tp) return a.tp < b.tp ? -1 : 1;
    int c = symbolic_y_compare(a.y, b.y);
    if (c == 0) return 0;
    return a.y_ascending ? c : -c;
}

// x appears no later than y when walking clockwise from `base`
// (cyclic order; base itself has rank 0).
bool at_or_before_from(const TourPos& base, const TourPos& x,
                       const TourPos& y) {
    bool xw = tour_cmp(x, base) < 0;   // x wrapped past the tour origin
    bool yw = tour_cmp(y, base) < 0;
    if (xw != yw) return yw;
    return tour_cmp(x, y) <= 0;
}

// x ∈ [lo, hi] walking clockwise from lo (closed cyclic interval).
bool in_closed_cyclic(const TourPos& x, const TourPos& lo,
                      const TourPos& hi) {
    return at_or_before_from(lo, x, hi);
}

// x ∈ [lo, hi) — half-open; lo itself is a member, hi is not.
bool in_half_open_cyclic(const TourPos& x, const TourPos& lo,
                         const TourPos& hi) {
    return at_or_before_from(lo, x, hi) && tour_cmp(x, hi) != 0;
}

// ════════════════════════════════════════════════════════════════
//  Fused-arc endpoint positions
// ════════════════════════════════════════════════════════════════

TourPos fused_arc_start(const Submap& S, const Polygon& C,
                        std::size_t arc_idx) {
    const Arc& a = S.arc(arc_idx);
    return tour_pos(C, a.first_edge, a.first_side,
                    S.arc_start_symbolic_y(arc_idx, C));
}

TourPos fused_arc_end(const Submap& S, const Polygon& C,
                      std::size_t arc_idx) {
    const Arc& a = S.arc(arc_idx);
    return tour_pos(C, a.last_edge, a.last_side,
                    S.arc_end_symbolic_y(arc_idx, C));
}

// Does the fused table arc contain the ∂C point (edge, side, y)?
// [C91 §2.4 tex 142]: a wrap-spanning arc is one structure — its edge
// coverage is checked per leg, and its tour interval is cyclic when it
// crosses C's start turnaround (the tour origin); the end turnaround is
// linear in tour order (LEFT n−1 → RIGHT n−1 are adjacent positions).
bool fused_arc_contains(const Submap& S, const Polygon& C,
                        std::size_t arc_idx,
                        std::size_t edge, Side side, SymbolicY y) {
    const Arc& a = S.arc(arc_idx);
    assert(S.start_vertex != NONE && S.end_vertex != NONE &&
           "[C91 §2.4(iii)]: S's C endpoints must be set");
    if (!a.covers(edge, side, S.start_vertex, S.end_vertex))
        return false;
    TourPos x = tour_pos(C, edge, side, y);
    TourPos s = fused_arc_start(S, C, arc_idx);
    TourPos e = fused_arc_end(S, C, arc_idx);
    if (tour_cmp(s, e) <= 0)
        return tour_cmp(s, x) <= 0 && tour_cmp(x, e) <= 0;
    // Crosses the tour origin (C's start turnaround).
    return tour_cmp(s, x) <= 0 || tour_cmp(x, e) <= 0;
}

} // namespace

// ════════════════════════════════════════════════════════════════
//  compute_arc_provenance — [C91 §3.2 tex 244]
// ════════════════════════════════════════════════════════════════

std::vector<ArcProvenance> compute_arc_provenance(
        const Submap& S, const Polygon& C,
        const Submap& S1, const Polygon& C1,
        const Submap& S2, const Polygon& C2) {
    const std::size_t n1e = C1.num_edges();
    assert(C.num_edges() == C1.num_edges() + C2.num_edges() &&
           "[C91 §3 tex 160]: C = C₁ ∪ C₂ shares one vertex");

    std::vector<ArcProvenance> out(S.num_arcs());
    for (std::size_t ai = 0; ai < S.num_arcs(); ++ai) {
        const Arc& a = S.arc(ai);
        assert(!a.dead &&
               "[C91 §3.1 tex 226]: fused submap is freshly built (no dead arcs)");
        // [C91 §2.4 tex 142]: a fused arc may double-back around an
        // endpoint of C — but never around BOTH: a double-wrap arc
        // would cover a whole ∂C side including the junction, where
        // chords were added at every companion ([C91 §3.1 tex 224]).
        assert(!(a.first_side == a.last_side && a.wraps()) &&
               "[C91 §3.2 tex 244]: no fused arc double-wraps (junction "
               "chords break both sides)");

        // [C91 §3.2 tex 244]: an arc of S "cannot overlap both ∂C₁ and
        // ∂C₂" — chords were added at every junction companion, so no
        // arc straddles the junction.  (A wrap arc's two legs flank one
        // C endpoint, so first/last edges agree on the side of the
        // junction.)
        const bool on_c1 = a.first_edge < n1e;
        assert((a.last_edge < n1e) == on_c1 &&
               "[C91 §3.2 tex 244]: fused arc must lie on one ∂Cᵢ only");

        const Submap&  Si = on_c1 ? S1 : S2;
        const Polygon& Ci = on_c1 ? C1 : C2;
        const std::size_t off = on_c1 ? 0 : n1e;

        const std::size_t edge_i = a.first_edge - off;
        const SymbolicY sy = S.arc_start_symbolic_y(ai, C);

        // [C91 §2.4 tex 144]: all Sᵢ arcs through the arc's start point.
        auto cands = Si.double_identify(edge_i, sy, Ci);
        assert(cands.count >= 1 &&
               "[C91 §2.4 tex 144]: ∂Cᵢ point must lie on some Sᵢ arc");

        // The containing arc covers the point AND continues clockwise
        // from it — reject candidates on the wrong ∂C side and those
        // ENDING exactly there (arcs tile ∂Cᵢ, so exactly one covers
        // forward).  Since every Sᵢ chord endpoint is still a chord
        // endpoint of S ([C91 §3.2 tex 244]), the surviving candidate
        // contains the whole fused arc.
        assert(Si.start_vertex != NONE && Si.end_vertex != NONE &&
               "[C91 §2.4(iii)]: Sᵢ's C endpoints must be set");
        std::size_t chosen = NONE;
        std::size_t covering_ender = NONE;
        for (std::size_t k = 0; k < cands.count; ++k) {
            std::size_t si_arc = cands.arcs[k];
            const Arc& sa = Si.arc(si_arc);
            if (!sa.covers(edge_i, a.first_side,
                           Si.start_vertex, Si.end_vertex)) continue;
            bool ends_here =
                sa.last_edge == edge_i && sa.last_side == a.first_side &&
                symbolic_y_equal(Si.arc_end_symbolic_y(si_arc, Ci), sy);
            if (ends_here) {
                covering_ender = si_arc;
                continue;
            }
            assert(chosen == NONE &&
                   "[C91 §2.4(iii)]: arcs tile ∂Cᵢ — exactly one covers "
                   "clockwise-forward from any point");
            chosen = si_arc;
        }
        // [C91 §2.2 tex 96]: a zero-length fused arc is a single ∂C
        // point with no forward extent — an Sᵢ arc ENDING exactly there
        // still contains it.
        if (chosen == NONE && a.edge_count == 0)
            chosen = covering_ender;
        assert(chosen != NONE &&
               "[C91 §3.2 tex 244]: every fused arc lies inside an Sᵢ arc");
        out[ai] = ArcProvenance{on_c1, chosen};
    }
    return out;
}

// ════════════════════════════════════════════════════════════════
//  fused_region_cycle — [C91 §3.2 tex 238]
// ════════════════════════════════════════════════════════════════

FusedRegionCycle fused_region_cycle(const Submap& S, const Polygon& C,
                                    std::size_t region,
                                    const std::vector<std::size_t>& arcs_of_region) {
    assert(region < S.num_nodes() && !S.node(region).dead);

    // [C91 §2.2 tex 96] Lemma 2.2: the clockwise ∂C order of a region's
    // arcs is its boundary order (counterclockwise around the region).
    struct Keyed { std::size_t arc; TourPos pos; };
    std::vector<Keyed> sorted;
    sorted.reserve(arcs_of_region.size());
    for (std::size_t ai : arcs_of_region) {
        assert(ai < S.num_arcs() && !S.arc(ai).dead &&
               S.arc(ai).region_node == region &&
               "[C91 §3.2]: inventory entries must be live arcs of the region");
        sorted.push_back({ai, fused_arc_start(S, C, ai)});
    }
    assert(!sorted.empty() && "[C91 §2.2 tex 96]: region has ≥ 1 arc");
    std::sort(sorted.begin(), sorted.end(),
              [](const Keyed& a, const Keyed& b) {
                  return tour_cmp(a.pos, b.pos) < 0;
              });
#ifndef NDEBUG
    // [C91 §2.2 tex 96]: consecutive cycle arcs are separated by an
    // exit chord whose far endpoint is a distinct ∂C position, so two
    // arcs of one region can never start at the same TourPos.  A tie
    // would silently permute the boundary cycle and corrupt the
    // nonconsecutiveness test and insert_chord's chain walk.
    for (std::size_t i = 1; i < sorted.size(); ++i)
        assert(tour_cmp(sorted[i - 1].pos, sorted[i].pos) != 0 &&
               "[C91 §2.2 tex 96]: region arcs must start at distinct "
               "∂C positions");
#endif

    // [C91 §2.4 tex 142]: a wrap-spanning arc is ONE arc-structure that
    // double-backs around C's endpoint — the table never splits it, so
    // the sorted arc-structures ARE the boundary cycle.
    FusedRegionCycle cycle;
    for (const Keyed& k : sorted) {
        assert(cycle.count < FusedRegionCycle::MAX_ARCS &&
               "[C91 §3.2 tex 238]: fused region arc count is bounded "
               "(≤ 2 runs of constant length)");
        // [C91 §2.2 tex 96]: "some arcs may be of zero length."  Under
        // [C91 §2.1] every P-edge is nonnull, so edge_count == 0 ⟺ the
        // arc has zero geometric extent.
        cycle.arcs[cycle.count++] =
            CycleArc{k.arc, S.arc(k.arc).edge_count == 0};
    }
    return cycle;
}

namespace {

// [C91 §2.1 tex 72]: is the ∂C point (edge, side) at vertex `vidx` one of
// a y-extremum's INSIDE-pair duplicates?  "one of these pairs, the one on
// the 'inside' of the turn, gives rise to a chord of null length" — such
// a point sees only its distance-0 sibling and can never source a chord
// to another arc, so Lemma 3.2 never shoots from it.  The inside face
// computation mirrors rebuild_submap's junction null-chord synthesis
// ([C91 §3.1 tex 224], fusion.cpp).
bool is_inside_companion(const Polygon& C, std::size_t edge, Side side,
                         std::size_t vidx) {
    if (vidx == 0 || vidx + 1 >= C.num_vertices()) return false;  // case 3
    const Point& u = C.vertex(vidx - 1);
    const Point& v = C.vertex(vidx);
    const Point& w = C.vertex(vidx + 1);
    if (!is_local_y_extremum(u, v, w)) return false;              // case 1
    const bool is_max = point_y_above(v, u) && point_y_above(v, w);

    // Branch x-order infinitesimally off v ([C91 §2 tex 47], exact
    // cross-multiplied slopes as in fusion.cpp's junction synthesis).
    const double t_prev = is_max ? (u.x - v.x) * (v.y - w.y)
                                 : (u.x - v.x) * (w.y - v.y);
    const double t_next = is_max ? (w.x - v.x) * (v.y - u.y)
                                 : (w.x - v.x) * (u.y - v.y);
    assert(t_prev != t_next &&
           "[C91 §2 tex 47]: extremum branches have distinct x-slopes");
    const bool prev_left = t_prev < t_next;

    auto minus_x_face = [&](std::size_t e) -> Side {
        const auto& ed = C.edge(e);
        bool asc = symbolic_y_less(symbolic_y_of(C.vertex(ed.start_idx)),
                                   symbolic_y_of(C.vertex(ed.end_idx)));
        return asc ? LEFT : RIGHT;
    };
    auto plus_x_face = [&](std::size_t e) -> Side {
        return minus_x_face(e) == LEFT ? RIGHT : LEFT;
    };

    // The wedge (inside of the turn) lies between the two branches: the
    // prev edge's inside face looks toward the next branch and vice
    // versa.
    Side inside_next = prev_left ? minus_x_face(vidx) : plus_x_face(vidx);
    Side inside_prev = prev_left ? plus_x_face(vidx - 1)
                                 : minus_x_face(vidx - 1);
    if (edge == vidx && side == inside_next) return true;
    if (edge == vidx - 1 && side == inside_prev) return true;
    return false;
}

} // namespace

// ════════════════════════════════════════════════════════════════
//  local_shoot_fused — [C91 §3.2 tex 244/246]
// ════════════════════════════════════════════════════════════════

RayHit local_shoot_fused(Point p, SymbolicY p_y, Side direction,
                         const FusedRegionCycle& cycle,
                         const FusedShootContext& ctx,
                         bool require_hit) {
    assert(ctx.S && ctx.C && ctx.C1 && ctx.C2 && ctx.ray1 && ctx.ray2 &&
           ctx.provenance);
    const Submap& S = *ctx.S;
    const Polygon& C = *ctx.C;
    const std::size_t n1e = ctx.C1->num_edges();

    // [C91 §2 tex 47]: the ray's perturbed height is p_y; carry the tag
    // in Point.index for the oracle's symbolic comparisons.
    p.index = p_y.tag;

    RayHit best;
    best.hit = false;

    for (std::size_t li = 0; li < cycle.count; ++li) {
        {
            const std::size_t ai = cycle.arcs[li].arc;
            const Arc& a = S.arc(ai);
            const ArcProvenance& pr = (*ctx.provenance)[ai];
            const std::size_t off = pr.on_c1 ? 0 : n1e;
            const RayShootingOracle& oracle =
                pr.on_c1 ? *ctx.ray1 : *ctx.ray2;

            // [C91 §3.2 tex 244]: shoot into the Sᵢ arc containing this
            // fused arc, restricted to the fused arc's span — one call
            // per arc-structure; a wrap-spanning target keeps its
            // double-backing side flags ([C91 §2.4 tex 142] /
            // [C91 §3.0(i) tex 169]).
            Subarc target;
            target.first_edge = a.first_edge - off;
            target.first_side = a.first_side;
            target.last_edge = a.last_edge - off;
            target.last_side = a.last_side;
            // [C91 §3.0(i) tex 169]: α' is "specified by its two
            // endpoints" — exact ys are frame-independent; only the
            // edges carry the −off shift into Cᵢ's frame.
            target.first_y = S.arc_start_symbolic_y(ai, C);
            target.last_y = S.arc_end_symbolic_y(ai, C);
            assert_subarc_clockwise(target);

            RayHit hit = oracle.shoot(p, direction, pr.arc_in_si, target);
            if (!hit.hit) continue;

            // [C91 §3.0(i) tex 169 + §2.1 tex 70]: a direct hit lies
            // forward of p; a wrapped (through-infinity) hit lies behind.
            double d = (direction == LEFT) ? (p.x - hit.x)
                                           : (hit.x - p.x);
            if (hit.wrapped)
                assert(d <= 0.0 &&
                       "[C91 §2.1 tex 70]: a wrapped hit lies at or behind "
                       "the source in the travel direction (d == 0 is the "
                       "source's own duplicate met after a full wrap)");
            else
                assert(d > 0.0 &&
                       "[C91 §3.0(i) tex 169]: a direct hit lies strictly "
                       "forward (same-position hits are only reachable "
                       "through the wrap)");

            hit.edge += off;    // back to C's frame

            // [C91 §3.0(i) tex 169]: the Subarc target is endpoint-
            // exact, so the oracle's report lies on the fused arc by
            // contract; [C91 Lemma 2.1 tex 77] pins the hit to a wall
            // of the region, so the arc's exact span must contain it.
            assert(fused_arc_contains(S, C, ai, hit.edge, hit.side, p_y) &&
                   "[C91 §3.0(i) tex 169 + Lemma 2.1]: the report lies "
                   "on α′ (the fused arc's exact span)");

            hit.hit_arc_idx = ai;   // fused-table attribution

            if (!best.hit) {
                best = hit;
                continue;
            }
            // [C91 §2.1 tex 70]: nearest = lexicographic (wrapped, d) —
            // every direct hit precedes every wrapped one; within a
            // class, smaller signed distance is nearer (for wrapped hits
            // d < 0, and the most negative is met first after the wrap).
            double bd = (direction == LEFT) ? (p.x - best.x)
                                            : (best.x - p.x);
            if (hit.wrapped != best.wrapped) {
                if (!hit.wrapped) best = hit;
            } else if (d < bd) {
                best = hit;
            } else if (d == bd) {
                // [C91 §2.1 tex 72]: double-boundary disambiguation —
                // same first-face rule as [C91 §3.1] local_shoot
                // (fusion.cpp): the ray strikes the wall whose free side
                // opposes its travel.  This holds at distance 0 too
                // (after a full wrap the source's own wall still faces
                // away from the ray and is never struck), so there is no
                // distance-0 inversion.
                [[maybe_unused]] auto opposes_ray = [&](const RayHit& h) {
                    assert(h.edge < C.num_edges());
                    const auto& e = C.edge(h.edge);
                    bool asc = symbolic_y_less(
                        symbolic_y_of(C.vertex(e.start_idx)),
                        symbolic_y_of(C.vertex(e.end_idx)));
                    Side minus_x_face = asc ? LEFT : RIGHT;
                    Side plus_x_face = asc ? RIGHT : LEFT;
                    return h.side == ((direction == RIGHT) ? minus_x_face
                                                           : plus_x_face);
                };
                assert(opposes_ray(hit) && opposes_ray(best) &&
                       "[C91 §3.0(i) tex 169]: a reported hit is a wall "
                       "opposing the ray's travel");
                if (hit.edge != best.edge) {
                    // Two opposing walls at one point of C: the shared
                    // vertex of the two edges, at the ray's exact
                    // symbolic level ([C91 §2 tex 47]: SoS admits no
                    // other coincidence).  At a y-extremum exactly one
                    // face is the inside-of-turn duplicate ([C91 §2.1
                    // tex 72]), which sees only its distance-0 sibling —
                    // the outer wall shields the wedge, so the non-inside
                    // face is struck.  At a non-extremum both faces name
                    // the same companion point; keep the first found.
                    const std::size_t vidx = std::max(hit.edge, best.edge);
                    assert(vidx == std::min(hit.edge, best.edge) + 1 &&
                           "[C91 §2 tex 47]: same-position walls lie on "
                           "edges sharing one vertex of C");
                    assert(symbolic_y_equal(
                               symbolic_y_of(C.vertex(vidx)), p_y) &&
                           "[C91 §2 tex 47]: the tie point is the shared "
                           "vertex at the ray's level");
                    const bool hit_inside =
                        is_inside_companion(C, hit.edge, hit.side, vidx);
                    const bool best_inside =
                        is_inside_companion(C, best.edge, best.side, vidx);
                    if (best_inside && !hit_inside) best = hit;
                }
            }
        }
    }

    // [C91 §3.1 tex 181 / §2.2 Lemma 2.1]: regions are closed under
    // visibility — shooting from a boundary point of the region hits.
    if (require_hit)
        assert(best.hit &&
               "[C91 §3.2 tex 244]: local shoot within a fused region must hit");
    return best;
}

// ════════════════════════════════════════════════════════════════
//  find_visible_point — [C91 §3.2 Lemma 3.2]
// ════════════════════════════════════════════════════════════════

namespace {

// Does chord (y, {e1,s1}—{e2,s2}) duplicate an existing chord of the
// region?  [C91 §2.1 tex 70]: a ∂C point's horizontal visibility is
// unique, so a site at an existing chord endpoint can only re-derive
// that same chord — which must not be added again ([C91 §2.2 tex 102]:
// the dual graph is a tree, no duplicate chords).
bool duplicates_region_chord(const Submap& S, std::size_t region,
                             SymbolicY y,
                             std::size_t e1, Side s1,
                             std::size_t e2, Side s2) {
    for (std::size_t ci : S.node(region).incident_chords) {
        const Chord& c = S.chord(ci);
        if (c.dead) continue;
        if (!symbolic_y_equal(c.symbolic_y(), y)) continue;
        bool fwd = c.left_edge == e1 && c.left_side == s1 &&
                   c.right_edge == e2 && c.right_side == s2;
        bool rev = c.left_edge == e2 && c.left_side == s2 &&
                   c.right_edge == e1 && c.right_side == s1;
        if (fwd || rev) return true;
    }
    return false;
}

struct SiteShot {
    bool success = false;   // hit lands on A₂ and is a NEW chord
    RayHit hit{};
};

// [C91 §2.1 tex 70]: a ∂C point's horizontal visibility is unique.  If
// the point already IS a chord endpoint of the region, its visibility is
// realized by that chord — a shot from it would just retrace the chord
// (as fusion's case (i) reads t off the stored chord, [C91 §3.1 tex 220]).
struct RealizedVisibility {
    bool realized = false;
    std::size_t other_edge = NONE;  // the chord's far endpoint
    Side other_side = LEFT;
};

RealizedVisibility region_chord_at(const Submap& S, std::size_t region,
                                   std::size_t edge, Side side,
                                   SymbolicY y) {
    for (std::size_t ci : S.node(region).incident_chords) {
        const Chord& c = S.chord(ci);
        if (c.dead) continue;
        if (!symbolic_y_equal(c.symbolic_y(), y)) continue;
        if (c.left_edge == edge && c.left_side == side)
            return {true, c.right_edge, c.right_side};
        if (c.right_edge == edge && c.right_side == side)
            return {true, c.left_edge, c.left_side};
    }
    return {};
}

// [C91 §3.2 tex 246]: "we can efficiently check, in time O(f(γ₂)),
// whether a given vertex of A₁ qualifies as the point v sought" — shoot
// from the ∂C point (edge_c, side, y) within the region and test whether
// the point it sees lies on A₂.  Sites whose visibility is already
// realized by an existing chord ([C91 §2.1 tex 70]) never qualify.
SiteShot try_site(std::size_t edge_c, Side side, SymbolicY y,
                  const Submap& S, const Polygon& C, std::size_t region,
                  std::size_t A2, const FusedRegionCycle& cycle,
                  const FusedShootContext& fctx) {
    if (region_chord_at(S, region, edge_c, side, y).realized)
        return SiteShot{};

    Point p{edge_x_at_y(C, edge_c, y), y.y, y.tag};
    Side dir = shooting_direction(edge_c, side, C);
    RayHit hit = local_shoot_fused(p, y, dir, cycle, fctx);

    SiteShot out;
    out.hit = hit;
    out.success = hit.hit_arc_idx == A2 &&
        !duplicates_region_chord(S, region, y, edge_c, side,
                                 hit.edge, hit.side);
    return out;
}

// One ∂ᾱ endpoint of an S_α chord, in both frames.
struct AlphaPoint {
    std::size_t edge_a = NONE;  // ᾱ-frame edge
    Side side = LEFT;
    std::size_t edge_c = NONE;  // C-frame edge
    double x = 0.0;
    TourPos pos_a{};            // ∂ᾱ tour position
};

// Everything Lemma 3.2's per-subarc search needs.
struct PieceSearchContext {
    const Submap* S;
    const Polygon* C;
    std::size_t region;
    std::size_t A1;              // region arc of S (table index)
    std::size_t A2;
    const FusedRegionCycle* cycle;
    const FusedShootContext* fctx;
    // The subarc α and its submap S_α ([C91 §3.0(ii) tex 170] cond (3)).
    const ArcPiece* piece;
    const Polygon* Ci;
    std::size_t edge_off;       // Cᵢ edge → C edge
    std::size_t lo;             // Cᵢ edge of ᾱ's edge 0
    Side s;                     // α's ∂C side (single, cond (2))
    // A₂'s exact span in C-frame tour positions.
    TourPos a2_start{};
    TourPos a2_end{};
    // α's endpoints c (cw-first) and d (cw-last) — ∂ᾱ and C frames.
    TourPos c_pos_a{}, d_pos_a{};   // ᾱ-frame ∂ᾱ tour
    TourPos c_pos_c{}, d_pos_c{};   // C-frame ∂C tour
};

// Fill in a VisiblePoint from a successful site shot.
VisiblePoint make_success(const PieceSearchContext& ctx,
                          std::size_t p_edge_c, Side p_side,
                          SymbolicY y, const RayHit& hit) {
    VisiblePoint vp;
    vp.found = true;
    assert(fused_arc_contains(*ctx.S, *ctx.C, ctx.A1, p_edge_c, p_side,
                              y) &&
           "[C91 §3.2]: the site must lie on A₁");
    vp.p_table_arc = ctx.A1;
    vp.p_edge = p_edge_c;
    vp.p_side = p_side;
    vp.p_x = edge_x_at_y(*ctx.C, p_edge_c, y);
    vp.y = y;
    vp.q_table_arc = hit.hit_arc_idx;
    vp.q_edge = hit.edge;
    vp.q_side = hit.side;
    vp.q_x = hit.x;
    return vp;
}

// [C91 §3.2 tex 250–262]: the basic test at an internal node of T, plus
// the branch decision.  Returns found=true on success; otherwise sets
// *next_node to the child to descend into.
VisiblePoint descend_step(const PieceSearchContext& ctx,
                          const TreeDecomposition& td,
                          std::size_t node_idx,
                          std::size_t* next_node) {
    const Submap& Sa = *ctx.piece->submap;
    const Polygon& Ca = *ctx.piece->curve;
    const TDNode& node = td.node(node_idx);
    const Chord& ab = Sa.chord(node.chord_idx);

    // [C91 §2.3]: an h(γ)-granular submap has no null-length chords —
    // a null chord's inner region is empty and degree-1, so contracting
    // it never exceeds the weight bound (granularity criterion (ii)
    // would fail).  The descent therefore never meets one.
    assert(!ab.is_null_length &&
           "[C91 §2.3]: granular S_α cannot contain null-length chords");

    const SymbolicY y_ab = ab.symbolic_y();

    // Chord endpoints in ᾱ and C frames.
    auto alpha_point = [&](bool left) -> AlphaPoint {
        AlphaPoint ap;
        ap.edge_a = left ? ab.left_edge : ab.right_edge;
        ap.side = left ? ab.left_side : ab.right_side;
        ap.edge_c = ctx.edge_off + ctx.lo + ap.edge_a;
        ap.x = edge_x_at_y(Ca, ap.edge_a, y_ab);
        ap.pos_a = tour_pos(Ca, ap.edge_a, ap.side, y_ab);
        return ap;
    };
    AlphaPoint el = alpha_point(true);
    AlphaPoint er = alpha_point(false);
    assert(tour_cmp(el.pos_a, er.pos_a) != 0 &&
           "[C91 §2.1 tex 70]: non-null chord endpoints are distinct ∂ᾱ points");

    // [C91 §2.5 tex 150]: label the endpoints so that a, c, b occur
    // clockwise with respect to D (the outside disk).  ∂C-clockwise is
    // counterclockwise w.r.t. D, so a is the endpoint reached FIRST when
    // walking ∂ᾱ clockwise from c.
    bool left_is_a = at_or_before_from(ctx.c_pos_a, el.pos_a, er.pos_a);
    const AlphaPoint& A = left_is_a ? el : er;
    const AlphaPoint& B = left_is_a ? er : el;
    bool a_is_left_endpoint = left_is_a;

    // Membership in α: α is the FULL side-s of ∂ᾱ (the piece is a
    // vertex-to-vertex subchain, [C91 §3.0(ii)(3) tex 170]).
    bool a_in = (A.side == ctx.s);
    bool b_in = (B.side == ctx.s);
    // [C91 §2.5 tex 150]: α (the cw d→c arc of ∂D) reaching b must pass
    // a first (a, c, b clockwise w.r.t. D).
    assert(!(b_in && !a_in) &&
           "[C91 §2.5]: b ∈ α requires a ∈ α (circle order)");

    // ── Shots from a and b (if in α): success check + a'/b'.
    // [C91 §3.2 tex 260]: "To compute a' and b' ... by local shooting in
    // the region of S incident upon A₁ and A₂ ... no shooting is needed
    // for a or b if the point in question does not lie in α."
    VisiblePoint fail;
    bool a_prime_exists = false, b_prime_exists = false;
    TourPos a_prime_pos{}, b_prime_pos{};
    std::size_t ap_edge = NONE, bp_edge = NONE;
    Side ap_side = LEFT, bp_side = LEFT;
    double ap_x = 0.0, bp_x = 0.0;

    auto endpoint_shot = [&](const AlphaPoint& from, const AlphaPoint& other,
                             bool other_in, bool* prime_exists,
                             TourPos* prime_pos, std::size_t* pe, Side* ps,
                             double* px, VisiblePoint* success) -> bool {
        Side dir = shooting_direction(from.edge_c, from.side, *ctx.C);
        // ab is a chord of V(ᾱ): the unique shooting direction at `from`
        // points along the chord toward the other endpoint.
        assert(((dir == RIGHT) == (other.x > from.x)) &&
               "[C91 §2.1 tex 70]: shooting direction at a chord endpoint "
               "points toward the chord's other endpoint");
        // [C91 §2.1 tex 70]: if `from` coincides with an existing chord
        // endpoint of the region (possible only at α's own endpoints c/d),
        // its visibility IS that chord — the ray retraces it, so the
        // first point of ab ∩ A is the chord's far endpoint (a point of
        // the neighboring arc ⊆ A), read off in O(1) exactly as fusion's
        // case (i) reads t ([C91 §3.1 tex 220]).  No new chord from here.
        {
            RealizedVisibility rv = region_chord_at(
                *ctx.S, ctx.region, from.edge_c, from.side, y_ab);
            if (rv.realized) {
                *prime_exists = true;
                *prime_pos = tour_pos(*ctx.C, rv.other_edge, rv.other_side,
                                      y_ab);
                *pe = rv.other_edge;
                *ps = rv.other_side;
                *px = edge_x_at_y(*ctx.C, rv.other_edge, y_ab);
                return false;
            }
        }
        SiteShot shot = try_site(from.edge_c, from.side, y_ab,
                                 *ctx.S, *ctx.C, ctx.region,
                                 ctx.A2, *ctx.cycle, *ctx.fctx);
        // The other chord endpoint is a ∂C point ON the ray, so a direct
        // crossing exists at d_other and the first hit — on ∂R by
        // [C91 §2.2 Lemma 2.1] — is direct (never through infinity).
        assert(!shot.hit.wrapped &&
               "[C91 §2.2 Lemma 2.1]: shot toward the chord's other "
               "endpoint cannot wrap");
        if (shot.success) {
            *success = make_success(ctx, from.edge_c, from.side, y_ab,
                                    shot.hit);
            return true;
        }
        double d_hit = (dir == LEFT) ? (from.x - shot.hit.x)
                                     : (shot.hit.x - from.x);
        double d_other = (dir == LEFT) ? (from.x - other.x)
                                       : (other.x - from.x);
        assert(d_hit <= d_other &&
               "[C91 §2.2 Lemma 2.1]: the other chord endpoint is on ∂C, "
               "so the first hit cannot lie beyond it");
        if (d_hit < d_other) {
            // Strict interior of ab: the hit is a point of A ([C91 §2.5]
            // — the open chord avoids ∂ᾱ, so any interior ∂C hit is on
            // ∂C \ ∂ᾱ ⊆ A).
            *prime_exists = true;
            *prime_pos = tour_pos(*ctx.C, shot.hit.edge, shot.hit.side, y_ab);
            *pe = shot.hit.edge;
            *ps = shot.hit.side;
            *px = shot.hit.x;
        } else if (!other_in) {
            // Clear segment, other endpoint on A: the first/last point of
            // ab ∩ A is the other endpoint itself.
            *prime_exists = true;
            *prime_pos = tour_pos(*ctx.C, other.edge_c, other.side, y_ab);
            *pe = other.edge_c;
            *ps = other.side;
            *px = other.x;
        }
        // else: clear segment with both endpoints in α ⟹ ab ∩ A = ∅.
        return false;
    };

    if (a_in) {
        VisiblePoint success;
        if (endpoint_shot(A, B, b_in, &a_prime_exists, &a_prime_pos,
                          &ap_edge, &ap_side, &ap_x, &success))
            return success;
    }
    if (b_in) {
        VisiblePoint success;
        if (endpoint_shot(B, A, a_in, &b_prime_exists, &b_prime_pos,
                          &bp_edge, &bp_side, &bp_x, &success))
            return success;
        // Lemma 2.4: a' is the FIRST and b' the LAST point of ab ∩ A —
        // if the set is nonempty both exist.
        assert(a_prime_exists == b_prime_exists &&
               "[C91 §2.5 Lemma 2.4]: a' and b' exist together when both "
               "a, b ∈ B₁ ∪ B₂");
    }

    // ── The two pieces α₁, α₂ of ∂ᾱ between a and b (tex 250).
    // piece₁ := (a → b) clockwise, piece₂ := (b → a) clockwise (open).
    // α ∩ pieceᵢ emptiness: α is the closed side-s tour interval [c, d];
    // (x, y) ∩ [c, d] = ∅ ⟺ [c, d] ⊆ [y, x] (closed complement).
    auto alpha_avoids_open = [&](const TourPos& x, const TourPos& y) {
        return in_closed_cyclic(ctx.c_pos_a, y, x) &&
               in_closed_cyclic(ctx.d_pos_a, y, x) &&
               at_or_before_from(y, ctx.c_pos_a, ctx.d_pos_a);
    };
    bool piece1_empty = alpha_avoids_open(A.pos_a, B.pos_a);
    bool piece2_empty = alpha_avoids_open(B.pos_a, A.pos_a);
    assert(!(piece1_empty && piece2_empty) &&
           "[C91 §3.2]: α has interior points, so some piece meets it");

    // Which ∂ᾱ piece contains c (with the a-tiebreak: c == a ⟹ walking
    // clockwise from a enters α immediately, so c⁺ ∈ piece₁).
    bool c_in_piece1 = in_half_open_cyclic(ctx.c_pos_a, A.pos_a, B.pos_a);

    // ── Decide the rejected piece (tex 250: "find some i ∈ {1,2} such
    // that α ∩ αᵢ is empty or has no point that sees A₂").
    bool reject_piece1;
    if (piece1_empty) {
        reject_piece1 = true;
    } else if (piece2_empty) {
        reject_piece1 = false;
    } else {
        // [C91 §2.5 Lemma 2.4]: both B₁, B₂ nonempty — classify the
        // shielding of A's pieces and locate A₂'s piece.
        //
        // [C91 §2.5]: "the third case, where the traversal from c to d
        // avoids both a and b, was eliminated earlier, since it
        // corresponds to a situation where one of the Bᵢ is empty" — α
        // is connected, so a ∉ α ∧ b ∉ α forces α inside one piece.
        assert(a_in &&
               "[C91 §2.5]: both B₁ and B₂ nonempty requires a ∈ B₁ ∪ B₂");
        bool a_eq_b = a_prime_exists && b_prime_exists &&
                      ap_edge == bp_edge && ap_side == bp_side &&
                      ap_x == bp_x;
        ShieldingResult sh = classify_shielding(a_prime_exists,
                                                b_prime_exists, a_eq_b);

        // A₂'s piece along A (indexed from c per shielding.h).  A runs
        // clockwise along ∂C from d to c (tex 258), so rank positions
        // from d; pieces from c are delimited by a' then b'
        // ([c,a'], [a',b'], [b',d] — Lemma 2.4).  A boundary exactly at
        // A₂'s start resolves upward (A₂ extends clockwise from it).
        auto at_or_before_d = [&](const TourPos& u, const TourPos& v) {
            return at_or_before_from(ctx.d_pos_c, u, v);
        };
        std::size_t piece_idx;
        if (!a_prime_exists) {
            assert(!b_prime_exists);
            piece_idx = 0;                       // single piece [c, d]
        } else if (!b_prime_exists) {
            piece_idx = at_or_before_d(a_prime_pos, ctx.a2_start) ? 0 : 1;
            // consistency: piece 1 = [d, a'] must contain A₂'s end too.
            if (piece_idx == 1)
                assert(at_or_before_d(ctx.a2_end, a_prime_pos) &&
                       "[C91 §2.5 Lemma 2.4]: A₂ must lie in one piece of A");
        } else if (a_eq_b) {
            piece_idx = at_or_before_d(a_prime_pos, ctx.a2_start) ? 0 : 1;
            if (piece_idx == 1)
                assert(at_or_before_d(ctx.a2_end, a_prime_pos));
        } else {
            // [C91 §2.5 Lemma 2.4]: along A from c: a' then b'.
            assert(at_or_before_d(b_prime_pos, a_prime_pos) &&
                   "[C91 §2.5 Lemma 2.4]: piece order along A is c, a', b', d");
            if (at_or_before_d(a_prime_pos, ctx.a2_start)) {
                piece_idx = 0;
            } else if (at_or_before_d(b_prime_pos, ctx.a2_start)) {
                piece_idx = 1;
                assert(at_or_before_d(ctx.a2_end, a_prime_pos));
            } else {
                piece_idx = 2;
                assert(at_or_before_d(ctx.a2_end, b_prime_pos));
            }
        }
        assert(piece_idx < sh.num_pieces &&
               "[C91 §2.5 Lemma 2.4]: A₂'s piece must be one of the 1–3 pieces");

        // B₁ is the piece of ∂D between a and b containing c
        // ([C91 §2.5 tex 150]: H₁ ∋ c; B₁ = closure(α ∩ H₁)).
        BoundarySide rejected_b = sh.shielded_from[piece_idx];
        reject_piece1 = (rejected_b == BoundarySide::B1) == c_in_piece1;
    }

    // ── Branch to the child covering the KEPT piece (tex 252: "we find
    // the child of the root that corresponds to α₂").
    //
    // The arc of S_α STARTING at a (clockwise) lies inside piece₁
    // (a → b); a count-1 slot instead holds the arc ENDING at a, which
    // lies inside piece₂.  Its region identifies the dual-tree component
    // of that piece; tree_decomposition.cpp's convention maps
    // left_child ↔ the component of chord.region[0].
    const Chord::AdjArcs& a_adj = a_is_left_endpoint ? ab.left_adj
                                                     : ab.right_adj;
    std::size_t probe_arc;
    bool probe_in_piece1;
    if (a_adj.count == 2) {
        // [s2-adjacency-convention]: identify the arc STARTING at a by
        // its full start position (edge + side + derived start-y, [C91
        // §2.4 tex 133]) — the y alone is ambiguous whenever the ENDER
        // slot's arc starts at the chord's other endpoint b (a leaf
        // pocket b→a has start-y == y_ab) or at another chord sourced
        // at the same vertex.  Same rule as fusion.cpp's
        // arc_starts_at_chord_slot.
        const std::size_t a_edge =
            a_is_left_endpoint ? ab.left_edge : ab.right_edge;
        const Side a_side = a_is_left_endpoint ? ab.left_side : ab.right_side;
        auto starts_at_a = [&](std::size_t ai) {
            const Arc& arc = Sa.arc(ai);
            return arc.first_edge == a_edge && arc.first_side == a_side &&
                   symbolic_y_equal(Sa.arc_start_symbolic_y(ai, Ca), y_ab);
        };
        bool first_starts = starts_at_a(a_adj.arcs[0]);
        probe_arc = first_starts ? a_adj.arcs[0] : a_adj.arcs[1];
        assert((first_starts || starts_at_a(a_adj.arcs[1])) &&
               "[C91 §2.4(ii)]: one adj slot must hold the arc starting "
               "at the chord endpoint");
        probe_in_piece1 = true;
    } else {
        // [C91 §2.2 tex 94]: vertex endpoint's single slot = before-arc
        // (ends at a) — it lies clockwise-before a, i.e. in piece₂.
        probe_arc = a_adj.arcs[0];
        probe_in_piece1 = false;
    }
    std::size_t probe_region = Sa.arc(probe_arc).region_node;
    assert((probe_region == ab.region[0] || probe_region == ab.region[1]) &&
           "[C91 §2.4(ii)]: adj arc belongs to one of the chord's regions");
    bool probe_child_is_left = (probe_region == ab.region[0]);

    bool keep_piece1 = !reject_piece1;
    bool keep_probe_side = (keep_piece1 == probe_in_piece1);
    bool go_left = (keep_probe_side == probe_child_is_left);
    *next_node = go_left ? node.left_child : node.right_child;
    assert(*next_node != NONE &&
           "[C91 §2.3]: internal TD node has two children");
    return fail;
}

// [C91 §3.2 tex 248–256]: search one subarc α of A₁ — descend T, then
// check the leaf region's vertices belonging to α.
VisiblePoint search_piece(const PieceSearchContext& ctx) {
    const Submap& Sa = *ctx.piece->submap;
    const Polygon& Ca = *ctx.piece->curve;

    const TreeDecomposition& td = Sa.tree_decomposition();
    assert(!td.empty() &&
           "[C91 §3.0(ii)(3) tex 170]: normal-form S_α carries its tree "
           "decomposition");

    std::size_t node_idx = td.root();
    while (td.node(node_idx).is_internal()) {
        std::size_t next = NONE;
        VisiblePoint vp = descend_step(ctx, td, node_idx, &next);
        if (vp.found) return vp;
        node_idx = next;
    }

    // ── Leaf: "we examine each vertex of the region associated with it
    // and, among those belonging to α, we check whether any of them can
    // see A₂" (tex 252).  O(h(γ₂)) vertices, O(f(γ₂)) per shot.
    std::size_t leaf_region = td.node(node_idx).region_idx;
    RegionArcs arcs = collect_region_arcs(Sa, leaf_region);
    for (std::size_t k = 0; k < arcs.count; ++k) {
        const Arc& ra = Sa.arc(arcs.arcs[k]);
        // Enumerate the arc's side-s legs ([C91 §2.4 tex 142]: a wrapped
        // arc covers per-leg edge ranges).
        assert(Sa.start_vertex != NONE && Sa.end_vertex != NONE &&
               "[C91 §2.4(iii)]: S_α's endpoints must be set");
        ArcLeg all_legs[3];
        std::size_t total_legs =
            ra.legs(Sa.start_vertex, Sa.end_vertex, all_legs);
        for (std::size_t g = 0; g < total_legs; ++g) {
            if (all_legs[g].side != ctx.s) continue;
            for (std::size_t ee = all_legs[g].lo; ee <= all_legs[g].hi;
                 ++ee) {
                // Both endpoint vertices of edge ee (ᾱ vertex indices
                // ee and ee+1).
                for (std::size_t vv = ee; vv <= ee + 1; ++vv) {
                    SymbolicY vy = symbolic_y_of(Ca.vertex(vv));
                    // [C91 §2.1 tex 72]: inside-pair duplicates see only
                    // their distance-0 sibling — never candidates.
                    {
                        std::size_t edge_c = ctx.edge_off + ctx.lo + ee;
                        std::size_t vidx_c = symbolic_y_equal(
                            vy, symbolic_y_of(
                                ctx.C->vertex(edge_c))) ? edge_c
                                                        : edge_c + 1;
                        if (is_inside_companion(*ctx.C, edge_c, ctx.s,
                                                vidx_c))
                            continue;
                    }
                    // The site must lie on THIS arc (a mid-edge chord
                    // endpoint may cut the edge short); vertices outside
                    // the arc belong to a neighboring region.
                    TourPos site = tour_pos(Ca, ee, ctx.s, vy);
                    TourPos a_start = tour_pos(
                        Ca, ra.first_edge, ra.first_side,
                        Sa.arc_start_symbolic_y(arcs.arcs[k], Ca));
                    TourPos a_end = tour_pos(
                        Ca, ra.last_edge, ra.last_side,
                        Sa.arc_end_symbolic_y(arcs.arcs[k], Ca));
                    if (!in_closed_cyclic(site, a_start, a_end)) continue;

                    std::size_t edge_c = ctx.edge_off + ctx.lo + ee;
                    SiteShot shot = try_site(edge_c, ctx.s, vy,
                                             *ctx.S, *ctx.C, ctx.region,
                                             ctx.A2, *ctx.cycle, *ctx.fctx);
                    if (shot.success)
                        return make_success(ctx, edge_c, ctx.s, vy,
                                            shot.hit);
                }
            }
        }
    }

    return VisiblePoint{};
}

} // namespace

VisiblePoint find_visible_point(const Submap& S, const Polygon& C,
                                std::size_t region,
                                std::size_t A1, std::size_t A2,
                                const FusedRegionCycle& cycle,
                                const std::vector<ArcProvenance>& provenance,
                                const ConformalityOracles& oracles) {
    assert(S.arc(A1).edge_count > 0 && S.arc(A2).edge_count > 0 &&
           "[C91 §2.1 tex 70/72]: zero-length arcs are single ∂C points "
           "whose visibility is already realized — never candidates");

    const std::size_t n1e = oracles.C1->num_edges();

    FusedShootContext fctx;
    fctx.S = &S;
    fctx.C = &C;
    fctx.C1 = oracles.C1;
    fctx.C2 = oracles.C2;
    fctx.ray1 = oracles.ray1;
    fctx.ray2 = oracles.ray2;
    fctx.provenance = &provenance;

    // A₂'s exact span (C frame) for the shielding test.
    TourPos a2_start = fused_arc_start(S, C, A2);
    TourPos a2_end = fused_arc_end(S, C, A2);

    // Is the ∂C point on A₁?  (Sites are only valid shooting locations
    // when they lie on the region's boundary.)
    auto on_A1 = [&](std::size_t edge_c, Side side, SymbolicY y) -> bool {
        return fused_arc_contains(S, C, A1, edge_c, side, y);
    };

    // [C91 §3.2 tex 244/248]: "we invoke the arc-cutter associated with
    // the arc of S₁ or S₂ containing A₁."  A₁ is one arc-structure —
    // wrap-spanning arcs included ([C91 §2.4 tex 142]) — inside one Sᵢ
    // arc, so a single cut covers it; the cutter subdivides any
    // double-backing ([C91 §3.0(ii)(2) tex 170]).
    {
        const ArcProvenance& pr = provenance[A1];

        const std::size_t off = pr.on_c1 ? 0 : n1e;
        const Polygon& Ci = pr.on_c1 ? *oracles.C1 : *oracles.C2;
        const ArcCuttingOracle& cutter = pr.on_c1 ? *oracles.cut1
                                                  : *oracles.cut2;
        const std::size_t g_bound = pr.on_c1 ? oracles.g_gamma1
                                             : oracles.g_gamma2;
        const std::size_t h_bound = pr.on_c1 ? oracles.h_gamma1
                                             : oracles.h_gamma2;

        const Arc& a1_struct = S.arc(A1);
        Subarc target;
        target.first_edge = a1_struct.first_edge - off;
        target.first_side = a1_struct.first_side;
        target.last_edge = a1_struct.last_edge - off;
        target.last_side = a1_struct.last_side;
        // [C91 §3.0(i)(ii) tex 169/170]: α' is "specified by its two
        // endpoints" — exact ys are frame-independent; only the edges
        // carry the −off shift into Cᵢ's frame.
        target.first_y = S.arc_start_symbolic_y(A1, C);
        target.last_y = S.arc_end_symbolic_y(A1, C);
        assert_subarc_clockwise(target);

        std::vector<ArcPiece> pieces = cutter.cut(pr.arc_in_si, target);
        // [C91 §3.0(ii) tex 170]: mandatory post-condition validation.
        assert_cut_postconditions(Ci, target, pieces.data(), pieces.size(),
                                  g_bound, h_bound);

        // [C91 §3.2 tex 248]: "Except for at most two single-edge
        // subarcs attached to the endpoints of A₁ ... We search each
        // subarc in turn, stopping as soon as we find a good vertex or
        // point, if ever."  Boundary pieces carry no submap, so they are
        // exempt from the tree search — but their O(1) vertices lying on
        // A₁ are still candidates and are checked directly in O(f(γ₂))
        // each ([C91 §3.2 tex 246]).  When A₁ spans at most two partial
        // edges the cut yields ONLY boundary pieces, and Lemma 3.3's
        // guaranteed vertex must remain reachable through them.
        for (const ArcPiece& piece : pieces) {
            if (piece.is_boundary_piece) {
                assert(piece.subarc.first_edge == piece.subarc.last_edge);
                const std::size_t e = piece.subarc.first_edge;
                const Side s = piece.subarc.first_side;
                for (std::size_t vv = e; vv <= e + 1; ++vv) {
                    SymbolicY vy = symbolic_y_of(Ci.vertex(vv));
                    if (!on_A1(off + e, s, vy)) continue;
                    // [C91 §2.1 tex 72]: inside-pair duplicates see only
                    // their distance-0 sibling — never candidates.
                    {
                        std::size_t vidx_c = symbolic_y_equal(
                            vy, symbolic_y_of(C.vertex(off + e)))
                                ? off + e : off + e + 1;
                        if (is_inside_companion(C, off + e, s, vidx_c))
                            continue;
                    }
                    SiteShot shot = try_site(off + e, s, vy, S, C, region,
                                             A2, cycle, fctx);
                    if (shot.success) {
                        VisiblePoint vp;
                        vp.found = true;
                        vp.p_table_arc = A1;   // on_A1 verified above
                        vp.p_edge = off + e;
                        vp.p_side = s;
                        vp.p_x = edge_x_at_y(C, off + e, vy);
                        vp.y = vy;
                        vp.q_table_arc = shot.hit.hit_arc_idx;
                        vp.q_edge = shot.hit.edge;
                        vp.q_side = shot.hit.side;
                        vp.q_x = shot.hit.x;
                        return vp;
                    }
                }
                continue;
            }

            const std::size_t lo = std::min(piece.subarc.first_edge,
                                            piece.subarc.last_edge);
            [[maybe_unused]] const std::size_t hi =
                std::max(piece.subarc.first_edge, piece.subarc.last_edge);
            assert(piece.subarc.first_side == piece.subarc.last_side &&
                   "[C91 §3.0(ii)(2) tex 170]: pieces do not double-back");
            const Side s = piece.subarc.first_side;

            // Pin the ᾱ ↔ Cᵢ frame correspondence: ᾱ's vertex 0 is Cᵢ's
            // vertex lo (SoS tags are global P-frame indices).
            assert(piece.curve->vertex(0).index == Ci.vertex(lo).index &&
                   "[C91 §3.0(ii)(3) tex 170]: ᾱ must be the ascending "
                   "vertex-to-vertex subchain of Cᵢ over the piece's range");
            assert(piece.curve->num_edges() == hi - lo + 1);

            PieceSearchContext pctx;
            pctx.S = &S;
            pctx.C = &C;
            pctx.region = region;
            pctx.A1 = A1;
            pctx.A2 = A2;
            pctx.cycle = &cycle;
            pctx.fctx = &fctx;
            pctx.piece = &piece;
            pctx.Ci = &Ci;
            pctx.edge_off = off;
            pctx.lo = lo;
            pctx.s = s;
            pctx.a2_start = a2_start;
            pctx.a2_end = a2_end;

            // α's endpoints: c = clockwise-first, d = clockwise-last.  In
            // ᾱ's frame α is the full side-s: LEFT runs vertex 0 → n̄,
            // RIGHT runs vertex n̄ → 0 (clockwise tour).
            const Polygon& Ca = *piece.curve;
            const std::size_t na = Ca.num_edges();
            if (s == LEFT) {
                pctx.c_pos_a = tour_pos(Ca, 0, LEFT,
                                        symbolic_y_of(Ca.vertex(0)));
                pctx.d_pos_a = tour_pos(Ca, na - 1, LEFT,
                                        symbolic_y_of(Ca.vertex(na)));
                pctx.c_pos_c = tour_pos(C, off + lo, LEFT,
                                        symbolic_y_of(Ca.vertex(0)));
                pctx.d_pos_c = tour_pos(C, off + lo + na - 1, LEFT,
                                        symbolic_y_of(Ca.vertex(na)));
            } else {
                pctx.c_pos_a = tour_pos(Ca, na - 1, RIGHT,
                                        symbolic_y_of(Ca.vertex(na)));
                pctx.d_pos_a = tour_pos(Ca, 0, RIGHT,
                                        symbolic_y_of(Ca.vertex(0)));
                pctx.c_pos_c = tour_pos(C, off + lo + na - 1, RIGHT,
                                        symbolic_y_of(Ca.vertex(na)));
                pctx.d_pos_c = tour_pos(C, off + lo, RIGHT,
                                        symbolic_y_of(Ca.vertex(0)));
            }

            VisiblePoint vp = search_piece(pctx);
            if (vp.found) return vp;
        }
    }

    return VisiblePoint{};
}

// ════════════════════════════════════════════════════════════════
//  restore_conformality — [C91 §3.2 Lemma 3.4]
// ════════════════════════════════════════════════════════════════

void restore_conformality(Submap& S, const Polygon& C,
                          const ConformalityOracles& oracles) {
    assert(oracles.S1 && oracles.S2 && oracles.C1 && oracles.C2 &&
           oracles.ray1 && oracles.ray2 && oracles.cut1 && oracles.cut2 &&
           "[C91 §3.0 tex 166–170]: §3.2 requires all per-submap oracles");

    // [C91 §3.2 tex 244]: provenance, computed once (before mutations).
    std::vector<ArcProvenance> provenance = compute_arc_provenance(
        S, C, *oracles.S1, *oracles.C1, *oracles.S2, *oracles.C2);

    // Region → live arcs inventory (maintained across insertions so
    // cycle collection stays O(1) per bounded region — Lemma 3.4).
    std::vector<std::vector<std::size_t>> region_arcs(S.num_nodes());
    for (std::size_t ai = 0; ai < S.num_arcs(); ++ai) {
        const Arc& a = S.arc(ai);
        if (a.dead) continue;
        region_arcs[a.region_node].push_back(ai);
    }

    std::vector<std::size_t> work;
    for (std::size_t r = 0; r < S.num_nodes(); ++r)
        if (!S.node(r).dead && !region_arcs[r].empty())
            work.push_back(r);

    while (!work.empty()) {
        std::size_t r = work.back();
        work.pop_back();
        assert(!S.node(r).dead && "[C91 §3.2]: regions are never removed here");

        FusedRegionCycle cycle = fused_region_cycle(S, C, r, region_arcs[r]);
        if (cycle.count <= 4) continue;

        // [C91 §3.2 tex 264]: "apply Lemma 3.2 to every pair of
        // nonconsecutive arcs until we find a chord."
        bool found = false;
        const std::size_t k = cycle.count;
        for (std::size_t i = 0; i < k && !found; ++i) {
            if (cycle.arcs[i].is_zero_length) continue;
            for (std::size_t j = 0; j < k && !found; ++j) {
                // (i) of Lemma 3.3: i − j ≢ −1, 0, 1 (mod k).
                if (j == i || j == (i + 1) % k || j == (i + k - 1) % k)
                    continue;
                if (cycle.arcs[j].is_zero_length) continue;

                VisiblePoint vp = find_visible_point(
                    S, C, r, cycle.arcs[i].arc, cycle.arcs[j].arc,
                    cycle, provenance, oracles);
                if (!vp.found) continue;

                // The cycle's table arcs, in boundary order, for
                // insert_chord.
                std::vector<std::size_t> flat;
                for (std::size_t li = 0; li < k; ++li)
                    flat.push_back(cycle.arcs[li].arc);

                Submap::ChordPointSpec p{vp.p_table_arc, vp.p_edge,
                                         vp.p_side, vp.p_x};
                Submap::ChordPointSpec q{vp.q_table_arc, vp.q_edge,
                                         vp.q_side, vp.q_x};
                auto res = S.insert_chord(p, q, vp.y, r,
                                          flat.data(), flat.size(), C);

                // The appended halves inherit their source's provenance
                // ([C91 §3.2 tex 264]: "region arcs can only become
                // smaller" — Lemma 3.2 stays applicable).
                provenance.resize(S.num_arcs());
                provenance[res.p_after_arc] = provenance[vp.p_table_arc];
                provenance[res.q_after_arc] = provenance[vp.q_table_arc];

                // Refresh the two regions' inventories.
                region_arcs.resize(S.num_nodes());
                region_arcs[r].clear();
                region_arcs[res.new_region].clear();
                flat.push_back(res.p_after_arc);
                flat.push_back(res.q_after_arc);
                for (std::size_t ai : flat)
                    region_arcs[S.arc(ai).region_node].push_back(ai);

                work.push_back(r);
                work.push_back(res.new_region);
                found = true;
            }
        }
        // [C91 §3.2 Lemma 3.3]: k > 4 guarantees a nonconsecutive pair
        // with a vertex of ∂C on one seeing the other; Lemma 3.2 finds
        // it.  (Zero-length arcs cannot carry the guarantee — their
        // point's visibility is its distance-0 companion, [C91 §2.1
        // tex 72] — so the guaranteed pair is among the non-null arcs.)
        assert(found &&
               "[C91 §3.2 Lemma 3.3]: a region with more than four arcs "
               "must admit a chord between nonconsecutive arcs");
    }

    // [C91 §3.2 tex 264]: "we iterate on this process until no region
    // has more than four arcs" — postcondition.  A region's boundary
    // alternates arcs and chords ([C91 §2.2 tex 96]), so ≤ 4
    // arc-structures ⟺ degree ≤ 4 = conformality ([C91 §2.3 tex 114]).
#ifndef NDEBUG
    for (std::size_t r = 0; r < S.num_nodes(); ++r) {
        if (S.node(r).dead || region_arcs.size() <= r) continue;
        if (region_arcs[r].empty()) continue;
        FusedRegionCycle cycle = fused_region_cycle(S, C, r, region_arcs[r]);
        assert(cycle.count <= 4 &&
               "[C91 §3.2 tex 264]: every region must end with ≤ 4 arcs");
    }
#endif
    assert(S.is_conformal() &&
           "[C91 §2.3 tex 114]: S must be conformal after §3.2");
    S.assert_tree_property();
}

} // namespace chazelle
