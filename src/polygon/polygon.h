#pragma once

// [C91 §2.1]: a simple (nonclosed) polygonal curve C with vertices
// v1, ..., vn and its directed edge table.
//
// [C91 §2.4 tex 133]: "Let P be the input polygonal curve (the one
// whose visibility map is sought) and let C be the subchain of P whose
// visibility map (or submap) we wish to represent. ... We assume that
// the edges of P are stored in a table (the input table) in the order
// in which they occur along the boundary of P. ... The input table is
// read-only: it is never to be modified or even copied."  A Polygon is
// therefore a VIEW (offset + length) into one shared immutable input
// table; constructing the merged curve C = C₁ ∪ C₂ ([C91 §3 tex 160])
// is O(1) index arithmetic, never a copy — Lemma 4.1's sublinear
// per-merge budget ([C91 §4.1 tex 347–350]) forbids anything else.

#include "point.h"
#include "edge.h"
#include "perturbation.h"
#include "../common.h"

#include <cassert>
#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

namespace chazelle {

class Polygon {
public:
    // Root curve: allocates the input table (done once per input
    // setup; [C91 §2.4 tex 133]).  O(n).
    explicit Polygon(std::vector<Point> vertices);

    // [C91 §3 tex 160]: the merged curve C = C₁ ∪ C₂, where "C₁ ∩ C₂
    // is a vertex of P" — C₁ and C₂ must be contiguous views into the
    // SAME input table sharing exactly that junction vertex.  O(1):
    // the table is never copied ([C91 §2.4 tex 133]) and the y-extreme
    // vertices combine from the two halves.
    Polygon(const Polygon& c1, const Polygon& c2);

    // Contiguous subchain view [first, first + count) of this curve
    // ([C91 §2.4 tex 133]: every curve is a subchain of P).  O(count)
    // for the y-extreme scan — used at input setup, not inside merges.
    Polygon subchain(std::size_t first, std::size_t count) const;

    std::size_t num_vertices() const noexcept { return len_; }
    std::size_t num_edges()    const noexcept { return len_ - 1; }

    const Point& vertex(std::size_t i) const noexcept {
        assert(i < len_);
        return table_->vertices[offset_ + i];
    }

    // Directed edge i runs from local vertex i to i + 1 ([C91 §2.1]).
    Edge edge(std::size_t i) const noexcept {
        assert(i + 1 < len_);
        return Edge{i, i + 1};
    }

    // [C91 §2.1 Fig 2.2.3]: one of C's two endpoints.
    bool is_endpoint(std::size_t vertex_index) const noexcept {
        assert(vertex_index < len_ && "[C91 §2.1]: invalid vertex");
        return vertex_index == 0 || vertex_index == len_ - 1;
    }

    // [C91 §2.1 Fig 2.2.2]: local y-extremum under SoS.  Endpoints are
    // case 3, not case 2 — never extrema.
    bool is_y_extremum(std::size_t vertex_index) const noexcept {
        assert(vertex_index < len_ && "[C91 §2.1]: invalid vertex");
        if (is_endpoint(vertex_index)) return false;
        return is_local_y_extremum(vertex(vertex_index - 1),
                                   vertex(vertex_index),
                                   vertex(vertex_index + 1));
    }

    // [C91 §2.2]: count of nonnull-length edges in [lo, hi].  O(1) via
    // the table-wide prefix sums; used by region_weight ("max nonnull
    // edges in any of its arcs").
    std::size_t count_nonnull_edges(std::size_t lo,
                                    std::size_t hi) const noexcept {
        assert(lo <= hi && "[C91 §2.2]: lo must not exceed hi");
        assert(hi < num_edges() && "[C91 §2.2]: hi must be within bounds");
        return table_->nonnull_prefix[offset_ + hi + 1] -
               table_->nonnull_prefix[offset_ + lo];
    }

    // [C91 §2.2 tex 106]: edge of null (zero) length — padding edges of
    // [C91 §4 tex 316] are the only source (input curves are simple).
    bool edge_is_null(std::size_t i) const noexcept {
        return count_nonnull_edges(i, i) == 0;
    }

    // Smallest nonnull edge ≥ i (local frame), or num_edges() if none.
    // O(1) via the table-wide successor array.  [C91 §2.2 tex 106]
    // makes every weight/size bound of the paper count NONNULL edges
    // only, so scans whose budgets are stated in region weights ([C91
    // §3.4 tex 306]: "checking all the O(γ) edges of the region") must
    // step over null runs in O(1), not edge by edge.
    std::size_t next_nonnull_edge(std::size_t i) const noexcept {
        assert(i < num_edges());
        std::size_t t = table_->next_nonnull[offset_ + i];
        return (t >= offset_ + num_edges()) ? num_edges() : t - offset_;
    }

    // [C91 §2.4 tex 133]: position of this view's first vertex in the
    // one input table of P — the paper's "pointers to the input table"
    // ([C91 §3.0 tex 169]) are global table positions; views translate
    // with this offset.  The up-phase ([C91 §4.1 tex 336]) needs it to
    // align curves with the chain grid of [C91 §4 tex 316].
    std::size_t table_offset() const noexcept { return offset_; }

    // [C91 §3.4 tex 306]: the ray-shooting structure's vertical line is
    // anchored "to the right of all the vertices of P", and its polar
    // segments are delimited by the curve's global y-extremes; the
    // extreme vertices (under the SoS order of [C91 §2 tex 47]) are
    // scanned once at input setup and combined in O(1) per merge.
    std::size_t max_y_vertex() const noexcept { return max_y_vertex_; }
    std::size_t min_y_vertex() const noexcept { return min_y_vertex_; }

private:
    // The one input table of P: immutable, shared by every subchain
    // view, never copied ([C91 §2.4 tex 133]).
    struct Table {
        std::vector<Point> vertices;
        // nonnull_prefix[i] = number of nonnull edges among the table's
        // edges [0, i).
        std::vector<std::size_t> nonnull_prefix;
        // next_nonnull[i] = smallest table edge ≥ i of nonzero length
        // (= #edges if none) — O(1) null-run skipping for scans whose
        // budgets count nonnull edges only ([C91 §2.2 tex 106]).
        std::vector<std::size_t> next_nonnull;
    };

    Polygon() = default;                       // internal (subchain)

    void find_y_extremes();                    // O(len_) scan

    std::shared_ptr<const Table> table_;
    std::size_t offset_ = 0;                   // first vertex in table
    std::size_t len_ = 0;                      // number of vertices
    std::size_t max_y_vertex_ = 0;             // local indices
    std::size_t min_y_vertex_ = 0;
};

// [C91 §2 tex 47]: SoS-aware interpolation parameter on an edge.
//
// Returns t ∈ [0, 1] such that the (perturbed) horizontal line at
// `target_y` crosses edge `edge_idx` at the convex combination
// (1-t)·vs + t·ve.  When `target_y`'s SymbolicY matches an endpoint's
// (same raw y and tag), the perturbed crossing IS that endpoint —
// returns 0.0 or 1.0 exactly.  Otherwise raw y's strictly bracket
// target_y (else SoS would have forced a tag-match), and standard
// linear interpolation is exact.
//
// Required by the paper's "distinct y" relaxation under SoS: with
// repeated raw y-coords a horizontal edge has `vs.y == ve.y` and
// only the symbolic short-circuit produces a defined answer.
inline double edge_t_at_y(const Polygon& C, std::size_t edge_idx,
                          SymbolicY target_y) {
    const auto e = C.edge(edge_idx);
    const Point& vs = C.vertex(e.start_idx);
    const Point& ve = C.vertex(e.end_idx);
    if (symbolic_y_equal(target_y, symbolic_y_of(vs))) return 0.0;
    if (symbolic_y_equal(target_y, symbolic_y_of(ve))) return 1.0;
    // [C91 §2 tex 47]: with no tag-match, the perturbed horizontal line
    // crosses the OPEN edge, i.e. target_y lies symbolically strictly
    // between the endpoint ys.  (A raw-y tie with an endpoint is fine —
    // the crossing is infinitesimally beside that vertex and the linear
    // interpolation below returns the exact limit point.)
    assert(((symbolic_y_less(symbolic_y_of(vs), target_y) &&
             symbolic_y_less(target_y, symbolic_y_of(ve))) ||
            (symbolic_y_less(symbolic_y_of(ve), target_y) &&
             symbolic_y_less(target_y, symbolic_y_of(vs)))) &&
           "[C91 §2 tex 47]: query must lie strictly inside the edge's "
           "perturbed y-range (no crossing otherwise)");
    // Raw-horizontal edge: its endpoints are consecutive vertices of C
    // carrying consecutive SoS indices, so no third tag can fall
    // symbolically strictly between them — unreachable given the assert
    // above.
    assert(vs.y != ve.y &&
           "[C91 §2 tex 47]: horizontal edge requires SoS tag-match at "
           "one endpoint; strictly-between target is unreachable");
    return (target_y.y - vs.y) / (ve.y - vs.y);
}

// [C91 §2 tex 47]: SoS-aware x of where the perturbed horizontal line
// at `target_y` crosses edge `edge_idx`.  Companion to `edge_t_at_y`.
inline double edge_x_at_y(const Polygon& C, std::size_t edge_idx,
                          SymbolicY target_y) {
    const auto e = C.edge(edge_idx);
    const Point& vs = C.vertex(e.start_idx);
    const Point& ve = C.vertex(e.end_idx);
    double t = edge_t_at_y(C, edge_idx, target_y);
    return vs.x + t * (ve.x - vs.x);
}

// [C91 §2 tex 47]: x of the crossing of the perturbed horizontal line
// at `sy` with edge e, if any — the non-asserting companion of
// edge_t_at_y (endpoint tag matches short-circuit to the vertex; a
// strict symbolic containment interpolates).
inline bool edge_crossing_x(const Polygon& C, std::size_t e, SymbolicY sy,
                            double* x) {
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
    // [C91 §2 tex 47]: a raw-horizontal edge cannot reach here — its
    // endpoints carry consecutive SoS tags, so no tag lies symbolically
    // strictly between them (same lemma as edge_t_at_y).
    assert(vs.y != ve.y &&
           "[C91 §2 tex 47]: horizontal edge admits no strictly-interior "
           "crossing");
    double t = (sy.y - vs.y) / (ve.y - vs.y);
    *x = vs.x + t * (ve.x - vs.x);
    return true;
}

// ── Perturbed x-offsets at raw-coincident points ([C91 §2 tex 47]) ─
//
// Under the SoS perturbation, points whose RAW positions coincide
// separate: a ray-level crossing of an edge moves off the shared
// point P by
//     perturbed x  =  x_P + (dx/dy)_edge · δ · dsign,
// where δ = |perturbed(sy) − perturbed(y_P)| is one shared
// infinitesimal (both levels are the same pair of perturbed values
// for every participant at P — their edges all END at P's vertex by
// simplicity, [C91 §2.1 tex 68]) and dsign is the symbolic side of
// sy relative to the vertex.  `perturbed_x_offset` returns the
// coefficient (dx/dy)·dsign — the amount and direction the point
// REALLY moves, per unit δ — which is the [10] first-nonvanishing-
// term evaluation of any x-comparison between such points (the same
// derivation as ray_contact_precedes' configuration (3)).  A wall
// sitting exactly AT the vertex (tag-matched level) does not move:
// offset 0.  A vertical edge's crossing stays at x_P: offset 0.  A
// null edge is tag-matched by construction (its foreign-level span
// is empty, [C91 §2 tex 47]).
inline double perturbed_x_offset(const Polygon& C, SymbolicY sy,
                                  std::size_t edge) {
    const auto& ed = C.edge(edge);
    const Point& a = C.vertex(ed.start_idx);
    const Point& b = C.vertex(ed.end_idx);
    if (symbolic_y_equal(sy, symbolic_y_of(a)) ||
        symbolic_y_equal(sy, symbolic_y_of(b)))
        return 0.0;                       // wall AT the cluster vertex
    const bool a_tied = (a.y == sy.y);
    const bool b_tied = (b.y == sy.y);
    // Neither endpoint raw-tied: the point is a MID-EDGE crossing.
    // Two distinct edges cannot share a mid-interior point
    // ([C91 §2.1 tex 68] simplicity), so a raw coincidence involving
    // it is the edge's own two walls at one crossing — neither moves
    // relative to the other (offset 0), and equal offsets keep the
    // wrapped convention.
    if (!a_tied && !b_tied) return 0.0;
    assert(a_tied != b_tied &&
           "[C91 §2 tex 47]: raw-horizontal edges admit no strict "
           "crossing (consecutive endpoint tags)");
    const Point& at = a_tied ? a : b;
    const Point& aw = a_tied ? b : a;
    if (aw.x == at.x) return 0.0;         // vertical: stays at x_P
    const double slope = (aw.x - at.x) / (aw.y - at.y);
    return symbolic_y_less(sy, symbolic_y_of(at)) ? -slope : slope;
}

// Sentinel for "the source's perturbed x-offset was not supplied"
// (real offsets are finite).
inline constexpr double SOURCE_OFFSET_NONE =
    std::numeric_limits<double>::infinity();
inline bool has_source_offset(double off) noexcept {
    return off != SOURCE_OFFSET_NONE;
}

// [C91 §2 tex 47]: does a contact on `hit_edge`, raw-coincident with
// the ray's source (whose own perturbed x-offset is
// `source_x_offset`, computed ONCE by the original caller in the
// source's own frame — the offset is plain geometry, identical in
// every subchain frame), lie STRICTLY FORWARD of the source along
// `dir` in the perturbed picture?  Equal offsets (the source's own
// companion wall, a collinear edge, or a null stack) keep the
// wrapped convention (met after a full wrap) — the historic behavior
// for a source's own duplicates.
inline bool perturbed_hit_forward(const Polygon& C, SymbolicY sy,
                                  Side dir, double source_x_offset,
                                  std::size_t hit_edge) {
    const double hit_off = perturbed_x_offset(C, sy, hit_edge);
    if (source_x_offset == hit_off) return false;
    return (dir == RIGHT) ? (hit_off > source_x_offset)
                          : (hit_off < source_x_offset);
}

// [C91 §2 tex 47]: at a y-extremum vertex v (branches to u = prev and
// w = next), is the PREV branch infinitesimally LEFT of the NEXT branch
// on the wedge side?  Each branch's x-offset per unit of |Δy| is the
// exact rational (q.x − v.x) / |q.y − v.y| with nonnegative
// denominator, compared cross-multiplied; a null-length branch
// ([C91 §4 tex 316] padding: same coordinates, SoS-separated) is
// vertical — offset 0/1 — and a raw-horizontal branch has offset
// ±1/0.  Equal offsets would overlap the branches, contradicting the
// simplicity of P ([C91 §2.1 tex 68]).
inline bool extremum_prev_branch_left(const Point& u, const Point& v,
                                      const Point& w) {
    assert(is_local_y_extremum(u, v, w) &&
           "[C91 §2.1 tex 72]: wedge sides exist at y-extrema only");
    double nu = u.x - v.x, du = u.y > v.y ? u.y - v.y : v.y - u.y;
    double nw = w.x - v.x, dw = w.y > v.y ? w.y - v.y : v.y - w.y;
    if (nu == 0.0 && du == 0.0) du = 1.0;      // null branch: vertical
    if (nw == 0.0 && dw == 0.0) dw = 1.0;
    if (nu != 0.0 && du == 0.0) { nu = (nu > 0.0) ? 1.0 : -1.0; }
    if (nw != 0.0 && dw == 0.0) { nw = (nw > 0.0) ? 1.0 : -1.0; }
    const double lhs = nu * dw;
    const double rhs = nw * du;
    assert(lhs != rhs &&
           "[C91 §2 tex 47]: extremum branches have distinct x-offsets "
           "(equal offsets ⟹ overlapping edges, non-simple P)");
    return lhs < rhs;
}

// [C91 §2.1 tex 72]: is the ∂C point (edge, side) at vertex `vidx` one of
// a y-extremum's INSIDE-pair duplicates?  "one of these pairs, the one on
// the 'inside' of the turn, gives rise to a chord of null length" — such
// a point sees only its distance-0 sibling and can never source a chord
// to another arc, so Lemma 3.2 never shoots from it, and a ray tie at
// the vertex resolves to the OUTSIDE wall (the outer wall shields the
// wedge).  The inside face computation mirrors rebuild_submap's junction
// null-chord synthesis ([C91 §3.1 tex 224], fusion.cpp).
inline bool is_inside_companion(const Polygon& C, std::size_t edge,
                                Side side, std::size_t vidx) {
    if (vidx == 0 || vidx + 1 >= C.num_vertices()) return false;  // case 3
    const Point& u = C.vertex(vidx - 1);
    const Point& v = C.vertex(vidx);
    const Point& w = C.vertex(vidx + 1);
    if (!is_local_y_extremum(u, v, w)) return false;              // case 1

    // Branch x-order infinitesimally off v ([C91 §2 tex 47]).
    const bool prev_left = extremum_prev_branch_left(u, v, w);

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

// [C91 §2 tex 47]: tie-break between two ray contacts at the SAME raw
// travel position (equal (wrapped, distance) in the metric of
// [C91 §2.1 tex 70], same ray).  Simplicity of P admits exactly two
// configurations:
//
//  (1) Two walls at one tag-matched vertex of C — the ray's symbolic y
//      IS the vertex's, and both contact edges are incident to it.  At
//      a y-extremum exactly one wall is the inside-of-turn duplicate
//      ([C91 §2.1 tex 72]), which sees only its distance-0 sibling —
//      the OUTSIDE wall is struck.  At a non-extremum both walls name
//      the same ∂C companion under its two edge labels ([C91 §3.0(i)
//      tex 169]: either edge "contains" the point); keep the old.
//
//  (2) The padding cluster ([C91 §4 tex 316]): a strict (non-endpoint)
//      crossing of the single real edge arriving at the cluster point
//      ties in RAW x with a vertex contact on the null run.  Under the
//      y-only perturbation ([C91 §2 tex 47]; x stays exact) the strict
//      crossing sits at x₀ − (dx/dy)·δ, δ > 0 the infinitesimal level
//      distance to the tied endpoint, so its side of the cluster
//      follows from the edge's slope sign.  A VERTICAL real edge
//      (dx = 0) keeps x₀ exactly, and the null run hangs inside the
//      double boundary's sliver at its end ([C91 §2.1 tex 70]: the
//      thin-polygon thickness is the smallest infinitesimal) — the
//      real wall shields the run from both directions.
//
// Precondition on callers: no ray is ever shot from an inside-pair
// duplicate toward its own wedge — such a point sees only its
// distance-0 sibling ([C91 §2.1 tex 72]), its visibility IS the null
// chord and is read off, never shot (fusion's case (i) read-off,
// §3.2's candidate-vertex filter, the naive V(C) construction's
// null-chord
// synthesis).  For every other source the wedge is shielded by its
// outside walls, which is what configuration (1) resolves to.
//
// Returns true iff the NEW contact strictly precedes the OLD one along
// the ray; false keeps the old (including equivalent-label ties).
inline bool ray_contact_precedes(const Polygon& C, SymbolicY sy, Side dir,
                                 std::size_t e_new, Side s_new,
                                 std::size_t e_old, Side s_old) {
    if (e_new == e_old && s_new == s_old) return false;   // same wall
    auto tag_matched_vertex = [&](std::size_t e) -> std::size_t {
        const auto& ed = C.edge(e);
        if (symbolic_y_equal(sy, symbolic_y_of(C.vertex(ed.start_idx))))
            return ed.start_idx;
        if (symbolic_y_equal(sy, symbolic_y_of(C.vertex(ed.end_idx))))
            return ed.end_idx;
        return NONE;
    };
    const std::size_t vn = tag_matched_vertex(e_new);
    const std::size_t vo = tag_matched_vertex(e_old);
    if (vn != NONE && vo != NONE) {
        // (1): both contacts at one vertex ([C91 §2 tex 47]: one vertex
        // per SoS tag).
        assert(vn == vo &&
               "[C91 §2 tex 47]: a symbolic y names a unique vertex");
        return is_inside_companion(C, e_old, s_old, vo) &&
               !is_inside_companion(C, e_new, s_new, vn);
    }
    if (vn == NONE && vo == NONE) {
        // (3): TWO strict crossings at one raw x — they flank a vertex
        // whose RAW y equals the ray's under a different SoS tag (its
        // two incident edges both cross the perturbed level beside it;
        // a coincidence on non-adjacent edges would put a raw point of
        // one edge on another, contradicting simplicity, [C91 §2.1
        // tex 68]).  The perturbed crossings sit at x_v + (dx/dy)ᵢ·δ
        // with ONE shared infinitesimal δ = perturbed(sy) −
        // perturbed(y_v) ([C91 §2 tex 47]), so the ray order follows
        // the slope order — flipped when sy lies symbolically below
        // the vertex, and by travel direction.  Equal slopes are
        // collinear edges through the vertex: one geometric wall under
        // two labels; keep the old.
        const std::size_t vidx = std::max(e_new, e_old);
        assert(vidx == std::min(e_new, e_old) + 1 &&
               "[C91 §2.1 tex 68]: coincident strict crossings flank "
               "one shared vertex");
        const Point& vb = C.vertex(vidx);
        assert(vb.y == sy.y && vb.index != sy.tag &&
               "[C91 §2 tex 47]: the flanked vertex shares the ray's "
               "raw level under a different tag");
        auto dxdy = [&](std::size_t e) -> double {
            const auto& ed = C.edge(e);
            const Point& a = C.vertex(ed.start_idx);
            const Point& b = C.vertex(ed.end_idx);
            assert(a.y != b.y &&
                   "[C91 §2 tex 47]: raw-horizontal edges admit no "
                   "strict crossing");
            return (b.x - a.x) / (b.y - a.y);
        };
        const double dsign =
            symbolic_y_less(sy, symbolic_y_of(vb)) ? -1.0 : 1.0;
        const double cn = dxdy(e_new) * dsign;
        const double co = dxdy(e_old) * dsign;
        if (cn == co) return false;             // collinear: same wall
        return (dir == RIGHT) ? (cn < co) : (cn > co);
    }
    // (2): exactly one strict crossing.
    const bool new_is_strict = (vn == NONE);
    const std::size_t es = new_is_strict ? e_new : e_old;
    const auto& ed = C.edge(es);
    const Point& vs = C.vertex(ed.start_idx);
    const Point& ve = C.vertex(ed.end_idx);
    // The strict edge's tied endpoint sits at the coincidence point (an
    // interior coincidence would again violate simplicity); a raw-
    // horizontal edge admits no strict crossing (edge_t_at_y).
    assert((vs.y == sy.y) != (ve.y == sy.y) &&
           "[C91 §2 tex 47]: exactly one endpoint of the strict edge "
           "lies at the tied raw level");
    const Point& at = (vs.y == sy.y) ? vs : ve;
    const Point& away = (vs.y == sy.y) ? ve : vs;
    assert(away.y != at.y);
    // Perturbed crossing x = at.x + (dx/dy)·(sy_level − at_level):
    // sign(correction) = sign(dx/dy) flipped when sy lies symbolically
    // below the tied endpoint.  Zero correction (vertical edge): the
    // strict wall shields the null run — strict wins both directions.
    const bool below = symbolic_y_less(sy, symbolic_y_of(at));
    int slope_sign = (away.x == at.x) ? 0
                    : (((away.x > at.x) == (away.y > at.y)) ? 1 : -1);
    int corr = below ? -slope_sign : slope_sign;
    // Strict contact precedes iff its corrected x lies toward the ray's
    // source: corr ≤ 0 for a RIGHT-traveling ray, corr ≥ 0 for LEFT.
    const bool strict_first = (dir == RIGHT) ? (corr <= 0) : (corr >= 0);
    return new_is_strict == strict_first;
}

} // namespace chazelle
