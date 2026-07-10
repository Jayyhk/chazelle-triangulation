// src/visibility/up_phase.cpp — [C91 §4.1 tex 336–357]: The Up-Phase.

#include "up_phase.h"

#include "../merge/granularity.h"   // enforce_granularity ([C91 §3.3])

#include <algorithm>
#include <utility>
#include <vector>

namespace chazelle {

namespace {

// ════════════════════════════════════════════════════════════════
//  Dyadic chain covers on the [C91 §4 tex 316] grid
// ════════════════════════════════════════════════════════════════

struct ChainRef {
    std::size_t grade;
    std::size_t index;
};

// Decompose the P-global vertex interval [lo, hi] (lo < hi) into
// chains of the [C91 §4 tex 316] grid, in ascending position order,
// with at most two chains per grade and every grade < `max_grade_excl`
// — the classic two-pointer dyadic decomposition, O(max_grade_excl)
// time.  Requires hi − lo ≤ 2^{max_grade_excl}.  This is Lemma 4.1's
// partition: "In O(λ) time we can partition D into j ≤ 2λ chains,
// D₁, ..., D_j, in grades less than λ, with at most two chains per
// grade" ([C91 §4.1 tex 349]).
std::vector<ChainRef> dyadic_cover(std::size_t lo, std::size_t hi,
                                   std::size_t max_grade_excl) {
    assert(lo < hi && "cover needs a nonempty interval");
    assert(max_grade_excl >= 1 &&
           hi - lo <= (std::size_t{1} << max_grade_excl) &&
           "[C91 §4.1 tex 349]: the interval must fit within one "
           "grade-max_grade_excl chain length");
    std::vector<ChainRef> front, back;
    for (std::size_t mu = 0; lo < hi; ++mu) {
        const std::size_t step = std::size_t{1} << mu;
        if (mu + 1 == max_grade_excl) {
            // After the lower grades, lo and hi are step-aligned and
            // the remainder is one or two top-grade chains.
            assert(lo % step == 0 && hi % step == 0 &&
                   hi - lo <= 2 * step &&
                   "dyadic remainder is 1–2 aligned top-grade chains");
            for (; lo < hi; lo += step)
                front.push_back({mu, lo / step});
            break;
        }
        if ((lo / step) % 2 == 1) {
            front.push_back({mu, lo / step});
            lo += step;
        }
        if ((hi / step) % 2 == 1) {
            back.push_back({mu, hi / step - 1});
            hi -= step;
        }
    }
    for (std::size_t k = back.size(); k-- > 0; ) front.push_back(back[k]);
#ifndef NDEBUG
    for (std::size_t k = 0; k + 1 < front.size(); ++k)
        assert(((front[k].index + 1) << front[k].grade) ==
                   (front[k + 1].index << front[k + 1].grade) &&
               "cover chains are contiguous in ascending order");
#endif
    return front;
}

// Is the grid chain (grade, index) free of nonnull edges?  [C91 §2.2
// tex 106]: only [C91 §4 tex 316]'s padding produces null edges.
bool chain_all_null(const Polygon& P, std::size_t grade,
                    std::size_t index) {
    const std::size_t lo = index << grade;
    return P.count_nonnull_edges(lo, lo + (std::size_t{1} << grade) - 1)
           == 0;
}

// [C91 §4.1 tex 352]: the chain cover of one vertex-to-vertex piece —
// dyadic chains with every NONNULL-containing chain in grade ≤
// cap = ⌈βλ⌉ ("a collection of O(λ) chains in grades at most ⌈βλ⌉";
// the nonnull reading of tex 350's edge counts admits ALL-NULL padding
// chains of any grade < λ, whose canonical submaps are chordless and
// h-granular for every h — see up_phase.h).
//
// The base cover has ≤ 2 chains per grade; an all-real chain within a
// piece of ≤ 2^{cap} nonnull edges necessarily has grade ≤ cap, so at
// most ONE chain — the one containing P's real/padding boundary — can
// violate the cap.  Splitting it sheds one all-null right half per
// level (padding is a SUFFIX of P, so the nonnull edges are a prefix
// of the mixed chain) until the left part reaches grade cap.  O(λ)
// chains total (UpPhaseArcCutter::g_bound).
std::vector<ChainRef> piece_cover(const Polygon& P, std::size_t lo,
                                  std::size_t hi, std::size_t lambda,
                                  std::size_t cap) {
    assert(cap < lambda &&
           "[C91 §4.1 tex 352]: ⌈βλ⌉ < λ (early grades are naive)");
    assert(P.count_nonnull_edges(lo, hi - 1) <=
               (std::size_t{1} << cap) &&
           "[C91 §4.1 tex 350]: a piece has at most 2^{⌈βλ⌉} (nonnull) "
           "edges — the input submap is 2^{⌈βλ⌉}-semigranular");
    std::vector<ChainRef> cover = dyadic_cover(lo, hi, lambda);
    std::vector<ChainRef> out;
    out.reserve(cover.size() + lambda + 2);
    for (const ChainRef& c : cover) {
        if (c.grade <= cap || chain_all_null(P, c.grade, c.index)) {
            out.push_back(c);
            continue;
        }
        // The single mixed over-cap chain.
        ChainRef cur = c;
        std::vector<ChainRef> rights;
        while (cur.grade > cap) {
            assert(!chain_all_null(P, cur.grade, cur.index) &&
                   "over-cap splitting only reaches mixed chains");
            ChainRef left{cur.grade - 1, 2 * cur.index};
            ChainRef right{cur.grade - 1, 2 * cur.index + 1};
            // [C91 §4 tex 316]: padding is a suffix of P ⟹ the mixed
            // chain's nonnull edges are a prefix ⟹ the left half keeps
            // them all whenever the right half has any null edge...
            // more simply, the LEFT half always contains nonnull edges.
            assert(!chain_all_null(P, left.grade, left.index) &&
                   "[C91 §4 tex 316]: the real prefix stays in the "
                   "left half");
            // An over-cap ALL-REAL left half would put > 2^{cap}
            // nonnull edges inside the piece — excluded by the
            // semigranularity assert above, so descending terminates
            // at grade cap.
            rights.push_back(right);
            cur = left;
        }
        out.push_back(cur);
        for (std::size_t k = rights.size(); k-- > 0; )
            out.push_back(rights[k]);
    }
#ifndef NDEBUG
    for (std::size_t k = 0; k + 1 < out.size(); ++k)
        assert(((out[k].index + 1) << out[k].grade) ==
                   (out[k + 1].index << out[k + 1].grade) &&
               "piece cover stays contiguous and ascending");
    for (const ChainRef& c : out)
        assert((c.grade <= cap || chain_all_null(P, c.grade, c.index)) &&
               "[C91 §4.1 tex 352]: nonnull-containing chains have "
               "grade ≤ ⌈βλ⌉");
#endif
    return out;
}

// ════════════════════════════════════════════════════════════════
//  The tex 350–352 subarc decomposition (shared by both oracles)
// ════════════════════════════════════════════════════════════════

// One piece of the decomposition of α' — either a single-segment
// boundary piece (a partial edge at one of α''s endpoints, tex 350:
// "single line segments (at most two of them)") or a grid chain.
struct DecompPiece {
    bool boundary = false;
    Subarc subarc{};            // Cᵢ-local edges, exact endpoint ys
    ChainRef chain{NONE, NONE}; // valid iff !boundary
};

// Traversal-start / traversal-end vertices of the ∂C face (edge, side)
// ([C91 §2.4(iii) tex 138]: LEFT follows the edge, RIGHT reverses it).
std::size_t face_trav_start(std::size_t edge, Side side) {
    return (side == LEFT) ? edge : edge + 1;
}
std::size_t face_trav_end(std::size_t edge, Side side) {
    return (side == LEFT) ? edge + 1 : edge;
}

// [C91 §4.1 tex 350–352]: decompose the subarc α' (Cᵢ-local, exact
// endpoint ys) into its ≤ 2 boundary pieces + O(λ) grid chains, in
// clockwise ∂Cᵢ traversal order.  O(λ) time: subarc_legs is O(1), the
// covers are O(λ), and null runs are never walked edge-by-edge.
std::vector<DecompPiece> decompose_subarc(const UpPhase& up,
                                          const Polygon& Ci,
                                          const Subarc& target,
                                          std::size_t lambda) {
    assert(target.first_y.tag != SOS_NONE &&
           target.last_y.tag != SOS_NONE &&
           "[C91 §3.0(i) tex 169]: α' is specified by its two exact "
           "endpoints");
    const Polygon& P = up.graded().curve();
    const std::size_t off = Ci.table_offset();
    const std::size_t cap = ceil_beta(lambda);

    ArcLeg legs[3];
    const std::size_t nl = subarc_legs(target, 0, Ci.num_vertices() - 1,
                                       legs);
    assert(nl >= 1 && nl <= 3 &&
           "[C91 §4.1 tex 350]: a subarc has a constant number of "
           "no-double-backing pieces (1–3 legs, [C91 §2.4 tex 142])");

    std::vector<DecompPiece> out;
    for (std::size_t g = 0; g < nl; ++g) {
        const Side s = legs[g].side;
        const std::size_t lo = legs[g].lo, hi = legs[g].hi;
        // Traversal-first/-last faces of the leg (LEFT ascends edges,
        // RIGHT descends).
        const std::size_t f_first = (s == LEFT) ? lo : hi;
        const std::size_t f_last  = (s == LEFT) ? hi : lo;

        // A leg boundary is interior (at one of Cᵢ's endpoint
        // turnarounds — a vertex) unless it is α''s own endpoint, which
        // may sit mid-face: tex 350's "at most two" single segments.
        const bool first_partial =
            (g == 0) &&
            !symbolic_y_equal(target.first_y,
                              symbolic_y_of(Ci.vertex(
                                  face_trav_start(f_first, s))));
        const bool last_partial =
            (g + 1 == nl) &&
            !symbolic_y_equal(target.last_y,
                              symbolic_y_of(Ci.vertex(
                                  face_trav_end(f_last, s))));

        // Single-face leg with both ends partial: one boundary piece.
        if (lo == hi && first_partial && last_partial) {
            DecompPiece bp;
            bp.boundary = true;
            bp.subarc = Subarc{lo, s, lo, s, target.first_y,
                               target.last_y};
            out.push_back(bp);
            continue;
        }

        // Leading boundary piece: from α''s first endpoint to the
        // first face's traversal-end vertex.
        if (first_partial) {
            DecompPiece bp;
            bp.boundary = true;
            bp.subarc = Subarc{f_first, s, f_first, s, target.first_y,
                               symbolic_y_of(Ci.vertex(
                                   face_trav_end(f_first, s)))};
            out.push_back(bp);
        }

        // The vertex-to-vertex remainder, as a P-global vertex range.
        std::size_t vlo = lo + (first_partial && s == LEFT ? 1 : 0) +
                          (last_partial && s == RIGHT ? 1 : 0);
        std::size_t vhi = hi + 1 - (last_partial && s == LEFT ? 1 : 0) -
                          (first_partial && s == RIGHT ? 1 : 0);
        if (vlo < vhi) {
            std::vector<ChainRef> cover = piece_cover(
                P, off + vlo, off + vhi, lambda, cap);
            // Emit in traversal order: ascending for LEFT, descending
            // for RIGHT.
            auto emit = [&](const ChainRef& c) {
                const std::size_t glo = c.index << c.grade;
                const std::size_t ghi = glo +
                                        (std::size_t{1} << c.grade) - 1;
                DecompPiece cp;
                cp.chain = c;
                const std::size_t elo = glo - off;   // Cᵢ-local edges
                const std::size_t ehi = ghi - off;
                if (s == LEFT) {
                    cp.subarc = Subarc{
                        elo, s, ehi, s,
                        symbolic_y_of(Ci.vertex(elo)),
                        symbolic_y_of(Ci.vertex(ehi + 1))};
                } else {
                    cp.subarc = Subarc{
                        ehi, s, elo, s,
                        symbolic_y_of(Ci.vertex(ehi + 1)),
                        symbolic_y_of(Ci.vertex(elo))};
                }
                out.push_back(cp);
            };
            if (s == LEFT)
                for (const ChainRef& c : cover) emit(c);
            else
                for (std::size_t k = cover.size(); k-- > 0; )
                    emit(cover[k]);
        }

        // Trailing boundary piece: from the last face's traversal-start
        // vertex to α''s last endpoint.
        if (last_partial) {
            DecompPiece bp;
            bp.boundary = true;
            bp.subarc = Subarc{f_last, s, f_last, s,
                               symbolic_y_of(Ci.vertex(
                                   face_trav_start(f_last, s))),
                               target.last_y};
            out.push_back(bp);
        }
    }

    // Subdivision identity α' = α₁ ∪ ... ∪ α_k ([C91 §3.0(ii) tex 170]).
    assert(!out.empty());
    assert(out.front().subarc.first_edge == target.first_edge &&
           out.front().subarc.first_side == target.first_side &&
           symbolic_y_equal(out.front().subarc.first_y, target.first_y) &&
           "decomposition starts at α'.first");
    assert(out.back().subarc.last_edge == target.last_edge &&
           out.back().subarc.last_side == target.last_side &&
           symbolic_y_equal(out.back().subarc.last_y, target.last_y) &&
           "decomposition ends at α'.last");
    return out;
}

// ════════════════════════════════════════════════════════════════
//  On-α' membership with vertex-label canonicalization
// ════════════════════════════════════════════════════════════════

// [C91 §3.0(i) tex 169]: a report names "the edge of P that contains"
// the point hit — a contact at a NON-extremum interior vertex of Cᵢ is
// contained in BOTH incident edges, so target membership tries the
// adjacent label too (same rewrite rule as rebuild_submap's
// canonicalization); on success the hit's edge is rewritten to the
// label lying on α', keeping the report consistent with the target.
// At a y-extremum the same-side labels name DISTINCT duplication
// companions ([C91 §2.1 tex 72]) and are never exchanged.
bool target_contains_hit(const Subarc& target, const Polygon& Ci,
                         RayHit* h, SymbolicY sy) {
    const std::size_t cs = 0, ce = Ci.num_vertices() - 1;
    if (subarc_contains_point(target, Ci, h->edge, h->side, sy, cs, ce))
        return true;
    const auto& ed = Ci.edge(h->edge);
    std::size_t vidx = NONE;
    if (symbolic_y_equal(sy, symbolic_y_of(Ci.vertex(ed.start_idx))))
        vidx = ed.start_idx;
    else if (symbolic_y_equal(sy, symbolic_y_of(Ci.vertex(ed.end_idx))))
        vidx = ed.end_idx;
    if (vidx == NONE || Ci.is_endpoint(vidx) || Ci.is_y_extremum(vidx))
        return false;
    const std::size_t alt = (h->edge == vidx) ? vidx - 1 : vidx;
    if (alt == h->edge) return false;
    if (subarc_contains_point(target, Ci, alt, h->side, sy, cs, ce)) {
        h->edge = alt;
        return true;
    }
    return false;
}

} // namespace

// ════════════════════════════════════════════════════════════════
//  UpPhaseRayShooter — [C91 §3.0(i) tex 169] per tex 351–352
// ════════════════════════════════════════════════════════════════

UpPhaseRayShooter::UpPhaseRayShooter(const UpPhase& up, const Submap& Si,
                                     const Polygon& Ci,
                                     std::size_t lambda)
    : up_(&up), Si_(&Si), Ci_(&Ci), lambda_(lambda) {
    assert(ceil_beta(lambda) < lambda &&
           "[C91 §4.1 tex 352]: ⌈βλ⌉ < λ (early grades are naive)");
}

RayHit UpPhaseRayShooter::shoot(Point p, Side direction,
                                std::size_t arc_idx,
                                const Subarc& target,
                                double source_x_offset) const {
    assert(arc_idx < Si_->num_arcs() && !Si_->arc(arc_idx).dead &&
           "[C91 §3.0 tex 166]: α is specified by its arc-structure");
    assert_subarc_clockwise(target);
    const SymbolicY sy{p.y, p.index};
    const std::size_t off = Ci_->table_offset();

    std::vector<DecompPiece> pieces =
        decompose_subarc(*up_, *Ci_, target, lambda_);

    // First contact with ᾱ' = the lexicographic (wrapped, distance)
    // minimum over the pieces' first contacts ([C91 §4.1 tex 351]: "we
    // shoot toward each of the O(λ) subarcs of its decomposition and
    // determine the closest hit"; the obstacles ᾱⱼ partition ᾱ', so
    // the minimum over per-piece first contacts IS the first contact
    // with ᾱ', [C91 §2.1 tex 70]).
    RayHit best;
    best.hit = false;
    double best_d = 0.0;
    auto offer = [&](const RayHit& h) {
        assert(h.hit);
        const double d = (direction == RIGHT) ? (h.x - p.x)
                                              : (p.x - h.x);
        // [C91 §2.1 tex 70/§2 tex 47]: d < 0 ⟹ behind (wrapped);
        // d > 0 ⟹ forward (direct); raw d == 0 goes either way — the
        // perturbed x-offsets decide (perturbed_hit_forward).
        assert((d < 0.0 ? h.wrapped : (d > 0.0 ? !h.wrapped : true)) &&
               "[C91 §2.1 tex 70]: wrap flag consistent with the signed "
               "travel distance");
        bool better;
        if (!best.hit)
            better = true;
        else if (h.wrapped != best.wrapped)
            better = !h.wrapped;
        else if (d != best_d)
            better = d < best_d;
        else
            // Coincident raw positions across pieces (a shared chain
            // vertex, or the padding cluster) — [C91 §2 tex 47]
            // (polygon.h::ray_contact_precedes), in Cᵢ's frame.
            better = ray_contact_precedes(*Ci_, sy, direction, h.edge,
                                          h.side, best.edge, best.side);
        if (better) {
            best = h;
            best_d = d;
        }
    };

    for (const DecompPiece& piece : pieces) {
        if (piece.boundary) {
            // [C91 §4.1 tex 352]: "Shooting toward a single-edge subarc
            // is trivial."  The obstacle is the partial-edge SEGMENT
            // between the piece's two endpoint ys — it blocks from both
            // ∂P sides, so containment is tested on the ray's level,
            // not the struck side.
            const std::size_t e = piece.subarc.first_edge;
            double x;
            if (!edge_crossing_x(*Ci_, e, sy, &x)) continue;
            const SymbolicY ya = piece.subarc.first_y;
            const SymbolicY yb = piece.subarc.last_y;
            const bool within =
                (symbolic_y_leq(ya, sy) && symbolic_y_leq(sy, yb)) ||
                (symbolic_y_leq(yb, sy) && symbolic_y_leq(sy, ya));
            if (!within) continue;
            const auto& ed = Ci_->edge(e);
            const bool asc = symbolic_y_less(
                symbolic_y_of(Ci_->vertex(ed.start_idx)),
                symbolic_y_of(Ci_->vertex(ed.end_idx)));
            const Side minus_x = asc ? LEFT : RIGHT;
            RayHit h;
            h.hit = true;
            h.x = x;
            h.y = p.y;
            h.edge = e;
            h.side = (direction == RIGHT)
                         ? minus_x
                         : (minus_x == LEFT ? RIGHT : LEFT);
            {
                const double bd =
                    (direction == RIGHT) ? (x - p.x) : (p.x - x);
                h.wrapped = (bd < 0.0) ||
                            (bd == 0.0 &&
                             !(has_source_offset(source_x_offset) &&
                               perturbed_hit_forward(*Ci_, sy, direction,
                                                   source_x_offset, e)));
            }
            offer(h);
        } else {
            // [C91 §4.1 tex 351]: "Shooting toward any other subarc
            // makes use of the shooting structure of a canonical submap
            // for a chain in grade μ ≤ ⌈βλ⌉" (all-null padding chains
            // may sit higher; their structures are trivial).  The
            // chain's curve IS the piece's ᾱⱼ, so the Lemma 3.6 report
            // is the tex-169 first contact with the obstacle ᾱⱼ.
            const RayShootingStructure& rs =
                up_->chain_structure(piece.chain.grade, piece.chain.index);
            // The source's perturbed x-offset is plain geometry,
            // identical in every piece frame — it threads through
            // unchanged.
            RayHit h = rs.shoot_toward_boundary(p, direction,
                                                source_x_offset);
            if (!h.hit) continue;
            // Chain-local → Cᵢ-local edge frame.
            h.edge = (piece.chain.index << piece.chain.grade) + h.edge
                     - off;
            offer(h);
        }
    }

    if (!best.hit) return RayHit{};
    // [C91 §3.0(i) tex 169 + §3.2 tex 246]: the report is the first
    // contact of the ray with the OBSTACLE ᾱ' — "the point hit" plus
    // "the name of the edge of P that contains it" — and "does not
    // tell us if the point hit is on the desired arc or is the
    // companion of a point of the arc.  We already discussed how local
    // checking can decide which way it is": whether the contact lies
    // on α' proper is the CALLER's local check ([C91 §3.1 tex 220]'s
    // case (i) orientation inference and [C91 §3.1 tex 222]'s explicit
    // "proper orientation" disqualification depend on receiving the
    // companion contacts).  The report is therefore NOT filtered to
    // α'; only the edge label is canonicalized toward the target when
    // the contact sits at a non-extremum vertex contained in both
    // incident edges (tex 169's label ambiguity).
    (void)target_contains_hit(target, *Ci_, &best, sy);
    return best;
}

// ════════════════════════════════════════════════════════════════
//  UpPhaseArcCutter — [C91 §3.0(ii) tex 170] per tex 352–353
// ════════════════════════════════════════════════════════════════

UpPhaseArcCutter::UpPhaseArcCutter(const UpPhase& up, const Submap& Si,
                                   const Polygon& Ci, std::size_t lambda)
    : up_(&up), Si_(&Si), Ci_(&Ci), lambda_(lambda) {
    assert(ceil_beta(lambda) < lambda &&
           "[C91 §4.1 tex 352]: ⌈βλ⌉ < λ (early grades are naive)");
}

std::vector<ArcPiece> UpPhaseArcCutter::cut(std::size_t arc_idx,
                                            const Subarc& target) const {
    assert(arc_idx < Si_->num_arcs() && !Si_->arc(arc_idx).dead &&
           "[C91 §3.0 tex 166]: α is specified by its arc-structure");
    assert_subarc_clockwise(target);
    const std::size_t off = Ci_->table_offset();
    const std::size_t h_gamma = h_bound(lambda_);

    std::vector<DecompPiece> pieces =
        decompose_subarc(*up_, *Ci_, target, lambda_);

    std::vector<ArcPiece> out;
    out.reserve(pieces.size());
    for (const DecompPiece& piece : pieces) {
        ArcPiece ap;
        ap.subarc = piece.subarc;
        if (piece.boundary) {
            // [C91 §3.0(ii)(3) tex 170]: "except for ᾱ₁ and ᾱ₂ (in the
            // case where these are single-edge pieces attached to the
            // points of C corresponding to the endpoints of α')".
            ap.is_boundary_piece = true;
        } else {
            const std::size_t mu = piece.chain.grade;
            const std::size_t idx = piece.chain.index;
            ap.curve = &up_->graded().chain(mu, idx);
            ap.submap = &up_->chain_submap(mu, idx);
            // [C91 §4.1 tex 352]: "conformal submaps of granularity at
            // most 2^{⌈β⌈βλ⌉⌉}" — a grade-μ chain's canonical submap is
            // 2^{⌈βμ⌉}-granular (tex 338), ≤ h for μ ≤ ⌈βλ⌉; an
            // all-null padding chain's canonical submap is chordless
            // (total weight 0) and therefore h-granular outright.
            if (mu <= ceil_beta(lambda_)) {
                ap.granularity = UpPhase::grade_gamma(mu);
            } else {
                assert(ap.submap->num_live_chords() == 0 &&
                       "[C91 §4 tex 316]: an over-cap cover chain is "
                       "all-null padding — its canonical submap is "
                       "chordless");
                ap.granularity = h_gamma;
            }
            assert(ap.granularity <= h_gamma &&
                   "[C91 §4.1 tex 353]: h(γ) ≤ 2^{⌈β⌈βλ⌉⌉}");
            // The piece's curve is the vertex-to-vertex subchain of Cᵢ
            // over its edge range ([C91 §3.0(ii)(3) tex 170]).
            assert(ap.curve->table_offset() == (idx << mu) &&
                   ap.curve->num_edges() == (std::size_t{1} << mu) &&
                   "[C91 §4 tex 316]: grid chain spans its aligned "
                   "vertex range");
            assert((idx << mu) >= off &&
                   "chain lies within Cᵢ's table range");
        }
        out.push_back(ap);
    }

    // [C91 §3.0(ii) tex 170]: "in O(g(γᵢ)) time, subdivides the subarc
    // α' into at most g(γᵢ) subarcs."
    assert(out.size() <= g_bound(lambda_) &&
           "[C91 §4.1 tex 353]: g(γ) = O(λ) pieces");
    return out;
}

// ════════════════════════════════════════════════════════════════
//  Null-tail extension — [C91 §4 tex 316] + [C91 §2.2 tex 106]
// ════════════════════════════════════════════════════════════════
//
// The padding is one all-null SUFFIX of P (pad_curve), so in any
// dyadic partition every part at or past the real/padding boundary is
// all-null.  Merging with an all-null part is trivial: the tail is
// geometrically the single point at C's end, occupied by the curve
// already, so V(C ∪ tail) has exactly V(C)'s nonnull-relevant
// structure — the tail's own vertices contribute only weight-0 chords,
// none of which survives γ ≥ 1 granularity.  A canonical submap of
// V(C ∪ tail) is therefore the canonical submap of V(C) with
//
//  · end_vertex moved past the tail (arc legs re-derive their spans
//    from it, [C91 §2.4 tex 142], so the arc through C's end
//    turnaround silently absorbs the tail's ∂-faces — its nonnull
//    counts are unchanged, [C91 §2.2 tex 106]);
//
//  · when C's end vertex becomes a y-EXTREMUM of the union (the last
//    real edge ascends into it; the null branch is symbolically below,
//    [C91 §2 tex 47]), the old endpoint-companion structure ([C91 §2.1
//    tex 72] case 3) becomes the 4-duplicate case 2: the companion on
//    the last real edge's INSIDE face turns into an inside duplicate,
//    and any chord anchored at the old outside companion of that side
//    transfers to the tail's first null edge's outer face — same
//    geometric sight, new ∂C name.  The turn-crossing arc splits into
//    the apex's zero-length cap and the new end-wrap arc through the
//    tail's turnaround.
//
// This keeps every actual merge() of Lemma 4.1 between curves with
// real geometry and a REAL junction vertex; no fusion ever walks the
// padding cluster.

namespace {

// Extend S — a canonical submap of V(curve) — over the all-null grid
// range adjacent to its end, up to P's vertex `tail_hi_vertex`
// (inclusive).  See the header comment above; the caller re-asserts
// canonicalness.
void extend_with_null_tail(Submap& S, const Polygon& P,
                           Polygon& curve, std::size_t tail_hi_vertex) {
    const std::size_t old_end_local = curve.num_vertices() - 1;
    const std::size_t glob_lo = curve.table_offset();
    assert(tail_hi_vertex > glob_lo + old_end_local &&
           "extension must add at least one tail vertex");
    assert(P.count_nonnull_edges(glob_lo + old_end_local,
                                 tail_hi_vertex - 1) == 0 &&
           "[C91 §4 tex 316]: the tail is all-null");
    assert(curve.num_vertices() >= 2);

    const std::size_t e_r = old_end_local - 1;   // last old edge
    const Point& v_end = curve.vertex(old_end_local);
    const Point& v_prev = curve.vertex(old_end_local - 1);
    // The null branch is symbolically BELOW v_end ([C91 §2 tex 47]:
    // larger tag → lower perturbed y), so v_end becomes a y-MAXIMUM of
    // the union exactly when the last old edge ascends into it.
    const bool becomes_max = point_y_above(v_end, v_prev);

    if (becomes_max) {
        // [C91 §2.1 tex 72]: case 3 → case 2 at v_end.  e_r ascends ⟹
        // its east (inside) face is RIGHT; the tail's first null edge
        // descends symbolically ⟹ its outer (east) face is LEFT.  A
        // chord anchored at the old east companion (e_r, RIGHT, y_end)
        // transfers to the outer face of the tail's first edge — the
        // same geometric sight under its new ∂C name — and the arcs
        // bounded there re-anchor with it: the arc ENDING there becomes
        // the apex's zero-length cap, the arc STARTING there absorbs
        // the tail detour (its legs re-derive from the moved
        // end_vertex, [C91 §2.4 tex 142]).
        const SymbolicY ye = symbolic_y_of(v_end);
        const std::size_t e_null = old_end_local;    // first tail edge

        // Identify the bounded arcs BEFORE any rewrite (the start/end
        // derivations read the pre-surgery labels).
        std::size_t arc_ending = NONE, arc_starting = NONE;
        for (std::size_t ai = 0; ai < S.num_arcs(); ++ai) {
            const Arc& a = S.arc(ai);
            if (a.dead) continue;
            if (a.last_edge == e_r && a.last_side == RIGHT &&
                symbolic_y_equal(S.arc_end_symbolic_y(ai, curve), ye)) {
                assert(arc_ending == NONE);
                arc_ending = ai;
            }
            if (a.first_edge == e_r && a.first_side == RIGHT &&
                symbolic_y_equal(S.arc_start_symbolic_y(ai, curve),
                                 ye)) {
                assert(arc_starting == NONE);
                arc_starting = ai;
            }
        }

        std::size_t moved_chord = NONE;
        for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
            Chord& c = S.chord(ci);
            if (c.dead) continue;
            if (!symbolic_y_equal(c.symbolic_y(), ye)) continue;
            if (c.left_edge == e_r && c.left_side == RIGHT) {
                assert(moved_chord == NONE &&
                       "[C91 §2.1 tex 70]: one chord per ∂C point");
                c.left_edge = e_null;
                c.left_side = LEFT;
                moved_chord = ci;
            } else if (c.right_edge == e_r && c.right_side == RIGHT) {
                assert(moved_chord == NONE &&
                       "[C91 §2.1 tex 70]: one chord per ∂C point");
                c.right_edge = e_null;
                c.right_side = LEFT;
                moved_chord = ci;
            }
        }

        if (moved_chord != NONE) {
            assert(arc_ending != NONE && arc_starting != NONE &&
                   "[C91 §2.2 tex 96]: a chord endpoint bounds an "
                   "ending and a starting arc");
            S.arc(arc_ending).last_edge = e_null;
            S.arc(arc_ending).last_side = LEFT;
            S.arc(arc_starting).first_edge = e_null;
            S.arc(arc_starting).first_side = LEFT;
            // The arc through the (new) end turnaround is the one
            // absorbing the tail detour.
            S.end_arc = arc_starting;
        } else {
            assert(arc_ending == NONE && arc_starting == NONE &&
                   "[C91 §2.2 tex 96]: arcs are bounded at chord "
                   "endpoints only");
        }
    }

    curve = P.subchain(glob_lo, tail_hi_vertex + 1 - glob_lo);
    S.end_vertex = curve.num_vertices() - 1;
    // [C91 §2.2 tex 106]: nonnull counts are unchanged by the all-null
    // tail; the per-leg refresh re-derives the (possibly re-anchored)
    // spans.
    S.refresh_arc_edge_counts(curve);
    S.build_tree_decomposition();
}

// The chordless canonical submap of an ALL-NULL portion ([C91 §4 tex
// 316]): total weight 0, so every chord of V is removable under any
// γ ≥ 1 ([C91 §2.3 tex 121]) — the canonical submap is the single
// closed arc over all of ∂C ([C91 §2.4 tex 142]).
Submap all_null_canonical(const Polygon& curve) {
    Submap s;
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = curve.num_vertices() - 1;
    Arc a{};
    a.first_edge = 0;
    a.last_edge = 0;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = 0;
    [[maybe_unused]] std::size_t ai = s.add_arc(a);
    assert(s.start_arc == ai && s.end_arc == ai);
    s.build_tree_decomposition();
    return s;
}

} // namespace

// ════════════════════════════════════════════════════════════════
//  Lemma 4.1 — [C91 §4.1 tex 347–356]
// ════════════════════════════════════════════════════════════════

UpPhase::PortionResult UpPhase::compute_canonical_portion(
        std::size_t a, std::size_t b) const {
    const Polygon& P = graded_.curve();
    assert(a < b && b < P.num_vertices() &&
           "[C91 Lemma 4.1 tex 347]: D = v_a, ..., v_b within P");

    // λ from 2^{λ−1} < b − a ≤ 2^λ (tex 347).
    std::size_t lambda = 0;
    while ((std::size_t{1} << lambda) < b - a) ++lambda;
    assert(lambda >= 1 &&
           (b - a) > (std::size_t{1} << (lambda - 1)) - 1 &&
           (b - a) <= (std::size_t{1} << lambda) &&
           "[C91 Lemma 4.1 tex 347]: 2^{λ−1} < b − a ≤ 2^λ");
    assert(ceil_beta(lambda) < lambda &&
           "[C91 §4.1 tex 352]: Lemma 4.1 requires ⌈βλ⌉ < λ — smaller "
           "portions are the naive early grades (tex 345)");

    // γ = 2^{⌈βλ⌉}, "the granularity of a canonical submap of V(D)"
    // (tex 349).
    const std::size_t gamma = grade_gamma(lambda);

    // tex 349: "In O(λ) time we can partition D into j ≤ 2λ chains,
    // D₁, ..., D_j, in grades less than λ, with at most two chains per
    // grade.  This implies that, for each i, a canonical submap Sᵢ of
    // V(Dᵢ) is available."
    std::vector<ChainRef> parts = dyadic_cover(a, b, lambda);
    assert(parts.size() >= 1 && parts.size() <= 2 * lambda &&
           "[C91 §4.1 tex 349]: j ≤ 2λ partition chains");

    // [C91 §4 tex 316]: the padding is one all-null SUFFIX of P, so the
    // all-null parts form a suffix of the partition.  They carry no
    // weight ([C91 §2.2 tex 106]) and are attached by the null-tail
    // extension instead of by merging — every actual merge below has a
    // REAL junction vertex.
    std::size_t n_real_parts = parts.size();
    while (n_real_parts > 0 &&
           chain_all_null(P, parts[n_real_parts - 1].grade,
                          parts[n_real_parts - 1].index))
        --n_real_parts;
#ifndef NDEBUG
    // [C91 §4 tex 316]: the padding is a suffix of P, so after
    // trimming the all-null suffix no all-null part remains.
    for (std::size_t k = 0; k < n_real_parts; ++k)
        assert(!chain_all_null(P, parts[k].grade, parts[k].index) &&
               "[C91 §4 tex 316]: all-null parts form a suffix of the "
               "dyadic partition");
#endif
    if (n_real_parts == 0) {
        // The whole portion lies inside the padding: its canonical
        // submap is the chordless closed arc (weight 0).
        PortionResult out{P.subchain(a, b - a + 1), Submap{}};
        out.submap = all_null_canonical(out.curve);
        assert(out.submap.is_granular(gamma, out.curve));
        return out;
    }

    // tex 349: "Since the granularity of canonical submaps grows
    // monotonically with the size of the underlying polygonal curve, we
    // can trivially reset the granularity of each Sᵢ to γ (Section
    // 3.3)."  Resetting operates on COPIES: the stored chain submaps
    // stay alive as this level's oracle structures (tex 341) and for
    // the down-phase ([C91 §4 tex 323]).
    struct Node {
        Polygon curve;
        Submap submap;
    };
    std::vector<Node> level;
    level.reserve(n_real_parts);
    for (std::size_t pi = 0; pi < n_real_parts; ++pi) {
        const ChainRef& part = parts[pi];
        const Polygon& c = graded_.chain(part.grade, part.index);
        Submap s = chain_submap(part.grade, part.index);   // copy
        assert(grade_gamma(part.grade) <= gamma &&
               "[C91 §4.1 tex 349]: canonical granularity grows "
               "monotonically with curve size");
        enforce_granularity(s, c, gamma);
        s.normalize(c);
        assert(s.is_granular(gamma, c) &&
               "[C91 §3.3 tex 276]: reset yields a γ-granular submap");
        level.push_back(Node{c, std::move(s)});
    }

    // tex 349–350: "Let us now merge these submaps two-by-two (D₁ with
    // D₂, D₃ with D₄, etc.).  More generally, we consider a perfectly
    // balanced binary tree whose leaves are in bijection with the Dᵢ
    // and we merge submaps bottom-up by following the tree pattern.
    // Application of Lemma 3.5 results in a canonical submap of V(D)."
    // There are O(log λ) levels (tex 355), an odd node riding up a
    // level unmerged.
    while (level.size() > 1) {
        std::vector<Node> next;
        next.reserve((level.size() + 1) / 2);
        for (std::size_t i = 0; i + 1 < level.size(); i += 2) {
            Node& L = level[i];
            Node& R = level[i + 1];

            // tex 350–353: the oracles, realized from the processed
            // lower grades' structures.
            UpPhaseRayShooter ray_l(*this, L.submap, L.curve, lambda);
            UpPhaseRayShooter ray_r(*this, R.submap, R.curve, lambda);
            UpPhaseArcCutter cut_l(*this, L.submap, L.curve, lambda);
            UpPhaseArcCutter cut_r(*this, R.submap, R.curve, lambda);

            MergeInput in;
            in.C1 = &L.curve;
            in.C2 = &R.curve;
            in.S1 = &L.submap;
            in.S2 = &R.submap;
            in.gamma1 = gamma;
            in.gamma2 = gamma;
            in.gamma = gamma;
            in.ray_shooter_1 = &ray_l;
            in.ray_shooter_2 = &ray_r;
            in.arc_cutter_1 = &cut_l;
            in.arc_cutter_2 = &cut_r;
            in.g_gamma1 = in.g_gamma2 = UpPhaseArcCutter::g_bound(lambda);
            in.h_gamma1 = in.h_gamma2 = UpPhaseArcCutter::h_bound(lambda);

            MergeResult res = merge(in);
            next.push_back(Node{res.C, std::move(res.S)});
        }
        if (level.size() % 2 == 1) next.push_back(std::move(level.back()));
        level = std::move(next);
    }

    PortionResult out{level[0].curve, std::move(level[0].submap)};

    if (n_real_parts < parts.size()) {
        // Attach the all-null suffix by extension (see above).
        extend_with_null_tail(out.submap, P, out.curve, b);
#ifdef CHAZELLE_EXPENSIVE_ASSERTS
        out.submap.check_invariants(out.curve);
#endif
    }

    assert(out.curve.num_vertices() == b - a + 1 &&
           out.curve.table_offset() == a &&
           "[C91 Lemma 4.1 tex 347]: the result covers exactly D");
    assert(out.submap.is_conformal() &&
           "[C91 Lemma 4.1 tex 350]: the portion's submap is conformal");
    assert(out.submap.is_granular(gamma, out.curve) &&
           "[C91 Lemma 4.1 tex 350]: the portion's submap is "
           "2^{⌈βλ⌉}-granular");
    assert(!out.submap.tree_decomposition().empty() &&
           "[C91 Lemma 4.1 tex 350]: normal form carries the tree "
           "decomposition");
    return out;
}

// ════════════════════════════════════════════════════════════════
//  The up-phase driver — [C91 §4.1 tex 340–345/356–357]
// ════════════════════════════════════════════════════════════════

UpPhase::UpPhase(std::vector<Point> vertices)
    : graded_(std::move(vertices)) {
    const std::size_t p = graded_.p();
    submaps_.resize(p + 1);
    structures_.resize(p + 1);

    // tex 340: "For λ = 0, 1, ..., p in that order, we process grade λ."
    for (std::size_t lam = 0; lam <= p; ++lam) {
        const std::size_t count = graded_.num_chains(lam);
        submaps_[lam].reserve(count);        // pointer stability for (ii)
        structures_[lam].reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            const Polygon& c = graded_.chain(lam, i);

            // (i) tex 341: a canonical submap of V(C) for each chain C.
            Submap s;
            if (lam <= NAIVE_MAX_GRADE) {
                // tex 345: "This work can be done naively for the
                // early grades."
                s = build_canonical_submap_naive(c);
            } else {
                // tex 346: "Lemma 4.1 can be called upon to compute a
                // canonical submap of the visibility map of each chain
                // in grade λ" (tex 356).
                PortionResult res = compute_canonical_portion(
                    i << lam, (i + 1) << lam);
                s = std::move(res.submap);
            }
            assert(s.is_conformal() &&
                   s.is_granular(grade_gamma(lam), c) &&
                   !s.tree_decomposition().empty() &&
                   "[C91 §4.1 tex 338]: canonical = 2^{⌈βλ⌉}-granular, "
                   "conformal, normal form");
            submaps_[lam].push_back(std::move(s));

            // (ii) tex 342: "We preprocess each canonical submap for
            // ray-shooting along the lines of Lemma 3.6, setting γ to
            // the value 2^{⌈βλ⌉}."
            structures_[lam].push_back(
                std::make_unique<RayShootingStructure>(
                    submaps_[lam][i], c, grade_gamma(lam)));
        }
    }

    // tex 323: the up-phase ends with "the submap for the whole
    // polygon" — the single grade-p chain is P itself.
    assert(submaps_[p].size() == 1 &&
           "[C91 §4 tex 319]: grade p has the single chain P");
}

} // namespace chazelle
