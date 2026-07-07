// tests/s3.2-tests.cpp — Tests for [C91 §3.2]: Restoring Conformality.
//
// Fixture: a six-peak "comb" C = C₁ ∪ C₂ with a hand-built fused submap
// that faithfully mimics a §3.1 fusion output (junction chords at the
// C₁ ∩ C₂ vertex, tent-pocket caps, zero-length arcs).  Its big outer
// region has SEVEN logical arcs — the [C91 §3.2 tex 238] situation —
// and the oracles are REAL geometric implementations, so the whole
// Lemma 3.2/3.3/3.4 pipeline runs against true visibility.

#include "merge/conformality.h"
#include "merge/fusion.h"
#include "polygon/polygon.h"
#include "submap/submap.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <functional>
#include <memory>
#include <sys/prctl.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace chazelle;

namespace {
bool assert_fires(std::function<void()> fn) {
    pid_t pid = fork();
    if (pid == 0) {
        prctl(PR_SET_DUMPABLE, 0);
        if (freopen("/dev/null", "w", stderr) == nullptr) std::_Exit(2);
        fn();
        std::_Exit(0);
    }
    int status = 0;
    if (waitpid(pid, &status, 0) < 0) return false;
    return WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
}
}

// ════════════════════════════════════════════════════════════════
//  The comb geometry
// ════════════════════════════════════════════════════════════════
//
//  C: 13 vertices, peaks P1..P6 at odd indices, valleys between.
//  Junction C₁ ∩ C₂ = vertex 7 (the tallest apex, P4).
//
//     idx:  0      1      2     3      4     5      6     7      8
//           (0,0) (2,20) (4,6) (6,24) (8,4) (10,22)(12,5)(14,26)(16,2)
//     idx:  9      10       11     12
//           (18,23)(20,4.5) (22,25)(24,1)

static Polygon make_C() {
    return Polygon({{0,0,0}, {2,20,1}, {4,6,2}, {6,24,3}, {8,4,4},
                    {10,22,5}, {12,5,6}, {14,26,7}, {16,2,8},
                    {18,23,9}, {20,4.5,10}, {22,25,11}, {24,1,12}});
}
static Polygon make_C1() {
    return Polygon({{0,0,0}, {2,20,1}, {4,6,2}, {6,24,3}, {8,4,4},
                    {10,22,5}, {12,5,6}, {14,26,7}});
}
static Polygon make_C2() {
    return Polygon({{14,26,7}, {16,2,8}, {18,23,9}, {20,4.5,10},
                    {22,25,11}, {24,1,12}});
}

// Chordless single-region normal-form submap: one region bounded by
// the single closed arc — all of ∂C, one arc-structure stored cut at
// C's start turnaround ([C91 §2.4 tex 142/138]).
static Submap make_chordless_normal(const Polygon& poly) {
    Submap s;
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = poly.num_vertices() - 1;
    Arc a{};
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count =
        2 * poly.count_nonnull_edges(0, poly.num_edges() - 1);
    std::size_t ai0 = s.add_arc(a);
    assert(s.start_arc == ai0 && s.end_arc == ai0);
    s.build_tree_decomposition();
    return s;
}

// ════════════════════════════════════════════════════════════════
//  The fused submap over the comb
// ════════════════════════════════════════════════════════════════
//
//  Chords (all genuine V(C) chords):
//    c0 v1w: (0.6,6)–(4,6)      under-P1 cap, sourced at v2
//    c1 v1e: (4,6)–(7.8,6)      under-P2 cap, sourced at v2
//    c2 v3w: (8.1̄,5)–(12,5)     under-P3 cap, sourced at v6
//    c3 v3e: (12,5)–(15.75,5)   under-P4 cap, sourced at v6
//    c4 v5w: (16.2̄,4.5)–(20,4.5) under-P5 cap, sourced at v10
//    c5 v5e: (20,4.5)–(23.708̄,4.5) under-P6 cap, sourced at v10
//    c6 cAPEX: junction outside-pair chord at the apex (zero length)
//    c7 n7:  junction inside-pair null-length chord ([C91 §2.1 tex 72])
//
//  Regions: 0 = R_big (7 arcs — not yet conformal), 1..6 = tent
//  pockets P1..P6, 7 = cAPEX's empty inner, 8 = n7's empty inner.
//
//  [C91 §2.4 tex 142]: R_big's two arcs through C's endpoint
//  turnarounds are single double-backing structures — E* around C's
//  end vertex (LEFT [7..11] + RIGHT [11]) and S* around its start
//  vertex (RIGHT [0] + LEFT [0..6]).
//
//  Table arc indices (canonical order [C91 §2.4(iii) tex 138]):
//    LEFT-starting:  0 Lz[7..7]∅, 1 E* (7,L)→(11,R)
//    RIGHT-starting: 2 P6[11..10], 3 Z5[9]∅, 4 P5[9..8], 5 R2[8..7],
//           6 P4a[7], 7 Nz[7]∅, 8 P4b[6], 9 Z3[5]∅, 10 P3[5..4],
//           11 R3[4..3], 12 P2[3..2], 13 Z1[1]∅, 14 P1[1..0],
//           15 S* (0,R)→(6,L)

struct CombFixture {
    Polygon C, C1, C2;
    Submap S1, S2;
    Submap S;
    CombFixture()
        : C(make_C()), C1(make_C1()), C2(make_C2()),
          S1(make_chordless_normal(C1)), S2(make_chordless_normal(C2)) {
        build();
    }

    void build() {
        for (int i = 0; i < 9; ++i) S.add_node();

        auto add = [&](std::size_t fe, Side fs, std::size_t le, Side ls,
                       std::size_t region, std::size_t count)
                -> std::size_t {
            Arc a{};
            a.first_edge = fe; a.first_side = fs;
            a.last_edge = le; a.last_side = ls;
            a.region_node = region;
            a.edge_count = count;
            return S.add_arc(a);
        };
        S.start_vertex = 0;
        S.end_vertex = 12;

        add(7, LEFT, 7, LEFT, 7, 0);      // 0  Lz  (cAPEX inner, ∅)
        add(7, LEFT, 11, RIGHT, 0, 6);    // 1  E*  (end wrap, ᾱ=[7,11])
        add(11, RIGHT, 10, RIGHT, 6, 2);  // 2  P6
        add(9, RIGHT, 9, RIGHT, 0, 0);    // 3  Z5  (v10 outside pair, ∅)
        add(9, RIGHT, 8, RIGHT, 5, 2);    // 4  P5
        add(8, RIGHT, 7, RIGHT, 0, 2);    // 5  R2
        add(7, RIGHT, 7, RIGHT, 4, 1);    // 6  P4a
        add(7, RIGHT, 7, RIGHT, 8, 0);    // 7  Nz  (n7 inner, ∅)
        add(6, RIGHT, 6, RIGHT, 4, 1);    // 8  P4b
        add(5, RIGHT, 5, RIGHT, 0, 0);    // 9  Z3  (v6 outside pair, ∅)
        add(5, RIGHT, 4, RIGHT, 3, 2);    // 10 P3
        add(4, RIGHT, 3, RIGHT, 0, 2);    // 11 R3
        add(3, RIGHT, 2, RIGHT, 2, 2);    // 12 P2
        add(1, RIGHT, 1, RIGHT, 0, 0);    // 13 Z1  (v2 outside pair, ∅)
        add(1, RIGHT, 0, RIGHT, 1, 2);    // 14 P1
        add(0, RIGHT, 6, LEFT, 0, 8);     // 15 S*  (start wrap, ᾱ=[0,6])

        auto chord = [&](std::size_t le, Side lsd, Chord::AdjArcs ladj,
                         std::size_t re, Side rsd, Chord::AdjArcs radj,
                         double y, std::size_t tag,
                         std::size_t r0, std::size_t r1,
                         bool null_len = false) {
            Chord c{};
            c.left_edge = le; c.left_side = lsd; c.left_adj = ladj;
            c.right_edge = re; c.right_side = rsd; c.right_adj = radj;
            c.y = y; c.y_tag = tag;
            c.region[0] = r0; c.region[1] = r1;
            c.is_null_length = null_len;
            S.add_chord(c);
        };
        // c0 v1w
        chord(0, RIGHT, {{14, 15}, 2}, 1, RIGHT, {{13}, 1}, 6.0, 2, 0, 1);
        // c1 v1e
        chord(2, RIGHT, {{12}, 1}, 3, RIGHT, {{11, 12}, 2}, 6.0, 2, 0, 2);
        // c2 v3w
        chord(4, RIGHT, {{10, 11}, 2}, 5, RIGHT, {{9}, 1}, 5.0, 6, 0, 3);
        // c3 v3e
        chord(6, RIGHT, {{8}, 1}, 7, RIGHT, {{5, 6}, 2}, 5.0, 6, 0, 4);
        // c4 v5w
        chord(8, RIGHT, {{4, 5}, 2}, 9, RIGHT, {{3}, 1}, 4.5, 10, 0, 5);
        // c5 v5e
        chord(10, RIGHT, {{2}, 1}, 11, RIGHT, {{1, 2}, 2}, 4.5, 10, 0, 6);
        // c6 cAPEX (junction outside pair; zero length, distinct ∂C points)
        chord(6, LEFT, {{15}, 1}, 7, LEFT, {{0}, 1}, 26.0, 7, 0, 7);
        // c7 n7 (junction inside pair, null length [C91 §2.1 tex 72])
        chord(7, RIGHT, {{6}, 1}, 7, RIGHT, {{7}, 1}, 26.0, 7, 4, 8, true);

        // [C91 §2.4(iii) tex 138]: add_arc auto-registered the wrap
        // arcs as the endpoint pointers.
        assert(S.start_arc == 15 && S.end_arc == 1);
    }
};

// ════════════════════════════════════════════════════════════════
//  Geometric oracles
// ════════════════════════════════════════════════════════════════

// True geometric ray shooter over Cᵢ, restricted to the target subarc.
// [C91 §2.1 tex 70]: a ray that misses everything wraps through the
// point at infinity; same-position crossings (the source's duplicates)
// are only reachable after a full wrap.
struct GeomRayShooter : RayShootingOracle {
    const Polygon* Ci;
    explicit GeomRayShooter(const Polygon* c) : Ci(c) {}

    static bool crossing_x(const Polygon& C, std::size_t e, SymbolicY sy,
                           double* x) {
        const auto& ed = C.edge(e);
        const Point& vs = C.vertex(ed.start_idx);
        const Point& ve = C.vertex(ed.end_idx);
        SymbolicY y0 = symbolic_y_of(vs);
        SymbolicY y1 = symbolic_y_of(ve);
        if (symbolic_y_equal(sy, y0)) { *x = vs.x; return true; }
        if (symbolic_y_equal(sy, y1)) { *x = ve.x; return true; }
        bool between =
            (symbolic_y_less(y0, sy) && symbolic_y_less(sy, y1)) ||
            (symbolic_y_less(y1, sy) && symbolic_y_less(sy, y0));
        if (!between) return false;
        double t = (sy.y - vs.y) / (ve.y - vs.y);
        *x = vs.x + t * (ve.x - vs.x);
        return true;
    }

    RayHit shoot(Point p, Side dir, std::size_t /*arc_idx*/,
                 const Subarc& target) const override {
        SymbolicY sy{p.y, p.index};
        // [C91 §2.4 tex 142]: wrap-spanning targets decompose per leg.
        ArcLeg legs[3];
        std::size_t nl = subarc_legs(target, 0, Ci->num_vertices() - 1,
                                     legs);

        RayHit best;
        best.hit = false;
        double best_d = 0.0;
        for (std::size_t g = 0; g < nl; ++g) {
            for (std::size_t e = legs[g].lo; e <= legs[g].hi; ++e) {
                double x;
                if (!crossing_x(*Ci, e, sy, &x)) continue;
                // A ray traveling RIGHT strikes −x-facing walls; LEFT,
                // +x-facing ([C91 §2.1 tex 72] double boundary).
                const auto& ed = Ci->edge(e);
                bool asc = symbolic_y_less(
                    symbolic_y_of(Ci->vertex(ed.start_idx)),
                    symbolic_y_of(Ci->vertex(ed.end_idx)));
                Side minus_x = asc ? LEFT : RIGHT;
                Side struck = (dir == RIGHT)
                    ? minus_x : (minus_x == LEFT ? RIGHT : LEFT);
                // [C91 §3.0(i) tex 169]: α' is endpoint-exact — skip
                // candidates off α' (wrong side OR beyond an endpoint
                // on a shared boundary edge) and keep scanning.
                if (!subarc_contains_point(target, *Ci, e, struck, sy,
                                           0, Ci->num_vertices() - 1))
                    continue;
                double d = (dir == RIGHT) ? (x - p.x) : (p.x - x);
                // d > 0: direct.  d ≤ 0: reachable only through the
                // wrap (d == 0 is the source's own duplicate position,
                // met after a full circuit).
                bool wrapped = (d <= 0.0);
                bool better;
                if (!best.hit) better = true;
                else if (wrapped != best.wrapped) better = !wrapped;
                else better = d < best_d;
                if (better) {
                    best.hit = true;
                    best.x = x;
                    best.y = p.y;
                    best.edge = e;
                    best.side = struck;
                    best.wrapped = wrapped;
                    best_d = d;
                }
            }
        }
        return best;
    }
};

// Geometric arc cutter over Cᵢ: subdivides the target into single-edge
// boundary pieces at mid-edge (partial) ends plus vertex-to-vertex
// middle chunks of at most `piece_len` edges, each with a chordless
// normal-form submap of V(ᾱ) — valid per [C91 §3.0(ii) tex 170] (a
// chordless submap is conformal and h-granular for h ≥ |ᾱ|).
struct GeomArcCutter : ArcCuttingOracle {
    const Polygon* Ci;
    const Polygon* C;                // full curve (for chord scan)
    const Submap* const* S_ptr;      // live fused submap
    std::size_t off;                 // Cᵢ edge → C edge
    std::size_t piece_len;
    mutable std::vector<std::unique_ptr<Polygon>> curves;
    mutable std::vector<std::unique_ptr<Submap>> submaps;

    GeomArcCutter(const Polygon* ci, const Polygon* c,
                  const Submap* const* sp, std::size_t o, std::size_t pl)
        : Ci(ci), C(c), S_ptr(sp), off(o), piece_len(pl) {}

    // Is there a live chord of the fused submap with a MID-EDGE endpoint
    // on (C-edge, side)?  Then a target ending on that edge may end
    // mid-edge, and the end edge conservatively becomes a boundary piece
    // ([C91 §3.0(ii)(3) tex 170]).  Conservative marking loses no
    // candidate: boundary pieces' vertices are still checked directly
    // (with the on-A₁ filter) by find_visible_point.
    bool mid_edge_endpoint(std::size_t edge_c, Side side) const {
        const Submap& S = **S_ptr;
        auto at_vertex = [&](SymbolicY y) {
            const auto& e = C->edge(edge_c);
            return symbolic_y_equal(y, symbolic_y_of(C->vertex(e.start_idx))) ||
                   symbolic_y_equal(y, symbolic_y_of(C->vertex(e.end_idx)));
        };
        for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
            const Chord& c = S.chord(ci);
            if (c.dead) continue;
            bool l = c.left_edge == edge_c && c.left_side == side &&
                     !at_vertex(c.symbolic_y());
            bool r = c.right_edge == edge_c && c.right_side == side &&
                     !at_vertex(c.symbolic_y());
            if (l || r) return true;
        }
        return false;
    }

    Submap build_piece_submap(const Polygon& alpha) const {
        return make_chordless_normal(alpha);
    }

    std::vector<ArcPiece> cut(std::size_t arc_idx,
                              const Subarc& target) const override {
        // [C91 §3.0(ii)(2) tex 170]: pieces must not double-back — a
        // wrap-spanning target is first split at Cᵢ's endpoint
        // turnaround(s) into its single-side legs ([C91 §2.4 tex 142]),
        // each cut recursively, pieces emitted in traversal order.
        if (target.first_side != target.last_side ||
            (target.first_side == LEFT &&
             target.first_edge > target.last_edge) ||
            (target.first_side == RIGHT &&
             target.first_edge < target.last_edge)) {
            ArcLeg legs[3];
            std::size_t nl = subarc_legs(target, 0,
                                         Ci->num_vertices() - 1, legs);
            std::vector<ArcPiece> out;
            for (std::size_t g = 0; g < nl; ++g) {
                Subarc leg_sub = (legs[g].side == LEFT)
                    ? Subarc{legs[g].lo, LEFT, legs[g].hi, LEFT}
                    : Subarc{legs[g].hi, RIGHT, legs[g].lo, RIGHT};
                auto pieces = cut(arc_idx, leg_sub);
                out.insert(out.end(), pieces.begin(), pieces.end());
            }
            return out;
        }
        const Side s = target.first_side;
        const std::size_t lo = std::min(target.first_edge, target.last_edge);
        const std::size_t hi = std::max(target.first_edge, target.last_edge);

        // Which geometric end (low/high edge) is the traversal
        // first/last: LEFT ascends, RIGHT descends.
        const std::size_t first_e = target.first_edge;
        const std::size_t last_e = target.last_edge;
        const bool first_partial = mid_edge_endpoint(off + first_e, s);
        const bool last_partial = mid_edge_endpoint(off + last_e, s);
        const bool low_partial = (s == LEFT) ? first_partial : last_partial;
        const bool high_partial = (s == LEFT) ? last_partial : first_partial;

        std::size_t mlo = lo + (low_partial ? 1 : 0);
        std::size_t mhi_p1 = hi + 1 - (high_partial ? 1 : 0);

        auto boundary_piece = [&](std::size_t e) {
            ArcPiece p;
            p.subarc = Subarc{e, s, e, s};
            p.is_boundary_piece = true;
            return p;
        };
        auto middle_piece = [&](std::size_t a, std::size_t b) {
            std::vector<Point> vs;
            for (std::size_t v = a; v <= b + 1; ++v)
                vs.push_back(Ci->vertex(v));
            curves.push_back(std::make_unique<Polygon>(std::move(vs)));
            submaps.push_back(std::make_unique<Submap>(
                build_piece_submap(*curves.back())));
            ArcPiece p;
            p.subarc = (s == LEFT) ? Subarc{a, s, b, s}
                                   : Subarc{b, s, a, s};
            p.curve = curves.back().get();
            p.submap = submaps.back().get();
            return p;
        };

        std::vector<std::size_t> chunk_lo, chunk_hi;
        for (std::size_t a = mlo; a < mhi_p1; a += piece_len) {
            chunk_lo.push_back(a);
            chunk_hi.push_back(std::min(a + piece_len - 1, mhi_p1 - 1));
        }

        std::vector<ArcPiece> out;
        if (lo == hi && low_partial && high_partial) {
            // Single-edge target partial at both ends: ONE boundary
            // piece covers it ([C91 §3.0(ii)(3) tex 170]).
            out.push_back(boundary_piece(lo));
            return out;
        }
        if (s == LEFT) {
            if (low_partial) out.push_back(boundary_piece(lo));
            for (std::size_t k = 0; k < chunk_lo.size(); ++k)
                out.push_back(middle_piece(chunk_lo[k], chunk_hi[k]));
            if (high_partial) out.push_back(boundary_piece(hi));
        } else {
            if (high_partial) out.push_back(boundary_piece(hi));
            for (std::size_t k = chunk_lo.size(); k-- > 0; )
                out.push_back(middle_piece(chunk_lo[k], chunk_hi[k]));
            if (low_partial) out.push_back(boundary_piece(lo));
        }
        return out;
    }
};

// Bundle: fixture + oracles wired into a ConformalityOracles.
struct CombRig {
    CombFixture fx;
    const Submap* S_live;
    GeomRayShooter ray1, ray2;
    GeomArcCutter cut1, cut2;
    ConformalityOracles oracles;

    explicit CombRig(std::size_t piece_len = 1)
        : S_live(&fx.S),
          ray1(&fx.C1), ray2(&fx.C2),
          cut1(&fx.C1, &fx.C, &S_live, 0, piece_len),
          cut2(&fx.C2, &fx.C, &S_live, fx.C1.num_edges(), piece_len) {
        oracles.S1 = &fx.S1;
        oracles.S2 = &fx.S2;
        oracles.C1 = &fx.C1;
        oracles.C2 = &fx.C2;
        oracles.ray1 = &ray1;
        oracles.ray2 = &ray2;
        oracles.cut1 = &cut1;
        oracles.cut2 = &cut2;
        oracles.g_gamma1 = 64;
        oracles.g_gamma2 = 64;
        oracles.h_gamma1 = 64;
        oracles.h_gamma2 = 64;
    }
};

// Live arcs of a region (inventory scan).
static std::vector<std::size_t> region_arcs_of(const Submap& S,
                                               std::size_t region) {
    std::vector<std::size_t> out;
    for (std::size_t ai = 0; ai < S.num_arcs(); ++ai)
        if (!S.arc(ai).dead && S.arc(ai).region_node == region)
            out.push_back(ai);
    return out;
}

// Brute-force geometric validation: the chord's endpoints are mutually
// visible with respect to C — no edge of C crosses its open interior.
static void assert_chord_geometrically_valid(const Polygon& C,
                                             const Chord& c) {
    SymbolicY y = c.symbolic_y();
    double x1 = edge_x_at_y(C, c.left_edge, y);
    double x2 = edge_x_at_y(C, c.right_edge, y);
    if (x1 == x2) return;   // zero-length (duplicate-pair) chord
    double lo = std::min(x1, x2);
    double hi = std::max(x1, x2);
    for (std::size_t e = 0; e < C.num_edges(); ++e) {
        double x;
        if (!GeomRayShooter::crossing_x(C, e, y, &x)) continue;
        assert(!(x > lo && x < hi) &&
               "[C91 §2.1 tex 70]: an added chord's open interior must "
               "not cross C (mutual visibility)");
    }
}

// ════════════════════════════════════════════════════════════════
//  1. Fixture sanity + arc_end_symbolic_y
// ════════════════════════════════════════════════════════════════

static void test_fixture_and_arc_end_y() {
    CombFixture fx;
    fx.S.check_invariants(fx.C);
    fx.S.assert_tree_property();
    assert(!fx.S.is_conformal() && "R_big has degree 7 — not conformal yet");

    // arc_end_symbolic_y: R2 (arc 5) ends at v3e's mid-edge endpoint —
    // the chord's y (5, tag 6), via c3's slot-0 reference.
    SymbolicY e5 = fx.S.arc_end_symbolic_y(5, fx.C);
    assert(symbolic_y_equal(e5, SymbolicY{5.0, 6}));
    // S* (arc 15) ends at cAPEX's vertex endpoint: the apex (26, tag 7)
    // — its LEFT leg, past C's start turnaround ([C91 §2.4 tex 142]).
    SymbolicY e15 = fx.S.arc_end_symbolic_y(15, fx.C);
    assert(symbolic_y_equal(e15, SymbolicY{26.0, 7}));
    // ...and STARTS at c0's mid-edge endpoint (6, tag 2) on its RIGHT leg.
    SymbolicY s15 = fx.S.arc_start_symbolic_y(15, fx.C);
    assert(symbolic_y_equal(s15, SymbolicY{6.0, 2}));
    // The null arc Nz (arc 7) ends at the null chord's y (26, tag 7).
    SymbolicY e7 = fx.S.arc_end_symbolic_y(7, fx.C);
    assert(symbolic_y_equal(e7, SymbolicY{26.0, 7}));

    std::printf("  [PASS] fixture_and_arc_end_y\n");
}

// ════════════════════════════════════════════════════════════════
//  2. compute_arc_provenance — [C91 §3.2 tex 244]
// ════════════════════════════════════════════════════════════════

static void test_arc_provenance() {
    CombFixture fx;
    auto prov = compute_arc_provenance(fx.S, fx.C, fx.S1, fx.C1,
                                       fx.S2, fx.C2);
    assert(prov.size() == 16);
    // Each chordless Sᵢ has ONE closed arc (index 0) covering all of
    // ∂Cᵢ ([C91 §2.4 tex 142/138]); every fused arc maps into it.
    auto expect = [&](std::size_t ai, bool on_c1, std::size_t si_arc) {
        assert(prov[ai].on_c1 == on_c1);
        assert(prov[ai].arc_in_si == si_arc);
    };
    expect(0, false, 0);   // Lz  ⊂ S₂ (edge 7 = C₂'s edge 0)
    expect(1, false, 0);   // E*  ⊂ S₂
    expect(5, false, 0);   // R2  ⊂ S₂
    expect(7, false, 0);   // Nz  ⊂ S₂
    expect(8, true, 0);    // P4b ⊂ S₁
    expect(11, true, 0);   // R3  ⊂ S₁
    expect(15, true, 0);   // S*  ⊂ S₁

    std::printf("  [PASS] arc_provenance\n");
}

// ════════════════════════════════════════════════════════════════
//  3. fused_region_cycle — [C91 §3.2 tex 238]
// ════════════════════════════════════════════════════════════════

static void test_fused_region_cycle() {
    CombFixture fx;

    // R_big: 7 arc-structures ([C91 §2.4 tex 142]: the wrap-spanning
    // arcs E* and S* are single structures) — the cycle is the sorted
    // table arcs directly.
    auto cycle = fused_region_cycle(fx.S, fx.C, 0, region_arcs_of(fx.S, 0));
    assert(cycle.count == 7 &&
           "[C91 §3.2 tex 238]: the fused big region has 7 arcs");
    // Boundary order (by clockwise tour position of each arc's start):
    // E*, Z5, R2, Z3, R3, Z1, S*.
    assert(cycle.arcs[0].arc == 1 && !cycle.arcs[0].is_zero_length);
    assert(cycle.arcs[1].arc == 3 && cycle.arcs[1].is_zero_length);
    assert(cycle.arcs[2].arc == 5 && !cycle.arcs[2].is_zero_length);
    assert(cycle.arcs[3].arc == 9 && cycle.arcs[3].is_zero_length);
    assert(cycle.arcs[4].arc == 11);
    assert(cycle.arcs[5].arc == 13 && cycle.arcs[5].is_zero_length);
    assert(cycle.arcs[6].arc == 15 && !cycle.arcs[6].is_zero_length);

    // Pocket P4 (region 4): two arcs (split by the junction null
    // chord n7 — exactly why [C91 §3.1 tex 224] adds junction chords).
    auto p4 = fused_region_cycle(fx.S, fx.C, 4, region_arcs_of(fx.S, 4));
    assert(p4.count == 2);

    // cAPEX's inner region: one zero-length arc.
    auto az = fused_region_cycle(fx.S, fx.C, 7, region_arcs_of(fx.S, 7));
    assert(az.count == 1 && az.arcs[0].is_zero_length);

    std::printf("  [PASS] fused_region_cycle\n");
}

// ════════════════════════════════════════════════════════════════
//  4. local_shoot_fused — [C91 §3.2 tex 244/246]
// ════════════════════════════════════════════════════════════════

static void test_local_shoot_fused() {
    CombRig rig;
    auto& fx = rig.fx;
    auto prov = compute_arc_provenance(fx.S, fx.C, fx.S1, fx.C1,
                                       fx.S2, fx.C2);
    auto cycle = fused_region_cycle(fx.S, fx.C, 0, region_arcs_of(fx.S, 0));

    FusedShootContext ctx;
    ctx.S = &fx.S; ctx.C = &fx.C; ctx.C1 = &fx.C1; ctx.C2 = &fx.C2;
    ctx.ray1 = &rig.ray1; ctx.ray2 = &rig.ray2;
    ctx.provenance = &prov;

    // (a) Direct: v4's under-west companion shoots LEFT and hits
    // edge 0's east wall at (0.4, 4) — on S*'s RIGHT leg (arc 15).
    RayHit h = local_shoot_fused(Point{8, 4, 4}, SymbolicY{4.0, 4}, LEFT,
                                 cycle, ctx);
    assert(h.hit && !h.wrapped);
    assert(h.hit_arc_idx == 15);
    assert(h.edge == 0 && h.side == RIGHT);
    assert(h.x == 0.4);

    // (b) Direct: v8's under-east companion shoots RIGHT and hits
    // edge 11's west wall at x = 22 + 2·23/24 — on E*'s RIGHT leg (arc 1).
    h = local_shoot_fused(Point{16, 2, 8}, SymbolicY{2.0, 8}, RIGHT,
                          cycle, ctx);
    assert(h.hit && !h.wrapped);
    assert(h.hit_arc_idx == 1);
    assert(h.edge == 11 && h.side == RIGHT);

    // (c) Wrapped ([C91 §2.1 tex 70]): P1's apex west-outside companion
    // has nothing to its west; the ray wraps and strikes edge 11's east
    // wall — on E*'s LEFT leg (arc 1).
    h = local_shoot_fused(Point{2, 20, 1}, SymbolicY{20.0, 1}, LEFT,
                          cycle, ctx);
    assert(h.hit && h.wrapped);
    assert(h.hit_arc_idx == 1);
    assert(h.edge == 11 && h.side == LEFT);

    std::printf("  [PASS] local_shoot_fused\n");
}

// ════════════════════════════════════════════════════════════════
//  5. insert_chord — [C91 §3.2 tex 264]
// ════════════════════════════════════════════════════════════════

static void test_insert_chord() {
    CombFixture fx;

    // Insert v4's under-west chord (0.4,4)–(8,4): p on R3 (arc 11, at
    // the polygon vertex v4), q on S* (arc 15, mid-edge on edge 0's
    // RIGHT leg — splitting a wrap-spanning structure,
    // [C91 §2.4 tex 142]).
    std::size_t flat[] = {1, 3, 5, 9, 11, 13, 15};
    Submap::ChordPointSpec p{11, 3, RIGHT, 8.0};
    Submap::ChordPointSpec q{15, 0, RIGHT, 0.4};
    auto res = fx.S.insert_chord(p, q, SymbolicY{4.0, 4}, 0, flat, 7, fx.C);

    assert(res.chord_idx == 8);
    assert(res.new_region == 9);
    assert(fx.S.num_nodes() == 10);
    assert(fx.S.num_chords() == 9);
    assert(fx.S.num_arcs() == 18);
    fx.S.assert_tree_property();

    const Chord& nc = fx.S.chord(8);
    // Slot order by ascending x: left = q (mid-edge → 2 adj arcs),
    // right = p (at vertex 4 → 1 adj arc, [C91 §2.2 tex 94]).
    assert(nc.left_edge == 0 && nc.left_side == RIGHT);
    assert(nc.left_adj.count == 2);
    assert(nc.left_adj.arcs[0] == 15 && nc.left_adj.arcs[1] == res.q_after_arc);
    assert(nc.right_edge == 3 && nc.right_side == RIGHT);
    assert(nc.right_adj.count == 1 && nc.right_adj.arcs[0] == 11);
    assert(nc.region[0] == 0 && nc.region[1] == 9);

    // Chain: the new region owns R3's after-half, Z1, S*'s before-half.
    assert(fx.S.arc(res.p_after_arc).region_node == 9);
    assert(fx.S.arc(13).region_node == 9);          // Z1
    assert(fx.S.arc(15).region_node == 9);          // S* before-half
    assert(fx.S.arc(res.q_after_arc).region_node == 0);
    assert(fx.S.arc(11).region_node == 0);          // R3 before-half

    // Splits: R3 = [4..3] split AT v4 (edge 3's end vertex) → edge 3
    // belongs entirely to the after-half ([C91 §2.4 tex 133]: the
    // arc-structure's pointers are the edges the arc spans): before
    // [4..4] (ends at v4, ec 1) + after [3..3] (starts at v4, ec 1).
    assert(fx.S.arc(11).last_edge == 4 && fx.S.arc(11).first_edge == 4);
    assert(fx.S.arc(11).edge_count == 1);
    assert(fx.S.arc(res.p_after_arc).first_edge == 3 &&
           fx.S.arc(res.p_after_arc).last_edge == 3);
    assert(fx.S.arc(res.p_after_arc).edge_count == 1);

    // Splitting S* at the mid-edge point on its RIGHT leg leaves the
    // before-half a plain RIGHT arc [0,0] and gives the double-backing
    // to the after-half (R(0)→L(6), ᾱ = [0,6], [C91 §2.4 tex 142]);
    // the start-turn pointer moves with it.
    assert(fx.S.arc(15).first_side == RIGHT &&
           fx.S.arc(15).last_side == RIGHT &&
           fx.S.arc(15).edge_count == 1);
    assert(fx.S.arc(res.q_after_arc).first_side == RIGHT &&
           fx.S.arc(res.q_after_arc).last_side == LEFT &&
           fx.S.arc(res.q_after_arc).edge_count == 8);
    assert(fx.S.start_arc == res.q_after_arc);

    // Adjacency re-pointing: c1 (v1e)'s slot 0 referenced R3 as its
    // ENDING arc — now the after-half; and c1 itself moved to region 9.
    assert(fx.S.chord(1).right_adj.arcs[0] == res.p_after_arc);
    assert(fx.S.chord(1).region[0] == 9);
    // c2 (v3w)'s slot 1 referenced R3 as its STARTING arc — unchanged.
    assert(fx.S.chord(2).left_adj.arcs[1] == 11);
    // c0 (v1w) moved to region 9 (its footprint went with the chain).
    assert(fx.S.chord(0).region[0] == 9);

    // Region cycles after the split: new region has 3 arcs; old region
    // 0 now has 6.
    auto c9 = fused_region_cycle(fx.S, fx.C, 9, region_arcs_of(fx.S, 9));
    assert(c9.count == 3);
    auto c0 = fused_region_cycle(fx.S, fx.C, 0, region_arcs_of(fx.S, 0));
    assert(c0.count == 6);

    std::printf("  [PASS] insert_chord\n");
}

// [C91 §2.1 tex 70/72]: on a curve where the horizontal circle at a
// vertex v's y meets C only at v, v's two opposite-side companions see
// each other THROUGH INFINITY — a genuine V(C) chord of zero geometric
// length (equal endpoint x) between distinct ∂C points, exactly like
// the fixture's cAPEX.  [C91 §3.2 tex 264/267]: when Lemma 3.3's
// guaranteed chord is such a chord, insert_chord must accept it; the
// slot order falls to the clockwise-∂C tie-break used by the fusion
// rebuild's canonicalization.
static void test_insert_chord_equal_x_wrap() {
    // y-monotone curve: every horizontal circle crosses C exactly once.
    Polygon C({{0, 0, 0}, {4, 10, 1}, {2, 20, 2}, {6, 30, 3}});

    Submap S;
    for (int i = 0; i < 3; ++i) S.add_node();   // 0=mid, 1=top, 2=bottom
    S.start_vertex = 0;
    S.end_vertex = 3;

    auto add = [&](std::size_t fe, Side fs, std::size_t le, Side ls,
                   std::size_t region, std::size_t count) {
        Arc a{};
        a.first_edge = fe; a.first_side = fs;
        a.last_edge = le; a.last_side = ls;
        a.region_node = region;
        a.edge_count = count;
        return S.add_arc(a);
    };
    // Chords: c0 = augmented ([C91 §2.2 tex 92]) wrap chord between the
    // two sides of edge 2 at y 25 (both endpoints at x 4); c1 = v1's
    // companion-pair wrap chord ([C91 §2.1 tex 72]: "possibly the same
    // chord for the two companions"), both endpoints at x 4.
    add(1, LEFT, 2, LEFT, 0, 2);      // 0 α1: v1L → mid-e2 (LEFT)
    add(2, LEFT, 2, RIGHT, 1, 1);     // 1 α2: end wrap around v3
    add(2, RIGHT, 1, RIGHT, 0, 2);    // 2 α3: mid-e2 (RIGHT) → v1R
    add(0, RIGHT, 0, LEFT, 2, 1);     // 3 α4: start wrap around v0
    assert(S.end_arc == 1 && S.start_arc == 3);

    auto chord = [&](std::size_t le, Side lsd, Chord::AdjArcs ladj,
                     std::size_t re, Side rsd, Chord::AdjArcs radj,
                     double y, std::size_t tag,
                     std::size_t r0, std::size_t r1) {
        Chord c{};
        c.left_edge = le; c.left_side = lsd; c.left_adj = ladj;
        c.right_edge = re; c.right_side = rsd; c.right_adj = radj;
        c.y = y; c.y_tag = tag;
        c.region[0] = r0; c.region[1] = r1;
        c.is_null_length = false;
        S.add_chord(c);
    };
    chord(2, LEFT, {{0, 1}, 2}, 2, RIGHT, {{1, 2}, 2}, 25.0, 99, 1, 0);
    chord(0, LEFT, {{3}, 1}, 1, RIGHT, {{2}, 1}, 10.0, 1, 2, 0);

    // Insert v2's companion-pair chord: p = v2's LEFT-side companion on
    // α1, q = its RIGHT-side companion on α3 — distinct ∂C points, both
    // at (2, 20).
    std::size_t cyc[] = {0, 2};
    Submap::ChordPointSpec p{0, 1, LEFT, 2.0};
    Submap::ChordPointSpec q{2, 1, RIGHT, 2.0};
    auto res = S.insert_chord(p, q, SymbolicY{20.0, 2}, 0, cyc, 2, C);

    assert(S.num_chords() == 3);
    S.assert_tree_property();

    const Chord& nc = S.chord(res.chord_idx);
    // Equal x ⟹ slot order by clockwise ∂C position: (1,LEFT) < (1,RIGHT).
    assert(nc.left_edge == 1 && nc.left_side == LEFT);
    assert(nc.right_edge == 1 && nc.right_side == RIGHT);
    // Vertex endpoints: one adj arc each ([C91 §2.2 tex 94]), the
    // before-halves.
    assert(nc.left_adj.count == 1 && nc.left_adj.arcs[0] == 0);
    assert(nc.right_adj.count == 1 && nc.right_adj.arcs[0] == 2);
    // The chord wraps ([C91 §2.1 tex 70]).
    assert(chord_runs_through_infinity(C, nc));

    // Vertex splits at v2 put each incident edge entirely on one half
    // ([C91 §2.4 tex 133]).
    assert(S.arc(0).first_edge == 1 && S.arc(0).last_edge == 1 &&
           S.arc(0).edge_count == 1);
    assert(S.arc(res.p_after_arc).first_edge == 2 &&
           S.arc(res.p_after_arc).last_edge == 2 &&
           S.arc(res.p_after_arc).edge_count == 1);
    assert(S.arc(2).first_edge == 2 && S.arc(2).last_edge == 2 &&
           S.arc(2).edge_count == 1);
    assert(S.arc(res.q_after_arc).first_edge == 1 &&
           S.arc(res.q_after_arc).last_edge == 1 &&
           S.arc(res.q_after_arc).edge_count == 1);

    // New region = p's after-half + q's before-half (the e2 sides);
    // region 0 keeps the e1 sides.
    assert(S.arc(res.p_after_arc).region_node == res.new_region);
    assert(S.arc(2).region_node == res.new_region);
    assert(S.arc(0).region_node == 0);
    assert(S.arc(res.q_after_arc).region_node == 0);

    // c0's footprint moved with the chain; its ending-arc slot was
    // repointed to the after-half, its starting-arc slot kept.
    assert(S.chord(0).region[0] == 1 &&
           S.chord(0).region[1] == res.new_region);
    assert(S.chord(0).left_adj.arcs[0] == res.p_after_arc);
    assert(S.chord(0).right_adj.arcs[1] == 2);
    // c1's ending-arc slot at v1R now ends at q's after-half.
    assert(S.chord(1).right_adj.arcs[0] == res.q_after_arc);
    assert(S.chord(1).region[0] == 2 && S.chord(1).region[1] == 0);

    std::printf("  [PASS] insert_chord_equal_x_wrap\n");
}

// ════════════════════════════════════════════════════════════════
//  6. find_visible_point — [C91 §3.2 Lemma 3.2]
// ════════════════════════════════════════════════════════════════

static void test_find_visible_point() {
    CombRig rig;
    auto& fx = rig.fx;
    auto prov = compute_arc_provenance(fx.S, fx.C, fx.S1, fx.C1,
                                       fx.S2, fx.C2);
    auto cycle = fused_region_cycle(fx.S, fx.C, 0, region_arcs_of(fx.S, 0));
    // cycle: 0=E*, 1=Z5, 2=R2, 3=Z3, 4=R3, 5=Z1, 6=S*

    // (a) R2 → E*: v8's under-east companion sees edge 11's west
    // wall (nonconsecutive: indices 2 and 0).  Found via a BOUNDARY
    // piece's vertex ([C91 §3.2 tex 246/248] — R2 spans two partial
    // edges, so the cut has no vertex-to-vertex middle).
    VisiblePoint vp = find_visible_point(fx.S, fx.C, 0, cycle.arcs[2].arc,
                                         cycle.arcs[0].arc, cycle, prov,
                                         rig.oracles);
    assert(vp.found);
    assert(symbolic_y_equal(vp.y, SymbolicY{2.0, 8}));
    assert(vp.p_table_arc == 5 && vp.p_edge == 8 && vp.p_side == RIGHT);
    assert(vp.q_table_arc == 1 && vp.q_edge == 11 && vp.q_side == RIGHT);

    // (b) R3 → S*: v4's under-west companion sees edge 0's east
    // wall on S*'s RIGHT leg (indices 4 and 6).
    vp = find_visible_point(fx.S, fx.C, 0, cycle.arcs[4].arc,
                            cycle.arcs[6].arc, cycle, prov, rig.oracles);
    assert(vp.found);
    assert(symbolic_y_equal(vp.y, SymbolicY{4.0, 4}));
    assert(vp.q_table_arc == 15 && vp.q_x == 0.4);

    // (c) S* → R2: no vertex of the outer arc sees the under-comb
    // strip between the teeth — Lemma 3.2 must come back empty (sound
    // failure; this exercises the arc-cutter's split of the
    // wrap-spanning target at C's start turnaround ([C91 §3.0(ii)(2)
    // tex 170]) and the full vertex scan over edges 0..6).
    vp = find_visible_point(fx.S, fx.C, 0, cycle.arcs[6].arc,
                            cycle.arcs[2].arc, cycle, prov, rig.oracles);
    assert(!vp.found);

    // (d) Zero-length arcs are never candidates ([C91 §2.1 tex 70/72]).
    assert(assert_fires([&]{
        find_visible_point(fx.S, fx.C, 0, cycle.arcs[1].arc,
                           cycle.arcs[4].arc, cycle, prov, rig.oracles);
    }));

    std::printf("  [PASS] find_visible_point\n");
}

// ════════════════════════════════════════════════════════════════
//  7. find_visible_point — tree-decomposition descent (Lemma 3.2's
//     binary search + the [C91 §2.5 Lemma 2.4] shielding test)
// ════════════════════════════════════════════════════════════════

// Fixture cutter for the wrap-spanning arc S* ([C91 §2.4 tex 142]):
// the target R(0)→L(6) is split at C's start turnaround
// ([C91 §3.0(ii)(2) tex 170]) into a single-edge boundary piece (the
// RIGHT leg, partial at v1w's mid-edge endpoint) and ONE
// vertex-to-vertex piece [0..6] on the LEFT whose S_α has a real V(ᾱ)
// chord (P2's apex east chord (6,24)–(13.81,24)), forcing the descent
// through an internal TD node.
struct DescentCutter : ArcCuttingOracle {
    const Polygon* C1;
    mutable std::vector<std::unique_ptr<Polygon>> curves;
    mutable std::vector<std::unique_ptr<Submap>> submaps;
    explicit DescentCutter(const Polygon* c1) : C1(c1) {}

    std::vector<ArcPiece> cut(std::size_t /*arc_idx*/,
                              const Subarc& target) const override {
        std::vector<ArcPiece> out;
        assert(target.first_edge == 0 && target.first_side == RIGHT &&
               target.last_edge == 6 && target.last_side == LEFT &&
               "descent test cuts the wrap-spanning S*");
        // RIGHT leg: single edge 0, partial at its traversal start
        // (v1w's mid-edge endpoint) → one boundary piece.
        {
            ArcPiece p;
            p.subarc = Subarc{0, RIGHT, 0, RIGHT};
            p.is_boundary_piece = true;
            out.push_back(p);
        }
        // LEFT leg: vertex-to-vertex [0..6] with a 2-region S_α.
        std::vector<Point> vs;
        for (std::size_t v = 0; v <= 7; ++v) vs.push_back(C1->vertex(v));
        curves.push_back(std::make_unique<Polygon>(std::move(vs)));

        auto sm = std::make_unique<Submap>();
        std::size_t r_out = sm->add_node();
        std::size_t r_pocket = sm->add_node();
        sm->start_vertex = 0;
        sm->end_vertex = 7;
        // [C91 §2.4 tex 142]: r_out's single arc runs from the chord's
        // mid-edge endpoint through BOTH of ᾱ's endpoint turnarounds
        // back to the apex — one DOUBLE-WRAP structure W_out
        // (first=(6,LEFT) > last=(2,LEFT); ᾱ-range = everything).
        Arc a{};
        a.first_edge = 3; a.last_edge = 6;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r_pocket; a.edge_count = 4;
        std::size_t l_pocket = sm->add_arc(a);        // apex → chord end
        a = {};
        a.first_edge = 6; a.last_edge = 2;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r_out; a.edge_count = 7;
        std::size_t w_out = sm->add_arc(a);           // chord end → apex

        // P2's apex east chord of V(ᾱ): (6,24) → (13.81, 24) on edge 6.
        Chord c{};
        c.region[0] = r_out; c.region[1] = r_pocket;
        c.left_edge = 3; c.left_side = LEFT;          // at the apex vertex
        c.right_edge = 6; c.right_side = LEFT;        // mid-edge
        c.y = 24.0; c.y_tag = 3;
        c.left_adj = {{w_out}, 1};
        c.right_adj = {{l_pocket, w_out}, 2};
        sm->add_chord(c);
        assert(sm->start_arc == w_out && sm->end_arc == w_out &&
               "[C91 §2.4 tex 142]: the double-wrap arc is both "
               "endpoint arcs");
        sm->build_tree_decomposition();
        submaps.push_back(std::move(sm));

        ArcPiece p;
        p.subarc = Subarc{0, LEFT, 6, LEFT};
        p.curve = curves.back().get();
        p.submap = submaps.back().get();
        out.push_back(p);
        return out;
    }
};

static void test_descent() {
    CombRig rig;   // rig's cut1 is unused here; we substitute DescentCutter
    auto& fx = rig.fx;
    DescentCutter dcut(&fx.C1);
    ConformalityOracles oracles = rig.oracles;
    oracles.cut1 = &dcut;
    oracles.g_gamma1 = 4;

    auto prov = compute_arc_provenance(fx.S, fx.C, fx.S1, fx.C1,
                                       fx.S2, fx.C2);
    auto cycle = fused_region_cycle(fx.S, fx.C, 0, region_arcs_of(fx.S, 0));

    // S* → R2: the descent enters S_α's tree at the apex chord,
    // rejects the pocket side (its α-portion cannot see A₂ — the
    // [C91 §2.5 Lemma 2.4] machinery), reaches the r_out leaf, scans its
    // vertices, and correctly finds nothing.
    VisiblePoint vp = find_visible_point(fx.S, fx.C, 0, cycle.arcs[6].arc,
                                         cycle.arcs[2].arc, cycle, prov,
                                         oracles);
    assert(!vp.found);

    // S* → R3: v2 (the valley between P1 and P2) is in the POCKET
    // side... it is not on r_out's leaf.  But R3 is seen from v4's own
    // arc, not from S*; the outer arc's vertices cannot see the
    // under-comb strip here either.
    vp = find_visible_point(fx.S, fx.C, 0, cycle.arcs[6].arc,
                            cycle.arcs[4].arc, cycle, prov, oracles);
    assert(!vp.found);

    std::printf("  [PASS] descent\n");
}

// ════════════════════════════════════════════════════════════════
//  8. restore_conformality — [C91 §3.2 Lemma 3.4] end-to-end
// ════════════════════════════════════════════════════════════════

static void run_restore_and_check(std::size_t piece_len) {
    CombRig rig(piece_len);
    auto& fx = rig.fx;
    const std::size_t chords_before = fx.S.num_chords();
    const std::size_t nodes_before = fx.S.num_nodes();

    restore_conformality(fx.S, fx.C, rig.oracles);

    // [C91 §2.3 tex 114]: conformal; tree property maintained.
    assert(fx.S.is_conformal());
    fx.S.assert_tree_property();

    // [C91 §3.2 tex 264]: every region has ≤ 4 (paper) arcs.
    for (std::size_t r = 0; r < fx.S.num_nodes(); ++r) {
        if (fx.S.node(r).dead) continue;
        auto inv = region_arcs_of(fx.S, r);
        if (inv.empty()) continue;
        auto cyc = fused_region_cycle(fx.S, fx.C, r, inv);
        assert(cyc.count <= 4);
    }

    // The 7-arc region needs ≥ 2 chords (7 → k₁+k₂ = 9, both ≥ 3); each
    // insertion splits one region.
    std::size_t added = fx.S.num_chords() - chords_before;
    assert(added >= 2 && added <= 3);
    assert(fx.S.num_nodes() - nodes_before == added);

    // Every added chord is a genuine visibility chord of C.
    for (std::size_t ci = chords_before; ci < fx.S.num_chords(); ++ci) {
        const Chord& c = fx.S.chord(ci);
        assert(!c.dead && !c.is_null_length);
        assert_chord_geometrically_valid(fx.C, c);
        // [C91 §2 tex 47]: sourced at a vertex of C.
        assert(c.y_tag < fx.C.num_vertices());
    }

    // Old chords untouched (restoration only ADDS, [C91 §3.2 tex 264]).
    for (std::size_t ci = 0; ci < chords_before; ++ci)
        assert(!fx.S.chord(ci).dead);
}

static void test_restore_conformality_e2e() {
    run_restore_and_check(1);    // single-edge pieces (leaf-only search)
    std::printf("  [PASS] restore_conformality_e2e (piece_len=1)\n");
    run_restore_and_check(64);   // whole-unit pieces (chordless S_α)
    std::printf("  [PASS] restore_conformality_e2e (piece_len=64)\n");
}

// ════════════════════════════════════════════════════════════════
//  9. Death tests
// ════════════════════════════════════════════════════════════════

static void test_insert_chord_asserts() {
    // (a) A chord endpoint coinciding with an existing chord endpoint
    // must fire ([C91 §2.1 tex 70]: visibility already realized).
    assert(assert_fires([]{
        CombFixture fx;
        std::size_t flat[] = {1, 3, 5, 9, 11, 13, 15};
        // p at v1w's existing mid-edge endpoint on (edge 0, RIGHT, y=6).
        Submap::ChordPointSpec p{15, 0, RIGHT, 0.6};
        Submap::ChordPointSpec q{11, 3, RIGHT, 7.0};
        fx.S.insert_chord(p, q, SymbolicY{6.0, 2}, 0, flat, 7, fx.C);
    }));

    // (b) p and q on the same arc must fire (Lemma 3.3's chords connect
    // nonconsecutive — hence distinct — arcs).
    assert(assert_fires([]{
        CombFixture fx;
        std::size_t flat[] = {1, 3, 5, 9, 11, 13, 15};
        Submap::ChordPointSpec p{11, 4, RIGHT, 8.2};
        Submap::ChordPointSpec q{11, 3, RIGHT, 7.0};
        fx.S.insert_chord(p, q, SymbolicY{4.7, 4}, 0, flat, 7, fx.C);
    }));

    std::printf("  [PASS] insert_chord_asserts\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::setbuf(stdout, nullptr);
    std::printf("[C91 §3.2 tests]:\n");
    test_fixture_and_arc_end_y();
    test_arc_provenance();
    test_fused_region_cycle();
    test_local_shoot_fused();
    test_insert_chord();
    test_insert_chord_equal_x_wrap();
    test_find_visible_point();
    test_descent();
    test_restore_conformality_e2e();
    test_insert_chord_asserts();
    std::printf("All §3.2 tests passed.\n");
    return 0;
}
