// tests/s3.4-tests.cpp — Tests for [C91 §3.4]: Implementing the Oracles.
//
// Covers: the [LT79] planar separator ((i)–(iii) of [C91 §3.4 tex
// 297–302] plus the δ = 2/3 decomposition of tex 304), the Lemma 3.6
// ray-shooting structure (S* faces, dual graph G, vertical line,
// queries validated against a brute-force first-contact reference), the
// [C91 §3.0(i)] oracle adapter, and the first true end-to-end merges —
// the comb of tests/s3.2-tests.cpp with its junction at C's global
// y-maximum, which exercises the through-infinity machinery of
// [C91 §2.1 tex 70] end to end (TODO items formerly blocked on §3.4).

#include "merge/conformality.h"
#include "merge/fusion.h"
#include "merge/granularity.h"
#include "merge/merge.h"
#include "merge/ray_shooting.h"
#include "polygon/polygon.h"
#include "submap/submap.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <functional>
#include <memory>
#include <vector>
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

// Deterministic PRNG (fixed seed; no libc rand).
struct Rng {
    unsigned long long s = 0x9e3779b97f4a7c15ull;
    Rng() = default;
    explicit Rng(unsigned long long seed) : s(seed ^ 0x9e3779b97f4a7c15ull) {}
    unsigned long long next() {
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        return s;
    }
    double uniform(double lo, double hi) {
        return lo + (hi - lo) * (double)(next() % 1000000ull) / 1000000.0;
    }
};
}

// ════════════════════════════════════════════════════════════════
//  Comb fixture (tests/s3.2-tests.cpp geometry)
// ════════════════════════════════════════════════════════════════

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

// One region bounded by the single closed arc — all of ∂C, one
// arc-structure stored cut at C's start turnaround
// ([C91 §2.4 tex 142/138]).
static Submap make_chordless_normal(const Polygon& poly) {
    Submap s;
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = poly.num_vertices() - 1;
    Arc a{};
    a.first_edge = 0; a.last_edge = 0;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 0;
    a.edge_count = poly.count_nonnull_edges(0, poly.num_edges() - 1);
    std::size_t ai0 = s.add_arc(a);
    assert(s.start_arc == ai0 && s.end_arc == ai0);
    s.build_tree_decomposition();
    return s;
}

// The fused comb submap of tests/s3.2-tests.cpp (8 chords, 18 arcs).
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

        add(7, LEFT, 7, LEFT, 7, 0);      // 0  Lz
        add(7, LEFT, 11, RIGHT, 0, 5);    // 1  E*  (end wrap, ᾱ=[7,11])
        add(11, RIGHT, 10, RIGHT, 6, 2);  // 2  P6
        add(9, RIGHT, 9, RIGHT, 0, 0);    // 3  Z5
        add(9, RIGHT, 8, RIGHT, 5, 2);    // 4  P5
        add(8, RIGHT, 7, RIGHT, 0, 2);    // 5  R2
        add(7, RIGHT, 7, RIGHT, 4, 1);    // 6  P4a
        add(7, RIGHT, 7, RIGHT, 8, 0);    // 7  Nz
        add(6, RIGHT, 6, RIGHT, 4, 1);    // 8  P4b
        add(5, RIGHT, 5, RIGHT, 0, 0);    // 9  Z3
        add(5, RIGHT, 4, RIGHT, 3, 2);    // 10 P3
        add(4, RIGHT, 3, RIGHT, 0, 2);    // 11 R3
        add(3, RIGHT, 2, RIGHT, 2, 2);    // 12 P2
        add(1, RIGHT, 1, RIGHT, 0, 0);    // 13 Z1
        add(1, RIGHT, 0, RIGHT, 1, 2);    // 14 P1
        add(0, RIGHT, 6, LEFT, 0, 7);     // 15 S*  (start wrap, ᾱ=[0,6])
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
        chord(0, RIGHT, {{14, 15}, 2}, 1, RIGHT, {{13}, 1}, 6.0, 2, 0, 1);
        chord(2, RIGHT, {{12}, 1}, 3, RIGHT, {{11, 12}, 2}, 6.0, 2, 0, 2);
        chord(4, RIGHT, {{10, 11}, 2}, 5, RIGHT, {{9}, 1}, 5.0, 6, 0, 3);
        chord(6, RIGHT, {{8}, 1}, 7, RIGHT, {{5, 6}, 2}, 5.0, 6, 0, 4);
        chord(8, RIGHT, {{4, 5}, 2}, 9, RIGHT, {{3}, 1}, 4.5, 10, 0, 5);
        chord(10, RIGHT, {{2}, 1}, 11, RIGHT, {{1, 2}, 2}, 4.5, 10, 0, 6);
        chord(6, LEFT, {{15}, 1}, 7, LEFT, {{0}, 1}, 26.0, 7, 0, 7);
        chord(7, RIGHT, {{6}, 1}, 7, RIGHT, {{7}, 1}, 26.0, 7, 4, 8, true);
        // [C91 §2.4(iii) tex 138]: add_arc auto-registered the wrap
        // arcs as the endpoint pointers ([C91 §2.4 tex 142]).
        assert(S.start_arc == 15 && S.end_arc == 1);
    }
};

// ════════════════════════════════════════════════════════════════
//  Geometric test oracles (tests/s3.2-tests.cpp conventions)
// ════════════════════════════════════════════════════════════════

// [C91 §3.0(i) tex 169]-style test shooter: nearest crossing with the
// target subarc only, in the (wrapped, distance) order of
// [C91 §2.1 tex 70].
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
                const auto& ed = Ci->edge(e);
                bool asc = symbolic_y_less(
                    symbolic_y_of(Ci->vertex(ed.start_idx)),
                    symbolic_y_of(Ci->vertex(ed.end_idx)));
                Side minus_x = asc ? LEFT : RIGHT;
                Side struck = (dir == RIGHT)
                    ? minus_x : (minus_x == LEFT ? RIGHT : LEFT);
                if (struck != legs[g].side) continue;
                double d = (dir == RIGHT) ? (x - p.x) : (p.x - x);
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

// [C91 §3.0(ii) tex 170] test cutter with no dependence on the live
// fused submap: BOTH end edges become single-edge boundary pieces
// (allowed at the first/last position whether or not α' stops
// mid-edge), and the interior is cut into vertex-to-vertex chunks of at
// most `piece_len` edges, each with a chordless normal-form submap of
// V(ᾱⱼ) (conformal and h-granular for h ≥ |ᾱⱼ|).
struct SafeArcCutter : ArcCuttingOracle {
    const Polygon* Ci;
    std::size_t piece_len;
    mutable std::vector<std::unique_ptr<Polygon>> curves;
    mutable std::vector<std::unique_ptr<Submap>> submaps;

    SafeArcCutter(const Polygon* ci, std::size_t pl)
        : Ci(ci), piece_len(pl) {}

    // One single-side leg; boundary pieces only at the ends flagged as
    // the TARGET's own endpoints ([C91 §3.0(ii)(3) tex 170]: boundary
    // pieces attach to the endpoints of α'; a leg end at one of Cᵢ's
    // turnarounds is a vertex end and becomes a middle piece).
    std::vector<ArcPiece> cut_leg(const Subarc& leg,
                                  bool target_first,
                                  bool target_last) const {
        const Side s = leg.first_side;
        const std::size_t lo = std::min(leg.first_edge, leg.last_edge);
        const std::size_t hi = std::max(leg.first_edge, leg.last_edge);
        const bool blo = (s == LEFT) ? target_first : target_last;
        const bool bhi = (s == LEFT) ? target_last : target_first;

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
                make_chordless_normal(*curves.back())));
            ArcPiece p;
            p.subarc = (s == LEFT) ? Subarc{a, s, b, s}
                                   : Subarc{b, s, a, s};
            p.curve = curves.back().get();
            p.submap = submaps.back().get();
            return p;
        };

        std::vector<ArcPiece> out;
        if (lo == hi && (blo || bhi)) {
            out.push_back(boundary_piece(lo));
            return out;
        }
        std::size_t mlo = lo + (blo ? 1 : 0);
        std::size_t mhi = hi - (bhi ? 1 : 0);
        std::vector<std::size_t> chunk_lo, chunk_hi;
        for (std::size_t a = mlo; a <= mhi; a += piece_len) {
            chunk_lo.push_back(a);
            chunk_hi.push_back(std::min(a + piece_len - 1, mhi));
        }
        if (s == LEFT) {
            if (blo) out.push_back(boundary_piece(lo));
            for (std::size_t k = 0; k < chunk_lo.size(); ++k)
                out.push_back(middle_piece(chunk_lo[k], chunk_hi[k]));
            if (bhi) out.push_back(boundary_piece(hi));
        } else {
            if (bhi) out.push_back(boundary_piece(hi));
            for (std::size_t k = chunk_lo.size(); k-- > 0; )
                out.push_back(middle_piece(chunk_lo[k], chunk_hi[k]));
            if (blo) out.push_back(boundary_piece(lo));
        }
        return out;
    }

    std::vector<ArcPiece> cut(std::size_t /*arc_idx*/,
                              const Subarc& target) const override {
        // [C91 §3.0(ii)(2) tex 170]: pieces must not double-back — a
        // wrap-spanning target is split at Cᵢ's endpoint turnaround(s)
        // into its single-side legs ([C91 §2.4 tex 142]); the turnaround
        // ends are vertex ends, so only the target's OWN endpoints may
        // produce boundary pieces.
        ArcLeg legs[3];
        std::size_t nl = subarc_legs(target, 0, Ci->num_vertices() - 1,
                                     legs);
        std::vector<ArcPiece> out;
        for (std::size_t g = 0; g < nl; ++g) {
            Subarc leg_sub = (legs[g].side == LEFT)
                ? Subarc{legs[g].lo, LEFT, legs[g].hi, LEFT}
                : Subarc{legs[g].hi, RIGHT, legs[g].lo, RIGHT};
            auto pieces = cut_leg(leg_sub, g == 0, g + 1 == nl);
            out.insert(out.end(), pieces.begin(), pieces.end());
        }
        return out;
    }
};

// ════════════════════════════════════════════════════════════════
//  Brute-force ray reference — first contact with ∂C, any side,
//  lexicographic (wrapped, distance) order ([C91 §2.1 tex 70])
// ════════════════════════════════════════════════════════════════

static RayHit brute_shoot(const Polygon& C, Point p, Side dir) {
    SymbolicY sy{p.y, p.index};
    RayHit best;
    double best_d = 0.0;
    for (std::size_t e = 0; e < C.num_edges(); ++e) {
        double x;
        if (!GeomRayShooter::crossing_x(C, e, sy, &x)) continue;
        const auto& ed = C.edge(e);
        bool asc = symbolic_y_less(symbolic_y_of(C.vertex(ed.start_idx)),
                                   symbolic_y_of(C.vertex(ed.end_idx)));
        Side minus_x = asc ? LEFT : RIGHT;
        Side struck = (dir == RIGHT) ? minus_x
                                     : (minus_x == LEFT ? RIGHT : LEFT);
        double d = (dir == RIGHT) ? (x - p.x) : (p.x - x);
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
    return best;
}

static void check_against_brute(const RayShootingStructure& rs,
                                const Polygon& C, Point p, Side dir) {
    RayHit got = rs.shoot_toward_boundary(p, dir);
    RayHit want = brute_shoot(C, p, dir);
    assert(got.hit == want.hit &&
           "[C91 Lemma 3.6 tex 310–312]: structure and brute force must "
           "agree on hit existence");
    if (!got.hit) return;
    assert(got.wrapped == want.wrapped && got.x == want.x &&
           "[C91 Lemma 3.6 tex 310–312]: structure must report the first "
           "contact in the wrap metric");
    // Coincident contacts (several edges meeting at the struck point)
    // may resolve to different edge names; anything else must match.
    if (got.edge != want.edge) {
        double xa, xb;
        assert(GeomRayShooter::crossing_x(C, got.edge, SymbolicY{p.y, p.index}, &xa) &&
               GeomRayShooter::crossing_x(C, want.edge, SymbolicY{p.y, p.index}, &xb) &&
               xa == xb && "differing edges only at a coincident contact");
    } else {
        assert(got.side == want.side);
    }
}

// ════════════════════════════════════════════════════════════════
//  Conformal γ-granular comb: run the real §3.2 + §3.3 pipeline on the
//  fused fixture, then normalize — a genuine normal-form input for the
//  §3.4 structure.
// ════════════════════════════════════════════════════════════════

struct GranularComb {
    CombFixture fx;
    GeomRayShooter ray1, ray2;
    SafeArcCutter cut1, cut2;
    std::size_t gamma = 0;

    GranularComb() : ray1(&fx.C1), ray2(&fx.C2),
                     cut1(&fx.C1, 64), cut2(&fx.C2, 64) {
        ConformalityOracles o;
        o.S1 = &fx.S1; o.S2 = &fx.S2;
        o.C1 = &fx.C1; o.C2 = &fx.C2;
        o.ray1 = &ray1; o.ray2 = &ray2;
        o.cut1 = &cut1; o.cut2 = &cut2;
        o.g_gamma1 = 64; o.g_gamma2 = 64;
        o.h_gamma1 = 64; o.h_gamma2 = 64;
        restore_conformality(fx.S, fx.C, o);

        // γ = the maximum region weight after conformalization ⟹ the
        // submap is γ-semigranular, [C91 §3.3 tex 276]'s precondition.
        for (std::size_t r = 0; r < fx.S.num_nodes(); ++r)
            if (!fx.S.node(r).dead)
                gamma = std::max(gamma, fx.S.region_weight(r));
        enforce_granularity(fx.S, fx.C, gamma);
        fx.S.normalize(fx.C);
    }
};

// ════════════════════════════════════════════════════════════════
//  Conformal-but-not-granular comb: stop after §3.2 + normalize.  A
//  normal-form conformal γ-semigranular submap (γ = max weight) with
//  MANY more regions than the granular collapse — enough that the
//  §3.4 dual graph G is large enough for the [LT79] separator to
//  RECURSE (nonempty D*, multiple D_i subsets), so the ray-shooting
//  queries genuinely route through the separator decomposition.
// ════════════════════════════════════════════════════════════════

struct ConformalComb {
    CombFixture fx;
    GeomRayShooter ray1, ray2;
    SafeArcCutter cut1, cut2;
    std::size_t gamma = 0;

    ConformalComb() : ray1(&fx.C1), ray2(&fx.C2),
                      cut1(&fx.C1, 64), cut2(&fx.C2, 64) {
        ConformalityOracles o;
        o.S1 = &fx.S1; o.S2 = &fx.S2;
        o.C1 = &fx.C1; o.C2 = &fx.C2;
        o.ray1 = &ray1; o.ray2 = &ray2;
        o.cut1 = &cut1; o.cut2 = &cut2;
        o.g_gamma1 = 64; o.g_gamma2 = 64;
        o.h_gamma1 = 64; o.h_gamma2 = 64;
        restore_conformality(fx.S, fx.C, o);
        for (std::size_t r = 0; r < fx.S.num_nodes(); ++r)
            if (!fx.S.node(r).dead)
                gamma = std::max(gamma, fx.S.region_weight(r));
        // [C91 §3.4 tex 286] takes a conformal γ-granular submap; the
        // ray-shooting structure needs only conformal + semigranular
        // (weights ≤ γ, the O(γ) region-scan bound of tex 306), so we
        // skip enforce_granularity to retain the finer subdivision.
        fx.S.normalize(fx.C);
    }
};

// ════════════════════════════════════════════════════════════════
//  1. μ = 1 structure (chordless submaps) — [C91 §3.4 tex 297]
// ════════════════════════════════════════════════════════════════
//
// The [LT79] planar separator itself ([C91 §3.4 tex 297–304]) is
// covered by tests/lt79-tests.cpp; here we test the Lemma 3.6
// ray-shooting structure that consumes it.

static void test_structure_trivial_mu1() {
    Polygon C = make_C();
    Submap S = make_chordless_normal(C);
    RayShootingStructure rs(S, C, C.num_edges());
    assert(rs.mu() == 1 &&
           "[C91 §3.4 tex 297]: the chordless submap has one face");

    Rng rng;
    for (std::size_t v = 0; v < C.num_vertices(); ++v) {
        check_against_brute(rs, C, C.vertex(v), LEFT);
        check_against_brute(rs, C, C.vertex(v), RIGHT);
    }
    for (int i = 0; i < 200; ++i) {
        Point p{rng.uniform(-3.0, 27.0), rng.uniform(-2.0, 28.0), SOS_NONE};
        check_against_brute(rs, C, p, (i & 1) ? LEFT : RIGHT);
    }
    // Above / below the whole curve: no contact at all.
    RayHit h = rs.shoot_toward_boundary(Point{5, 100, SOS_NONE}, LEFT);
    assert(!h.hit && "[C91 §2.1 tex 70]: no contact above C's y-range");

    std::printf("  [PASS] structure_trivial_mu1\n");
}

// ════════════════════════════════════════════════════════════════
//  3. The γ-granular comb structure — faces, G, vertical line
// ════════════════════════════════════════════════════════════════

static void test_structure_comb() {
    GranularComb gc;
    const Submap& S = gc.fx.S;
    const Polygon& C = gc.fx.C;

    RayShootingStructure rs(S, C, gc.gamma);

    // Faces ↔ nonempty regions ([C91 §3.4 tex 286]): every live region
    // except null-chord interiors gets a face.
    std::size_t live = S.num_live_nodes();
    std::size_t nulls = 0;
    for (std::size_t ci = 0; ci < S.num_chords(); ++ci)
        if (!S.chord(ci).dead && S.chord(ci).is_null_length) ++nulls;
    assert(rs.mu() == live - nulls &&
           "[C91 §3.4 tex 286]: empty regions have no associated faces");
    assert(rs.mu() > 1 && "the granular comb retains structure");

    // The junction outside-pair chord (cAPEX) survives granularity (its
    // contraction weight 12 exceeds γ) and is the ONLY chord through
    // infinity, so the vertical line has exactly one crossing
    // ([C91 §3.4 tex 306]); the region above it is the polar cap.
    assert(rs.vertical_line().size() == 1 &&
           "[C91 §3.4 tex 306]: exactly the apex outside-pair chord "
           "crosses the vertical line");
    {
        const auto& lc = rs.vertical_line()[0];
        assert(symbolic_y_equal(lc.y, SymbolicY{26.0, 7}));
        assert(S.region_weight(lc.region_above) == 0 &&
               "the region above the topmost crossing is the polar cap");
    }

    // Queries vs brute force — vertices (their own SoS tags), random
    // free points, both directions, including wrapped answers.
    Rng rng;
    for (std::size_t v = 0; v < C.num_vertices(); ++v) {
        check_against_brute(rs, C, C.vertex(v), LEFT);
        check_against_brute(rs, C, C.vertex(v), RIGHT);
    }
    for (int i = 0; i < 500; ++i) {
        Point p{rng.uniform(-3.0, 27.0), rng.uniform(-2.0, 28.0), SOS_NONE};
        check_against_brute(rs, C, p, (i & 1) ? LEFT : RIGHT);
    }

    std::printf("  [PASS] structure_comb (mu=%zu, gamma=%zu)\n",
                rs.mu(), gc.gamma);
}

// ════════════════════════════════════════════════════════════════
//  3b. Structure over a larger submap that forces the [LT79]
//      separator to RECURSE — [C91 §3.4 tex 304] within §3.4
// ════════════════════════════════════════════════════════════════

static void test_structure_separator_recursion() {
    ConformalComb cc;
    const Submap& S = cc.fx.S;
    const Polygon& C = cc.fx.C;

    S.check_invariants(C);
    assert(S.is_conformal());
    assert(S.is_semigranular(cc.gamma));

    RayShootingStructure rs(S, C, cc.gamma);

    // The finer subdivision yields μ ≫ 1, so the decomposition
    // ([C91 §3.4 tex 304]) is not a single leaf: the root piece
    // (μ³ > μ²) is split by planar_separator, producing a NONEMPTY D*
    // and MULTIPLE D_i subsets.  This is what makes the query route
    // through the separator (D* scan → double identification →
    // per-D_i scan, tex 306–308) rather than a trivial whole-graph
    // scan.
    assert(rs.mu() > 4 &&
           "the conformal comb keeps many regions (μ > μ^{2/3} leaf "
           "threshold)");
    assert(rs.decomposition().dstar_size >= 1 &&
           "[C91 §3.4 tex 304]: a non-leaf μ forces a nonempty D*");
    assert(rs.decomposition().num_subsets >= 2 &&
           "[C91 §3.4 tex 304]: the separator partitions G into ≥ 2 D_i");

    // Every query still matches the brute-force first contact,
    // including those whose ray-shoot lands in a non-D* region (the
    // tex 308 D_i path) and the wrapped answers past the apex crossing.
    Rng rng(0xC0FFEEu);
    for (std::size_t v = 0; v < C.num_vertices(); ++v) {
        check_against_brute(rs, C, C.vertex(v), LEFT);
        check_against_brute(rs, C, C.vertex(v), RIGHT);
    }
    for (int i = 0; i < 3000; ++i) {
        Point p{rng.uniform(-3.0, 27.0), rng.uniform(-2.0, 28.0), SOS_NONE};
        check_against_brute(rs, C, p, (i & 1) ? LEFT : RIGHT);
    }

    std::printf("  [PASS] structure_separator_recursion "
                "(mu=%zu, |D*|=%zu, subsets=%zu)\n",
                rs.mu(), rs.decomposition().dstar_size,
                rs.decomposition().num_subsets);
}

// ════════════════════════════════════════════════════════════════
//  3c. Wrap-straddling regions: single double-backing arc-structures
//      ([C91 §2.4 tex 142]) — per-region structure count == degree ≤ 4
//      ([C91 §2.3 tex 114]); collect_region_arcs and region_weight
//      reach the wrap arcs through chord adjacency alone
// ════════════════════════════════════════════════════════════════

// Ground truth: a region's live arcs by full-table scan.
static std::vector<std::size_t> region_table_arcs(const Submap& S,
                                                  std::size_t r) {
    std::vector<std::size_t> out;
    for (std::size_t ai = 0; ai < S.num_arcs(); ++ai)
        if (!S.arc(ai).dead && S.arc(ai).region_node == r)
            out.push_back(ai);
    return out;
}

static void test_collect_region_arcs_wrap_straddle() {
    ConformalComb cc;
    const Submap& S = cc.fx.S;

    // [C91 §2.4 tex 142]: the conformal comb keeps regions whose arcs
    // pass around C's endpoint turnarounds — each such arc is ONE
    // double-backing structure, so every region has exactly `degree`
    // arc-structures, at most 4 ([C91 §2.3 tex 114] / [C91 §3.2 tex
    // 264]).
    bool saw_wrap_region = false;
    for (std::size_t r = 0; r < S.num_nodes(); ++r) {
        if (S.node(r).dead) continue;
        auto want = region_table_arcs(S, r);
        if (want.empty()) continue;

        assert(want.size() <= 4 &&
               "[C91 §2.3 tex 114]: conformal region has ≤ 4 "
               "arc-structures — wrap-spanning arcs are never split");
        for (std::size_t ai : want)
            if (S.arc(ai).wraps()) saw_wrap_region = true;

        // collect_region_arcs must gather EVERY arc-structure of the
        // region — wrap arcs are reached as their chords' before-arcs
        // ([C91 §2.2 tex 96] alternation), with no endpoint special
        // cases.
        RegionArcs got = collect_region_arcs(S, r);
        assert(got.count == want.size() &&
               "[C91 §3.1 tex 181]: collect_region_arcs must gather "
               "every arc-structure of the region");
        for (std::size_t ai : want) {
            bool found = false;
            for (std::size_t g : got) if (g == ai) { found = true; break; }
            assert(found &&
                   "[C91 §3.1 tex 181]: collect_region_arcs reaches "
                   "every arc-structure through chord adjacency");
        }
    }
    assert(saw_wrap_region &&
           "[C91 §2.4 tex 142]: the conformal comb must contain a "
           "region with a double-backing arc (else this test proves "
           "nothing about the single-structure representation)");

    std::printf("  [PASS] collect_region_arcs_wrap_straddle\n");
}

static void test_region_weight_wrap_straddle() {
    ConformalComb cc;
    const Submap& S = cc.fx.S;

    // region_weight gathers through chord adjacency ([C91 §2.4(ii)
    // tex 137]); wrap-spanning arcs are single structures whose
    // edge_count covers the union of their legs ([C91 §2.4 tex 142] /
    // [C91 §2.2 tex 106]).  Referee it against the true max edge_count
    // over ALL the region's arc-structures.
    bool saw_wrap = false;
    for (std::size_t r = 0; r < S.num_nodes(); ++r) {
        if (S.node(r).dead) continue;
        auto arcs = region_table_arcs(S, r);
        std::size_t truth = 0;
        for (std::size_t ai : arcs) {
            truth = std::max(truth, S.arc(ai).edge_count);
            if (S.arc(ai).wraps()) saw_wrap = true;
        }
        assert(S.region_weight(r) == truth &&
               "[C91 §2.2 tex 106]: region_weight must equal the max "
               "edge_count over ALL the region's arc-structures — "
               "including double-backing ones");
    }
    assert(saw_wrap &&
           "the conformal comb must contain a double-backing arc");

    std::printf("  [PASS] region_weight_wrap_straddle\n");
}

// ════════════════════════════════════════════════════════════════
//  3d. A merged submap with NO through-infinity chord — exercises
//      build_vertical_line's region_infinity_ branch and the query's
//      tex-308 no-D*-hit vertical-line path ([C91 §3.4 tex 306–308])
// ════════════════════════════════════════════════════════════════

static void test_structure_no_wrapped_chords() {
    // The merged junction (4,5) is a LOCAL max but NOT the global max
    // (9), so its outside pair is a finite chord — no chord of the
    // fusion runs through infinity ([C91 §2.1 tex 70]).  With μ > 1 and
    // an empty wrapped-chord set, build_vertical_line takes the
    // region_infinity_ branch, and queries that strike no D* region
    // route through it (the tex-308 vertical-line path), rather than the
    // vertical-line binary search the comb exercises.
    Polygon C1({{0,0,0}, {2,3,1}, {4,5,2}});
    Polygon C2({{4,5,2}, {6,2,3}, {8,9,4}});
    Submap S1 = make_chordless_normal(C1);
    Submap S2 = make_chordless_normal(C2);
    const std::size_t g = std::max(C1.num_edges(), C2.num_edges());
    SubmapRayShooter r1(S1, C1, g), r2(S2, C2, g);
    SafeArcCutter c1(&C1, 64), c2(&C2, 64);

    MergeInput in;
    in.C1 = &C1; in.C2 = &C2; in.S1 = &S1; in.S2 = &S2;
    in.gamma1 = g; in.gamma2 = g; in.gamma = g;
    in.ray_shooter_1 = &r1; in.ray_shooter_2 = &r2;
    in.arc_cutter_1 = &c1; in.arc_cutter_2 = &c2;
    in.g_gamma1 = 64; in.g_gamma2 = 64; in.h_gamma1 = 64; in.h_gamma2 = 64;
    MergeResult res = merge(in);

    std::size_t gamma = 0;
    for (std::size_t r = 0; r < res.S.num_nodes(); ++r)
        if (!res.S.node(r).dead)
            gamma = std::max(gamma, res.S.region_weight(r));

    RayShootingStructure rs(res.S, res.C, gamma);
    assert(rs.mu() > 1 && "the merge must retain structure");
    assert(rs.vertical_line().empty() &&
           "[C91 §2.1 tex 70]: no chord wraps ⟹ empty vertical line");
    assert(rs.region_at_infinity() != NONE &&
           rs.region_at_infinity() < res.S.num_nodes() &&
           !res.S.node(rs.region_at_infinity()).dead &&
           "[C91 §3.4 tex 306]: a live polar region is elected when no "
           "chord crosses the vertical line");

    Rng rng(0x5EEDu);
    for (std::size_t v = 0; v < res.C.num_vertices(); ++v) {
        check_against_brute(rs, res.C, res.C.vertex(v), LEFT);
        check_against_brute(rs, res.C, res.C.vertex(v), RIGHT);
    }
    for (int i = 0; i < 2000; ++i) {
        Point p{rng.uniform(-2.0, 10.0), rng.uniform(-2.0, 11.0), SOS_NONE};
        check_against_brute(rs, res.C, p, (i & 1) ? LEFT : RIGHT);
    }

    std::printf("  [PASS] structure_no_wrapped_chords (mu=%zu, "
                "region_infinity path)\n", rs.mu());
}

// ════════════════════════════════════════════════════════════════
//  3e. Regression: merge across a LOCAL-MINIMUM junction — the mid-edge
//      current point p must carry its borrowed SoS tag ([C91 §2 tex 47])
// ════════════════════════════════════════════════════════════════

static void test_local_min_junction() {
    // Junction (4,2) is a LOCAL min of C (neighbours 3 and 4) but not the
    // GLOBAL min (v0 = 1).  In the case (ii) startup, c₀ = (1,2) on ∂C₁
    // sees ∂C₂ only at the junction vertex, which lies exactly at c₀'s y.
    // c₀ has no vertex tag of its own and borrows a₀'s (the junction's);
    // if that tag is dropped (Point.index = NONE) the oracle shoots at a
    // perturbed y that slips past the junction, tripping the [C91 §3.1
    // tex 181] "local shoot must hit" invariant.  This merge must
    // complete and the result must be geometrically correct.
    Polygon C1({{0,1,0}, {2,3,1}, {4,2,2}});
    Polygon C2({{4,2,2}, {6,4,3}, {8,7,4}});
    Submap S1 = make_chordless_normal(C1);
    Submap S2 = make_chordless_normal(C2);
    const std::size_t g = std::max(C1.num_edges(), C2.num_edges());
    SubmapRayShooter r1(S1, C1, g), r2(S2, C2, g);
    SafeArcCutter c1(&C1, 64), c2(&C2, 64);

    MergeInput in;
    in.C1 = &C1; in.C2 = &C2; in.S1 = &S1; in.S2 = &S2;
    in.gamma1 = g; in.gamma2 = g; in.gamma = g;
    in.ray_shooter_1 = &r1; in.ray_shooter_2 = &r2;
    in.arc_cutter_1 = &c1; in.arc_cutter_2 = &c2;
    in.g_gamma1 = 64; in.g_gamma2 = 64; in.h_gamma1 = 64; in.h_gamma2 = 64;
    MergeResult res = merge(in);           // must not abort

    res.S.check_invariants(res.C);
    assert(res.S.is_conformal());
    assert(res.S.is_granular(g, res.C));

    std::size_t gamma = 0;
    for (std::size_t r = 0; r < res.S.num_nodes(); ++r)
        if (!res.S.node(r).dead)
            gamma = std::max(gamma, res.S.region_weight(r));
    RayShootingStructure rs(res.S, res.C, gamma);

    Rng rng(0xA11CEu);
    for (std::size_t v = 0; v < res.C.num_vertices(); ++v) {
        check_against_brute(rs, res.C, res.C.vertex(v), LEFT);
        check_against_brute(rs, res.C, res.C.vertex(v), RIGHT);
    }
    for (int i = 0; i < 2000; ++i) {
        Point p{rng.uniform(-2.0, 10.0), rng.uniform(-1.0, 9.0), SOS_NONE};
        check_against_brute(rs, res.C, p, (i & 1) ? LEFT : RIGHT);
    }

    std::printf("  [PASS] local_min_junction\n");
}

// ════════════════════════════════════════════════════════════════
//  4. [C91 §3.0(i)] adapter — subarc filtering
// ════════════════════════════════════════════════════════════════

static void test_oracle_adapter() {
    Polygon C1 = make_C1();
    Submap S1 = make_chordless_normal(C1);
    SubmapRayShooter shooter(S1, C1, C1.num_edges());

    // The chordless S₁ has ONE closed arc (index 0) covering all of
    // ∂C₁ ([C91 §2.4 tex 142/138]); subarc targets select within it.
    // v4 = (8,4): shooting LEFT hits edge 0's east wall at (0.4, 4).
    Point p = C1.vertex(4);
    Subarc right_whole{6, RIGHT, 0, RIGHT};
    RayHit h = shooter.shoot(p, LEFT, 0, right_whole);
    assert(h.hit && !h.wrapped && h.edge == 0 && h.side == RIGHT &&
           h.x == 0.4);

    // Same shot toward the LEFT side only: the first contact is on the
    // RIGHT side, so the report is empty ([C91 §3.0(i) tex 169]: the
    // report concerns α' only).
    Subarc left_whole{0, LEFT, 6, LEFT};
    h = shooter.shoot(p, LEFT, 0, left_whole);
    assert(!h.hit);

    // Restricting to a right-side range that excludes edge 0 also
    // filters the hit out.
    Subarc right_partial{6, RIGHT, 2, RIGHT};
    h = shooter.shoot(p, LEFT, 0, right_partial);
    assert(!h.hit);

    // P1's apex west companion: nothing to its west — the ray wraps
    // ([C91 §2.1 tex 70]) and strikes edge 6's east wall.
    Point apex_w = C1.vertex(1);
    h = shooter.shoot(apex_w, LEFT, 0, right_whole);
    assert(h.hit && h.wrapped && h.edge == 6 && h.side == RIGHT);

    std::printf("  [PASS] oracle_adapter\n");
}

// ════════════════════════════════════════════════════════════════
//  5. End-to-end merge — the comb, junction at C's global maximum
// ════════════════════════════════════════════════════════════════

// The junction (14,26) is the global y-maximum of C₁, so pass 1's a₀ is
// an inside duplicate ([C91 §3.1 tex 191]) and pass 2's startup shot
// wraps through infinity, discovering the apex outside-pair chord —
// exactly the through-infinity fusion scenario formerly blocked on
// [C91 §3.4].

static void test_fusion_wrapped_junction() {
    Polygon C1 = make_C1(), C2 = make_C2();
    Submap S1 = make_chordless_normal(C1);
    Submap S2 = make_chordless_normal(C2);
    SubmapRayShooter ray1(S1, C1, 7), ray2(S2, C2, 5);

    // Pass 2 walks C₂ (junction at its first vertex); its startup shot
    // travels east from the apex outside duplicate, wraps, and comes
    // back to the C₁-side duplicate: the discovered chord is the apex
    // outside pair (6,LEFT)–(7,LEFT) at y = {26, 7} in C's frame.
    FusionState st2;
    st2.junction_at_end = false;
    fuse_submaps(st2, S2, C2, S1, C1, ray2, ray1);
    bool found_apex = false;
    for (const auto& dc : st2.chords) {
        if (!symbolic_y_equal(dc.y, SymbolicY{26.0, 7})) continue;
        // Walker = C₂ (edge 0 ↦ C-edge 7); target = C₁ (edge 6).
        std::size_t le = dc.left_on_walker ? dc.left_edge + 7 : dc.left_edge;
        std::size_t re = dc.right_on_walker ? dc.right_edge + 7
                                            : dc.right_edge;
        if ((le == 6 && re == 7) || (le == 7 && re == 6)) found_apex = true;
    }
    assert(found_apex &&
           "[C91 §2.1 tex 70] + [C91 §3.1 tex 191]: the wrapped startup "
           "discovers the apex outside-pair chord");

    std::printf("  [PASS] fusion_wrapped_junction\n");
}

static MergeResult run_comb_merge(std::size_t gamma) {
    Polygon C1 = make_C1(), C2 = make_C2();
    static std::vector<std::unique_ptr<Polygon>> keep_p;
    static std::vector<std::unique_ptr<Submap>> keep_s;
    keep_p.push_back(std::make_unique<Polygon>(C1));
    keep_p.push_back(std::make_unique<Polygon>(C2));
    keep_s.push_back(std::make_unique<Submap>(
        make_chordless_normal(*keep_p[keep_p.size() - 2])));
    keep_s.push_back(std::make_unique<Submap>(
        make_chordless_normal(*keep_p[keep_p.size() - 1])));
    const Polygon& c1 = *keep_p[keep_p.size() - 2];
    const Polygon& c2 = *keep_p[keep_p.size() - 1];
    Submap& s1 = *keep_s[keep_s.size() - 2];
    Submap& s2 = *keep_s[keep_s.size() - 1];

    static std::vector<std::unique_ptr<SubmapRayShooter>> keep_r;
    static std::vector<std::unique_ptr<SafeArcCutter>> keep_c;
    keep_r.push_back(std::make_unique<SubmapRayShooter>(s1, c1, 7));
    keep_r.push_back(std::make_unique<SubmapRayShooter>(s2, c2, 7));
    keep_c.push_back(std::make_unique<SafeArcCutter>(&c1, 64));
    keep_c.push_back(std::make_unique<SafeArcCutter>(&c2, 64));

    MergeInput in;
    in.C1 = &c1; in.C2 = &c2;
    in.S1 = &s1; in.S2 = &s2;
    in.gamma1 = 7; in.gamma2 = 7; in.gamma = gamma;
    in.ray_shooter_1 = keep_r[keep_r.size() - 2].get();
    in.ray_shooter_2 = keep_r[keep_r.size() - 1].get();
    in.arc_cutter_1 = keep_c[keep_c.size() - 2].get();
    in.arc_cutter_2 = keep_c[keep_c.size() - 1].get();
    in.g_gamma1 = 64; in.g_gamma2 = 64;
    in.h_gamma1 = 64; in.h_gamma2 = 64;
    return merge(in);
}

static void test_merge_comb_e2e() {
    // γ = 12 ≥ the total weight: γ-granularity then means NO chords at
    // all (any surviving chord would be removable — see TODO.md's §3.3
    // deviation note), so the merge must collapse to the chordless
    // normal form of V(C).
    {
        MergeResult r = run_comb_merge(12);
        r.S.check_invariants(r.C);
        assert(r.S.is_conformal());
        assert(r.S.is_granular(12, r.C));
        assert(r.S.num_live_chords() == 0 &&
               "γ ≥ total weight ⟹ γ-granular ⟺ chordless");
        assert(r.S.num_live_nodes() == 1);
    }

    // γ = 7: [C91 Lemma 3.5 tex 279] — merge() asserts conformality,
    // γ-granularity, and normal form internally; re-verify externally
    // and check that real structure survives (the merge is not a
    // trivial collapse).
    {
        MergeResult r = run_comb_merge(7);
        r.S.check_invariants(r.C);
        assert(r.S.is_conformal() &&
               "[C91 Lemma 3.5 tex 279]: merge output is conformal");
        assert(r.S.is_granular(7, r.C) &&
               "[C91 Lemma 3.5 tex 279]: merge output is γ-granular");
        assert(!r.S.tree_decomposition().empty() &&
               "[C91 §2.4(iv)]: normal form includes the tree "
               "decomposition");
        assert(r.S.num_live_chords() >= 1 &&
               "γ = 7 < total weight forces surviving exit chords");
        // Every region weight within γ ([C91 §2.3]).
        for (std::size_t rg = 0; rg < r.S.num_nodes(); ++rg)
            if (!r.S.node(rg).dead)
                assert(r.S.region_weight(rg) <= 7);
    }

    std::printf("  [PASS] merge_comb_e2e\n");
}

// ════════════════════════════════════════════════════════════════
//  6. Death tests — structure preconditions
// ════════════════════════════════════════════════════════════════

static void test_structure_preconditions() {
    // Non-conformal input must fire ([C91 §3.4 tex 286]).
    assert(assert_fires([]{
        CombFixture fx;               // R_big has degree 7
        RayShootingStructure rs(fx.S, fx.C, 100);
        (void)rs;
    }));
    // γ below a region weight must fire (semigranularity, tex 286/306).
    assert(assert_fires([]{
        Polygon C = make_C();
        Submap S = make_chordless_normal(C);
        RayShootingStructure rs(S, C, 1);
        (void)rs;
    }));
    std::printf("  [PASS] structure_preconditions\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::setbuf(stdout, nullptr);
    std::printf("[C91 §3.4 tests]:\n");
    test_structure_trivial_mu1();
    test_structure_comb();
    test_structure_separator_recursion();
    test_collect_region_arcs_wrap_straddle();
    test_region_weight_wrap_straddle();
    test_structure_no_wrapped_chords();
    test_local_min_junction();
    test_oracle_adapter();
    test_fusion_wrapped_junction();
    test_merge_comb_e2e();
    test_structure_preconditions();
    std::printf("All §3.4 tests passed.\n");
    return 0;
}
