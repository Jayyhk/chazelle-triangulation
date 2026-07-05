// tests/s3.3-tests.cpp — Tests for [C91 §3.3]: Maintaining Granularity.
//
// Covers: the [C91 §2.3 tex 121] criterion-(ii) enforcement pass
// (tex 276), remove_chord's [C91 §2.2 tex 96/108] junction gluing at
// vertex endpoints and across null-length chords, Submap::normalize
// (tex 276 "We can now put S in normal form"), and the full
// §3.2 → §3.3 pipeline on the six-peak comb with REAL geometric
// oracles.

#include "merge/conformality.h"
#include "merge/granularity.h"
#include "polygon/polygon.h"
#include "submap/submap.h"

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

// Brute-force region weight: max nonnull-edge count over the region's
// live arcs ([C91 §2.2 tex 106]), by full table scan.
std::size_t brute_region_weight(const Submap& s, std::size_t region) {
    std::size_t w = 0;
    for (std::size_t ai = 0; ai < s.num_arcs(); ++ai) {
        const Arc& a = s.arc(ai);
        if (a.dead || a.region_node != region) continue;
        if (a.edge_count > w) w = a.edge_count;
    }
    return w;
}
}

// ════════════════════════════════════════════════════════════════
//  The six-peak comb fixture ([C91 §3.2] test geometry) + oracles.
//  Mirrors tests/s3.2-tests.cpp — see there for the full derivation.
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
        add(7, LEFT, 11, RIGHT, 0, 5);    // 1  E*  (end wrap, ᾱ=[7,11])
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
        // arcs as the endpoint pointers ([C91 §2.4 tex 142]).
        assert(S.start_arc == 15 && S.end_arc == 1);
    }
};

// True geometric ray shooter over Cᵢ (see s3.2-tests.cpp).
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

// Geometric arc cutter over Cᵢ (see s3.2-tests.cpp).
struct GeomArcCutter : ArcCuttingOracle {
    const Polygon* Ci;
    const Polygon* C;
    const Submap* const* S_ptr;
    std::size_t off;
    std::size_t piece_len;
    mutable std::vector<std::unique_ptr<Polygon>> curves;
    mutable std::vector<std::unique_ptr<Submap>> submaps;

    GeomArcCutter(const Polygon* ci, const Polygon* c,
                  const Submap* const* sp, std::size_t o, std::size_t pl)
        : Ci(ci), C(c), S_ptr(sp), off(o), piece_len(pl) {}

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
                make_chordless_normal(*curves.back())));
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

// ════════════════════════════════════════════════════════════════
//  Simple 3-region fixture (see s2.2-tests) for deterministic
//  enforcement checks: r0 —c0(v1)— r1 —c1(v2)— r2 on a 5-vertex chain.
// ════════════════════════════════════════════════════════════════

static Polygon chain_polygon() {
    return Polygon({{0,0,0}, {1,1,1}, {2,2,2}, {3,3,3}, {4,4,4}});
}

// [C91 §2.4 tex 142]: r2's arc double-backs around C's end vertex (E:
// L(2)→R(2), ᾱ=[2,3]) and r0's around its start vertex (S: R(0)→L(0),
// ᾱ=[0,0]) — single structures, never split.
static Submap build_3region_submap() {
    Submap s;
    std::size_t r0 = s.add_node();
    std::size_t r1 = s.add_node();
    std::size_t r2 = s.add_node();
    s.start_vertex = 0;
    s.end_vertex = 4;

    Arc a{};
    auto add = [&](std::size_t fe, Side fs, std::size_t le, Side ls,
                   std::size_t region, std::size_t count) {
        a = {};
        a.first_edge = fe; a.last_edge = le;
        a.first_side = fs; a.last_side = ls;
        a.region_node = region; a.edge_count = count;
        return s.add_arc(a);
    };
    std::size_t a1 = add(1, LEFT, 1, LEFT, r1, 1);
    std::size_t aE = add(2, LEFT, 2, RIGHT, r2, 2);   // end wrap
    std::size_t a4 = add(1, RIGHT, 1, RIGHT, r1, 1);
    std::size_t aS = add(0, RIGHT, 0, LEFT, r0, 1);   // start wrap

    // Vertex-endpoint chords: ONE adj arc per endpoint (before-arc,
    // [C91 §2.2 tex 94 + §2.4(ii)] 1-slot convention).
    Chord c0{};
    c0.region[0] = r0; c0.region[1] = r1;
    c0.left_adj = {{aS}, 1}; c0.right_adj = {{a4}, 1};
    c0.left_edge = 0; c0.right_edge = 0; c0.y = 1.0; c0.y_tag = 1;
    s.add_chord(c0);

    Chord c1{};
    c1.region[0] = r1; c1.region[1] = r2;
    c1.left_adj = {{a1}, 1}; c1.right_adj = {{aE}, 1};
    c1.left_edge = 1; c1.right_edge = 1; c1.y = 2.0; c1.y_tag = 2;
    s.add_chord(c1);

    assert(s.start_arc == aS && s.end_arc == aE &&
           "[C91 §2.4(iii) tex 138]: endpoint arcs auto-registered");
    return s;
}

// ════════════════════════════════════════════════════════════════
//  1. Deterministic enforcement on the 3-region chain
// ════════════════════════════════════════════════════════════════

static void test_enforce_3region_partial() {
    // γ = 2: c0 is removable (deg(r0) = 1 < 3; contraction glues
    // S+a1 and a4+(the merged arc) into ONE start-wrap arc R(1)→L(1)
    // with ᾱ=[0,1] (ec 2 ≤ γ, [C91 §2.4 tex 142]).  After it, c1's
    // contraction would close ∂C (ec 4 > γ) — kept.  [C91 §2.3
    // tex 121] criterion (ii) holds.
    Submap s = build_3region_submap();
    Polygon poly = chain_polygon();
    s.check_invariants(poly);

    enforce_granularity(s, poly, 2);

    assert(s.num_live_nodes() == 2 && s.num_live_chords() == 1);
    assert(!s.chord(1).dead && "c1's contraction weight 4 > γ = 2: kept");

    // Live arcs: the merged start-wrap arc (ᾱ=[0,1], ec 2) and E
    // (ᾱ=[2,3], ec 2) — both single double-backing structures
    // ([C91 §2.4 tex 142]).
    assert(s.num_live_arcs() == 2);
    std::size_t two_edge = 0;
    for (std::size_t ai = 0; ai < s.num_arcs(); ++ai)
        if (!s.arc(ai).dead && s.arc(ai).edge_count == 2) ++two_edge;
    assert(two_edge == 2 && "S+a1+a4 glued into one wrap arc; E untouched");

    // [C91 §3.3 tex 276]: "We can now put S in normal form."
    s.normalize(poly);
    s.check_invariants(poly);
    assert(!s.tree_decomposition().empty());
    assert(s.is_conformal());
    assert(s.is_granular(2, poly) &&
           "[C91 Lemma 3.5]: output must be γ-granular");

    // [C91 §2.2 tex 96/108]: junction gluing keeps every live arc
    // reachable, so region_weight matches brute force for every live
    // region after removals + normalize.
    for (std::size_t r = 0; r < s.num_nodes(); ++r) {
        if (s.node(r).dead) continue;
        assert(s.region_weight(r) == brute_region_weight(s, r));
    }

    // double_identify works on the re-normalized table: v2's y on
    // edge 1 passes through the glued arcs.
    auto res = s.double_identify(1, SymbolicY{2.0, 2}, poly);
    assert(res.count >= 1 && res.count <= 6);

    std::printf("  [PASS] enforce_3region_partial\n");
}

static void test_enforce_3region_full() {
    // γ = 4: both chords removable → single chordless region with the
    // two full-side arcs [0..3] / [3..0] (ec 4 each).
    Submap s = build_3region_submap();
    Polygon poly = chain_polygon();

    enforce_granularity(s, poly, 4);

    assert(s.num_live_nodes() == 1 && s.num_live_chords() == 0);
    // [C91 §2.2 tex 94 / §2.4 tex 142]: removing the last chord closes
    // ∂C into the single closed arc covering all of C.
    assert(s.num_live_arcs() == 1 &&
           "full contraction leaves the single closed arc");
    for (std::size_t ai = 0; ai < s.num_arcs(); ++ai)
        if (!s.arc(ai).dead) {
            assert(s.arc(ai).first_side == LEFT &&
                   s.arc(ai).last_side == RIGHT &&
                   s.arc(ai).first_edge == 0 && s.arc(ai).last_edge == 0);
            assert(s.arc(ai).edge_count == 4);
        }

    s.normalize(poly);
    s.check_invariants(poly);
    // [C91 §2.3 tex 123]: no exit chord + criterion (i) ⟹ γ-granular.
    assert(s.is_granular(4, poly));
    assert(s.region_weight(0) == 4 &&
           s.region_weight(0) == brute_region_weight(s, 0));

    std::printf("  [PASS] enforce_3region_full\n");
}

// ════════════════════════════════════════════════════════════════
//  2. remove_chord glue chains on the comb (vertex + null-length)
// ════════════════════════════════════════════════════════════════

static void test_comb_glue_chains() {
    CombFixture fx;
    fx.S.check_invariants(fx.C);

    // ── Remove c3 (v3e): LEFT endpoint at v6's east duplicate (vertex
    // junction) glues P4b(8) with the zero-length Z3(9) — a
    // zero-length mate contributes no edges, so P4b keeps ec 1.
    // RIGHT endpoint (15.75, 5) is mid-edge: glues R2(5) + P4a(6) →
    // [8..7], ec 2.
    fx.S.remove_chord(3, fx.C);
    fx.S.assert_tree_property();
    assert(fx.S.node(4).dead && "pocket P4 merged into the outer region");
    assert(!fx.S.arc(8).dead && fx.S.arc(8).edge_count == 1 &&
           "P4b + zero-length Z3: ec unchanged");
    assert(fx.S.arc(9).dead && "Z3 tombstoned by the vertex glue");
    assert(!fx.S.arc(5).dead && fx.S.arc(5).edge_count == 2 &&
           "R2 + P4a glued at the mid-edge endpoint (ec 2)");
    assert(fx.S.arc(6).dead && "P4a tombstoned by the mid-edge glue");

    // ── Remove n7 (null-length, [C91 §2.2 tex 108]): the outer
    // before-arc (glued R2, arc 5), the inner null arc Nz(7), and the
    // outer after-arc P4b(8) fuse into one arc [8..6], ec 3.
    fx.S.remove_chord(7, fx.C);
    fx.S.assert_tree_property();
    assert(fx.S.node(8).dead && "n7's empty inner region merged");
    assert(!fx.S.arc(5).dead && fx.S.arc(5).edge_count == 3 &&
           "[C91 §2.2 tex 108]: null chord removal fuses "
           "before + null + after into one arc");
    assert(fx.S.arc(7).dead && fx.S.arc(8).dead);

    // Weights readable through chord adjacency + C-endpoint pointers.
    assert(fx.S.region_weight(0) == brute_region_weight(fx.S, 0));

    std::printf("  [PASS] comb_glue_chains\n");
}

// ════════════════════════════════════════════════════════════════
//  3. Precondition death tests
// ════════════════════════════════════════════════════════════════

static void test_precondition_asserts() {
    // (a) [C91 §3.3 tex 276]: §3.3 runs on the CONFORMAL output of
    // §3.2 — the raw comb fusion (region 0 has degree 7) must fire.
    assert(assert_fires([]{
        CombFixture fx;
        enforce_granularity(fx.S, fx.C, 100);
    }));

    // (b) [C91 §3.3 tex 276]: S must be γ-semigranular on entry
    // (γ ≥ γ₂); γ = 1 < max region weight 2 must fire.
    assert(assert_fires([]{
        Submap s = build_3region_submap();
        Polygon poly = chain_polygon();
        enforce_granularity(s, poly, 1);
    }));

    // (c) Submap::normalize requires conformality ([C91 §2.4(iv)]).
    assert(assert_fires([]{
        CombFixture fx;
        fx.S.normalize(fx.C);
    }));

    std::printf("  [PASS] precondition_asserts\n");
}

// ════════════════════════════════════════════════════════════════
//  4. Full §3.2 → §3.3 pipeline on the comb (real geometric oracles)
// ════════════════════════════════════════════════════════════════

static void run_pipeline_and_check(std::size_t gamma,
                                   bool expect_full_contraction) {
    CombRig rig;
    auto& fx = rig.fx;

    // [C91 §3.2]: restore conformality (adds 2–3 chords, appends split
    // arc halves — table leaves canonical order).
    restore_conformality(fx.S, fx.C, rig.oracles);
    assert(fx.S.is_conformal());

    // [C91 §3.3 tex 276]: enforce γ-granularity, then normal form.
    enforce_granularity(fx.S, fx.C, gamma);
    fx.S.normalize(fx.C);

    // [C91 Lemma 3.5 tex 279]: normal-form γ-granular conformal submap.
    fx.S.check_invariants(fx.C);
    assert(fx.S.is_conformal());
    assert(fx.S.is_granular(gamma, fx.C));
    assert(!fx.S.tree_decomposition().empty());

    if (expect_full_contraction) {
        // γ ≥ |C| edges: every contraction weight ≤ γ, and a finite
        // tree always has an edge incident upon a degree-<3 node, so
        // criterion (ii) forces contraction down to a single region.
        assert(fx.S.num_live_chords() == 0 &&
               fx.S.num_live_nodes() == 1 &&
               "γ ≥ total weight forces full contraction");
    }

    // region_weight is sound for every region ([C91 §2.2 tex 96/108]
    // junction gluing + the C-endpoint arc pointers).
    for (std::size_t r = 0; r < fx.S.num_nodes(); ++r) {
        if (fx.S.node(r).dead) continue;
        assert(fx.S.region_weight(r) == brute_region_weight(fx.S, r));
        assert(fx.S.region_weight(r) <= gamma &&
               "[C91 §2.3 tex 120]: criterion (i) after §3.3");
    }

    // double_identify works on the re-normalized table ([C91 §2.4
    // tex 144]); every returned arc passes through the query point.
    auto res = fx.S.double_identify(3, SymbolicY{6.0, 2}, fx.C);
    assert(res.count >= 1 && res.count <= 6);
    for (std::size_t ai : res) {
        const Arc& a = fx.S.arc(ai);
        auto [lo, hi] = a.underlying_edge_range(0, 12);
        assert(lo <= 3 && 3 <= hi &&
               "double_identify result must contain the query edge");
    }
}

static void test_pipeline_gammas() {
    run_pipeline_and_check(7, false);
    std::printf("  [PASS] pipeline (gamma=7)\n");
    run_pipeline_and_check(9, false);
    std::printf("  [PASS] pipeline (gamma=9)\n");
    run_pipeline_and_check(100, true);
    std::printf("  [PASS] pipeline (gamma=100)\n");
}

// ════════════════════════════════════════════════════════════════
//  5. Normal form repairs the §3.2-broken table
// ════════════════════════════════════════════════════════════════

static void test_normalize_repairs_table() {
    CombRig rig;
    auto& fx = rig.fx;
    std::size_t arcs_before = fx.S.num_arcs();

    restore_conformality(fx.S, fx.C, rig.oracles);
    assert(fx.S.num_arcs() > arcs_before &&
           "§3.2 appends split arc halves");

    // insert_chord broke canonical order; double_identify must fail
    // fast until normalize re-sorts ([C91 §2.4 tex 144]).
    {
        // Copy the fixture state into the child (CombRig is not
        // trivially shareable across fork, so rebuild inside).
        assert(assert_fires([]{
            CombRig r2;
            restore_conformality(r2.fx.S, r2.fx.C, r2.oracles);
            (void)r2.fx.S.double_identify(3, SymbolicY{6.0, 2}, r2.fx.C);
        }));
    }

    fx.S.normalize(fx.C);
    fx.S.check_invariants(fx.C);
    auto res = fx.S.double_identify(3, SymbolicY{6.0, 2}, fx.C);
    assert(res.count >= 1);

    // normalize is idempotent.
    fx.S.normalize(fx.C);
    fx.S.check_invariants(fx.C);

    std::printf("  [PASS] normalize_repairs_table\n");
}

// ════════════════════════════════════════════════════════════════
//  6. The degree-drop re-check (paper's tex 276 once-only claim)
// ════════════════════════════════════════════════════════════════

// [C91 §3.3 tex 276] claims "chords need be processed only once since
// the removals cannot make any chord removable if it was not already
// so before."  That holds for the weight half of the test but not the
// degree half: contracting a leaf edge lowers the surviving node's
// degree.  This test builds the counterexample shape — the §3.2-added
// hub chords ordered BEFORE the leaf pocket chords — and checks that
// enforcement still reaches the (unique, γ ≥ total weight) fully
// contracted state, which requires re-examining hub chords after their
// pockets are removed.
static void test_degree_drop_recheck() {
    CombRig rig;
    auto& fx = rig.fx;
    restore_conformality(fx.S, fx.C, rig.oracles);
    fx.S.normalize(fx.C);   // compacted canonical state for cloning

    std::size_t nch = fx.S.num_chords();
    assert(nch >= 10 && "comb + §3.2 chords: 8 fused + ≥2 added");

    // Clone with the chord table REVERSED: the §3.2-added chords (the
    // region-0 hub splitters) come first, the leaf pockets last, so
    // the initial in-turn pass reaches every hub chord while both of
    // its endpoints still have degree ≥ 3.
    Submap clone;
    for (std::size_t r = 0; r < fx.S.num_nodes(); ++r) {
        clone.add_node();
        clone.node(r).dead = fx.S.node(r).dead;
    }
    clone.start_vertex = fx.S.start_vertex;
    clone.end_vertex = fx.S.end_vertex;
    for (std::size_t ai = 0; ai < fx.S.num_arcs(); ++ai)
        clone.add_arc(fx.S.arc(ai));
    for (std::size_t ci = nch; ci-- > 0; ) {
        Chord c = fx.S.chord(ci);
        clone.add_chord(c);
    }
    clone.start_arc = fx.S.start_arc;
    clone.end_arc = fx.S.end_arc;
    clone.start_vertex = fx.S.start_vertex;
    clone.end_vertex = fx.S.end_vertex;
    clone.check_invariants(fx.C);

    // γ ≥ total weight: the unique criterion-(ii) fixed point is the
    // single-region submap (any remaining edge would have a leaf and
    // contraction weight ≤ γ).  Without the re-check, a hub chord
    // skipped at degree (3,3)+ would survive and the postcondition
    // assert inside enforce_granularity would fire.
    enforce_granularity(clone, fx.C, 100);
    assert(clone.num_live_chords() == 0 && clone.num_live_nodes() == 1 &&
           "[C91 Lemma 3.5]: γ ≥ total weight must fully contract");

    clone.normalize(fx.C);
    assert(clone.is_granular(100, fx.C));

    std::printf("  [PASS] degree_drop_recheck\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::setbuf(stdout, nullptr);
    std::printf("[C91 §3.3 tests]:\n");
    test_enforce_3region_partial();
    test_enforce_3region_full();
    test_comb_glue_chains();
    test_precondition_asserts();
    test_pipeline_gammas();
    test_normalize_repairs_table();
    test_degree_drop_recheck();
    std::printf("All §3.3 tests passed.\n");
    return 0;
}
