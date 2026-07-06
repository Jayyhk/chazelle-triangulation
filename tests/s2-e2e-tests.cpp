// tests/s2-e2e-tests.cpp — Exhaustive property-based tests for [C91 §2.0–2.4].
//
// Strategy: generate random polygonal curves, build submaps with
// visibility chords, then verify every operation against brute-force
// reference implementations.  Runs thousands of iterations to shake
// out corner cases.

#include "polygon/polygon.h"
#include "polygon/perturbation.h"
#include "submap/submap.h"
#include "submap/tree_decomposition.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <random>
#include <set>
#include <vector>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════════
//  [C91 §0]. Random polygon generator
// ════════════════════════════════════════════════════════════════════

// Generate a random nonclosed polygonal curve with n vertices.
// All y-coordinates are distinct (SoS-safe).  x-coordinates random.
// Vertex indices are 0..n-1.
static Polygon random_polygon(std::mt19937& rng, std::size_t n) {
    assert(n >= 2);
    std::vector<double> ys(n);
    // Distinct y: use consecutive integers + small perturbation.
    for (std::size_t i = 0; i < n; ++i)
        ys[i] = static_cast<double>(i) * 10.0
               + std::uniform_real_distribution<double>(0.1, 9.9)(rng);

    // Shuffle y-values to get non-monotone curves.
    std::shuffle(ys.begin(), ys.end(), rng);

    std::vector<Point> pts(n);
    for (std::size_t i = 0; i < n; ++i) {
        pts[i].x = std::uniform_real_distribution<double>(-100, 100)(rng);
        pts[i].y = ys[i];
        pts[i].index = i;
    }
    return Polygon(std::move(pts));
}

// ════════════════════════════════════════════════════════════════════
//  [C91 §0]b. Brute-force helpers
// ════════════════════════════════════════════════════════════════════

// Brute-force: find all arcs that contain a given (edge_idx, y) point.
// Linear scan of all arcs — O(m).
[[maybe_unused]]
static std::vector<std::size_t> brute_double_identify(
        const Submap& s, std::size_t edge_idx, SymbolicY /*y*/,
        const Polygon& /*poly*/) {
    std::vector<std::size_t> result;
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        const auto& a = s.arc(i);

        // [C91 §2.4 tex 142]: coverage per leg — a wrap-spanning arc is
        // ONE structure whose legs flank C's endpoint turnaround(s).
        // For brute-force, accept any arc covering this edge on either
        // side — the optimized double_identify does y-disambiguation
        // via arc_start_symbolic_y, so its result is a SUBSET of this
        // conservative set.
        if (a.covers(edge_idx, LEFT, s.start_vertex, s.end_vertex) ||
            a.covers(edge_idx, RIGHT, s.start_vertex, s.end_vertex))
            result.push_back(i);
    }
    return result;
}

// Brute-force region weight: linear scan of all arcs for a given region.
static std::size_t brute_region_weight(const Submap& s,
                                        std::size_t node_idx) {
    std::size_t max_ec = 0;
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead) continue;
        if (s.arc(i).region_node == node_idx) {
            max_ec = std::max(max_ec, s.arc(i).edge_count);
        }
    }
    return max_ec;
}

// Brute-force is_conformal: check all node degrees.
static bool brute_is_conformal(const Submap& s) {
    for (std::size_t i = 0; i < s.num_nodes(); ++i) {
        if (s.node(i).dead) continue;
        if (s.node(i).degree() > 4) return false;
    }
    return true;
}

// Brute-force is_semigranular.
static bool brute_is_semigranular(const Submap& s, std::size_t gamma) {
    for (std::size_t i = 0; i < s.num_nodes(); ++i) {
        if (s.node(i).dead) continue;
        if (brute_region_weight(s, i) > gamma) return false;
    }
    return true;
}

// ════════════════════════════════════════════════════════════════════
//  [C91 §1]. Submap builder: construct a valid submap from a polygon
//       by placing chords at interior vertices (y-extrema and
//       non-extrema) to create a multi-region submap.
// ════════════════════════════════════════════════════════════════════

// Information about a chord to be placed.
struct ChordSpec {
    double y;
    std::size_t y_tag;
    std::size_t left_edge;  // edge on LEFT side
    std::size_t right_edge; // edge on RIGHT side
    Side left_side;
    Side right_side;
    bool is_null_length;
};

// Find the edge containing vertex v on a given side.
// For LEFT side: edge v-1 ends at v, edge v starts at v.
// Convention: the edge "before" v in ∂C traversal direction.
//   LEFT (ascending): edge v-1 (the edge ending at v).
//   RIGHT (descending): edge v (the edge starting at v).
[[maybe_unused]]
static std::size_t edge_before_vertex(std::size_t v, Side side,
                                       std::size_t num_edges) {
    if (side == LEFT) {
        // LEFT traversal is ascending.  Edge v-1 ends at vertex v.
        assert(v > 0);
        return v - 1;
    } else {
        // RIGHT traversal is descending.  Edge v starts at vertex v.
        assert(v < num_edges);
        return v;
    }
}

// Build a submap by placing horizontal chords at every interior
// vertex of the polygon.  This creates a visibility-map-like
// structure.  Each interior vertex v gets a chord at y=v.y
// connecting the LEFT side (edge v-1) to the RIGHT side (edge v
// or wherever the horizontal ray hits).
//
// For simplicity, we place chords only at polygon vertices (not at
// arbitrary y-coordinates).  This means both chord endpoints are
// always polygon vertices, so remove_chord won't merge arcs (the
// paper's non-vertex endpoint case).  We test that case separately.
static Submap build_vertex_chord_submap(const Polygon& poly,
                                         std::size_t num_chords_max,
                                         std::mt19937& rng) {
    std::size_t n = poly.num_vertices();
    std::size_t ne = poly.num_edges();

    // Collect interior vertices (not endpoints).
    std::vector<std::size_t> interior_verts;
    for (std::size_t v = 1; v + 1 < n; ++v)
        interior_verts.push_back(v);

    // Shuffle and take up to num_chords_max.
    std::shuffle(interior_verts.begin(), interior_verts.end(), rng);
    std::size_t num_chords = std::min(num_chords_max,
                                       interior_verts.size());
    interior_verts.resize(num_chords);

    // Sort by y to build arcs in ∂C order.
    // Each chord at vertex v splits the LEFT side at edge v-1/v boundary
    // and the RIGHT side at edge v/v-1 boundary.
    // Sort vertices by edge index (= vertex index - 1 for LEFT).
    std::sort(interior_verts.begin(), interior_verts.end());

    // Build the submap:
    // Regions: chord i and chord i+1 define a region between them.
    // With k chords at sorted positions, we get k+1 regions.
    // Region 0: from C start to first chord.
    // Region j: from chord j-1 to chord j.
    // Region k: from last chord to C end.

    Submap s;
    std::size_t num_regions = num_chords + 1;
    for (std::size_t i = 0; i < num_regions; ++i)
        s.add_node();

    // [C91 §2.4(iii) tex 138]: endpoints first, so add_arc classifies
    // the wrap arcs.
    s.start_vertex = 0;
    s.end_vertex = n - 1;

    if (num_chords == 0) {
        // Chordless: one region bounded by the single closed arc
        // ([C91 §2.4 tex 142/138]).
        Arc a{};
        a.first_edge = 0; a.last_edge = 0;
        a.first_side = LEFT; a.last_side = RIGHT;
        a.region_node = 0;
        a.edge_count = 2 * poly.count_nonnull_edges(0, ne - 1);
        std::size_t ai = s.add_arc(a);
        assert(s.start_arc == ai && s.end_arc == ai);
        return s;
    }

    // [C91 §2.4 tex 142]: region 0's arc double-backs around C's start
    // vertex (S) and region k's around C's end vertex (E) — single
    // structures; regions 1..k-1 get one plain arc per ∂C side.
    //
    // Table order ([C91 §2.4(iii) tex 138]): LEFT-starting arcs
    // ascending first_edge (middles, then E), then RIGHT-starting
    // descending (middles, then S — the start-turn arc, last).
    std::vector<std::size_t> left_arc(num_regions, NONE);
    std::vector<std::size_t> right_arc(num_regions, NONE);
    for (std::size_t i = 1; i < num_chords; ++i) {
        std::size_t lo = interior_verts[i - 1];       // starts at v_{i-1}
        std::size_t hi = interior_verts[i] - 1;       // ends at v_i
        Arc a{};
        a.first_edge = lo; a.last_edge = hi;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = i;
        a.edge_count = poly.count_nonnull_edges(lo, hi);
        left_arc[i] = s.add_arc(a);
    }
    {
        // E: region k — from the last chord's LEFT endpoint around C's
        // end vertex to its RIGHT endpoint; ᾱ = [v_{k-1}, ne-1].
        std::size_t v = interior_verts[num_chords - 1];
        Arc a{};
        a.first_edge = v; a.last_edge = v;
        a.first_side = LEFT; a.last_side = RIGHT;
        a.region_node = num_chords;
        a.edge_count = 2 * poly.count_nonnull_edges(v, ne - 1);
        right_arc[num_chords] = left_arc[num_chords] = s.add_arc(a);
    }
    for (std::size_t i = num_chords; i-- > 1; ) {
        // RIGHT arc of region i: [v_i − 1 .. v_{i-1}] descending.
        std::size_t hi = interior_verts[i] - 1;
        std::size_t lo = interior_verts[i - 1];
        Arc a{};
        a.first_edge = hi; a.last_edge = lo;
        a.first_side = RIGHT; a.last_side = RIGHT;
        a.region_node = i;
        a.edge_count = poly.count_nonnull_edges(lo, hi);
        right_arc[i] = s.add_arc(a);
    }
    {
        // S: region 0 — from the first chord's RIGHT endpoint around
        // C's start vertex to its LEFT endpoint; ᾱ = [0, v_0 − 1].
        std::size_t v = interior_verts[0];
        Arc a{};
        a.first_edge = v - 1; a.last_edge = v - 1;
        a.first_side = RIGHT; a.last_side = LEFT;
        a.region_node = 0;
        a.edge_count = 2 * poly.count_nonnull_edges(0, v - 1);
        right_arc[0] = left_arc[0] = s.add_arc(a);
    }

    // Add chords.
    for (std::size_t i = 0; i < num_chords; ++i) {
        std::size_t v = interior_verts[i];
        Chord c{};
        c.region[0] = i;
        c.region[1] = i + 1;
        c.left_edge = v - 1;
        c.right_edge = v; // RIGHT side edge
        c.left_side = LEFT;
        c.right_side = RIGHT;
        c.y = poly.vertex(v).y;
        c.y_tag = v;
        c.is_null_length = false;

        // Adjacent arcs: ONE per polygon-vertex endpoint — the
        // before-arc in ∂C traversal order ([C91 §2.2 tex 94 +
        // §2.4(ii)] 1-slot convention).  Region i's LEFT-side arc ends
        // at the chord's LEFT endpoint; region i+1's RIGHT-side arc at
        // its RIGHT endpoint — for the boundary regions these are the
        // wrap arcs S and E ([C91 §2.4 tex 142]).
        c.left_adj = {{left_arc[i]}, 1};
        c.right_adj = {{right_arc[i + 1]}, 1};

        s.add_chord(c);
    }

    // [C91 §2.4(iii) tex 138]: add_arc auto-registered the endpoint
    // pointers from the wrap classes.
    assert(s.start_arc == right_arc[0] &&
           s.end_arc == left_arc[num_chords]);

    return s;
}

// Build a simpler 2-region submap with a non-vertex chord for testing
// arc merging in remove_chord.
static Submap build_nonvertex_chord_submap(const Polygon& poly,
                                            std::mt19937& rng) {
    std::size_t ne = poly.num_edges();
    std::size_t n = poly.num_vertices();
    if (ne < 2) {
        // Too small for a non-vertex chord.  Chordless submap: the
        // single closed arc ([C91 §2.4 tex 142/138]).
        Submap s;
        s.add_node();
        s.start_vertex = 0; s.end_vertex = 1;
        Arc a{};
        a.first_edge = 0; a.last_edge = 0;
        a.first_side = LEFT; a.last_side = RIGHT;
        a.region_node = 0; a.edge_count = 2 * poly.count_nonnull_edges(0, 0);
        s.add_arc(a);
        return s;
    }

    // Pick a random edge for the chord.
    std::size_t chord_edge = std::uniform_int_distribution<std::size_t>(
        0, ne - 1)(rng);

    // Chord y is midpoint of the edge (not a vertex).
    const auto& e = poly.edge(chord_edge);
    double y_mid = (poly.vertex(e.start_idx).y +
                    poly.vertex(e.end_idx).y) / 2.0;
    // Use a tag that doesn't match any vertex.
    std::size_t y_tag = n + 100;

    Submap s;
    s.add_node(); // r0
    s.add_node(); // r1
    s.start_vertex = 0; s.end_vertex = n - 1;

    // [C91 §2.4 tex 142]: r1's arc double-backs around C's end vertex
    // (E), r0's around C's start vertex (S) — one structure each.
    Arc a{};
    a.first_edge = chord_edge; a.last_edge = chord_edge;
    a.first_side = LEFT; a.last_side = RIGHT;
    a.region_node = 1;
    a.edge_count = 2 * poly.count_nonnull_edges(chord_edge, ne - 1);
    std::size_t E = s.add_arc(a);

    a = {};
    a.first_edge = chord_edge; a.last_edge = chord_edge;
    a.first_side = RIGHT; a.last_side = LEFT;
    a.region_node = 0;
    a.edge_count = 2 * poly.count_nonnull_edges(0, chord_edge);
    std::size_t S = s.add_arc(a);

    // Chord: non-vertex endpoint on both sides (same edge); mid-edge
    // endpoints record {before, after} ([C91 §2.4(ii) tex 137]).
    Chord c{};
    c.region[0] = 0; c.region[1] = 1;
    c.left_edge = chord_edge; c.right_edge = chord_edge;
    c.left_side = LEFT; c.right_side = RIGHT;
    c.y = y_mid; c.y_tag = y_tag;
    c.is_null_length = false;
    c.left_adj = {{S, E}, 2};
    c.right_adj = {{E, S}, 2};
    s.add_chord(c);

    assert(s.start_arc == S && s.end_arc == E);

    return s;
}

// ════════════════════════════════════════════════════════════════════
//  Test 1: SoS perturbation properties
// ════════════════════════════════════════════════════════════════════

static void test_sos_properties(std::mt19937& rng, int iters) {
    std::printf("  test_sos_properties (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(3, 50)(rng);
        auto poly = random_polygon(rng, n);

        // All vertices must have distinct symbolic y.
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                SymbolicY yi = symbolic_y_of(poly.vertex(i));
                SymbolicY yj = symbolic_y_of(poly.vertex(j));
                assert(!symbolic_y_equal(yi, yj) &&
                       "SoS: all vertices must have distinct symbolic y");
                // Total order: exactly one of <, > holds.
                assert((symbolic_y_less(yi, yj) || symbolic_y_less(yj, yi)) &&
                       "SoS: total order violated");
            }
        }

        // Transitivity: sort and verify.
        std::vector<SymbolicY> sorted(n);
        for (std::size_t i = 0; i < n; ++i)
            sorted[i] = symbolic_y_of(poly.vertex(i));
        std::sort(sorted.begin(), sorted.end(), symbolic_y_less);
        for (std::size_t i = 0; i + 1 < n; ++i) {
            assert(symbolic_y_less(sorted[i], sorted[i + 1]) &&
                   "SoS: sorted order not strictly ascending");
        }

        // y-extremum classification: endpoints never extrema.
        assert(!poly.is_y_extremum(0));
        assert(!poly.is_y_extremum(n - 1));

        // count_nonnull_edges: range queries consistent.
        for (std::size_t lo = 0; lo < poly.num_edges(); ++lo) {
            std::size_t total = 0;
            for (std::size_t hi = lo; hi < poly.num_edges(); ++hi) {
                const auto& e = poly.edge(hi);
                bool nonnull = (poly.vertex(e.start_idx).x !=
                                poly.vertex(e.end_idx).x) ||
                               (poly.vertex(e.start_idx).y !=
                                poly.vertex(e.end_idx).y);
                total += nonnull ? 1 : 0;
                assert(poly.count_nonnull_edges(lo, hi) == total &&
                       "count_nonnull_edges mismatch with brute force");
            }
        }
    }
    std::printf("  [PASS] test_sos_properties\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 2: Submap construction invariants
// ════════════════════════════════════════════════════════════════════

static void test_submap_invariants(std::mt19937& rng, int iters) {
    std::printf("  test_submap_invariants (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 30)(rng);
        auto poly = random_polygon(rng, n);

        std::size_t max_chords = std::min(n - 2,
            std::uniform_int_distribution<std::size_t>(0, 10)(rng));
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        // Tree property.
        s.assert_tree_property();

        // Full invariant check.
        s.check_invariants();

        // Conformality (vertex chords create degree ≤ 2 nodes).
        assert(s.is_conformal() &&
               "vertex-chord submap should be conformal (linear chain)");

        // region_weight matches brute force.
        for (std::size_t i = 0; i < s.num_nodes(); ++i) {
            if (s.node(i).dead) continue;
            std::size_t opt = s.region_weight(i);
            std::size_t brute = brute_region_weight(s, i);
            assert(opt == brute &&
                   "region_weight mismatch with brute force");
        }

        // is_conformal matches brute force.
        assert(s.is_conformal() == brute_is_conformal(s));

        // is_semigranular matches brute force for various gamma.
        for (std::size_t g = 0; g <= n + 5; ++g) {
            assert(s.is_semigranular(g) == brute_is_semigranular(s, g) &&
                   "is_semigranular mismatch");
        }
    }
    std::printf("  [PASS] test_submap_invariants\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 3: remove_chord correctness
// ════════════════════════════════════════════════════════════════════

static void test_remove_chord(std::mt19937& rng, int iters) {
    std::printf("  test_remove_chord (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);

        std::size_t max_chords = std::uniform_int_distribution<std::size_t>(
            1, std::min(n - 2, std::size_t(8)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        // Remove chords in random order.
        std::vector<std::size_t> chord_indices;
        for (std::size_t i = 0; i < s.num_chords(); ++i)
            chord_indices.push_back(i);
        std::shuffle(chord_indices.begin(), chord_indices.end(), rng);

        for (std::size_t ci : chord_indices) {
            if (s.chord(ci).dead) continue;

            std::size_t pre_nodes = s.num_live_nodes();
            std::size_t pre_chords = s.num_live_chords();

            s.remove_chord(ci, poly);

            // Tree property preserved.
            s.assert_tree_property();

            // Exactly one node and one chord removed.
            assert(s.num_live_nodes() == pre_nodes - 1);
            assert(s.num_live_chords() == pre_chords - 1);

            // [C91 §2.2 tex 96/108]: junction gluing keeps every live
            // arc reachable (chord before-arc slots + the C-endpoint
            // arc pointers), so region_weight is exact even between
            // cascaded removals.
            for (std::size_t r = 0; r < s.num_nodes(); ++r) {
                if (s.node(r).dead) continue;
                assert(s.region_weight(r) == brute_region_weight(s, r));
            }
        }

        // After all removals: single region, no chords.
        assert(s.num_live_nodes() == 1);
        assert(s.num_live_chords() == 0);
        s.assert_tree_property();

        // Removing every vertex chord — the last removal closes ∂C into
        // the single closed arc ([C91 §2.2 tex 96] / [C91 §2.4 tex 142]).
        assert(s.num_live_arcs() == 1);

        // [C91 §3.3 tex 276]: put S in normal form.
        s.normalize(poly);
        assert(s.num_nodes() == 1);
        assert(s.num_chords() == 0);
        s.check_invariants(poly);
        assert(s.region_weight(0) == brute_region_weight(s, 0));
        assert(!s.tree_decomposition().empty());
    }
    std::printf("  [PASS] test_remove_chord\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 4: remove_chord with arc merging (non-vertex endpoints)
// ════════════════════════════════════════════════════════════════════

static void test_remove_chord_arc_merge(std::mt19937& rng, int iters) {
    std::printf("  test_remove_chord_arc_merge (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);

        auto s = build_nonvertex_chord_submap(poly, rng);

        if (s.num_live_chords() == 0) continue; // skip trivial

        // Removing the only chord glues at both mid-edge endpoints and
        // closes ∂C into the single closed arc ([C91 §2.2 tex 94] /
        // [C91 §2.4 tex 142]).
        s.remove_chord(0, poly);
        s.assert_tree_property();
        assert(s.num_live_nodes() == 1);
        assert(s.num_live_chords() == 0);
        assert(s.num_live_arcs() == 1 &&
               "last-chord removal closes ∂C into one arc");

        // region_weight matches brute force.
        for (std::size_t i = 0; i < s.num_nodes(); ++i) {
            if (s.node(i).dead) continue;
            assert(s.region_weight(i) == brute_region_weight(s, i));
        }

        // Compact and verify.
        s.compact();
        assert(s.num_arcs() == 1);
        s.check_invariants();
    }
    std::printf("  [PASS] test_remove_chord_arc_merge\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 4b: remove_chord — leaf region with a single arc adjacent at
//  BOTH non-vertex endpoints ([C91 §2.2 tex 94] + [C91 §2.3 tex 121]).
//  Gluing at both endpoints chains three arcs into one.
// ════════════════════════════════════════════════════════════════════

static void test_remove_chord_shared_arc(std::mt19937& rng, int iters) {
    std::printf("  test_remove_chord_shared_arc (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t ne = poly.num_edges();

        // Random split edges e1 < e2 for the chord's two endpoints.
        std::size_t e1 = std::uniform_int_distribution<std::size_t>(
            0, ne - 2)(rng);
        std::size_t e2 = std::uniform_int_distribution<std::size_t>(
            e1 + 1, ne - 1)(rng);
        bool left_side_leaf = std::uniform_int_distribution<int>(0, 1)(rng);

        Submap s;
        std::size_t r0 = s.add_node();
        std::size_t r1 = s.add_node();
        s.start_vertex = 0; s.end_vertex = n - 1;

        // [C91 §2.4 tex 142]: r0's single arc runs the LONG way around —
        // through BOTH of C's endpoint turnarounds — one DOUBLE-WRAP
        // structure W ("an arc might wrap around both sides of C");
        // r1's leaf arc A spans between the chord's endpoints.
        Arc a{};
        std::size_t A = NONE, W = NONE;
        if (left_side_leaf) {
            // A = LEFT [e1, e2] (r1); W = LEFT [e2, ne-1] + all of
            // RIGHT + LEFT [0, e1]: first=(e2,LEFT) > last=(e1,LEFT).
            a = {}; a.first_edge = e1; a.last_edge = e2;
            a.first_side = LEFT; a.last_side = LEFT;
            a.region_node = r1; a.edge_count = poly.count_nonnull_edges(e1, e2);
            A = s.add_arc(a);
            a = {}; a.first_edge = e2; a.last_edge = e1;
            a.first_side = LEFT; a.last_side = LEFT;
            a.region_node = r0;
            // [C91 §2.2 tex 106]: nonnull ∂C edges per leg.
            a.edge_count = poly.count_nonnull_edges(e2, ne - 1) +
                           poly.count_nonnull_edges(0, ne - 1) +
                           poly.count_nonnull_edges(0, e1);
            W = s.add_arc(a);
        } else {
            // Mirror on the RIGHT side: A = RIGHT [e2, e1] (r1);
            // W = RIGHT [0, e1] + all of LEFT + RIGHT [e2, ne-1]:
            // first=(e1,RIGHT) < last=(e2,RIGHT).
            a = {}; a.first_edge = e2; a.last_edge = e1;
            a.first_side = RIGHT; a.last_side = RIGHT;
            a.region_node = r1; a.edge_count = poly.count_nonnull_edges(e1, e2);
            A = s.add_arc(a);
            a = {}; a.first_edge = e1; a.last_edge = e2;
            a.first_side = RIGHT; a.last_side = RIGHT;
            a.region_node = r0;
            // [C91 §2.2 tex 106]: nonnull ∂C edges per leg.
            a.edge_count = poly.count_nonnull_edges(0, e1) +
                           poly.count_nonnull_edges(0, ne - 1) +
                           poly.count_nonnull_edges(e2, ne - 1);
            W = s.add_arc(a);
        }
        assert(s.start_arc == W && s.end_arc == W &&
               "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

        // Chord: both endpoints mid-edge (tag matches no vertex); W ends
        // at the first endpoint and starts at the second, A the reverse
        // ([C91 §2.4(ii) tex 137] {before, after} slots).
        Chord c{};
        c.region[0] = r0; c.region[1] = r1;
        if (left_side_leaf) {
            c.left_edge = e1;  c.left_side = LEFT;
            c.right_edge = e2; c.right_side = LEFT;
        } else {
            c.left_edge = e2;  c.left_side = RIGHT;
            c.right_edge = e1; c.right_side = RIGHT;
        }
        const auto& ce = poly.edge(c.left_edge);
        c.y = (poly.vertex(ce.start_idx).y + poly.vertex(ce.end_idx).y) / 2.0;
        c.y_tag = n + 100;
        c.left_adj  = {{W, A}, 2};
        c.right_adj = {{A, W}, 2};
        s.add_chord(c);

        // [C91 §2.3 tex 121]: predicted contraction weight == actual.
        std::size_t simulated = s.simulated_contraction_weight(0, poly);

        std::size_t survivor = s.remove_chord(0, poly);
        s.assert_tree_property();
        assert(s.num_live_nodes() == 1);
        assert(s.num_live_arcs() == 1 &&
               "removing the last chord closes ∂C into the single "
               "closed arc ([C91 §2.4 tex 142])");
        assert(simulated == brute_region_weight(s, survivor) &&
               "simulated_contraction_weight != actual for shared-arc chord");
        assert(s.region_weight(survivor) == brute_region_weight(s, survivor));

        s.compact();
        s.check_invariants(poly);
    }
    std::printf("  [PASS] test_remove_chord_shared_arc\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 5: simulated_contraction_weight vs actual removal
// ════════════════════════════════════════════════════════════════════

static void test_simulated_contraction_weight(std::mt19937& rng, int iters) {
    std::printf("  test_simulated_contraction_weight (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords = std::uniform_int_distribution<std::size_t>(
            1, std::min(n - 2, std::size_t(6)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        // For each chord, compare simulated weight vs actual removal.
        for (std::size_t ci = 0; ci < s.num_chords(); ++ci) {
            if (s.chord(ci).dead) continue;

            std::size_t simulated = s.simulated_contraction_weight(ci, poly);

            // Actually remove in a copy.
            Submap copy = s;
            std::size_t survivor = copy.remove_chord(ci, poly);
            std::size_t actual = brute_region_weight(copy, survivor);

            assert(simulated == actual &&
                   "simulated_contraction_weight != actual post-removal weight");
        }
    }
    std::printf("  [PASS] test_simulated_contraction_weight\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 6: double_identify — optimized vs brute-force
// ════════════════════════════════════════════════════════════════════

static void test_double_identify(std::mt19937& rng, int iters) {
    std::printf("  test_double_identify (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords = std::uniform_int_distribution<std::size_t>(
            0, std::min(n - 2, std::size_t(6)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        // Must compact before double_identify.
        s.compact();

        // Query at every polygon vertex.
        for (std::size_t v = 0; v < n; ++v) {
            SymbolicY qy = symbolic_y_of(poly.vertex(v));

            // Find which edges this vertex is on.
            // Vertex v is the start of edge v (if v < ne) and end of edge v-1.
            if (v < poly.num_edges()) {
                auto result = s.double_identify(v, qy, poly);

                // Every returned arc must actually contain this edge.
                for (std::size_t k = 0; k < result.count; ++k) {
                    std::size_t ai = result.arcs[k];
                    assert(ai < s.num_arcs());
                    assert(!s.arc(ai).dead);
                }

                // Result count must be ≤ 6.
                assert(result.count <= 6);

                // Result must be non-empty (vertex is on this edge).
                assert(result.count > 0 &&
                       "double_identify should find at least one arc "
                       "at a polygon vertex on its edge");
            }
        }

        // Query at random points on random edges.
        for (int q = 0; q < 20; ++q) {
            std::size_t e = std::uniform_int_distribution<std::size_t>(
                0, poly.num_edges() - 1)(rng);
            // Random y between edge endpoints.
            double y0 = poly.vertex(poly.edge(e).start_idx).y;
            double y1 = poly.vertex(poly.edge(e).end_idx).y;
            double qy_val = y0 + std::uniform_real_distribution<double>(
                0.01, 0.99)(rng) * (y1 - y0);
            // Tag that doesn't match any vertex.
            SymbolicY qy{qy_val, n + 200 + static_cast<std::size_t>(q)};

            auto result = s.double_identify(e, qy, poly);
            assert(result.count <= 6);

            // Every returned arc must cover this edge on some ∂C side
            // ([C91 §2.4 tex 142]: per-leg coverage).
            for (std::size_t k = 0; k < result.count; ++k) {
                std::size_t ai = result.arcs[k];
                const auto& a = s.arc(ai);
                assert((a.covers(e, LEFT, s.start_vertex, s.end_vertex) ||
                        a.covers(e, RIGHT, s.start_vertex, s.end_vertex)) &&
                       "double_identify returned arc that doesn't "
                       "contain the queried edge");
            }

            // Should find at least one arc (the edge must be in some arc).
            // (For compacted submaps, every edge is in at least one arc per side.)
            assert(result.count >= 1 &&
                   "double_identify should find at least one arc "
                   "for any edge in the polygon");
        }

        // Query for an edge beyond the polygon — should return 0.
        {
            auto result = s.double_identify(poly.num_edges() + 10,
                                             {0.0, 0}, poly);
            assert(result.count == 0 &&
                   "double_identify should return 0 for out-of-range edge");
        }
    }
    std::printf("  [PASS] test_double_identify\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 7: compact + check_invariants
// ════════════════════════════════════════════════════════════════════

static void test_compact(std::mt19937& rng, int iters) {
    std::printf("  test_compact (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords = std::uniform_int_distribution<std::size_t>(
            1, std::min(n - 2, std::size_t(6)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        // Remove some chords.
        std::size_t to_remove = std::uniform_int_distribution<std::size_t>(
            0, s.num_live_chords())(rng);
        std::vector<std::size_t> removable;
        for (std::size_t ci = 0; ci < s.num_chords(); ++ci)
            if (!s.chord(ci).dead) removable.push_back(ci);
        std::shuffle(removable.begin(), removable.end(), rng);
        for (std::size_t i = 0; i < to_remove && i < removable.size(); ++i)
            s.remove_chord(removable[i], poly);

        std::size_t live_nodes = s.num_live_nodes();
        std::size_t live_chords = s.num_live_chords();
        std::size_t live_arcs = s.num_live_arcs();

        // Collect pre-compact weights.
        std::vector<std::pair<std::size_t, std::size_t>> pre_weights;
        for (std::size_t i = 0; i < s.num_nodes(); ++i) {
            if (!s.node(i).dead)
                pre_weights.push_back({i, brute_region_weight(s, i)});
        }

        // Compact.
        s.compact();

        // Post-compact: raw counts = live counts.
        assert(s.num_nodes() == live_nodes);
        assert(s.num_chords() == live_chords);
        assert(s.num_arcs() == live_arcs);

        // No dead entries.
        for (std::size_t i = 0; i < s.num_nodes(); ++i)
            assert(!s.node(i).dead);
        for (std::size_t i = 0; i < s.num_chords(); ++i)
            assert(!s.chord(i).dead);
        for (std::size_t i = 0; i < s.num_arcs(); ++i)
            assert(!s.arc(i).dead);

        // Full invariant check.
        s.check_invariants();

        // [C91 §3.3]: junction gluing + C-endpoint arc pointers keep
        // region_weight exact after cascaded removals.
        for (std::size_t i = 0; i < s.num_nodes(); ++i)
            assert(s.region_weight(i) == brute_region_weight(s, i));
        for (std::size_t i = 0; i < s.num_arcs(); ++i) {
            assert(s.arc(i).region_node < s.num_nodes());
        }
    }
    std::printf("  [PASS] test_compact\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 8: Tree decomposition properties
// ════════════════════════════════════════════════════════════════════

static std::size_t td_depth(const TreeDecomposition& td, std::size_t idx) {
    if (td.node(idx).is_leaf()) return 0;
    std::size_t ld = td_depth(td, td.node(idx).left_child);
    std::size_t rd = td_depth(td, td.node(idx).right_child);
    return 1 + std::max(ld, rd);
}

static std::size_t td_count_leaves(const TreeDecomposition& td,
                                    std::size_t idx) {
    if (td.node(idx).is_leaf()) return 1;
    return td_count_leaves(td, td.node(idx).left_child)
         + td_count_leaves(td, td.node(idx).right_child);
}

static std::size_t td_count_internal(const TreeDecomposition& td,
                                      std::size_t idx) {
    if (td.node(idx).is_leaf()) return 0;
    return 1 + td_count_internal(td, td.node(idx).left_child)
             + td_count_internal(td, td.node(idx).right_child);
}

static void test_tree_decomposition(std::mt19937& rng, int iters) {
    std::printf("  test_tree_decomposition (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 30)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords = std::uniform_int_distribution<std::size_t>(
            0, std::min(n - 2, std::size_t(10)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);
        s.compact();

        assert(s.is_conformal());
        s.build_tree_decomposition();
        const auto& td = s.tree_decomposition();

        assert(!td.empty());

        std::size_t num_regions = s.num_nodes();
        std::size_t num_chords_live = s.num_chords();

        // Bijection: internal nodes ↔ chords, leaves ↔ regions.
        std::size_t internals = td_count_internal(td, td.root());
        std::size_t leaves = td_count_leaves(td, td.root());
        assert(internals == num_chords_live &&
               "tree decomposition internals must biject with chords");
        assert(leaves == num_regions &&
               "tree decomposition leaves must biject with regions");

        // Total tree decomposition nodes = chords + regions.
        assert(td.size() == num_chords_live + num_regions);

        // Depth bound: O(log m).  Specifically ≤ ceil(log2(m+1)) * 2.
        // (Paper says O(log m); generous bound for small m.)
        if (num_regions > 1) {
            std::size_t depth = td_depth(td, td.root());
            // For n regions, depth ≤ ~2 * log2(n) should hold.
            double log_bound = 2.0 * std::ceil(std::log2(
                static_cast<double>(num_regions) + 1.0)) + 2.0;
            assert(depth <= static_cast<std::size_t>(log_bound) &&
                   "tree decomposition depth exceeds O(log m) bound");
        }

        // Verify parent pointers.
        assert(td.node(td.root()).parent == NONE);
        for (std::size_t i = 0; i < td.size(); ++i) {
            const auto& node = td.node(i);
            if (node.is_internal()) {
                assert(td.node(node.left_child).parent == i);
                assert(td.node(node.right_child).parent == i);
            }
        }

        // Verify leaf region indices are valid and unique.
        std::set<std::size_t> leaf_regions;
        for (std::size_t i = 0; i < td.size(); ++i) {
            if (td.node(i).is_leaf()) {
                assert(td.node(i).region_idx < s.num_nodes());
                assert(leaf_regions.insert(td.node(i).region_idx).second &&
                       "duplicate region in tree decomposition leaves");
            }
        }
        assert(leaf_regions.size() == num_regions);

        // Verify internal chord indices are valid and unique.
        std::set<std::size_t> internal_chords;
        for (std::size_t i = 0; i < td.size(); ++i) {
            if (td.node(i).is_internal()) {
                assert(td.node(i).chord_idx < s.num_chords());
                assert(internal_chords.insert(td.node(i).chord_idx).second &&
                       "duplicate chord in tree decomposition internals");
            }
        }
        assert(internal_chords.size() == num_chords_live);
    }
    std::printf("  [PASS] test_tree_decomposition\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 9: is_granular consistency
// ════════════════════════════════════════════════════════════════════

static void test_granularity(std::mt19937& rng, int iters) {
    std::printf("  test_granularity (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords = std::uniform_int_distribution<std::size_t>(
            0, std::min(n - 2, std::size_t(6)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        for (std::size_t gamma = 0; gamma <= n + 5; ++gamma) {
            bool semi = s.is_semigranular(gamma);
            bool gran = s.is_granular(gamma, poly);

            // Granular implies semigranular.
            if (gran) {
                assert(semi && "granular must imply semigranular");
            }

            // Verify semigranular against brute force.
            assert(semi == brute_is_semigranular(s, gamma));

            // If no chords and semigranular, must be granular.
            if (s.num_live_chords() == 0 && semi) {
                assert(gran && "no-chord semigranular must be granular");
            }

            // Verify granularity condition (ii) directly.
            if (semi && s.num_live_chords() > 0) {
                bool condition_ii = true;
                for (std::size_t ci = 0; ci < s.num_chords(); ++ci) {
                    if (s.chord(ci).dead) continue;
                    const auto& c = s.chord(ci);
                    std::size_t d0 = s.node(c.region[0]).degree();
                    std::size_t d1 = s.node(c.region[1]).degree();
                    if (d0 >= 3 && d1 >= 3) continue;
                    if (s.simulated_contraction_weight(ci, poly) <= gamma) {
                        condition_ii = false;
                        break;
                    }
                }
                assert(gran == condition_ii &&
                       "is_granular mismatch with direct condition (ii) check");
            }
        }
    }
    std::printf("  [PASS] test_granularity\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 10: simulated_contraction_weight with non-vertex chords
// ════════════════════════════════════════════════════════════════════

static void test_simulated_contraction_nonvertex(std::mt19937& rng,
                                                  int iters) {
    std::printf("  test_simulated_contraction_nonvertex (%d iters)...\n",
                iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);
        auto s = build_nonvertex_chord_submap(poly, rng);

        if (s.num_live_chords() == 0) continue;

        std::size_t simulated = s.simulated_contraction_weight(0, poly);

        // Actually remove and measure.
        Submap copy = s;
        std::size_t survivor = copy.remove_chord(0, poly);
        std::size_t actual = brute_region_weight(copy, survivor);

        assert(simulated == actual &&
               "simulated_contraction_weight != actual for non-vertex chord");
    }
    std::printf("  [PASS] test_simulated_contraction_nonvertex\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 11: Full pipeline — build, remove some, compact, tree decomposition, verify
// ════════════════════════════════════════════════════════════════════

static void test_full_pipeline(std::mt19937& rng, int iters) {
    std::printf("  test_full_pipeline (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(5, 30)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords = std::uniform_int_distribution<std::size_t>(
            1, std::min(n - 2, std::size_t(10)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        // Step 1: verify initial state.
        s.assert_tree_property();
        s.check_invariants();
        assert(s.is_conformal());

        // Step 2: remove random subset of chords.
        std::size_t to_remove = std::uniform_int_distribution<std::size_t>(
            0, s.num_live_chords() / 2 + 1)(rng);
        std::vector<std::size_t> removable;
        for (std::size_t ci = 0; ci < s.num_chords(); ++ci)
            if (!s.chord(ci).dead) removable.push_back(ci);
        std::shuffle(removable.begin(), removable.end(), rng);
        for (std::size_t i = 0; i < to_remove && i < removable.size(); ++i) {
            if (s.chord(removable[i]).dead) continue;
            s.remove_chord(removable[i], poly);
            s.assert_tree_property();
        }

        // Step 3: compact.
        s.compact();
        s.check_invariants();
        assert(s.is_conformal());

        // Step 4: normal form ([C91 §3.3 tex 276]) + weight cross-checks.
        s.normalize(poly);
        s.check_invariants(poly);
        for (std::size_t i = 0; i < s.num_nodes(); ++i)
            assert(s.region_weight(i) == brute_region_weight(s, i));
        for (std::size_t i = 0; i < s.num_arcs(); ++i) {
            assert(s.arc(i).region_node < s.num_nodes());
        }

        // Step 5: build tree decomposition.
        s.build_tree_decomposition();
        const auto& td = s.tree_decomposition();
        assert(!td.empty());
        assert(td_count_leaves(td, td.root()) == s.num_nodes());
        assert(td_count_internal(td, td.root()) == s.num_chords());

        // Step 6: double_identify at every vertex.
        for (std::size_t v = 0; v < poly.num_vertices(); ++v) {
            if (v >= poly.num_edges()) continue;
            SymbolicY qy = symbolic_y_of(poly.vertex(v));
            auto result = s.double_identify(v, qy, poly);
            assert(result.count >= 1 && result.count <= 6);
        }

        // Step 7: granularity cross-checks vs brute force ([C91 §2.3
        // tex 120] criterion (i) via is_semigranular).
        {
            std::size_t max_w = 0;
            for (std::size_t i = 0; i < s.num_nodes(); ++i)
                max_w = std::max(max_w, brute_region_weight(s, i));
            assert(s.is_semigranular(max_w));
            if (max_w > 0)
                assert(!s.is_semigranular(max_w - 1));
        }
    }
    std::printf("  [PASS] test_full_pipeline\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 12: Edge cases — minimal polygons
// ════════════════════════════════════════════════════════════════════

static void test_edge_cases() {
    std::printf("  test_edge_cases...\n"); std::fflush(stdout);

    // 2-vertex polygon (1 edge): single region, no chords possible.
    {
        Polygon p({{0, 0, 0}, {1, 1, 1}});
        assert(p.num_edges() == 1);
        assert(p.is_endpoint(0));
        assert(p.is_endpoint(1));

        Submap s;
        s.add_node();
        s.start_vertex = 0; s.end_vertex = 1;
        // The single closed arc ([C91 §2.4 tex 142/138]).
        Arc a{};
        a.first_edge = 0; a.last_edge = 0;
        a.first_side = LEFT; a.last_side = RIGHT;
        a.region_node = 0; a.edge_count = 1;
        std::size_t ai0 = s.add_arc(a);
        assert(s.start_arc == ai0 && s.end_arc == ai0);

        s.check_invariants();
        assert(s.region_weight(0) == 1);
        assert(s.is_conformal());
        assert(s.is_semigranular(1));
        assert(s.is_granular(1, p));

        s.build_tree_decomposition();
        assert(s.tree_decomposition().size() == 1);

        // [C91 §2.4 tex 142/144]: the closed arc covers the point on
        // both ∂C sides but is ONE arc-structure — reported once.
        auto r = s.double_identify(0, {0.0, 0}, p);
        assert(r.count == 1);
    }

    // 3-vertex polygon (2 edges): one interior vertex.
    {
        Polygon p({{0, 0, 0}, {1, 5, 1}, {2, 2, 2}});
        assert(p.num_edges() == 2);
        assert(p.is_y_extremum(1)); // y-max

        // Single region bounded by the closed arc ([C91 §2.4 tex 142]).
        Submap s;
        s.add_node();
        s.start_vertex = 0; s.end_vertex = 2;
        Arc a{};
        a.first_edge = 0; a.last_edge = 0;
        a.first_side = LEFT; a.last_side = RIGHT;
        a.region_node = 0; a.edge_count = 2;
        std::size_t ai0 = s.add_arc(a);
        assert(s.start_arc == ai0 && s.end_arc == ai0);

        s.check_invariants();
        assert(s.region_weight(0) == 2);
    }

    // Zero-length edge polygon.
    {
        Polygon p({{0, 0, 0}, {0, 0, 1}, {1, 1, 2}});
        assert(p.num_edges() == 2);
        assert(p.count_nonnull_edges(0, 0) == 0); // edge 0 is zero-length
        assert(p.count_nonnull_edges(1, 1) == 1);
        assert(p.count_nonnull_edges(0, 1) == 1);
    }

    std::printf("  [PASS] test_edge_cases\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 13: Stress test — many chords, many removals
// ════════════════════════════════════════════════════════════════════

static void test_stress(std::mt19937& rng, int iters) {
    std::printf("  test_stress (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(10, 50)(rng);
        auto poly = random_polygon(rng, n);
        // Use all interior vertices as chords.
        auto s = build_vertex_chord_submap(poly, n - 2, rng);

        s.assert_tree_property();
        s.check_invariants();
        assert(s.is_conformal());

        // Remove all chords in random order.
        std::vector<std::size_t> all_chords;
        for (std::size_t ci = 0; ci < s.num_chords(); ++ci)
            all_chords.push_back(ci);
        std::shuffle(all_chords.begin(), all_chords.end(), rng);

        for (std::size_t ci : all_chords) {
            if (s.chord(ci).dead) continue;

            // Verify simulated weight before removal.
            std::size_t sim_w = s.simulated_contraction_weight(ci, poly);
            Submap copy = s;
            std::size_t survivor = copy.remove_chord(ci, poly);
            std::size_t actual_w = brute_region_weight(copy, survivor);
            assert(sim_w == actual_w);

            s.remove_chord(ci, poly);
            s.assert_tree_property();
        }

        assert(s.num_live_nodes() == 1);
        assert(s.num_live_chords() == 0);

        s.compact();
        s.check_invariants();

        // Build tree decomposition on single-region.
        s.build_tree_decomposition();
        assert(s.tree_decomposition().size() == 1);
    }
    std::printf("  [PASS] test_stress\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 14: Null-length chord at y-extremum
// ════════════════════════════════════════════════════════════════════

static void test_null_length_chord(std::mt19937& rng, int iters) {
    std::printf("  test_null_length_chord (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        // W-shape polygon guaranteed to have y-extrema.
        std::size_t n = std::uniform_int_distribution<std::size_t>(5, 15)(rng);
        auto poly = random_polygon(rng, n);

        // Find y-extrema.
        std::vector<std::size_t> extrema;
        for (std::size_t v = 1; v + 1 < n; ++v) {
            if (poly.is_y_extremum(v))
                extrema.push_back(v);
        }
        if (extrema.empty()) continue;

        // Pick one extremum and build a submap with a null-length chord there.
        std::size_t ext_v = extrema[std::uniform_int_distribution<std::size_t>(
            0, extrema.size() - 1)(rng)];

        Submap s;
        std::size_t r0 = s.add_node(); // main region
        std::size_t r1 = s.add_node(); // empty region inside the null-length chord
        s.start_vertex = 0; s.end_vertex = n - 1;

        // [C91 §2.4 tex 142]: r0's single arc W runs from the null
        // chord's position at the extremum the long way around ∂C —
        // through BOTH of C's endpoint turnarounds — back to it: one
        // DOUBLE-WRAP structure, first=(ext_v,LEFT) > last=(ext_v-1,LEFT).
        Arc a{};
        // N: r1's zero-length arc inside the null-length chord.
        a.first_edge = ext_v; a.last_edge = ext_v;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r1;
        a.edge_count = 0;
        std::size_t N = s.add_arc(a);

        std::size_t ne = poly.num_edges();
        a = {};
        a.first_edge = ext_v; a.last_edge = ext_v - 1;
        a.first_side = LEFT; a.last_side = LEFT;
        a.region_node = r0;
        // [C91 §2.2 tex 106]: nonnull ∂C edges per leg.
        a.edge_count = poly.count_nonnull_edges(ext_v, ne - 1) +
                       poly.count_nonnull_edges(0, ne - 1) +
                       poly.count_nonnull_edges(0, ext_v - 1);
        std::size_t W = s.add_arc(a);
        assert(s.start_arc == W && s.end_arc == W &&
               "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

        // Null-length chord at the extremum — symbolic y is the
        // vertex's own (y, SoS index) ([C91 §2 tex 47]); the outer slot
        // holds W (ends at the first duplicate), the other the inner
        // null arc ([C91 §2.1 tex 72]).
        Chord c{};
        c.region[0] = r0; c.region[1] = r1;
        c.left_edge = ext_v - 1; c.right_edge = ext_v - 1;
        c.left_side = LEFT; c.right_side = LEFT;
        c.y = poly.vertex(ext_v).y; c.y_tag = ext_v;
        c.is_null_length = true;
        c.left_adj = {{W}, 1};
        c.right_adj = {{N}, 1};
        s.add_chord(c);

        // Invariants.
        s.assert_tree_property();
        assert(s.region_weight(r1) == 0); // empty region inside null-length chord
        // W is the null chord's outer before-arc slot ([C91 §2.2 tex
        // 96] alternation) — region_weight is exact.
        assert(brute_region_weight(s, r0) > 0);
        assert(s.region_weight(r0) == brute_region_weight(s, r0));
        assert(s.is_conformal());

        // Remove the null-length chord → single region.
        s.remove_chord(0, poly);
        s.assert_tree_property();
        assert(s.num_live_nodes() == 1);

        s.compact();
        s.check_invariants();
    }
    std::printf("  [PASS] test_null_length_chord\n");
}

// ════════════════════════════════════════════════════════════════════
//  Test 15: Repeated compact is idempotent
// ════════════════════════════════════════════════════════════════════

static void test_compact_idempotent(std::mt19937& rng, int iters) {
    std::printf("  test_compact_idempotent (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 15)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords = std::uniform_int_distribution<std::size_t>(
            1, std::min(n - 2, std::size_t(5)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        // Remove some chords.
        for (std::size_t ci = 0; ci < s.num_chords(); ++ci) {
            if (s.chord(ci).dead) continue;
            if (std::uniform_int_distribution<int>(0, 1)(rng))
                s.remove_chord(ci, poly);
        }

        // First compact.
        s.compact();
        std::size_t n1 = s.num_nodes();
        std::size_t c1 = s.num_chords();
        std::size_t a1 = s.num_arcs();
        s.check_invariants();

        // Second compact: should be no-op.
        s.compact();
        assert(s.num_nodes() == n1);
        assert(s.num_chords() == c1);
        assert(s.num_arcs() == a1);
        s.check_invariants();
    }
    std::printf("  [PASS] test_compact_idempotent\n");
}

// ════════════════════════════════════════════════════════════════════

int main(int argc, char** argv) {
    // Allow overriding seed from command line.
    unsigned seed = 42;
    if (argc > 1) seed = static_cast<unsigned>(std::atoi(argv[1]));
    std::mt19937 rng(seed);

    std::printf("[C91 §2.0–2.4 e2e tests (seed=%u)]:\n", seed);
    std::fflush(stdout);

    std::setbuf(stdout, nullptr); // unbuffered for crash diagnostics
    test_edge_cases();
    test_sos_properties(rng, 500);
    test_submap_invariants(rng, 500);
    test_remove_chord(rng, 500);
    test_remove_chord_arc_merge(rng, 500);
    test_remove_chord_shared_arc(rng, 500);
    test_simulated_contraction_weight(rng, 500);
    test_simulated_contraction_nonvertex(rng, 500);
    test_double_identify(rng, 300);
    test_compact(rng, 500);
    test_tree_decomposition(rng, 300);
    test_granularity(rng, 300);
    test_null_length_chord(rng, 300);
    test_compact_idempotent(rng, 300);
    test_full_pipeline(rng, 300);
    test_stress(rng, 100);

    std::printf("\nAll §2.0–2.4 e2e tests passed (seed=%u).\n", seed);
    return 0;
}
