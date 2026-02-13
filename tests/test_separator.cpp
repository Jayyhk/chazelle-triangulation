/// Tests for PlanarGraph construction + Lipton-Tarjan separator.
///
///   - Graph construction, edge orbits (CW/CCW), face traversal
///   - find_separator balance guarantee
///   - iterated_separator piece-size guarantee

#include "separator/planar_graph.h"
#include "separator/separator.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <numeric>
#include <unordered_set>
#include <vector>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════════
//  PlanarGraph construction
// ════════════════════════════════════════════════════════════════════

static void test_single_edge() {
    std::printf("Test 1: single edge graph\n");

    PlanarGraph g;
    g.add_vertex(1.0);
    g.add_vertex(1.0);
    g.add_edge(0, 1);

    assert(g.num_vertices() == 2);
    assert(g.num_edges() == 1);

    // Both vertices should point to edge 0.
    assert(g.vertex(0).some_edge == 0);
    assert(g.vertex(1).some_edge == 0);

    // Edge endpoints.
    assert(g.edge(0).endpoint[0] == 0);
    assert(g.edge(0).endpoint[1] == 1);

    // CW/CCW of the only edge at each vertex loops to itself.
    assert(g.next_cw(0, 0) == 0);
    assert(g.next_ccw(0, 0) == 0);
    assert(g.next_cw(0, 1) == 0);

    std::printf("  PASS\n\n");
}

static void test_triangle_graph() {
    std::printf("Test 2: triangle graph (3 vertices, 3 edges)\n");

    PlanarGraph g;
    g.add_vertex(1.0);
    g.add_vertex(1.0);
    g.add_vertex(1.0);

    g.add_edge(0, 1);
    g.add_edge(1, 2);
    g.add_edge(2, 0);

    assert(g.num_vertices() == 3);
    assert(g.num_edges() == 3);

    // Each vertex should have exactly 2 edges in its orbit.
    for (std::size_t v = 0; v < 3; ++v) {
        std::vector<std::size_t> incident;
        g.for_each_edge_cw(v, [&](std::size_t ei) {
            incident.push_back(ei);
        });
        assert(incident.size() == 2);
    }

    std::printf("  PASS\n\n");
}

static void test_edge_orbit_order() {
    std::printf("Test 3: CW orbit consistency\n");

    // Build a star: central vertex 0 connected to 1,2,3,4.
    PlanarGraph g;
    g.add_vertex(1.0); // center = 0
    for (int i = 1; i <= 4; ++i)
        g.add_vertex(1.0);

    // Insert edges 0-1, 0-2, 0-3, 0-4.
    for (int i = 1; i <= 4; ++i)
        g.add_edge(0, static_cast<std::size_t>(i));

    // Check that walking CW from vertex 0 visits all 4 edges.
    std::vector<std::size_t> orbit;
    g.for_each_edge_cw(0, [&](std::size_t ei) {
        orbit.push_back(ei);
    });
    assert(orbit.size() == 4);

    // The orbit should be a cyclic permutation.  Verify no duplicates.
    std::unordered_set<std::size_t> unique(orbit.begin(), orbit.end());
    assert(unique.size() == 4);

    std::printf("  PASS\n\n");
}

static void test_face_traversal() {
    std::printf("Test 4: face traversal on triangle\n");

    PlanarGraph g;
    g.add_vertex(1.0);
    g.add_vertex(1.0);
    g.add_vertex(1.0);

    g.add_edge(0, 1);
    g.add_edge(1, 2);
    g.add_edge(2, 0);

    // Trace a face starting from edge 0, heading towards vertex 1.
    auto [e1, d1] = g.face_next(0, 1);
    auto [e2, d2] = g.face_next(e1, d1);
    auto [e3, d3] = g.face_next(e2, d2);

    // After 3 steps we should return to the start (triangle face).
    assert(e3 == 0);
    assert(d3 == 1);

    std::printf("  PASS\n\n");
}

static void test_deleted_edge_skip() {
    std::printf("Test 5: for_each_edge_cw skips deleted edges\n");

    PlanarGraph g;
    g.add_vertex(1.0);
    g.add_vertex(1.0);
    g.add_vertex(1.0);

    g.add_edge(0, 1);
    g.add_edge(0, 2);

    // Delete edge 0.
    g.edge(0).deleted = true;

    std::vector<std::size_t> orbit;
    g.for_each_edge_cw(0, [&](std::size_t ei) {
        orbit.push_back(ei);
    });
    // Only edge 1 should remain.
    assert(orbit.size() == 1);
    assert(orbit[0] == 1);

    std::printf("  PASS\n\n");
}

// ════════════════════════════════════════════════════════════════════
//  PlanarGraph: larger constructions
// ════════════════════════════════════════════════════════════════════

/// Build a cycle graph C_n as a planar graph.
static PlanarGraph build_cycle(std::size_t n) {
    PlanarGraph g;
    for (std::size_t i = 0; i < n; ++i)
        g.add_vertex(1.0);
    for (std::size_t i = 0; i < n; ++i)
        g.add_edge(i, (i + 1) % n);
    return g;
}

/// Build a grid graph of rows × cols vertices.
static PlanarGraph build_grid(std::size_t rows, std::size_t cols) {
    PlanarGraph g;
    for (std::size_t i = 0; i < rows * cols; ++i)
        g.add_vertex(1.0);
    auto idx = [cols](std::size_t r, std::size_t c) { return r * cols + c; };
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t c = 0; c < cols; ++c) {
            if (c + 1 < cols) g.add_edge(idx(r, c), idx(r, c + 1));
            if (r + 1 < rows) g.add_edge(idx(r, c), idx(r + 1, c));
        }
    }
    return g;
}

static void test_cycle_degree() {
    std::printf("Test 6: cycle graph degree = 2\n");

    auto g = build_cycle(10);
    for (std::size_t v = 0; v < 10; ++v) {
        std::size_t deg = 0;
        g.for_each_edge_cw(v, [&](std::size_t) { ++deg; });
        assert(deg == 2);
    }

    std::printf("  PASS\n\n");
}

static void test_grid_graph() {
    std::printf("Test 7: 4x4 grid graph construction\n");

    auto g = build_grid(4, 4);
    assert(g.num_vertices() == 16);
    // Horizontal edges: 4 rows × 3 = 12.
    // Vertical edges: 3 × 4 cols = 12.
    assert(g.num_edges() == 24);

    // Corner vertex (0,0) should have degree 2.
    std::size_t deg_corner = 0;
    g.for_each_edge_cw(0, [&](std::size_t) { ++deg_corner; });
    assert(deg_corner == 2);

    // Interior vertex (1,1) = index 5 should have degree 4.
    std::size_t deg_interior = 0;
    g.for_each_edge_cw(5, [&](std::size_t) { ++deg_interior; });
    assert(deg_interior == 4);

    std::printf("  PASS\n\n");
}

// ════════════════════════════════════════════════════════════════════
//  Lipton-Tarjan separator
// ════════════════════════════════════════════════════════════════════

static void test_separator_small_triangle() {
    std::printf("Test 8: find_separator on triangle graph\n");

    PlanarGraph g;
    for (int i = 0; i < 3; ++i) g.add_vertex(1.0);
    g.add_edge(0, 1);
    g.add_edge(1, 2);
    g.add_edge(2, 0);

    auto result = find_separator(g);

    // Verify partition: A ∪ B ∪ C = {0,1,2}.
    std::size_t total = result.A.size() + result.B.size() + result.C.size();
    assert(total == 3);

    std::printf("  A=%zu B=%zu C=%zu\n",
                result.A.size(), result.B.size(), result.C.size());
    std::printf("  PASS\n\n");
}

static void test_separator_cycle() {
    std::printf("Test 9: find_separator on cycle C_12\n");

    auto g = build_cycle(12);
    auto result = find_separator(g);

    std::size_t total = result.A.size() + result.B.size() + result.C.size();
    assert(total == 12);

    // Disjointness: no vertex in two sets.
    std::unordered_set<std::size_t> all;
    for (auto v : result.A) { assert(all.insert(v).second); }
    for (auto v : result.B) { assert(all.insert(v).second); }
    for (auto v : result.C) { assert(all.insert(v).second); }

    // Balance: neither A nor B has > 2/3 of total cost.
    // With uniform cost 1.0, this means each should have ≤ 8 vertices.
    assert(result.A.size() <= 8);
    assert(result.B.size() <= 8);

    // Separator bound: |C| ≤ 2√(2n) ≈ 2√24 ≈ 9.8
    assert(result.C.size() <= 10);

    std::printf("  A=%zu B=%zu C=%zu\n",
                result.A.size(), result.B.size(), result.C.size());
    std::printf("  PASS\n\n");
}

static void test_separator_grid() {
    std::printf("Test 10: find_separator on 5x5 grid\n");

    auto g = build_grid(5, 5);
    auto result = find_separator(g);

    std::size_t n = 25;
    std::size_t total = result.A.size() + result.B.size() + result.C.size();
    assert(total == n);

    // Balance check.
    double cost_A = static_cast<double>(result.A.size());
    double cost_B = static_cast<double>(result.B.size());
    assert(cost_A <= 2.0 / 3.0 * static_cast<double>(n) + 1);
    assert(cost_B <= 2.0 / 3.0 * static_cast<double>(n) + 1);

    // Separator size bound: |C| ≤ 2√(2·25) ≈ 14.14
    assert(result.C.size() <= 15);

    std::printf("  A=%zu B=%zu C=%zu\n",
                result.A.size(), result.B.size(), result.C.size());
    std::printf("  PASS\n\n");
}

static void test_separator_disjointness() {
    std::printf("Test 11: separator partition disjointness (larger graph)\n");

    auto g = build_grid(6, 6);
    auto result = find_separator(g);

    std::unordered_set<std::size_t> setA(result.A.begin(), result.A.end());
    std::unordered_set<std::size_t> setB(result.B.begin(), result.B.end());
    std::unordered_set<std::size_t> setC(result.C.begin(), result.C.end());

    // No overlaps.
    for (auto v : result.A) assert(setB.count(v) == 0 && setC.count(v) == 0);
    for (auto v : result.B) assert(setA.count(v) == 0 && setC.count(v) == 0);
    for (auto v : result.C) assert(setA.count(v) == 0 && setB.count(v) == 0);

    // Complete: all vertices accounted for.
    std::size_t total = setA.size() + setB.size() + setC.size();
    assert(total == 36);

    std::printf("  PASS\n\n");
}

static void test_separator_no_cross_edges() {
    std::printf("Test 12: no edges between A and B (separator property)\n");

    auto g = build_grid(5, 5);
    auto result = find_separator(g);

    std::unordered_set<std::size_t> setA(result.A.begin(), result.A.end());
    std::unordered_set<std::size_t> setB(result.B.begin(), result.B.end());

    // Check that no edge connects a vertex in A to a vertex in B directly.
    for (std::size_t ei = 0; ei < g.num_edges(); ++ei) {
        auto& e = g.edge(ei);
        if (e.deleted) continue;
        bool u_in_A = setA.count(e.endpoint[0]) > 0;
        bool v_in_B = setB.count(e.endpoint[1]) > 0;
        bool u_in_B = setB.count(e.endpoint[0]) > 0;
        bool v_in_A = setA.count(e.endpoint[1]) > 0;
        assert(!(u_in_A && v_in_B));
        assert(!(u_in_B && v_in_A));
    }

    std::printf("  PASS\n\n");
}

// ════════════════════════════════════════════════════════════════════
//  Iterated separator
// ════════════════════════════════════════════════════════════════════

static void test_iterated_separator_small() {
    std::printf("Test 13: iterated_separator on C_10 (max piece size 4)\n");

    auto g = build_cycle(10);
    auto result = iterated_separator(g, 4);

    // Every piece should have ≤ 4 vertices.
    for (auto& piece : result.pieces) {
        assert(piece.size() <= 4);
    }

    // All vertices accounted for.
    std::size_t accounted = 0;
    for (auto& piece : result.pieces) accounted += piece.size();
    // Plus separator nodes.
    for (std::size_t v = 0; v < 10; ++v) {
        if (result.separator_nodes[v]) ++accounted;
    }
    assert(accounted == 10);

    std::printf("  %zu pieces, ", result.pieces.size());
    std::size_t sep_count = 0;
    for (auto b : result.separator_nodes) if (b) ++sep_count;
    std::printf("%zu separator nodes\n", sep_count);
    std::printf("  PASS\n\n");
}

static void test_iterated_separator_grid() {
    std::printf("Test 14: iterated_separator on 6x6 grid (max piece size 10)\n");

    auto g = build_grid(6, 6);
    auto result = iterated_separator(g, 10);

    // Every piece should have ≤ 10 vertices.
    for (auto& piece : result.pieces) {
        assert(piece.size() <= 10);
    }

    // All 36 vertices accounted for.
    std::size_t accounted = 0;
    for (auto& piece : result.pieces) accounted += piece.size();
    for (std::size_t v = 0; v < 36; ++v) {
        if (result.separator_nodes[v]) ++accounted;
    }
    assert(accounted == 36);

    // vertex_piece should be consistent.
    for (std::size_t pi = 0; pi < result.pieces.size(); ++pi) {
        for (std::size_t v : result.pieces[pi]) {
            assert(result.vertex_piece[v] == pi);
        }
    }
    // Separator nodes should have vertex_piece == NONE.
    for (std::size_t v = 0; v < 36; ++v) {
        if (result.separator_nodes[v]) {
            assert(result.vertex_piece[v] == NONE);
        }
    }

    std::printf("  %zu pieces\n", result.pieces.size());
    std::printf("  PASS\n\n");
}

// ════════════════════════════════════════════════════════════════════
//  main
// ════════════════════════════════════════════════════════════════════

int main() {
    std::printf("=== PlanarGraph + Separator tests ===\n\n");

    test_single_edge();
    test_triangle_graph();
    test_edge_orbit_order();
    test_face_traversal();
    test_deleted_edge_skip();
    test_cycle_degree();
    test_grid_graph();
    test_separator_small_triangle();
    test_separator_cycle();
    test_separator_grid();
    test_separator_disjointness();
    test_separator_no_cross_edges();
    test_iterated_separator_small();
    test_iterated_separator_grid();

    std::printf("All %d tests passed.\n", 14);
    return 0;
}
