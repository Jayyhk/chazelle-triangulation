#include "algorithm/polygon/perturbation.h"
#include "algorithm/polygon/polygon.h"
#include "algorithm/submap/submap.h"
#include "algorithm/submap/tree_decomposition.h"

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

static Polygon random_polygon(std::mt19937& rng, std::size_t n) {
    assert(n >= 2);
    std::vector<Exact> ys(n);

    for (std::size_t i = 0; i < n; ++i)
        ys[i] =
            static_cast<Exact>(i) * 10.0 + std::uniform_real_distribution<double>(0.1, 9.9)(rng);

    std::shuffle(ys.begin(), ys.end(), rng);

    std::vector<Point> pts(n);
    for (std::size_t i = 0; i < n; ++i) {
        pts[i].x = std::uniform_real_distribution<double>(-100, 100)(rng);
        pts[i].y = ys[i];
        pts[i].index = i;
    }
    return Polygon(std::move(pts));
}

[[maybe_unused]]
static std::vector<std::size_t> brute_double_identify(const Submap& s, std::size_t edge_idx,
                                                      const SymbolicY&, const Polygon&) {
    std::vector<std::size_t> result;
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead)
            continue;
        const auto& a = s.arc(i);

        if (a.covers(edge_idx, LEFT, s.start_vertex, s.end_vertex) ||
            a.covers(edge_idx, RIGHT, s.start_vertex, s.end_vertex))
            result.push_back(i);
    }
    return result;
}

static std::size_t brute_region_weight(const Submap& s, std::size_t node_idx) {
    std::size_t max_ec = 0;
    for (std::size_t i = 0; i < s.num_arcs(); ++i) {
        if (s.arc(i).dead)
            continue;
        if (s.arc(i).region_node == node_idx) {
            max_ec = std::max(max_ec, s.arc(i).edge_count);
        }
    }
    return max_ec;
}

static bool brute_is_conformal(const Submap& s) {
    for (std::size_t i = 0; i < s.num_nodes(); ++i) {
        if (s.node(i).dead)
            continue;
        if (s.node(i).degree() > 4)
            return false;
    }
    return true;
}

static bool brute_is_semigranular(const Submap& s, std::size_t granularity) {
    for (std::size_t i = 0; i < s.num_nodes(); ++i) {
        if (s.node(i).dead)
            continue;
        if (brute_region_weight(s, i) > granularity)
            return false;
    }
    return true;
}

struct ChordSpec {
    Exact y;
    std::size_t y_tag;
    std::size_t left_edge;
    std::size_t right_edge;
    Side left_side;
    Side right_side;
    bool is_null_length;
};

[[maybe_unused]]
static std::size_t edge_before_vertex(std::size_t v, Side side, std::size_t num_edges) {
    if (side == LEFT) {
        assert(v > 0);
        return v - 1;
    } else {
        assert(v < num_edges);
        return v;
    }
}

static Submap build_vertex_chord_submap(const Polygon& poly, std::size_t num_chords_max,
                                        std::mt19937& rng) {
    std::size_t n = poly.num_vertices();
    std::size_t ne = poly.num_edges();

    std::vector<std::size_t> interior_verts;
    for (std::size_t v = 1; v + 1 < n; ++v)
        interior_verts.push_back(v);

    std::shuffle(interior_verts.begin(), interior_verts.end(), rng);
    std::size_t num_chords = std::min(num_chords_max, interior_verts.size());
    interior_verts.resize(num_chords);

    std::sort(interior_verts.begin(), interior_verts.end());

    Submap s;
    std::size_t num_regions = num_chords + 1;
    for (std::size_t i = 0; i < num_regions; ++i)
        s.add_node();

    s.start_vertex = 0;
    s.end_vertex = n - 1;

    if (num_chords == 0) {
        Arc a{};
        a.first_edge = 0;
        a.last_edge = 0;
        a.first_side = LEFT;
        a.last_side = RIGHT;
        a.region_node = 0;
        a.edge_count = 2 * poly.count_nonnull_edges(0, ne - 1);
        std::size_t ai = s.add_arc(a);
        assert(s.start_arc == ai && s.end_arc == ai);
        return s;
    }

    std::vector<std::size_t> left_arc(num_regions, NONE);
    std::vector<std::size_t> right_arc(num_regions, NONE);
    for (std::size_t i = 1; i < num_chords; ++i) {
        std::size_t lo = interior_verts[i - 1];
        std::size_t hi = interior_verts[i] - 1;
        Arc a{};
        a.first_edge = lo;
        a.last_edge = hi;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = i;
        a.edge_count = poly.count_nonnull_edges(lo, hi);
        left_arc[i] = s.add_arc(a);
    }
    {
        std::size_t v = interior_verts[num_chords - 1];
        Arc a{};
        a.first_edge = v;
        a.last_edge = v;
        a.first_side = LEFT;
        a.last_side = RIGHT;
        a.region_node = num_chords;
        a.edge_count = 2 * poly.count_nonnull_edges(v, ne - 1);
        right_arc[num_chords] = left_arc[num_chords] = s.add_arc(a);
    }
    for (std::size_t i = num_chords; i-- > 1;) {
        std::size_t hi = interior_verts[i] - 1;
        std::size_t lo = interior_verts[i - 1];
        Arc a{};
        a.first_edge = hi;
        a.last_edge = lo;
        a.first_side = RIGHT;
        a.last_side = RIGHT;
        a.region_node = i;
        a.edge_count = poly.count_nonnull_edges(lo, hi);
        right_arc[i] = s.add_arc(a);
    }
    {
        std::size_t v = interior_verts[0];
        Arc a{};
        a.first_edge = v - 1;
        a.last_edge = v - 1;
        a.first_side = RIGHT;
        a.last_side = LEFT;
        a.region_node = 0;
        a.edge_count = 2 * poly.count_nonnull_edges(0, v - 1);
        right_arc[0] = left_arc[0] = s.add_arc(a);
    }

    for (std::size_t i = 0; i < num_chords; ++i) {
        std::size_t v = interior_verts[i];
        Chord c{};
        c.region[0] = i;
        c.region[1] = i + 1;
        c.left_edge = v - 1;
        c.right_edge = v;
        c.left_side = LEFT;
        c.right_side = RIGHT;
        c.y = poly.vertex(v).y;
        c.y_tag = v;
        c.is_null_length = false;

        c.left_adj = {{left_arc[i]}, 1};
        c.right_adj = {{right_arc[i + 1]}, 1};

        s.add_chord(c);
    }

    assert(s.start_arc == right_arc[0] && s.end_arc == left_arc[num_chords]);

    return s;
}

static Submap build_nonvertex_chord_submap(const Polygon& poly, std::mt19937& rng) {
    std::size_t ne = poly.num_edges();
    std::size_t n = poly.num_vertices();
    if (ne < 2) {
        Submap s;
        s.add_node();
        s.start_vertex = 0;
        s.end_vertex = 1;
        Arc a{};
        a.first_edge = 0;
        a.last_edge = 0;
        a.first_side = LEFT;
        a.last_side = RIGHT;
        a.region_node = 0;
        a.edge_count = 2 * poly.count_nonnull_edges(0, 0);
        s.add_arc(a);
        return s;
    }

    std::size_t chord_edge = std::uniform_int_distribution<std::size_t>(0, ne - 1)(rng);

    const auto& e = poly.edge(chord_edge);
    Exact y_mid = (poly.vertex(e.start_idx).y + poly.vertex(e.end_idx).y) / 2.0;

    std::size_t y_tag = n + 100;

    Submap s;
    s.add_node();
    s.add_node();
    s.start_vertex = 0;
    s.end_vertex = n - 1;

    Arc a{};
    a.first_edge = chord_edge;
    a.last_edge = chord_edge;
    a.first_side = LEFT;
    a.last_side = RIGHT;
    a.region_node = 1;
    a.edge_count = 2 * poly.count_nonnull_edges(chord_edge, ne - 1);
    std::size_t E = s.add_arc(a);

    a = {};
    a.first_edge = chord_edge;
    a.last_edge = chord_edge;
    a.first_side = RIGHT;
    a.last_side = LEFT;
    a.region_node = 0;
    a.edge_count = 2 * poly.count_nonnull_edges(0, chord_edge);
    std::size_t submap = s.add_arc(a);

    Chord c{};
    c.region[0] = 0;
    c.region[1] = 1;
    c.left_edge = chord_edge;
    c.right_edge = chord_edge;
    c.left_side = LEFT;
    c.right_side = RIGHT;
    c.y = y_mid;
    c.y_tag = y_tag;
    c.is_null_length = false;
    c.left_adj = {{submap, E}, 2};
    c.right_adj = {{E, submap}, 2};
    s.add_chord(c);

    assert(s.start_arc == submap && s.end_arc == E);

    return s;
}

static void test_sos_properties(std::mt19937& rng, int iters) {
    std::printf("  test_sos_properties (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(3, 50)(rng);
        auto poly = random_polygon(rng, n);

        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                SymbolicY yi = symbolic_y_of(poly.vertex(i));
                SymbolicY yj = symbolic_y_of(poly.vertex(j));
                assert(!symbolic_y_equal(yi, yj) &&
                       "SoS: all vertices must have distinct symbolic y");

                assert((symbolic_y_less(yi, yj) || symbolic_y_less(yj, yi)) &&
                       "SoS: total order violated");
            }
        }

        std::vector<SymbolicY> sorted(n);
        for (std::size_t i = 0; i < n; ++i)
            sorted[i] = symbolic_y_of(poly.vertex(i));
        std::sort(sorted.begin(), sorted.end(), symbolic_y_less);
        for (std::size_t i = 0; i + 1 < n; ++i) {
            assert(symbolic_y_less(sorted[i], sorted[i + 1]) &&
                   "SoS: sorted order not strictly ascending");
        }

        assert(!poly.is_y_extremum(0));
        assert(!poly.is_y_extremum(n - 1));

        for (std::size_t lo = 0; lo < poly.num_edges(); ++lo) {
            std::size_t total = 0;
            for (std::size_t hi = lo; hi < poly.num_edges(); ++hi) {
                const auto& e = poly.edge(hi);
                bool nonnull = (poly.vertex(e.start_idx).x != poly.vertex(e.end_idx).x) ||
                               (poly.vertex(e.start_idx).y != poly.vertex(e.end_idx).y);
                total += nonnull ? 1 : 0;
                assert(poly.count_nonnull_edges(lo, hi) == total &&
                       "count_nonnull_edges mismatch with brute force");
            }
        }
    }
    std::printf("  [PASS] test_sos_properties\n");
}

static void test_submap_invariants(std::mt19937& rng, int iters) {
    std::printf("  test_submap_invariants (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 30)(rng);
        auto poly = random_polygon(rng, n);

        std::size_t max_chords =
            std::min(n - 2, std::uniform_int_distribution<std::size_t>(0, 10)(rng));
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        s.assert_tree_property();

        s.check_invariants();

        assert(s.is_conformal() && "vertex-chord submap should be conformal (linear chain)");

        for (std::size_t i = 0; i < s.num_nodes(); ++i) {
            if (s.node(i).dead)
                continue;
            std::size_t opt = s.region_weight(i);
            std::size_t brute = brute_region_weight(s, i);
            assert(opt == brute && "region_weight mismatch with brute force");
        }

        assert(s.is_conformal() == brute_is_conformal(s));

        for (std::size_t g = 0; g <= n + 5; ++g) {
            assert(s.is_semigranular(g) == brute_is_semigranular(s, g) &&
                   "is_semigranular mismatch");
        }
    }
    std::printf("  [PASS] test_submap_invariants\n");
}

static void test_remove_chord(std::mt19937& rng, int iters) {
    std::printf("  test_remove_chord (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);

        std::size_t max_chords =
            std::uniform_int_distribution<std::size_t>(1, std::min(n - 2, std::size_t(8)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        std::vector<std::size_t> chord_indices;
        chord_indices.reserve(s.num_chords());
        for (std::size_t i = 0; i < s.num_chords(); ++i)
            chord_indices.push_back(i);
        std::shuffle(chord_indices.begin(), chord_indices.end(), rng);

        for (std::size_t ci : chord_indices) {
            if (s.chord(ci).dead)
                continue;

            std::size_t pre_nodes = s.num_live_nodes();
            std::size_t pre_chords = s.num_live_chords();

            s.remove_chord(ci, poly);

            s.assert_tree_property();

            assert(s.num_live_nodes() == pre_nodes - 1);
            assert(s.num_live_chords() == pre_chords - 1);

            for (std::size_t r = 0; r < s.num_nodes(); ++r) {
                if (s.node(r).dead)
                    continue;
                assert(s.region_weight(r) == brute_region_weight(s, r));
            }
        }

        assert(s.num_live_nodes() == 1);
        assert(s.num_live_chords() == 0);
        s.assert_tree_property();

        assert(s.num_live_arcs() == 1);

        s.normalize(poly);
        assert(s.num_nodes() == 1);
        assert(s.num_chords() == 0);
        s.check_invariants(poly);
        assert(s.region_weight(0) == brute_region_weight(s, 0));
        assert(!s.tree_decomposition().empty());
    }
    std::printf("  [PASS] test_remove_chord\n");
}

static void test_remove_chord_arc_merge(std::mt19937& rng, int iters) {
    std::printf("  test_remove_chord_arc_merge (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);

        auto s = build_nonvertex_chord_submap(poly, rng);

        if (s.num_live_chords() == 0)
            continue;

        s.remove_chord(0, poly);
        s.assert_tree_property();
        assert(s.num_live_nodes() == 1);
        assert(s.num_live_chords() == 0);
        assert(s.num_live_arcs() == 1 && "last-chord removal closes ∂C into one arc");

        for (std::size_t i = 0; i < s.num_nodes(); ++i) {
            if (s.node(i).dead)
                continue;
            assert(s.region_weight(i) == brute_region_weight(s, i));
        }

        s.compact();
        assert(s.num_arcs() == 1);
        s.check_invariants();
    }
    std::printf("  [PASS] test_remove_chord_arc_merge\n");
}

static void test_remove_chord_shared_arc(std::mt19937& rng, int iters) {
    std::printf("  test_remove_chord_shared_arc (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t ne = poly.num_edges();

        std::size_t e1 = std::uniform_int_distribution<std::size_t>(0, ne - 2)(rng);
        std::size_t e2 = std::uniform_int_distribution<std::size_t>(e1 + 1, ne - 1)(rng);
        bool left_side_leaf = std::uniform_int_distribution<int>(0, 1)(rng);

        Submap s;
        std::size_t r0 = s.add_node();
        std::size_t r1 = s.add_node();
        s.start_vertex = 0;
        s.end_vertex = n - 1;

        Arc a{};
        std::size_t A = NONE, W = NONE;
        if (left_side_leaf) {
            a = {};
            a.first_edge = e1;
            a.last_edge = e2;
            a.first_side = LEFT;
            a.last_side = LEFT;
            a.region_node = r1;
            a.edge_count = poly.count_nonnull_edges(e1, e2);
            A = s.add_arc(a);
            a = {};
            a.first_edge = e2;
            a.last_edge = e1;
            a.first_side = LEFT;
            a.last_side = LEFT;
            a.region_node = r0;

            a.edge_count = poly.count_nonnull_edges(e2, ne - 1) +
                           poly.count_nonnull_edges(0, ne - 1) + poly.count_nonnull_edges(0, e1);
            W = s.add_arc(a);
        } else {
            a = {};
            a.first_edge = e2;
            a.last_edge = e1;
            a.first_side = RIGHT;
            a.last_side = RIGHT;
            a.region_node = r1;
            a.edge_count = poly.count_nonnull_edges(e1, e2);
            A = s.add_arc(a);
            a = {};
            a.first_edge = e1;
            a.last_edge = e2;
            a.first_side = RIGHT;
            a.last_side = RIGHT;
            a.region_node = r0;

            a.edge_count = poly.count_nonnull_edges(0, e1) + poly.count_nonnull_edges(0, ne - 1) +
                           poly.count_nonnull_edges(e2, ne - 1);
            W = s.add_arc(a);
        }
        assert(s.start_arc == W && s.end_arc == W &&
               "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

        Chord c{};
        c.region[0] = r0;
        c.region[1] = r1;
        if (left_side_leaf) {
            c.left_edge = e1;
            c.left_side = LEFT;
            c.right_edge = e2;
            c.right_side = LEFT;
        } else {
            c.left_edge = e2;
            c.left_side = RIGHT;
            c.right_edge = e1;
            c.right_side = RIGHT;
        }
        const auto& ce = poly.edge(c.left_edge);
        c.y = (poly.vertex(ce.start_idx).y + poly.vertex(ce.end_idx).y) / 2.0;
        c.y_tag = n + 100;
        c.left_adj = {{W, A}, 2};
        c.right_adj = {{A, W}, 2};
        s.add_chord(c);

        std::size_t simulated = s.simulated_contraction_weight(0, poly);

        std::size_t survivor = s.remove_chord(0, poly);
        s.assert_tree_property();
        assert(s.num_live_nodes() == 1);
        assert(s.num_live_arcs() == 1 && "removing the last chord closes ∂C into the single "
                                         "closed arc ([C91 §2.4 tex 142])");
        assert(simulated == brute_region_weight(s, survivor) &&
               "simulated_contraction_weight != actual for shared-arc chord");
        assert(s.region_weight(survivor) == brute_region_weight(s, survivor));

        s.compact();
        s.check_invariants(poly);
    }
    std::printf("  [PASS] test_remove_chord_shared_arc\n");
}

static void test_simulated_contraction_weight(std::mt19937& rng, int iters) {
    std::printf("  test_simulated_contraction_weight (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords =
            std::uniform_int_distribution<std::size_t>(1, std::min(n - 2, std::size_t(6)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        for (std::size_t ci = 0; ci < s.num_chords(); ++ci) {
            if (s.chord(ci).dead)
                continue;

            std::size_t simulated = s.simulated_contraction_weight(ci, poly);

            Submap copy = s;
            std::size_t survivor = copy.remove_chord(ci, poly);
            std::size_t actual = brute_region_weight(copy, survivor);

            assert(simulated == actual &&
                   "simulated_contraction_weight != actual post-removal weight");
        }
    }
    std::printf("  [PASS] test_simulated_contraction_weight\n");
}

static void test_double_identify(std::mt19937& rng, int iters) {
    std::printf("  test_double_identify (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords =
            std::uniform_int_distribution<std::size_t>(0, std::min(n - 2, std::size_t(6)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        s.compact();

        for (std::size_t v = 0; v < n; ++v) {
            SymbolicY qy = symbolic_y_of(poly.vertex(v));

            if (v < poly.num_edges()) {
                auto result = s.double_identify(v, qy, poly);

                for (std::size_t k = 0; k < result.count; ++k) {
                    std::size_t ai = result.arcs[k];
                    assert(ai < s.num_arcs());
                    assert(!s.arc(ai).dead);
                    const auto& a = s.arc(ai);
                    assert((a.covers(v, LEFT, s.start_vertex, s.end_vertex) ||
                            a.covers(v, RIGHT, s.start_vertex, s.end_vertex)) &&
                           "double_identify returned arc that doesn't "
                           "contain the queried edge");
                }

                assert(result.count <= 6);

                assert(result.count > 0 && "double_identify should find at least one arc "
                                           "at a polygon vertex on its edge");
            }
        }

        for (int q = 0; q < 20; ++q) {
            std::size_t e =
                std::uniform_int_distribution<std::size_t>(0, poly.num_edges() - 1)(rng);

            Exact y0 = poly.vertex(poly.edge(e).start_idx).y;
            Exact y1 = poly.vertex(poly.edge(e).end_idx).y;
            Exact qy_val = y0 + std::uniform_real_distribution<double>(0.01, 0.99)(rng) * (y1 - y0);

            SymbolicY qy{qy_val, n + 200 + static_cast<std::size_t>(q)};

            auto result = s.double_identify(e, qy, poly);
            assert(result.count <= 6);

            for (std::size_t k = 0; k < result.count; ++k) {
                std::size_t ai = result.arcs[k];
                const auto& a = s.arc(ai);
                assert((a.covers(e, LEFT, s.start_vertex, s.end_vertex) ||
                        a.covers(e, RIGHT, s.start_vertex, s.end_vertex)) &&
                       "double_identify returned arc that doesn't "
                       "contain the queried edge");
            }

            assert(result.count >= 1 && "double_identify should find at least one arc "
                                        "for any edge in the polygon");
        }
    }
    std::printf("  [PASS] test_double_identify\n");
}

static void test_compact(std::mt19937& rng, int iters) {
    std::printf("  test_compact (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords =
            std::uniform_int_distribution<std::size_t>(1, std::min(n - 2, std::size_t(6)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        std::size_t to_remove =
            std::uniform_int_distribution<std::size_t>(0, s.num_live_chords())(rng);
        std::vector<std::size_t> removable;
        for (std::size_t ci = 0; ci < s.num_chords(); ++ci)
            if (!s.chord(ci).dead)
                removable.push_back(ci);
        std::shuffle(removable.begin(), removable.end(), rng);
        for (std::size_t i = 0; i < to_remove && i < removable.size(); ++i)
            s.remove_chord(removable[i], poly);

        std::size_t live_nodes = s.num_live_nodes();
        std::size_t live_chords = s.num_live_chords();
        std::size_t live_arcs = s.num_live_arcs();

        std::vector<std::pair<std::size_t, std::size_t>> pre_weights;
        for (std::size_t i = 0; i < s.num_nodes(); ++i) {
            if (!s.node(i).dead)
                pre_weights.push_back({i, brute_region_weight(s, i)});
        }

        s.compact();

        assert(s.num_nodes() == live_nodes);
        assert(s.num_chords() == live_chords);
        assert(s.num_arcs() == live_arcs);

        for (std::size_t i = 0; i < s.num_nodes(); ++i)
            assert(!s.node(i).dead);
        for (std::size_t i = 0; i < s.num_chords(); ++i)
            assert(!s.chord(i).dead);
        for (std::size_t i = 0; i < s.num_arcs(); ++i)
            assert(!s.arc(i).dead);

        s.check_invariants();

        for (std::size_t i = 0; i < s.num_nodes(); ++i)
            assert(s.region_weight(i) == brute_region_weight(s, i));
        for (std::size_t i = 0; i < s.num_arcs(); ++i) {
            assert(s.arc(i).region_node < s.num_nodes());
        }
    }
    std::printf("  [PASS] test_compact\n");
}

static std::size_t td_depth(const TreeDecomposition& td, std::size_t idx) {
    if (td.node(idx).is_leaf())
        return 0;
    std::size_t ld = td_depth(td, td.node(idx).left_child);
    std::size_t rd = td_depth(td, td.node(idx).right_child);
    return 1 + std::max(ld, rd);
}

static std::size_t td_count_leaves(const TreeDecomposition& td, std::size_t idx) {
    if (td.node(idx).is_leaf())
        return 1;
    return td_count_leaves(td, td.node(idx).left_child) +
           td_count_leaves(td, td.node(idx).right_child);
}

static std::size_t td_count_internal(const TreeDecomposition& td, std::size_t idx) {
    if (td.node(idx).is_leaf())
        return 0;
    return 1 + td_count_internal(td, td.node(idx).left_child) +
           td_count_internal(td, td.node(idx).right_child);
}

static void test_tree_decomposition(std::mt19937& rng, int iters) {
    std::printf("  test_tree_decomposition (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 30)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords =
            std::uniform_int_distribution<std::size_t>(0, std::min(n - 2, std::size_t(10)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);
        s.compact();

        assert(s.is_conformal());
        s.build_tree_decomposition();
        const auto& td = s.tree_decomposition();

        assert(!td.empty());

        std::size_t num_regions = s.num_nodes();
        std::size_t num_chords_live = s.num_chords();

        std::size_t internals = td_count_internal(td, td.root());
        std::size_t leaves = td_count_leaves(td, td.root());
        assert(internals == num_chords_live &&
               "tree decomposition internals must biject with chords");
        assert(leaves == num_regions && "tree decomposition leaves must biject with regions");

        assert(td.size() == num_chords_live + num_regions);

        if (num_regions > 1) {
            std::size_t depth = td_depth(td, td.root());
            std::size_t depth_bound = 1;
            for (std::size_t remaining = num_regions; remaining > 1; ++depth_bound)
                remaining -= remaining / 4 + (remaining % 4 != 0);
            assert(depth <= depth_bound && "tree decomposition depth exceeds O(log m) bound");
        }

        assert(td.node(td.root()).parent == NONE);
        for (std::size_t i = 0; i < td.size(); ++i) {
            const auto& node = td.node(i);
            if (node.is_internal()) {
                assert(td.node(node.left_child).parent == i);
                assert(td.node(node.right_child).parent == i);
            }
        }

        std::set<std::size_t> leaf_regions;
        for (std::size_t i = 0; i < td.size(); ++i) {
            if (td.node(i).is_leaf()) {
                assert(td.node(i).region_idx < s.num_nodes());
                assert(leaf_regions.insert(td.node(i).region_idx).second &&
                       "duplicate region in tree decomposition leaves");
            }
        }
        assert(leaf_regions.size() == num_regions);

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

static void test_granularity(std::mt19937& rng, int iters) {
    std::printf("  test_granularity (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords =
            std::uniform_int_distribution<std::size_t>(0, std::min(n - 2, std::size_t(6)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        for (std::size_t granularity = 0; granularity <= n + 5; ++granularity) {
            bool semi = s.is_semigranular(granularity);
            bool gran = s.is_granular(granularity, poly);

            if (gran) {
                assert(semi && "granular must imply semigranular");
            }

            assert(semi == brute_is_semigranular(s, granularity));

            if (s.num_live_chords() == 0 && semi) {
                assert(gran && "no-chord semigranular must be granular");
            }

            if (semi && s.num_live_chords() > 0) {
                bool condition_ii = true;
                for (std::size_t ci = 0; ci < s.num_chords(); ++ci) {
                    if (s.chord(ci).dead)
                        continue;
                    const auto& c = s.chord(ci);
                    std::size_t d0 = s.node(c.region[0]).degree();
                    std::size_t d1 = s.node(c.region[1]).degree();
                    if (d0 >= 3 && d1 >= 3)
                        continue;
                    if (s.simulated_contraction_weight(ci, poly) <= granularity) {
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

static void test_simulated_contraction_nonvertex(std::mt19937& rng, int iters) {
    std::printf("  test_simulated_contraction_nonvertex (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 20)(rng);
        auto poly = random_polygon(rng, n);
        auto s = build_nonvertex_chord_submap(poly, rng);

        if (s.num_live_chords() == 0)
            continue;

        std::size_t simulated = s.simulated_contraction_weight(0, poly);

        Submap copy = s;
        std::size_t survivor = copy.remove_chord(0, poly);
        std::size_t actual = brute_region_weight(copy, survivor);

        assert(simulated == actual &&
               "simulated_contraction_weight != actual for non-vertex chord");
    }
    std::printf("  [PASS] test_simulated_contraction_nonvertex\n");
}

static void test_full_pipeline(std::mt19937& rng, int iters) {
    std::printf("  test_full_pipeline (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(5, 30)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords =
            std::uniform_int_distribution<std::size_t>(1, std::min(n - 2, std::size_t(10)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        s.assert_tree_property();
        s.check_invariants();
        assert(s.is_conformal());

        std::size_t to_remove =
            std::uniform_int_distribution<std::size_t>(0, s.num_live_chords() / 2 + 1)(rng);
        std::vector<std::size_t> removable;
        for (std::size_t ci = 0; ci < s.num_chords(); ++ci)
            if (!s.chord(ci).dead)
                removable.push_back(ci);
        std::shuffle(removable.begin(), removable.end(), rng);
        for (std::size_t i = 0; i < to_remove && i < removable.size(); ++i) {
            if (s.chord(removable[i]).dead)
                continue;
            s.remove_chord(removable[i], poly);
            s.assert_tree_property();
        }

        s.compact();
        s.check_invariants();
        assert(s.is_conformal());

        s.normalize(poly);
        s.check_invariants(poly);
        for (std::size_t i = 0; i < s.num_nodes(); ++i)
            assert(s.region_weight(i) == brute_region_weight(s, i));
        for (std::size_t i = 0; i < s.num_arcs(); ++i) {
            assert(s.arc(i).region_node < s.num_nodes());
        }

        s.build_tree_decomposition();
        const auto& td = s.tree_decomposition();
        assert(!td.empty());
        assert(td_count_leaves(td, td.root()) == s.num_nodes());
        assert(td_count_internal(td, td.root()) == s.num_chords());

        for (std::size_t v = 0; v < poly.num_vertices(); ++v) {
            if (v >= poly.num_edges())
                continue;
            SymbolicY qy = symbolic_y_of(poly.vertex(v));
            auto result = s.double_identify(v, qy, poly);
            assert(result.count >= 1 && result.count <= 6);
        }

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

static void test_edge_cases() {
    std::printf("  test_edge_cases...\n");
    std::fflush(stdout);

    {
        Polygon p({{0, 0, 0}, {1, 1, 1}});
        assert(p.num_edges() == 1);
        assert(p.is_endpoint(0));
        assert(p.is_endpoint(1));

        Submap s;
        s.add_node();
        s.start_vertex = 0;
        s.end_vertex = 1;

        Arc a{};
        a.first_edge = 0;
        a.last_edge = 0;
        a.first_side = LEFT;
        a.last_side = RIGHT;
        a.region_node = 0;
        a.edge_count = 1;
        std::size_t ai0 = s.add_arc(a);
        assert(s.start_arc == ai0 && s.end_arc == ai0);

        s.check_invariants();
        assert(s.region_weight(0) == 1);
        assert(s.is_conformal());
        assert(s.is_semigranular(1));
        assert(s.is_granular(1, p));

        s.build_tree_decomposition();
        assert(s.tree_decomposition().size() == 1);

        auto r = s.double_identify(0, {0.0, 0}, p);
        assert(r.count == 1);
    }

    {
        Polygon p({{0, 0, 0}, {1, 5, 1}, {2, 2, 2}});
        assert(p.num_edges() == 2);
        assert(p.is_y_extremum(1));

        Submap s;
        s.add_node();
        s.start_vertex = 0;
        s.end_vertex = 2;
        Arc a{};
        a.first_edge = 0;
        a.last_edge = 0;
        a.first_side = LEFT;
        a.last_side = RIGHT;
        a.region_node = 0;
        a.edge_count = 2;
        std::size_t ai0 = s.add_arc(a);
        assert(s.start_arc == ai0 && s.end_arc == ai0);

        s.check_invariants();
        assert(s.region_weight(0) == 2);
    }

    {
        Polygon p({{0, 0, 0}, {0, 0, 1}, {1, 1, 2}});
        assert(p.num_edges() == 2);
        assert(p.count_nonnull_edges(0, 0) == 0);
        assert(p.count_nonnull_edges(1, 1) == 1);
        assert(p.count_nonnull_edges(0, 1) == 1);
    }

    std::printf("  [PASS] test_edge_cases\n");
}

static void test_stress(std::mt19937& rng, int iters) {
    std::printf("  test_stress (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(10, 50)(rng);
        auto poly = random_polygon(rng, n);

        auto s = build_vertex_chord_submap(poly, n - 2, rng);

        s.assert_tree_property();
        s.check_invariants();
        assert(s.is_conformal());

        std::vector<std::size_t> all_chords;
        all_chords.reserve(s.num_chords());
        for (std::size_t ci = 0; ci < s.num_chords(); ++ci)
            all_chords.push_back(ci);
        std::shuffle(all_chords.begin(), all_chords.end(), rng);

        for (std::size_t ci : all_chords) {
            if (s.chord(ci).dead)
                continue;

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

        s.build_tree_decomposition();
        assert(s.tree_decomposition().size() == 1);
    }
    std::printf("  [PASS] test_stress\n");
}

static void test_null_length_chord(std::mt19937& rng, int iters) {
    std::printf("  test_null_length_chord (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(5, 15)(rng);
        auto poly = random_polygon(rng, n);

        std::vector<std::size_t> extrema;
        for (std::size_t v = 1; v + 1 < n; ++v) {
            if (poly.is_y_extremum(v))
                extrema.push_back(v);
        }
        if (extrema.empty())
            continue;

        std::size_t ext_v =
            extrema[std::uniform_int_distribution<std::size_t>(0, extrema.size() - 1)(rng)];

        Submap s;
        std::size_t r0 = s.add_node();
        std::size_t r1 = s.add_node();
        s.start_vertex = 0;
        s.end_vertex = n - 1;

        Arc a{};

        a.first_edge = ext_v;
        a.last_edge = ext_v;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = r1;
        a.edge_count = 0;
        std::size_t N = s.add_arc(a);

        std::size_t ne = poly.num_edges();
        a = {};
        a.first_edge = ext_v;
        a.last_edge = ext_v - 1;
        a.first_side = LEFT;
        a.last_side = LEFT;
        a.region_node = r0;

        a.edge_count = poly.count_nonnull_edges(ext_v, ne - 1) +
                       poly.count_nonnull_edges(0, ne - 1) + poly.count_nonnull_edges(0, ext_v - 1);
        std::size_t W = s.add_arc(a);
        assert(s.start_arc == W && s.end_arc == W &&
               "[C91 §2.4 tex 142]: the double-wrap arc is both endpoint arcs");

        Chord c{};
        c.region[0] = r0;
        c.region[1] = r1;
        c.left_edge = ext_v - 1;
        c.right_edge = ext_v - 1;
        c.left_side = LEFT;
        c.right_side = LEFT;
        c.y = poly.vertex(ext_v).y;
        c.y_tag = ext_v;
        c.is_null_length = true;
        c.left_adj = {{W}, 1};
        c.right_adj = {{N}, 1};
        s.add_chord(c);

        s.assert_tree_property();
        assert(s.region_weight(r1) == 0);

        assert(brute_region_weight(s, r0) > 0);
        assert(s.region_weight(r0) == brute_region_weight(s, r0));
        assert(s.is_conformal());

        s.remove_chord(0, poly);
        s.assert_tree_property();
        assert(s.num_live_nodes() == 1);

        s.compact();
        s.check_invariants();
    }
    std::printf("  [PASS] test_null_length_chord\n");
}

static void test_compact_idempotent(std::mt19937& rng, int iters) {
    std::printf("  test_compact_idempotent (%d iters)...\n", iters);
    for (int iter = 0; iter < iters; ++iter) {
        std::size_t n = std::uniform_int_distribution<std::size_t>(4, 15)(rng);
        auto poly = random_polygon(rng, n);
        std::size_t max_chords =
            std::uniform_int_distribution<std::size_t>(1, std::min(n - 2, std::size_t(5)))(rng);
        auto s = build_vertex_chord_submap(poly, max_chords, rng);

        for (std::size_t ci = 0; ci < s.num_chords(); ++ci) {
            if (s.chord(ci).dead)
                continue;
            if (std::uniform_int_distribution<int>(0, 1)(rng))
                s.remove_chord(ci, poly);
        }

        s.compact();
        std::size_t n1 = s.num_nodes();
        std::size_t c1 = s.num_chords();
        std::size_t a1 = s.num_arcs();
        s.check_invariants();

        s.compact();
        assert(s.num_nodes() == n1);
        assert(s.num_chords() == c1);
        assert(s.num_arcs() == a1);
        s.check_invariants();
    }
    std::printf("  [PASS] test_compact_idempotent\n");
}

int main(int argc, char** argv) {
    unsigned seed = 42;
    if (argc > 1)
        seed = static_cast<unsigned>(std::atoi(argv[1]));
    std::mt19937 rng(seed);

    std::printf("[C91 §2.0–2.4 e2e tests (seed=%u)]:\n", seed);
    std::fflush(stdout);

    std::setbuf(stdout, nullptr);
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
