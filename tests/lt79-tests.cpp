#include "algorithm/separator/planar_separator.h"
#include "support/random.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <queue>
#include <vector>

using namespace chazelle;
using chazelle::test::DeterministicRandomGenerator;

static EmbeddedPlanarGraph grid_graph(std::size_t k) {
    auto id = [&](std::size_t r, std::size_t c) { return r * k + c; };
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    std::vector<std::vector<std::size_t>> he(k, std::vector<std::size_t>(k)),
        ve(k, std::vector<std::size_t>(k));
    for (std::size_t r = 0; r < k; ++r)
        for (std::size_t c = 0; c + 1 < k; ++c) {
            he[r][c] = edges.size();
            edges.push_back({id(r, c), id(r, c + 1)});
        }
    for (std::size_t r = 0; r + 1 < k; ++r)
        for (std::size_t c = 0; c < k; ++c) {
            ve[r][c] = edges.size();
            edges.push_back({id(r, c), id(r + 1, c)});
        }

    std::vector<std::vector<std::size_t>> rot(k * k);
    for (std::size_t r = 0; r < k; ++r)
        for (std::size_t c = 0; c < k; ++c) {
            auto& v = rot[id(r, c)];
            if (c + 1 < k)
                v.push_back(he[r][c]);
            if (r + 1 < k)
                v.push_back(ve[r][c]);
            if (c > 0)
                v.push_back(he[r][c - 1]);
            if (r > 0)
                v.push_back(ve[r - 1][c]);
        }
    return EmbeddedPlanarGraph(k * k, std::move(edges), rot);
}

static EmbeddedPlanarGraph path_graph(std::size_t n) {
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    std::vector<std::vector<std::size_t>> rot(n);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        edges.push_back({i, i + 1});
        rot[i].push_back(i);
        rot[i + 1].push_back(i);
    }
    return EmbeddedPlanarGraph(n, std::move(edges), rot);
}

static EmbeddedPlanarGraph cycle_graph(std::size_t n) {
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    edges.reserve(n);
    std::vector<std::vector<std::size_t>> rot(n);
    for (std::size_t i = 0; i < n; ++i)
        edges.push_back({i, (i + 1) % n});
    for (std::size_t i = 0; i < n; ++i) {
        rot[i].push_back(i);
        rot[i].push_back((i + n - 1) % n);
    }
    return EmbeddedPlanarGraph(n, std::move(edges), rot);
}

static EmbeddedPlanarGraph random_tree(std::size_t n, unsigned seed) {
    DeterministicRandomGenerator rng(seed);
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    std::vector<std::vector<std::size_t>> rot(n);
    for (std::size_t i = 1; i < n; ++i) {
        std::size_t j = rng.next() % i;
        rot[i].push_back(edges.size());
        rot[j].push_back(edges.size());
        edges.push_back({j, i});
    }
    return EmbeddedPlanarGraph(n, std::move(edges), rot);
}

static EmbeddedPlanarGraph wheel_graph(std::size_t n) {
    assert(n >= 4);
    std::size_t m = n - 1;
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    edges.reserve(2 * m);
    std::vector<std::vector<std::size_t>> rot(n);
    for (std::size_t i = 0; i < m; ++i)
        edges.push_back({1 + i, 1 + (i + 1) % m});
    for (std::size_t i = 0; i < m; ++i)
        edges.push_back({0, 1 + i});
    for (std::size_t i = 0; i < m; ++i)
        rot[0].push_back(m + i);
    for (std::size_t i = 0; i < m; ++i) {
        rot[1 + i].push_back(i);
        rot[1 + i].push_back(m + i);
        rot[1 + i].push_back((i + m - 1) % m);
    }
    return EmbeddedPlanarGraph(n, std::move(edges), rot);
}

static EmbeddedPlanarGraph diag_grid_graph(std::size_t k, unsigned seed) {
    DeterministicRandomGenerator rng(seed);
    auto id = [&](std::size_t r, std::size_t c) { return r * k + c; };
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    std::vector<std::vector<std::pair<double, std::size_t>>> ang(k * k);
    auto add = [&](std::size_t a, std::size_t b, double dxa, double dya) {
        std::size_t e = edges.size();
        edges.push_back({a, b});
        ang[a].push_back({std::atan2(dya, dxa), e});
        ang[b].push_back({std::atan2(-dya, -dxa), e});
    };
    for (std::size_t r = 0; r < k; ++r)
        for (std::size_t c = 0; c < k; ++c) {
            if (c + 1 < k)
                add(id(r, c), id(r, c + 1), 1, 0);
            if (r + 1 < k)
                add(id(r, c), id(r + 1, c), 0, 1);
            if (r + 1 < k && c + 1 < k) {
                if (rng.next() & 1)
                    add(id(r, c), id(r + 1, c + 1), 1, 1);
                else
                    add(id(r, c + 1), id(r + 1, c), -1, 1);
            }
        }
    std::vector<std::vector<std::size_t>> rot(k * k);
    for (std::size_t v = 0; v < k * k; ++v) {
        std::sort(ang[v].begin(), ang[v].end());
        for (auto& p : ang[v])
            rot[v].push_back(p.second);
    }
    return EmbeddedPlanarGraph(k * k, std::move(edges), rot);
}

static void verify_partition(const EmbeddedPlanarGraph& g, const std::vector<SepPart>& part) {
    const std::size_t n = g.num_vertices();
    assert(part.size() == n);

    std::size_t na = 0, nb = 0, nd = 0;
    for (SepPart p : part) {
        if (p == SepPart::A)
            ++na;
        else if (p == SepPart::B)
            ++nb;
        else
            ++nd;
    }
    assert(na + nb + nd == n);

    assert(3 * na <= 2 * n && "(ii): |A| ≤ 2n/3");
    assert(3 * nb <= 2 * n && "(ii): |B| ≤ 2n/3");

    assert(nd * nd <= 8 * n && "(iii): |D| ≤ √(8n)");

    for (std::size_t e = 0; e < g.num_edges(); ++e) {
        SepPart pu = part[g.edge_u(e)], pv = part[g.edge_v(e)];
        assert(
            !((pu == SepPart::A && pv == SepPart::B) || (pu == SepPart::B && pv == SepPart::A)) &&
            "(i): no edge joins A and B");
    }
}

static void verify_decomposition(const EmbeddedPlanarGraph& g, const SeparatorDecomposition& d) {
    const std::size_t mu = g.num_vertices();
    assert(d.subset.size() == mu);

    std::vector<std::size_t> sizes(d.num_subsets, 0);
    std::size_t dstar = 0;
    for (std::size_t v = 0; v < mu; ++v) {
        if (d.subset[v] == NONE)
            ++dstar;
        else {
            assert(d.subset[v] < d.num_subsets);
            ++sizes[d.subset[v]];
        }
    }
    assert(dstar == d.dstar_size);
    __extension__ typedef unsigned __int128 u128;
    for (std::size_t s : sizes)
        assert((u128)s * s * s <= (u128)mu * mu && "[C91 §3.4 tex 304]: each |D_i| ≤ μ^{2/3}");

    assert((long double)dstar <= 19.0L * powl((long double)mu, 2.0L / 3.0L) + 1.0L &&
           "[C91 §3.4 tex 304]: |D*| = O(μ^{2/3})");

    std::vector<bool> seen(mu, false);
    for (std::size_t v0 = 0; v0 < mu; ++v0) {
        if (d.subset[v0] == NONE || seen[v0])
            continue;
        seen[v0] = true;
        std::queue<std::size_t> q;
        q.push(v0);
        while (!q.empty()) {
            std::size_t v = q.front();
            q.pop();
            for (std::size_t h : g.incident_halves(v)) {
                std::size_t w = g.half_to(h);
                if (d.subset[w] == NONE || seen[w])
                    continue;
                assert(d.subset[w] == d.subset[v0] &&
                       "[C91 §3.4 tex 304]: subsets are separated by D*");
                seen[w] = true;
                q.push(w);
            }
        }
    }
}

static void run_family(const EmbeddedPlanarGraph& g) {
    verify_partition(g, planar_separator(g));
    verify_decomposition(g, build_separator_decomposition(g));
}

static void test_root_inside_cycle() {
    const std::vector<std::pair<std::size_t, std::size_t>> edges = {
        {0, 12}, {0, 10}, {8, 0}, {6, 0},  {0, 9}, {2, 0},  {5, 0},
        {1, 11}, {12, 1}, {2, 3}, {4, 10}, {6, 7}, {12, 11}};
    const std::vector<std::vector<std::size_t>> rotations = {{0, 1, 2, 3, 4, 5, 6},
                                                             {7, 8},
                                                             {9, 5},
                                                             {9},
                                                             {10},
                                                             {6},
                                                             {11, 3},
                                                             {11},
                                                             {2},
                                                             {4},
                                                             {10, 1},
                                                             {7, 12},
                                                             {0, 12, 8}};
    run_family(EmbeddedPlanarGraph(13, edges, rotations));
    std::printf("  [PASS] root_inside_cycle\n");
}

static void test_grids() {
    for (std::size_t k : {1u, 2u, 3u, 4u, 5u, 7u, 10u, 16u, 25u, 31u})
        run_family(grid_graph(k));
    std::printf("  [PASS] grids\n");
}

static void test_paths_and_cycles() {
    for (std::size_t n : {1u, 2u, 3u, 5u, 17u, 100u, 501u})
        run_family(path_graph(n));
    for (std::size_t n : {3u, 4u, 9u, 64u, 200u})
        run_family(cycle_graph(n));
    std::printf("  [PASS] paths_and_cycles\n");
}

static void test_random_trees() {
    for (unsigned seed = 0; seed < 50; ++seed)
        run_family(random_tree(200 + 13 * seed, seed));
    for (unsigned seed = 0; seed < 8; ++seed)
        run_family(random_tree(2000 + 313 * seed, seed ^ 0xbeefu));
    std::printf("  [PASS] random_trees\n");
}

static void test_wheels() {
    for (std::size_t n : {4u, 5u, 10u, 100u, 333u})
        run_family(wheel_graph(n));
    std::printf("  [PASS] wheels\n");
}

static void test_diag_grids() {
    for (unsigned seed = 0; seed < 30; ++seed)
        run_family(diag_grid_graph(6 + seed % 20, seed));
    std::printf("  [PASS] diag_grids\n");
}

static void test_disconnected() {
    {
        std::vector<std::pair<std::size_t, std::size_t>> edges;
        std::vector<std::vector<std::size_t>> rot(30);
        for (std::size_t c = 0; c < 6; ++c)
            for (std::size_t i = 0; i + 1 < 5; ++i) {
                std::size_t a = 5 * c + i;
                rot[a].push_back(edges.size());
                rot[a + 1].push_back(edges.size());
                edges.push_back({a, a + 1});
            }
        EmbeddedPlanarGraph g(30, std::move(edges), rot);
        auto part = planar_separator(g);
        verify_partition(g, part);
        std::size_t nd = 0;
        for (SepPart p : part)
            if (p == SepPart::D)
                ++nd;
        assert(nd == 0 && "[LT79 tex 360–365]: light components pack with D = ∅");
        verify_decomposition(g, build_separator_decomposition(g));
    }
    {
        std::vector<std::pair<std::size_t, std::size_t>> edges;
        std::vector<std::vector<std::size_t>> rot(60);
        for (std::size_t i = 0; i + 1 < 50; ++i) {
            rot[i].push_back(edges.size());
            rot[i + 1].push_back(edges.size());
            edges.push_back({i, i + 1});
        }
        for (std::size_t i = 50; i + 1 < 60; ++i) {
            rot[i].push_back(edges.size());
            rot[i + 1].push_back(edges.size());
            edges.push_back({i, i + 1});
        }
        EmbeddedPlanarGraph g(60, std::move(edges), rot);
        auto part = planar_separator(g);
        verify_partition(g, part);
        for (std::size_t v = 50; v < 60; ++v)
            assert(part[v] != SepPart::D && "[LT79 tex 378–380]: D lies in one component");
        verify_decomposition(g, build_separator_decomposition(g));
    }
    {
        std::vector<std::pair<std::size_t, std::size_t>> edges;
        std::vector<std::vector<std::size_t>> rot(59);
        auto add = [&](std::size_t a, std::size_t b) {
            rot[a].push_back(edges.size());
            rot[b].push_back(edges.size());
            edges.push_back({a, b});
        };
        std::size_t next = 1;
        std::vector<std::size_t> l1v, l3v;
        for (int i = 0; i < 14; ++i) {
            l1v.push_back(next);
            add(0, next++);
        }
        std::size_t x = next++;
        add(l1v[0], x);
        for (int i = 0; i < 14; ++i) {
            l3v.push_back(next);
            add(x, next++);
        }
        std::size_t y = next++;
        add(l3v[0], y);
        for (int i = 0; i < 13; ++i)
            add(y, next++);
        assert(next == 44);
        for (std::size_t c = 0; c < 3; ++c) {
            for (std::size_t i = 0; i + 1 < 5; ++i)
                add(next + i, next + i + 1);
            next += 5;
        }
        assert(next == 59);
        EmbeddedPlanarGraph g(59, std::move(edges), rot);
        auto part = planar_separator(g);
        verify_partition(g, part);
        verify_decomposition(g, build_separator_decomposition(g));
    }
    std::printf("  [PASS] disconnected\n");
}

static void test_induced_subgraph() {
    auto g = grid_graph(4);
    std::vector<std::size_t> nodes = {0, 1, 2, 4, 5, 6, 9, 10};
    std::vector<std::size_t> ids;
    EmbeddedPlanarGraph sub = g.induced(nodes, &ids);
    assert(sub.num_vertices() == nodes.size());
    assert(ids == nodes);

    std::size_t expect = 0;
    for (std::size_t e = 0; e < g.num_edges(); ++e) {
        bool u_in = std::find(nodes.begin(), nodes.end(), g.edge_u(e)) != nodes.end();
        bool v_in = std::find(nodes.begin(), nodes.end(), g.edge_v(e)) != nodes.end();
        if (u_in && v_in)
            ++expect;
    }
    assert(sub.num_edges() == expect);
    run_family(sub);
    std::printf("  [PASS] induced_subgraph\n");
}

int main() {
    std::setbuf(stdout, nullptr);
    std::printf("[LT79 tests]:\n");
    test_grids();
    test_paths_and_cycles();
    test_random_trees();
    test_wheels();
    test_diag_grids();
    test_disconnected();
    test_induced_subgraph();
    test_root_inside_cycle();
    std::printf("All LT79 tests passed.\n");
    return 0;
}
