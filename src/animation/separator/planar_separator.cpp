#include "planar_separator.h"
#include "../polygon/exact.h"
#include "../trace.h"

#include <algorithm>

namespace chazelle::animation {

EmbeddedPlanarGraph::EmbeddedPlanarGraph(std::size_t num_vertices,
                                         std::vector<std::pair<std::size_t, std::size_t>> edges,
                                         const std::vector<std::vector<std::size_t>>& rotations)
    : num_vertices_(num_vertices), edges_(std::move(edges)) {
    assert(rotations.size() == num_vertices_ &&
           "[LT79 §3 tex 524–528]: one rotation list per vertex");

    nxt_.assign(2 * edges_.size(), NONE);
    prv_.assign(2 * edges_.size(), NONE);
    rot_.assign(num_vertices_, {});

    std::vector<std::size_t> seen(edges_.size(), 0);

    for (std::size_t v = 0; v < num_vertices_; ++v) {
        rot_[v].reserve(rotations[v].size());
        for (std::size_t e : rotations[v]) {
            assert(e < edges_.size() && "rotation refers to unknown edge");
            assert((edges_[e].first == v || edges_[e].second == v) &&
                   "rotation lists an edge not incident to its vertex");
            assert(edges_[e].first != edges_[e].second &&
                   "[LT79]: loops are excluded from the representation");
            std::size_t h = (edges_[e].first == v) ? 2 * e : 2 * e + 1;
            rot_[v].push_back(h);
            ++seen[e];
        }

        const auto& r = rot_[v];
        for (std::size_t i = 0; i < r.size(); ++i) {
            std::size_t j = (i + 1 == r.size()) ? 0 : i + 1;
            nxt_[r[i]] = r[j];
            prv_[r[j]] = r[i];
        }
    }

    for (std::size_t e = 0; e < edges_.size(); ++e) {
        assert(seen[e] == 2 && "[LT79 §3 tex 524–528]: each edge appears exactly once in "
                               "each endpoint's rotation");
        (void)e;
    }
    (void)seen;
}

EmbeddedPlanarGraph EmbeddedPlanarGraph::induced(const std::vector<std::size_t>& nodes,
                                                 std::vector<std::size_t>* old_index) const {
    std::vector<std::size_t> new_index(num_vertices_, NONE);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        assert(nodes[i] < num_vertices_);
        assert(new_index[nodes[i]] == NONE && "duplicate node in induced()");
        new_index[nodes[i]] = i;
    }

    std::vector<std::pair<std::size_t, std::size_t>> new_edges;
    std::vector<std::size_t> edge_map(edges_.size(), NONE);
    std::vector<std::vector<std::size_t>> new_rot(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        std::size_t v = nodes[i];
        for (std::size_t h : rot_[v]) {
            std::size_t w = half_to(h);
            if (new_index[w] == NONE)
                continue;
            std::size_t e = h >> 1;
            if (edge_map[e] == NONE) {
                edge_map[e] = new_edges.size();
                new_edges.emplace_back(new_index[edges_[e].first], new_index[edges_[e].second]);
            }
            new_rot[i].push_back(edge_map[e]);
        }
    }

    if (old_index)
        *old_index = nodes;
    return EmbeddedPlanarGraph(nodes.size(), std::move(new_edges), new_rot);
}

namespace {

struct HalfEdgeGraph {
    std::vector<std::size_t> to;
    std::vector<std::size_t> nxt, prv;
    std::vector<std::size_t> first_half;
    std::size_t n = 0;

    std::size_t num_edges() const {
        return to.size() / 2;
    }
    std::size_t origin(std::size_t h) const {
        return to[h ^ 1];
    }

    std::size_t add_vertex() {
        first_half.push_back(NONE);
        return n++;
    }

    std::size_t raw_edge(std::size_t u, std::size_t v) {
        assert(u != v && "[LT79]: loops never arise");
        (void)u;
        to.push_back(v);
        to.push_back(u);
        nxt.push_back(NONE);
        nxt.push_back(NONE);
        prv.push_back(NONE);
        prv.push_back(NONE);
        return to.size() - 2;
    }

    void insert_before(std::size_t h, std::size_t pos) {
        std::size_t p = prv[pos];
        nxt[p] = h;
        prv[h] = p;
        nxt[h] = pos;
        prv[pos] = h;
    }
    void insert_after(std::size_t h, std::size_t pos) {
        std::size_t nx = nxt[pos];
        nxt[pos] = h;
        prv[h] = pos;
        nxt[h] = nx;
        prv[nx] = h;
    }

    template <typename F> void for_out_halves(std::size_t v, F&& f) const {
        std::size_t h0 = first_half[v];
        if (h0 == NONE)
            return;
        std::size_t h = h0;
        do {
            f(h);
            h = nxt[h];
        } while (h != h0);
    }

    std::size_t face_next(std::size_t h) const {
        return nxt[h ^ 1];
    }

    std::size_t count_faces() const {
        std::vector<bool> vis(to.size(), false);
        std::size_t faces = 0;
        for (std::size_t h = 0; h < to.size(); ++h) {
            if (vis[h])
                continue;
            ++faces;
            std::size_t x = h;
            do {
                vis[x] = true;
                x = face_next(x);
            } while (x != h);
        }
        return faces;
    }

    void assert_euler() const {
#ifndef NDEBUG
        assert(n >= 1 && !to.empty());
        assert(n + count_faces() == num_edges() + 2 &&
               "[LT79 §3 tex 534–538]: rotation system must be a planar "
               "(sphere) embedding");
#endif
    }
};

struct ReducedPlanarGraph {
    HalfEdgeGraph g;
    std::size_t apex = 0;
    std::vector<std::size_t> orig;
    std::vector<std::size_t> parent;
    std::vector<std::size_t> parent_edge;
    std::vector<std::size_t> depth;
    std::vector<std::size_t> desc;
    std::vector<std::size_t> bfs_order;
    std::size_t total = 0;

    std::size_t cost(std::size_t v) const {
        return v == apex ? 0 : 1;
    }
    bool is_tree_edge(std::size_t e) const {
        std::size_t u = g.origin(2 * e), w = g.to[2 * e];
        return parent_edge[u] == e || parent_edge[w] == e;
    }

    std::size_t tree_half(std::size_t a, std::size_t b) const {
        std::size_t te = (parent[a] == b) ? parent_edge[a] : parent_edge[b];
        assert(te != NONE && (parent[a] == b || parent[b] == a) &&
               "tree_half requires a tree edge");
        return (g.origin(2 * te) == a) ? 2 * te : 2 * te + 1;
    }
};

struct SeparatorCycle {
    std::vector<bool> on_v;
    std::vector<bool> on_e;
    std::vector<std::size_t> cyc_next, cyc_prev;
    std::vector<std::size_t> cyc_half;
    std::size_t nontree = NONE;
    std::size_t inside_half = NONE;

    std::size_t inside_cost = 0;
    std::size_t cycle_count = 0;
    std::size_t lca = NONE;

    std::size_t traversal_half(const HalfEdgeGraph& g) const {
        std::size_t a = g.origin(2 * nontree);
        if (cyc_half[a] != NONE && (cyc_half[a] >> 1) == nontree)
            return cyc_half[a];
        std::size_t b = g.to[2 * nontree];
        assert(cyc_half[b] != NONE && (cyc_half[b] >> 1) == nontree);
        return cyc_half[b];
    }
};

template <typename F>
void rotation_between(const HalfEdgeGraph& g, std::size_t from_excl, std::size_t to_excl, F&& f) {
    if (from_excl == to_excl)
        return;
    for (std::size_t x = g.nxt[from_excl]; x != to_excl; x = g.nxt[x]) {
        assert(x != from_excl && "rotation walk must terminate");
        f(x);
    }
}

std::size_t hanging_cost(const ReducedPlanarGraph& s, std::size_t c, std::size_t w) {
    if (s.parent[w] == c)
        return s.desc[w];
    assert(s.parent[c] == w && "[LT79 §3 tex 606–609]: a crossing tree edge is a parent edge");
    return s.total - s.desc[c];
}

void triangulate(HalfEdgeGraph& g) {
    std::vector<bool> vis(g.to.size(), false);
    std::vector<std::size_t> multiplicity(g.n, 0);
    for (std::size_t h0 = 0; h0 < g.to.size(); ++h0) {
        if (vis[h0])
            continue;
        std::vector<std::size_t> walk;
        std::size_t x = h0;
        do {
            vis[x] = true;
            walk.push_back(x);
            x = g.face_next(x);
        } while (x != h0);

        std::size_t k = walk.size();
        if (k <= 3)
            continue;
        std::vector<std::size_t> nx(k), pv(k);
        std::vector<bool> live(k, true);
        std::vector<std::size_t> candidates;
        std::size_t distinct = 0;
        for (std::size_t i = 0; i < k; ++i) {
            nx[i] = (i + 1) % k;
            pv[i] = (i + k - 1) % k;
            if (multiplicity[g.origin(walk[i])]++ == 0)
                ++distinct;
            candidates.push_back(i);
        }
        assert(distinct >= 3 && "[LT79 Step 7 tex 595]: a simple connected graph with "
                                "at least three vertices admits triangular faces");
        std::size_t len = k;
        while (len > 3 && !candidates.empty()) {
            std::size_t i1 = candidates.back();
            candidates.pop_back();
            if (!live[i1])
                continue;
            std::size_t i2 = nx[i1];
            std::size_t h1 = walk[i1], h2 = walk[i2];
            std::size_t a = g.origin(h1), c = g.to[h2];
            const std::size_t b = g.to[h1];
            if (a == c || (distinct == 3 && multiplicity[b] == 1))
                continue;

            std::size_t p = g.raw_edge(a, c);
            g.insert_before(p, h1);
            g.insert_after(p ^ 1, h2 ^ 1);
            assert(g.face_next(h2) == (p ^ 1) && g.face_next(p ^ 1) == h1 &&
                   "ear cut must close the triangle (h1, h2, twin(p))");
            vis.push_back(true);
            vis.push_back(true);

            walk[i1] = p;
            nx[i1] = nx[i2];
            pv[nx[i2]] = i1;
            live[i2] = false;
            if (--multiplicity[b] == 0)
                --distinct;
            --len;
            candidates.push_back(pv[i1]);
            candidates.push_back(i1);
        }
        assert(len == 3 && distinct == 3 &&
               "[LT79 Step 7 tex 595–597]: every face must become a triangle");
        for (std::size_t h : walk)
            multiplicity[g.origin(h)] = 0;
    }
    g.assert_euler();
#ifndef NDEBUG
    for (std::size_t h = 0; h < g.to.size(); ++h) {
        const std::size_t h2 = g.face_next(h), h3 = g.face_next(h2);
        assert(h2 != h && h3 != h && g.face_next(h3) == h &&
               "[LT79 Step 7 tex 595–597]: three edges per face");
    }
#endif
}

void build_initial_cycle(const ReducedPlanarGraph& s, std::size_t e, SeparatorCycle& cs) {
    const HalfEdgeGraph& g = s.g;
    std::size_t v1 = g.origin(2 * e), w1 = g.to[2 * e];

    cs.on_v.assign(g.n, false);
    cs.on_e.assign(g.num_edges(), false);
    cs.cyc_next.assign(g.n, NONE);
    cs.cyc_prev.assign(g.n, NONE);
    cs.cyc_half.assign(g.n, NONE);

    std::vector<std::size_t> pv, pw;
    for (std::size_t x = v1; x != NONE; x = s.parent[x])
        pv.push_back(x);
    for (std::size_t x = w1; x != NONE; x = s.parent[x])
        pw.push_back(x);
    while (pv.size() >= 2 && pw.size() >= 2 && pv[pv.size() - 2] == pw[pw.size() - 2]) {
        pv.pop_back();
        pw.pop_back();
    }
    assert(pv.back() == pw.back() && "tree paths must meet at the LCA");
    cs.lca = pv.back();

    std::vector<std::size_t> cyc;
    cyc.push_back(v1);
    for (std::size_t i = 0; i + 1 < pw.size(); ++i)
        cyc.push_back(pw[i]);
    if (cs.lca != v1)
        cyc.push_back(cs.lca);
    for (std::size_t i = pv.size() - 1; i-- > 1;)
        cyc.push_back(pv[i]);

    cs.cycle_count = cyc.size();
    cs.nontree = e;
    for (std::size_t i = 0; i < cyc.size(); ++i) {
        std::size_t a = cyc[i];
        std::size_t b = cyc[(i + 1) % cyc.size()];
        cs.on_v[a] = true;
        cs.cyc_next[a] = b;
        cs.cyc_prev[b] = a;
        std::size_t half_ab = (i == 0) ? 2 * e : s.tree_half(a, b);
        assert(g.origin(half_ab) == a && g.to[half_ab] == b);
        cs.on_e[half_ab >> 1] = true;
        cs.cyc_half[a] = half_ab;
    }

    std::size_t cost_alpha = 0, cost_beta = 0, cycle_cost = 0;
    for (std::size_t i = 0; i < cyc.size(); ++i) {
        std::size_t c = cyc[i];
        cycle_cost += s.cost(c);
        std::size_t out = cs.cyc_half[c];
        std::size_t in = cs.cyc_half[cyc[(i + cyc.size() - 1) % cyc.size()]];
        assert(g.to[in] == c);
        auto tally = [&](std::size_t x, std::size_t& acc) {
            std::size_t w = g.to[x];
            if (cs.on_v[w])
                return;
            if (!s.is_tree_edge(x >> 1))
                return;
            acc += hanging_cost(s, c, w);
        };
        rotation_between(g, out, in ^ 1, [&](std::size_t x) { tally(x, cost_alpha); });
        rotation_between(g, in ^ 1, out, [&](std::size_t x) { tally(x, cost_beta); });
    }

    assert(cost_alpha + cost_beta + cycle_cost == s.total &&
           "[LT79 §3 tex 605–610]: the two sides plus the cycle must "
           "account for the whole graph");

    if (cost_alpha >= cost_beta) {
        cs.inside_cost = cost_alpha;
        cs.inside_half = 2 * e + 1;
    } else {
        cs.inside_cost = cost_beta;
        cs.inside_half = 2 * e;
    }
}

struct CycleSideScan {
    const ReducedPlanarGraph* s = nullptr;
    const SeparatorCycle* cs = nullptr;
    const std::vector<bool>* chain_mark = nullptr;
    std::size_t arc_from = NONE, arc_to = NONE;
    const std::vector<std::size_t>* chain = nullptr;
    std::size_t new_half = NONE;
    bool new_at_start = false;

    bool scan_alpha = false;

    bool done = false;
    std::size_t cost = 0;

    std::size_t cur = NONE;
    int seg = 0;
    std::size_t cpos = 0;
    std::size_t x = NONE, xstop = NONE;

    std::size_t out_half() const {
        if (seg == 0) {
            if (cur != arc_to)
                return cs->cyc_half[cur];
            if (chain->empty())
                return new_half;
            return new_at_start ? new_half : s->tree_half(cur, chain->front());
        }
        if (cpos + 1 < chain->size())
            return s->tree_half(cur, (*chain)[cpos + 1]);
        return new_at_start ? s->tree_half(cur, arc_from) : new_half;
    }
    std::size_t in_half() const {
        if (seg == 0) {
            if (cur == arc_from) {
                if (chain->empty())
                    return new_half;
                return new_at_start ? s->tree_half(chain->back(), arc_from) : new_half;
            }
            return cs->cyc_half[cs->cyc_prev[cur]];
        }
        if (cpos == 0)
            return new_at_start ? new_half : s->tree_half(arc_to, cur);
        return s->tree_half((*chain)[cpos - 1], cur);
    }

    void setup_vertex() {
        std::size_t o = out_half(), i = in_half();
        assert(s->g.origin(o) == cur && s->g.to[i] == cur);
        if (scan_alpha) {
            x = s->g.nxt[o];
            xstop = i ^ 1;
        } else {
            x = s->g.nxt[i ^ 1];
            xstop = o;
        }
        if (x == xstop)
            x = NONE;
    }

    void begin() {
        cur = arc_from;
        seg = 0;
        cpos = 0;
        setup_vertex();
    }

    bool advance_vertex() {
        if (seg == 0) {
            if (cur == arc_to) {
                if (chain->empty())
                    return false;
                seg = 1;
                cpos = 0;
                cur = chain->front();
            } else {
                cur = cs->cyc_next[cur];
            }
        } else {
            if (cpos + 1 >= chain->size())
                return false;
            ++cpos;
            cur = (*chain)[cpos];
        }
        setup_vertex();
        return true;
    }

    bool step() {
        if (done)
            return false;
        while (x == NONE) {
            if (!advance_vertex()) {
                done = true;
                return false;
            }
        }
        std::size_t w = s->g.to[x];
        if (!cs->on_v[w] && !(*chain_mark)[w] && s->is_tree_edge(x >> 1))
            cost += hanging_cost(*s, cur, w);
        x = s->g.nxt[x];
        if (x == xstop)
            x = NONE;
        return true;
    }
};

struct TreePathToCycle {
    std::size_t z = NONE;
    std::vector<std::size_t> chain;
};

TreePathToCycle find_path_to_cycle(const ReducedPlanarGraph& s, const SeparatorCycle& cs,
                                   std::size_t y, std::vector<bool>& markY,
                                   std::vector<bool>& markL) {
    TreePathToCycle r;
    if (cs.on_v[y]) {
        r.z = y;
        return r;
    }
    std::vector<std::size_t> listY{y}, listL{cs.lca};
    markY[y] = true;
    assert(cs.on_v[cs.lca] && "cycle lca must be maintained");

    bool y_alive = true, l_alive = true;
    std::size_t meetY = NONE, meetL = NONE;
    while (true) {
        if (y_alive) {
            std::size_t p = s.parent[listY.back()];
            if (p == NONE) {
                y_alive = false;
            } else if (cs.on_v[p]) {
                r.z = p;
                r.chain = std::move(listY);
                for (std::size_t v : r.chain)
                    markY[v] = false;
                for (std::size_t v : listL)
                    markL[v] = false;
                return r;
            } else if (markL[p]) {
                meetY = NONE;
                meetL = p;
                break;
            } else {
                listY.push_back(p);
                markY[p] = true;
            }
        }
        if (l_alive) {
            std::size_t q = s.parent[listL.back()];
            if (q == NONE) {
                l_alive = false;
            } else if (markY[q]) {
                meetY = q;
                meetL = NONE;
                break;
            } else {
                assert(!cs.on_v[q] && "the cycle lca's ancestors lie off the cycle");
                listL.push_back(q);
                markL[q] = true;
            }
        }
        assert((y_alive || l_alive) && "both ascents reach the root, so they must meet");
    }

    r.z = cs.lca;
    std::size_t lprime = (meetL != NONE) ? meetL : meetY;
    for (std::size_t v : listY) {
        r.chain.push_back(v);
        if (v == lprime)
            break;
    }
    if (r.chain.empty() || r.chain.back() != lprime) {
        assert(meetL != NONE);

        r.chain.assign(listY.begin(), listY.end());
        r.chain.push_back(lprime);
    }

    std::size_t cut = listL.size();
    for (std::size_t i = 0; i < listL.size(); ++i)
        if (listL[i] == lprime) {
            cut = i;
            break;
        }

    for (std::size_t i = cut; i-- > 1;)
        r.chain.push_back(listL[i]);

    for (std::size_t v : listY)
        markY[v] = false;
    for (std::size_t v : listL)
        markL[v] = false;
    return r;
}

void improve_cycle(const ReducedPlanarGraph& s, SeparatorCycle& cs, std::size_t total_units) {
    const HalfEdgeGraph& g = s.g;
    std::vector<bool> markY(g.n, false), markL(g.n, false);
    std::vector<bool> chain_mark(g.n, false);

#ifndef NDEBUG
    std::size_t guard = 2 * g.count_faces() + 8;
#endif

    while (3 * cs.inside_cost > 2 * total_units) {
        assert(guard > 0 && "[LT79 §3 tex 656–658]: Step 9 must terminate within the "
                            "face budget");
#ifndef NDEBUG
        --guard;
#endif

        std::size_t f1 = cs.inside_half;
        std::size_t f2 = g.face_next(f1);
        std::size_t f3 = g.face_next(f2);

        assert(f2 != f1 && f3 != f1 && g.face_next(f3) == f1 &&
               "[LT79 Step 9 tex 620]: locate the inside triangle "
               "guaranteed by Step 7");

        std::size_t vstar = g.origin(f1), wstar = g.to[f1];
        std::size_t y = g.to[f2];
        assert((f1 >> 1) == cs.nontree);
        assert(y != vstar && y != wstar &&
               "triangles of a loopless multigraph have distinct corners");

        bool e2_on = cs.on_e[f2 >> 1];
        bool e3_on = cs.on_e[f3 >> 1];
        assert(!(e2_on && e3_on) && "[LT79 tex 228–229] case 1: the face cannot BE the cycle "
                                    "while vertices remain inside");

        if (e2_on || e3_on) {
            std::size_t drop = e2_on ? wstar : vstar;
            std::size_t keep = e2_on ? vstar : wstar;
            std::size_t tri = e2_on ? f3 : f2;
            std::size_t new_e = tri >> 1;
            std::size_t dropped_cycle_e = e2_on ? (f2 >> 1) : (f3 >> 1);

            if (g.origin(cs.traversal_half(g)) == keep) {
                assert(cs.cyc_next[keep] == drop && cs.cyc_next[drop] == y);
                cs.cyc_next[keep] = y;
                cs.cyc_prev[y] = keep;
                cs.cyc_half[keep] = (g.origin(2 * new_e) == keep) ? 2 * new_e : 2 * new_e + 1;
            } else {
                assert(cs.cyc_next[y] == drop && cs.cyc_next[drop] == keep);
                cs.cyc_next[y] = keep;
                cs.cyc_prev[keep] = y;
                cs.cyc_half[y] = (g.origin(2 * new_e) == y) ? 2 * new_e : 2 * new_e + 1;
            }
            cs.on_v[drop] = false;
            cs.cyc_next[drop] = cs.cyc_prev[drop] = cs.cyc_half[drop] = NONE;
            cs.on_e[cs.nontree] = false;
            cs.on_e[dropped_cycle_e] = false;
            cs.on_e[new_e] = true;
            if (drop == cs.lca)
                cs.lca = y;
            cs.cycle_count -= 1;
            cs.nontree = new_e;
            cs.inside_half = tri ^ 1;
            continue;
        }

        bool e2_tree = s.is_tree_edge(f2 >> 1);
        bool e3_tree = s.is_tree_edge(f3 >> 1);
        assert(!(e2_tree && e3_tree) && "[LT79 tex 237–238] case 3a: two tree edges would close a "
                                        "tree cycle");

        if (e2_tree || e3_tree) {
            assert(!cs.on_v[y] && "a tree edge from the cycle to an on-cycle vertex "
                                  "would close a tree cycle");
            std::size_t attach = e2_tree ? wstar : vstar;
            std::size_t other = e2_tree ? vstar : wstar;
            std::size_t tri = e2_tree ? f3 : f2;
            std::size_t new_e = tri >> 1;
            std::size_t tree_e = e2_tree ? (f2 >> 1) : (f3 >> 1);

            if (g.origin(cs.traversal_half(g)) == other) {
                assert(cs.cyc_next[other] == attach);
                cs.cyc_next[other] = y;
                cs.cyc_prev[y] = other;
                cs.cyc_next[y] = attach;
                cs.cyc_prev[attach] = y;
                cs.cyc_half[other] = (g.origin(2 * new_e) == other) ? 2 * new_e : 2 * new_e + 1;
                cs.cyc_half[y] = s.tree_half(y, attach);
            } else {
                assert(cs.cyc_next[attach] == other);
                cs.cyc_next[attach] = y;
                cs.cyc_prev[y] = attach;
                cs.cyc_next[y] = other;
                cs.cyc_prev[other] = y;
                cs.cyc_half[attach] = s.tree_half(attach, y);
                cs.cyc_half[y] = (g.origin(2 * new_e) == y) ? 2 * new_e : 2 * new_e + 1;
            }
            cs.on_v[y] = true;
            cs.on_e[cs.nontree] = false;
            cs.on_e[tree_e] = true;
            cs.on_e[new_e] = true;
            if (attach == cs.lca && s.parent[cs.lca] == y)
                cs.lca = y;
            cs.cycle_count += 1;
            assert(cs.inside_cost >= s.cost(y));
            cs.inside_cost -= s.cost(y);
            cs.nontree = new_e;
            cs.inside_half = tri ^ 1;
            continue;
        }

        TreePathToCycle pr = find_path_to_cycle(s, cs, y, markY, markL);
        assert(pr.z < cs.on_v.size() && cs.on_v[pr.z] && "[LT79 tex 628]: path ends on the cycle");
        assert((pr.chain.empty() ? pr.z == y : pr.chain.front() == y) &&
               "[LT79 tex 626]: path starts at the triangle's third vertex");
        std::size_t path_cost = 0;
        for (std::size_t i = 0; i < pr.chain.size(); ++i) {
            std::size_t v = pr.chain[i];
            [[maybe_unused]] std::size_t next = i + 1 < pr.chain.size() ? pr.chain[i + 1] : pr.z;
            assert(!cs.on_v[v] && !chain_mark[v] &&
                   "[LT79 tex 661–663]: path has distinct off-cycle vertices");
            assert((s.parent[v] == next || s.parent[next] == v) &&
                   "[LT79 tex 626–628]: each path step is a tree edge");
            chain_mark[v] = true;
            path_cost += s.cost(v);
        }

        std::size_t T_old = cs.traversal_half(g);
        std::size_t A = g.origin(T_old), B = g.to[T_old];
        assert((A == vstar && B == wstar) || (A == wstar && B == vstar));
        std::size_t half_A_to_y = (A == vstar) ? (f3 ^ 1) : f2;
        std::size_t half_y_to_B = (B == wstar) ? (f2 ^ 1) : f3;
        std::size_t tri_A = (A == vstar) ? f3 : f2;
        std::size_t tri_B = (B == wstar) ? f2 : f3;
        assert(g.origin(half_A_to_y) == A && g.to[half_A_to_y] == y);
        assert(g.origin(half_y_to_B) == y && g.to[half_y_to_B] == B);

        std::vector<std::size_t> chain2(pr.chain.rbegin(), pr.chain.rend());
        CycleSideScan sa, sb;
        sa.s = &s;
        sa.cs = &cs;
        sa.chain_mark = &chain_mark;
        sa.arc_from = pr.z;
        sa.arc_to = A;
        sa.chain = &pr.chain;
        sa.new_half = half_A_to_y;
        sa.new_at_start = true;
        sa.scan_alpha = (tri_A == half_A_to_y);
        sb.s = &s;
        sb.cs = &cs;
        sb.chain_mark = &chain_mark;
        sb.arc_from = B;
        sb.arc_to = pr.z;
        sb.chain = &chain2;
        sb.new_half = half_y_to_B;
        sb.new_at_start = false;
        sb.scan_alpha = (tri_B == half_y_to_B);
        sa.begin();
        sb.begin();

        while (true) {
            if (!sa.step())
                break;
            if (!sb.step())
                break;
        }
        std::size_t inA, inB;
        if (sa.done) {
            inA = sa.cost;
            assert(cs.inside_cost >= inA + path_cost);
            inB = cs.inside_cost - inA - path_cost;
        } else {
            assert(sb.done);
            inB = sb.cost;
            assert(cs.inside_cost >= inB + path_cost);
            inA = cs.inside_cost - inB - path_cost;
        }

        bool pickA = (inA >= inB);
        std::size_t new_e = pickA ? (half_A_to_y >> 1) : (half_y_to_B >> 1);
        std::size_t tri = pickA ? tri_A : tri_B;

        bool lca_dropped = false;
        auto discard = [&](std::size_t first, std::size_t last) {
            std::size_t v = first;
            while (true) {
                std::size_t nx = cs.cyc_next[v];
                cs.on_e[cs.cyc_half[v] >> 1] = false;
                if (v == cs.lca)
                    lca_dropped = true;
                cs.on_v[v] = false;
                cs.cyc_next[v] = cs.cyc_prev[v] = cs.cyc_half[v] = NONE;
                cs.cycle_count -= 1;
                if (v == last)
                    break;
                v = nx;
            }
        };
        cs.on_e[cs.nontree] = false;
        if (pickA) {
            if (B != pr.z) {
                std::size_t before_z = cs.cyc_prev[pr.z];
                discard(B, before_z);
            }
        } else {
            if (A != pr.z) {
                std::size_t after_z = cs.cyc_next[pr.z];
                cs.on_e[cs.cyc_half[pr.z] >> 1] = false;
                discard(after_z, A);
            }
        }

        auto link = [&](std::size_t a, std::size_t b, std::size_t half) {
            cs.cyc_next[a] = b;
            cs.cyc_prev[b] = a;
            cs.cyc_half[a] = half;
            assert(g.origin(half) == a && g.to[half] == b);
            cs.on_e[half >> 1] = true;
        };
        if (pickA) {
            if (pr.chain.empty()) {
                link(A, pr.z, half_A_to_y);
            } else {
                link(A, pr.chain[0], half_A_to_y);
                cs.on_v[pr.chain[0]] = true;
                cs.cycle_count += 1;
                for (std::size_t i = 0; i + 1 < pr.chain.size(); ++i) {
                    link(pr.chain[i], pr.chain[i + 1], s.tree_half(pr.chain[i], pr.chain[i + 1]));
                    cs.on_v[pr.chain[i + 1]] = true;
                    cs.cycle_count += 1;
                }
                link(pr.chain.back(), pr.z, s.tree_half(pr.chain.back(), pr.z));
            }
        } else {
            if (chain2.empty()) {
                link(pr.z, B, half_y_to_B);
            } else {
                std::size_t prev = pr.z;
                for (std::size_t idx = 0; idx < chain2.size(); ++idx) {
                    std::size_t v = chain2[idx];
                    link(prev, v, s.tree_half(prev, v));
                    cs.on_v[v] = true;
                    cs.cycle_count += 1;
                    prev = v;
                }
                link(prev, B, half_y_to_B);
            }
        }

        std::size_t chain_min = NONE;
        for (std::size_t v : pr.chain)
            if (chain_min == NONE || s.depth[v] < s.depth[chain_min])
                chain_min = v;
        if (chain_min != NONE && s.depth[chain_min] < s.depth[pr.z]) {
            cs.lca = chain_min;
        } else if (lca_dropped) {
            cs.lca = pr.z;
        }

        for (std::size_t v : pr.chain)
            chain_mark[v] = false;

        cs.inside_cost = pickA ? inA : inB;
        cs.nontree = new_e;
        cs.inside_half = tri ^ 1;
    }
}

std::vector<std::uint8_t> classify_sides(const ReducedPlanarGraph& s, const SeparatorCycle& cs) {
    const HalfEdgeGraph& g = s.g;
    std::vector<std::uint8_t> side(g.n, 3);

    std::size_t T = cs.traversal_half(g);
    bool inside_is_alpha = (cs.inside_half == (T ^ 1));

    for (std::size_t c = 0; c < g.n; ++c) {
        if (!cs.on_v[c])
            continue;
        side[c] = 2;
        std::size_t out = cs.cyc_half[c];
        std::size_t in = cs.cyc_half[cs.cyc_prev[c]];
        rotation_between(g, out, in ^ 1, [&](std::size_t x) {
            std::size_t w = g.to[x];
            if (cs.on_v[w] || !s.is_tree_edge(x >> 1))
                return;
            side[w] = inside_is_alpha ? 0 : 1;
        });
        rotation_between(g, in ^ 1, out, [&](std::size_t x) {
            std::size_t w = g.to[x];
            if (cs.on_v[w] || !s.is_tree_edge(x >> 1))
                return;
            side[w] = inside_is_alpha ? 1 : 0;
        });
    }

    if (!cs.on_v[s.apex] && side[s.apex] == 3) {
        std::size_t pl = s.parent[cs.lca];
        assert(pl != NONE && side[pl] <= 1 &&
               "[LT79 tex 606–609]: the above-lca component is classified "
               "through the lca's parent edge");
        side[s.apex] = side[pl];
    }

    for (std::size_t v : s.bfs_order) {
        if (side[v] != 3)
            continue;
        assert(s.parent[v] != NONE && side[s.parent[v]] <= 1 &&
               "hanging subtrees inherit their crossing edge's side");
        side[v] = side[s.parent[v]];
    }

#ifndef NDEBUG
    std::size_t check_inside = 0;
    for (std::size_t v = 0; v < g.n; ++v)
        if (side[v] == 0)
            check_inside += s.cost(v);
    assert(check_inside == cs.inside_cost &&
           "[LT79 §3 tex 616–637]: incremental inside-cost bookkeeping "
           "must agree with a direct recount");
#endif
    return side;
}

ReducedPlanarGraph build_reduced_planar_graph(const EmbeddedPlanarGraph& g,
                                              const std::vector<std::size_t>& level, std::size_t l0,
                                              std::size_t l2) {
    ReducedPlanarGraph s;
    auto in_cluster = [&](std::size_t v) { return level[v] != NONE && level[v] <= l0; };
    auto in_middle = [&](std::size_t v) {
        return level[v] != NONE && level[v] > l0 && level[v] < l2;
    };

    std::vector<std::size_t> id(g.num_vertices(), NONE);
    s.g.add_vertex();
    s.orig.push_back(NONE);
    for (std::size_t v = 0; v < g.num_vertices(); ++v) {
        if (!in_middle(v))
            continue;
        id[v] = s.g.add_vertex();
        s.orig.push_back(v);
    }
    s.total = s.g.n - 1;
    assert(s.total >= 1 && "the shrink branch requires a nonempty middle");

    std::vector<bool> table(g.num_vertices(), false);
    std::vector<std::size_t> kept_half(g.num_vertices(), NONE);
    std::vector<std::size_t> apex_order;
    std::vector<bool> half_seen(2 * g.num_edges(), false);

    for (std::size_t v = 0; v < g.num_vertices(); ++v) {
        if (!in_cluster(v))
            continue;
        for (std::size_t h0 : g.incident_halves(v)) {
            if (half_seen[h0])
                continue;
            if (in_cluster(g.half_to(h0)))
                continue;

            std::size_t h = h0;
            do {
                assert(!half_seen[h]);
                half_seen[h] = true;
                std::size_t w = g.half_to(h);
                if (in_middle(w) && !table[w]) {
                    table[w] = true;
                    kept_half[w] = h;
                    apex_order.push_back(w);
                }

                std::size_t j = g.rot_next(h);
                while (in_cluster(g.half_to(j)))
                    j = g.rot_next(j ^ 1);
                h = j;
            } while (h != h0);
        }
    }
    assert(!apex_order.empty() && "the middle is BFS-reachable only through the cluster");

    std::vector<std::size_t> apex_edge(g.num_vertices(), NONE);
    std::vector<std::vector<std::size_t>> rot(s.g.n);
    for (std::size_t w : apex_order) {
        std::size_t half = s.g.raw_edge(0, id[w]);
        apex_edge[w] = half;
        rot[0].push_back(half);
    }
    std::vector<std::size_t> piece_edge_half(g.num_edges(), NONE);
    for (std::size_t v = 0; v < g.num_vertices(); ++v) {
        if (!in_middle(v))
            continue;
        for (std::size_t h : g.incident_halves(v)) {
            std::size_t w = g.half_to(h);
            std::size_t e = h >> 1;
            if (in_middle(w)) {
                if (piece_edge_half[e] == NONE)
                    piece_edge_half[e] = s.g.raw_edge(id[v], id[w]);
                std::size_t mine = (s.g.origin(piece_edge_half[e]) == id[v])
                                       ? piece_edge_half[e]
                                       : (piece_edge_half[e] ^ 1);
                rot[id[v]].push_back(mine);
            } else if (in_cluster(w)) {
                if (kept_half[v] == (h ^ 1))
                    rot[id[v]].push_back(apex_edge[v] ^ 1);
            }
        }
    }
    for (std::size_t v = 0; v < s.g.n; ++v) {
        assert(!rot[v].empty() && "no isolated vertices in the shrunken graph");
        s.g.first_half[v] = rot[v].front();
        for (std::size_t i = 0; i < rot[v].size(); ++i) {
            std::size_t j = (i + 1 == rot[v].size()) ? 0 : i + 1;
            s.g.nxt[rot[v][i]] = rot[v][j];
            s.g.prv[rot[v][j]] = rot[v][i];
        }
    }
    s.g.assert_euler();

    s.parent.assign(s.g.n, NONE);
    s.parent_edge.assign(s.g.n, NONE);
    s.depth.assign(s.g.n, NONE);
    s.desc.assign(s.g.n, 0);
    s.bfs_order.clear();
    s.bfs_order.push_back(0);
    s.depth[0] = 0;
    for (std::size_t qi = 0; qi < s.bfs_order.size(); ++qi) {
        std::size_t v = s.bfs_order[qi];
        s.g.for_out_halves(v, [&](std::size_t h) {
            std::size_t w = s.g.to[h];
            if (s.depth[w] != NONE)
                return;
            s.depth[w] = s.depth[v] + 1;
            s.parent[w] = v;
            s.parent_edge[w] = h >> 1;
            s.bfs_order.push_back(w);
        });
    }
    assert(s.bfs_order.size() == s.g.n && "[LT79 §3 tex 592–595]: the shrunken graph is connected "
                                          "through the apex");

    for (std::size_t v = 0; v < s.g.n; ++v) {
        assert(s.depth[v] <= l2 - l0 - 1 &&
               "[LT79 tex 303–304]: shrunken BFS radius ≤ l2 − l0 − 1");
        (void)v;
    }
    for (std::size_t qi = s.bfs_order.size(); qi-- > 0;) {
        std::size_t v = s.bfs_order[qi];
        s.desc[v] += s.cost(v);
        if (s.parent[v] != NONE)
            s.desc[s.parent[v]] += s.desc[v];
    }

    triangulate(s.g);
    return s;
}

void partition_connected_graph(const EmbeddedPlanarGraph& g, const std::vector<std::size_t>& comp,
                               std::vector<SepPart>& part) {
    const std::size_t total_n = g.num_vertices();
    const std::size_t n_comp = comp.size();

    std::vector<std::size_t> level(g.num_vertices(), NONE);
    std::vector<std::size_t> order;
    order.reserve(n_comp);
    level[comp[0]] = 0;
    order.push_back(comp[0]);
    std::size_t r = 0;
    for (std::size_t qi = 0; qi < order.size(); ++qi) {
        std::size_t v = order[qi];
        r = std::max(r, level[v]);
        for (std::size_t h : g.incident_halves(v)) {
            std::size_t w = g.half_to(h);
            if (level[w] != NONE)
                continue;
            level[w] = level[v] + 1;
            order.push_back(w);
        }
    }
    assert(order.size() == n_comp && "component must be BFS-connected");

    std::vector<std::size_t> L(r + 2, 0);
    for (std::size_t v : order)
        ++L[level[v]];

    std::size_t l1 = NONE, k = 0;
    {
        std::size_t pref = 0;
        for (std::size_t l = 0; l <= r; ++l) {
            pref += L[l];
            if (2 * pref > total_n) {
                l1 = l;
                k = pref;
                break;
            }
        }
    }
    assert(l1 != NONE && "[LT79 tex 331–334]: a component of cost > 2/3 > 1/2 always "
                         "yields l1");

    auto sq_le = [](std::size_t a, std::size_t four_b) { return a * a <= four_b; };
    std::size_t l0 = NONE;
    for (std::size_t l = l1 + 1; l-- > 0;)
        if (sq_le(L[l] + 2 * (l1 - l), 4 * k)) {
            l0 = l;
            break;
        }
    std::size_t l2 = NONE;
    for (std::size_t l = l1 + 1; l <= r + 1; ++l)
        if (sq_le(L[l] + 2 * (l - l1 - 1), 4 * (n_comp - k))) {
            l2 = l;
            break;
        }
    assert(l0 != NONE && l2 != NONE &&
           "[LT79 tex 345–356]: suitable levels l0 and l2 always exist");

    std::size_t below = 0, middle = 0, above = 0;
    for (std::size_t l = 0; l < l0; ++l)
        below += L[l];
    for (std::size_t l = l0 + 1; l < l2; ++l)
        middle += L[l];
    for (std::size_t l = l2 + 1; l <= r; ++l)
        above += L[l];

    auto band_of = [&](std::size_t v) -> int {
        std::size_t l = level[v];
        if (l == l0 || l == l2)
            return 4;
        if (l < l0)
            return 0;
        if (l < l2)
            return 1;
        return 2;
    };

    if (3 * middle <= 2 * total_n) {
        int star_band = (middle >= below && middle >= above) ? 1 : (below >= above ? 0 : 2);
        const std::size_t a_star = star_band == 0 ? below : star_band == 1 ? middle : above;
        const std::size_t b_star = below + middle + above - a_star;
        const bool a_is_pair = b_star > a_star;
        for (std::size_t v : order) {
            int b = band_of(v);
            part[v] =
                (b == 4) ? SepPart::D : ((b == star_band) != a_is_pair ? SepPart::A : SepPart::B);
        }
        return;
    }

    ReducedPlanarGraph s = build_reduced_planar_graph(g, level, l0, l2);

    std::size_t e0 = NONE;
    for (std::size_t e = 0; e < s.g.num_edges(); ++e)
        if (!s.is_tree_edge(e)) {
            e0 = e;
            break;
        }
    assert(e0 != NONE && "triangulated shrunken graph with ≥ 3 vertices has a nontree "
                         "edge");

    SeparatorCycle cs;
    build_initial_cycle(s, e0, cs);
    improve_cycle(s, cs, total_n);

    std::vector<std::uint8_t> side = classify_sides(s, cs);

    std::size_t in_cost = 0, out_cost = 0;
    for (std::size_t v = 1; v < s.g.n; ++v) {
        if (side[v] == 0)
            in_cost += 1;
        else if (side[v] == 1)
            out_cost += 1;
    }
    assert(in_cost == cs.inside_cost);
    assert(3 * in_cost <= 2 * total_n && "[LT79 tex 616–618]: Step 9 exits with inside cost ≤ 2/3");
    assert(3 * out_cost <= 2 * total_n &&
           "[LT79 lem:radius tex 197–200]: the final cycle bounds both "
           "sides by 2/3");
    std::uint8_t a_side = (in_cost >= out_cost) ? 0 : 1;

    std::size_t cyc_real = 0;
    for (std::size_t v = 1; v < s.g.n; ++v)
        if (side[v] == 2)
            ++cyc_real;
    assert(cyc_real + 1 >= cs.cycle_count);
    assert(cyc_real <= 2 * (l2 - l0 - 1) &&
           "[LT79 lem:radius tex 206–208]: cycle length ≤ 2r + 1 with "
           "r = l2 − l0 − 1; the root is on the cycle (≤ 2r non-root "
           "vertices) or off it (the paths meet below the root, ≤ 2r − 1 "
           "total), so the non-apex count is at most 2r either way");

    for (std::size_t v : order) {
        int b = band_of(v);
        part[v] = (b == 4) ? SepPart::D : SepPart::B;
    }

    for (std::size_t sv = 1; sv < s.g.n; ++sv) {
        std::size_t v = s.orig[sv];
        if (side[sv] == 2)
            part[v] = SepPart::D;
        else
            part[v] = (side[sv] == a_side) ? SepPart::A : SepPart::B;
    }
}

}

std::vector<SepPart> planar_separator(const EmbeddedPlanarGraph& g) {
    const std::size_t n = g.num_vertices();
    assert(n >= 1);
    std::vector<SepPart> part(n, SepPart::B);

    std::vector<std::size_t> comp_id(n, NONE);
    std::vector<std::vector<std::size_t>> comps;
    for (std::size_t v0 = 0; v0 < n; ++v0) {
        if (comp_id[v0] != NONE)
            continue;
        comps.emplace_back();
        auto& c = comps.back();
        comp_id[v0] = comps.size() - 1;
        c.push_back(v0);
        for (std::size_t qi = 0; qi < c.size(); ++qi)
            for (std::size_t h : g.incident_halves(c[qi])) {
                std::size_t w = g.half_to(h);
                if (comp_id[w] != NONE)
                    continue;
                comp_id[w] = comps.size() - 1;
                c.push_back(w);
            }
    }

    std::size_t big = 0;
    for (std::size_t i = 1; i < comps.size(); ++i)
        if (comps[i].size() > comps[big].size())
            big = i;

    if (3 * comps[big].size() > 2 * n) {
        partition_connected_graph(g, comps[big], part);
    } else if (3 * comps[big].size() > n) {
        for (std::size_t v : comps[big])
            part[v] = SepPart::A;
    } else {
        std::size_t acc = 0;
        for (const auto& c : comps) {
            if (3 * acc > n)
                break;
            for (std::size_t v : c)
                part[v] = SepPart::A;
            acc += c.size();
        }
    }

#ifndef NDEBUG
    std::size_t na = 0, nb = 0, nd = 0;
    for (std::size_t v = 0; v < n; ++v) {
        if (part[v] == SepPart::A)
            ++na;
        else if (part[v] == SepPart::B)
            ++nb;
        else
            ++nd;
    }
    assert(3 * na <= 2 * n && "[C91 §3.4 tex 300] (ii): |A| ≤ 2μ/3");
    assert(3 * nb <= 2 * n && "[C91 §3.4 tex 300] (ii): |B| ≤ 2μ/3");
    assert(nd * nd <= 8 * n && "[C91 §3.4 tex 301] (iii): |D| ≤ √(8μ)");
    for (std::size_t e = 0; e < g.num_edges(); ++e) {
        SepPart pu = part[g.edge_u(e)], pv = part[g.edge_v(e)];
        assert(
            !((pu == SepPart::A && pv == SepPart::B) || (pu == SepPart::B && pv == SepPart::A)) &&
            "[C91 §3.4 tex 299] (i): no edge joins A and B");
    }
#endif
    return part;
}

SeparatorDecomposition build_separator_decomposition(const EmbeddedPlanarGraph& g,
                                                     std::size_t trace_structure) {
    const std::size_t mu = g.num_vertices();
    SeparatorDecomposition out;
    out.subset.assign(mu, NONE);

    auto is_leaf = [&](std::size_t sz) {
        __extension__ typedef unsigned __int128 u128;
        const u128 s = sz, m = mu;
        return s * s * s <= m * m;
    };

    struct Item {
        EmbeddedPlanarGraph graph;
        std::vector<std::size_t> ids;
    };
    std::vector<Item> stack;
    {
        std::vector<std::size_t> all(mu);
        for (std::size_t v = 0; v < mu; ++v)
            all[v] = v;
        std::vector<std::size_t> ids;
        stack.push_back(Item{g.induced(all, &ids), std::move(ids)});
    }

    while (!stack.empty()) {
        Item it = std::move(stack.back());
        stack.pop_back();
        std::size_t sz = it.graph.num_vertices();
        if (sz == 0)
            continue;
        if (is_leaf(sz)) {
            if (auto* trace = AnimationTrace::current(); trace && trace_structure != NONE)
                trace->indices("separator_leaf",
                               {{"structure", trace_structure}, {"subset", out.num_subsets}},
                               "faces", it.ids);
            for (std::size_t v = 0; v < sz; ++v)
                out.subset[it.ids[v]] = out.num_subsets;
            ++out.num_subsets;
            continue;
        }
        std::vector<SepPart> part = planar_separator(it.graph);
        if (auto* trace = AnimationTrace::current(); trace && trace_structure != NONE) {
            std::vector<std::size_t> parts;
            parts.reserve(part.size());
            for (auto value : part)
                parts.push_back(static_cast<std::size_t>(value));
            trace->indices("separator_partition", {{"structure", trace_structure}}, "faces", it.ids,
                           "parts", parts);
        }
        std::vector<std::size_t> a_nodes, b_nodes;
        for (std::size_t v = 0; v < sz; ++v) {
            if (part[v] == SepPart::D) {
                out.subset[it.ids[v]] = NONE;
                ++out.dstar_size;
            } else if (part[v] == SepPart::A) {
                a_nodes.push_back(v);
            } else {
                b_nodes.push_back(v);
            }
        }
        for (auto* nodes : {&a_nodes, &b_nodes}) {
            if (nodes->empty())
                continue;
            std::vector<std::size_t> sub_ids;
            EmbeddedPlanarGraph sub = it.graph.induced(*nodes, &sub_ids);
            std::vector<std::size_t> top_ids(nodes->size());
            for (std::size_t i = 0; i < nodes->size(); ++i)
                top_ids[i] = it.ids[(*nodes)[i]];
            stack.push_back(Item{std::move(sub), std::move(top_ids)});
        }
    }

#ifndef NDEBUG
    {
        const Exact count{out.dstar_size}, size{mu};
        assert(count * count * count <= 19 * 19 * 19 * size * size &&
               "[C91 §3.4 tex 304]: |D*| = O(μ^{2/3}) with the derived "
               "constant");
    }

    {
        std::vector<std::size_t> sizes(out.num_subsets, 0);
        for (std::size_t v = 0; v < mu; ++v)
            if (out.subset[v] != NONE)
                ++sizes[out.subset[v]];
        for (std::size_t szv : sizes) {
            __extension__ typedef unsigned __int128 u128;
            assert((u128)szv * szv * szv <= (u128)mu * mu &&
                   "[C91 §3.4 tex 304]: each |D_i| ≤ μ^{2/3}");
        }
    }
#endif
    return out;
}

}
