// src/merge/planar_separator.cpp
//
// [LT79 §3 tex 520–668]: the Partitioning Algorithm (Steps 1–10), plus
// the recursive decomposition of [C91 §3.4 tex 304].  All step numbers
// and inequalities cite papers/lipton-tarjan1979-transcribed.tex.
//
// Costs are the unit costs of [LT79 Corollary sqrt-sep tex 382–393]
// (1/n per vertex, apex of the shrunken graph 0 per [LT79 tex 301]), so
// every cost comparison is an exact integer test on vertex counts:
// cost(X) ≤ 2/3 ⟺ 3·|X| ≤ 2·n.

#include "planar_separator.h"

#include <algorithm>
#include <cmath>

namespace chazelle {

// ════════════════════════════════════════════════════════════════
//  EmbeddedPlanarGraph
// ════════════════════════════════════════════════════════════════

EmbeddedPlanarGraph::EmbeddedPlanarGraph(
        std::size_t num_vertices,
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
        // [LT79 §3 tex 526–527]: the cyclic clockwise/counter-clockwise
        // rotation links around each endpoint.
        const auto& r = rot_[v];
        for (std::size_t i = 0; i < r.size(); ++i) {
            std::size_t j = (i + 1 == r.size()) ? 0 : i + 1;
            nxt_[r[i]] = r[j];
            prv_[r[j]] = r[i];
        }
    }

    for (std::size_t e = 0; e < edges_.size(); ++e) {
        assert(seen[e] == 2 &&
               "[LT79 §3 tex 524–528]: each edge appears exactly once in "
               "each endpoint's rotation");
        (void)e;
    }
    (void)seen;
}

EmbeddedPlanarGraph EmbeddedPlanarGraph::induced(
        const std::vector<std::size_t>& nodes,
        std::vector<std::size_t>* old_index) const {
    std::vector<std::size_t> new_index(num_vertices_, NONE);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        assert(nodes[i] < num_vertices_);
        assert(new_index[nodes[i]] == NONE && "duplicate node in induced()");
        new_index[nodes[i]] = i;
    }

    // [C91 §3.4 tex 304]: keep "the edges of G that join only nodes of
    // A".  Rotation order preserved — an induced subgraph of an embedded
    // planar graph inherits its embedding.
    std::vector<std::pair<std::size_t, std::size_t>> new_edges;
    std::vector<std::size_t> edge_map(edges_.size(), NONE);
    std::vector<std::vector<std::size_t>> new_rot(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        std::size_t v = nodes[i];
        for (std::size_t h : rot_[v]) {
            std::size_t w = half_to(h);
            if (new_index[w] == NONE) continue;
            std::size_t e = h >> 1;
            if (edge_map[e] == NONE) {
                edge_map[e] = new_edges.size();
                new_edges.emplace_back(new_index[edges_[e].first],
                                       new_index[edges_[e].second]);
            }
            new_rot[i].push_back(edge_map[e]);
        }
    }

    if (old_index) *old_index = nodes;
    return EmbeddedPlanarGraph(nodes.size(), std::move(new_edges), new_rot);
}

// ════════════════════════════════════════════════════════════════
//  Mutable half-edge graph for [LT79] Steps 6–9
// ════════════════════════════════════════════════════════════════

namespace {

struct HalfEdgeGraph {
    std::vector<std::size_t> to;         // per half: destination
    std::vector<std::size_t> nxt, prv;   // cw rotation links around origin
    std::vector<std::size_t> first_half; // per vertex: any outgoing half
    std::size_t n = 0;

    std::size_t num_edges() const { return to.size() / 2; }
    std::size_t origin(std::size_t h) const { return to[h ^ 1]; }

    std::size_t add_vertex() {
        first_half.push_back(NONE);
        return n++;
    }

    std::size_t raw_edge(std::size_t u, std::size_t v) {
        assert(u != v && "[LT79]: loops never arise");
        (void)u;
        to.push_back(v); to.push_back(u);
        nxt.push_back(NONE); nxt.push_back(NONE);
        prv.push_back(NONE); prv.push_back(NONE);
        return to.size() - 2;            // the u→v half
    }

    void insert_before(std::size_t h, std::size_t pos) {
        std::size_t p = prv[pos];
        nxt[p] = h; prv[h] = p;
        nxt[h] = pos; prv[pos] = h;
    }
    void insert_after(std::size_t h, std::size_t pos) {
        std::size_t nx = nxt[pos];
        nxt[pos] = h; prv[h] = pos;
        nxt[h] = nx; prv[nx] = h;
    }

    template <typename F>
    void for_out_halves(std::size_t v, F&& f) const {
        std::size_t h0 = first_half[v];
        if (h0 == NONE) return;
        std::size_t h = h0;
        do { f(h); h = nxt[h]; } while (h != h0);
    }

    // The next half of the face lying on a fixed side of h.
    std::size_t face_next(std::size_t h) const { return nxt[h ^ 1]; }

    std::size_t count_faces() const {
        std::vector<bool> vis(to.size(), false);
        std::size_t faces = 0;
        for (std::size_t h = 0; h < to.size(); ++h) {
            if (vis[h]) continue;
            ++faces;
            std::size_t x = h;
            do { vis[x] = true; x = face_next(x); } while (x != h);
        }
        return faces;
    }

    // V − E + F = 2 ⟺ the rotation system is a sphere embedding — the
    // structural content of [LT79 §3 tex 534–538] Step 1.
    void assert_euler() const {
#ifndef NDEBUG
        assert(n >= 1 && !to.empty());
        assert(n + count_faces() == num_edges() + 2 &&
               "[LT79 §3 tex 534–538]: rotation system must be a planar "
               "(sphere) embedding");
#endif
    }
};

// The shrunken middle graph of [LT79 §3 tex 575–586] Step 6, with the
// Step 7 (tex 592–597) spanning tree and descendant costs.
struct Shrunken {
    HalfEdgeGraph g;
    std::size_t apex = 0;                 // vertex 0; cost 0 [LT79 tex 301]
    std::vector<std::size_t> orig;        // shrunken vertex → piece vertex
    std::vector<std::size_t> parent;      // BFS parent (NONE at apex)
    std::vector<std::size_t> parent_edge; // BFS parent edge id
    std::vector<std::size_t> depth;
    std::vector<std::size_t> desc;        // descendant cost incl. self
    std::vector<std::size_t> bfs_order;
    std::size_t total = 0;                // total cost = # middle vertices

    std::size_t cost(std::size_t v) const { return v == apex ? 0 : 1; }
    bool is_tree_edge(std::size_t e) const {
        std::size_t u = g.origin(2 * e), w = g.to[2 * e];
        return parent_edge[u] == e || parent_edge[w] == e;
    }
    // Half of tree edge (a, b) directed a→b.
    std::size_t tree_half(std::size_t a, std::size_t b) const {
        std::size_t te = (parent[a] == b) ? parent_edge[a] : parent_edge[b];
        assert(te != NONE && (parent[a] == b || parent[b] == a) &&
               "tree_half requires a tree edge");
        return (g.origin(2 * te) == a) ? 2 * te : 2 * te + 1;
    }
};

// ── Current fundamental cycle for [LT79 §3 tex 603–641] Steps 8–9 ──
//
// The cycle is a doubly-linked traversal (cyc_next / cyc_half).  At a
// cycle vertex c with traversal halves o(c) = cyc_half[c] (outgoing)
// and i(c) = cyc_half[cyc_prev[c]] (incoming), the rotation splits into
//   side-α(c): halves strictly between o(c) and twin(i(c)) cw,
//   side-β(c): halves strictly between twin(i(c)) and o(c) cw;
// each is one geometric side of the Jordan curve, consistently along
// the whole traversal.  The face of half h occupies the rotation corner
// (rot-prev(h), h) at origin(h); hence F(twin(T)) lies on side-α and
// F(T) on side-β, where T is the traversal half of the nontree edge.
struct CycleState {
    std::vector<bool> on_v;
    std::vector<bool> on_e;
    std::vector<std::size_t> cyc_next, cyc_prev;
    std::vector<std::size_t> cyc_half;   // outgoing traversal half
    std::size_t nontree = NONE;
    std::size_t inside_half = NONE;      // half of `nontree` whose face
                                         // (face_next trace) lies inside
    std::size_t inside_cost = 0;         // strict-inside cost
    std::size_t cycle_count = 0;         // # cycle vertices
    std::size_t lca = NONE;              // min-depth cycle vertex

    // Traversal half of the nontree edge.
    std::size_t traversal_half(const HalfEdgeGraph& g) const {
        std::size_t a = g.origin(2 * nontree);
        if (cyc_half[a] != NONE && (cyc_half[a] >> 1) == nontree)
            return cyc_half[a];
        std::size_t b = g.to[2 * nontree];
        assert(cyc_half[b] != NONE && (cyc_half[b] >> 1) == nontree);
        return cyc_half[b];
    }
};

// Visit the halves strictly between `from_excl` and `to_excl` in cw
// (nxt) order at their common origin.
template <typename F>
void rotation_between(const HalfEdgeGraph& g, std::size_t from_excl,
                      std::size_t to_excl, F&& f) {
    if (from_excl == to_excl) return;
    for (std::size_t x = g.nxt[from_excl]; x != to_excl; x = g.nxt[x]) {
        assert(x != from_excl && "rotation walk must terminate");
        f(x);
    }
}

// [LT79 §3 tex 606–609]: cost associated with a tree edge crossing the
// boundary from vertex c to hanging vertex w — "the descendant cost of
// w if c is the parent of w, and the cost of all vertices minus the
// descendant cost of c if w is the parent of c."
std::size_t hanging_cost(const Shrunken& s, std::size_t c, std::size_t w) {
    if (s.parent[w] == c) return s.desc[w];
    assert(s.parent[c] == w &&
           "[LT79 §3 tex 606–609]: a crossing tree edge is a parent edge");
    return s.total - s.desc[c];
}

// ════════════════════════════════════════════════════════════════
//  Step 7 (part): triangulation — [LT79 §3 tex 595–597]
// ════════════════════════════════════════════════════════════════
//
// "Make all faces of the new graph into triangles by scanning the
// boundary of each face and adding (nontree) edges as necessary."
// The paper leaves the details open (tex 652–653: "We urge readers to
// fill in the details"); we cut ears: any face-walk position with
// w_{i−1} ≠ w_{i+1} admits the chord (w_{i−1}, w_{i+1}).  Chords may
// duplicate existing edges (parallel edges are harmless — a face of a
// loopless multigraph with 3 sides has 3 distinct vertices, which is
// all Lemma lem:radius's case analysis [LT79 tex 222–267] needs).  A
// face whose walk alternates between two vertices admits no loop-free
// chord and is left as a residual parallel bundle; Step 9 re-anchors
// the cycle across such faces (they enclose no vertices, so this
// preserves both the cost bookkeeping and the [LT79 tex 656–667]
// termination argument).
void triangulate(HalfEdgeGraph& g) {
    std::vector<bool> vis(g.to.size(), false);
    for (std::size_t h0 = 0; h0 < g.to.size(); ++h0) {
        if (vis[h0]) continue;
        std::vector<std::size_t> walk;
        std::size_t x = h0;
        do {
            vis[x] = true;
            walk.push_back(x);
            x = g.face_next(x);
        } while (x != h0);

        // Doubly-linked positions over `walk` + a cursor that backs up
        // one slot after each cut: amortized O(|walk|).
        std::size_t k = walk.size();
        if (k <= 3) continue;
        std::vector<std::size_t> nx(k), pv(k);
        for (std::size_t i = 0; i < k; ++i) {
            nx[i] = (i + 1) % k;
            pv[i] = (i + k - 1) % k;
        }
        std::size_t len = k;
        std::size_t cur = 0;
        std::size_t since_progress = 0;
        while (len > 3 && since_progress <= len) {
            std::size_t i1 = cur, i2 = nx[cur];
            std::size_t h1 = walk[i1], h2 = walk[i2];
            std::size_t a = g.origin(h1), c = g.to[h2];
            if (a == c) {                     // loop chord — skip
                cur = nx[cur];
                ++since_progress;
                continue;
            }
            // Chord a→c cutting off triangle (h1, h2, twin(p)): p goes
            // before h1 in rot(a), twin(p) after twin(h2) in rot(c);
            // then face_next(h2) = twin(p), face_next(twin(p)) = h1 and
            // the reduced face continues at the old face_next(h2).
            std::size_t p = g.raw_edge(a, c);
            g.insert_before(p, h1);
            g.insert_after(p ^ 1, h2 ^ 1);
            assert(g.face_next(h2) == (p ^ 1) && g.face_next(p ^ 1) == h1 &&
                   "ear cut must close the triangle (h1, h2, twin(p))");
            vis.push_back(true); vis.push_back(true);
            // Replace h1, h2 by p in the walk.
            walk[i1] = p;
            nx[i1] = nx[i2]; pv[nx[i2]] = i1;
            --len;
            cur = pv[i1];
            since_progress = 0;
        }
        // len > 3 here ⟹ fully alternating residual bundle face
        // (every candidate chord was a loop) — see note above.
    }
    g.assert_euler();
}

// ════════════════════════════════════════════════════════════════
//  Step 8 — [LT79 §3 tex 603–614]
// ════════════════════════════════════════════════════════════════

void build_initial_cycle(const Shrunken& s, std::size_t e, CycleState& cs) {
    const HalfEdgeGraph& g = s.g;
    std::size_t v1 = g.origin(2 * e), w1 = g.to[2 * e];

    cs.on_v.assign(g.n, false);
    cs.on_e.assign(g.num_edges(), false);
    cs.cyc_next.assign(g.n, NONE);
    cs.cyc_prev.assign(g.n, NONE);
    cs.cyc_half.assign(g.n, NONE);

    // [LT79 §3 tex 603–605]: "Locate the corresponding cycle by
    // following parent pointers from v1 and w1."
    std::vector<std::size_t> pv, pw;
    for (std::size_t x = v1; x != NONE; x = s.parent[x]) pv.push_back(x);
    for (std::size_t x = w1; x != NONE; x = s.parent[x]) pw.push_back(x);
    while (pv.size() >= 2 && pw.size() >= 2 &&
           pv[pv.size() - 2] == pw[pw.size() - 2]) {
        pv.pop_back();
        pw.pop_back();
    }
    assert(pv.back() == pw.back() && "tree paths must meet at the LCA");
    cs.lca = pv.back();

    // Traversal: v1 →(nontree)→ w1 →tree→ lca →tree→ … → v1.  The lca
    // appears exactly once: it is skipped from pw/pv when it coincides
    // with an endpoint of the nontree edge (one chain empty).
    std::vector<std::size_t> cyc;
    cyc.push_back(v1);
    for (std::size_t i = 0; i + 1 < pw.size(); ++i) cyc.push_back(pw[i]);
    if (cs.lca != v1) cyc.push_back(cs.lca);
    for (std::size_t i = pv.size() - 1; i-- > 1; ) cyc.push_back(pv[i]);

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

    // [LT79 §3 tex 605–610]: "Compute the cost on each side of this
    // cycle by scanning the tree edges incident on either side."
    std::size_t cost_alpha = 0, cost_beta = 0, cycle_cost = 0;
    for (std::size_t i = 0; i < cyc.size(); ++i) {
        std::size_t c = cyc[i];
        cycle_cost += s.cost(c);
        std::size_t out = cs.cyc_half[c];
        std::size_t in = cs.cyc_half[cyc[(i + cyc.size() - 1) % cyc.size()]];
        assert(g.to[in] == c);
        auto tally = [&](std::size_t x, std::size_t& acc) {
            std::size_t w = g.to[x];
            if (cs.on_v[w]) return;
            if (!s.is_tree_edge(x >> 1)) return;
            acc += hanging_cost(s, c, w);
        };
        rotation_between(g, out, in ^ 1,
                         [&](std::size_t x) { tally(x, cost_alpha); });
        rotation_between(g, in ^ 1, out,
                         [&](std::size_t x) { tally(x, cost_beta); });
    }
    // Every non-cycle vertex hangs under exactly one crossing tree edge.
    assert(cost_alpha + cost_beta + cycle_cost == s.total &&
           "[LT79 §3 tex 605–610]: the two sides plus the cycle must "
           "account for the whole graph");

    // "Determine which side of the cycle has greater cost and call it
    // the 'inside'" (tex 609–610).  F(twin(T)) lies on side-α for the
    // traversal half T = 2e of the nontree edge (corner rule).
    if (cost_alpha >= cost_beta) {
        cs.inside_cost = cost_alpha;
        cs.inside_half = 2 * e + 1;
    } else {
        cs.inside_cost = cost_beta;
        cs.inside_half = 2 * e;
    }
}

// ════════════════════════════════════════════════════════════════
//  Step 9 — [LT79 §3 tex 616–641]
// ════════════════════════════════════════════════════════════════

// Lazy boundary scan of one case-3c candidate cycle, emitting hanging
// tree-edge costs on the candidate's inside side.  Used for [LT79 §3
// tex 629–632]: "Scan the tree edges inside the (y, w_i) cycle,
// alternately scanning an edge in one cycle and an edge in the other
// cycle.  Stop scanning when all edges inside one of the cycles have
// been scanned."
//
// Candidate boundary: [arc_from ..old cyc_next.. arc_to], then either
// the new nontree edge into chain.front() (`new_at_start`, the
// candidate through A = origin of the old traversal) or a tree edge
// (the candidate through B), then the chain, then back to arc_from
// (tree edge resp. new nontree edge).  The candidate's inside face is
// F(twin(tri_half)); by the corner rule the inside lies on side-α of
// the candidate traversal iff tri_half equals the candidate's
// traversal half (`scan_alpha`), else on side-β.
struct SideScan {
    const Shrunken* s = nullptr;
    const CycleState* cs = nullptr;
    const std::vector<bool>* chain_mark = nullptr;
    std::size_t arc_from = NONE, arc_to = NONE;
    const std::vector<std::size_t>* chain = nullptr;
    std::size_t new_half = NONE;   // traversal half of the new nontree edge
    bool new_at_start = false;
    // Which rotation side of the candidate traversal is its inside —
    // side-α iff the candidate's triangle half equals its traversal
    // half (corner rule).
    bool scan_alpha = false;

    bool done = false;
    std::size_t cost = 0;

    std::size_t cur = NONE;
    int seg = 0;                   // 0 = arc, 1 = chain
    std::size_t cpos = 0;
    std::size_t x = NONE, xstop = NONE;

    std::size_t out_half() const {
        if (seg == 0) {
            if (cur != arc_to) return cs->cyc_half[cur];
            if (chain->empty()) return new_half;
            return new_at_start ? new_half
                                : s->tree_half(cur, chain->front());
        }
        if (cpos + 1 < chain->size())
            return s->tree_half(cur, (*chain)[cpos + 1]);
        return new_at_start ? s->tree_half(cur, arc_from) : new_half;
    }
    std::size_t in_half() const {
        if (seg == 0) {
            if (cur == arc_from) {
                if (chain->empty()) return new_half;
                return new_at_start ? s->tree_half(chain->back(), arc_from)
                                    : new_half;
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
            // side-α: strictly between o and twin(i).
            x = s->g.nxt[o];
            xstop = i ^ 1;
        } else {
            // side-β: strictly between twin(i) and o.
            x = s->g.nxt[i ^ 1];
            xstop = o;
        }
        if (x == xstop) x = NONE;
    }

    void begin() {
        cur = arc_from;
        seg = 0;
        cpos = 0;
        setup_vertex();
    }

    // Advance to the next boundary vertex; false when the walk closes.
    bool advance_vertex() {
        if (seg == 0) {
            if (cur == arc_to) {
                if (chain->empty()) return false;
                seg = 1;
                cpos = 0;
                cur = chain->front();
            } else {
                cur = cs->cyc_next[cur];
            }
        } else {
            if (cpos + 1 >= chain->size()) return false;
            ++cpos;
            cur = (*chain)[cpos];
        }
        setup_vertex();
        return true;
    }

    // Process one scanned half; false when the whole walk is complete.
    bool step() {
        if (done) return false;
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
        if (x == xstop) x = NONE;
        return true;
    }
};

// Tree path from y down/around to the current cycle.  Normally the
// parent ascent from y reaches the cycle directly at z ([LT79 §3 tex
// 626–628]: "determine the tree path from y to the (v_i, w_i) cycle by
// following parent pointers from y; let z be the vertex on the cycle
// reached during this search").  When the BFS root lies strictly
// inside the current cycle, the ascent can reach the root without
// touching the cycle; the fundamental-cycle path then detours through
// l' = lca(y, L) where L is the cycle's minimum-depth vertex, giving
// z = L and chain = y…l'…(child of L).  We find the meeting point by
// ascending from y and from L alternately (the same alternation device
// as tex 629–632), so the work stays proportional to the chain length
// and the [LT79 tex 656–667] O(n) accounting is preserved (every chain
// vertex joins the next boundary).
struct PathResult {
    std::size_t z = NONE;
    std::vector<std::size_t> chain;   // y first; excludes z; may be empty
};

PathResult find_path_to_cycle(const Shrunken& s, const CycleState& cs,
                              std::size_t y, std::vector<bool>& markY,
                              std::vector<bool>& markL) {
    PathResult r;
    if (cs.on_v[y]) {                 // y already on the cycle: z = y
        r.z = y;
        return r;
    }
    std::vector<std::size_t> listY{y}, listL{cs.lca};
    markY[y] = true;
    assert(cs.on_v[cs.lca] && "cycle lca must be maintained");

    bool y_alive = true, l_alive = true;
    std::size_t meetY = NONE, meetL = NONE;  // (index into listY, vertex)
    while (true) {
        if (y_alive) {
            std::size_t p = s.parent[listY.back()];
            if (p == NONE) {
                y_alive = false;
            } else if (cs.on_v[p]) {
                r.z = p;                       // normal case
                r.chain = std::move(listY);
                for (std::size_t v : r.chain) markY[v] = false;
                for (std::size_t v : listL) markL[v] = false;
                return r;
            } else if (markL[p]) {
                meetY = NONE; meetL = p;       // met on L's path
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
                meetY = q; meetL = NONE;       // met on y's path
                break;
            } else {
                assert(!cs.on_v[q] &&
                       "the cycle lca's ancestors lie off the cycle");
                listL.push_back(q);
                markL[q] = true;
            }
        }
        assert((y_alive || l_alive) &&
               "both ascents reach the root, so they must meet");
    }

    // Root-inside detour: z = L; chain = y…l' then l'…(child of L).
    r.z = cs.lca;
    std::size_t lprime = (meetL != NONE) ? meetL : meetY;
    for (std::size_t v : listY) {
        r.chain.push_back(v);
        if (v == lprime) break;
    }
    if (r.chain.empty() || r.chain.back() != lprime) {
        // l' lies on L's list; append the ascent part of y fully, then
        // descend: listY entirely above? No — l' on listL means the y
        // ascent stopped when its next parent was marked L-side.
        assert(meetL != NONE);
        // chain currently = all of listY (none equal to lprime).
        r.chain.assign(listY.begin(), listY.end());
        r.chain.push_back(lprime);
    }
    // Descend from l' to L's child: the suffix of listL below l'.
    std::size_t cut = listL.size();
    for (std::size_t i = 0; i < listL.size(); ++i)
        if (listL[i] == lprime) { cut = i; break; }
    // Elements listL[1..cut) are strictly between L and l'.
    for (std::size_t i = cut; i-- > 1; )
        r.chain.push_back(listL[i]);

    for (std::size_t v : listY) markY[v] = false;
    for (std::size_t v : listL) markL[v] = false;
    return r;
}

// One Step-9 improvement pass; returns when 3·inside ≤ 2·total_units.
void improve_cycle(const Shrunken& s, CycleState& cs,
                   std::size_t total_units) {
    const HalfEdgeGraph& g = s.g;
    std::vector<bool> markY(g.n, false), markL(g.n, false);
    std::vector<bool> chain_mark(g.n, false);

    // Termination: each iteration removes ≥ 1 face from the inside
    // ([LT79 §3 tex 656–658]).
    std::size_t guard = 2 * g.count_faces() + 8;

    while (3 * cs.inside_cost > 2 * total_units) {
        assert(guard-- > 0 &&
               "[LT79 §3 tex 656–658]: Step 9 must terminate within the "
               "face budget");

        std::size_t f1 = cs.inside_half;
        std::size_t f2 = g.face_next(f1);
        std::size_t f3 = g.face_next(f2);

        if (f3 == f1 || (g.face_next(f3) != f1)) {
            // Inside face is a 2-gon or a residual parallel bundle
            // (vertex set {v, w}; see triangulate()).  Re-anchor the
            // cycle on the first non-tree parallel edge across the
            // face: same tree path, same inside cost, one less face
            // inside.
            std::size_t pick = NONE;
            for (std::size_t h = f2; h != f1; h = g.face_next(h)) {
                std::size_t e = h >> 1;
                if (e == cs.nontree) continue;
                if (s.is_tree_edge(e)) continue;
                pick = h;
                break;
            }
            assert(pick != NONE &&
                   "a bundle face inside a positive-cost cycle carries a "
                   "non-tree parallel edge");
            std::size_t c0 = g.origin(cs.traversal_half(g));
            cs.on_e[cs.nontree] = false;
            cs.nontree = pick >> 1;
            cs.on_e[cs.nontree] = true;
            cs.cyc_half[c0] =
                (g.origin(2 * cs.nontree) == c0) ? 2 * cs.nontree
                                                 : 2 * cs.nontree + 1;
            cs.inside_half = pick ^ 1;
            continue;
        }

        // Triangle (v*, w*, y) inside, adjacent to the nontree edge
        // ([LT79 §3 tex 620–621]).
        std::size_t vstar = g.origin(f1), wstar = g.to[f1];
        std::size_t y = g.to[f2];
        assert((f1 >> 1) == cs.nontree);
        assert(y != vstar && y != wstar &&
               "triangles of a loopless multigraph have distinct corners");

        bool e2_on = cs.on_e[f2 >> 1];
        bool e3_on = cs.on_e[f3 >> 1];
        assert(!(e2_on && e3_on) &&
               "[LT79 tex 228–229] case 1: the face cannot BE the cycle "
               "while vertices remain inside");

        if (e2_on || e3_on) {
            // ── Case 2 [LT79 tex 231–233 / 620–624] ──  Drop the
            // nontree endpoint between the two cycle edges; the inside
            // is unchanged and one face leaves it.
            std::size_t drop = e2_on ? wstar : vstar;
            std::size_t keep = e2_on ? vstar : wstar;
            std::size_t tri = e2_on ? f3 : f2;    // half of the new edge
            std::size_t new_e = tri >> 1;
            std::size_t dropped_cycle_e = e2_on ? (f2 >> 1) : (f3 >> 1);

            // The old traversal direction across the nontree edge (its
            // stored half's origin) fixes where the relink goes; the
            // cyc_next test alone is ambiguous on a 2-vertex cycle.
            if (g.origin(cs.traversal_half(g)) == keep) {
                // … → keep →(nontree)→ drop → y → …  becomes keep → y.
                assert(cs.cyc_next[keep] == drop && cs.cyc_next[drop] == y);
                cs.cyc_next[keep] = y;
                cs.cyc_prev[y] = keep;
                cs.cyc_half[keep] =
                    (g.origin(2 * new_e) == keep) ? 2 * new_e : 2 * new_e + 1;
            } else {
                // … → y → drop →(nontree)→ keep → …  becomes y → keep.
                assert(cs.cyc_next[y] == drop && cs.cyc_next[drop] == keep);
                cs.cyc_next[y] = keep;
                cs.cyc_prev[keep] = y;
                cs.cyc_half[y] =
                    (g.origin(2 * new_e) == y) ? 2 * new_e : 2 * new_e + 1;
            }
            cs.on_v[drop] = false;
            cs.cyc_next[drop] = cs.cyc_prev[drop] = cs.cyc_half[drop] = NONE;
            cs.on_e[cs.nontree] = false;
            cs.on_e[dropped_cycle_e] = false;
            cs.on_e[new_e] = true;
            if (drop == cs.lca) cs.lca = y;
            cs.cycle_count -= 1;
            cs.nontree = new_e;
            cs.inside_half = tri ^ 1;
            continue;
        }

        bool e2_tree = s.is_tree_edge(f2 >> 1);
        bool e3_tree = s.is_tree_edge(f3 >> 1);
        assert(!(e2_tree && e3_tree) &&
               "[LT79 tex 237–238] case 3a: two tree edges would close a "
               "tree cycle");

        if (e2_tree || e3_tree) {
            // ── Case 3b [LT79 tex 240–251 / 620–624] ──  Insert y on
            // the cycle; inside loses y and one face.
            assert(!cs.on_v[y] &&
                   "a tree edge from the cycle to an on-cycle vertex "
                   "would close a tree cycle");
            std::size_t attach = e2_tree ? wstar : vstar;
            std::size_t other = e2_tree ? vstar : wstar;
            std::size_t tri = e2_tree ? f3 : f2;   // half of new nontree
            std::size_t new_e = tri >> 1;
            std::size_t tree_e = e2_tree ? (f2 >> 1) : (f3 >> 1);

            // As in case 2, the old traversal direction across the
            // nontree edge fixes the insertion side (cyc_next alone is
            // ambiguous on a 2-vertex cycle).
            if (g.origin(cs.traversal_half(g)) == other) {
                // … → other →(nontree)→ attach → …  becomes
                // other →(new_e)→ y →(tree_e)→ attach.
                assert(cs.cyc_next[other] == attach);
                cs.cyc_next[other] = y;
                cs.cyc_prev[y] = other;
                cs.cyc_next[y] = attach;
                cs.cyc_prev[attach] = y;
                cs.cyc_half[other] =
                    (g.origin(2 * new_e) == other) ? 2 * new_e
                                                   : 2 * new_e + 1;
                cs.cyc_half[y] = s.tree_half(y, attach);
            } else {
                // … → attach →(nontree)→ other → …  becomes
                // attach →(tree_e)→ y →(new_e)→ other.
                assert(cs.cyc_next[attach] == other);
                cs.cyc_next[attach] = y;
                cs.cyc_prev[y] = attach;
                cs.cyc_next[y] = other;
                cs.cyc_prev[other] = y;
                cs.cyc_half[attach] = s.tree_half(attach, y);
                cs.cyc_half[y] =
                    (g.origin(2 * new_e) == y) ? 2 * new_e : 2 * new_e + 1;
            }
            cs.on_v[y] = true;
            cs.on_e[cs.nontree] = false;
            cs.on_e[tree_e] = true;
            cs.on_e[new_e] = true;
            if (attach == cs.lca && s.parent[cs.lca] == y) cs.lca = y;
            cs.cycle_count += 1;
            assert(cs.inside_cost >= s.cost(y));
            cs.inside_cost -= s.cost(y);
            cs.nontree = new_e;
            cs.inside_half = tri ^ 1;
            continue;
        }

        // ── Case 3c [LT79 tex 253–265 / 626–635] ──  Both (v*, y) and
        // (y, w*) are nontree; pick the candidate cycle with more cost
        // inside.
        PathResult pr = find_path_to_cycle(s, cs, y, markY, markL);
        std::size_t path_cost = 0;
        for (std::size_t v : pr.chain) {
            chain_mark[v] = true;
            path_cost += s.cost(v);
        }

        // Orient by the STORED traversal of the old nontree edge:
        // A →(nontree)→ B with A = origin(T), B = to(T); the tree part
        // of the cycle runs B → … → A in cyc_next order and z lies on
        // it.  Candidate through A keeps the arc [z ..cyc_next.. A] and
        // closes A →(edge (A,y))→ y →chain→ z; candidate through B
        // keeps [B ..cyc_next.. z] and closes z →chain→ y →(edge
        // (y,B))→ B.  The triangle halves are f2 = (w*→y), f3 = (y→v*),
        // so the candidate edges/halves depend on whether A is v* or w*.
        std::size_t T_old = cs.traversal_half(g);
        std::size_t A = g.origin(T_old), B = g.to[T_old];
        assert((A == vstar && B == wstar) || (A == wstar && B == vstar));
        std::size_t half_A_to_y = (A == vstar) ? (f3 ^ 1) : f2;
        std::size_t half_y_to_B = (B == wstar) ? (f2 ^ 1) : f3;
        std::size_t tri_A = (A == vstar) ? f3 : f2;   // triangle half of (A,y)
        std::size_t tri_B = (B == wstar) ? f2 : f3;   // triangle half of (y,B)
        assert(g.origin(half_A_to_y) == A && g.to[half_A_to_y] == y);
        assert(g.origin(half_y_to_B) == y && g.to[half_y_to_B] == B);

        // Scan side per the corner rule: the candidate's inside face is
        // F(twin(tri_half)); with candidate traversal half T_cand, the
        // inside is side-α iff twin(tri_half) == twin(T_cand).
        std::vector<std::size_t> chain2(pr.chain.rbegin(), pr.chain.rend());
        SideScan sa, sb;
        sa.s = &s; sa.cs = &cs; sa.chain_mark = &chain_mark;
        sa.arc_from = pr.z; sa.arc_to = A; sa.chain = &pr.chain;
        sa.new_half = half_A_to_y;
        sa.new_at_start = true;
        sa.scan_alpha = (tri_A == half_A_to_y);
        sb.s = &s; sb.cs = &cs; sb.chain_mark = &chain_mark;
        sb.arc_from = B; sb.arc_to = pr.z; sb.chain = &chain2;
        sb.new_half = half_y_to_B;
        sb.new_at_start = false;
        sb.scan_alpha = (tri_B == half_y_to_B);
        sa.begin();
        sb.begin();

        // Alternate scanning ([LT79 tex 629–632]).
        while (true) {
            if (!sa.step()) break;
            if (!sb.step()) break;
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

        // "Let (v_{i+1}, w_{i+1}) be the edge among (v_i, y) and
        // (y, w_i) whose cycle has more cost inside it" (tex 634–635).
        bool pickA = (inA >= inB);
        std::size_t new_e = pickA ? (half_A_to_y >> 1) : (half_y_to_B >> 1);
        std::size_t tri = pickA ? tri_A : tri_B;

        // Discard the old-cycle arc belonging to the other candidate,
        // unmarking as we go (each such vertex is strictly outside all
        // subsequent cycles — the accounting of [LT79 tex 656–667]).
        bool lca_dropped = false;
        auto discard = [&](std::size_t first, std::size_t last) {
            std::size_t v = first;
            while (true) {
                std::size_t nx = cs.cyc_next[v];
                cs.on_e[cs.cyc_half[v] >> 1] = false;
                if (v == cs.lca) lca_dropped = true;
                cs.on_v[v] = false;
                cs.cyc_next[v] = cs.cyc_prev[v] = cs.cyc_half[v] = NONE;
                cs.cycle_count -= 1;
                if (v == last) break;
                v = nx;
            }
        };
        cs.on_e[cs.nontree] = false;
        if (pickA) {
            // Keep [z .. A]; discard B..(vertex before z) and the edge
            // into z.
            if (B != pr.z) {
                std::size_t before_z = cs.cyc_prev[pr.z];
                discard(B, before_z);
            }
        } else {
            // Keep [B .. z]; discard (vertex after z)..A and the edge
            // out of z.
            if (A != pr.z) {
                std::size_t after_z = cs.cyc_next[pr.z];
                cs.on_e[cs.cyc_half[pr.z] >> 1] = false;
                discard(after_z, A);
            }
        }

        // Splice in the chain (+ new nontree edge) and remark.
        auto link = [&](std::size_t a, std::size_t b, std::size_t half) {
            cs.cyc_next[a] = b;
            cs.cyc_prev[b] = a;
            cs.cyc_half[a] = half;
            assert(g.origin(half) == a && g.to[half] == b);
            cs.on_e[half >> 1] = true;
        };
        if (pickA) {
            // A → y → chain … → z (traversal direction preserved on the
            // kept arc).
            if (pr.chain.empty()) {
                link(A, pr.z, half_A_to_y);              // A → y = z
            } else {
                link(A, pr.chain[0], half_A_to_y);       // A → y
                cs.on_v[pr.chain[0]] = true;
                cs.cycle_count += 1;
                for (std::size_t i = 0; i + 1 < pr.chain.size(); ++i) {
                    link(pr.chain[i], pr.chain[i + 1],
                         s.tree_half(pr.chain[i], pr.chain[i + 1]));
                    cs.on_v[pr.chain[i + 1]] = true;
                    cs.cycle_count += 1;
                }
                link(pr.chain.back(), pr.z,
                     s.tree_half(pr.chain.back(), pr.z));
            }
        } else {
            // z → chain(reversed) … → y → B.
            if (chain2.empty()) {
                link(pr.z, B, half_y_to_B);              // y = z → B
            } else {
                std::size_t prev = pr.z;
                for (std::size_t idx = 0; idx < chain2.size(); ++idx) {
                    std::size_t v = chain2[idx];
                    link(prev, v, s.tree_half(prev, v));
                    cs.on_v[v] = true;
                    cs.cycle_count += 1;
                    prev = v;
                }
                link(prev, B, half_y_to_B);              // y → B
            }
        }

        // Maintain the cycle's minimum-depth vertex.  In the
        // root-detour case (find_path_to_cycle) the chain climbs above
        // z = old lca, and the new minimum is the chain's shallowest
        // vertex; otherwise the minimum is the old lca if it survived,
        // else z (the top of the kept arc).
        std::size_t chain_min = NONE;
        for (std::size_t v : pr.chain)
            if (chain_min == NONE || s.depth[v] < s.depth[chain_min])
                chain_min = v;
        if (chain_min != NONE && s.depth[chain_min] < s.depth[pr.z]) {
            cs.lca = chain_min;
        } else if (lca_dropped) {
            cs.lca = pr.z;
        }

        for (std::size_t v : pr.chain) chain_mark[v] = false;

        cs.inside_cost = pickA ? inA : inB;
        cs.nontree = new_e;
        cs.inside_half = tri ^ 1;
    }
}

// Classify every shrunken vertex against the final cycle: 0 = the side
// containing the inside face, 1 = other side, 2 = on the cycle.
std::vector<std::uint8_t> classify_sides(const Shrunken& s,
                                         const CycleState& cs) {
    const HalfEdgeGraph& g = s.g;
    std::vector<std::uint8_t> side(g.n, 3);

    // Which rotation side is the inside?  inside_half == twin(T) ⟺ the
    // inside face lies on side-α (corner rule).
    std::size_t T = cs.traversal_half(g);
    bool inside_is_alpha = (cs.inside_half == (T ^ 1));

    for (std::size_t c = 0; c < g.n; ++c) {
        if (!cs.on_v[c]) continue;
        side[c] = 2;
        std::size_t out = cs.cyc_half[c];
        std::size_t in = cs.cyc_half[cs.cyc_prev[c]];
        rotation_between(g, out, in ^ 1, [&](std::size_t x) {
            std::size_t w = g.to[x];
            if (cs.on_v[w] || !s.is_tree_edge(x >> 1)) return;
            side[w] = inside_is_alpha ? 0 : 1;
        });
        rotation_between(g, in ^ 1, out, [&](std::size_t x) {
            std::size_t w = g.to[x];
            if (cs.on_v[w] || !s.is_tree_edge(x >> 1)) return;
            side[w] = inside_is_alpha ? 1 : 0;
        });
    }
    // The component of the tree above the cycle's minimum-depth vertex
    // (containing the apex when the apex is off the cycle) hangs under
    // the single crossing edge (lca, parent(lca)); the scan classified
    // parent(lca), and the apex shares its side.
    if (!cs.on_v[s.apex] && side[s.apex] == 3) {
        std::size_t pl = s.parent[cs.lca];
        assert(pl != NONE && side[pl] <= 1 &&
               "[LT79 tex 606–609]: the above-lca component is classified "
               "through the lca's parent edge");
        side[s.apex] = side[pl];
    }
    // Flood down the tree in BFS order (parents precede children).
    for (std::size_t v : s.bfs_order) {
        if (side[v] != 3) continue;
        assert(s.parent[v] != NONE && side[s.parent[v]] <= 1 &&
               "hanging subtrees inherit their crossing edge's side");
        side[v] = side[s.parent[v]];
    }

#ifndef NDEBUG
    std::size_t check_inside = 0;
    for (std::size_t v = 0; v < g.n; ++v)
        if (side[v] == 0) check_inside += s.cost(v);
    assert(check_inside == cs.inside_cost &&
           "[LT79 §3 tex 616–637]: incremental inside-cost bookkeeping "
           "must agree with a direct recount");
#endif
    return side;
}

// ════════════════════════════════════════════════════════════════
//  Step 6 — [LT79 §3 tex 575–586]
// ════════════════════════════════════════════════════════════════

// Shrink levels ≤ l0 to the apex, delete levels ≥ l2.  `level` is
// per piece-vertex (NONE off this component), l0 < l2.
Shrunken build_shrunken(const EmbeddedPlanarGraph& g,
                        const std::vector<std::size_t>& level,
                        std::size_t l0, std::size_t l2) {
    Shrunken s;
    auto in_cluster = [&](std::size_t v) {
        return level[v] != NONE && level[v] <= l0;
    };
    auto in_middle = [&](std::size_t v) {
        return level[v] != NONE && level[v] > l0 && level[v] < l2;
    };

    // Map middle vertices; apex = 0.
    std::vector<std::size_t> id(g.num_vertices(), NONE);
    s.g.add_vertex();                    // apex
    s.orig.push_back(NONE);
    for (std::size_t v = 0; v < g.num_vertices(); ++v) {
        if (!in_middle(v)) continue;
        id[v] = s.g.add_vertex();
        s.orig.push_back(v);
    }
    s.total = s.g.n - 1;
    assert(s.total >= 1 && "the shrink branch requires a nonempty middle");

    // [LT79 §3 tex 576–584]: "Construct a Boolean table with one entry
    // per vertex … Scan the edges incident to this tree clockwise
    // around the tree.  When scanning an edge (v, w) with v in the
    // tree, check the table entry for w.  If it is true, delete edge
    // (v, w).  If it is false, change it to true, construct an edge
    // (x, w), and delete edge (v, w)."  The clockwise scan is the
    // boundary walk of the cluster in the rotation system; a cluster
    // whose complement is disconnected in the embedding has several
    // boundary cycles — each is walked and their edge lists
    // concatenated (each enclosed pocket occupies its own wedge of the
    // apex rotation; the Euler assertion below validates the result).
    std::vector<bool> table(g.num_vertices(), false);
    std::vector<std::size_t> kept_half(g.num_vertices(), NONE);
    std::vector<std::size_t> apex_order;         // middle vertices in cw order
    std::vector<bool> half_seen(2 * g.num_edges(), false);

    for (std::size_t v = 0; v < g.num_vertices(); ++v) {
        if (!in_cluster(v)) continue;
        for (std::size_t h0 : g.incident_halves(v)) {
            if (half_seen[h0]) continue;
            if (in_cluster(g.half_to(h0))) continue;
            // Walk this boundary cycle.
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
                // Next boundary half: rotate clockwise around the
                // cluster, diving along cluster-internal edges.
                std::size_t j = g.rot_next(h);
                while (in_cluster(g.half_to(j)))
                    j = g.rot_next(j ^ 1);
                h = j;
            } while (h != h0);
        }
    }
    assert(!apex_order.empty() &&
           "the middle is BFS-reachable only through the cluster");

    // Build the half-edge structure.  Apex edges first, in scan order.
    std::vector<std::size_t> apex_edge(g.num_vertices(), NONE);
    std::vector<std::vector<std::size_t>> rot(s.g.n);
    for (std::size_t w : apex_order) {
        std::size_t half = s.g.raw_edge(0, id[w]);
        apex_edge[w] = half;             // the apex→w half
        rot[0].push_back(half);
    }
    std::vector<std::size_t> piece_edge_half(g.num_edges(), NONE);
    for (std::size_t v = 0; v < g.num_vertices(); ++v) {
        if (!in_middle(v)) continue;
        for (std::size_t h : g.incident_halves(v)) {
            std::size_t w = g.half_to(h);
            std::size_t e = h >> 1;
            if (in_middle(w)) {
                if (piece_edge_half[e] == NONE)
                    piece_edge_half[e] = s.g.raw_edge(id[v], id[w]);
                std::size_t mine = (s.g.origin(piece_edge_half[e]) == id[v])
                    ? piece_edge_half[e] : (piece_edge_half[e] ^ 1);
                rot[id[v]].push_back(mine);
            } else if (in_cluster(w)) {
                // Kept exactly when this is the half the scan retained.
                if (kept_half[v] == (h ^ 1))
                    rot[id[v]].push_back(apex_edge[v] ^ 1);
            }
            // deleted levels: dropped
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

    // [LT79 §3 tex 592–595] Step 7: BFS tree from the apex + descendant
    // costs.
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
            if (s.depth[w] != NONE) return;
            s.depth[w] = s.depth[v] + 1;
            s.parent[w] = v;
            s.parent_edge[w] = h >> 1;
            s.bfs_order.push_back(w);
        });
    }
    assert(s.bfs_order.size() == s.g.n &&
           "[LT79 §3 tex 592–595]: the shrunken graph is connected "
           "through the apex");
    // Radius bound of [LT79 tex 303–304]: "The new graph has a spanning
    // tree of radius l2 − l1 − 1" (with l1 := l0 here).
    for (std::size_t v = 0; v < s.g.n; ++v) {
        assert(s.depth[v] <= l2 - l0 - 1 &&
               "[LT79 tex 303–304]: shrunken BFS radius ≤ l2 − l0 − 1");
        (void)v;
    }
    for (std::size_t qi = s.bfs_order.size(); qi-- > 0; ) {
        std::size_t v = s.bfs_order[qi];
        s.desc[v] += s.cost(v);
        if (s.parent[v] != NONE) s.desc[s.parent[v]] += s.desc[v];
    }

    triangulate(s.g);
    return s;
}

// ════════════════════════════════════════════════════════════════
//  Connected-component partition — [LT79 thm:main tex 325–356]
// ════════════════════════════════════════════════════════════════

void lt_connected(const EmbeddedPlanarGraph& g,
                  const std::vector<std::size_t>& comp,
                  std::vector<SepPart>& part) {
    const std::size_t total_n = g.num_vertices();
    const std::size_t n_comp = comp.size();

    // Step 3 (tex 550–552): BFS levels within the component.
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
            if (level[w] != NONE) continue;
            level[w] = level[v] + 1;
            order.push_back(w);
        }
    }
    assert(order.size() == n_comp && "component must be BFS-connected");

    std::vector<std::size_t> L(r + 2, 0);        // L(r+1) = 0 [LT79 tex 278]
    for (std::size_t v : order) ++L[level[v]];

    // Step 4 (tex 558–561): l1 = first level with cost(0..l1) ≥ 1/2.
    std::size_t l1 = NONE, k = 0;
    {
        std::size_t pref = 0;
        for (std::size_t l = 0; l <= r; ++l) {
            pref += L[l];
            if (2 * pref >= total_n) { l1 = l; k = pref; break; }
        }
    }
    assert(l1 != NONE &&
           "[LT79 tex 331–334]: a component of cost > 2/3 > 1/2 always "
           "yields l1");

    // Step 5 (tex 567–569): highest l0 ≤ l1 with L(l0) + 2(l1−l0) ≤ 2√k,
    // lowest l2 ≥ l1+1 with L(l2) + 2(l2−l1−1) ≤ 2√(n_comp − k).
    auto sq_le = [](std::size_t a, std::size_t four_b) {
        return a * a <= four_b;          // a ≤ 2√b ⟺ a² ≤ 4b
    };
    std::size_t l0 = NONE;
    for (std::size_t l = l1 + 1; l-- > 0; )
        if (sq_le(L[l] + 2 * (l1 - l), 4 * k)) { l0 = l; break; }
    std::size_t l2 = NONE;
    for (std::size_t l = l1 + 1; l <= r + 1; ++l)
        if (sq_le(L[l] + 2 * (l - l1 - 1), 4 * (n_comp - k))) {
            l2 = l;
            break;
        }
    assert(l0 != NONE && l2 != NONE &&
           "[LT79 tex 345–356]: suitable levels l0 and l2 always exist");

    // Counts of the level bands.
    std::size_t below = 0, middle = 0, above = 0;
    for (std::size_t l = 0; l < l0; ++l) below += L[l];
    for (std::size_t l = l0 + 1; l < l2; ++l) middle += L[l];
    for (std::size_t l = l2 + 1; l <= r; ++l) above += L[l];

    auto band_of = [&](std::size_t v) -> int {
        std::size_t l = level[v];
        if (l == l0 || l == l2) return 4;        // separator levels
        if (l < l0) return 0;
        if (l < l2) return 1;
        return 2;
    };

    if (3 * middle <= 2 * total_n) {
        // [LT79 lem:levels tex 296–298]: "let A be the most costly part
        // of the three, let B be the remaining two parts, and let C be
        // the set of vertices on levels l1 [= l0] and l2."  That is the
        // component-internal partition (A*, B*).  [LT79 thm:main tex
        // 372–376] extends it to the whole graph: "Let A be the set
        // among A* and B* with greater cost ... and let B be the
        // remaining vertices of G" — without the max rule, B (= the two
        // lighter bands plus every other component) can exceed 2n/3
        // when the bands are balanced and the graph is disconnected.
        // Taking A = B* is edge-safe: the three bands are pairwise
        // non-adjacent, being separated by the deleted levels l0/l2
        // (BFS levels of adjacent vertices differ by at most 1).
        int star_band = (middle >= below && middle >= above) ? 1
                      : (below >= above ? 0 : 2);
        const std::size_t a_star =
            star_band == 0 ? below : star_band == 1 ? middle : above;
        const std::size_t b_star = below + middle + above - a_star;
        const bool a_is_pair = b_star > a_star;  // A = B* (two bands)
        for (std::size_t v : order) {
            int b = band_of(v);
            part[v] = (b == 4) ? SepPart::D
                    : ((b == star_band) != a_is_pair ? SepPart::A
                                                     : SepPart::B);
        }
        return;
    }

    // [LT79 lem:levels tex 300–313] + Steps 6–10: shrink and find the
    // cycle separator.
    Shrunken s = build_shrunken(g, level, l0, l2);

    // Step 8 (tex 603): "Choose any nontree edge."
    std::size_t e0 = NONE;
    for (std::size_t e = 0; e < s.g.num_edges(); ++e)
        if (!s.is_tree_edge(e)) { e0 = e; break; }
    assert(e0 != NONE &&
           "triangulated shrunken graph with ≥ 3 vertices has a nontree "
           "edge");

    CycleState cs;
    build_initial_cycle(s, e0, cs);
    improve_cycle(s, cs, total_n);

    std::vector<std::uint8_t> side = classify_sides(s, cs);

    // Cost on each side (real vertices only; the apex has cost 0 and is
    // ignored when mapping back — [LT79 tex 306–310]: "Let A be the set
    // among A* and B* having greater cost, let C consist of the
    // vertices on levels l1 [= l0] and l2 … plus the vertices in C*
    // minus the root of the tree, and let B contain the remaining
    // vertices").
    std::size_t in_cost = 0, out_cost = 0;
    for (std::size_t v = 1; v < s.g.n; ++v) {
        if (side[v] == 0) in_cost += 1;
        else if (side[v] == 1) out_cost += 1;
    }
    assert(in_cost == cs.inside_cost);
    assert(3 * in_cost <= 2 * total_n &&
           "[LT79 tex 616–618]: Step 9 exits with inside cost ≤ 2/3");
    assert(3 * out_cost <= 2 * total_n &&
           "[LT79 lem:radius tex 197–200]: the final cycle bounds both "
           "sides by 2/3");
    std::uint8_t a_side = (in_cost >= out_cost) ? 0 : 1;

    // Cycle-size bound [LT79 tex 341–342]:
    // |C| ≤ L(l0) + L(l2) + 2(l2 − l0 − 1) ≤ 2√k + 2√(n−k) ≤ 2√(2n).
    std::size_t cyc_real = 0;
    for (std::size_t v = 1; v < s.g.n; ++v)
        if (side[v] == 2) ++cyc_real;
    assert(cyc_real + 1 >= cs.cycle_count);  // apex on cycle or not
    assert(cyc_real <= 2 * (l2 - l0 - 1) &&
           "[LT79 lem:radius tex 206–208]: cycle length ≤ 2r + 1 with "
           "r = l2 − l0 − 1; the root is on the cycle (≤ 2r non-root "
           "vertices) or off it (the paths meet below the root, ≤ 2r − 1 "
           "total), so the non-apex count is at most 2r either way");

    for (std::size_t v : order) {
        int b = band_of(v);
        part[v] = (b == 4) ? SepPart::D : SepPart::B;
    }
    // Middle vertices take their shrunken side ([LT79 tex 306–310]).
    for (std::size_t sv = 1; sv < s.g.n; ++sv) {
        std::size_t v = s.orig[sv];
        if (side[sv] == 2) part[v] = SepPart::D;
        else part[v] = (side[sv] == a_side) ? SepPart::A : SepPart::B;
    }
}

} // namespace

// ════════════════════════════════════════════════════════════════
//  planar_separator — [LT79 thm:main tex 358–380 + Cor sqrt-sep]
// ════════════════════════════════════════════════════════════════

std::vector<SepPart> planar_separator(const EmbeddedPlanarGraph& g) {
    const std::size_t n = g.num_vertices();
    assert(n >= 1);
    std::vector<SepPart> part(n, SepPart::B);

    // Step 2 (tex 541–544): connected components.
    std::vector<std::size_t> comp_id(n, NONE);
    std::vector<std::vector<std::size_t>> comps;
    for (std::size_t v0 = 0; v0 < n; ++v0) {
        if (comp_id[v0] != NONE) continue;
        comps.emplace_back();
        auto& c = comps.back();
        comp_id[v0] = comps.size() - 1;
        c.push_back(v0);
        for (std::size_t qi = 0; qi < c.size(); ++qi)
            for (std::size_t h : g.incident_halves(c[qi])) {
                std::size_t w = g.half_to(h);
                if (comp_id[w] != NONE) continue;
                comp_id[w] = comps.size() - 1;
                c.push_back(w);
            }
    }

    std::size_t big = 0;
    for (std::size_t i = 1; i < comps.size(); ++i)
        if (comps[i].size() > comps[big].size()) big = i;

    if (3 * comps[big].size() > 2 * n) {
        // [LT79 tex 372–376]: recurse into the heavy component; A = the
        // heavier of its two sides, B = everything else.
        lt_connected(g, comps[big], part);
    } else if (3 * comps[big].size() > n) {
        // [LT79 tex 367–370]: a component with cost in (1/3, 2/3].
        for (std::size_t v : comps[big]) part[v] = SepPart::A;
    } else {
        // [LT79 tex 358–365]: all components ≤ 1/3; take a prefix of
        // components until its cost exceeds 1/3.
        std::size_t acc = 0;
        for (const auto& c : comps) {
            if (3 * acc > n) break;
            for (std::size_t v : c) part[v] = SepPart::A;
            acc += c.size();
        }
    }

    // Postconditions ([C91 §3.4 tex 297–302] (i)–(iii)).
#ifndef NDEBUG
    std::size_t na = 0, nb = 0, nd = 0;
    for (std::size_t v = 0; v < n; ++v) {
        if (part[v] == SepPart::A) ++na;
        else if (part[v] == SepPart::B) ++nb;
        else ++nd;
    }
    assert(3 * na <= 2 * n && "[C91 §3.4 tex 300] (ii): |A| ≤ 2μ/3");
    assert(3 * nb <= 2 * n && "[C91 §3.4 tex 300] (ii): |B| ≤ 2μ/3");
    assert(nd * nd <= 8 * n && "[C91 §3.4 tex 301] (iii): |D| ≤ √(8μ)");
    for (std::size_t e = 0; e < g.num_edges(); ++e) {
        SepPart pu = part[g.edge_u(e)], pv = part[g.edge_v(e)];
        assert(!((pu == SepPart::A && pv == SepPart::B) ||
                 (pu == SepPart::B && pv == SepPart::A)) &&
               "[C91 §3.4 tex 299] (i): no edge joins A and B");
    }
#endif
    return part;
}

// ════════════════════════════════════════════════════════════════
//  build_separator_decomposition — [C91 §3.4 tex 304]
// ════════════════════════════════════════════════════════════════

SeparatorDecomposition build_separator_decomposition(
        const EmbeddedPlanarGraph& g) {
    const std::size_t mu = g.num_vertices();
    SeparatorDecomposition out;
    out.subset.assign(mu, NONE);

    // Stop threshold: piece ≤ μ^δ with δ = 2/3 ⟺ piece³ ≤ μ².
    // 128-bit exact: sz³ overflows 64 bits for sz ≥ 2 642 246, a
    // realistic node count.
    auto is_leaf = [&](std::size_t sz) {
        __extension__ typedef unsigned __int128 u128;
        const u128 s = sz, m = mu;
        return s * s * s <= m * m;
    };

    struct Item {
        EmbeddedPlanarGraph graph;
        std::vector<std::size_t> ids;    // piece vertex → top-level node
    };
    std::vector<Item> stack;
    {
        std::vector<std::size_t> all(mu);
        for (std::size_t v = 0; v < mu; ++v) all[v] = v;
        std::vector<std::size_t> ids;
        stack.push_back(Item{g.induced(all, &ids), std::move(ids)});
    }

    while (!stack.empty()) {
        Item it = std::move(stack.back());
        stack.pop_back();
        std::size_t sz = it.graph.num_vertices();
        if (sz == 0) continue;
        if (is_leaf(sz)) {
            // [C91 §3.4 tex 304]: one subset D_i per remaining piece,
            // "each of size at most μ^{2/3}".
            for (std::size_t v = 0; v < sz; ++v)
                out.subset[it.ids[v]] = out.num_subsets;
            ++out.num_subsets;
            continue;
        }
        std::vector<SepPart> part = planar_separator(it.graph);
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
            if (nodes->empty()) continue;
            std::vector<std::size_t> sub_ids;
            EmbeddedPlanarGraph sub = it.graph.induced(*nodes, &sub_ids);
            std::vector<std::size_t> top_ids(nodes->size());
            for (std::size_t i = 0; i < nodes->size(); ++i)
                top_ids[i] = it.ids[(*nodes)[i]];
            stack.push_back(Item{std::move(sub), std::move(top_ids)});
        }
    }

#ifndef NDEBUG
    // [C91 §3.4 tex 304]: |D*| = O(μ^δ), δ = 2/3.  Concrete constant:
    // recursion nodes of size class j (sizes in [r·(3/2)^j, r·(3/2)^{j+1})
    // with r = μ^{2/3}) are vertex-disjoint, so there are at most
    // μ/(r·(3/2)^j) of them; each contributes ≤ √(8·size).  Summing the
    // geometric series: |D*| ≤ √8·√(3/2)/(1 − √(2/3)) · μ^{2/3}
    // < 19·μ^{2/3}.
    {
        long double bound = 19.0L;
        long double m23 = 1.0L;
        // μ^{2/3} via integer-safe cube comparison is awkward; a long
        // double estimate with a +1 slack is exact enough for the
        // assertion's purpose.
        m23 = powl((long double)mu, 2.0L / 3.0L);
        assert((long double)out.dstar_size <= bound * m23 + 1.0L &&
               "[C91 §3.4 tex 304]: |D*| = O(μ^{2/3}) with the derived "
               "constant");
    }
    // Every subset obeys the μ^{2/3} bound (checked structurally by
    // is_leaf, re-checked here).
    {
        std::vector<std::size_t> sizes(out.num_subsets, 0);
        for (std::size_t v = 0; v < mu; ++v)
            if (out.subset[v] != NONE) ++sizes[out.subset[v]];
        for (std::size_t szv : sizes) {
            // 128-bit exact, as in is_leaf: sz³ overflows 64 bits at
            // realistic node counts.
            __extension__ typedef unsigned __int128 u128;
            assert((u128)szv * szv * szv <= (u128)mu * mu &&
                   "[C91 §3.4 tex 304]: each |D_i| ≤ μ^{2/3}");
        }
    }
#endif
    return out;
}

} // namespace chazelle
