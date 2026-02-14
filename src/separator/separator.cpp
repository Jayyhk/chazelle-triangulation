#include "separator/separator.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <queue>
#include <vector>

namespace chazelle {

// ════════════════════════════════════════════════════════════════════
//  Internal helpers
// ════════════════════════════════════════════════════════════════════

namespace {

// ── Step 2: Connected components via BFS ────────────────────────

struct Component {
    std::vector<std::size_t> vertices;
    double total_cost = 0.0;
};

std::vector<Component>
find_components(const PlanarGraph& graph) {
    const std::size_t nv = graph.num_vertices();
    std::vector<bool> visited(nv, false);
    std::vector<Component> comps;

    for (std::size_t start = 0; start < nv; ++start) {
        if (visited[start] || graph.vertex(start).deleted) continue;
        Component comp;
        std::queue<std::size_t> q;
        q.push(start);
        visited[start] = true;
        while (!q.empty()) {
            std::size_t v = q.front();
            q.pop();
            comp.vertices.push_back(v);
            comp.total_cost += graph.vertex(v).cost;
            graph.for_each_edge_cw(v, [&](std::size_t ei) {
                std::size_t w = graph.edge(ei).other(v);
                if (!visited[w] && !graph.vertex(w).deleted) {
                    visited[w] = true;
                    q.push(w);
                }
            });
        }
        comps.push_back(std::move(comp));
    }
    return comps;
}

// ── Step 3: BFS spanning tree ───────────────────────────────────

/// Run BFS from `root`, setting level / parent / parent_edge / is_tree
/// on the graph.  Returns the maximum level.
std::size_t bfs_spanning_tree(PlanarGraph& graph,
                              const std::vector<std::size_t>& verts,
                              std::size_t root) {
    // Reset BFS fields for relevant vertices.
    for (std::size_t v : verts) {
        auto& vx = graph.vertex(v);
        vx.level       = NONE;
        vx.parent      = NONE;
        vx.parent_edge = NONE;
        vx.desc_cost   = vx.cost;
    }
    // Clear tree marks on edges incident to relevant vertices only.
    for (std::size_t v : verts) {
        graph.for_each_edge_cw(v, [&](std::size_t ei) {
            graph.edge(ei).is_tree = false;
        });
    }

    std::queue<std::size_t> q;
    graph.vertex(root).level = 0;
    q.push(root);
    std::size_t max_level = 0;

    while (!q.empty()) {
        std::size_t v = q.front();
        q.pop();
        std::size_t lv = graph.vertex(v).level;
        graph.for_each_edge_cw(v, [&](std::size_t ei) {
            std::size_t w = graph.edge(ei).other(v);
            if (graph.vertex(w).level != NONE || graph.vertex(w).deleted)
                return;
            graph.vertex(w).level       = lv + 1;
            graph.vertex(w).parent      = v;
            graph.vertex(w).parent_edge = ei;
            graph.edge(ei).is_tree      = true;
            max_level = std::max(max_level, lv + 1);
            q.push(w);
        });
    }
    return max_level;
}

// ── Step 7 helper: triangulate all faces ────────────────────────

/// Triangulate every face of the embedded graph by fan-triangulation.
/// Added edges are marked as non-tree.
void triangulate_faces(PlanarGraph& graph,
                       const std::vector<std::size_t>& verts) {
    // Collect all directed half-edges, then trace faces.
    // A directed half-edge is (edge_index, destination_vertex).
    // We track visited half-edges to avoid re-processing faces.
    struct HalfEdge {
        std::size_t edge_idx;
        std::size_t dest;
    };

    // Use a set of (edge_idx * 2 + side) to track visited directed edges.
    const std::size_t ne = graph.num_edges();
    std::vector<bool> visited(ne * 2, false);

    auto half_id = [&](std::size_t ei, std::size_t dest) -> std::size_t {
        int s = graph.edge(ei).side_of(dest);
        return ei * 2 + static_cast<std::size_t>(s);
    };

    for (std::size_t v : verts) {
        if (graph.vertex(v).deleted) continue;
        graph.for_each_edge_cw(v, [&](std::size_t start_ei) {
            // Directed edge: towards v.
            std::size_t hid = half_id(start_ei, v);
            if (visited[hid]) return;

            // Trace the face boundary.
            std::vector<HalfEdge> face;
            std::size_t ei = start_ei;
            std::size_t dest = v;
            do {
                std::size_t h = half_id(ei, dest);
                if (h < visited.size()) visited[h] = true;
                face.push_back({ei, dest});
                auto [nei, ndest] = graph.face_next(ei, dest);
                ei   = nei;
                dest = ndest;
            } while (ei != start_ei || dest != v);

            // If face has > 3 edges, triangulate via fan from face[0].dest.
            if (face.size() <= 3) return;

            std::size_t anchor = face[0].dest;
            for (std::size_t i = 2; i + 1 < face.size(); ++i) {
                std::size_t w = face[i].dest;
                // Add non-tree edge (anchor, w).
                // Find where to splice at anchor: after face[0].edge_idx.
                // Find where to splice at w: after face[i].edge_idx.
                std::size_t new_ei =
                    graph.add_edge(anchor, w, face[0].edge_idx,
                                   face[i].edge_idx);
                graph.edge(new_ei).is_tree = false;
                // Update visited capacity if edges grew.
                visited.resize(graph.num_edges() * 2, false);
            }
        });
    }
}

// ── Step 8–9 helpers: fundamental cycle operations ──────────────

/// Find the LCA of u and v in the BFS tree, and return the cycle
/// as a list of vertices from u to LCA to v (inclusive).
std::vector<std::size_t>
fundamental_cycle_vertices(const PlanarGraph& graph,
                           std::size_t u, std::size_t v) {
    // Walk both up to the same level.
    std::vector<std::size_t> path_u, path_v;
    std::size_t a = u, b = v;

    // Bring to same level.
    while (graph.vertex(a).level > graph.vertex(b).level) {
        path_u.push_back(a);
        a = graph.vertex(a).parent;
    }
    while (graph.vertex(b).level > graph.vertex(a).level) {
        path_v.push_back(b);
        b = graph.vertex(b).parent;
    }
    // Walk both up until meeting.
    while (a != b) {
        path_u.push_back(a);
        path_v.push_back(b);
        a = graph.vertex(a).parent;
        b = graph.vertex(b).parent;
    }
    path_u.push_back(a); // LCA

    // Cycle: path_u (u→LCA) + reverse of path_v (LCA→v).
    std::vector<std::size_t> cycle = std::move(path_u);
    for (auto it = path_v.rbegin(); it != path_v.rend(); ++it) {
        cycle.push_back(*it);
    }
    return cycle;
}

/// Compute the cost of vertices strictly inside a fundamental cycle.
/// Uses the BFS tree's descendant costs to avoid a full BFS.
///
/// For a fundamental cycle in a BFS tree formed by non-tree edge (u,v),
/// the cycle goes u→LCA→v along tree paths.  The "inside" vertices are
/// exactly the descendants of cycle vertices that hang off the interior
/// side.  We compute inside cost as:
///   sum over each cycle vertex w: desc_cost(w) − cost(w)
///     − sum of desc_cost(child) for children of w that are ON the cycle
/// This gives us the total cost of all descendants of the cycle that are
/// NOT on the cycle itself.
///
/// Returns: the raw (non-cycle descendant) cost.  Caller decides how
///          to interpret interior vs exterior.
double compute_one_side_cost(const PlanarGraph& graph,
                             const std::vector<std::size_t>& cycle,
                             const std::vector<bool>& on_cycle) {
    double desc_cost_sum = 0.0;
    for (std::size_t v : cycle) {
        double on_cycle_child_cost = 0.0;
        graph.for_each_edge_cw(v, [&](std::size_t ei) {
            if (!graph.edge(ei).is_tree) return;
            std::size_t w = graph.edge(ei).other(v);
            if (graph.vertex(w).parent == v && on_cycle[w]) {
                on_cycle_child_cost += graph.vertex(w).desc_cost;
            }
        });
        double off = graph.vertex(v).desc_cost - graph.vertex(v).cost
                   - on_cycle_child_cost;
        if (off > 0.0) desc_cost_sum += off;
    }
    return desc_cost_sum;
}

// ── Compute descendant costs bottom-up ──────────────────────────

void compute_descendant_costs(PlanarGraph& graph,
                              const std::vector<std::size_t>& verts,
                              std::size_t max_level) {
    // Bucket vertices by level, then process bottom-up.
    std::vector<std::vector<std::size_t>> by_level(max_level + 1);
    for (std::size_t v : verts) {
        if (graph.vertex(v).level != NONE) {
            by_level[graph.vertex(v).level].push_back(v);
        }
    }
    for (auto it = by_level.rbegin(); it != by_level.rend(); ++it) {
        for (std::size_t v : *it) {
            auto& vx = graph.vertex(v);
            vx.desc_cost = vx.cost;
            // Add descendant costs of tree children.
            graph.for_each_edge_cw(v, [&](std::size_t ei) {
                if (!graph.edge(ei).is_tree) return;
                std::size_t w = graph.edge(ei).other(v);
                if (graph.vertex(w).parent == v) {
                    vx.desc_cost += graph.vertex(w).desc_cost;
                }
            });
        }
    }
}

} // anonymous namespace

// ════════════════════════════════════════════════════════════════════
//  find_separator — Lipton-Tarjan 10-step algorithm
// ════════════════════════════════════════════════════════════════════

SeparatorResult find_separator(PlanarGraph& graph) {
    const std::size_t nv = graph.num_vertices();
    SeparatorResult result;

    // Collect non-deleted vertices.
    std::vector<std::size_t> all_verts;
    all_verts.reserve(nv);
    double total_cost = 0.0;
    for (std::size_t i = 0; i < nv; ++i) {
        if (!graph.vertex(i).deleted) {
            all_verts.push_back(i);
            total_cost += graph.vertex(i).cost;
        }
    }

    if (all_verts.empty()) return result;

    // Normalise costs to sum to 1.
    if (total_cost > 0.0) {
        for (std::size_t v : all_verts) {
            graph.vertex(v).cost /= total_cost;
        }
    }

    // ── Step 2: connected components ────────────────────────────
    auto comps = find_components(graph);

    // Find the heaviest component.
    std::size_t heavy_idx = 0;
    for (std::size_t i = 1; i < comps.size(); ++i) {
        if (comps[i].total_cost > comps[heavy_idx].total_cost)
            heavy_idx = i;
    }

    if (comps[heavy_idx].total_cost <= 2.0 / 3.0) {
        // Trivial separator: C = ∅, greedily pack A until cost > 1/3.
        double a_cost = 0.0;
        bool filling_a = true;
        for (std::size_t ci = 0; ci < comps.size(); ++ci) {
            auto& target = (filling_a ? result.A : result.B);
            for (std::size_t v : comps[ci].vertices) {
                target.push_back(v);
            }
            a_cost += comps[ci].total_cost;
            if (filling_a && a_cost > 1.0 / 3.0) {
                filling_a = false;
            }
        }
        // Restore original costs before returning.
        if (total_cost > 0.0) {
            for (std::size_t v : all_verts) {
                graph.vertex(v).cost *= total_cost;
            }
        }
        return result;
    }

    // Work on the heaviest component.
    auto& comp = comps[heavy_idx];

    // Put all other components into B.
    for (std::size_t ci = 0; ci < comps.size(); ++ci) {
        if (ci == heavy_idx) continue;
        for (std::size_t v : comps[ci].vertices)
            result.B.push_back(v);
    }

    // ── Step 3: BFS spanning tree ───────────────────────────────
    std::size_t root = comp.vertices[0];
    std::size_t max_level = bfs_spanning_tree(graph, comp.vertices, root);

    // Compute level sizes and cumulative costs.
    std::vector<double> level_cost(max_level + 1, 0.0);
    std::vector<std::size_t> level_size(max_level + 1, 0);
    for (std::size_t v : comp.vertices) {
        std::size_t lv = graph.vertex(v).level;
        level_cost[lv] += graph.vertex(v).cost;
        level_size[lv]++;
    }

    // ── Step 4: find median level l₁ ────────────────────────────
    std::size_t l1 = 0;
    {
        double cum = 0.0;
        for (std::size_t l = 0; l <= max_level; ++l) {
            cum += level_cost[l];
            if (cum > 0.5 * comp.total_cost) {
                l1 = l;
                break;
            }
        }
    }

    double k_cost = 0.0; // cost of levels 0..l1
    std::size_t k_size = 0;
    for (std::size_t l = 0; l <= l1; ++l) {
        k_cost += level_cost[l];
        k_size += level_size[l];
    }

    double n_comp = static_cast<double>(comp.vertices.size());

    // ── Step 5: find l₀ ≤ l₁ and l₂ ≥ l₁+1 ────────────────────
    // l₀: highest level ≤ l₁ such that L(l₀) + 2(l₁ − l₀) ≤ 2√k
    std::size_t l0 = 0;
    {
        double bound = 2.0 * std::sqrt(static_cast<double>(k_size));
        for (std::size_t l = l1; ; --l) {
            double val = static_cast<double>(level_size[l])
                       + 2.0 * static_cast<double>(l1 - l);
            if (val <= bound) {
                l0 = l;
                break;
            }
            // The paper proves l₀ must exist (proof by contradiction).
            assert(l > 0 && "L-T Step 5: level l0 must exist");
        }
    }

    // l₂: lowest level ≥ l₁+1 such that L(l₂) + 2(l₂ − l₁ − 1) ≤ 2√(n−k)
    std::size_t l2 = max_level;
    {
        double n_minus_k = n_comp - static_cast<double>(k_size);
        double bound = 2.0 * std::sqrt(n_minus_k);
        bool found_l2 = false;
        for (std::size_t l = l1 + 1; l <= max_level; ++l) {
            double val = static_cast<double>(level_size[l])
                       + 2.0 * static_cast<double>(l - l1 - 1);
            if (val <= bound) {
                l2 = l;
                found_l2 = true;
                break;
            }
        }
        // The paper proves l₂ must exist (proof by contradiction).
        assert(found_l2 && "L-T Step 5: level l2 must exist");
    }

    // Check if the middle region (levels l₀+1 .. l₂−1) is small enough.
    double middle_cost = 0.0;
    std::size_t middle_size = 0;
    for (std::size_t l = l0 + 1; l < l2; ++l) {
        middle_cost += level_cost[l];
        middle_size += level_size[l];
    }

    // Separator = levels l₀ and l₂.
    std::vector<std::size_t> sep_levels;
    for (std::size_t v : comp.vertices) {
        std::size_t lv = graph.vertex(v).level;
        if (lv == l0 || lv == l2) {
            sep_levels.push_back(v);
        }
    }

    if (middle_cost <= 2.0 / 3.0 * comp.total_cost) {
        // The level separator suffices — no cycle needed.
        // C = levels l₀ ∪ l₂
        // Split top / middle / bottom among A and B.
        result.C.insert(result.C.end(), sep_levels.begin(), sep_levels.end());
        for (std::size_t v : comp.vertices) {
            std::size_t lv = graph.vertex(v).level;
            if (lv == l0 || lv == l2) continue; // already in C
            if (lv < l0)        result.A.push_back(v);
            else if (lv > l2)   result.A.push_back(v);
            else                result.B.push_back(v); // l0 < lv < l2
        }
        // Balance: if A is too heavy, swap.
        double a_cost = 0.0, b_cost = 0.0;
        for (std::size_t v : result.A) a_cost += graph.vertex(v).cost;
        for (std::size_t v : result.B) b_cost += graph.vertex(v).cost;
        if (a_cost > 2.0 / 3.0 + 1e-9) {
            std::swap(result.A, result.B);
        }

        // Restore original costs before returning.
        if (total_cost > 0.0) {
            for (std::size_t v : all_verts) {
                graph.vertex(v).cost *= total_cost;
            }
        }
        return result;
    }

    // ── Steps 6–9: need a cycle separator in the middle region ──

    // Step 6: Per Lipton-Tarjan §3 Step 6, contract levels 0..l₀
    // into a single supervertex x and delete levels ≥ l₂.  Edges
    // from the contracted subtree to middle vertices (l₀+1..l₂-1)
    // are redirected through x via an Euler tour of the subtree,
    // preserving the planar embedding.  The BFS tree from x then has
    // depth ≤ l₂ − l₀ − 1, bounding fundamental cycles.
    std::vector<std::pair<std::size_t, bool>> step6_saved;  // vertex states
    std::vector<std::pair<std::size_t, bool>> step6_edge_saved; // edge states

    // In a connected BFS tree every level 0..max_level has vertices.
    // The paper proves l₀ exists; level l₀ is therefore non-empty.
    std::size_t super_x = NONE;
    for (std::size_t v : comp.vertices) {
        if (graph.vertex(v).level == l0) {
            super_x = v;
            break;
        }
    }
    assert(super_x != NONE && "BFS level l0 must have vertices");

    // 6a. Soft-delete vertices at levels ≥ l₂.
    for (std::size_t v : comp.vertices) {
        if (graph.vertex(v).level >= l2) {
            step6_saved.emplace_back(v, graph.vertex(v).deleted);
            graph.vertex(v).deleted = true;
        }
    }

    // 6b. Euler tour of the BFS tree at levels 0..l₀, per the paper:
    //     "Construct a Boolean table … Scan the edges incident to
    //      this tree clockwise around the tree."
    //     For each edge (v,w) with v in the tree and table[w] false,
    //     record w as needing a redirect edge from x.
    struct RedirectInfo {
        std::size_t w;          // middle vertex
        std::size_t pos_at_w;   // edge CW-before (v,w) at w
    };
    std::vector<RedirectInfo> redirects;
    std::vector<bool> table(nv, false);
    for (std::size_t v : comp.vertices) {
        if (graph.vertex(v).level <= l0)
            table[v] = true;  // in contracted set
    }

    auto tree_scan = [&](auto&& self, std::size_t v) -> void {
        graph.for_each_edge_cw(v, [&](std::size_t ei) {
            std::size_t w = graph.edge(ei).other(v);
            if (graph.vertex(w).deleted) return;  // l₂+ vertex
            if (table[w]) {
                // In contracted set or already connected to x.
                // Recurse on tree children within levels 0..l₀.
                if (graph.vertex(w).level <= l0 &&
                    graph.edge(ei).is_tree &&
                    graph.vertex(w).parent == v)
                {
                    self(self, w);
                }
            } else {
                // First encounter with middle vertex w.
                table[w] = true;
                int ws = graph.edge(ei).side_of(w);
                redirects.push_back({w, graph.edge(ei).ccw[ws]});
            }
        });
    };
    tree_scan(tree_scan, root);

    // 6c. Soft-delete ALL edges of vertices at levels 0..l₀.
    for (std::size_t v : comp.vertices) {
        if (graph.vertex(v).level > l0) continue;
        graph.for_each_edge_cw(v, [&](std::size_t ei) {
            step6_edge_saved.emplace_back(ei, graph.edge(ei).deleted);
            graph.edge(ei).deleted = true;
        });
    }

    // 6d. Add redirect edges from super_x to each middle vertex
    //     discovered during the Euler tour, in tour order.
    std::size_t edges_before_redirect = graph.num_edges();
    {
        std::size_t last_x_edge = NONE;
        for (auto& [w, pos_at_w] : redirects) {
            std::size_t new_ei = graph.add_edge(
                super_x, w, last_x_edge, pos_at_w);
            last_x_edge = new_ei;
        }
    }

    // 6e. Soft-delete 0..l₀ vertices except super_x; zero its cost
    //     per the paper ("shrink to a single vertex of cost zero").
    double saved_super_x_cost = graph.vertex(super_x).cost;
    graph.vertex(super_x).cost = 0.0;
    for (std::size_t v : comp.vertices) {
        std::size_t lv = graph.vertex(v).level;
        if ((lv < l0) || (lv == l0 && v != super_x)) {
            step6_saved.emplace_back(v, graph.vertex(v).deleted);
            graph.vertex(v).deleted = true;
        }
    }

    // Collect middle vertices (including super_x).
    std::vector<std::size_t> middle_verts;
    middle_verts.push_back(super_x);
    for (std::size_t v : comp.vertices) {
        std::size_t lv = graph.vertex(v).level;
        if (lv > l0 && lv < l2) {
            middle_verts.push_back(v);
        }
    }

    // Step 7: Rebuild BFS from super_x (depth ≤ l₂ − l₀ − 1),
    // compute descendant costs, and triangulate faces.
    std::size_t middle_max_level =
        bfs_spanning_tree(graph, middle_verts, super_x);
    compute_descendant_costs(graph, middle_verts, middle_max_level);
    triangulate_faces(graph, middle_verts);

    // Steps 8–9: Find a fundamental cycle that separates the middle.
    // Per Step 8, choose any nontree edge in the contracted graph.
    // After the middle-region BFS, levels have been reassigned, so
    // we identify valid edges by checking for non-deleted neighbors.
    std::size_t cycle_edge = NONE;
    for (std::size_t v : middle_verts) {
        graph.for_each_edge_cw(v, [&](std::size_t ei) {
            if (cycle_edge != NONE) return;
            if (graph.edge(ei).is_tree) return;
            std::size_t w = graph.edge(ei).other(v);
            if (!graph.vertex(w).deleted) {
                cycle_edge = ei;
            }
        });
        if (cycle_edge != NONE) break;
    }

    std::vector<std::size_t> cycle_verts;
    if (cycle_edge != NONE) {
        std::size_t cu = graph.edge(cycle_edge).endpoint[0];
        std::size_t cv = graph.edge(cycle_edge).endpoint[1];
        cycle_verts = fundamental_cycle_vertices(graph, cu, cv);

        // Step 9: improve the cycle.
        //
        // Lipton-Tarjan Step 9: incremental O(V) cycle improvement.
        // Each iteration absorbs one interior vertex into the cycle,
        // removing one face from the inside.  Cost update is O(1)
        // per iteration.  Total O(V) iterations × O(1) = O(V).

        // Mark cycle vertices in a vector<bool> for O(1) lookup.
        std::vector<bool> on_cycle(nv, false);
        for (std::size_t v : cycle_verts) on_cycle[v] = true;

        // Compute initial cycle cost and raw descendant cost.
        double cycle_cost_val = 0.0;
        for (std::size_t v : cycle_verts) {
            cycle_cost_val += graph.vertex(v).cost;
        }
        double one_side_desc =
            compute_one_side_cost(graph, cycle_verts, on_cycle);
        double other_side_desc =
            comp.total_cost - cycle_cost_val - one_side_desc;
        double inside_cost =
            std::max(one_side_desc, other_side_desc);

        // Determine interior vs exterior by marking one side via
        // BFS from one face of the initial non-tree edge.  O(V).
        std::vector<bool> is_interior(nv, false);
        double interior_cost = one_side_desc;
        {
            std::size_t cu0 = graph.edge(cycle_edge).endpoint[0];
            std::size_t cv0 = graph.edge(cycle_edge).endpoint[1];
            auto [ne_t, nd_t] = graph.face_next(cycle_edge, cv0);
            std::size_t y_probe = nd_t;

            // BFS from y_probe, not crossing cycle.
            if (!on_cycle[y_probe] && !graph.vertex(y_probe).deleted) {
                std::queue<std::size_t> bfs_q;
                is_interior[y_probe] = true;
                bfs_q.push(y_probe);
                while (!bfs_q.empty()) {
                    std::size_t bv = bfs_q.front(); bfs_q.pop();
                    graph.for_each_edge_cw(bv, [&](std::size_t ei) {
                        std::size_t w = graph.edge(ei).other(bv);
                        if (on_cycle[w] || is_interior[w] ||
                            graph.vertex(w).deleted) return;
                        is_interior[w] = true;
                        bfs_q.push(w);
                    });
                }
            }
            // Verify by comparing marked cost with one_side_desc.
            double marked_cost = 0.0;
            for (std::size_t v : middle_verts) {
                if (is_interior[v]) marked_cost += graph.vertex(v).cost;
            }
            // If marked set matches the OTHER side, flip.
            if (std::abs(marked_cost - other_side_desc) <
                std::abs(marked_cost - one_side_desc)) {
                // Guessed the wrong side — re-mark from the other face.
                is_interior.assign(nv, false);
                auto [ne_t2, nd_t2] = graph.face_next(cycle_edge, cu0);
                std::size_t y_probe2 = nd_t2;
                if (!on_cycle[y_probe2] &&
                    !graph.vertex(y_probe2).deleted) {
                    std::queue<std::size_t> bfs_q;
                    is_interior[y_probe2] = true;
                    bfs_q.push(y_probe2);
                    while (!bfs_q.empty()) {
                        std::size_t bv = bfs_q.front(); bfs_q.pop();
                        graph.for_each_edge_cw(bv, [&](std::size_t ei) {
                            std::size_t w = graph.edge(ei).other(bv);
                            if (on_cycle[w] || is_interior[w] ||
                                graph.vertex(w).deleted) return;
                            is_interior[w] = true;
                            bfs_q.push(w);
                        });
                    }
                }
                // Swap interior/exterior costs to match new marking.
                interior_cost = other_side_desc;
            }
        }

        double exterior_cost = comp.total_cost - cycle_cost_val
                             - interior_cost;

        std::size_t current_edge = cycle_edge;
        int max_iters = static_cast<int>(comp.vertices.size());

        while (inside_cost > 2.0 / 3.0 * comp.total_cost
               && max_iters-- > 0)
        {
            // Find the triangle on the inside of current_edge.
            std::size_t cu2 = graph.edge(current_edge).endpoint[0];
            std::size_t cv2 = graph.edge(current_edge).endpoint[1];

            auto [ne1, nd1] = graph.face_next(current_edge, cv2);
            std::size_t y1 = nd1;
            auto [ne2, nd2] = graph.face_next(current_edge, cu2);
            std::size_t y2 = nd2;

            // Pick the interior vertex.
            std::size_t y = NONE;
            if (is_interior[y1]) {
                y = y1;
            } else if (is_interior[y2]) {
                y = y2;
            }
            if (y == NONE) break; // No interior vertex to absorb.

            // Find a non-tree edge among (cu2,y) and (cv2,y).
            std::size_t nontree_edge = NONE;
            graph.for_each_edge_cw(y, [&](std::size_t ei) {
                if (nontree_edge != NONE) return;
                if (graph.edge(ei).is_tree) return;
                std::size_t w = graph.edge(ei).other(y);
                if (w == cu2 || w == cv2) {
                    nontree_edge = ei;
                }
            });

            if (nontree_edge == NONE) break;

            // §L-T Step 9: Absorb vertex y into the cycle.
            // Incremental cost update: O(1).
            on_cycle[y] = true;
            is_interior[y] = false;
            cycle_cost_val += graph.vertex(y).cost;
            interior_cost -= graph.vertex(y).cost;
            if (interior_cost < 0.0) interior_cost = 0.0;
            exterior_cost = comp.total_cost - cycle_cost_val
                          - interior_cost;
            inside_cost = std::max(interior_cost, exterior_cost);

            current_edge = nontree_edge;
        }

        // Reconstruct cycle_verts from on_cycle marking.
        cycle_verts.clear();
        for (std::size_t v : middle_verts) {
            if (on_cycle[v]) cycle_verts.push_back(v);
        }
    }

    // ── Step 10: assemble partition ─────────────────────────────
    // Clean up redirect edges and triangulation edges.
    for (std::size_t ei = edges_before_redirect; ei < graph.num_edges(); ++ei) {
        graph.edge(ei).deleted = true;
    }

    // Restore edges soft-deleted during contraction.
    for (auto& [ei, was_deleted] : step6_edge_saved) {
        graph.edge(ei).deleted = was_deleted;
    }

    // Restore vertices deleted in Step 6 and super_x cost.
    graph.vertex(super_x).cost = saved_super_x_cost;
    for (auto& [v, was_deleted] : step6_saved) {
        graph.vertex(v).deleted = was_deleted;
    }

    // C = separator levels (l₀, l₂) + cycle vertices
    std::vector<bool> is_sep(nv, false);
    for (std::size_t v : sep_levels) {
        if (!is_sep[v]) { is_sep[v] = true; result.C.push_back(v); }
    }
    for (std::size_t v : cycle_verts) {
        if (!is_sep[v]) { is_sep[v] = true; result.C.push_back(v); }
    }

    // Partition remaining vertices into A (inside cycle) and B (outside).
    // BFS from tree root, not crossing separator.
    std::size_t tree_root = root;
    std::vector<bool> reachable(nv, false);
    {
        std::queue<std::size_t> q;
        if (!is_sep[tree_root]) {
            reachable[tree_root] = true;
            q.push(tree_root);
        }
        while (!q.empty()) {
            std::size_t v = q.front();
            q.pop();
            graph.for_each_edge_cw(v, [&](std::size_t ei) {
                std::size_t w = graph.edge(ei).other(v);
                if (reachable[w] || graph.vertex(w).deleted) return;
                if (is_sep[w]) return;
                reachable[w] = true;
                q.push(w);
            });
        }
    }

    for (std::size_t v : comp.vertices) {
        if (is_sep[v]) continue;
        if (reachable[v])
            result.A.push_back(v);
        else
            result.B.push_back(v);
    }

    // Ensure balance: if one side > 2/3, swap labels.
    double a_cost = 0.0, b_cost = 0.0;
    for (std::size_t v : result.A) a_cost += graph.vertex(v).cost;
    for (std::size_t v : result.B) b_cost += graph.vertex(v).cost;
    if (a_cost > 2.0 / 3.0 + 1e-9 && b_cost <= 2.0 / 3.0 + 1e-9) {
        std::swap(result.A, result.B);
    }

    // Restore original costs.
    for (std::size_t v : all_verts) {
        graph.vertex(v).cost *= total_cost;
    }

    return result;
}

// ════════════════════════════════════════════════════════════════════
//  iterated_separator
// ════════════════════════════════════════════════════════════════════

IteratedSeparatorResult
iterated_separator(PlanarGraph& graph, std::size_t max_piece_size) {
    const std::size_t nv = graph.num_vertices();
    IteratedSeparatorResult result;
    result.separator_nodes.assign(nv, false);
    result.vertex_piece.assign(nv, NONE);

    // Start with all non-deleted vertices as a single piece.
    std::vector<std::size_t> initial;
    for (std::size_t i = 0; i < nv; ++i) {
        if (!graph.vertex(i).deleted) initial.push_back(i);
    }

    // DFS-based recursion: at each call only `piece` vertices are
    // non-deleted (guaranteed by the caller).  Soft-delete A or B
    // in O(|piece|) per call, not O(nv).  Total work per recursion
    // level is O(μ) (pieces are disjoint) → O(μ log μ) overall.
    //
    // We use a vector<bool> to track which vertices are part of
    // the current piece for O(1) membership tests.
    std::vector<bool> in_piece(nv, false);

    std::function<void(std::vector<std::size_t>)> recurse =
        [&](std::vector<std::size_t> piece) {
        if (piece.size() <= max_piece_size) {
            // Small enough — record as a final piece.
            std::size_t pidx = result.pieces.size();
            for (std::size_t v : piece) {
                result.vertex_piece[v] = pidx;
            }
            result.pieces.push_back(std::move(piece));
            return;
        }

        // Mark piece vertices in the shared boolean vector.
        for (std::size_t v : piece) in_piece[v] = true;

        // Reset costs to uniform.
        for (std::size_t v : piece) {
            graph.vertex(v).cost = 1.0;
        }

        auto sep = find_separator(graph);

        // Clear the piece marking.
        for (std::size_t v : piece) in_piece[v] = false;

        // Mark separator vertices.
        for (std::size_t v : sep.C) {
            result.separator_nodes[v] = true;
            graph.vertex(v).deleted = true; // remove from both sides
        }

        // Build sub-pieces (excluding separator vertices).
        std::vector<std::size_t> a_piece, b_piece;
        for (std::size_t v : sep.A) {
            if (!result.separator_nodes[v]) a_piece.push_back(v);
        }
        for (std::size_t v : sep.B) {
            if (!result.separator_nodes[v]) b_piece.push_back(v);
        }

        // Recurse on A: soft-delete B in O(|B|), then restore.
        if (!a_piece.empty()) {
            for (std::size_t v : b_piece)
                graph.vertex(v).deleted = true;
            recurse(std::move(a_piece));
            for (std::size_t v : b_piece)
                graph.vertex(v).deleted = false;
        }

        // Recurse on B: soft-delete A in O(|A|), then restore.
        // (A was already consumed by move, but sep.A still has the IDs.)
        if (!b_piece.empty()) {
            for (std::size_t v : sep.A) {
                if (!result.separator_nodes[v])
                    graph.vertex(v).deleted = true;
            }
            recurse(std::move(b_piece));
            for (std::size_t v : sep.A) {
                if (!result.separator_nodes[v])
                    graph.vertex(v).deleted = false;
            }
        }

        // Restore separator vertices (un-delete for parent's use).
        for (std::size_t v : sep.C) {
            graph.vertex(v).deleted = false;
        }
    };

    recurse(std::move(initial));

    return result;
}

} // namespace chazelle
