#include "separator/separator.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <queue>
#include <unordered_set>
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
///   sum over each cycle vertex w (except LCA): desc_cost(w) − cost(w)
///     − sum of desc_cost(child) for children of w that are ON the cycle
/// This gives us the total cost of all descendants of the cycle that are
/// NOT on the cycle itself, which is exactly the interior cost.
double compute_inside_cost(const PlanarGraph& graph,
                           const std::vector<std::size_t>& cycle,
                           const std::unordered_set<std::size_t>& cycle_set,
                           double total_cost) {
    // Cost of cycle vertices themselves.
    double cycle_cost = 0.0;
    for (std::size_t v : cycle) {
        cycle_cost += graph.vertex(v).cost;
    }

    // For each cycle vertex, sum the descendant costs of tree children
    // that are also on the cycle.  The difference between desc_cost and
    // (cost + on-cycle-children desc_cost) gives the off-cycle descendants,
    // which are the interior vertices.
    //
    // In a fundamental cycle u→LCA→v, every cycle vertex except the LCA
    // has exactly one cycle-child (the next vertex on the path toward u
    // or v).  The LCA has two cycle-children (one from each path).
    double interior_desc_cost = 0.0;

    for (std::size_t v : cycle) {
        double on_cycle_child_cost = 0.0;
        graph.for_each_edge_cw(v, [&](std::size_t ei) {
            if (!graph.edge(ei).is_tree) return;
            std::size_t w = graph.edge(ei).other(v);
            if (graph.vertex(w).parent == v && cycle_set.count(w)) {
                on_cycle_child_cost += graph.vertex(w).desc_cost;
            }
        });
        // Interior descendants of v = desc_cost(v) − cost(v) − on_cycle_child_cost
        double off = graph.vertex(v).desc_cost - graph.vertex(v).cost
                   - on_cycle_child_cost;
        if (off > 0.0) interior_desc_cost += off;
    }

    // interior_desc_cost is one side.  The other side is:
    double other_side = total_cost - cycle_cost - interior_desc_cost;

    // "Inside" = the heavier side (the side we need to reduce).
    return std::max(interior_desc_cost, other_side);
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
            if (l == 0) { l0 = 0; break; }
        }
    }

    // l₂: lowest level ≥ l₁+1 such that L(l₂) + 2(l₂ − l₁ − 1) ≤ 2√(n−k)
    std::size_t l2 = max_level;
    {
        double n_minus_k = n_comp - static_cast<double>(k_size);
        double bound = 2.0 * std::sqrt(n_minus_k);
        for (std::size_t l = l1 + 1; l <= max_level; ++l) {
            double val = static_cast<double>(level_size[l])
                       + 2.0 * static_cast<double>(l - l1 - 1);
            if (val <= bound) {
                l2 = l;
                break;
            }
        }
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

    // Step 6: Lipton-Tarjan graph shrinking.
    // Contract levels 0..l₀ into a single supervertex x, and delete
    // levels ≥ l₂.  Then rebuild BFS tree in the shrunken graph rooted
    // at x, which has depth ≤ l₂ − l₀ − 1.  This bounds the
    // fundamental cycle length to ≤ 2(l₂ − l₀ − 1) + 1 ≤ 2√(2n).
    //
    // Implementation: pick an arbitrary level-l₀ vertex as the
    // supervertex x (it inherits all edges from levels 0..l₀ to
    // levels l₀+1..l₂−1).  Soft-delete all other vertices outside
    // [l₀, l₂−1].  Then rebuild BFS from x on the middle region.
    std::vector<std::pair<std::size_t, bool>> step6_saved;

    // Find a representative vertex at level l₀ to serve as supervertex x.
    std::size_t super_x = NONE;
    for (std::size_t v : comp.vertices) {
        if (graph.vertex(v).level == l0) {
            super_x = v;
            break;
        }
    }
    if (super_x == NONE) super_x = root; // fallback

    // Soft-delete all vertices outside the middle region [l₀, l₂−1],
    // except the supervertex x which stays at level l₀.
    for (std::size_t v : comp.vertices) {
        std::size_t lv = graph.vertex(v).level;
        if (lv < l0 || lv >= l2) {
            step6_saved.emplace_back(v, graph.vertex(v).deleted);
            graph.vertex(v).deleted = true;
        } else if (lv == l0 && v != super_x) {
            // Other l₀ vertices are also contracted into super_x.
            step6_saved.emplace_back(v, graph.vertex(v).deleted);
            graph.vertex(v).deleted = true;
        }
    }

    // Collect middle vertices (including the supervertex).
    std::vector<std::size_t> middle_verts;
    middle_verts.push_back(super_x);
    for (std::size_t v : comp.vertices) {
        std::size_t lv = graph.vertex(v).level;
        if (lv > l0 && lv < l2) {
            middle_verts.push_back(v);
        }
    }

    // Step 6 continued: Rebuild BFS spanning tree from super_x on the
    // middle region.  This gives a new tree with depth ≤ l₂ − l₀ − 1,
    // which is critical for bounding separator size.
    std::size_t middle_max_level =
        bfs_spanning_tree(graph, middle_verts, super_x);

    // Step 7: compute descendant costs and triangulate.
    compute_descendant_costs(graph, middle_verts, middle_max_level);
    std::size_t edges_before_triangulation = graph.num_edges();
    triangulate_faces(graph, middle_verts);

    // Steps 8–9: Find a fundamental cycle that separates the middle.
    // Find any non-tree edge among middle vertices.
    std::size_t cycle_edge = NONE;
    for (std::size_t v : middle_verts) {
        graph.for_each_edge_cw(v, [&](std::size_t ei) {
            if (cycle_edge != NONE) return;
            if (graph.edge(ei).is_tree) return;
            std::size_t w = graph.edge(ei).other(v);
            std::size_t lw = graph.vertex(w).level;
            if (lw > l0 && lw < l2) {
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
        // Each iteration finds the triangle on the inside of the
        // current non-tree edge and replaces the edge with one that
        // removes at least one face from the inside.
        //
        // Full implementation of the three sub-cases:
        std::unordered_set<std::size_t> cycle_set(
            cycle_verts.begin(), cycle_verts.end());

        double inside_cost = compute_inside_cost(
            graph, cycle_verts, cycle_set, comp.total_cost);

        std::size_t current_edge = cycle_edge;
        int max_iters = static_cast<int>(comp.vertices.size());

        while (inside_cost > 2.0 / 3.0 * comp.total_cost
               && max_iters-- > 0)
        {
            // Find the triangle on the inside of current_edge.
            std::size_t cu2 = graph.edge(current_edge).endpoint[0];
            std::size_t cv2 = graph.edge(current_edge).endpoint[1];

            // Two faces border this edge.  Try both and pick the one
            // whose third vertex is "inside."
            auto [ne1, nd1] = graph.face_next(current_edge, cv2);
            std::size_t y1 = nd1;
            auto [ne2, nd2] = graph.face_next(current_edge, cu2);
            std::size_t y2 = nd2;

            // Pick y that is inside the cycle (not on the cycle and
            // not outside).
            std::size_t y = NONE;

            if (cycle_set.find(y1) == cycle_set.end()) {
                // y1 might be inside.
                y = y1;
            }
            if (y == NONE && cycle_set.find(y2) == cycle_set.end()) {
                y = y2;
            }
            if (y == NONE) break; // both triangle vertices are on cycle

            // Find the two edges from the triangle: (cu2, y) and (y, cv2)
            // (or (cv2, y) and (y, cu2) depending on orientation).
            // Find a non-tree edge among the triangle's edges.
            std::size_t nontree_edge = NONE;
            graph.for_each_edge_cw(y, [&](std::size_t ei) {
                if (nontree_edge != NONE) return;
                if (graph.edge(ei).is_tree) return;
                std::size_t w = graph.edge(ei).other(y);
                if (w == cu2 || w == cv2) {
                    nontree_edge = ei;
                }
            });

            if (nontree_edge == NONE) break; // shouldn't happen

            // Replace current_edge with nontree_edge.
            current_edge = nontree_edge;

            // Recompute cycle and inside cost.
            std::size_t nu = graph.edge(current_edge).endpoint[0];
            std::size_t nv2 = graph.edge(current_edge).endpoint[1];
            cycle_verts = fundamental_cycle_vertices(graph, nu, nv2);
            cycle_set.clear();
            cycle_set.insert(cycle_verts.begin(), cycle_verts.end());
            inside_cost = compute_inside_cost(
                graph, cycle_verts, cycle_set, comp.total_cost);
        }
    }

    // ── Step 10: assemble partition ─────────────────────────────
    // Clean up edges added by triangulate_faces.
    for (std::size_t ei = edges_before_triangulation; ei < graph.num_edges(); ++ei) {
        graph.edge(ei).deleted = true;
    }

    // Restore vertices deleted in Step 6.
    for (auto& [v, was_deleted] : step6_saved) {
        graph.vertex(v).deleted = was_deleted;
    }

    // C = separator levels (l₀, l₂) + cycle vertices
    std::unordered_set<std::size_t> sep_set;
    for (std::size_t v : sep_levels) sep_set.insert(v);
    for (std::size_t v : cycle_verts) sep_set.insert(v);

    for (std::size_t v : sep_set) {
        result.C.push_back(v);
    }

    // Partition remaining vertices into A (inside cycle) and B (outside).
    // BFS from tree root, not crossing separator.
    std::size_t tree_root = root;
    std::vector<bool> reachable(nv, false);
    {
        std::queue<std::size_t> q;
        if (sep_set.find(tree_root) == sep_set.end()) {
            reachable[tree_root] = true;
            q.push(tree_root);
        }
        while (!q.empty()) {
            std::size_t v = q.front();
            q.pop();
            graph.for_each_edge_cw(v, [&](std::size_t ei) {
                std::size_t w = graph.edge(ei).other(v);
                if (reachable[w] || graph.vertex(w).deleted) return;
                if (sep_set.count(w)) return;
                reachable[w] = true;
                q.push(w);
            });
        }
    }

    for (std::size_t v : comp.vertices) {
        if (sep_set.count(v)) continue;
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

    // Work queue of pieces to (potentially) split.
    std::queue<std::vector<std::size_t>> work;
    work.push(std::move(initial));

    while (!work.empty()) {
        auto piece = std::move(work.front());
        work.pop();

        if (piece.size() <= max_piece_size) {
            // Small enough — record as a final piece.
            std::size_t pidx = result.pieces.size();
            for (std::size_t v : piece) {
                result.vertex_piece[v] = pidx;
            }
            result.pieces.push_back(std::move(piece));
            continue;
        }

        // Build a subgraph for this piece and find a separator.
        // Soft-delete ALL vertices not in this piece so that
        // find_separator only operates on piece vertices.
        std::unordered_set<std::size_t> in_piece(piece.begin(), piece.end());

        // Save and soft-delete every non-piece vertex that's currently active.
        std::vector<std::pair<std::size_t, bool>> saved;
        for (std::size_t v = 0; v < nv; ++v) {
            if (!in_piece.count(v) && !graph.vertex(v).deleted) {
                saved.emplace_back(v, false);
                graph.vertex(v).deleted = true;
            }
        }

        // Reset costs to uniform.
        for (std::size_t v : piece) {
            graph.vertex(v).cost = 1.0;
        }

        auto sep = find_separator(graph);

        // Restore.
        for (auto& [v, _] : saved) {
            graph.vertex(v).deleted = false;
        }

        // Mark separator vertices.
        for (std::size_t v : sep.C) {
            result.separator_nodes[v] = true;
        }

        // Recurse on A and B (exclude separator vertices).
        if (!sep.A.empty()) {
            std::vector<std::size_t> a_piece;
            for (std::size_t v : sep.A) {
                if (!result.separator_nodes[v]) a_piece.push_back(v);
            }
            if (!a_piece.empty()) work.push(std::move(a_piece));
        }
        if (!sep.B.empty()) {
            std::vector<std::size_t> b_piece;
            for (std::size_t v : sep.B) {
                if (!result.separator_nodes[v]) b_piece.push_back(v);
            }
            if (!b_piece.empty()) work.push(std::move(b_piece));
        }
    }

    return result;
}

} // namespace chazelle
