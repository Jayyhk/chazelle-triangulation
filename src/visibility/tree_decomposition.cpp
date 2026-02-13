#include "visibility/tree_decomposition.h"
#include "visibility/submap.h"

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <queue>

namespace chazelle {

void TreeDecomposition::build(const Submap& submap) {
    nodes_.clear();
    root_ = TD_NONE;
    max_depth_ = 0;

    // Collect all live nodes and chords.
    std::vector<std::size_t> all_nodes;
    for (std::size_t i = 0; i < submap.num_nodes(); ++i) {
        if (!submap.node(i).deleted) all_nodes.push_back(i);
    }

    std::vector<std::size_t> all_chords;
    for (std::size_t i = 0; i < submap.num_chords(); ++i) {
        auto& c = submap.chord(i);
        if (c.region[0] != NONE && c.region[1] != NONE &&
            !submap.node(c.region[0]).deleted &&
            !submap.node(c.region[1]).deleted) {
            all_chords.push_back(i);
        }
    }

    if (all_nodes.empty()) return;

    // If only one region, create a single leaf.
    if (all_chords.empty()) {
        TDNode leaf;
        leaf.region_idx = all_nodes[0];
        leaf.depth = 0;
        nodes_.push_back(leaf);
        root_ = 0;
        return;
    }

    root_ = decompose(submap, all_nodes, all_chords, TD_NONE, 0);
}

std::size_t TreeDecomposition::decompose(
    const Submap& submap,
    const std::vector<std::size_t>& subtree_nodes,
    const std::vector<std::size_t>& subtree_chords,
    std::size_t parent_td,
    std::size_t depth) {

    max_depth_ = std::max(max_depth_, depth);

    // Base case: no chords → single region → leaf.
    if (subtree_chords.empty()) {
        assert(subtree_nodes.size() == 1);
        TDNode leaf;
        leaf.region_idx = subtree_nodes[0];
        leaf.parent = parent_td;
        leaf.depth = depth;
        std::size_t idx = nodes_.size();
        nodes_.push_back(leaf);
        return idx;
    }

    // Find the centroid EDGE of the subtree (§2.3).
    // This is an edge (chord) whose removal splits the tree into two
    // subtrees each with ≤ ¾ of the total edges.  Such an edge always
    // exists in any tree.

    // Build adjacency for the subtree using hash maps keyed by actual
    // subtree node IDs.  This ensures O(|subtree|) allocation per
    // recursion level instead of O(submap.num_nodes()), giving
    // O(r log r) total across all levels (geometric series).
    std::unordered_set<std::size_t> node_set(subtree_nodes.begin(),
                                              subtree_nodes.end());

    // Find subtree sizes by rooting at an arbitrary node.
    std::size_t arb_root = subtree_nodes[0];
    std::size_t total_edges = subtree_chords.size();

    // Build local adjacency map: node → [(chord_idx, neighbor_node)]
    // Keyed only by subtree nodes — O(|subtree|) space.
    struct AdjEntry { std::size_t chord_idx; std::size_t neighbor; };
    std::unordered_map<std::size_t, std::vector<AdjEntry>> adj;
    adj.reserve(subtree_nodes.size());
    for (std::size_t ci : subtree_chords) {
        auto& c = submap.chord(ci);
        std::size_t u = c.region[0], v = c.region[1];
        adj[u].push_back({ci, v});
        adj[v].push_back({ci, u});
    }

    // BFS from arb_root.  All maps keyed by subtree node IDs only.
    std::unordered_map<std::size_t, std::size_t> parent_node;
    std::unordered_map<std::size_t, std::size_t> parent_chord;
    std::unordered_map<std::size_t, std::size_t> subtree_size;
    parent_node.reserve(subtree_nodes.size());
    parent_chord.reserve(subtree_nodes.size());
    subtree_size.reserve(subtree_nodes.size());

    std::vector<std::size_t> bfs_order;
    bfs_order.reserve(subtree_nodes.size());

    {
        std::queue<std::size_t> q;
        q.push(arb_root);
        parent_node[arb_root] = arb_root; // sentinel
        while (!q.empty()) {
            std::size_t u = q.front(); q.pop();
            bfs_order.push_back(u);
            for (auto& [ci, v] : adj[u]) {
                if (!parent_node.count(v) && node_set.count(v)) {
                    parent_node[v] = u;
                    parent_chord[v] = ci;
                    q.push(v);
                }
            }
        }
    }

    // Compute subtree sizes (count of edges in subtree rooted at v).
    // Process in reverse BFS order (leaves first).
    for (auto it = bfs_order.rbegin(); it != bfs_order.rend(); ++it) {
        std::size_t v = *it;
        subtree_size[v] = 0;
        for (auto& [ci, w] : adj[v]) {
            auto pit = parent_node.find(w);
            if (pit != parent_node.end() && pit->second == v) {
                subtree_size[v] += subtree_size[w] + 1;
            }
        }
    }

    // Find centroid edge: edge (parent→child) where both sides ≤ ¾.
    std::size_t best_chord = subtree_chords[0];
    std::size_t best_max = total_edges; // worst case
    std::size_t best_child = NONE;

    for (auto& v : bfs_order) {
        if (v == arb_root) continue;
        std::size_t child_side = subtree_size[v];
        std::size_t other_side = total_edges - 1 - child_side;
        std::size_t mx = std::max(child_side, other_side);
        if (mx < best_max) {
            best_max = mx;
            best_chord = parent_chord[v];
            best_child = v;
        }
    }

    // Create the internal TD node for this centroid chord.
    TDNode internal;
    internal.chord_idx = best_chord;
    internal.parent = parent_td;
    internal.depth = depth;
    std::size_t td_idx = nodes_.size();
    nodes_.push_back(internal);

    // Split subtree into two halves.
    auto& cc = submap.chord(best_chord);
    std::size_t side_a_root, side_b_root;
    if (best_child != NONE) {
        side_b_root = best_child;
        side_a_root = parent_node[best_child];
    } else {
        side_a_root = cc.region[0];
        side_b_root = cc.region[1];
    }

    // BFS from each side to collect nodes and chords.
    auto collect_side = [&](std::size_t root_node) {
        std::vector<std::size_t> side_nodes;
        std::vector<std::size_t> side_chords;
        std::unordered_set<std::size_t> visited;
        std::queue<std::size_t> q;
        q.push(root_node);
        visited.insert(root_node);
        while (!q.empty()) {
            std::size_t u = q.front(); q.pop();
            side_nodes.push_back(u);
            for (auto& [ci, v] : adj[u]) {
                if (ci == best_chord) continue;
                if (!visited.count(v) && node_set.count(v)) {
                    visited.insert(v);
                    q.push(v);
                    side_chords.push_back(ci);
                }
            }
        }
        return std::make_pair(side_nodes, side_chords);
    };

    auto [nodes_a, chords_a] = collect_side(side_a_root);
    auto [nodes_b, chords_b] = collect_side(side_b_root);

    // No cleanup needed — local hash maps go out of scope.

    // Recurse on both sides.
    std::size_t left  = decompose(submap, nodes_a, chords_a, td_idx, depth + 1);
    std::size_t right = decompose(submap, nodes_b, chords_b, td_idx, depth + 1);

    nodes_[td_idx].left_child  = left;
    nodes_[td_idx].right_child = right;

    return td_idx;
}

} // namespace chazelle
