#include "visibility/tree_decomposition.h"
#include "visibility/submap.h"

#include <algorithm>
#include <cassert>
#include <queue>
#include <unordered_map>
#include <vector>

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

    // §2.3: O(m log m + 1) deterministic tree decomposition.
    //
    // Use compact local indexing with vectors instead of hash maps.
    // At each recursion level, we map the subtree's submap node IDs
    // into a dense [0, N) range, so all lookups are O(1) array
    // accesses.  Total work across all levels of each recursion depth
    // is O(|subtree|), giving O(r log r) overall (¾ shrinkage).

    const std::size_t N = subtree_nodes.size();
    const std::size_t total_edges = subtree_chords.size();

    // §2.3: "By using straightforward tree-labeling techniques we can
    // find the centroid node, and from there, the first edge to be
    // removed, in linear time.  Proceeding recursively gives us a
    // simple O(m log m + 1)-time algorithm for computing the tree
    // decomposition."
    //
    // Build a mapping: submap node ID → compact local index [0, N).
    // Use a hash map for O(1) expected-time lookups, giving O(N) per
    // recursion level and O(r log r) total — matching the paper's
    // O(m log m + 1) bound.
    std::unordered_map<std::size_t, std::size_t> id_map;
    id_map.reserve(N);
    for (std::size_t i = 0; i < N; ++i)
        id_map[subtree_nodes[i]] = i;

    auto compact_id = [&](std::size_t sid) -> std::size_t {
        auto it = id_map.find(sid);
        assert(it != id_map.end());
        return it->second;
    };

    // Build local adjacency: adj[local_id] → [(chord_idx, neighbor_local_id)]
    struct AdjEntry { std::size_t chord_idx; std::size_t neighbor; };
    std::vector<std::vector<AdjEntry>> adj(N);
    for (std::size_t ci : subtree_chords) {
        auto& c = submap.chord(ci);
        std::size_t lu = compact_id(c.region[0]);
        std::size_t lv = compact_id(c.region[1]);
        adj[lu].push_back({ci, lv});
        adj[lv].push_back({ci, lu});
    }

    // BFS from local node 0.
    std::vector<std::size_t> parent_local(N, NONE);
    std::vector<std::size_t> parent_chord_local(N, NONE);
    std::vector<std::size_t> sub_size(N, 0);
    std::vector<std::size_t> bfs_order;
    bfs_order.reserve(N);

    {
        std::queue<std::size_t> q;
        q.push(0);
        parent_local[0] = 0; // sentinel: root is its own parent
        while (!q.empty()) {
            std::size_t u = q.front(); q.pop();
            bfs_order.push_back(u);
            for (auto& [ci, v] : adj[u]) {
                if (parent_local[v] == NONE) {
                    parent_local[v] = u;
                    parent_chord_local[v] = ci;
                    q.push(v);
                }
            }
        }
    }

    // Compute subtree edge counts in reverse BFS order (leaves first).
    for (auto it = bfs_order.rbegin(); it != bfs_order.rend(); ++it) {
        std::size_t v = *it;
        sub_size[v] = 0;
        for (auto& [ci, w] : adj[v]) {
            if (parent_local[w] == v) {
                sub_size[v] += sub_size[w] + 1;
            }
        }
    }

    // §2.3: "pick the centroid of the submap's tree and observe that
    // there exists at least one incident edge whose removal leaves two
    // subtrees, each with a number of edges at most three-quarters the
    // original number."
    //
    // Step 1: Find the centroid NODE — the node whose removal minimises
    // the maximum component edge-count.
    std::size_t centroid = 0;
    std::size_t centroid_max_comp = total_edges;

    for (std::size_t v : bfs_order) {
        std::size_t max_comp = 0;
        for (auto& [ci, w] : adj[v]) {
            if (parent_local[w] == v) {
                max_comp = std::max(max_comp, sub_size[w]);
            }
        }
        if (v != 0) {
            max_comp = std::max(max_comp,
                                total_edges - sub_size[v] - 1);
        }
        if (max_comp < centroid_max_comp) {
            centroid_max_comp = max_comp;
            centroid = v;
        }
    }

    // Step 2: Among the centroid's incident edges, pick the one whose
    // removal gives the most balanced split.  The ¾ guarantee holds
    // because the submap is conformal (degree ≤ 4).
    std::size_t best_chord = subtree_chords[0];
    std::size_t best_max = total_edges;
    std::size_t best_child = NONE;

    for (auto& [ci, w] : adj[centroid]) {
        std::size_t child_side_edges;
        std::size_t chord_for_edge;
        std::size_t child_node;

        if (parent_local[w] == centroid) {
            // w is a child of centroid.
            child_side_edges = sub_size[w];
            chord_for_edge = parent_chord_local[w];
            child_node = w;
        } else if (centroid != 0) {
            // w is the parent of centroid.
            child_side_edges = sub_size[centroid];
            chord_for_edge = parent_chord_local[centroid];
            child_node = centroid;
        } else {
            continue;
        }

        std::size_t other_side = total_edges - 1 - child_side_edges;
        std::size_t mx = std::max(child_side_edges, other_side);
        if (mx < best_max) {
            best_max = mx;
            best_chord = chord_for_edge;
            best_child = child_node;
        }
    }

    // Create the internal TD node for this centroid chord.
    TDNode internal;
    internal.chord_idx = best_chord;
    internal.parent = parent_td;
    internal.depth = depth;
    std::size_t td_idx = nodes_.size();
    nodes_.push_back(internal);

    // Split subtree into two halves around best_chord.
    std::size_t side_a_root, side_b_root; // local IDs
    if (best_child != NONE) {
        side_b_root = best_child;
        side_a_root = parent_local[best_child];
    } else {
        auto& cc = submap.chord(best_chord);
        side_a_root = compact_id(cc.region[0]);
        side_b_root = compact_id(cc.region[1]);
    }

    // BFS from each side to collect submap node IDs and chord indices.
    std::vector<bool> visited(N, false);

    auto collect_side = [&](std::size_t root_local) {
        std::vector<std::size_t> side_nodes;
        std::vector<std::size_t> side_chords;
        std::queue<std::size_t> q;
        q.push(root_local);
        visited[root_local] = true;
        while (!q.empty()) {
            std::size_t u = q.front(); q.pop();
            side_nodes.push_back(subtree_nodes[u]); // map back to submap ID
            for (auto& [ci, v] : adj[u]) {
                if (ci == best_chord) continue;
                if (!visited[v]) {
                    visited[v] = true;
                    q.push(v);
                    side_chords.push_back(ci);
                }
            }
        }
        return std::make_pair(side_nodes, side_chords);
    };

    auto [nodes_a, chords_a] = collect_side(side_a_root);
    auto [nodes_b, chords_b] = collect_side(side_b_root);

    // Recurse on both sides.
    std::size_t left  = decompose(submap, nodes_a, chords_a, td_idx, depth + 1);
    std::size_t right = decompose(submap, nodes_b, chords_b, td_idx, depth + 1);

    nodes_[td_idx].left_child  = left;
    nodes_[td_idx].right_child = right;

    return td_idx;
}

} // namespace chazelle
