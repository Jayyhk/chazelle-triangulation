/// src/submap/tree_decomposition.cpp

#include "tree_decomposition.h"
#include "submap.h"

#include <algorithm>
#include <cassert>

namespace chazelle {

void TreeDecomposition::build(const Submap& submap) {
    nodes_.clear();
    root_ = NONE;

    // [C91 §2.4(iv)]: "If the submap is conformal, then its tree
    // decomposition should be available."
    assert(submap.is_conformal() &&
           "§2.4(iv): tree decomposition requires conformal submap");

    // [C91 §2.3] (tex 114): tree decomposition requires ≥ 1 region.
    assert(submap.num_nodes() >= 1 &&
           "§2.3: submap must have at least one region");

    // Precondition: submap must be compacted (no dead entries).
    // build() indexes into nodes_/chords_ directly — dead entries
    // would produce a wrong decomposition.
    assert(submap.num_live_nodes() == submap.num_nodes() &&
           submap.num_live_chords() == submap.num_chords() &&
           "§2.4(iv): tree decomposition requires compacted submap");

    // Collect all chords and regions.
    std::vector<std::size_t> all_chords;
    for (std::size_t i = 0; i < submap.num_chords(); ++i) {
        all_chords.push_back(i);
    }
    std::vector<std::size_t> all_regions;
    for (std::size_t i = 0; i < submap.num_nodes(); ++i) {
        all_regions.push_back(i);
    }

    if (all_chords.empty()) {
        // No chords → single leaf (one region).
        // The assert above guarantees num_nodes() >= 1, so
        // all_regions is never empty here.
        assert(!all_regions.empty() &&
               "§2.3: no-chord submap must have exactly one region");
        TDNode leaf;
        leaf.region_idx = all_regions[0];
        root_ = nodes_.size();
        nodes_.push_back(leaf);
        return;
    }

    // [C91 §2.3] (tex 116): Allocate reusable buffer once — O(m).
    std::vector<std::size_t> region_local_buf(submap.num_nodes(), NONE);
    root_ = decompose(submap, all_chords, all_regions, NONE, region_local_buf);

    // [C91 §2.3]: "The internal nodes (resp. leaves) are in bijection
    // with the exit chords (resp. regions) of the submap."
    std::size_t num_internal = 0, num_leaves = 0;
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i].is_internal()) ++num_internal;
        else ++num_leaves;
    }
    assert(num_internal == all_chords.size() &&
           "§2.3: internal TD nodes must biject with chords");
    assert(num_leaves == all_regions.size() &&
           "§2.3: TD leaves must biject with regions");
}

std::size_t TreeDecomposition::decompose(
        const Submap& submap,
        const std::vector<std::size_t>& chords,
        const std::vector<std::size_t>& regions,
        std::size_t parent_td,
        std::vector<std::size_t>& region_local) {

    if (chords.empty()) {
        // Leaf: single region.
        assert(regions.size() == 1 &&
               "§2.3: leaf of tree decomposition must be a single region");
        TDNode leaf;
        leaf.region_idx = regions[0];
        leaf.parent = parent_td;
        std::size_t idx = nodes_.size();
        nodes_.push_back(leaf);
        return idx;
    }

    // [C91 §2.3]: "straightforward tree-labeling techniques."
    // Build local adjacency using the shared region_local buffer.
    // Write O(n_local) entries, clean up after — O(n_local) per call,
    // O(m) per level, O(m log m) total.

    // Map global region indices to local [0, regions.size()).
    for (std::size_t i = 0; i < regions.size(); ++i)
        region_local[regions[i]] = i;

    std::size_t n = regions.size();
    std::size_t m = chords.size();
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> adj(n);
    // adj[local_region] = list of (chord_index_in_chords_vec, local_other_region)

    for (std::size_t ci = 0; ci < m; ++ci) {
        const auto& c = submap.chord(chords[ci]);
        std::size_t lr0 = region_local[c.region[0]];
        std::size_t lr1 = region_local[c.region[1]];
        adj[lr0].push_back({ci, lr1});
        adj[lr1].push_back({ci, lr0});
    }

    // Find centroid: root an arbitrary tree, compute subtree sizes,
    // pick the node where max subtree ≤ n/2 (which guarantees ≤ 3/4
    // of edges for bounded-degree trees).
    std::vector<std::size_t> subtree_size(n, 1);
    std::vector<std::size_t> parent_local(n, NONE);
    std::vector<std::size_t> order;
    order.reserve(n);

    // BFS to establish parent pointers and processing order.
    std::vector<bool> visited(n, false);
    order.push_back(0);
    visited[0] = true;
    for (std::size_t qi = 0; qi < order.size(); ++qi) {
        std::size_t u = order[qi];
        for (auto [ci, v] : adj[u]) {
            if (!visited[v]) {
                visited[v] = true;
                parent_local[v] = u;
                order.push_back(v);
            }
        }
    }

    // Compute subtree sizes bottom-up.
    for (std::size_t i = order.size(); i-- > 0; ) {
        std::size_t u = order[i];
        if (parent_local[u] != NONE)
            subtree_size[parent_local[u]] += subtree_size[u];
    }

    // Find centroid: node where max(child subtree, n - subtree) ≤ n/2.
    std::size_t centroid = 0;
    std::size_t best_max = n;
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t max_sub = n - subtree_size[i]; // parent side
        for (auto [ci, v] : adj[i]) {
            if (parent_local[v] == i) // v is a child
                max_sub = std::max(max_sub, subtree_size[v]);
        }
        if (max_sub < best_max) {
            best_max = max_sub;
            centroid = i;
        }
    }

    // [C91 §2.3]: "there exists at least one incident edge whose
    // removal leaves two subtrees, each with at most three-quarters
    // the original number."
    // Pick the edge that gives the best split.
    std::size_t best_chord_local = 0;
    std::size_t best_split = n;
    for (auto [ci, v] : adj[centroid]) {
        // Removing this edge splits into: subtree containing v, and the rest.
        std::size_t v_size;
        if (parent_local[v] == centroid) {
            v_size = subtree_size[v];
        } else {
            v_size = n - subtree_size[centroid];
        }
        std::size_t split = std::max(v_size, n - v_size);
        if (split < best_split) {
            best_split = split;
            best_chord_local = ci;
        }
    }

    // [C91 §2.3] (tex 114): "each with a number of edges at most
    // three-quarters the original number."  Paper bound is on edges
    // (n−1 total), so ≤ 3(n−1)/4 edges ⟹ ≤ ⌈3(n−1)/4⌉ + 1 regions.
    assert(best_split <= (3 * (n - 1) + 3) / 4 &&
           "§2.3: centroid split must give ≤ 3/4 edges on each side");

    // The chosen chord becomes an internal node of the decomposition.
    std::size_t chosen_chord = chords[best_chord_local];
    TDNode internal;
    internal.chord_idx = chosen_chord;
    internal.parent = parent_td;
    std::size_t td_idx = nodes_.size();
    nodes_.push_back(internal);

    // Split regions and chords into two sides.
    const auto& cc = submap.chord(chosen_chord);
    std::size_t side0 = region_local[cc.region[0]];
    [[maybe_unused]] std::size_t side1 = region_local[cc.region[1]];

    // BFS from side0, not crossing chosen chord.
    std::vector<bool> in_left(n, false);
    {
        std::vector<std::size_t> q;
        q.push_back(side0);
        in_left[side0] = true;
        for (std::size_t qi = 0; qi < q.size(); ++qi) {
            std::size_t u = q[qi];
            for (auto [ci, v] : adj[u]) {
                if (ci != best_chord_local && !in_left[v]) {
                    in_left[v] = true;
                    q.push_back(v);
                }
            }
        }
    }

    std::vector<std::size_t> left_regions, right_regions;
    std::vector<std::size_t> left_chords, right_chords;

    for (std::size_t i = 0; i < n; ++i) {
        if (in_left[i])
            left_regions.push_back(regions[i]);
        else
            right_regions.push_back(regions[i]);
    }
    for (std::size_t ci = 0; ci < m; ++ci) {
        if (ci == best_chord_local) continue;
        const auto& c = submap.chord(chords[ci]);
        std::size_t lr0 = region_local[c.region[0]];
        if (in_left[lr0])
            left_chords.push_back(chords[ci]);
        else
            right_chords.push_back(chords[ci]);
    }

    // Recurse on both sides.
    std::size_t left_td = decompose(submap, left_chords, left_regions, td_idx, region_local);
    std::size_t right_td = decompose(submap, right_chords, right_regions, td_idx, region_local);

    nodes_[td_idx].left_child = left_td;
    nodes_[td_idx].right_child = right_td;

    // No explicit cleanup needed: each recursive call overwrites
    // region_local entries for its own regions before reading them,
    // so stale values from sibling/parent calls are never observed.

    return td_idx;
}

} // namespace chazelle
