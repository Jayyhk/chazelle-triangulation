// src/submap/tree_decomposition.cpp

#include "tree_decomposition.h"
#include "submap.h"

#include <algorithm>
#include <cassert>

namespace chazelle {

void TreeDecomposition::build(const Submap& submap) {
    nodes_.clear();
    root_ = NONE;

    assert(submap.is_conformal() &&
           "[C91 §2.4(iv)]: tree decomposition requires conformal submap");
    assert(submap.num_nodes() >= 1 &&
           "[C91 §2.3]: submap must have ≥ 1 region");
    // Compacted — build() indexes directly into nodes_/chords_, so dead
    // entries would corrupt the decomposition.
    assert(submap.num_live_nodes() == submap.num_nodes() &&
           submap.num_live_chords() == submap.num_chords() &&
           "[C91 §2.4(iv)]: tree decomposition requires compacted submap");

    std::vector<std::size_t> all_chords;
    for (std::size_t i = 0; i < submap.num_chords(); ++i)
        all_chords.push_back(i);
    std::vector<std::size_t> all_regions;
    for (std::size_t i = 0; i < submap.num_nodes(); ++i)
        all_regions.push_back(i);

    if (all_chords.empty()) {
        // Single-region submap → one leaf.
        assert(!all_regions.empty() && "[C91 §2.3]: no-chord submap has 1 region");
        TDNode leaf;
        leaf.region_idx = all_regions[0];
        root_ = nodes_.size();
        nodes_.push_back(leaf);
        return;
    }

    // [C91 §2.3 tex 116]: reusable buffer, O(m).
    std::vector<std::size_t> region_local_buf(submap.num_nodes(), NONE);
    root_ = decompose(submap, all_chords, all_regions, NONE, region_local_buf);

    // [C91 §2.3]: internal nodes ↔ chords; leaves ↔ regions (bijection).
    std::size_t num_internal = 0, num_leaves = 0;
    for (const auto& n : nodes_) {
        if (n.is_internal()) ++num_internal; else ++num_leaves;
    }
    assert(num_internal == all_chords.size() &&
           "[C91 §2.3]: tree decomposition internals must biject with chords");
    assert(num_leaves == all_regions.size() &&
           "[C91 §2.3]: tree decomposition leaves must biject with regions");
}

std::size_t TreeDecomposition::decompose(
        const Submap& submap,
        const std::vector<std::size_t>& chords,
        const std::vector<std::size_t>& regions,
        std::size_t parent_idx,
        std::vector<std::size_t>& region_local) {

    if (chords.empty()) {
        assert(regions.size() == 1 && "[C91 §2.3]: leaf is a single region");
        TDNode leaf;
        leaf.region_idx = regions[0];
        leaf.parent = parent_idx;
        std::size_t idx = nodes_.size();
        nodes_.push_back(leaf);
        return idx;
    }

    // Map global region indices → local [0, regions.size()) via the shared
    // buffer.  Per call: O(n_local) writes, O(m log m) total.
    for (std::size_t i = 0; i < regions.size(); ++i)
        region_local[regions[i]] = i;

    std::size_t n = regions.size();
    std::size_t m = chords.size();
    // adj[local_region] = (chord_idx_in_chords_vec, local_other_region) list.
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> adj(n);

    for (std::size_t ci = 0; ci < m; ++ci) {
        const auto& c = submap.chord(chords[ci]);
        std::size_t lr0 = region_local[c.region[0]];
        std::size_t lr1 = region_local[c.region[1]];
        adj[lr0].push_back({ci, lr1});
        adj[lr1].push_back({ci, lr0});
    }

    // Centroid: BFS to root the tree, get subtree sizes bottom-up, pick
    // the node minimizing max(child subtree, n − subtree).
    std::vector<std::size_t> subtree_size(n, 1);
    std::vector<std::size_t> parent_local(n, NONE);
    std::vector<std::size_t> order;
    order.reserve(n);
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
    for (std::size_t i = order.size(); i-- > 0; ) {
        std::size_t u = order[i];
        if (parent_local[u] != NONE)
            subtree_size[parent_local[u]] += subtree_size[u];
    }

    std::size_t centroid = 0;
    std::size_t best_max = n;
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t max_sub = n - subtree_size[i]; // parent side
        for (auto [ci, v] : adj[i])
            if (parent_local[v] == i)
                max_sub = std::max(max_sub, subtree_size[v]);
        if (max_sub < best_max) { best_max = max_sub; centroid = i; }
    }

    // [C91 §2.3]: some edge incident on the centroid splits ≤ 3/4 / 3/4.
    std::size_t best_chord_local = 0;
    std::size_t best_split = n;
    for (auto [ci, v] : adj[centroid]) {
        // Removing this edge splits into: subtree containing v + the rest.
        std::size_t v_size = (parent_local[v] == centroid)
            ? subtree_size[v]
            : n - subtree_size[centroid];
        std::size_t split = std::max(v_size, n - v_size);
        if (split < best_split) { best_split = split; best_chord_local = ci; }
    }

    // [C91 §2.3 tex 114]: paper bound is on edges (n−1 total): ≤ ⌊3(n−1)/4⌋
    // edges per side ⟹ ≤ ⌊3(n−1)/4⌋ + 1 regions.
    assert(best_split <= 3 * (n - 1) / 4 + 1 &&
           "[C91 §2.3]: centroid split must give ≤ 3/4 edges per side");

    std::size_t chosen_chord = chords[best_chord_local];
    TDNode internal;
    internal.chord_idx = chosen_chord;
    internal.parent = parent_idx;
    std::size_t new_idx = nodes_.size();
    nodes_.push_back(internal);

    // BFS from one chord side, not crossing the chosen chord.
    const auto& cc = submap.chord(chosen_chord);
    std::size_t side0 = region_local[cc.region[0]];
    [[maybe_unused]] std::size_t side1 = region_local[cc.region[1]];
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

    std::size_t left_idx = decompose(submap, left_chords, left_regions, new_idx, region_local);
    std::size_t right_idx = decompose(submap, right_chords, right_regions, new_idx, region_local);
    nodes_[new_idx].left_child  = left_idx;
    nodes_[new_idx].right_child = right_idx;

    // No cleanup of region_local: each recursive call rewrites entries
    // for its own regions before reading, so stale values are never seen.
    return new_idx;
}

} // namespace chazelle
