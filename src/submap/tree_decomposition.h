#pragma once

/// [C91 §2.3]: Hierarchical tree decomposition of a conformal submap.
///
/// "The idea is to pick the centroid of the submap's tree and observe
/// that there exists at least one incident edge whose removal leaves
/// two subtrees, each with a number of edges at most three-quarters
/// the original number.  Associating the removed edge with the root
/// of a binary tree and recursing in this fashion with respect to
/// the root's two children provides a tree decomposition."
///
/// Internal nodes ↔ exit chords.  Leaves ↔ regions.
/// Depth = O(log m).  Built in O(m log m) time.

#include "../common.h"

#include <cstddef>
#include <vector>

namespace chazelle {

struct TDNode {
    /// For internal nodes: index of the chord in the Submap.
    std::size_t chord_idx = NONE;

    /// For leaf nodes: index of the region in the Submap.
    std::size_t region_idx = NONE;

    std::size_t parent      = NONE;
    std::size_t left_child  = NONE;
    std::size_t right_child = NONE;

    bool is_leaf()     const noexcept { return chord_idx == NONE; }
    bool is_internal() const noexcept { return chord_idx != NONE; }
};

class TreeDecomposition {
public:
    /// Build the decomposition from a conformal submap.
    /// Requires forward declaration — Submap is defined in submap.h.
    void build(const class Submap& submap);

    std::size_t root() const noexcept { return root_; }

    const TDNode& node(std::size_t i) const noexcept {
        return nodes_[i];
    }

    std::size_t size() const noexcept { return nodes_.size(); }
    bool empty()       const noexcept { return nodes_.empty(); }

private:
    std::vector<TDNode> nodes_;
    std::size_t root_ = NONE;

    /// Recursive centroid decomposition on a subtree.
    /// [C91 §2.3] (tex 116): O(m log m) via reusable buffer.
    /// region_local_buf: pre-allocated buffer of size submap.num_nodes(),
    ///   shared across all recursive calls.  Used entries are cleaned
    ///   up after each call (O(n_local) per call, O(m) per level).
    std::size_t decompose(
        const class Submap& submap,
        const std::vector<std::size_t>& chords_in_subtree,
        const std::vector<std::size_t>& regions_in_subtree,
        std::size_t parent_td,
        std::vector<std::size_t>& region_local_buf);
};

} // namespace chazelle
