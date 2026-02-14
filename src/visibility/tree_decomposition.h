#pragma once

/// TreeDecomposition — hierarchical centroid decomposition of a conformal
/// submap tree.
///
/// Per Chazelle §2.3:
/// - Internal nodes correspond to exit chords (centroid edges).
/// - Leaves correspond to regions of the submap.
/// - Depth is O(log r) where r = number of regions.
/// - Each level splits the tree via a centroid edge whose removal yields
///   two subtrees, each with ≤ ¾ of the original edge count.
///
/// Used by conformality restoration (§3.2) for binary search over tree
/// levels: at each level, test the centroid chord using ray-shooting
/// and the topological lemma (Lemma 2.4) to narrow down the search.

#include <cstddef>
#include <limits>
#include <vector>

namespace chazelle {

/// Sentinel for tree decomposition indices.  Same underlying value as NONE
/// but used in TD-specific code for readability.
inline constexpr std::size_t TD_NONE = std::numeric_limits<std::size_t>::max();
static_assert(TD_NONE == std::numeric_limits<std::size_t>::max(),
              "TD_NONE must be the max sentinel value");

/// A node in the tree decomposition hierarchy.
struct TDNode {
    /// If this is an internal node: the chord index in the submap that
    /// was used as the centroid edge for the split.
    std::size_t chord_idx = TD_NONE;

    /// If this is a leaf: the region (submap node) index.
    std::size_t region_idx = TD_NONE;

    bool is_leaf() const { return chord_idx == TD_NONE; }
    bool is_internal() const { return chord_idx != TD_NONE; }

    /// Parent in the decomposition tree (TD_NONE for root).
    std::size_t parent = TD_NONE;

    /// Children (at most 2 for binary decomposition).
    std::size_t left_child  = TD_NONE;
    std::size_t right_child = TD_NONE;

    /// Depth in the decomposition tree.
    std::size_t depth = 0;
};

class Submap; // forward declaration

class TreeDecomposition {
public:
    TreeDecomposition() = default;

    /// Build the tree decomposition from a conformal submap.
    /// Time: O(r log r) where r = number of regions (deterministic).
    void build(const Submap& submap);

    /// Root of the decomposition tree.
    std::size_t root() const { return root_; }

    /// Access a node.
    const TDNode& node(std::size_t i) const { return nodes_[i]; }

    /// Total number of nodes (internal + leaves).
    std::size_t size() const { return nodes_.size(); }

    /// Maximum depth.
    std::size_t max_depth() const { return max_depth_; }

private:
    std::vector<TDNode> nodes_;
    std::size_t root_ = TD_NONE;
    std::size_t max_depth_ = 0;

    /// Recursive centroid decomposition.
    /// @param submap         The submap being decomposed.
    /// @param subtree_nodes  The set of submap-tree node indices in this subtree.
    /// @param subtree_chords The set of chord indices forming the subtree.
    /// @param parent_td      Parent index in the decomposition tree.
    /// @param depth          Current depth.
    /// @param node_id_map    Dense vector (size ≥ num_nodes) mapping
    ///                       submap node ID → compact local index.
    ///                       Shared across all recursion levels for O(1)
    ///                       lookup per node (§2.3: O(m log m+1) total).
    /// @return Index of the created TD node.
    std::size_t decompose(
        const Submap& submap,
        const std::vector<std::size_t>& subtree_nodes,
        const std::vector<std::size_t>& subtree_chords,
        std::size_t parent_td,
        std::size_t depth,
        std::vector<std::size_t>& node_id_map);
};

} // namespace chazelle
