#pragma once

// [C91 §2.3]: Hierarchical centroid decomposition of a conformal submap's
// tree.  "Pick the centroid; some incident edge's removal leaves two
// subtrees each with ≤ 3/4 the edges.  Associate the removed edge with
// the root, recurse on its children."
//
// Internal nodes ↔ exit chords; leaves ↔ regions.  Depth O(log m),
// built in O(m log m).

#include "../common.h"

#include <cstddef>
#include <vector>

namespace chazelle {

struct TDNode {
    // Internal node: index of the chord in the Submap.  Leaf: NONE.
    std::size_t chord_idx = NONE;
    // Leaf: index of the region in the Submap.  Internal: NONE.
    std::size_t region_idx = NONE;

    std::size_t parent      = NONE;
    std::size_t left_child  = NONE;
    std::size_t right_child = NONE;

    bool is_leaf()     const noexcept { return chord_idx == NONE; }
    bool is_internal() const noexcept { return chord_idx != NONE; }
};

class TreeDecomposition {
public:
    // Submap forward-declared; defined in submap.h.
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

    // [§2.3 tex 116]: O(m log m) via a shared region_local_buf
    // (pre-allocated to submap.num_nodes(); cleaned per call).
    std::size_t decompose(
        const class Submap& submap,
        const std::vector<std::size_t>& chords_in_subtree,
        const std::vector<std::size_t>& regions_in_subtree,
        std::size_t parent_td,
        std::vector<std::size_t>& region_local_buf);
};

} // namespace chazelle
