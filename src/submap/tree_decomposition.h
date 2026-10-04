#pragma once

#include "../common.h"

#include <cstddef>
#include <vector>

namespace chazelle {

struct TreeDecompositionNode {
    std::size_t chord_idx = NONE;

    std::size_t region_idx = NONE;

    std::size_t parent = NONE;
    std::size_t left_child = NONE;
    std::size_t right_child = NONE;

    bool is_leaf() const noexcept {
        return chord_idx == NONE;
    }
    bool is_internal() const noexcept {
        return chord_idx != NONE;
    }
};

class TreeDecomposition {
public:
    void build(const class Submap& submap);

    std::size_t root() const noexcept {
        return root_;
    }

    const TreeDecompositionNode& node(std::size_t i) const noexcept {
        return nodes_[i];
    }

    std::size_t size() const noexcept {
        return nodes_.size();
    }
    bool empty() const noexcept {
        return nodes_.empty();
    }

private:
    std::vector<TreeDecompositionNode> nodes_;
    std::size_t root_ = NONE;
};

}
