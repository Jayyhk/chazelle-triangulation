#include "tree_decomposition.h"
#include "submap.h"

#include <algorithm>
#include <cassert>
#include <utility>

namespace chazelle {

namespace {

struct ComponentVertex {
    std::size_t region = NONE;
    std::size_t begin = 0;
    std::size_t end = 0;
    std::size_t maximum_end = 0;
    std::size_t size = 1;
    std::size_t height = 1;
    std::size_t parent = NONE;
    std::size_t left = NONE;
    std::size_t right = NONE;
};

class ComponentOrder {
public:
    explicit ComponentOrder(std::vector<ComponentVertex> vertices)
        : vertices_(std::move(vertices)) {}

    std::size_t size(std::size_t root) const {
        return root == NONE ? 0 : vertices_[root].size;
    }

    std::size_t balanced(std::size_t first, std::size_t last) {
        if (first == last)
            return NONE;
        const std::size_t middle = first + (last - first) / 2;
        attach_left(middle, balanced(first, middle));
        attach_right(middle, balanced(middle + 1, last));
        update(middle);
        return middle;
    }

    std::size_t region(std::size_t vertex) const {
        return vertices_[vertex].region;
    }

    std::size_t select(std::size_t root, std::size_t position) const {
        assert(position < size(root));
        while (true) {
            const std::size_t left_size = size(vertices_[root].left);
            if (position == left_size)
                return root;
            if (position < left_size) {
                root = vertices_[root].left;
            } else {
                position -= left_size + 1;
                root = vertices_[root].right;
            }
        }
    }

    std::size_t rank(std::size_t vertex) const {
        std::size_t result = size(vertices_[vertex].left);
        while (vertices_[vertex].parent != NONE) {
            const std::size_t parent = vertices_[vertex].parent;
            if (vertices_[parent].right == vertex)
                result += size(vertices_[parent].left) + 1;
            vertex = parent;
        }
        return result;
    }

    std::size_t subtree_size(std::size_t root, std::size_t vertex) const {
        return lower_bound(root, vertices_[vertex].end) - rank(vertex);
    }

    std::size_t centroid(std::size_t root) const {
        const std::size_t half = size(root) / 2;
        const std::size_t median_key = vertices_[select(root, half)].begin;
        std::size_t first = 0;
        std::size_t last = half + 1;
        while (first < last) {
            const std::size_t middle = first + (last - first) / 2;
            const std::size_t ancestor = ancestor_before(root, middle + 1, median_key);
            assert(ancestor != NONE && "[C91 §2.3 tex 114]: a component contains its root");
            if (subtree_size(root, ancestor) > half)
                first = middle + 1;
            else
                last = middle;
        }
        assert(first > 0 && "[C91 §2.3 tex 114]: the root subtree contains the component");
        return ancestor_before(root, first, median_key);
    }

    std::pair<std::size_t, std::size_t> split(std::size_t root, std::size_t count) {
        assert(count <= size(root));
        if (root == NONE)
            return {NONE, NONE};
        const std::size_t left = detach_left(root);
        const std::size_t right = detach_right(root);
        vertices_[root].parent = NONE;
        const std::size_t left_size = size(left);
        if (count <= left_size) {
            const auto [before, after] = split(left, count);
            return {before, join_with_vertex(after, root, right)};
        }
        const auto [before, after] = split(right, count - left_size - 1);
        return {join_with_vertex(left, root, before), after};
    }

    std::size_t join(std::size_t left, std::size_t right) {
        if (left == NONE)
            return right;
        if (right == NONE)
            return left;
        const auto [before, last] = split(left, size(left) - 1);
        assert(size(last) == 1);
        return join_with_vertex(before, last, right);
    }

private:
    std::vector<ComponentVertex> vertices_;

    std::size_t height(std::size_t root) const {
        return root == NONE ? 0 : vertices_[root].height;
    }

    void attach_left(std::size_t root, std::size_t child) {
        vertices_[root].left = child;
        if (child != NONE)
            vertices_[child].parent = root;
    }

    void attach_right(std::size_t root, std::size_t child) {
        vertices_[root].right = child;
        if (child != NONE)
            vertices_[child].parent = root;
    }

    std::size_t detach_left(std::size_t root) {
        const std::size_t child = vertices_[root].left;
        vertices_[root].left = NONE;
        if (child != NONE)
            vertices_[child].parent = NONE;
        return child;
    }

    std::size_t detach_right(std::size_t root) {
        const std::size_t child = vertices_[root].right;
        vertices_[root].right = NONE;
        if (child != NONE)
            vertices_[child].parent = NONE;
        return child;
    }

    void update(std::size_t root) {
        ComponentVertex& vertex = vertices_[root];
        vertex.size = size(vertex.left) + size(vertex.right) + 1;
        vertex.height = std::max(height(vertex.left), height(vertex.right)) + 1;
        vertex.maximum_end = vertex.end;
        if (vertex.left != NONE)
            vertex.maximum_end = std::max(vertex.maximum_end, vertices_[vertex.left].maximum_end);
        if (vertex.right != NONE)
            vertex.maximum_end = std::max(vertex.maximum_end, vertices_[vertex.right].maximum_end);
    }

    std::size_t rotate_left(std::size_t root) {
        const std::size_t right = detach_right(root);
        attach_right(root, detach_left(right));
        attach_left(right, root);
        update(root);
        update(right);
        return right;
    }

    std::size_t rotate_right(std::size_t root) {
        const std::size_t left = detach_left(root);
        attach_left(root, detach_right(left));
        attach_right(left, root);
        update(root);
        update(left);
        return left;
    }

    std::size_t rebalance(std::size_t root) {
        update(root);
        if (height(vertices_[root].left) > height(vertices_[root].right) + 1) {
            const std::size_t left = vertices_[root].left;
            if (height(vertices_[left].right) > height(vertices_[left].left))
                attach_left(root, rotate_left(detach_left(root)));
            root = rotate_right(root);
        } else if (height(vertices_[root].right) > height(vertices_[root].left) + 1) {
            const std::size_t right = vertices_[root].right;
            if (height(vertices_[right].left) > height(vertices_[right].right))
                attach_right(root, rotate_right(detach_right(root)));
            root = rotate_left(root);
        }
        assert(height(vertices_[root].left) <= height(vertices_[root].right) + 1 &&
               height(vertices_[root].right) <= height(vertices_[root].left) + 1);
        vertices_[root].parent = NONE;
        return root;
    }

    std::size_t join_with_vertex(std::size_t left, std::size_t vertex, std::size_t right) {
        if (height(left) > height(right) + 1) {
            attach_right(left, join_with_vertex(detach_right(left), vertex, right));
            return rebalance(left);
        }
        if (height(right) > height(left) + 1) {
            attach_left(right, join_with_vertex(left, vertex, detach_left(right)));
            return rebalance(right);
        }
        attach_left(vertex, left);
        attach_right(vertex, right);
        update(vertex);
        vertices_[vertex].parent = NONE;
        return vertex;
    }

    std::size_t lower_bound(std::size_t root, std::size_t key) const {
        std::size_t count = 0;
        while (root != NONE) {
            if (vertices_[root].begin < key) {
                count += size(vertices_[root].left) + 1;
                root = vertices_[root].right;
            } else {
                root = vertices_[root].left;
            }
        }
        return count;
    }

    std::size_t ancestor_before(std::size_t root, std::size_t prefix, std::size_t key) const {
        if (root == NONE || prefix == 0 || vertices_[root].maximum_end <= key)
            return NONE;
        const std::size_t left_size = size(vertices_[root].left);
        if (prefix > left_size + 1) {
            const std::size_t found =
                ancestor_before(vertices_[root].right, prefix - left_size - 1, key);
            if (found != NONE)
                return found;
        }
        if (prefix > left_size && vertices_[root].end > key)
            return root;
        return ancestor_before(vertices_[root].left, std::min(prefix, left_size), key);
    }
};

class DecompositionBuilder {
public:
    explicit DecompositionBuilder(const Submap& submap)
        : submap_(submap), parents_(submap.num_nodes(), NONE), positions_(submap.num_nodes(), NONE),
          removed_(submap.num_chords(), false), order_(make_order()) {
        nodes.reserve(2 * submap.num_chords() + 1);
    }

    std::vector<TreeDecompositionNode> nodes;

    std::size_t build() {
        return decompose(order_.balanced(0, submap_.num_nodes()), NONE);
    }

private:
    const Submap& submap_;
    std::vector<std::size_t> parents_;
    std::vector<std::size_t> positions_;
    std::vector<bool> removed_;
    ComponentOrder order_;

    std::vector<ComponentVertex> make_order() {
        std::vector<ComponentVertex> vertices;
        vertices.reserve(submap_.num_nodes());
        std::vector<std::size_t> pending{0};
        while (!pending.empty()) {
            const std::size_t region = pending.back();
            pending.pop_back();
            assert(positions_[region] == NONE && "[C91 §2.2 tex 110]: the dual graph is a tree");
            positions_[region] = vertices.size();
            ComponentVertex vertex;
            vertex.region = region;
            vertex.begin = vertices.size();
            vertex.end = vertex.begin + 1;
            vertices.push_back(vertex);
            for (std::size_t chord_index : submap_.node(region).incident_chords) {
                const Chord& chord = submap_.chord(chord_index);
                const std::size_t neighbor =
                    chord.region[0] == region ? chord.region[1] : chord.region[0];
                if (neighbor == parents_[region])
                    continue;
                assert(parents_[neighbor] == NONE && neighbor != 0 &&
                       "[C91 §2.2 tex 110]: the dual graph has no cycle");
                parents_[neighbor] = region;
                pending.push_back(neighbor);
            }
        }
        assert(vertices.size() == submap_.num_nodes() &&
               "[C91 §2.2 tex 110]: the dual graph is connected");
        for (std::size_t i = vertices.size(); i-- > 0;) {
            const std::size_t parent = parents_[vertices[i].region];
            if (parent != NONE)
                vertices[positions_[parent]].end =
                    std::max(vertices[positions_[parent]].end, vertices[i].end);
        }
        return vertices;
    }

    std::size_t decompose(std::size_t root, std::size_t parent) {
        const std::size_t count = order_.size(root);
        const std::size_t index = nodes.size();
        nodes.push_back({.parent = parent});
        if (count == 1) {
            nodes[index].region_idx = order_.region(root);
            return index;
        }
        const std::size_t centroid = order_.centroid(root);
        const std::size_t region = order_.region(centroid);
        std::size_t largest_branch = 0;
        std::size_t chosen = NONE;
        for (std::size_t chord_index : submap_.node(region).incident_chords) {
            if (removed_[chord_index])
                continue;
            const Chord& chord = submap_.chord(chord_index);
            const std::size_t neighbor =
                chord.region[0] == region ? chord.region[1] : chord.region[0];
            const std::size_t branch = parents_[neighbor] == region
                                           ? order_.subtree_size(root, positions_[neighbor])
                                           : count - order_.subtree_size(root, centroid);
            assert(branch <= count / 2 &&
                   "[C91 §2.3 tex 114]: every centroid branch has at most half the vertices");
            if (branch > largest_branch) {
                largest_branch = branch;
                chosen = chord_index;
            }
        }
        assert(chosen != NONE);
        const std::size_t edge_count = count - 1;
        [[maybe_unused]] const std::size_t maximum_side_edges = edge_count - (edge_count + 3) / 4;
        assert(
            largest_branch - 1 <= maximum_side_edges &&
            count - largest_branch - 1 <= maximum_side_edges &&
            "[C91 §2.3 tex 114]: a centroid edge leaves at most three quarters of the edges on either side");
        const Chord& chord = submap_.chord(chosen);
        const std::size_t child =
            parents_[chord.region[0]] == chord.region[1] ? chord.region[0] : chord.region[1];
        const std::size_t child_vertex = positions_[child];
        const std::size_t begin = order_.rank(child_vertex);
        const std::size_t subtree_count = order_.subtree_size(root, child_vertex);
        const auto [before, rest] = order_.split(root, begin);
        const auto [inside, after] = order_.split(rest, subtree_count);
        const std::size_t outside = order_.join(before, after);
        removed_[chosen] = true;
        nodes[index].chord_idx = chosen;
        const bool first_is_inside = child == chord.region[0];
        const std::size_t left = decompose(first_is_inside ? inside : outside, index);
        const std::size_t right = decompose(first_is_inside ? outside : inside, index);
        nodes[index].left_child = left;
        nodes[index].right_child = right;
        return index;
    }
};

}

void TreeDecomposition::build(const Submap& submap) {
    assert(submap.is_conformal() && "[C91 §2.4(iv)]: tree decomposition requires conformal submap");
    assert(submap.num_nodes() >= 1 && "[C91 §2.3]: submap must have at least one region");
    assert(submap.num_live_nodes() == submap.num_nodes() &&
           submap.num_live_chords() == submap.num_chords() &&
           "[C91 §2.4(iv)]: tree decomposition requires compacted submap");
    assert(submap.num_chords() + 1 == submap.num_nodes() &&
           "[C91 §2.2 tex 110]: the dual graph is a tree");
    DecompositionBuilder builder(submap);
    root_ = builder.build();
    nodes_ = std::move(builder.nodes);
    assert(nodes_.size() == 2 * submap.num_chords() + 1 &&
           "[C91 §2.3]: each chord is an internal node and each region is a leaf");
}

}
