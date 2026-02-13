#include "triangulate/monotone.h"
#include "geometry/perturbation.h"

namespace chazelle {

void triangulate_monotone(std::vector<VertexNode>& nodes,
                          std::size_t first,
                          std::size_t last,
                          std::vector<Triangle>& out) {
    // `last` defines the sub-polygon end for the caller (algo2);
    // the actual walk uses the circular linked list rooted at `first`.
    // We use `last` as a stop sentinel to avoid walking past the
    // sub-polygon boundary.
    (void)last;

    // Count vertices in this sub-polygon.
    std::size_t num_verts = 0;
    {
        std::size_t v = first;
        do {
            ++num_verts;
            v = nodes[v].next;
        } while (v != first);
    }

    if (num_verts < 3) return;

    // Determine start vertex.
    // Find the y-max and y-min of the sub-polygon.
    std::size_t y_max = first, y_min = first;
    {
        std::size_t v = first;
        do {
            Point pv{nodes[v].x, nodes[v].y, nodes[v].index};
            Point pmax{nodes[y_max].x, nodes[y_max].y, nodes[y_max].index};
            Point pmin{nodes[y_min].x, nodes[y_min].y, nodes[y_min].index};
            if (perturbed_y_compare(pv, pmax) > 0) y_max = v;
            if (perturbed_y_compare(pv, pmin) < 0) y_min = v;
            v = nodes[v].next;
        } while (v != first);
    }

    // Determine if the monotone chain is on the right or left.
    // If y_max → next → … → y_min goes in decreasing y, chain is on right.
    // Start = topmost if chain on right, bottommost if chain on left.
    //
    // Check: walk from y_max via next.  If we hit y_min before looping,
    // the next-chain is the monotone side.  If y_max.next has lower y
    // than y_max.prev (other than the edge vertex), chain is on the
    // "next" side (right).
    //
    // Determine start vertex: y_max if the base edge (y_max, y_min)
    // goes via next (i.e., they're adjacent via the base edge),
    // y_min otherwise.
    std::size_t start;
    if (nodes[y_max].next == y_min || nodes[y_min].prev == y_max) {
        // Base edge goes y_max → y_min via next.  The monotone chain
        // wraps the other way (via prev), placing it on the left side.
        // Per Algorithm 3: "topmost if chain on right, bottommost if
        // on left."  Chain on left → start = bottommost.
        start = y_min;
    } else {
        // Base edge goes y_min → y_max.  Chain on the right.
        start = y_max;
    }

    // Walk from next(start), removing convex vertices.
    // Use a stack to avoid re-scanning reflex vertices:
    // push reflex vertices, pop when a convex vertex enables ear removal.
    // Each vertex is pushed and popped at most once, giving O(n) total.
    std::vector<std::size_t> stack;
    stack.push_back(start);
    stack.push_back(nodes[start].next);

    std::size_t current = nodes[nodes[start].next].next;
    std::size_t iterations_left = 2 * num_verts; // safety bound (each vertex pushed + popped at most once)

    while (num_verts >= 3 && iterations_left-- > 0) {
        std::size_t top = stack.back();

        Point pp{nodes[stack[stack.size() >= 2 ? stack.size() - 2 : 0]].x,
                 nodes[stack[stack.size() >= 2 ? stack.size() - 2 : 0]].y,
                 nodes[stack[stack.size() >= 2 ? stack.size() - 2 : 0]].index};
        Point pc{nodes[top].x, nodes[top].y, nodes[top].index};
        Point pn{nodes[current].x, nodes[current].y, nodes[current].index};

        // Check if top of stack forms a convex ear with its neighbors.
        if (stack.size() >= 2 && orient2d(pp, pc, pn) >= 0.0) {
            // Emit triangle and pop the ear vertex.
            std::size_t prev_on_stack = stack[stack.size() - 2];
            out.push_back({nodes[prev_on_stack].index,
                           nodes[top].index,
                           nodes[current].index});

            // Remove top from linked list.
            nodes[prev_on_stack].next = current;
            nodes[current].prev = prev_on_stack;
            stack.pop_back();
            --num_verts;
            // Don't advance current — check if the new top is also convex.
        } else {
            // Reflex or base case — push current and advance.
            stack.push_back(current);
            current = nodes[current].next;
        }
    }
}

} // namespace chazelle
