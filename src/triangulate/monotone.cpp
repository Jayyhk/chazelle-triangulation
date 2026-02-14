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

    // Walk from next(start), removing convex vertices using backtracking.
    // Fournier-Montuno Algorithm 3:
    //   u = next(start)
    //   while n >= 3:
    //     if angle(prev(u), u, next(u)) is convex:
    //       output triangle
    //       remove u
    //       u = prev(u)
    //     else:
    //       u = next(u)

    std::size_t u = nodes[start].next;
    std::size_t safety_counter = 2 * num_verts * num_verts; // Upper bound for safety

    while (num_verts >= 3 && safety_counter-- > 0) {
        std::size_t p = nodes[u].prev;
        std::size_t n = nodes[u].next;

        Point prev_pt{nodes[p].x, nodes[p].y, nodes[p].index};
        Point curr_pt{nodes[u].x, nodes[u].y, nodes[u].index};
        Point next_pt{nodes[n].x, nodes[n].y, nodes[n].index};

        // Determine convexity.
        // For a simple polygon, "convex" means the internal angle is < 180 degrees.
        // Since we traverse in CCW order, this corresponds to a specific orientation.
        // orient2d returns > 0 for CCW turn (left turn), which is convex for CCW traversal.
        // orient2d returns < 0 for CW turn (right turn), which is reflex.
        // == 0 overlaps, treat as convex to remove degenerate vertices.
        if (orient2d(prev_pt, curr_pt, next_pt) >= 0.0) {
            // Convex (or collinear). Output triangle.
            out.push_back({nodes[p].index, nodes[u].index, nodes[n].index});

            // Remove u from the linked list.
            nodes[p].next = n;
            nodes[n].prev = p;
            
            // Decrement vertex count.
            --num_verts;

            // Backtrack: u = prev(u).
            // Algorithm 3 explicitly handles the case where the removed vertex was the start:
            // "if current = first then current = next(first) else current = save;"
            if (u == first) {
                u = n; // next(first)
                // Note: we do not update 'first' because it is a value argument
                // and serves only as a check for this specific boundary condition.
            } else {
                u = p; // save (prev(current))
            }
        } else {
            // Reflex. Advance.
            u = n;
        }
    }
}

} // namespace chazelle
