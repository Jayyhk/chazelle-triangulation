#include "triangulate/triangulate.h"

namespace chazelle {

namespace {

/// The `diagonal` function (Fournier-Montuno §4.2):
/// If the vertex's trapezoid is Class B (top and bottom are not polygon-
/// adjacent), return the bottom vertex index and clear the trapezoid.
/// Otherwise return FM_NONE.
std::size_t diagonal(std::vector<VertexNode>& nodes,
                     std::vector<Trapezoid>& traps,
                     std::size_t v) {
    // Per FM Algorithm 1: type-2 vertices (local y-minima) can have
    // two trapezoids.  Check both slots.
    for (int slot = 0; slot < 2; ++slot) {
        std::size_t& ti_ref = (slot == 0) ? nodes[v].trapezoid_idx
                                          : nodes[v].trapezoid_idx2;
        std::size_t ti = ti_ref;
        if (ti == FM_NONE) continue;

        const Trapezoid& trap = traps[ti];
        std::size_t top = trap.top_vertex;
        std::size_t bot = trap.bottom_vertex;

        // Determine which vertex is `v` and which is the "other."
        std::size_t other = (top == v) ? bot : top;

        // Class A check: are top and bottom adjacent in the polygon?
        if (nodes[top].next == bot || nodes[bot].next == top) {
            // Class A — no diagonal needed from this slot.
            continue;
        }

        // Class B — insert diagonal.  Clear this trapezoid slot.
        ti_ref = FM_NONE;
        return other;
    }
    return FM_NONE;
}

/// Recursive Algorithm 2 (Fournier-Montuno §4.2).
///
/// Walk the sub-polygon linked list looking for a Class B diagonal.
/// If found, split the polygon into two sub-polygons along the diagonal
/// and recurse on each.  If no diagonal is found, the polygon is
/// unimonotone — triangulate it directly via Algorithm 3.
///
/// O(n) total: each vertex is visited exactly once per recursion level.
/// The do...while(current!=first) loop terminates naturally when the
/// walk returns to the starting vertex.  No done-flags needed —
/// diagonal() is idempotent (clears trapezoid_idx on first call,
/// returns FM_NONE on subsequent calls).
void algo2(std::vector<VertexNode>& nodes,
           std::vector<Trapezoid>& traps,
           std::size_t first,
           [[maybe_unused]] std::size_t last,
           std::vector<Triangle>& out) {
    std::size_t current = first;

    // Walk the polygon looking for a Class B diagonal.
    do {
        std::size_t bottom = diagonal(nodes, traps, current);
        if (bottom != FM_NONE) {
            // Split the linked list into two disjoint sub-polygons
            // at the diagonal (current, bottom).  O(1) per split.

            std::size_t save_current_next = nodes[current].next;
            std::size_t save_current_prev = nodes[current].prev;
            std::size_t save_bottom_next = nodes[bottom].next;
            std::size_t save_bottom_prev = nodes[bottom].prev;

            // --- Sub-polygon 1: bottom → … → current → bottom ---
            nodes[current].next = bottom;
            nodes[bottom].prev  = current;

            algo2(nodes, traps, bottom, current, out);

            // --- Sub-polygon 2: current → save_next → … → bottom → current ---
            nodes[current].next = save_current_next;
            nodes[save_current_next].prev = current;
            nodes[bottom].prev = save_bottom_prev;
            nodes[save_bottom_prev].next = bottom;
            nodes[bottom].next = current;
            nodes[current].prev = bottom;

            algo2(nodes, traps, current, bottom, out);

            // Restore original links for correctness of parent calls.
            nodes[current].next = save_current_next;
            nodes[current].prev = save_current_prev;
            nodes[bottom].next  = save_bottom_next;
            nodes[bottom].prev  = save_bottom_prev;
            nodes[save_current_next].prev = current;
            nodes[save_current_prev].next = current;
            nodes[save_bottom_next].prev  = bottom;
            nodes[save_bottom_prev].next  = bottom;
            return;
        }

        current = nodes[current].next;
    } while (current != first);

    // No Class B trapezoid found — polygon is unimonotone.
    triangulate_monotone(nodes, first, last, out);
}

} // anonymous namespace

std::vector<Triangle>
fm_triangulate(std::vector<VertexNode>& nodes,
               std::vector<Trapezoid>& traps,
               std::size_t first,
               std::size_t last) {
    std::vector<Triangle> triangles;
    triangles.reserve(nodes.size()); // n−2 triangles at most

    algo2(nodes, traps, first, last, triangles);
    return triangles;
}

} // namespace chazelle
