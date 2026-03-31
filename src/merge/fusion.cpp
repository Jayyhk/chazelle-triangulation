/// src/merge/fusion.cpp

#include "fusion.h"

#include <algorithm>

namespace chazelle {

std::vector<FusionVertex> build_fusion_sequence(const Submap& S,
                                                 const Polygon& C) {
    std::size_t n = C.num_vertices();
    assert(n >= 2);

    // [C91 §3.1]: The junction vertex is C₁ ∩ C₂ = last vertex of C₁.
    std::size_t junction_v = n - 1;
    std::size_t junction_edge = junction_v - 1; // last edge of C₁

    // [C91 §3.1]: "Let a_{m+1} and a₀ be the companion vertices,
    // as they appear next to each other clockwise around ∂C₁,
    // resulting from the duplication of the vertex C₁ ∩ C₂."
    //
    // Clockwise ∂C₁: LEFT ascending → junction → RIGHT descending.
    // a_{m+1} = LEFT companion at junction (end of LEFT traversal).
    // a₀ = RIGHT companion at junction (start of RIGHT traversal).
    //
    // Both have the junction vertex's y-coordinate and SoS tag.
    SymbolicY junction_y = symbolic_y_of(C.vertex(junction_v));

    FusionVertex a_0;
    a_0.y = junction_y;
    a_0.edge = junction_edge;
    a_0.side = RIGHT;  // RIGHT companion: start of clockwise tour
    a_0.chord_idx = NONE;
    a_0.is_left_endpoint = false;
    a_0.is_companion = true;

    FusionVertex a_m1;
    a_m1.y = junction_y;
    a_m1.edge = junction_edge;
    a_m1.side = LEFT;   // LEFT companion: end of clockwise tour
    a_m1.chord_idx = NONE;
    a_m1.is_left_endpoint = false;
    a_m1.is_companion = true;

    // [C91 §3.1]: "Let a₁, a₂, ..., aₘ be the canonical vertex
    // enumeration of S₁.  Recall that this enumerates the exit chord
    // endpoints in S₁ as we encounter them going clockwise around ∂C₁."
    //
    // Each chord has two endpoints (LEFT and RIGHT).  Collect them all.
    std::vector<FusionVertex> seq;
    seq.reserve(2 * S.num_live_chords() + 2);

    for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
        if (S.chord(ci).dead) continue;
        const auto& c = S.chord(ci);

        FusionVertex left_ep;
        left_ep.y = c.symbolic_y();
        left_ep.edge = c.left_edge;
        left_ep.side = c.left_side;
        left_ep.chord_idx = ci;
        left_ep.is_left_endpoint = true;
        left_ep.is_companion = false;
        seq.push_back(left_ep);

        FusionVertex right_ep;
        right_ep.y = c.symbolic_y();
        right_ep.edge = c.right_edge;
        right_ep.side = c.right_side;
        right_ep.chord_idx = ci;
        right_ep.is_left_endpoint = false;
        right_ep.is_companion = false;
        seq.push_back(right_ep);
    }

    // Sort in clockwise ∂C₁ order starting from a₀.
    //
    // Clockwise ∂C₁ from a₀:
    //   RIGHT side: junction_edge descending to edge 0
    //   LEFT side: edge 0 ascending to junction_edge
    //
    // Assign a traversal position for sorting:
    //   RIGHT vertices: pos = (junction_edge - edge)  [0 at junction, increasing toward start]
    //   LEFT vertices: pos = (junction_edge + 1) + edge  [continues after RIGHT]
    //
    // Within the same edge, sort by y along the ∂C traversal direction.
    std::size_t right_half_len = junction_edge + 1; // edges 0..junction_edge

    auto traversal_pos = [&](const FusionVertex& v) -> std::size_t {
        if (v.side == RIGHT) {
            return (junction_edge - v.edge);
        } else {
            return right_half_len + v.edge;
        }
    };

    std::sort(seq.begin(), seq.end(),
        [&](const FusionVertex& a, const FusionVertex& b) {
            std::size_t pa = traversal_pos(a);
            std::size_t pb = traversal_pos(b);
            if (pa != pb) return pa < pb;
            // Same edge: sort by y along ∂C traversal direction.
            // RIGHT side descends in edge index but the y-order within
            // an edge depends on the edge's geometric direction.
            // LEFT side: ∂C traversal follows edge direction (ascending).
            // RIGHT side: ∂C traversal is reversed.
            // For vertices on the same edge, use symbolic y:
            // on LEFT, earlier in traversal = lower y on ascending edges
            // (but depends on edge direction — just use y for now,
            // ties broken by SoS).
            return symbolic_y_less(a.y, b.y);
        });

    // [C91 §3.1]: "we assume that [a₀ and a_{m+1}] are not [in the
    // sequence] and therefore add them to the sequence."
    // Prepend a₀, append a_{m+1}.
    std::vector<FusionVertex> result;
    result.reserve(seq.size() + 2);
    result.push_back(a_0);
    result.insert(result.end(), seq.begin(), seq.end());
    result.push_back(a_m1);

    return result;
}

} // namespace chazelle
