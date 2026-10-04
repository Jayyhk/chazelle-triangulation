#include "boundary_geometry.h"

#include <algorithm>

namespace chazelle {

RegionArcs collect_region_arcs(const Submap& submap, std::size_t region) {
    assert(region < submap.num_nodes() && !submap.node(region).dead);

    RegionArcs out;

    assert(submap.node(region).degree() <= 4 &&
           "[C91 §2.3]: conformal regions MUST have degree ≤ 4");
    auto check_adj = [&](const Chord::AdjArcs& adj) {
        for (std::size_t k = 0; k < adj.count; ++k) {
            std::size_t ai = adj.arcs[k];
            assert(ai < submap.num_arcs() && !submap.arc(ai).dead);
            if (submap.arc(ai).region_node == region) {
                bool dup = false;
                for (std::size_t i = 0; i < out.count; ++i)
                    if (out.arcs[i] == ai) {
                        dup = true;
                        break;
                    }
                if (!dup)
                    out.push(ai);
            }
        }
    };

    for (std::size_t ci : submap.node(region).incident_chords) {
        assert(ci < submap.num_chords());
        assert(!submap.chord(ci).dead && "[C91 §2.4]: normal-form (compacted) submap must have no "
                                         "dead chords in incident_chords");
        check_adj(submap.chord(ci).left_adj);
        check_adj(submap.chord(ci).right_adj);
    }

    if (submap.node(region).incident_chords.empty()) {
        assert(submap.num_live_chords() == 0 &&
               "[C91 §2.2 tex 102]: a chord-free region exists only in "
               "the chordless (single-region) submap");
        assert(submap.start_arc != NONE && submap.start_arc == submap.end_arc &&
               submap.start_arc < submap.num_arcs() && !submap.arc(submap.start_arc).dead &&
               submap.arc(submap.start_arc).region_node == region &&
               "[C91 §2.4(iii) tex 138]: the chordless submap's closed "
               "arc is the endpoint arc");
        out.push(submap.start_arc);
    }

    assert(out.count == std::max<std::size_t>(submap.node(region).degree(), 1) &&
           "[C91 §2.2 tex 96]: arc count == region degree (≥ 1)");

    return out;
}

Side shooting_direction(std::size_t edge, Side side, const Polygon& curve) {
    assert(edge < curve.num_edges());
    const auto& e = curve.edge(edge);
    SymbolicY start_y = symbolic_y_of(curve.vertex(e.start_idx));
    SymbolicY end_y = symbolic_y_of(curve.vertex(e.end_idx));
    bool edge_ascending = symbolic_y_less(start_y, end_y);

    if (side == LEFT)
        return edge_ascending ? LEFT : RIGHT;
    else
        return edge_ascending ? RIGHT : LEFT;
}

bool chord_runs_through_infinity(const Polygon& curve, const Chord& c) {
    assert(!c.is_null_length && "[C91 §2.2 tex 96]: null-length chords occupy a single point "
                                "and run through nothing");
    Exact left_x = edge_x_at_y(curve, c.left_edge, c.symbolic_y());
    Exact right_x = edge_x_at_y(curve, c.right_edge, c.symbolic_y());
    assert(left_x <= right_x && "chord slots are ordered by ascending x");
    if (left_x == right_x) {
        const Exact left_x_offset = perturbed_x_offset(curve, c.symbolic_y(), c.left_edge);
        const Exact right_x_offset = perturbed_x_offset(curve, c.symbolic_y(), c.right_edge);
        if (left_x_offset == right_x_offset)
            return true;
        const bool left_is_west = left_x_offset < right_x_offset;
        const std::size_t west_edge = left_is_west ? c.left_edge : c.right_edge;
        const Side west_side = left_is_west ? c.left_side : c.right_side;
        return shooting_direction(west_edge, west_side, curve) == LEFT;
    }
    return shooting_direction(c.left_edge, c.left_side, curve) == LEFT;
}

bool arc_starts_at_chord_slot(const Submap& submap, const Polygon& curve, const Chord& c,
                              bool left_slot, std::size_t arc_idx) {
    const Arc& a = submap.arc(arc_idx);
    std::size_t endpoint_edge = left_slot ? c.left_edge : c.right_edge;
    Side endpoint_side = left_slot ? c.left_side : c.right_side;
    return a.first_edge == endpoint_edge && a.first_side == endpoint_side &&
           symbolic_y_equal(submap.arc_start_symbolic_y(arc_idx, curve), c.symbolic_y());
}

}
