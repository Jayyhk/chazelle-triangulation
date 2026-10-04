#include "../submap/chord_inventory.h"
#include "fusion.h"

#include <utility>

namespace chazelle {

void rebuild_submap(Submap& submap, const Polygon& curve, const Submap& first_submap,
                    const Polygon& first_curve, const Submap& second_submap,
                    [[maybe_unused]] const Polygon& second_curve, const FusionState& state1,
                    const FusionState& state2) {
    assert(submap.num_nodes() == 0 && submap.num_arcs() == 0 && submap.num_chords() == 0 &&
           "[C91 §3.1 tex 226]: rebuild_submap requires a fresh output submap");
    assert(curve.num_edges() == first_curve.num_edges() + second_curve.num_edges() &&
           "[C91 §3 tex 160]: C = C₁ ∪ C₂ shares one vertex");

    const std::size_t n_c1_edges = first_curve.num_edges();

    using Pending = PendingChord;
    std::vector<Pending> pending;
    pending.reserve(state1.chords.size() + state2.chords.size() + first_submap.num_live_chords() +
                    second_submap.num_live_chords());

    auto edge_in_c = [&](std::size_t edge, bool on_first_curve) {
        return on_first_curve ? edge : (n_c1_edges + edge);
    };

    auto ingest_discovered = [&](const FusionState& st, bool pass_starts_on_first_curve) {
        for (const auto& dc : st.chords) {
            pending.push_back(
                {dc.y,
                 edge_in_c(dc.left_edge, dc.left_on_first_curve == pass_starts_on_first_curve),
                 dc.left_side,
                 edge_in_c(dc.right_edge, dc.right_on_first_curve == pass_starts_on_first_curve),
                 dc.right_side, false});
        }
    };
    ingest_discovered(state1, true);
    ingest_discovered(state2, false);

    const std::size_t junction_vidx_in_c = first_curve.num_vertices() - 1;
    [[maybe_unused]] const std::size_t junction_tag = first_curve.vertex(junction_vidx_in_c).index;
    assert(junction_tag == second_curve.vertex(0).index &&
           "[C91 §3 tex 160]: junction shared between C₁ and C₂");

    assert(junction_vidx_in_c >= 1 && junction_vidx_in_c + 1 < curve.num_vertices() &&
           "[C91 §3 tex 160]: junction is interior to C");

    assert(state1.invalidated_first_chords.size() == first_submap.num_chords() &&
           state2.invalidated_second_chords.size() == first_submap.num_chords() &&
           state1.invalidated_second_chords.size() == second_submap.num_chords() &&
           state2.invalidated_first_chords.size() == second_submap.num_chords() &&
           "[C91 §3.1 tex 224]: both passes account for every input chord");
    auto bit_at = [](const std::vector<bool>& b, std::size_t i) {
        assert(i < b.size());
        return b[i];
    };
    auto s1_invalid = [&](std::size_t i) {
        return bit_at(state1.invalidated_first_chords, i) ||
               bit_at(state2.invalidated_second_chords, i);
    };
    auto s2_invalid = [&](std::size_t i) {
        return bit_at(state1.invalidated_second_chords, i) ||
               bit_at(state2.invalidated_first_chords, i);
    };

    auto ingest_old = [&](const Submap& submap, bool on_first_curve, auto is_invalid) {
        for (std::size_t ci = 0; ci < submap.num_chords(); ++ci) {
            const Chord& c = submap.chord(ci);
            if (c.dead)
                continue;

            assert(!(c.is_null_length && c.y_tag == junction_tag) &&
                   "[C91 §2.1 tex 72 case 3]: Sᵢ cannot contain a "
                   "null-length chord sourced at its own C-endpoint");
            if (is_invalid(ci))
                continue;
            pending.push_back({c.symbolic_y(), edge_in_c(c.left_edge, on_first_curve), c.left_side,
                               edge_in_c(c.right_edge, on_first_curve), c.right_side,
                               c.is_null_length});
        }
    };
    ingest_old(first_submap, true, s1_invalid);
    ingest_old(second_submap, false, s2_invalid);

    if (is_local_y_extremum(curve.vertex(junction_vidx_in_c - 1), curve.vertex(junction_vidx_in_c),
                            curve.vertex(junction_vidx_in_c + 1))) {
        const Point& u = curve.vertex(junction_vidx_in_c - 1);
        const Point& v = curve.vertex(junction_vidx_in_c);
        const Point& w = curve.vertex(junction_vidx_in_c + 1);
        const bool is_max = point_y_above(v, u) && point_y_above(v, w);

        const bool prev_left = curve.previous_branch_left(junction_vidx_in_c);

        const Side minus_x_face = is_max ? RIGHT : LEFT;
        const Side plus_x_face = is_max ? LEFT : RIGHT;
        const Side inside = prev_left ? minus_x_face : plus_x_face;

        Pending p;
        p.y = symbolic_y_of(v);
        p.left_edge_c = p.right_edge_c = n_c1_edges;
        p.left_side = p.right_side = inside;
        p.is_null_length = true;
        pending.push_back(p);
    }

    build_submap_from_chords(submap, curve, std::move(pending));
}

}
