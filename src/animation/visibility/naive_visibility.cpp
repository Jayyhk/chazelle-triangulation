#include "naive_visibility.h"
#include "../merge/granularity.h"
#include "../submap/boundary_geometry.h"
#include "../submap/chord_inventory.h"
#include "../trace.h"
#include "chain.h"

#include <vector>

namespace chazelle::animation {

RayHit naive_first_contact(const Polygon& curve, const Point& p, const SymbolicY& sy, Side dir,
                           std::size_t source_edge) {
    const SourceOffset source_x_offset =
        (source_edge == NONE) ? SourceOffset{}
                              : SourceOffset{perturbed_x_offset(curve, sy, source_edge)};
    RayHit best;
    best.hit = false;
    Exact nearest_distance = 0.0;
    for (std::size_t e = 0; e < curve.num_edges(); ++e) {
        Exact x;
        if (!edge_crossing_x(curve, e, sy, &x))
            continue;

        const auto& ed = curve.edge(e);
        bool asc = symbolic_y_less(symbolic_y_of(curve.vertex(ed.start_idx)),
                                   symbolic_y_of(curve.vertex(ed.end_idx)));
        Side minus_x = asc ? LEFT : RIGHT;
        Side struck = (dir == RIGHT) ? minus_x : (minus_x == LEFT ? RIGHT : LEFT);

        Exact d = (dir == RIGHT) ? (x - p.x) : (p.x - x);

        bool wrapped =
            (d < 0.0) || (d == 0.0 && !perturbed_hit_forward(curve, sy, dir, source_x_offset, e));
        bool better;
        if (!best.hit)
            better = true;
        else if (wrapped != best.wrapped)
            better = !wrapped;
        else if (d != nearest_distance)
            better = d < nearest_distance;
        else
            better = ray_contact_precedes(curve, sy, dir, e, struck, best.edge, best.side);
        if (better) {
            best.hit = true;
            best.x = x;
            best.y = sy.y;
            best.edge = e;
            best.side = struck;
            best.wrapped = wrapped;
            nearest_distance = d;
        }
    }
    return best;
}

namespace {

std::vector<PendingChord> full_visibility_chords(const Polygon& curve) {
    std::vector<PendingChord> chords;
    const std::size_t n = curve.num_vertices();

    auto append_chord = [&](PendingChord pending) {
        chords.push_back(pending);
        if (auto* trace = AnimationTrace::current()) {
            canonicalize_chord(pending, curve);
            Chord chord;
            chord.y = pending.y.y;
            chord.y_tag = pending.y.tag;
            chord.left_edge = pending.left_edge_c;
            chord.left_side = pending.left_side;
            chord.right_edge = pending.right_edge_c;
            chord.right_side = pending.right_side;
            chord.is_null_length = pending.is_null_length;
            trace->chord("visibility_chord", curve, chord);
        }
    };

    auto shoot_from = [&](std::size_t edge, Side side, std::size_t vidx) {
        assert(!is_inside_companion(curve, edge, side, vidx) &&
               "[C91 §2.1 tex 72]: inside-pair duplicates get the null "
               "chord, never a shot (ray_contact_precedes precondition)");
        const SymbolicY vy = symbolic_y_of(curve.vertex(vidx));
        const Side dir = shooting_direction(edge, side, curve);
        Point p{curve.vertex(vidx).x, vy.y, vy.tag};
        RayHit h = naive_first_contact(curve, p, vy, dir, edge);
        if (auto* trace = AnimationTrace::current())
            trace->ray(curve, p, dir, h);

        assert(h.hit && "[C91 §2.1 tex 70]: a chord ray always hits C again");
        PendingChord visibility_chord;
        visibility_chord.y = vy;
        visibility_chord.left_edge_c = edge;
        visibility_chord.left_side = side;
        visibility_chord.right_edge_c = h.edge;
        visibility_chord.right_side = h.side;
        visibility_chord.is_null_length = false;
        append_chord(visibility_chord);
    };

    for (std::size_t v = 0; v < n; ++v) {
        if (curve.is_endpoint(v)) {
            const std::size_t e = (v == 0) ? 0 : curve.num_edges() - 1;
            shoot_from(e, LEFT, v);
            shoot_from(e, RIGHT, v);
        } else if (curve.is_y_extremum(v)) {
            const bool next_left_inside = is_inside_companion(curve, v, LEFT, v);
            assert(next_left_inside != is_inside_companion(curve, v, RIGHT, v) &&
                   "[C91 §2.1 tex 72]: exactly one face of the next edge "
                   "is the side facing the inside companion");
            const Side inside_next = next_left_inside ? LEFT : RIGHT;
            PendingChord null_chord;
            null_chord.y = symbolic_y_of(curve.vertex(v));
            null_chord.left_edge_c = null_chord.right_edge_c = v;
            null_chord.left_side = null_chord.right_side = inside_next;
            null_chord.is_null_length = true;
            append_chord(null_chord);

            const bool prev_left_inside = is_inside_companion(curve, v - 1, LEFT, v);
            assert(prev_left_inside != is_inside_companion(curve, v - 1, RIGHT, v) &&
                   "[C91 §2.1 tex 72]: exactly one face of the previous "
                   "edge is the side facing the inside companion");
            const Side outside_prev = prev_left_inside ? RIGHT : LEFT;
            const Side outside_next = next_left_inside ? RIGHT : LEFT;
            shoot_from(v - 1, outside_prev, v);
            shoot_from(v, outside_next, v);
        } else {
            assert(shooting_direction(v - 1, LEFT, curve) == shooting_direction(v, LEFT, curve) &&
                   shooting_direction(v - 1, RIGHT, curve) == shooting_direction(v, RIGHT, curve) &&
                   "[C91 §2.1 tex 72]: the chord direction is continuous "
                   "through a non-extremum vertex");
            shoot_from(v - 1, LEFT, v);
            shoot_from(v - 1, RIGHT, v);
        }
    }
    return chords;
}

}

Submap build_full_visibility_map(const Polygon& curve) {
    Submap submap;
    build_submap_from_chords(submap, curve, full_visibility_chords(curve));

    assert(submap.is_conformal() && "[C91 §2.1 tex 70]: the full V(C) is conformal (trapezoids)");

    assert(submap.is_semigranular(1) && "[C91 §2.1 tex 72]: full V(C) arcs span at most one edge");

    submap.build_tree_decomposition();
    return submap;
}

std::size_t canonical_granularity(std::size_t num_edges) {
    assert(num_edges >= 1 && "[C91 §2.1]: a curve has ≥ 1 edge");

    std::size_t k = 0;
    while ((std::size_t{1} << k) < num_edges)
        ++k;
    return std::size_t{1} << ceil_beta(k);
}

Submap build_canonical_submap_naive(const Polygon& curve) {
    Submap submap = build_full_visibility_map(curve);
    const std::size_t granularity = canonical_granularity(curve.num_edges());

    enforce_granularity(submap, curve, granularity);

    submap.normalize(curve);

    assert(submap.is_conformal() && submap.is_granular(granularity, curve) &&
           !submap.tree_decomposition().empty() &&
           "[C91 §4.1 tex 327]: canonical = 2^{⌈β⌈log m⌉⌉}-granular, "
           "conformal, normal form");
    return submap;
}

}
