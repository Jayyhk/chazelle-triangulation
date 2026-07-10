// src/visibility/naive_visibility.cpp — [C91 §4.1 tex 345] naive base
// case: full V(C) per [C91 §2.1 tex 68–72], canonicalized via §3.3.

#include "naive_visibility.h"

#include "../merge/fusion.h"       // PendingChord, build_submap_from_chords,
                                   // shooting_direction
#include "../merge/granularity.h"  // enforce_granularity
#include "chain.h"                 // ceil_beta ([C91 §4 tex 323])

#include <vector>

namespace chazelle {

// ════════════════════════════════════════════════════════════════
//  naive_first_contact — [C91 §2.1 tex 70]
// ════════════════════════════════════════════════════════════════

RayHit naive_first_contact(const Polygon& C, Point p, SymbolicY sy,
                           Side dir, std::size_t source_edge) {
    const double src_offset = (source_edge == NONE)
                                ? SOURCE_OFFSET_NONE
                                : perturbed_x_offset(C, sy, source_edge);
    RayHit best;
    best.hit = false;
    double best_d = 0.0;
    for (std::size_t e = 0; e < C.num_edges(); ++e) {
        double x;
        if (!edge_crossing_x(C, e, sy, &x)) continue;
        // [C91 §2.1 tex 72]: traveling RIGHT the ray arrives from −x
        // and strikes the −x-facing wall; LEFT symmetrically.
        const auto& ed = C.edge(e);
        bool asc = symbolic_y_less(symbolic_y_of(C.vertex(ed.start_idx)),
                                   symbolic_y_of(C.vertex(ed.end_idx)));
        Side minus_x = asc ? LEFT : RIGHT;
        Side struck = (dir == RIGHT) ? minus_x
                                     : (minus_x == LEFT ? RIGHT : LEFT);
        // [C91 §2.1 tex 70]: wrap metric — direct contacts (d > 0)
        // precede all through-infinity ones (d ≤ 0; d == 0 is the
        // source's own position, met last after a full wrap); within a
        // class the ray order is ascending d.
        double d = (dir == RIGHT) ? (x - p.x) : (p.x - x);
        // [C91 §2 tex 47]: raw d = 0 shares the source's raw point;
        // with the source's edge known, comparing the perturbed
        // x-offsets decides forward vs met-after-wrap
        // (polygon.h::perturbed_hit_forward).
        bool wrapped = (d < 0.0) ||
                       (d == 0.0 &&
                        !(has_source_offset(src_offset) &&
                          perturbed_hit_forward(C, sy, dir, src_offset, e)));
        bool better;
        if (!best.hit)
            better = true;
        else if (wrapped != best.wrapped)
            better = !wrapped;
        else if (d != best_d)
            better = d < best_d;
        else
            // Same raw travel position: [C91 §2 tex 47] tie-break.
            better = ray_contact_precedes(C, sy, dir, e, struck,
                                          best.edge, best.side);
        if (better) {
            best.hit = true;
            best.x = x;
            best.y = sy.y;
            best.edge = e;
            best.side = struck;
            best.wrapped = wrapped;
            best_d = d;
        }
    }
    return best;
}

// ════════════════════════════════════════════════════════════════
//  Full V(C) — [C91 §2.1 tex 68–72]
// ════════════════════════════════════════════════════════════════

namespace {

// [C91 §2.1 tex 68–72]: the chord inventory of the full V(C).  "We have
// now ensured that each vertex of ∂C is incident upon exactly one chord
// of the visibility map" (tex 72): each companion vertex of ∂C shoots
// its unique chord direction, EXCEPT the extremum inside pairs, whose
// chord is the synthesized null-length one ("the pair on the 'inside'
// of the turn gives rise to a chord of null length").  Chords are
// discovered from both companion endpoints where both are vertices of
// ∂C; build_submap_from_chords deduplicates.
std::vector<PendingChord> full_visibility_chords(const Polygon& C) {
    std::vector<PendingChord> chords;
    const std::size_t n = C.num_vertices();

    // One chord ray from the ∂C companion at vertex `vidx`, lying on
    // (edge, side); direction is the point's unique chord direction
    // ([C91 §2.1 tex 72]).
    auto shoot_from = [&](std::size_t edge, Side side, std::size_t vidx) {
        assert(!is_inside_companion(C, edge, side, vidx) &&
               "[C91 §2.1 tex 72]: inside-pair duplicates get the null "
               "chord, never a shot (ray_contact_precedes precondition)");
        const SymbolicY vy = symbolic_y_of(C.vertex(vidx));
        const Side dir = shooting_direction(edge, side, C);
        Point p{C.vertex(vidx).x, vy.y, vy.tag};
        RayHit h = naive_first_contact(C, p, vy, dir, edge);
        // [C91 §2.1 tex 70]: "it would actually wrap around in the
        // spherical plane until it hits C again" — at worst the ray
        // returns to the source vertex's opposite wall (d = 0).
        assert(h.hit &&
               "[C91 §2.1 tex 70]: a chord ray always hits C again");
        PendingChord pc;
        pc.y = vy;
        pc.left_edge_c = edge;
        pc.left_side = side;
        pc.right_edge_c = h.edge;
        pc.right_side = h.side;
        pc.is_null_length = false;
        chords.push_back(pc);
    };

    for (std::size_t v = 0; v < n; ++v) {
        if (C.is_endpoint(v)) {
            // [C91 §2.1 tex 72] case 3: "for each endpoint of C we
            // create two companion vertices ... next to each other as
            // well as on both sides of C" — one chord each.
            const std::size_t e = (v == 0) ? 0 : C.num_edges() - 1;
            shoot_from(e, LEFT, v);
            shoot_from(e, RIGHT, v);
        } else if (C.is_y_extremum(v)) {
            // [C91 §2.1 tex 72] case 2: four ∂C vertices.  Inside pair:
            // the null-length chord, recorded on the NEXT edge's wedge
            // face — the same convention as rebuild_submap's junction
            // synthesis ([C91 §3.1 tex 224], fusion.cpp).
            const bool next_left_inside = is_inside_companion(C, v, LEFT, v);
            assert(next_left_inside != is_inside_companion(C, v, RIGHT, v) &&
                   "[C91 §2.1 tex 72]: exactly one face of the next edge "
                   "is the wedge face");
            const Side inside_next = next_left_inside ? LEFT : RIGHT;
            PendingChord nc;
            nc.y = symbolic_y_of(C.vertex(v));
            nc.left_edge_c = nc.right_edge_c = v;
            nc.left_side = nc.right_side = inside_next;
            nc.is_null_length = true;
            chords.push_back(nc);

            // Outside pair: one chord ray each ("two horizontal
            // segments from each vertex, one in each direction",
            // tex 70 — realized at the outside duplicates).
            const bool prev_left_inside =
                is_inside_companion(C, v - 1, LEFT, v);
            assert(prev_left_inside !=
                       is_inside_companion(C, v - 1, RIGHT, v) &&
                   "[C91 §2.1 tex 72]: exactly one face of the previous "
                   "edge is the wedge face");
            const Side outside_prev = prev_left_inside ? RIGHT : LEFT;
            const Side outside_next = next_left_inside ? RIGHT : LEFT;
            shoot_from(v - 1, outside_prev, v);
            shoot_from(v, outside_next, v);
        } else {
            // [C91 §2.1 tex 72] case 1: two companions, one per ∂C
            // side; both incident edges contain each companion — use
            // the LOWER edge (the label canonicalization of
            // build_submap_from_chords).  The chord direction is
            // side-determined and agrees across the two incident edges
            // (both pass y-monotonically through v).
            assert(shooting_direction(v - 1, LEFT, C) ==
                       shooting_direction(v, LEFT, C) &&
                   shooting_direction(v - 1, RIGHT, C) ==
                       shooting_direction(v, RIGHT, C) &&
                   "[C91 §2.1 tex 72]: the chord direction is continuous "
                   "through a non-extremum vertex");
            shoot_from(v - 1, LEFT, v);
            shoot_from(v - 1, RIGHT, v);
        }
    }
    return chords;
}

} // namespace

Submap build_full_visibility_map(const Polygon& C) {
    Submap S;
    build_submap_from_chords(S, C, full_visibility_chords(C));

    // [C91 §2.1 tex 70]: V(C)'s regions are "triangles or trapezoids
    // with two horizontal segments and two nonhorizontal segments" —
    // each horizontal side holds the collinear chords of a single
    // vertex's level ([C91 §2 tex 47]: same symbolic y ⟹ same source
    // vertex, ≤ 2 chords + the null slit), so every region has degree
    // ≤ 4: V(C) is conformal ([C91 §2.3 tex 114]).
    assert(S.is_conformal() &&
           "[C91 §2.1 tex 70]: the full V(C) is conformal (trapezoids)");

    // [C91 §2.1 tex 72]: every vertex of ∂C is a chord endpoint, so no
    // arc of V(C) spans a vertex interiorly — every arc covers at most
    // one edge per leg: V(C) is 1-semigranular ([C91 §2.2 tex 106]).
    assert(S.is_semigranular(1) &&
           "[C91 §2.1 tex 72]: full V(C) arcs span at most one edge");

    S.build_tree_decomposition();
    return S;
}

// ════════════════════════════════════════════════════════════════
//  Canonical submaps — [C91 §4.1 tex 338]
// ════════════════════════════════════════════════════════════════

std::size_t canonical_granularity(std::size_t num_edges) {
    assert(num_edges >= 1 && "[C91 §2.1]: a curve has ≥ 1 edge");
    // ⌈log₂ m⌉ (smallest k with 2^k ≥ m), then 2^{⌈βk⌉} in exact
    // integer arithmetic ([C91 §4 tex 323]: β = 1/5 exact rational).
    std::size_t k = 0;
    while ((std::size_t{1} << k) < num_edges) ++k;
    return std::size_t{1} << ceil_beta(k);
}

Submap build_canonical_submap_naive(const Polygon& C) {
    Submap S = build_full_visibility_map(C);
    const std::size_t gamma = canonical_granularity(C.num_edges());

    // [C91 §3.3 tex 276]: "γ-granularity, for any γ ≥ γ₂, can be
    // enforced ... in time linear in the size of the submap tree."
    // Preconditions hold: V(C) is conformal and 1-semigranular ⊆
    // γ-semigranular.
    enforce_granularity(S, C, gamma);

    // [C91 §3.3 tex 276]: "We can now put S in normal form, which
    // includes computing its tree decomposition."
    S.normalize(C);

    assert(S.is_conformal() && S.is_granular(gamma, C) &&
           !S.tree_decomposition().empty() &&
           "[C91 §4.1 tex 338]: canonical = 2^{⌈β⌈log m⌉⌉}-granular, "
           "conformal, normal form");
    return S;
}

} // namespace chazelle
