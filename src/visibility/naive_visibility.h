#pragma once

// [C91 §4.1 tex 345]: "This work can be done naively for the early
// grades" — the base of the up-phase.  The full visibility map V(C) is
// computed by brute force straight from its definition ([C91 §2.1 tex
// 68–72]: "extending two horizontal segments from each vertex v_i, one
// in each direction.  Each segment, which we call a chord, is extended
// until it meets another point of C.  If it were to go to infinity in
// the Cartesian plane, then it would actually wrap around in the
// spherical plane until it hits C again"), then brought to the
// canonical granularity of [C91 §4.1 tex 338] by the §3.3 machinery.
//
// Costs are polynomial in the curve size (O(m²) chord shots) — used
// only on the O(1)-size chains of the early grades, so the up-phase's
// linear budget ([C91 §4.1 tex 357]) is unaffected — and on test
// reference paths.

#include "../polygon/polygon.h"
#include "../submap/submap.h"
#include "../merge/oracle.h"

#include <cstddef>

namespace chazelle {

// [C91 §2.1 tex 70]: first contact of the horizontal ray from p (whose
// perturbed level is `sy`; p.x is the source x) traveling `dir`, over
// the WHOLE curve C, in the lexicographic (wrapped, distance) metric;
// coincident raw positions resolve per [C91 §2 tex 47]
// (polygon.h::ray_contact_precedes).  O(m) per call — the naive
// reference shooter.
RayHit naive_first_contact(const Polygon& C, Point p, SymbolicY sy,
                           Side dir, std::size_t source_edge = NONE);
// (source_edge is C-frame; its perturbed x-offset is derived
// internally.)

// [C91 §2.1 tex 68–72]: the full visibility map V(C), in canonical
// arc-sequence order with its tree decomposition.  Every vertex of ∂C
// is incident upon exactly one chord (tex 72); extremum inside pairs
// contribute their null-length chords.  O(m²).
Submap build_full_visibility_map(const Polygon& C);

// [C91 §4.1 tex 338]: the canonical granularity 2^{⌈β⌈log m⌉⌉} of a
// curve with m = `num_edges` contiguous edges of P.  For a chain in
// grade λ (m = 2^λ) this is 2^{⌈βλ⌉}.
std::size_t canonical_granularity(std::size_t num_edges);

// [C91 §4.1 tex 338/345]: a canonical submap of V(C) — "2^{⌈β⌈log m⌉⌉}-
// granular, conformal, and represented in normal form" — computed
// naively: full V(C), γ-granularity enforcement ([C91 §3.3 tex 276]),
// normal form.  O(m²).
Submap build_canonical_submap_naive(const Polygon& C);

} // namespace chazelle
