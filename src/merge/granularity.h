#pragma once

// [C91 §3.3 tex 274–280]: Maintaining Granularity.
//
// "Since by making S conformal we did not remove any exit chord, it is
// still the case that, as observed in the proof of Lemma 3.2, no arc
// has more than γ₂ edges.  Therefore, S is conformal and
// γ₂-semigranular.  We must now check whether the second criterion for
// γ₂ granularity holds.  This criterion says that contracting any edge
// of the submap tree that is incident upon at least one node of degree
// less than 3 produces a new node whose weight exceeds γ₂.  This is
// very easy to enforce: if an edge does not pass the test, we just
// contract it by removing its corresponding exit chord (and those
// endpoints that are not vertices of ∂C)." (tex 276)

#include "../polygon/polygon.h"
#include "../submap/submap.h"

#include <cstddef>

namespace chazelle {

// [C91 §3.3 tex 276]: enforce γ-granularity (γ ≥ γ₂) on the conformal,
// γ-semigranular fusion S by contracting every removable exit chord —
// an edge of the submap tree is removable iff it is incident upon a
// node of degree < 3 AND contracting it yields a merged node of weight
// ≤ γ ([C91 §2.3 tex 121] criterion (ii) negated).
//
// "We process each exit chord in turn and check whether it is
// removable" (tex 276).  The paper adds "chords need be processed only
// once since the removals cannot make any chord removable if it was
// not already so before"; that claim holds for the weight half of the
// test (weights only grow under contraction) but NOT for the degree
// half: contracting a leaf edge (degree-1 node) DECREASES the merged
// node's degree by one (d_u + d_v − 2 = d_u − 1), so a chord whose two
// endpoints both had degree 3 when it was processed can later become
// incident upon a degree-2 node and turn removable.  (Example: two
// adjacent degree-3 hubs, each with two leaf pockets of tiny weight,
// γ ≥ total weight — the hub–hub chord is skipped first, every pocket
// chord is then contracted, and both hubs end at degree 1 with the
// hub–hub chord removable.)  Since Lemma 3.5 (tex 279) asserts the
// output IS γ-granular, we re-examine the ≤ 4 chords incident upon the
// merged node after each removal.  Each removal enqueues O(1) rechecks
// and each check is O(1), so the total stays "linear in the size of
// the submap tree" (tex 276).
//
// Weights are maintained in a local per-region table seeded by one
// O(M) sweep of the arc-sequence table, updated per removal with the
// contraction simulation's result — region weights are NOT recomputed
// from chord adjacency mid-cascade.
//
// Preconditions (asserted): S conformal ([C91 §2.3 tex 114]; §3.2 has
// run) and γ-semigranular (tex 276).
// Postconditions (asserted): S conformal ("the removal keeps the
// submap conformal", tex 276), γ-semigranular ("this will not cause a
// violation of the first criterion", tex 276), and criterion (ii)
// holds for every remaining chord ([C91 §2.3 tex 121]).
//
// S is NOT left in normal form — Submap::normalize owns tex 276's "We
// can now put S in normal form" step.
//
// Time: O(M) where M = size of the submap tree (tex 276).
void enforce_granularity(Submap& S, const Polygon& C, std::size_t gamma);

} // namespace chazelle
