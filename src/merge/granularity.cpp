// src/merge/granularity.cpp — [C91 §3.3 tex 274–280]: Maintaining
// Granularity.

#include "granularity.h"

#include <algorithm>
#include <cassert>
#include <vector>

namespace chazelle {

void enforce_granularity(Submap& S, const Polygon& C, std::size_t gamma) {
    assert(gamma >= 1 && "[C91 §2.3]: granularity parameter must be positive");

    // [C91 §3.3 tex 276]: "S is conformal and γ₂-semigranular" on entry
    // (§3.2 restored conformality without removing exit chords, so no
    // arc exceeds γ₂ ≤ γ edges).
    assert(S.is_conformal() &&
           "[C91 §3.3 tex 276]: S must be conformal after §3.2");
    S.assert_tree_property();

    // ── Per-region weights, one O(M) sweep ([C91 §2.2 tex 106]:
    // weight = max nonnull-edge count over the region's arcs; 0 for
    // empty regions).  Maintained locally across the removal cascade —
    // the chord-adjacency walk of Submap::region_weight is O(1) per
    // query but this table makes each update O(1) too.
    std::vector<std::size_t> weight(S.num_nodes(), 0);
    for (std::size_t ai = 0; ai < S.num_arcs(); ++ai) {
        const Arc& a = S.arc(ai);
        if (a.dead) continue;
        assert(a.region_node < weight.size());
        weight[a.region_node] = std::max(weight[a.region_node], a.edge_count);
    }

#ifndef NDEBUG
    // [C91 §3.3 tex 276]: γ₂-semigranular and γ ≥ γ₂ ⟹ γ-semigranular.
    for (std::size_t r = 0; r < S.num_nodes(); ++r) {
        if (S.node(r).dead) continue;
        assert(weight[r] <= gamma &&
               "[C91 §3.3 tex 276]: S must be γ-semigranular on entry");
    }
#endif

    // [C91 §2.3 tex 121] criterion (ii), negated: removable ⟺ incident
    // upon a node of degree < 3 AND merged weight ≤ γ.
    auto contraction_weight = [&](std::size_t ci) -> std::size_t {
        const Chord& c = S.chord(ci);
        return S.simulated_contraction_weight(ci, C,
                                              weight[c.region[0]],
                                              weight[c.region[1]]);
    };
    auto degree_condition = [&](std::size_t ci) -> bool {
        const Chord& c = S.chord(ci);
        return S.node(c.region[0]).degree() < 3 ||
               S.node(c.region[1]).degree() < 3;
    };

    // ── [C91 §3.3 tex 276]: "We process each exit chord in turn and
    // check whether it is removable."  Worklist: the initial pass in
    // table order, plus O(1) re-checks of the merged node's incident
    // chords after each removal (see granularity.h on why the paper's
    // once-only claim needs the degree-drop re-check).
    std::vector<std::size_t> work;
    work.reserve(S.num_chords());
    for (std::size_t ci = S.num_chords(); ci-- > 0; )
        work.push_back(ci);        // pop order = table order

    std::size_t removals = 0;
    while (!work.empty()) {
        std::size_t ci = work.back();
        work.pop_back();
        if (S.chord(ci).dead) continue;
        if (!degree_condition(ci)) continue;

        std::size_t merged = contraction_weight(ci);
        if (merged > gamma) continue;   // passes criterion (ii): keep

        // [C91 §3.3 tex 276]: "we just contract it by removing its
        // corresponding exit chord (and those endpoints that are not
        // vertices of ∂C)."
        std::size_t survivor = S.remove_chord(ci, C);
        ++removals;
        assert(removals <= S.num_chords() && "each chord is removed once");

        // Merged node's weight is exactly the simulated contraction
        // weight (simulated_contraction_weight mirrors remove_chord's
        // glue semantics).  [C91 §3.3 tex 276]: "this will not cause a
        // violation of the first criterion" — merged ≤ γ by the test.
        assert(survivor < weight.size());
        weight[survivor] = merged;

        // [C91 §3.3 tex 276]: "the removal keeps the submap conformal"
        // — the merged degree is d_u + d_v − 2 ≤ 4 since min(d_u, d_v)
        // < 3 and max(d_u, d_v) ≤ 4.
        assert(S.node(survivor).degree() <= 4 &&
               "[C91 §3.3 tex 276]: contraction must preserve conformality");

        // Re-check the merged node's ≤ 4 incident chords: a leaf-edge
        // contraction lowers this node's degree, which can newly
        // satisfy the degree condition (weights only grow, so the
        // weight half of removability is monotone and needs no
        // re-check on its own).
        for (std::size_t nc : S.node(survivor).incident_chords) {
            assert(!S.chord(nc).dead &&
                   "incident_chords holds live chords only");
            work.push_back(nc);
        }
    }

#ifndef NDEBUG
    // ── Postconditions ([C91 §3.3 tex 276] / Lemma 3.5 tex 279) ─────
    assert(S.is_conformal() &&
           "[C91 §3.3 tex 276]: removals keep the submap conformal");
    S.assert_tree_property();

    // Local weight table matches the actual table state.
    {
        std::vector<std::size_t> check(S.num_nodes(), 0);
        for (std::size_t ai = 0; ai < S.num_arcs(); ++ai) {
            const Arc& a = S.arc(ai);
            if (a.dead) continue;
            check[a.region_node] =
                std::max(check[a.region_node], a.edge_count);
        }
        for (std::size_t r = 0; r < S.num_nodes(); ++r) {
            if (S.node(r).dead) continue;
            assert(check[r] == weight[r] &&
                   "[C91 §3.3]: maintained weights must match the table");
            // [C91 §2.3 tex 120] criterion (i).
            assert(check[r] <= gamma &&
                   "[C91 §3.3 tex 276]: S stays γ-semigranular");
        }
    }

    // [C91 §2.3 tex 121] criterion (ii): no chord is still removable.
    for (std::size_t ci = 0; ci < S.num_chords(); ++ci) {
        if (S.chord(ci).dead) continue;
        if (!degree_condition(ci)) continue;
        assert(contraction_weight(ci) > gamma &&
               "[C91 §2.3 tex 121]: contracting any edge incident upon a "
               "degree-<3 node must exceed γ (Lemma 3.5 postcondition)");
    }
#endif
}

} // namespace chazelle
