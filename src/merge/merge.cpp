// src/merge/merge.cpp

#include "merge.h"
#include "conformality.h"
#include "fusion.h"
#include "granularity.h"

namespace chazelle {

// [C91 §3 tex 160]: Build C = C₁ ∪ C₂.  C₁ and C₂ share exactly one
// vertex (C₁'s last = C₂'s first) and are contiguous subchains of the
// one input table of P, so C is the O(1) union view — [C91 §2.4 tex
// 133]: the input table "is never to be modified or even copied", and
// Lemma 4.1's sublinear per-merge budget ([C91 §4.1 tex 347–350])
// forbids a Θ(n₁ + n₂) copy here.
static Polygon build_merged_curve(const Polygon& C1, const Polygon& C2) {
    return Polygon(C1, C2);
}

// [C91 §3.1]: Stage 1 — fuse S₁ and S₂ via discovered chords from the
// fusion vertex traversal.
//
// [C91 §3.1 tex 179]: "By symmetry, we may limit our discussion to the
// problem of fusing S₁ into S₂ ... The idea is then to repeat the work
// described below with respect to S₂ (i.e., fusing S₂ into S₁), and
// set up a new submap S based on the information collected."  Pass 2
// runs the same algorithm with the arguments swapped and the junction
// at the walked curve's FIRST vertex (fusion.h).  [C91 §3.1 tex 226]:
// rebuild_submap then sets S up in normal form from the chord
// inventory (without its tree decomposition — §3.2's job).
static void fuse(Submap& S, const MergeInput& in, const Polygon& C) {
    FusionState st1;
    st1.junction_at_end = true;
    fuse_submaps(st1, *in.S1, *in.C1, *in.S2, *in.C2,
                 *in.ray_shooter_1, *in.ray_shooter_2);

    FusionState st2;
    st2.junction_at_end = false;
    fuse_submaps(st2, *in.S2, *in.C2, *in.S1, *in.C1,
                 *in.ray_shooter_2, *in.ray_shooter_1);

    rebuild_submap(S, C, *in.S1, *in.C1, *in.S2, *in.C2, st1, st2);
}

// [C91 §3.2]: Stage 2 — reduce every region to at most four arcs by
// adding visibility chords (Lemmas 3.2–3.4; see conformality.h).
static void restore_conformality(Submap& S, const MergeInput& in,
                                        const Polygon& C) {
    ConformalityOracles o;
    o.S1 = in.S1;
    o.S2 = in.S2;
    o.C1 = in.C1;
    o.C2 = in.C2;
    o.ray1 = in.ray_shooter_1;
    o.ray2 = in.ray_shooter_2;
    o.cut1 = in.arc_cutter_1;
    o.cut2 = in.arc_cutter_2;
    o.g_gamma1 = in.g_gamma1;
    o.g_gamma2 = in.g_gamma2;
    o.h_gamma1 = in.h_gamma1;
    o.h_gamma2 = in.h_gamma2;
    restore_conformality(S, C, o);
}

// [C91 §3.3 tex 274–280]: Stage 3 — enforce γ-granularity by removing
// every removable exit chord (granularity.h), then "put S in normal
// form, which includes computing its tree decomposition" (tex 276).
static void maintain_granularity(Submap& S, const MergeInput& in,
                                 const Polygon& C) {
    // [C91 §3.3 tex 276]: "γ-granularity, for any γ ≥ γ₂, can be
    // enforced in this nondeterministic fashion in time linear in the
    // size of the submap tree."
    enforce_granularity(S, C, in.gamma);

    // [C91 §3.3 tex 276]: "We can now put S in normal form."
    S.normalize(C);

    // [C91 §3.3 Lemma 3.5 tex 279]: the result is a normal-form
    // γ-granular conformal submap of V(C).
    assert(S.is_conformal() &&
           "[C91 Lemma 3.5 tex 279]: merge output must be conformal");
    assert(S.is_granular(in.gamma, C) &&
           "[C91 Lemma 3.5 tex 279]: merge output must be γ-granular");
    assert(!S.tree_decomposition().empty() &&
           "[C91 §2.4(iv) tex 139]: normal form includes the tree "
           "decomposition");
}

MergeResult merge(const MergeInput& in) {
    assert_merge_preconditions(in);

    MergeResult result(build_merged_curve(*in.C1, *in.C2));

    // [C91 §3]: "The merge proceeds in three stages."
    fuse(result.S, in, result.C);
    restore_conformality(result.S, in, result.C);
    maintain_granularity(result.S, in, result.C);

    return result;
}

} // namespace chazelle
