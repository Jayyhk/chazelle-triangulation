// src/merge/merge.cpp

#include "merge.h"
#include "conformality.h"
#include "granularity.h"

namespace chazelle {

// [C91 §3]: Build C = C₁ ∪ C₂.  C₁ and C₂ share exactly one vertex
// (C₁'s last = C₂'s first); the merged curve drops the duplicate.
static Polygon build_merged_curve(const Polygon& C1, const Polygon& C2) {
    std::vector<Point> vertices;
    vertices.reserve(C1.num_vertices() + C2.num_vertices() - 1);

    for (std::size_t i = 0; i < C1.num_vertices(); ++i)
        vertices.push_back(C1.vertex(i));

    // Skip C₂'s first vertex (= C₁'s last vertex).
    for (std::size_t i = 1; i < C2.num_vertices(); ++i)
        vertices.push_back(C2.vertex(i));

    return Polygon(std::move(vertices));
}

// [C91 §3.1]: Stage 1 — fuse S₁ and S₂ via discovered chords from the
// fusion vertex traversal.
//
// TODO: wire up the three calls
//   1. st1.junction_at_end = true;
//      fuse_submaps(st1, S₁, C₁, S₂, C₂, ray1, ray2)   — S₁ into S₂
//   2. st2.junction_at_end = false;                     — [C91 §3.1 tex 179]
//      fuse_submaps(st2, S₂, C₂, S₁, C₁, ray2, ray1)   — by symmetry
//   3. rebuild_submap(out_S, C, S₁, C₁, S₂, C₂, st1, st2)
//
// Remaining blocker:
//   (a) [C91 §3.4]: no real ray-shooting oracle exists yet — fuse_submaps
//       calls local_shoot, which trips [C91 §3.1 tex 181]'s "local shoot
//       must hit" assertion under the no-op oracle stubs used in tests.
static void fuse(Submap& /*S*/, const MergeInput& /*in*/,
                  const Polygon& /*C*/) {
}

// [C91 §3.2]: Stage 2 — reduce every region to at most four arcs by
// adding visibility chords (Lemmas 3.2–3.4; see conformality.h).  On the
// pre-fusion empty submap (stage 1 is still blocked on [C91 §3.4]'s
// oracle) this is a no-op — no arcs, no regions to cut.
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
    if (S.num_arcs() == 0 && S.num_nodes() == 0) return;  // fuse() TODO gate
    restore_conformality(S, C, o);
}

// [C91 §3.3 tex 274–280]: Stage 3 — enforce γ-granularity by removing
// every removable exit chord (granularity.h), then "put S in normal
// form, which includes computing its tree decomposition" (tex 276).
// On the pre-fusion empty submap (stage 1 is still blocked on
// [C91 §3.4]'s oracle) this is a no-op — no tree to contract.
static void maintain_granularity(Submap& S, const MergeInput& in,
                                 const Polygon& C) {
    if (S.num_arcs() == 0 && S.num_nodes() == 0) return;  // fuse() TODO gate

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
