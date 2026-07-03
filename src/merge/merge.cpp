// src/merge/merge.cpp

#include "merge.h"

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

// [C91 §3.2]: Stage 2 — restore conformality by cutting up >4-arc
// regions using the arc-cutting oracle (binary search through tree
// decompositions + ray-shooters).  TODO: implement.
static void restore_conformality(Submap& /*S*/, const MergeInput& /*in*/,
                                  const Polygon& /*C*/) {
}

// [C91 §3.3]: Stage 3 — bring S to the desired granularity by removing
// chords, then put S in normal form.  TODO: implement.
static void enforce_granularity(Submap& /*S*/, const MergeInput& /*in*/,
                                 const Polygon& /*C*/) {
}

MergeResult merge(const MergeInput& in) {
    assert_merge_preconditions(in);

    MergeResult result(build_merged_curve(*in.C1, *in.C2));

    // [C91 §3]: "The merge proceeds in three stages."
    fuse(result.S, in, result.C);
    restore_conformality(result.S, in, result.C);
    enforce_granularity(result.S, in, result.C);

    return result;
}

} // namespace chazelle
