/// src/merge/merge.cpp

#include "merge.h"

namespace chazelle {

/// [C91 §3]: Build the merged curve C = C₁ ∪ C₂.
///
/// C₁ and C₂ share exactly one vertex (C₁'s last = C₂'s first).
/// The merged curve concatenates their vertices, deduplicating the
/// shared vertex.
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

/// [C91 §3.1]: "First we find which points of ∂C can be seen by the
/// endpoints of the exit chords of Sᵢ (i=1,2) and by the companion
/// vertices resulting from the duplication of C₁ ∩ C₂; this gives us
/// chords which we use to create a submap S of V(C) called the
/// fusion of S₁ and S₂."
///
/// TODO: (§3.1) implement.
static void fuse(Submap& /*S*/, const MergeInput& /*in*/,
                  const Polygon& /*C*/) {
}

/// [C91 §3.2]: "In the second stage we ensure that the submap is
/// conformal, which might involve adding new chords to cut up regions
/// with more than four arcs.  This is done by calling upon the
/// arc-cutting oracle... Finding new chords to cut up big regions is
/// carried out by binary search through the appropriate tree
/// decompositions, using the ray-shooting oracles along the way."
///
/// TODO: (§3.2) implement.
static void restore_conformality(Submap& /*S*/, const MergeInput& /*in*/,
                                  const Polygon& /*C*/) {
}

/// [C91 §3.3]: "In the third stage, finally, we bring the submap S
/// to the desired granularity by removing chords if necessary."
/// Then: "We can now put S in normal form."
///
/// TODO: (§3.3) implement.
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
