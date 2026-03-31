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

MergeResult merge(const MergeInput& in) {
    assert_merge_preconditions(in);

    MergeResult result(build_merged_curve(*in.C1, *in.C2));

    // TODO(§3.1): Fusion — discover chords, build submap S of V(C).
    // TODO(§3.2): Conformality — add chords until degree ≤ 4.
    // TODO(§3.3): Granularity — remove chords to enforce γ-granularity,
    //             then put S in normal form.

    return result;
}

} // namespace chazelle
