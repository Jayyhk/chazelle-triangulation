#include "chazelle/up_phase.h"
#include "chazelle/merge.h"

#include <cassert>
#include <cmath>

namespace chazelle {

/// Build a trivial canonical submap for a single edge.
/// A single edge has a trivial visibility map: one region, no chords,
/// one arc (the edge itself on each side of the double boundary).
static CanonicalSubmap build_grade0_submap(
    const Polygon& polygon,
    std::size_t chain_index) {

    auto range = polygon.chain_range(0, chain_index);

    CanonicalSubmap cs;
    cs.grade = 0;
    cs.chain_index = chain_index;
    cs.start_vertex = range.start_vertex;
    cs.end_vertex = range.end_vertex;

    // The visibility map of a single edge has:
    //   - 1 region (the single trapezoid / degenerate region)
    //   - 0 chords
    //   - 2 arcs (left side + right side of the edge)
    std::size_t region = cs.submap.add_node();

    // Left-side arc.
    ArcStructure left_arc;
    left_arc.first_edge = range.start_vertex; // edge index = start vertex
    left_arc.last_edge  = range.start_vertex;
    left_arc.first_side = Side::LEFT;
    left_arc.last_side  = Side::LEFT;
    left_arc.region_node = region;
    left_arc.edge_count = 1;
    std::size_t lai = cs.submap.add_arc(left_arc);
    cs.submap.node(region).arcs.push_back(lai);

    // Right-side arc.
    ArcStructure right_arc;
    right_arc.first_edge = range.start_vertex;
    right_arc.last_edge  = range.start_vertex;
    right_arc.first_side = Side::RIGHT;
    right_arc.last_side  = Side::RIGHT;
    right_arc.region_node = region;
    right_arc.edge_count = 1;
    std::size_t rai = cs.submap.add_arc(right_arc);
    cs.submap.node(region).arcs.push_back(rai);

    // §2.3: The arc-sequence table records which arcs pass through
    // the endpoints of C.  For a grade-0 chain (single edge), both
    // endpoints are covered by both arcs.
    cs.submap.start_arc = lai;
    cs.submap.end_arc   = rai;

    // Store chain vertex range and polygon pointer for §2.2 / §2.4.
    cs.submap.set_chain_info(range.start_vertex, range.end_vertex,
                             &polygon);

    cs.submap.recompute_weight(region);

    // Build ray-shooting oracle (trivial for one region).
    cs.oracle.build(cs.submap, polygon, 1);

    return cs;
}

void up_phase(const Polygon& polygon, GradeStorage& storage) {
    std::size_t p = polygon.num_grades();

    storage.init(p + 1, polygon.num_vertices());

    // Grade 0: one chain per edge, trivial submaps.
    std::size_t num_chains_0 = polygon.num_edges(); // 2^p chains
    for (std::size_t j = 0; j < num_chains_0; ++j) {
        storage.store(0, j, build_grade0_submap(polygon, j));
    }

    // Grades 1 through p: merge pairs of chains from the previous grade.
    for (std::size_t lambda = 1; lambda <= p; ++lambda) {
        std::size_t num_chains_lambda = polygon.num_chains(lambda);
        for (std::size_t j = 0; j < num_chains_lambda; ++j) {
            auto merged = merge_submaps(lambda, j, storage, polygon);
            storage.store(lambda, j, std::move(merged));
        }
    }
}

} // namespace chazelle
