#include "trapezoids.h"
#include "../submap/boundary_geometry.h"
#include "../visibility/down_phase.h"

#include <cassert>
#include <utility>

namespace chazelle {

namespace {

struct InteriorRegion {
    std::size_t top = NONE;
    std::size_t bottom = NONE;
    std::size_t left_edge = NONE;
    std::size_t right_edge = NONE;
    std::size_t above = NONE;
    std::size_t below = NONE;
    bool visited = false;
};

bool edge_incident_to_vertex(std::size_t edge, std::size_t vertex, std::size_t vertex_count) {
    return edge == vertex || (edge + 1) % vertex_count == vertex;
}

}

TrapezoidDecomposition extract_trapezoids(const Polygon& curve, const Submap& visibility_map,
                                          std::span<const std::size_t> original_vertex_positions) {
    const std::size_t vertex_count = original_vertex_positions.size();
    assert(vertex_count >= 3 && original_vertex_positions.front() == 0 &&
           original_vertex_positions.back() + 1 == curve.num_vertices() &&
           "[FM84 Algorithm 1 input]: the vertices describe a simple polygon boundary");
    assert(visibility_map.start_vertex == 0 &&
           visibility_map.end_vertex + 1 == curve.num_vertices() && visibility_map.is_conformal() &&
           visibility_map.is_semigranular(1) &&
           "[C91 §2.1 tex 68–72]: extraction requires the full visibility map of C");

    TrapezoidDecomposition result;
    result.vertex_trapezoids.resize(vertex_count);
    result.trapezoids.reserve(2 * vertex_count);
    std::vector<std::size_t> original_vertices(curve.num_vertices(), NONE);
    std::vector<std::size_t> original_edges(curve.num_edges(), NONE);
    Exact twice_area = 0;
    for (std::size_t i = 0; i < vertex_count; ++i) {
        const std::size_t first = original_vertex_positions[i];
        const std::size_t next = original_vertex_positions[(i + 1) % vertex_count];
        original_vertices[first] = i;
        const Point& a = curve.vertex(first);
        const Point& b = curve.vertex(next);
        assert((a.x != b.x || a.y != b.y) &&
               "[FM84 Algorithm 1 input]: polygon edges have distinct endpoints");
        twice_area += a.x * b.y - a.y * b.x;
        if (i + 1 == vertex_count)
            continue;
        assert(first < next && "[C91 §4 tex 316]: subdivision preserves boundary order");
        for (std::size_t edge = first; edge < next; ++edge) {
            original_edges[edge] = i;
            ++result.work.boundary_edges;
            [[maybe_unused]] const Point& p = curve.vertex(edge);
            [[maybe_unused]] const Point& q = curve.vertex(edge + 1);
            assert((p.x != q.x || p.y != q.y) &&
                   (p.x - a.x) * (b.y - a.y) == (p.y - a.y) * (b.x - a.x) &&
                   (p.x - a.x) * (p.x - b.x) + (p.y - a.y) * (p.y - b.y) <= 0 &&
                   "[C91 §4 tex 316]: padding vertices subdivide their original edge");
        }
    }
    assert(twice_area != 0 && "[FM84 Algorithm 1 input]: a simple polygon has nonzero area");
    const Side inside = twice_area > 0 ? LEFT : RIGHT;
    const std::size_t closing_edge = vertex_count - 1;
    const SymbolicY closing_start = symbolic_y_of(curve.vertex(curve.num_vertices() - 1));
    const SymbolicY closing_end = symbolic_y_of(curve.vertex(0));
    const SymbolicY closing_low =
        symbolic_y_less(closing_start, closing_end) ? closing_start : closing_end;
    const SymbolicY closing_high =
        symbolic_y_less(closing_start, closing_end) ? closing_end : closing_start;

    std::vector<InteriorRegion> regions(visibility_map.num_nodes());
    for (std::size_t region = 0; region < visibility_map.num_nodes(); ++region) {
        if (visibility_map.node(region).dead)
            continue;
        ++result.work.regions;
        InteriorRegion& cell = regions[region];
        for (const std::size_t chord_index : visibility_map.node(region).incident_chords) {
            const Chord& chord = visibility_map.chord(chord_index);
            assert(!chord.dead);
            const std::size_t vertex = curve.local_index_of_tag(chord.y_tag);
            assert(vertex != NONE &&
                   symbolic_y_equal(chord.symbolic_y(), symbolic_y_of(curve.vertex(vertex))) &&
                   "[C91 §2.1 tex 68]: every visibility chord originates at an input vertex");
            if (cell.top == NONE || point_y_above(curve.vertex(vertex), curve.vertex(cell.top)))
                cell.top = vertex;
            if (cell.bottom == NONE ||
                point_y_below(curve.vertex(vertex), curve.vertex(cell.bottom)))
                cell.bottom = vertex;
        }
        assert(cell.top != NONE && cell.bottom != NONE);
        if (cell.top == cell.bottom)
            continue;
        const SymbolicY low = symbolic_y_of(curve.vertex(cell.bottom));
        const SymbolicY high = symbolic_y_of(curve.vertex(cell.top));
        for (const std::size_t arc_index : collect_region_arcs(visibility_map, region)) {
            ++result.work.arcs;
            const Arc& arc = visibility_map.arc(arc_index);
            if (arc.edge_count == 0)
                continue;
            const SymbolicY start = visibility_map.arc_start_symbolic_y(arc_index, curve);
            const SymbolicY end = visibility_map.arc_end_symbolic_y(arc_index, curve);
            const SymbolicY arc_low = symbolic_y_less(start, end) ? start : end;
            const SymbolicY arc_high = symbolic_y_less(start, end) ? end : start;
            ArcSideRange ranges[3];
            const std::size_t count = arc.side_ranges(0, curve.num_vertices() - 1, ranges);
            [[maybe_unused]] std::size_t range_edges = 0;
            for (std::size_t range = 0; range < count; ++range)
                range_edges += ranges[range].last_edge - ranges[range].first_edge + 1;
            assert(range_edges <= 3 &&
                   "[C91 §2.1 tex 68–72]: a full-map arc spans one edge and endpoint contacts");
            for (std::size_t range = 0; range < count; ++range) {
                if (ranges[range].side != inside)
                    continue;
                for (std::size_t edge = ranges[range].first_edge; edge <= ranges[range].last_edge;
                     ++edge) {
                    ++result.work.arc_edges;
                    const Point& a = curve.vertex(edge);
                    const Point& b = curve.vertex(edge + 1);
                    const bool ascending = point_y_below(a, b);
                    const SymbolicY edge_low = symbolic_y_of(ascending ? a : b);
                    const SymbolicY edge_high = symbolic_y_of(ascending ? b : a);
                    const SymbolicY overlap_low =
                        symbolic_y_less(arc_low, edge_low) ? edge_low : arc_low;
                    const SymbolicY overlap_high =
                        symbolic_y_less(arc_high, edge_high) ? arc_high : edge_high;
                    if (!symbolic_y_less(overlap_low, overlap_high))
                        continue;
                    assert(symbolic_y_equal(overlap_low, low) &&
                           symbolic_y_equal(overlap_high, high) &&
                           "[C91 §2.1 tex 68]: a full-map side spans the trapezoid's height");
                    std::size_t& side_edge =
                        ((inside == LEFT) == ascending) ? cell.right_edge : cell.left_edge;
                    assert(
                        side_edge == NONE &&
                        "[C91 §2.1 tex 68]: a trapezoid has one edge on each nonhorizontal side");
                    side_edge = original_edges[edge];
                }
            }
        }
        if (cell.left_edge == NONE && cell.right_edge == NONE)
            continue;
        if (cell.left_edge == NONE || cell.right_edge == NONE) {
            assert(symbolic_y_leq(closing_low, low) && symbolic_y_leq(high, closing_high) &&
                   "[FM84 Algorithm 4a]: the closing edge bounds every cut interior region");
            (cell.left_edge == NONE ? cell.left_edge : cell.right_edge) = closing_edge;
        }
        assert(cell.left_edge != cell.right_edge &&
               "[FM84 Algorithm 1 output]: an interior trapezoid has two distinct side edges");
    }

    for (std::size_t chord_index = 0; chord_index < visibility_map.num_chords(); ++chord_index) {
        const Chord& chord = visibility_map.chord(chord_index);
        if (chord.dead || chord.is_null_length)
            continue;
        ++result.work.chords;
        InteriorRegion& first = regions[chord.region[0]];
        InteriorRegion& second = regions[chord.region[1]];
        if (first.left_edge == NONE || second.left_edge == NONE ||
            first.left_edge != second.left_edge || first.right_edge != second.right_edge)
            continue;
        const bool first_above =
            symbolic_y_equal(symbolic_y_of(curve.vertex(first.bottom)), chord.symbolic_y());
        InteriorRegion& above = first_above ? first : second;
        InteriorRegion& below = first_above ? second : first;
        assert(
            symbolic_y_equal(symbolic_y_of(curve.vertex(above.bottom)), chord.symbolic_y()) &&
            symbolic_y_equal(symbolic_y_of(curve.vertex(below.top)), chord.symbolic_y()) &&
            above.below == NONE && below.above == NONE &&
            "[FM84 Algorithm 4a]: consecutive pieces with the same side edges form one trapezoid");
        above.below = chord.region[first_above ? 1 : 0];
        below.above = chord.region[first_above ? 0 : 1];
    }

    for (std::size_t region = 0; region < regions.size(); ++region) {
        const InteriorRegion& first = regions[region];
        if (first.left_edge == NONE || first.above != NONE)
            continue;
        std::size_t last = region;
        for (;;) {
            InteriorRegion& cell = regions[last];
            assert(!cell.visited && "[FM84 Theorem 5]: merging trapezoid pieces visits each once");
            cell.visited = true;
            ++result.work.joined_regions;
            if (cell.below == NONE)
                break;
            last = cell.below;
        }
        const std::size_t top = original_vertices[first.top];
        const std::size_t bottom = original_vertices[regions[last].bottom];
        assert(top != NONE && bottom != NONE &&
               "[FM84 §2 tex 144–148]: trapezoids begin and end at original polygon vertices");
        const std::size_t trapezoid_index = result.trapezoids.size();
        result.trapezoids.push_back({top, bottom, first.left_edge, first.right_edge});
        VertexTrapezoids& vertex = result.vertex_trapezoids[top];
        assert(vertex.count < vertex.trapezoids.size() &&
               "[FM84 Algorithm 1 output]: a vertex points to at most two trapezoids");
        vertex.trapezoids[vertex.count++] = trapezoid_index;
        if (vertex.count == 2) {
            const Trapezoid& a = result.trapezoids[vertex.trapezoids[0]];
            if (!edge_incident_to_vertex(a.right_edge, top, vertex_count))
                std::swap(vertex.trapezoids[0], vertex.trapezoids[1]);
            assert(edge_incident_to_vertex(result.trapezoids[vertex.trapezoids[0]].right_edge, top,
                                           vertex_count) &&
                   edge_incident_to_vertex(result.trapezoids[vertex.trapezoids[1]].left_edge, top,
                                           vertex_count) &&
                   "[FM84 Algorithm 1 output]: two trapezoids are stored leftmost first");
        }
    }
    for ([[maybe_unused]] const InteriorRegion& cell : regions)
        assert((cell.left_edge == NONE || cell.visited) &&
               "[FM84 Theorem 5]: all interior regions belong to an output trapezoid");
    return result;
}

TrapezoidDecomposition compute_trapezoid_decomposition(const std::vector<Point>& vertices) {
    assert(vertices.size() >= 3 &&
           "[FM84 Algorithm 1 input]: a polygon needs at least three vertices");
    const UpPhase up(vertices);
    const DownPhaseResult down = compute_visibility_map(up);
    return extract_trapezoids(up.graded().curve(), down.visibility_map,
                              up.graded().original_vertex_positions());
}

}
