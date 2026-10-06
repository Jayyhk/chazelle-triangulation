#include "region_boundary.h"
#include "../submap/boundary_geometry.h"
#include "../submap/chord_inventory.h"

#include <algorithm>
#include <bit>
#include <cassert>
#include <utility>

namespace chazelle::animation {

namespace {

Exact absolute(const Exact& value) {
    return value < 0 ? -value : value;
}

Point primary_vertex(const Point& vertex) {
    return {vertex.x, vertex.y - Exact{vertex.index} * Exact::infinitesimal(3), vertex.index};
}

Polygon concatenate(const std::vector<RegionBoundaryPiece>& pieces) {
    assert(!pieces.empty());
    Polygon result = pieces[0].curve;
    for (std::size_t i = 1; i < pieces.size(); ++i)
        result = Polygon(result, pieces[i].curve);
    return result;
}

Submap transport(const Polygon& old_curve, const Submap& old_submap, const Polygon& new_curve,
                 bool reversed) {
    std::vector<PendingChord> chords;
    chords.reserve(old_submap.num_chords());
    for (std::size_t i = 0; i < old_submap.num_chords(); ++i) {
        const Chord& old = old_submap.chord(i);
        const std::size_t source = old_curve.local_index_of_tag(old.y_tag);
        assert(source != NONE && "[C91 §4.2 tex 367]: canonical pieces retain their input levels");
        auto edge = [&](std::size_t index) {
            return reversed ? old_curve.num_edges() - 1 - index : index;
        };
        auto side = [&](Side value) { return reversed ? (value == LEFT ? RIGHT : LEFT) : value; };
        const std::size_t new_source = reversed ? old_curve.num_vertices() - 1 - source : source;
        chords.push_back({symbolic_y_of(new_curve.vertex(new_source)), edge(old.left_edge),
                          side(old.left_side), edge(old.right_edge), side(old.right_side),
                          old.is_null_length});
    }
    Submap result;
    build_submap_from_chords(result, new_curve, chords);
    result.build_tree_decomposition();
    return result;
}

}

RegionBoundary::RegionBoundary(std::vector<RegionBoundaryPiece> boundary_pieces, std::size_t offset)
    : pieces(std::move(boundary_pieces)), curve(concatenate(pieces)), original_offset(offset) {
    assert(curve.vertex(0).index != curve.vertex(curve.num_vertices() - 1).index &&
           "[C91 §4.2 tex 367]: R* is punctured to make it nonclosed");
}

OriginalBoundaryLocation RegionBoundary::original_location(std::size_t edge, Side side) const {
    if (side != LEFT)
        return {};
    for (const RegionBoundaryPiece& piece : pieces) {
        if (edge < piece.curve.num_edges()) {
            if (piece.first_original_vertex == NONE)
                return piece.edge_locations[edge];
            return {(piece.reversed ? piece.first_original_vertex - 1 - edge
                                    : piece.first_original_vertex + edge) -
                        original_offset,
                    piece.original_side, piece.original_arc};
        }
        edge -= piece.curve.num_edges();
    }
    assert(false);
    return {};
}

SymbolicY RegionBoundary::original_level(std::size_t tag, const Polygon& original) const {
    for (const RegionBoundaryPiece& piece : pieces) {
        const std::size_t vertex = piece.curve.local_index_of_tag(tag);
        if (vertex == NONE)
            continue;
        if (piece.first_original_vertex == NONE)
            return piece.vertex_levels[vertex];
        const std::size_t index = piece.reversed ? piece.first_original_vertex - vertex
                                                 : piece.first_original_vertex + vertex;
        return symbolic_y_of(original.vertex(index));
    }
    assert(false && "[C91 §4.2 tex 377]: every extracted chord retains its original level");
    return {};
}

RegionBoundaryGeometry::RegionBoundaryGeometry(const UpPhase& up_phase, const Polygon& original)
    : up_phase_(&up_phase), first_original_vertex_(original.table_offset()),
      original_edge_count_(original.num_edges()),
      next_tag_(2 * up_phase.graded().curve().num_vertices()) {
    assert(std::has_single_bit(original_edge_count_) &&
           first_original_vertex_ % original_edge_count_ == 0 &&
           "[C91 Lemma 4.2 tex 364]: C is an aligned graded chain");
    const Exact width = Exact::infinitesimal(2);
    const Exact tilt = Exact::infinitesimal(1);
    std::array<std::vector<Point>, 2> sides;
    for (auto& vertices : sides)
        vertices.reserve(original.num_vertices());
    for (std::size_t i = 0; i < original.num_vertices(); ++i) {
        const Point vertex = primary_vertex(original.vertex(i));
        const Point& before = original.vertex(i == 0 ? 0 : i - 1);
        const Point& after = original.vertex(i + 1 == original.num_vertices() ? i : i + 1);
        const Point& unperturbed = original.vertex(i);
        Exact dx = i == 0 ? after.x - unperturbed.x : unperturbed.x - before.x;
        Exact dy = i == 0 ? after.y - unperturbed.y : unperturbed.y - before.y;
        const Exact length = absolute(dx) + absolute(dy);
        const Exact squared_length = dx * dx + dy * dy;
        assert(squared_length > 0 &&
               "[C91 §2.1 tex 68]: the curve has distinct consecutive vertices");
        Exact offset_x = -dy * length / squared_length;
        Exact offset_y = dx * length / squared_length;
        if (i > 0 && i + 1 < original.num_vertices()) {
            const Exact next_dx = after.x - unperturbed.x;
            const Exact next_dy = after.y - unperturbed.y;
            const Exact next_length = absolute(next_dx) + absolute(next_dy);
            const Exact determinant = dx * next_dy - dy * next_dx;
            if (determinant != 0) {
                offset_x = (length * next_dx - next_length * dx) / determinant;
                offset_y = (length * next_dy - next_length * dy) / determinant;
            } else {
                assert(dx * next_dx + dy * next_dy > 0 &&
                       "[C91 §2.1 tex 68]: consecutive collinear edges cannot overlap");
            }
        }
        Exact cap_x = 0;
        Exact cap_y = 0;
        if (i == 0 || i + 1 == original.num_vertices()) {
            const Exact sign = i == 0 ? -1 : 1;
            cap_x = sign * dx * length / squared_length;
            cap_y = sign * dy * length / squared_length;
        }
        for (const Side side : {LEFT, RIGHT}) {
            const Exact sign = side == LEFT ? 1 : -1;
            const std::size_t tag =
                2 * (first_original_vertex_ + i) + static_cast<std::size_t>(side);
            sides[static_cast<std::size_t>(side)].push_back(
                {vertex.x + width * (cap_x + sign * offset_x),
                 vertex.y + width * (cap_y + sign * offset_y) - Exact{tag} * tilt, tag});
        }
    }
    for (const Side side : {LEFT, RIGHT})
        sides_[static_cast<std::size_t>(side)] = std::make_unique<Polygon>(
            Polygon::symbolic_curve(std::move(sides[static_cast<std::size_t>(side)])));
}

CanonicalBoundaryChain& RegionBoundaryGeometry::canonical_chain(std::size_t grade,
                                                                std::size_t index, Side side,
                                                                bool reversed) const {
    auto& variant = chains_[2 * static_cast<std::size_t>(side) + (reversed ? 1 : 0)];
    if (variant.empty()) {
        variant.resize(static_cast<std::size_t>(std::bit_width(original_edge_count_)));
        for (std::size_t k = 0; k < variant.size(); ++k)
            variant[k].resize(original_edge_count_ >> k);
    }
    assert(grade < variant.size() && index >= (first_original_vertex_ >> grade) &&
           index - (first_original_vertex_ >> grade) < variant[grade].size() &&
           "[C91 Lemma 4.2 tex 364]: auxiliary chain preprocessing is confined to C");
    auto& cached = variant[grade][index - (first_original_vertex_ >> grade)];
    if (!cached) {
        const Polygon& original = up_phase_->graded().chain(grade, index);
        Polygon curve =
            (*sides_[static_cast<std::size_t>(side)])
                .subchain((index << grade) - first_original_vertex_, original.num_vertices(),
                          original.min_y_vertex(), original.max_y_vertex());
        if (reversed)
            curve = curve.reversed();
        Submap submap = transport(original, up_phase_->chain_submap(grade, index), curve, reversed);
        cached = std::make_unique<CanonicalBoundaryChain>(
            CanonicalBoundaryChain{std::move(curve), std::move(submap), nullptr});
    }
    return *cached;
}

Point RegionBoundaryGeometry::arc_endpoint(const Polygon& original, std::size_t global_edge,
                                           Side side, const SymbolicY& level, bool first) {
    const Point a = primary_vertex(original.vertex(global_edge));
    const Point b = primary_vertex(original.vertex(global_edge + 1));
    Exact position;
    if (level.tag == a.index && level.y == original.vertex(global_edge).y)
        position = 0;
    else if (level.tag == b.index && level.y == original.vertex(global_edge + 1).y)
        position = 1;
    else
        position = (level.y - Exact{level.tag} * Exact::infinitesimal(3) - a.y) / (b.y - a.y);
    const Exact inset = Exact::infinitesimal(0);
    position += ((side == LEFT) == first) ? inset : -inset;
    const Polygon& displaced = *sides_[static_cast<std::size_t>(side)];
    const Point& start = displaced.vertex(global_edge - first_original_vertex_);
    const Point& end = displaced.vertex(global_edge + 1 - first_original_vertex_);
    assert(next_tag_ < SOS_NONE);
    return {start.x + position * (end.x - start.x), start.y + position * (end.y - start.y),
            next_tag_++};
}

UpPhase::PortionResult RegionBoundaryGeometry::canonical_original_piece(std::size_t first,
                                                                        std::size_t last, Side side,
                                                                        bool reversed) const {
    assert(first_original_vertex_ <= first && first < last &&
           last - first_original_vertex_ <= original_edge_count_);
    UpPhase::PortionResult original =
        last - first <= 2
            ? UpPhase::PortionResult{up_phase_->graded().curve().subchain(first, last - first + 1),
                                     {}}
            : up_phase_->compute_canonical_portion(first, last);
    if (last - first <= 2)
        original.submap = build_canonical_submap_naive(original.curve);
    Polygon curve = (*sides_[static_cast<std::size_t>(side)])
                        .subchain(first - first_original_vertex_, last - first + 1,
                                  original.curve.min_y_vertex(), original.curve.max_y_vertex());
    if (reversed)
        curve = curve.reversed();
    Submap submap = transport(original.curve, original.submap, curve, reversed);
    return {std::move(curve), std::move(submap)};
}

UpPhase::PortionResult
RegionBoundaryGeometry::canonical_piece(const RegionBoundaryPiece& piece) const {
    if (piece.canonical_submap)
        return {piece.curve, *piece.canonical_submap};
    if (piece.first_original_vertex == NONE)
        return {piece.curve, build_canonical_submap_naive(piece.curve)};
    const std::size_t edges = piece.curve.num_edges();
    const std::size_t first =
        piece.reversed ? piece.first_original_vertex - edges : piece.first_original_vertex;
    const std::size_t last = first + edges;
    return canonical_original_piece(first, last, piece.original_side, piece.reversed);
}

RegionBoundary RegionBoundaryGeometry::boundary(const Polygon& curve, const Submap& submap,
                                                std::size_t region) {
    const Polygon& original = up_phase_->graded().curve();
    const std::size_t table_offset = curve.table_offset();
    RegionArcs arcs = collect_region_arcs(submap, region);
    for (std::size_t i = 1; i < arcs.count; ++i)
        for (std::size_t j = i; j > 0 && arcs.arcs[j] < arcs.arcs[j - 1]; --j)
            std::swap(arcs.arcs[j], arcs.arcs[j - 1]);
    assert(!submap.node(region).incident_chords.empty());
    std::vector<RegionBoundaryPiece> pieces;
    struct ArcEnd {
        Point point;
        SymbolicY level;
    };
    std::vector<ArcEnd> starts;
    std::vector<ArcEnd> ends;
    std::vector<std::vector<RegionBoundaryPiece>> arc_pieces;
    auto segment = [&](const Point& a, const Point& b, SymbolicY first_level, SymbolicY last_level,
                       OriginalBoundaryLocation location) {
        return RegionBoundaryPiece{Polygon::symbolic_curve({a, b}),
                                   NONE,
                                   false,
                                   LEFT,
                                   NONE,
                                   {location},
                                   {std::move(first_level), std::move(last_level)}};
    };
    for (std::size_t arc : arcs) {
        ArcSideRange ranges[3];
        const std::size_t count = submap.arc(arc).side_ranges(0, curve.num_vertices() - 1, ranges);
        std::vector<RegionBoundaryPiece> current;
        for (std::size_t r = 0; r < count; ++r) {
            const Side side = ranges[r].side;
            const bool reverse = side == RIGHT;
            std::size_t first_edge = reverse ? ranges[r].last_edge : ranges[r].first_edge;
            std::size_t last_edge = reverse ? ranges[r].first_edge : ranges[r].last_edge;
            const SymbolicY first_level =
                r == 0 ? submap.arc_start_symbolic_y(arc, curve)
                       : symbolic_y_of(curve.vertex(reverse ? first_edge + 1 : first_edge));
            const SymbolicY last_level =
                r + 1 == count ? submap.arc_end_symbolic_y(arc, curve)
                               : symbolic_y_of(curve.vertex(reverse ? last_edge : last_edge + 1));
            if (first_edge != last_edge &&
                symbolic_y_equal(first_level, symbolic_y_of(curve.vertex(
                                                  reverse ? first_edge : first_edge + 1))))
                first_edge = reverse ? first_edge - 1 : first_edge + 1;
            if (first_edge != last_edge &&
                symbolic_y_equal(last_level,
                                 symbolic_y_of(curve.vertex(reverse ? last_edge + 1 : last_edge))))
                last_edge = reverse ? last_edge + 1 : last_edge - 1;
            const Point first =
                arc_endpoint(original, table_offset + first_edge, side, first_level, true);
            const Point last =
                arc_endpoint(original, table_offset + last_edge, side, last_level, false);
            if (r == 0)
                starts.push_back({first, first_level});
            else {
                const Point& previous =
                    current.back().curve.vertex(current.back().curve.num_vertices() - 1);
                current.push_back(segment(previous, first, first_level, first_level, {}));
            }
            const OriginalBoundaryLocation first_location{first_edge, side, arc};
            if (first_edge == last_edge) {
                current.push_back(segment(first, last, first_level, last_level, first_location));
            } else {
                const std::size_t first_vertex =
                    table_offset + (reverse ? first_edge : first_edge + 1);
                const std::size_t last_vertex =
                    table_offset + (reverse ? last_edge + 1 : last_edge);
                const Polygon& boundary_side = *sides_[static_cast<std::size_t>(side)];
                current.push_back(segment(
                    first, boundary_side.vertex(first_vertex - first_original_vertex_), first_level,
                    symbolic_y_of(original.vertex(first_vertex)), first_location));
                if (first_vertex != last_vertex) {
                    auto canonical = canonical_original_piece(std::min(first_vertex, last_vertex),
                                                              std::max(first_vertex, last_vertex),
                                                              side, reverse);
                    current.push_back(
                        {std::move(canonical.curve), first_vertex, reverse, side, arc, {}, {}});
                    current.back().canonical_submap =
                        std::make_shared<Submap>(std::move(canonical.submap));
                }
                current.push_back(
                    segment(boundary_side.vertex(last_vertex - first_original_vertex_), last,
                            symbolic_y_of(original.vertex(last_vertex)), last_level,
                            {last_edge, side, arc}));
            }
            if (r + 1 == count)
                ends.push_back({last, last_level});
        }
        arc_pieces.push_back(std::move(current));
    }
    std::vector<RegionBoundaryPiece> connectors;
    for (std::size_t i = 0; i < arcs.count; ++i) {
        const std::size_t next = (i + 1) % arcs.count;
        const Chord* exit = nullptr;
        bool from_left = false;
        for (std::size_t chord_index : submap.node(region).incident_chords) {
            const Chord& chord = submap.chord(chord_index);
            if (!symbolic_y_equal(chord.symbolic_y(), ends[i].level) ||
                !symbolic_y_equal(chord.symbolic_y(), starts[next].level))
                continue;
            const Arc& before = submap.arc(arcs.arcs[i]);
            const Arc& after = submap.arc(arcs.arcs[next]);
            auto matches = [&](std::size_t edge, Side side, std::size_t slot_edge, Side slot_side) {
                if (edge == slot_edge && side == slot_side)
                    return true;
                const std::size_t vertex = curve.local_index_of_tag(chord.y_tag);
                if (vertex == NONE ||
                    !symbolic_y_equal(chord.symbolic_y(), symbolic_y_of(curve.vertex(vertex))))
                    return false;
                if (chord.is_null_length)
                    return is_inside_companion(curve, edge, side, vertex) &&
                           is_inside_companion(curve, slot_edge, slot_side, vertex);
                return side == slot_side && (edge == vertex || edge + 1 == vertex) &&
                       (slot_edge == vertex || slot_edge + 1 == vertex);
            };
            if (matches(before.last_edge, before.last_side, chord.left_edge, chord.left_side) &&
                matches(after.first_edge, after.first_side, chord.right_edge, chord.right_side)) {
                exit = &chord;
                from_left = true;
                break;
            }
            if (matches(before.last_edge, before.last_side, chord.right_edge, chord.right_side) &&
                matches(after.first_edge, after.first_side, chord.left_edge, chord.left_side)) {
                exit = &chord;
                break;
            }
        }
        assert(exit &&
               "[C91 §2.2 tex 96, §4.2 tex 370]: successive arcs are joined by an exit chord");
        const Exact period = 4 / Exact::infinitesimal(3);
        const Exact horizontal_shift =
            !exit->is_null_length && chord_runs_through_infinity(curve, *exit)
                ? (from_left ? -period : period)
                : Exact{0};
        std::vector<Point> points{ends[i].point, starts[next].point};
        std::vector<SymbolicY> levels(points.size(), exit->symbolic_y());
        connectors.push_back(
            {Polygon::symbolic_curve(std::move(points), horizontal_shift), NONE, false, LEFT, NONE,
             std::vector<OriginalBoundaryLocation>(levels.size() - 1), std::move(levels)});
    }
    RegionBoundaryPiece& cut = connectors.back();
    const Point& a = cut.curve.vertex(0);
    const Point& b = cut.curve.vertex(1);
    const Exact inset = Exact::infinitesimal(0);
    const Exact dx = cut.curve.edge_horizontal_delta(0);
    const Point first{a.x + 2 * inset * dx, a.y + 2 * inset * (b.y - a.y), next_tag_++};
    const Point last{a.x + inset * dx, a.y + inset * (b.y - a.y), next_tag_++};
    pieces.push_back({Polygon::symbolic_curve({first, b}, cut.curve.edge_horizontal_shift(0)),
                      NONE,
                      false,
                      LEFT,
                      NONE,
                      {{}},
                      {cut.vertex_levels[0], cut.vertex_levels[1]}});
    for (std::size_t i = 0; i < arcs.count; ++i) {
        for (auto& piece : arc_pieces[i])
            pieces.push_back(std::move(piece));
        if (i + 1 < arcs.count)
            pieces.push_back(std::move(connectors[i]));
    }
    pieces.push_back(segment(a, last, cut.vertex_levels[0], cut.vertex_levels[0], {}));
    return RegionBoundary(std::move(pieces), table_offset);
}

}
