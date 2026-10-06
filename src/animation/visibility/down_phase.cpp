#include "down_phase.h"
#include "../merge/conformality.h"
#include "../merge/granularity.h"
#include "../submap/boundary_geometry.h"
#include "../submap/chord_inventory.h"
#include "../trace.h"
#include "region_oracles.h"

#include <algorithm>
#include <bit>
#include <cassert>
#include <memory>
#include <utility>
#include <vector>

namespace chazelle::animation {

namespace {

Subarc whole_arc(const Polygon& curve, const Submap& submap, std::size_t arc) {
    const Arc& a = submap.arc(arc);
    return {a.first_edge,
            a.first_side,
            a.last_edge,
            a.last_side,
            submap.arc_start_symbolic_y(arc, curve),
            submap.arc_end_symbolic_y(arc, curve)};
}

bool contains(const Polygon& curve, const Submap& submap, std::size_t arc,
              const ChordEndpoint& endpoint, bool include_companion = true) {
    const Subarc boundary = whole_arc(curve, submap, arc);
    if (subarc_contains_point(boundary, curve, endpoint.edge_c, endpoint.side, endpoint.y, 0,
                              curve.num_vertices() - 1))
        return true;
    if (!include_companion)
        return false;
    const std::size_t vertex = curve.local_index_of_tag(endpoint.y.tag);
    if (vertex == NONE || !symbolic_y_equal(endpoint.y, symbolic_y_of(curve.vertex(vertex))) ||
        (endpoint.edge_c != vertex && endpoint.edge_c + 1 != vertex))
        return false;
    auto same_companion = [&](std::size_t edge, Side side, const SymbolicY& level) {
        if (!symbolic_y_equal(level, endpoint.y) || (edge != vertex && edge + 1 != vertex))
            return false;
        if (!curve.is_y_extremum(vertex))
            return side == endpoint.side;
        return is_inside_companion(curve, edge, side, vertex) &&
               is_inside_companion(curve, endpoint.edge_c, endpoint.side, vertex);
    };
    return same_companion(boundary.first_edge, boundary.first_side, boundary.first_y) ||
           same_companion(boundary.last_edge, boundary.last_side, boundary.last_y);
}

std::size_t endpoint_owner(const Polygon& curve, const Submap& submap, std::size_t original_arc,
                           const ChordEndpoint& endpoint) {
    const std::size_t region = submap.arc(original_arc).region_node;
    for (bool include_companion : {false, true}) {
        if (contains(curve, submap, original_arc, endpoint, include_companion))
            return original_arc;
        for (std::size_t arc : collect_region_arcs(submap, region))
            if (contains(curve, submap, arc, endpoint, include_companion))
                return arc;
        for (std::size_t chord_index : submap.node(region).incident_chords) {
            const Chord& chord = submap.chord(chord_index);
            for (const Chord::AdjArcs* adjacent : {&chord.left_adj, &chord.right_adj})
                for (std::size_t i = 0; i < adjacent->count; ++i)
                    if (contains(curve, submap, adjacent->arcs[i], endpoint, include_companion))
                        return adjacent->arcs[i];
        }
    }
    assert(false && "[C91 §4.2 tex 377]: extracted endpoints lie on the previous boundary arcs");
    return NONE;
}

bool equal_chord(const PendingChord& first, const PendingChord& second) {
    return symbolic_y_equal(first.y, second.y) && first.left_edge_c == second.left_edge_c &&
           first.left_side == second.left_side && first.right_edge_c == second.right_edge_c &&
           first.right_side == second.right_side && first.is_null_length == second.is_null_length;
}

class RefinementInventory {
public:
    RefinementInventory(const Polygon& curve, const Submap& coarse)
        : curve_(curve), coarse_(coarse), endpoints_(coarse.num_arcs()) {}

    void add(PendingChord chord, std::size_t first_arc, std::size_t second_arc, bool inherited) {
        canonicalize_chord(chord, curve_);
        const std::size_t index = chords_.size();
        chords_.push_back(std::move(chord));
        const PendingChord& stored = chords_.back();
        ChordEndpoint first{stored.left_edge_c, stored.left_side, stored.y, index, true};
        ChordEndpoint second{stored.right_edge_c, stored.right_side, stored.y, index, false};
        if (!inherited) {
            const bool first_matches = contains(curve_, coarse_, first_arc, first);
            const bool second_matches = contains(curve_, coarse_, second_arc, second);
            if (!first_matches && !second_matches && contains(curve_, coarse_, second_arc, first) &&
                contains(curve_, coarse_, first_arc, second))
                std::swap(first_arc, second_arc);
            first_arc = endpoint_owner(curve_, coarse_, first_arc, first);
            second_arc = endpoint_owner(curve_, coarse_, second_arc, second);
        }
        endpoints_[first_arc].push_back(first);
        endpoints_[second_arc].push_back(second);
    }

    Submap build() {
        auto precedes = [&](const auto& a, const auto& b) {
            return chord_endpoint_precedes(curve_, chords_, a, b);
        };
        for (auto& list : endpoints_)
            std::sort(list.begin(), list.end(), precedes);
        const Arc& start = coarse_.arc(coarse_.start_arc);
        const SymbolicY start_y = coarse_.arc_start_symbolic_y(coarse_.start_arc, curve_);
        auto position = [&](std::size_t edge, Side side) {
            return side == LEFT ? edge : 2 * curve_.num_edges() - 1 - edge;
        };
        const std::size_t start_position = position(start.first_edge, start.first_side);
        auto head = [&](const ChordEndpoint& endpoint) {
            const std::size_t key = position(endpoint.edge_c, endpoint.side);
            if (key != start_position)
                return key < start_position;
            const bool increasing =
                (endpoint.side == LEFT) ==
                point_y_below(curve_.vertex(endpoint.edge_c), curve_.vertex(endpoint.edge_c + 1));
            return increasing ? symbolic_y_less(endpoint.y, start_y)
                              : symbolic_y_greater(endpoint.y, start_y);
        };
        std::vector<ChordEndpoint> ordered;
        ordered.reserve(2 * chords_.size());
        auto append = [&](std::size_t arc, bool beginning) {
            for (const auto& endpoint : endpoints_[arc])
                if (arc != coarse_.start_arc || head(endpoint) == beginning)
                    ordered.push_back(endpoint);
        };
        append(coarse_.start_arc, true);
        for (std::size_t arc = 0; arc < coarse_.num_arcs(); ++arc)
            append(arc, false);
        assert(ordered.size() == 2 * chords_.size());
        std::vector<bool> duplicate(chords_.size(), false);
        for (std::size_t first = 0; first < ordered.size();) {
            std::size_t last = first + 1;
            while (last < ordered.size() && ordered[last].edge_c == ordered[first].edge_c &&
                   ordered[last].side == ordered[first].side &&
                   symbolic_y_equal(ordered[last].y, ordered[first].y))
                ++last;
            assert(
                last - first <= 128 &&
                "[C91 §4.2 tex 372]: a boundary junction creates only a constant number of incidences");
            std::sort(ordered.begin() + static_cast<std::ptrdiff_t>(first),
                      ordered.begin() + static_cast<std::ptrdiff_t>(last), precedes);
            for (std::size_t i = first; i < last; ++i) {
                const ChordEndpoint& endpoint = ordered[i];
                if (!endpoint.is_left_slot)
                    continue;
                for (std::size_t j = first; j < i; ++j) {
                    const ChordEndpoint& other = ordered[j];
                    if (other.is_left_slot &&
                        equal_chord(chords_[other.pending_idx], chords_[endpoint.pending_idx]))
                        duplicate[std::max(other.pending_idx, endpoint.pending_idx)] = true;
                }
            }
            first = last;
        }
        std::vector<std::size_t> remap(chords_.size(), NONE);
        std::vector<PendingChord> unique;
        unique.reserve(chords_.size());
        for (std::size_t i = 0; i < chords_.size(); ++i) {
            if (duplicate[i])
                continue;
            remap[i] = unique.size();
            unique.push_back(std::move(chords_[i]));
        }
        std::size_t retained = 0;
        for (ChordEndpoint endpoint : ordered) {
            if (duplicate[endpoint.pending_idx])
                continue;
            endpoint.pending_idx = remap[endpoint.pending_idx];
            ordered[retained++] = std::move(endpoint);
        }
        ordered.resize(retained);
        Submap result;
        build_submap_from_ordered_chords(result, curve_, unique, ordered);
        return result;
    }

private:
    const Polygon& curve_;
    const Submap& coarse_;
    std::vector<PendingChord> chords_;
    std::vector<std::vector<ChordEndpoint>> endpoints_;
};

void inherit(RefinementInventory& inventory, const Submap& original,
             const Submap* retained = nullptr) {
    for (std::size_t i = 0; i < original.num_chords(); ++i) {
        if (retained && retained->chord(i).dead)
            continue;
        const Chord& chord = original.chord(i);
        inventory.add({chord.symbolic_y(), chord.left_edge, chord.left_side, chord.right_edge,
                       chord.right_side, chord.is_null_length},
                      chord.left_adj.arcs[0], chord.right_adj.arcs[0], true);
    }
}

Submap restore_granularity(Submap submap, const Polygon& curve, std::size_t granularity) {
    Submap before = submap;
    if (auto* trace = AnimationTrace::current())
        trace->copy_submap(submap, before, curve);
    enforce_granularity(submap, curve, granularity);
    RefinementInventory inventory(curve, before);
    inherit(inventory, before, &submap);
    Submap result = inventory.build();
    result.build_tree_decomposition();
    assert(result.is_conformal() && result.is_granular(granularity, curve) &&
           "[C91 §3.3 tex 276, §4.2 tex 377]: removal restores granularity in linear work");
    return result;
}

Submap extract_visibility_map(const Polygon& curve, const Submap& augmented) {
    RefinementInventory inventory(curve, augmented);
    for (std::size_t i = 0; i < augmented.num_chords(); ++i) {
        const Chord& chord = augmented.chord(i);
        const std::size_t vertex = curve.local_index_of_tag(chord.y_tag);
        if (vertex == NONE ||
            !symbolic_y_equal(chord.symbolic_y(), symbolic_y_of(curve.vertex(vertex))))
            continue;
        if (chord.left_edge != vertex && chord.left_edge + 1 != vertex &&
            chord.right_edge != vertex && chord.right_edge + 1 != vertex)
            continue;
        inventory.add({chord.symbolic_y(), chord.left_edge, chord.left_side, chord.right_edge,
                       chord.right_side, chord.is_null_length},
                      chord.left_adj.arcs[0], chord.right_adj.arcs[0], true);
    }
    Submap result = inventory.build();
    result.build_tree_decomposition();
    assert(result.is_conformal() && result.is_semigranular(1) &&
           "[C91 §2.1 tex 70–72]: V(C) contains exactly the original vertex chords");
    return result;
}

Submap normalize_refinement(const Polygon& curve, const Submap& coarse, const Submap& refined) {
    std::vector<std::size_t> owners(refined.num_nodes(), NONE);
    std::vector<std::size_t> pending;
    auto seed = [&](std::size_t region, std::size_t owner) {
        if (owners[region] != NONE) {
            assert(owners[region] == owner &&
                   "[C91 §4.2 tex 377]: refining preserves the old region boundaries");
            return;
        }
        owners[region] = owner;
        pending.push_back(region);
    };
    if (coarse.num_chords() == 0)
        seed(0, 0);
    for (std::size_t i = 0; i < coarse.num_chords(); ++i) {
        if (coarse.chord(i).is_null_length) {
            for (std::size_t side = 0; side < 2; ++side)
                seed(refined.chord(i).region[side], coarse.chord(i).region[side]);
        } else {
            std::size_t old_below, old_above, new_below, new_above;
            coarse.chord_regions_below_above(i, curve, &old_below, &old_above);
            refined.chord_regions_below_above(i, curve, &new_below, &new_above);
            seed(new_below, old_below);
            seed(new_above, old_above);
        }
    }
    for (std::size_t i = 0; i < pending.size(); ++i) {
        const std::size_t region = pending[i];
        for (std::size_t chord_index : refined.node(region).incident_chords) {
            if (chord_index < coarse.num_chords())
                continue;
            const Chord& chord = refined.chord(chord_index);
            const std::size_t neighbor =
                chord.region[0] == region ? chord.region[1] : chord.region[0];
            seed(neighbor, owners[region]);
        }
    }
    assert(pending.size() == refined.num_nodes());
    RefinementInventory inventory(curve, coarse);
    inherit(inventory, coarse);
    for (std::size_t i = coarse.num_chords(); i < refined.num_chords(); ++i) {
        const Chord& chord = refined.chord(i);
        const RegionArcs arcs = collect_region_arcs(coarse, owners[chord.region[0]]);
        PendingChord value{chord.symbolic_y(), chord.left_edge,  chord.left_side,
                           chord.right_edge,   chord.right_side, chord.is_null_length};
        canonicalize_chord(value, curve);
        const ChordEndpoint first{value.left_edge_c, value.left_side, value.y, 0, true};
        const ChordEndpoint second{value.right_edge_c, value.right_side, value.y, 0, false};
        const std::size_t first_arc = endpoint_owner(curve, coarse, arcs.arcs[0], first);
        const std::size_t second_arc = endpoint_owner(curve, coarse, arcs.arcs[0], second);
        inventory.add(std::move(value), first_arc, second_arc, false);
    }
    return inventory.build();
}

void complete_boundary_levels(RefinementInventory& inventory, const Polygon& curve,
                              const Submap& coarse, const RegionArcs& arcs,
                              const UpPhaseRayShooter& rays,
                              const std::vector<std::size_t>& vertices) {
    for (std::size_t vertex : vertices) {
        const SymbolicY level = symbolic_y_of(curve.vertex(vertex));
        for (std::size_t source_arc : arcs) {
            for (std::size_t edge :
                 {vertex > 0 ? vertex - 1 : NONE, vertex < curve.num_edges() ? vertex : NONE}) {
                if (edge == NONE)
                    continue;
                for (Side side : {LEFT, RIGHT}) {
                    const ChordEndpoint source{edge, side, level, 0, true};
                    if (!contains(curve, coarse, source_arc, source))
                        continue;
                    if (is_inside_companion(curve, edge, side, vertex)) {
                        const Side inside =
                            is_inside_companion(curve, vertex, LEFT, vertex) ? LEFT : RIGHT;
                        inventory.add({level, vertex, inside, vertex, inside, true}, source_arc,
                                      source_arc, false);
                        continue;
                    }
                    const Point& origin = curve.vertex(vertex);
                    const Side direction = shooting_direction(edge, side, curve);
                    const SourceOffset offset = perturbed_x_offset(curve, level, edge);
                    RayHit best;
                    Exact best_distance;
                    for (std::size_t arc : arcs) {
                        RayHit hit = rays.shoot(origin, direction, arc,
                                                whole_arc(curve, coarse, arc), offset);
                        if (!hit.hit)
                            continue;
                        const Exact distance =
                            direction == RIGHT ? hit.x - origin.x : origin.x - hit.x;
                        if (!best.hit ||
                            (hit.wrapped != best.wrapped ? !hit.wrapped
                             : distance != best_distance
                                 ? distance < best_distance
                                 : ray_contact_precedes(curve, level, direction, hit.edge, hit.side,
                                                        best.edge, best.side))) {
                            best = hit;
                            best.hit_arc_idx = arc;
                            best_distance = distance;
                        }
                    }
                    assert(best.hit &&
                           "[C91 §4.2 tex 367–372]: boundary vertex chords stay in the region");
                    inventory.add({level, edge, side, best.edge, best.side, false}, source_arc,
                                  best.hit_arc_idx, false);
                }
            }
        }
    }
}

Submap refine(const UpPhase& up_phase, const Polygon& curve, const Submap& submap,
              std::size_t grade, std::unique_ptr<RegionBoundaryGeometry>& geometry,
              std::size_t& boundary_count) {
    const std::size_t next_grade = ceil_beta(grade);
    const std::size_t granularity = UpPhase::grade_granularity(next_grade);
    RefinementInventory inventory(curve, submap);
    inherit(inventory, submap);
    UpPhaseRayShooter boundary_rays(up_phase, submap, curve, grade);
    for (std::size_t region = 0; region < submap.num_nodes(); ++region) {
        const RegionArcs arcs = collect_region_arcs(submap, region);
        bool smaller = true;
        for (std::size_t arc : arcs)
            smaller = smaller && submap.arc(arc).edge_count <= granularity;
        if (smaller)
            continue;
        if (!geometry)
            geometry = std::make_unique<RegionBoundaryGeometry>(up_phase, curve);
        const RegionBoundary boundary = geometry->boundary(curve, submap, region);
        const auto canonical = canonical_region_boundary(*geometry, boundary, next_grade);
        ++boundary_count;
        std::vector<std::size_t> boundary_vertices{0, curve.num_vertices() - 1};
        for (std::size_t chord : submap.node(region).incident_chords) {
            const std::size_t vertex = curve.local_index_of_tag(submap.chord(chord).y_tag);
            if (vertex != NONE &&
                symbolic_y_equal(submap.chord(chord).symbolic_y(),
                                 symbolic_y_of(curve.vertex(vertex))) &&
                std::find(boundary_vertices.begin(), boundary_vertices.end(), vertex) ==
                    boundary_vertices.end())
                boundary_vertices.push_back(vertex);
        }
        complete_boundary_levels(inventory, curve, submap, arcs, boundary_rays, boundary_vertices);
        for (std::size_t i = 0; i < canonical.submap.num_chords(); ++i) {
            const Chord& chord = canonical.submap.chord(i);
            auto first = boundary.original_location(chord.left_edge, chord.left_side);
            auto second = boundary.original_location(chord.right_edge, chord.right_side);
            if (first.arc == NONE || second.arc == NONE)
                continue;
            const SymbolicY level = boundary.original_level(chord.y_tag, up_phase.graded().curve());
            const std::size_t vertex = curve.local_index_of_tag(level.tag);
            if (vertex != NONE && symbolic_y_equal(level, symbolic_y_of(curve.vertex(vertex))) &&
                std::find(boundary_vertices.begin(), boundary_vertices.end(), vertex) !=
                    boundary_vertices.end()) {
                if (first.edge == vertex || first.edge + 1 == vertex || second.edge == vertex ||
                    second.edge + 1 == vertex)
                    continue;
                const Exact first_x = edge_x_at_y(curve, first.edge, level);
                const Exact second_x = edge_x_at_y(curve, second.edge, level);
                const Exact vertex_x = curve.vertex(vertex).x;
                const Side direction = shooting_direction(first.edge, first.side, curve);
                const Exact first_offset = perturbed_x_offset(curve, level, first.edge);
                const Exact second_offset = perturbed_x_offset(curve, level, second.edge);
                const bool forward = direction == RIGHT
                                         ? (second_x > first_x ||
                                            (second_x == first_x && second_offset > first_offset))
                                         : (second_x < first_x ||
                                            (second_x == first_x && second_offset < first_offset));
                const bool blocks = direction == RIGHT
                                        ? (forward ? first_x < vertex_x && vertex_x < second_x
                                                   : first_x < vertex_x || vertex_x < second_x)
                                        : (forward ? second_x < vertex_x && vertex_x < first_x
                                                   : vertex_x < first_x || second_x < vertex_x);
                if (blocks)
                    continue;
            }
            const bool null = vertex != NONE &&
                              symbolic_y_equal(level, symbolic_y_of(curve.vertex(vertex))) &&
                              (first.edge == vertex || first.edge + 1 == vertex) &&
                              (second.edge == vertex || second.edge + 1 == vertex) &&
                              is_inside_companion(curve, first.edge, first.side, vertex) &&
                              is_inside_companion(curve, second.edge, second.side, vertex);
            if (first.edge == second.edge && first.side == second.side && !null)
                continue;
            PendingChord extracted{level, first.edge, first.side, second.edge, second.side, null};
            if (null) {
                extracted.left_edge_c = extracted.right_edge_c = vertex;
                extracted.left_side = extracted.right_side =
                    is_inside_companion(curve, vertex, LEFT, vertex) ? LEFT : RIGHT;
            }
            inventory.add(std::move(extracted), first.arc, second.arc, false);
        }
    }
    Submap result = inventory.build();
    UpPhaseRayShooter rays(up_phase, result, curve, grade);
    UpPhaseArcCutter cutter(up_phase, result, curve, grade);
    restore_conformality(result, curve, rays, cutter, UpPhaseArcCutter::piece_count_bound(grade),
                         UpPhaseArcCutter::piece_granularity_bound(grade));
    assert(result.is_semigranular(granularity) &&
           "[C91 §4.2 tex 377]: extracted arcs retain the smaller granularity");
    result = normalize_refinement(curve, submap, result);
    return restore_granularity(std::move(result), curve, granularity);
}

}

Submap refine_visibility_submap(const UpPhase& up_phase, const Polygon& curve, const Submap& submap,
                                std::size_t grade) {
    assert(std::has_single_bit(curve.num_edges()) &&
           curve.table_offset() % curve.num_edges() == 0 &&
           grade < static_cast<std::size_t>(std::bit_width(curve.num_edges())) &&
           "[C91 Lemma 4.2 tex 364]: C is a chain in grade l >= lambda");
    assert(grade >= 2 && submap.is_conformal() &&
           submap.is_granular(UpPhase::grade_granularity(grade), curve) &&
           !submap.tree_decomposition().empty() &&
           "[C91 Lemma 4.2 tex 364]: refinement starts with a normal-form granular submap");
    std::unique_ptr<RegionBoundaryGeometry> geometry;
    std::size_t boundary_count = 0;
    return refine(up_phase, curve, submap, grade, geometry, boundary_count);
}

DownPhaseResult compute_visibility_map(const UpPhase& up_phase, std::size_t grade,
                                       std::size_t chain_index) {
    animation_checkpoint("down_phase", grade);
    const Polygon& curve = up_phase.graded().chain(grade, chain_index);
    Submap submap = up_phase.chain_submap(grade, chain_index);
    if (auto* trace = AnimationTrace::current())
        trace->copy_submap(up_phase.chain_submap(grade, chain_index), submap, curve);
    std::size_t parameter = grade;
    std::unique_ptr<RegionBoundaryGeometry> geometry;
    DownPhaseResult result;
    while (UpPhase::grade_granularity(parameter) > REGION_COMPLETION_GRANULARITY) {
        animation_checkpoint("refine", parameter);
        submap = refine(up_phase, curve, submap, parameter, geometry, result.region_boundaries);
        const std::size_t next = ceil_beta(parameter);
        assert(next < parameter &&
               "[C91 §4.2 tex 379]: the inductive parameter strictly decreases");
        if (auto* trace = AnimationTrace::current())
            trace->submap("refined", curve, submap, UpPhase::grade_granularity(next));
        parameter = next;
        ++result.refinement_rounds;
    }
    animation_checkpoint("bounded_regions", parameter);
    RegionCompletion completion =
        complete_bounded_regions(curve, submap, UpPhase::grade_granularity(parameter));
    result.completion_work = completion.work;
    result.visibility_map = extract_visibility_map(curve, completion.visibility_map);
    if (auto* trace = AnimationTrace::current())
        trace->submap("visibility", curve, result.visibility_map, 1);
    return result;
}

DownPhaseResult compute_visibility_map(const UpPhase& up_phase) {
    return compute_visibility_map(up_phase, up_phase.graded().maximum_grade(), 0);
}

}
