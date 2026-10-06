#include "region_oracles.h"
#include "../merge/granularity.h"

#include <algorithm>
#include <cassert>
#include <utility>

namespace chazelle::animation {

RegionBoundaryOracles::RegionBoundaryOracles(const RegionBoundaryGeometry& geometry,
                                             const RegionBoundary& boundary,
                                             const Polygon& input_curve, std::size_t input_offset,
                                             std::size_t grade)
    : geometry_(&geometry), boundary_(&boundary), input_curve_(&input_curve),
      input_offset_(input_offset), grade_(grade) {}

CanonicalBoundaryChain& RegionBoundaryOracles::single_edge(std::size_t edge) const {
    for (const auto& existing : edges_)
        if (existing.edge == edge)
            return *existing.canonical;
    Polygon curve = input_curve_->subchain(edge, 2);
    Submap submap = build_canonical_submap_naive(curve);
    edges_.push_back({edge, std::make_unique<CanonicalBoundaryChain>(CanonicalBoundaryChain{
                                std::move(curve), std::move(submap), nullptr})});
    assert(edges_.size() <= 128 && "[C91 §4.2 tex 372]: only a constant number of edges are added");
    return *edges_.back().canonical;
}

std::vector<RegionBoundaryOracles::Piece>
RegionBoundaryOracles::partition(const Subarc& target) const {
    const Polygon& curve = *input_curve_;
    ArcSideRange ranges[3];
    const std::size_t range_count = subarc_side_ranges(target, 0, curve.num_vertices() - 1, ranges);
    std::vector<Piece> result;
    for (std::size_t r = 0; r < range_count; ++r) {
        const Side side = ranges[r].side;
        const std::size_t first_traversed =
            side == LEFT ? ranges[r].first_edge : ranges[r].last_edge;
        const std::size_t last_traversed =
            side == LEFT ? ranges[r].last_edge : ranges[r].first_edge;
        const bool first_partial =
            r == 0 && !symbolic_y_equal(target.first_y,
                                        symbolic_y_of(curve.vertex(
                                            side == LEFT ? first_traversed : first_traversed + 1)));
        const bool last_partial =
            r + 1 == range_count &&
            !symbolic_y_equal(
                target.last_y,
                symbolic_y_of(curve.vertex(side == LEFT ? last_traversed + 1 : last_traversed)));
        if (first_traversed == last_traversed && (first_partial || last_partial)) {
            result.push_back(
                {{first_traversed, side, last_traversed, side,
                  r == 0 ? target.first_y
                         : symbolic_y_of(
                               curve.vertex(side == LEFT ? first_traversed : first_traversed + 1)),
                  r + 1 == range_count ? target.last_y
                                       : symbolic_y_of(curve.vertex(
                                             side == LEFT ? last_traversed + 1 : last_traversed))},
                 nullptr,
                 0});
            continue;
        }
        if (first_partial)
            result.push_back({{first_traversed, side, first_traversed, side, target.first_y,
                               symbolic_y_of(curve.vertex(side == LEFT ? first_traversed + 1
                                                                       : first_traversed))},
                              nullptr,
                              0});
        const std::size_t lo =
            ranges[r].first_edge +
            ((first_partial && side == LEFT) || (last_partial && side == RIGHT) ? 1 : 0);
        const std::size_t end =
            ranges[r].last_edge + 1 -
            ((last_partial && side == LEFT) || (first_partial && side == RIGHT) ? 1 : 0);
        std::vector<Piece> full;
        std::size_t boundary_offset = 0;
        for (const RegionBoundaryPiece& boundary_piece : boundary_->pieces) {
            const std::size_t boundary_end = boundary_offset + boundary_piece.curve.num_edges();
            const std::size_t start_global = std::max(input_offset_ + lo, boundary_offset);
            const std::size_t end_global = std::min(input_offset_ + end, boundary_end);
            if (start_global < end_global) {
                if (boundary_piece.first_original_vertex == NONE) {
                    for (std::size_t global_edge = start_global; global_edge < end_global;
                         ++global_edge) {
                        const std::size_t edge = global_edge - input_offset_;
                        full.push_back(
                            {{edge, side, edge, side,
                              symbolic_y_of(curve.vertex(side == LEFT ? edge : edge + 1)),
                              symbolic_y_of(curve.vertex(side == LEFT ? edge + 1 : edge))},
                             &single_edge(edge),
                             1});
                    }
                } else {
                    std::size_t first =
                        boundary_piece.reversed
                            ? boundary_piece.first_original_vertex - (end_global - boundary_offset)
                            : boundary_piece.first_original_vertex + start_global - boundary_offset;
                    const std::size_t last = first + end_global - start_global;
                    assert(
                        last - first <= UpPhase::grade_granularity(grade_) &&
                        "[C91 §4.1 tex 340, §4.2 tex 372]: input subarcs have at most gamma edges");
                    struct Chain {
                        std::size_t first;
                        std::size_t grade;
                    };
                    std::vector<Chain> chains;
                    while (first < last) {
                        std::size_t grade = 0;
                        while (grade < ceil_beta(grade_) &&
                               first % (std::size_t{1} << (grade + 1)) == 0 &&
                               (std::size_t{1} << (grade + 1)) <= last - first)
                            ++grade;
                        chains.push_back({first, grade});
                        first += std::size_t{1} << grade;
                    }
                    if (boundary_piece.reversed)
                        std::reverse(chains.begin(), chains.end());
                    for (const Chain& reference : chains) {
                        const std::size_t length = std::size_t{1} << reference.grade;
                        const std::size_t global_edge =
                            boundary_piece.reversed
                                ? boundary_offset + boundary_piece.first_original_vertex -
                                      reference.first - length
                                : boundary_offset + reference.first -
                                      boundary_piece.first_original_vertex;
                        const std::size_t edge = global_edge - input_offset_;
                        CanonicalBoundaryChain& canonical = geometry_->canonical_chain(
                            reference.grade, reference.first >> reference.grade,
                            boundary_piece.original_side, boundary_piece.reversed);
                        full.push_back(
                            {{side == LEFT ? edge : edge + length - 1, side,
                              side == LEFT ? edge + length - 1 : edge, side,
                              symbolic_y_of(curve.vertex(side == LEFT ? edge : edge + length)),
                              symbolic_y_of(curve.vertex(side == LEFT ? edge + length : edge))},
                             &canonical,
                             UpPhase::grade_granularity(reference.grade)});
                    }
                }
            }
            boundary_offset = boundary_end;
        }
        if (side == RIGHT)
            std::reverse(full.begin(), full.end());
        result.insert(result.end(), full.begin(), full.end());
        if (last_partial)
            result.push_back(
                {{last_traversed, side, last_traversed, side,
                  symbolic_y_of(curve.vertex(side == LEFT ? last_traversed : last_traversed + 1)),
                  target.last_y},
                 nullptr,
                 0});
    }
    assert(
        !result.empty() && result.size() <= piece_count_bound(grade_) &&
        "[C91 §4.1 tex 344, §4.2 tex 372]: original chains and added edges require O(lambda) pieces");
    return result;
}

std::vector<ArcPiece> RegionBoundaryOracles::cut([[maybe_unused]] std::size_t arc,
                                                 const Subarc& target) const {
    std::vector<ArcPiece> result;
    for (const Piece& piece : partition(target))
        result.push_back({piece.subarc, piece.chain ? &piece.chain->submap : nullptr,
                          piece.chain ? &piece.chain->curve : nullptr, piece.chain == nullptr,
                          piece.granularity});
    return result;
}

RayHit RegionBoundaryOracles::shoot(Point origin, Side direction, [[maybe_unused]] std::size_t arc,
                                    const Subarc& target, SourceOffset offset) const {
    const Polygon& curve = *input_curve_;
    const SymbolicY level = symbolic_y_of(origin);
    RayHit best;
    Exact best_distance;
    for (const Piece& piece : partition(target)) {
        RayHit hit;
        if (piece.chain) {
            if (!piece.chain->rays)
                piece.chain->rays = std::make_unique<RayShootingStructure>(
                    piece.chain->submap, piece.chain->curve, piece.granularity);
            hit = piece.chain->rays->shoot_toward_boundary(origin, direction, offset);
            if (hit.hit)
                hit.edge += std::min(piece.subarc.first_edge, piece.subarc.last_edge);
        } else {
            Exact x;
            const std::size_t edge = piece.subarc.first_edge;
            if (edge_crossing_x(curve, edge, level, &x) &&
                ((symbolic_y_leq(piece.subarc.first_y, level) &&
                  symbolic_y_leq(level, piece.subarc.last_y)) ||
                 (symbolic_y_leq(piece.subarc.last_y, level) &&
                  symbolic_y_leq(level, piece.subarc.first_y)))) {
                const Side west =
                    point_y_below(curve.vertex(edge), curve.vertex(edge + 1)) ? LEFT : RIGHT;
                const Exact distance = direction == RIGHT ? x - origin.x : origin.x - x;
                hit = {true,
                       x,
                       level.y,
                       edge,
                       direction == RIGHT ? west : (west == LEFT ? RIGHT : LEFT),
                       distance < 0 ||
                           (distance == 0 &&
                            !perturbed_hit_forward(curve, level, direction, offset, edge)),
                       NONE};
            }
        }
        if (!hit.hit)
            continue;
        const Exact distance = direction == RIGHT ? hit.x - origin.x : origin.x - hit.x;
        if (!best.hit ||
            (hit.wrapped != best.wrapped ? !hit.wrapped
             : distance != best_distance ? distance < best_distance
                                         : ray_contact_precedes(curve, level, direction, hit.edge,
                                                                hit.side, best.edge, best.side))) {
            best = hit;
            best_distance = distance;
        }
    }
    if (best.hit && !subarc_contains_point(target, curve, best.edge, best.side, level, 0,
                                           curve.num_vertices() - 1)) {
        const std::size_t vertex = curve.local_index_of_tag(level.tag);
        if (vertex != NONE && !curve.is_endpoint(vertex) && !curve.is_y_extremum(vertex)) {
            const std::size_t adjacent = best.edge == vertex ? vertex - 1 : vertex;
            if (subarc_contains_point(target, curve, adjacent, best.side, level, 0,
                                      curve.num_vertices() - 1))
                best.edge = adjacent;
        }
    }
    return best;
}

UpPhase::PortionResult canonical_region_boundary(const RegionBoundaryGeometry& geometry,
                                                 const RegionBoundary& boundary,
                                                 std::size_t grade) {
    const std::size_t granularity = UpPhase::grade_granularity(grade);
    auto result = geometry.canonical_piece(boundary.pieces[0]);
    enforce_granularity(result.submap, result.curve, granularity);
    result.submap.normalize(result.curve);
    std::size_t offset = result.curve.num_edges();
    for (std::size_t i = 1; i < boundary.pieces.size(); ++i) {
        auto next = geometry.canonical_piece(boundary.pieces[i]);
        enforce_granularity(next.submap, next.curve, granularity);
        next.submap.normalize(next.curve);
        RegionBoundaryOracles first_oracles(geometry, boundary, result.curve, 0, grade);
        RegionBoundaryOracles second_oracles(geometry, boundary, next.curve, offset, grade);
        const MergeInput input{&result.curve,
                               &next.curve,
                               &result.submap,
                               &next.submap,
                               granularity,
                               granularity,
                               granularity,
                               &first_oracles,
                               &second_oracles,
                               &first_oracles,
                               &second_oracles,
                               RegionBoundaryOracles::piece_count_bound(grade),
                               RegionBoundaryOracles::piece_count_bound(grade),
                               UpPhaseArcCutter::piece_granularity_bound(grade),
                               UpPhaseArcCutter::piece_granularity_bound(grade)};
        MergeResult merged = merge(input);
        result = {std::move(merged.curve), std::move(merged.submap)};
        offset = result.curve.num_edges();
    }
    assert(
        result.submap.is_conformal() && result.submap.is_granular(granularity, result.curve) &&
        "[C91 §4.2 tex 367–372]: constant-many merges produce the smaller-granularity map of R*");
    return result;
}

}
