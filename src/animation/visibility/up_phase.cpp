#include "up_phase.h"
#include "../merge/granularity.h"
#include "../trace.h"

#include <algorithm>
#include <utility>
#include <vector>

namespace chazelle::animation {

namespace {

struct ChainReference {
    std::size_t grade;
    std::size_t index;
};

std::vector<ChainReference> partition_into_aligned_chains(std::size_t first_vertex,
                                                          std::size_t last_vertex,
                                                          std::size_t grade_limit_exclusive) {
    assert(first_vertex < last_vertex && "a chain partition requires a nonempty interval");
    assert(grade_limit_exclusive >= 1 &&
           last_vertex - first_vertex <= (std::size_t{1} << grade_limit_exclusive) &&
           "[C91 §4.1 tex 339]: the interval must fit within one "
           "chain of the limiting grade");
    std::vector<ChainReference> prefix_chains, suffix_chains;
    for (std::size_t piece_grade = 0; first_vertex < last_vertex; ++piece_grade) {
        const std::size_t chain_length = std::size_t{1} << piece_grade;
        if (piece_grade + 1 == grade_limit_exclusive) {
            assert(
                first_vertex % chain_length == 0 && last_vertex % chain_length == 0 &&
                last_vertex - first_vertex <= 2 * chain_length &&
                "the remaining interval contains one or two chains in the highest permitted grade");
            for (; first_vertex < last_vertex; first_vertex += chain_length)
                prefix_chains.push_back({piece_grade, first_vertex / chain_length});
            break;
        }
        if ((first_vertex / chain_length) % 2 == 1) {
            prefix_chains.push_back({piece_grade, first_vertex / chain_length});
            first_vertex += chain_length;
        }
        if ((last_vertex / chain_length) % 2 == 1) {
            suffix_chains.push_back({piece_grade, last_vertex / chain_length - 1});
            last_vertex -= chain_length;
        }
    }
    for (std::size_t index = suffix_chains.size(); index-- > 0;)
        prefix_chains.push_back(suffix_chains[index]);
#ifndef NDEBUG
    for (std::size_t index = 0; index + 1 < prefix_chains.size(); ++index)
        assert(((prefix_chains[index].index + 1) << prefix_chains[index].grade) ==
                   (prefix_chains[index + 1].index << prefix_chains[index + 1].grade) &&
               "partition chains are contiguous in ascending order");
#endif
    return prefix_chains;
}

std::vector<ChainReference>
partition_piece_into_chains(std::size_t first_vertex, std::size_t last_vertex,
                            std::size_t merge_grade,
                            [[maybe_unused]] std::size_t maximum_piece_grade) {
    assert(maximum_piece_grade < merge_grade && last_vertex > first_vertex &&
           last_vertex - first_vertex <= (std::size_t{1} << maximum_piece_grade) &&
           "[C91 §4.1 proof tex 341]: piece length ≤ 2^{⌈βλ⌉}");
    std::vector<ChainReference> chains =
        partition_into_aligned_chains(first_vertex, last_vertex, merge_grade);
    for ([[maybe_unused]] const ChainReference& chain : chains)
        assert(chain.grade <= maximum_piece_grade &&
               "[C91 §4.1 proof tex 341]: every chain grade ≤ ⌈βλ⌉");
    return chains;
}

struct SubarcPiece {
    bool is_endpoint_piece = false;
    Subarc subarc{};
    ChainReference chain{NONE, NONE};
};

std::size_t boundary_edge_start_vertex(std::size_t edge, Side side) {
    return (side == LEFT) ? edge : edge + 1;
}
std::size_t boundary_edge_end_vertex(std::size_t edge, Side side) {
    return (side == LEFT) ? edge + 1 : edge;
}

std::vector<SubarcPiece> decompose_subarc(const Polygon& input_curve, const Subarc& target,
                                          std::size_t merge_grade) {
    assert(target.first_y.tag != SOS_NONE && target.last_y.tag != SOS_NONE &&
           "[C91 §3.0(i) tex 169]: α' is specified by its two exact "
           "endpoints");
    const std::size_t input_table_offset = input_curve.table_offset();
    const std::size_t maximum_piece_grade = ceil_beta(merge_grade);

    ArcSideRange ranges[3];
    const std::size_t range_count =
        subarc_side_ranges(target, 0, input_curve.num_vertices() - 1, ranges);
    assert(range_count >= 1 && range_count <= 3 &&
           "[C91 §4.1 tex 341]: a subarc has a constant number of "
           "pieces lying on one side (1–3 pieces, [C91 §2.4 tex 142])");

    std::vector<SubarcPiece> result;
    for (std::size_t range_index = 0; range_index < range_count; ++range_index) {
        const Side side = ranges[range_index].side;
        const std::size_t first_edge = ranges[range_index].first_edge;
        const std::size_t last_edge = ranges[range_index].last_edge;

        const std::size_t first_traversed_edge = (side == LEFT) ? first_edge : last_edge;
        const std::size_t last_traversed_edge = (side == LEFT) ? last_edge : first_edge;

        const bool first_partial =
            (range_index == 0) &&
            !symbolic_y_equal(target.first_y,
                              symbolic_y_of(input_curve.vertex(
                                  boundary_edge_start_vertex(first_traversed_edge, side))));
        const bool last_partial =
            (range_index + 1 == range_count) &&
            !symbolic_y_equal(target.last_y,
                              symbolic_y_of(input_curve.vertex(
                                  boundary_edge_end_vertex(last_traversed_edge, side))));

        if (first_edge == last_edge && first_partial && last_partial) {
            SubarcPiece boundary_piece;
            boundary_piece.is_endpoint_piece = true;
            boundary_piece.subarc =
                Subarc{first_edge, side, first_edge, side, target.first_y, target.last_y};
            result.push_back(boundary_piece);
            continue;
        }

        if (first_partial) {
            SubarcPiece boundary_piece;
            boundary_piece.is_endpoint_piece = true;
            boundary_piece.subarc =
                Subarc{first_traversed_edge,
                       side,
                       first_traversed_edge,
                       side,
                       target.first_y,
                       symbolic_y_of(input_curve.vertex(
                           boundary_edge_end_vertex(first_traversed_edge, side)))};
            result.push_back(boundary_piece);
        }

        std::size_t first_full_vertex = first_edge + (first_partial && side == LEFT ? 1 : 0) +
                                        (last_partial && side == RIGHT ? 1 : 0);
        std::size_t last_full_vertex = last_edge + 1 - (last_partial && side == LEFT ? 1 : 0) -
                                       (first_partial && side == RIGHT ? 1 : 0);
        if (first_full_vertex < last_full_vertex) {
            std::vector<ChainReference> chains = partition_piece_into_chains(
                input_table_offset + first_full_vertex, input_table_offset + last_full_vertex,
                merge_grade, maximum_piece_grade);

            auto append_chain = [&](const ChainReference& chain) {
                const std::size_t first_global_edge = chain.index << chain.grade;
                const std::size_t last_global_edge =
                    first_global_edge + (std::size_t{1} << chain.grade) - 1;
                SubarcPiece chain_piece;
                chain_piece.chain = chain;
                const std::size_t first_local_edge = first_global_edge - input_table_offset;
                const std::size_t last_local_edge = last_global_edge - input_table_offset;
                if (side == LEFT) {
                    chain_piece.subarc =
                        Subarc{first_local_edge,
                               side,
                               last_local_edge,
                               side,
                               symbolic_y_of(input_curve.vertex(first_local_edge)),
                               symbolic_y_of(input_curve.vertex(last_local_edge + 1))};
                } else {
                    chain_piece.subarc =
                        Subarc{last_local_edge,
                               side,
                               first_local_edge,
                               side,
                               symbolic_y_of(input_curve.vertex(last_local_edge + 1)),
                               symbolic_y_of(input_curve.vertex(first_local_edge))};
                }
                result.push_back(chain_piece);
            };
            if (side == LEFT)
                for (const ChainReference& chain : chains)
                    append_chain(chain);
            else
                for (std::size_t index = chains.size(); index-- > 0;)
                    append_chain(chains[index]);
        }

        if (last_partial) {
            SubarcPiece boundary_piece;
            boundary_piece.is_endpoint_piece = true;
            boundary_piece.subarc =
                Subarc{last_traversed_edge,
                       side,
                       last_traversed_edge,
                       side,
                       symbolic_y_of(input_curve.vertex(
                           boundary_edge_start_vertex(last_traversed_edge, side))),
                       target.last_y};
            result.push_back(boundary_piece);
        }
    }

    assert(!result.empty());
    assert(result.front().subarc.first_edge == target.first_edge &&
           result.front().subarc.first_side == target.first_side &&
           symbolic_y_equal(result.front().subarc.first_y, target.first_y) &&
           "decomposition starts at α'.first");
    assert(result.back().subarc.last_edge == target.last_edge &&
           result.back().subarc.last_side == target.last_side &&
           symbolic_y_equal(result.back().subarc.last_y, target.last_y) &&
           "decomposition ends at α'.last");
    return result;
}

RayHit shoot_toward_single_edge_subarc(const Polygon& curve, const Subarc& subarc,
                                       const Point& origin, Side direction,
                                       const SourceOffset& source_x_offset) {
    const SymbolicY ray_y{origin.y, origin.index};
    const std::size_t edge_index = subarc.first_edge;
    Exact x;
    if (!edge_crossing_x(curve, edge_index, ray_y, &x))
        return RayHit{};
    const SymbolicY first_y = subarc.first_y;
    const SymbolicY last_y = subarc.last_y;
    const bool within = (symbolic_y_leq(first_y, ray_y) && symbolic_y_leq(ray_y, last_y)) ||
                        (symbolic_y_leq(last_y, ray_y) && symbolic_y_leq(ray_y, first_y));
    if (!within)
        return RayHit{};
    const auto& edge = curve.edge(edge_index);
    const bool edge_ascending = symbolic_y_less(symbolic_y_of(curve.vertex(edge.start_idx)),
                                                symbolic_y_of(curve.vertex(edge.end_idx)));
    const Side side_facing_left = edge_ascending ? LEFT : RIGHT;
    RayHit hit;
    hit.hit = true;
    hit.x = x;
    hit.y = origin.y;
    hit.edge = edge_index;
    hit.side = (direction == RIGHT) ? side_facing_left : (side_facing_left == LEFT ? RIGHT : LEFT);
    const Exact distance = (direction == RIGHT) ? (x - origin.x) : (origin.x - x);
    hit.wrapped = (distance < 0.0) ||
                  (distance == 0.0 &&
                   !perturbed_hit_forward(curve, ray_y, direction, source_x_offset, edge_index));
    return hit;
}

void align_hit_edge_with_subarc(const Subarc& target, const Polygon& input_curve, RayHit& hit,
                                const SymbolicY& ray_y) {
    const std::size_t first_curve_vertex = 0;
    const std::size_t last_curve_vertex = input_curve.num_vertices() - 1;
    if (subarc_contains_point(target, input_curve, hit.edge, hit.side, ray_y, first_curve_vertex,
                              last_curve_vertex))
        return;
    const auto& edge = input_curve.edge(hit.edge);
    std::size_t vertex_index = NONE;
    if (symbolic_y_equal(ray_y, symbolic_y_of(input_curve.vertex(edge.start_idx))))
        vertex_index = edge.start_idx;
    else if (symbolic_y_equal(ray_y, symbolic_y_of(input_curve.vertex(edge.end_idx))))
        vertex_index = edge.end_idx;
    if (vertex_index == NONE || input_curve.is_endpoint(vertex_index) ||
        input_curve.is_y_extremum(vertex_index))
        return;
    const std::size_t adjacent_edge = (hit.edge == vertex_index) ? vertex_index - 1 : vertex_index;
    if (adjacent_edge == hit.edge)
        return;
    if (subarc_contains_point(target, input_curve, adjacent_edge, hit.side, ray_y,
                              first_curve_vertex, last_curve_vertex))
        hit.edge = adjacent_edge;
}

}

UpPhaseRayShooter::UpPhaseRayShooter(const UpPhase& up, const Submap& input_submap,
                                     const Polygon& input_curve, std::size_t merge_grade)
    : up_phase_(&up), input_submap_(&input_submap), input_curve_(&input_curve),
      merge_grade_(merge_grade) {
    assert(ceil_beta(merge_grade) < merge_grade &&
           "[C91 §4.1 tex 341]: ⌈βλ⌉ < λ (early grades are naive)");
}

RayHit UpPhaseRayShooter::shoot(Point origin, Side direction, [[maybe_unused]] std::size_t arc_idx,
                                const Subarc& target, SourceOffset source_x_offset) const {
    assert(arc_idx < input_submap_->num_arcs() && !input_submap_->arc(arc_idx).dead &&
           "[C91 §3.0 tex 166]: α is specified by its arc-structure");
    assert_subarc_clockwise(target);
    const SymbolicY ray_y{origin.y, origin.index};
    const std::size_t input_table_offset = input_curve_->table_offset();

    std::vector<SubarcPiece> pieces = decompose_subarc(*input_curve_, target, merge_grade_);

    RayHit nearest_hit;
    nearest_hit.hit = false;
    Exact nearest_distance = 0.0;
    auto offer = [&](const RayHit& hit) {
        assert(hit.hit);
        const Exact distance = (direction == RIGHT) ? (hit.x - origin.x) : (origin.x - hit.x);

        assert((distance < 0.0 ? hit.wrapped : (distance > 0.0 ? !hit.wrapped : true)) &&
               "[C91 §2.1 tex 70]: wrap flag consistent with the signed "
               "travel distance");
        bool better;
        if (!nearest_hit.hit)
            better = true;
        else if (hit.wrapped != nearest_hit.wrapped)
            better = !hit.wrapped;
        else if (distance != nearest_distance)
            better = distance < nearest_distance;
        else
            better = ray_contact_precedes(*input_curve_, ray_y, direction, hit.edge, hit.side,
                                          nearest_hit.edge, nearest_hit.side);
        if (better) {
            nearest_hit = hit;
            nearest_distance = distance;
        }
    };

    for (const SubarcPiece& piece : pieces) {
        if (piece.is_endpoint_piece) {
            RayHit hit = shoot_toward_single_edge_subarc(*input_curve_, piece.subarc, origin,
                                                         direction, source_x_offset);
            if (hit.hit)
                offer(hit);
        } else {
            const RayShootingStructure& structure =
                up_phase_->chain_structure(piece.chain.grade, piece.chain.index);

            RayHit hit = structure.shoot_toward_boundary(origin, direction, source_x_offset);
            if (!hit.hit)
                continue;

            hit.edge = (piece.chain.index << piece.chain.grade) + hit.edge - input_table_offset;
            offer(hit);
        }
    }

    if (!nearest_hit.hit)
        return RayHit{};

    align_hit_edge_with_subarc(target, *input_curve_, nearest_hit, ray_y);
    return nearest_hit;
}

UpPhaseArcCutter::UpPhaseArcCutter(const UpPhase& up, const Submap& input_submap,
                                   const Polygon& input_curve, std::size_t merge_grade)
    : up_phase_(&up), input_submap_(&input_submap), input_curve_(&input_curve),
      merge_grade_(merge_grade) {
    assert(ceil_beta(merge_grade) < merge_grade &&
           "[C91 §4.1 tex 341]: ⌈βλ⌉ < λ (early grades are naive)");
}

std::vector<ArcPiece> UpPhaseArcCutter::cut([[maybe_unused]] std::size_t arc_idx,
                                            const Subarc& target) const {
    assert(arc_idx < input_submap_->num_arcs() && !input_submap_->arc(arc_idx).dead &&
           "[C91 §3.0 tex 166]: α is specified by its arc-structure");
    assert_subarc_clockwise(target);
    [[maybe_unused]] const std::size_t input_table_offset = input_curve_->table_offset();
    [[maybe_unused]] const std::size_t maximum_piece_granularity =
        piece_granularity_bound(merge_grade_);

    std::vector<SubarcPiece> pieces = decompose_subarc(*input_curve_, target, merge_grade_);

    std::vector<ArcPiece> result;
    result.reserve(pieces.size());
    for (const SubarcPiece& piece : pieces) {
        ArcPiece arc_piece;
        arc_piece.subarc = piece.subarc;
        if (piece.is_endpoint_piece) {
            arc_piece.is_boundary_piece = true;
        } else {
            const std::size_t piece_grade = piece.chain.grade;
            const std::size_t chain_index = piece.chain.index;
            arc_piece.curve = &up_phase_->graded().chain(piece_grade, chain_index);
            arc_piece.submap = &up_phase_->chain_submap(piece_grade, chain_index);

            assert(piece_grade <= ceil_beta(merge_grade_) &&
                   "[C91 §4.1 proof tex 341]: chain grade ≤ ⌈βλ⌉");
            arc_piece.granularity = UpPhase::grade_granularity(piece_grade);
            assert(arc_piece.granularity <= maximum_piece_granularity &&
                   "[C91 §4.1 tex 346]: h(γ) ≤ 2^{⌈β⌈βλ⌉⌉}");

            assert(arc_piece.curve->table_offset() == (chain_index << piece_grade) &&
                   arc_piece.curve->num_edges() == (std::size_t{1} << piece_grade) &&
                   "[C91 §4 tex 316]: chain spans its aligned "
                   "vertex range");
            assert((chain_index << piece_grade) >= input_table_offset &&
                   "chain lies within Cᵢ's table range");
        }
        result.push_back(arc_piece);
    }

    assert(result.size() <= piece_count_bound(merge_grade_) &&
           "[C91 §4.1 tex 344]: g(γ) = O(λ) pieces");
    return result;
}

UpPhase::PortionResult UpPhase::compute_canonical_portion(std::size_t first_vertex,
                                                          std::size_t last_vertex) const {
    [[maybe_unused]] const Polygon& input_curve = graded_.curve();
    assert(first_vertex < last_vertex && last_vertex < input_curve.num_vertices() &&
           "[C91 Lemma 4.1 tex 336]: D = v_a, ..., v_b within P");

    std::size_t merge_grade = 0;
    while ((std::size_t{1} << merge_grade) < last_vertex - first_vertex)
        ++merge_grade;
    assert(merge_grade >= 1 &&
           (last_vertex - first_vertex) > (std::size_t{1} << (merge_grade - 1)) &&
           (last_vertex - first_vertex) <= (std::size_t{1} << merge_grade) &&
           "[C91 Lemma 4.1 tex 336]: 2^{λ−1} < b − a ≤ 2^λ");
    assert(ceil_beta(merge_grade) < merge_grade &&
           "[C91 §4.1 tex 341]: Lemma 4.1 requires ⌈βλ⌉ < λ — smaller "
           "portions are the naive early grades (tex 333)");

    const std::size_t granularity = grade_granularity(merge_grade);

    std::vector<ChainReference> parts =
        partition_into_aligned_chains(first_vertex, last_vertex, merge_grade);
    assert(parts.size() >= 1 && parts.size() <= 2 * merge_grade &&
           "[C91 §4.1 tex 339]: j ≤ 2λ partition chains");

    std::vector<PortionResult> merge_level;
    merge_level.reserve(parts.size());
    for (const ChainReference& part : parts)
        merge_level.push_back(reset_chain_granularity(part.grade, part.index, granularity));

    while (merge_level.size() > 1) {
        std::vector<PortionResult> next_level;
        next_level.reserve((merge_level.size() + 1) / 2);
        for (std::size_t i = 0; i + 1 < merge_level.size(); i += 2)
            next_level.push_back(
                merge_portions(merge_level[i], merge_level[i + 1], merge_grade, granularity));
        if (merge_level.size() % 2 == 1)
            next_level.push_back(std::move(merge_level.back()));
        merge_level = std::move(next_level);
    }

    PortionResult result{merge_level[0].curve, std::move(merge_level[0].submap)};

    assert(result.curve.num_vertices() == last_vertex - first_vertex + 1 &&
           result.curve.table_offset() == first_vertex &&
           "[C91 Lemma 4.1 tex 336]: the result covers exactly D");
    assert(result.submap.is_conformal() &&
           "[C91 Lemma 4.1 tex 341]: the portion's submap is conformal");
    assert(result.submap.is_granular(granularity, result.curve) &&
           "[C91 Lemma 4.1 tex 341]: the portion's submap is "
           "2^{⌈βλ⌉}-granular");
    assert(!result.submap.tree_decomposition().empty() &&
           "[C91 Lemma 4.1 tex 341]: normal form carries the tree "
           "decomposition");
    return result;
}

UpPhase::PortionResult UpPhase::reset_chain_granularity(std::size_t grade, std::size_t index,
                                                        std::size_t granularity) const {
    const Polygon& curve = graded_.chain(grade, index);
    Submap submap = chain_submap(grade, index);
    if (auto* trace = AnimationTrace::current())
        trace->copy_submap(chain_submap(grade, index), submap, curve);
    assert(grade_granularity(grade) <= granularity &&
           "[C91 §4.1 tex 339]: canonical granularity grows "
           "monotonically with curve size");
    enforce_granularity(submap, curve, granularity);
    submap.normalize(curve);
    assert(submap.is_granular(granularity, curve) &&
           "[C91 §3.3 tex 276]: reset yields a γ-granular submap");
    return PortionResult{curve, std::move(submap)};
}

UpPhase::PortionResult UpPhase::merge_portions(const PortionResult& left,
                                               const PortionResult& right, std::size_t merge_grade,
                                               std::size_t granularity) const {
    UpPhaseRayShooter left_ray_shooter(*this, left.submap, left.curve, merge_grade);
    UpPhaseRayShooter right_ray_shooter(*this, right.submap, right.curve, merge_grade);
    UpPhaseArcCutter left_arc_cutter(*this, left.submap, left.curve, merge_grade);
    UpPhaseArcCutter right_arc_cutter(*this, right.submap, right.curve, merge_grade);

    const std::size_t maximum_piece_count = UpPhaseArcCutter::piece_count_bound(merge_grade);
    const std::size_t maximum_piece_granularity =
        UpPhaseArcCutter::piece_granularity_bound(merge_grade);
    const MergeInput merge_input{
        .first_curve = &left.curve,
        .second_curve = &right.curve,
        .first_submap = &left.submap,
        .second_submap = &right.submap,
        .first_granularity = granularity,
        .second_granularity = granularity,
        .granularity = granularity,
        .first_ray_shooter = &left_ray_shooter,
        .second_ray_shooter = &right_ray_shooter,
        .first_arc_cutter = &left_arc_cutter,
        .second_arc_cutter = &right_arc_cutter,
        .first_piece_count_bound = maximum_piece_count,
        .second_piece_count_bound = maximum_piece_count,
        .first_piece_granularity_bound = maximum_piece_granularity,
        .second_piece_granularity_bound = maximum_piece_granularity,
    };
    MergeResult result = merge(merge_input);
    return PortionResult{result.curve, std::move(result.submap)};
}

UpPhase::UpPhase(std::vector<Point> vertices) : graded_(std::move(vertices)) {
    animation_checkpoint("up_phase");
    if (auto* trace = AnimationTrace::current())
        trace->boundary(graded_.curve());
    const std::size_t maximum_grade = graded_.maximum_grade();
    submaps_.resize(maximum_grade + 1);
    structures_.resize(maximum_grade + 1);

    for (std::size_t grade = 0; grade <= maximum_grade; ++grade) {
        animation_checkpoint("grade", grade);
        const std::size_t chain_count = graded_.num_chains(grade);
        submaps_[grade].reserve(chain_count);
        structures_[grade].reserve(chain_count);
        for (std::size_t chain_index = 0; chain_index < chain_count; ++chain_index) {
            const Polygon& curve = graded_.chain(grade, chain_index);
            if (auto* trace = AnimationTrace::current())
                trace->chain(grade, chain_index, curve);

            Submap submap;
            if (grade <= NAIVE_MAX_GRADE) {
                submap = build_canonical_submap_naive(curve);
            } else {
                PortionResult result =
                    compute_canonical_portion(chain_index << grade, (chain_index + 1) << grade);
                submap = std::move(result.submap);
            }
            assert(submap.is_conformal() && submap.is_granular(grade_granularity(grade), curve) &&
                   !submap.tree_decomposition().empty() &&
                   "[C91 §4.1 tex 327]: canonical = 2^{⌈βλ⌉}-granular, "
                   "conformal, normal form");
            if (auto* trace = AnimationTrace::current())
                trace->submap("canonical", curve, submap, grade_granularity(grade));
            submaps_[grade].push_back(std::move(submap));

            structures_[grade].push_back(std::make_unique<RayShootingStructure>(
                submaps_[grade][chain_index], curve, grade_granularity(grade)));
        }
    }

    assert(submaps_[maximum_grade].size() == 1 &&
           "[C91 §4 tex 319]: grade p has the single chain P");
}

}
