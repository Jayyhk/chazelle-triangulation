#pragma once

#include "../merge/merge.h"
#include "../merge/oracle.h"
#include "../merge/ray_shooting.h"
#include "../polygon/polygon.h"
#include "../submap/submap.h"
#include "chain.h"
#include "naive_visibility.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace chazelle::animation {

class UpPhase {
public:
    explicit UpPhase(std::vector<Point> vertices);

    const GradedCurve& graded() const noexcept {
        return graded_;
    }

    static std::size_t grade_granularity(std::size_t merge_grade) noexcept {
        return std::size_t{1} << ceil_beta(merge_grade);
    }

    static constexpr std::size_t NAIVE_MAX_GRADE = 1;

    const Submap& chain_submap(std::size_t grade, std::size_t i) const {
        assert(grade < submaps_.size() && i < submaps_[grade].size() &&
               "[C91 §4.1 tex 333]: grade must be processed already");
        return submaps_[grade][i];
    }

    const RayShootingStructure& chain_structure(std::size_t grade, std::size_t i) const {
        assert(grade < structures_.size() && i < structures_[grade].size() &&
               "[C91 §4.1 tex 333]: grade must be processed already");
        return *structures_[grade][i];
    }

    struct PortionResult {
        Polygon curve;
        Submap submap;
    };
    PortionResult compute_canonical_portion(std::size_t first_vertex,
                                            std::size_t last_vertex) const;

private:
    PortionResult reset_chain_granularity(std::size_t grade, std::size_t index,
                                          std::size_t granularity) const;
    PortionResult merge_portions(const PortionResult& left, const PortionResult& right,
                                 std::size_t merge_grade, std::size_t granularity) const;

    GradedCurve graded_;

    std::vector<std::vector<Submap>> submaps_;
    std::vector<std::vector<std::unique_ptr<RayShootingStructure>>> structures_;
};

class UpPhaseRayShooter final : public RayShootingOracle {
public:
    UpPhaseRayShooter(const UpPhase& up, const Submap& input_submap, const Polygon& input_curve,
                      std::size_t merge_grade);

    RayHit shoot(Point origin, Side direction, std::size_t arc_idx, const Subarc& target,
                 SourceOffset source_x_offset = SOURCE_OFFSET_NONE) const override;

private:
    const UpPhase* up_phase_;
    const Submap* input_submap_;
    const Polygon* input_curve_;
    std::size_t merge_grade_;
};

class UpPhaseArcCutter final : public ArcCuttingOracle {
public:
    UpPhaseArcCutter(const UpPhase& up, const Submap& input_submap, const Polygon& input_curve,
                     std::size_t merge_grade);

    std::vector<ArcPiece> cut(std::size_t arc_idx, const Subarc& target) const override;

    static std::size_t piece_count_bound(std::size_t merge_grade) noexcept {
        return 2 + 6 * merge_grade;
    }

    static std::size_t piece_granularity_bound(std::size_t merge_grade) noexcept {
        return UpPhase::grade_granularity(ceil_beta(merge_grade));
    }

private:
    const UpPhase* up_phase_;
    const Submap* input_submap_;
    const Polygon* input_curve_;
    std::size_t merge_grade_;
};

}
