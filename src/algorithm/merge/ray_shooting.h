#pragma once

#include "../polygon/perturbation.h"
#include "../polygon/point.h"
#include "../polygon/polygon.h"
#include "../separator/planar_separator.h"
#include "../submap/submap.h"
#include "oracle.h"

#include <cstddef>
#include <vector>

namespace chazelle {

class RayShootingStructure {
public:
    RayShootingStructure(const Submap& submap, const Polygon& curve, std::size_t granularity);

    RayHit shoot_toward_boundary(const Point& p, Side direction,
                                 const SourceOffset& source_x_offset = SOURCE_OFFSET_NONE) const;

    std::size_t num_faces() const noexcept {
        return face_count_;
    }

    std::size_t face_of_region(std::size_t r) const {
        return face_of_region_[r];
    }
    std::size_t region_of_face(std::size_t f) const {
        return region_of_face_[f];
    }

    const std::vector<std::pair<std::size_t, std::size_t>>& dual_edges() const noexcept {
        return dual_edges_;
    }
    const SeparatorDecomposition& decomposition() const noexcept {
        return separator_decomposition_;
    }

    struct LineCrossing {
        std::size_t chord = NONE;
        SymbolicY y{};
        std::size_t region_below = NONE;
        std::size_t region_above = NONE;
    };
    const std::vector<LineCrossing>& vertical_line() const noexcept {
        return vertical_line_crossings_;
    }

    std::size_t region_at_infinity() const noexcept {
        return region_at_infinity_;
    }

private:
    const Submap* submap_;
    const Polygon* curve_;
    std::size_t granularity_;

    std::size_t face_count_ = 0;
    std::vector<std::size_t> face_of_region_;
    std::vector<std::size_t> region_of_face_;
    std::vector<std::vector<std::size_t>> arcs_of_region_;
    std::vector<std::pair<std::size_t, std::size_t>> dual_edges_;
    SeparatorDecomposition separator_decomposition_;
    std::vector<std::size_t> separator_faces_;

    std::vector<std::vector<std::size_t>> subset_faces_;
    std::vector<LineCrossing> vertical_line_crossings_;
    std::size_t region_at_infinity_ = NONE;

    struct BoundaryInterval {
        std::size_t lo_edge = NONE, hi_edge = NONE;
        SymbolicY lo_y{}, hi_y{};
        std::size_t region = NONE;
    };
    std::vector<BoundaryInterval> left_intervals_, right_intervals_;

    void regions_at_boundary(std::size_t edge, Side side, const SymbolicY& y,
                             std::vector<std::size_t>& out) const;

    void build_faces();
    void build_dual_graph_and_decomposition();
    void build_vertical_line();
};

}
