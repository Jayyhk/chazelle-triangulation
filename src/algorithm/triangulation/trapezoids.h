#pragma once

#include "../polygon/polygon.h"
#include "../submap/submap.h"

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace chazelle {

struct Trapezoid {
    std::size_t top_vertex = NONE;
    std::size_t bottom_vertex = NONE;
    std::size_t left_edge = NONE;
    std::size_t right_edge = NONE;
};

struct VertexTrapezoids {
    std::array<std::size_t, 2> trapezoids = {NONE, NONE};
    std::size_t count = 0;
};

struct TrapezoidExtractionWork {
    std::size_t boundary_edges = 0;
    std::size_t regions = 0;
    std::size_t arcs = 0;
    std::size_t arc_edges = 0;
    std::size_t chords = 0;
    std::size_t joined_regions = 0;
};

struct TrapezoidDecomposition {
    std::vector<Trapezoid> trapezoids;
    std::vector<VertexTrapezoids> vertex_trapezoids;
    TrapezoidExtractionWork work;
};

TrapezoidDecomposition extract_trapezoids(const Polygon& curve, const Submap& visibility_map,
                                          std::span<const std::size_t> original_vertex_positions);

TrapezoidDecomposition compute_trapezoid_decomposition(const std::vector<Point>& vertices);

}
