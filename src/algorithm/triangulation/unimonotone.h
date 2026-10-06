#pragma once

#include "trapezoids.h"

#include <cstddef>
#include <span>
#include <vector>

namespace chazelle {

struct TrapezoidDiagonal {
    std::size_t top_vertex = NONE;
    std::size_t bottom_vertex = NONE;
    std::size_t trapezoid = NONE;
};

struct UnimonotonePolygon {
    std::vector<std::size_t> vertices;
    std::size_t top_vertex = NONE;
    std::size_t bottom_vertex = NONE;
};

struct UnimonotoneDecompositionWork {
    std::size_t vertices_initialized = 0;
    std::size_t vertex_visits = 0;
    std::size_t trapezoids_examined = 0;
    std::size_t emitted_vertices = 0;
    std::size_t maximum_stack_size = 0;
};

struct UnimonotoneDecomposition {
    std::vector<UnimonotonePolygon> polygons;
    std::vector<TrapezoidDiagonal> diagonals;
    UnimonotoneDecompositionWork work;
};

UnimonotoneDecomposition decompose_unimonotone(std::span<const Point> vertices,
                                               const TrapezoidDecomposition& trapezoids);

UnimonotoneDecomposition compute_unimonotone_decomposition(const std::vector<Point>& vertices);

}
