#pragma once

#include "unimonotone.h"

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace chazelle {

struct Triangle {
    std::array<std::size_t, 3> vertices;
};

struct TriangulationWork {
    std::size_t polygon_vertices = 0;
    std::size_t convexity_tests = 0;
    std::size_t forward_steps = 0;
    std::size_t backward_steps = 0;
    std::size_t removed_vertices = 0;
    std::size_t triangle_incidents = 0;
};

struct Triangulation {
    std::vector<Triangle> triangles;
    std::vector<std::vector<std::size_t>> vertex_triangles;
    TriangulationWork work;
};

Triangulation triangulate_unimonotone(std::span<const Point> vertices,
                                      const UnimonotoneDecomposition& decomposition);

Triangulation triangulate_trapezoids(std::span<const Point> vertices,
                                     const TrapezoidDecomposition& trapezoids);

Triangulation triangulate_polygon(const std::vector<Point>& vertices);

}
