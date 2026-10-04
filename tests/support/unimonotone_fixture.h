#pragma once

#include "triangulation/trapezoids.h"

#include <cstddef>
#include <vector>

namespace chazelle::test {

struct TrapezoidizedPolygon {
    std::vector<Point> vertices;
    TrapezoidDecomposition trapezoids;
};

inline TrapezoidizedPolygon alternating_chains(std::size_t chain_size) {
    TrapezoidizedPolygon polygon;
    const std::size_t count = 2 * chain_size + 2;
    polygon.vertices.reserve(count);
    polygon.vertices.push_back({0, 0, 0});
    for (std::size_t vertex = 1; vertex <= chain_size; ++vertex)
        polygon.vertices.push_back({-1, Exact(2 * vertex - 1), vertex});
    polygon.vertices.push_back({0, Exact(2 * chain_size + 1), chain_size + 1});
    for (std::size_t vertex = chain_size; vertex > 0; --vertex)
        polygon.vertices.push_back({1, Exact(2 * vertex), polygon.vertices.size()});
    polygon.trapezoids.vertex_trapezoids.resize(count);
    const auto vertex_at_height = [count, chain_size](std::size_t height) {
        if (height == 0)
            return std::size_t{0};
        if (height == 2 * chain_size + 1)
            return chain_size + 1;
        return height % 2 == 0 ? count - height / 2 : (height + 1) / 2;
    };
    polygon.trapezoids.trapezoids.reserve(count - 1);
    for (std::size_t height = 1; height < count; ++height) {
        const std::size_t top = vertex_at_height(height);
        const std::size_t bottom = vertex_at_height(height - 1);
        const std::size_t index = polygon.trapezoids.trapezoids.size();
        polygon.trapezoids.trapezoids.push_back(
            {top, bottom, height / 2, count - 1 - (height - 1) / 2});
        polygon.trapezoids.vertex_trapezoids[top] = {{{index, NONE}}, 1};
    }
    return polygon;
}

}
