#pragma once

#include "polygon/polygon.h"

#include <algorithm>
#include <cassert>
#include <vector>

namespace chazelle::test {

inline std::vector<std::vector<Point>> polygon_fixtures() {
    return {
        {{0, 0, 0}, {4, 0, 1}, {2, 3, 2}},
        {{0, 0, 0}, {4, 0, 1}, {4, 4, 2}, {0, 4, 3}},
        {{0, 0, 0}, {2, 0, 1}, {4, 0, 2}, {4, 2, 3}, {4, 4, 4}, {0, 4, 5}},
        {{0, 0, 0}, {6, 0, 1}, {6, 6, 2}, {4, 6, 3}, {4, 2, 4}, {2, 2, 5}, {2, 6, 6}, {0, 6, 7}},
        {{0, -3, 0},
         {0, 2, 1},
         {1, 8, 2},
         {9, 9, 3},
         {8, 3, 4},
         {8, -2, 5},
         {6, 0, 6},
         {4, 4, 7},
         {2, 1, 8}},
        {{0, 8, 0}, {6, 6, 1}, {4, 4, 2}, {1, 2, 3}, {0, 0, 4}},
        {{0, 0, 0}, {6, 0, 1}, {6, 1, 2}, {5, 1, 3}, {5, 2, 4}, {4, 2, 5}, {4, 3, 6}, {0, 3, 7}}};
}

inline std::vector<Point> boundary_order(std::vector<Point> vertices, bool reverse,
                                         std::size_t first, std::size_t first_tag = 41) {
    assert(first < vertices.size() && vertices.size() < SOS_NONE - first_tag);
    if (reverse)
        std::reverse(vertices.begin(), vertices.end());
    std::rotate(vertices.begin(), vertices.begin() + static_cast<std::ptrdiff_t>(first),
                vertices.end());
    for (std::size_t vertex = 0; vertex < vertices.size(); ++vertex)
        vertices[vertex].index = first_tag + vertex;
    return vertices;
}

}
