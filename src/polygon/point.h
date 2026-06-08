#pragma once

// [C91 §2.1]: Vertex of a simple polygonal curve + ∂C side marker.

#include <cstddef>

namespace chazelle {

// [C91 §2.1]: "left and right sides of a snake" — which side of ∂C a point is on.
enum Side : unsigned char { LEFT, RIGHT };

// 2D vertex; `index` is the SoS tie-break tag ([C91 §2 tex 47,] [10]/[31]).
struct Point {
    double x = 0.0;
    double y = 0.0;
    std::size_t index = 0;
};

} // namespace chazelle
