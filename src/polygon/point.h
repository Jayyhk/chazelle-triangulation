#pragma once

/// [C91 §2]: Vertex of a simple polygonal curve.

#include <cstddef>

namespace chazelle {

/// A 2D vertex with an index used as the SoS tie-break tag.
/// [C91 §2]: "no two distinct vertices have the same y-coordinate"
/// — enforced symbolically via Edelsbrunner [10] / Yap [31].
struct Point {
    double x = 0.0;
    double y = 0.0;
    std::size_t index = 0;
};

} // namespace chazelle
