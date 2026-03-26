#pragma once

/// [C91 §2.1]: Vertex of a simple polygonal curve and ∂C side marker.

#include <cstddef>

namespace chazelle {

/// [C91 §2.1]: "the left and right sides of a snake."
/// Which side of the double boundary ∂C a point lies on.
enum Side : unsigned char { LEFT, RIGHT };

/// A 2D vertex with an index used as the SoS tie-break tag.
/// [C91 §2]: "no two distinct vertices have the same y-coordinate"
/// — enforced symbolically via Edelsbrunner [10] / Yap [31].
struct Point {
    double x = 0.0;
    double y = 0.0;
    std::size_t index = 0;
};

} // namespace chazelle
