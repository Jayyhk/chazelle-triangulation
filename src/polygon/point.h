#pragma once

// [C91 §2.1]: Vertex of a simple polygonal curve + ∂C side marker.

#include <cstddef>

namespace chazelle {

// [C91 §2.1]: "left and right sides of a snake" — which side of ∂C a point is on.
enum Side : unsigned char { LEFT, RIGHT };

// 2D vertex.
struct Point {
    double x = 0.0;
    double y = 0.0;
    // Dual role, one value: the vertex's position in P's input table
    // AND its SoS tie-break tag ([C91 §2 tex 47], [10]/[31] — larger
    // tag ⟹ lower perturbed y).  These coincide by construction:
    // tags are assigned in table order, so raw-y ties between
    // distinct vertices, padding clones ([C91 §4 tex 316]), and
    // duplicated companions all order by table position.  When a
    // Point is a derived location (a ray source or crossing rather
    // than a table vertex), `index` carries the SoS tag of the LEVEL
    // it rides, not a table position.
    std::size_t index = 0;
};

} // namespace chazelle
