#pragma once

#include <cstddef>

namespace chazelle {

/// Which side of the double boundary an arc or point lies on.
/// The double boundary ∂C is a conceptual thickening of curve C —
/// it is never constructed as a data structure. This flag tracks the
/// side implicitly (Chazelle §2.3: "the notion of double boundary
/// need not be encoded explicitly").
enum class Side : unsigned char { LEFT, RIGHT };

/// A 2D point with an index into the polygon's vertex array.
/// The index is used by symbolic y-perturbation to break ties.
struct Point {
    double x = 0.0;
    double y = 0.0;
    std::size_t index = 0; ///< Position in the polygon's input table.
};

} // namespace chazelle
