#pragma once

/// [C91 §2.4]: Arc-structure.
///
/// "Each arc is represented by a separate arc-structure."
/// Points directly into the input table (Polygon's edge array).
///
/// (ii): "each arc-structure has a pointer to the node of the tree
/// whose corresponding region is incident upon the arc."

#include "../polygon/point.h"
#include "../polygon/perturbation.h"
#include "../common.h"

#include <cstddef>

namespace chazelle {

struct Arc {
    /// Index of the first edge (e₁) in the input table.
    std::size_t first_edge = NONE;
    /// Which side of the double boundary for e₁.
    Side first_side = LEFT;

    /// Index of the last edge (eₜ) in the input table.
    /// For a single-edge arc (t=1): last_edge == first_edge.
    std::size_t last_edge = NONE;
    /// Which side of the double boundary for eₜ.
    Side last_side = LEFT;

    /// [C91 §2.4(ii)]: "a pointer to the node of the tree whose
    /// corresponding region is incident upon the arc."
    std::size_t region_node = NONE;

    /// [C91 §2.2]: Cached nonnull edge count for weight computation.
    std::size_t edge_count = 0;

    /// [C91 §2.4]: Symbolic y-coordinate of the arc's starting
    /// position on ∂C.  Used for disambiguation in double_identify
    /// when multiple arcs share the same edge ("sorting the endpoints
    /// falling within the same edges by considering y-coordinates").
    double key_y = 0.0;
    std::size_t key_y_tag = NONE;

    SymbolicY key_symbolic_y() const noexcept {
        return {key_y, key_y_tag};
    }
};

} // namespace chazelle
