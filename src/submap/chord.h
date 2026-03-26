#pragma once

/// [C91 §2.4]: Normal-form chord representation.
///
/// (ii): "Each edge of the tree contains a record describing the
/// corresponding chord as well as pointers to the arc-structures
/// of the two, three, or four arcs adjacent to it."

#include "../polygon/point.h"
#include "../polygon/perturbation.h"
#include "../common.h"

#include <array>
#include <cstddef>

namespace chazelle {

struct Chord {
    /// The two regions (tree nodes) this chord separates.
    /// [C91 §2.4(i)]: "standard edge/node adjacency."
    std::size_t region[2] = {NONE, NONE};

    /// [C91 §2.4(ii)]: "pointers to the arc-structures of the two,
    /// three, or four arcs adjacent to it."
    std::array<std::size_t, 4> adj_arcs = {NONE, NONE, NONE, NONE};
    std::size_t num_adj_arcs = 0;

    /// Chord endpoint positions on ∂C.
    /// Stored as edge indices into the input table + ∂C side flags.
    std::size_t left_edge  = NONE;
    std::size_t right_edge = NONE;
    Side left_side  = LEFT;
    Side right_side = RIGHT;

    /// Symbolic y-coordinate of the chord (horizontal segment).
    double y     = 0.0;
    std::size_t y_tag = NONE;

    SymbolicY symbolic_y() const noexcept { return {y, y_tag}; }

    /// True if this is a null-length chord (both endpoints at the
    /// same ∂C position).
    bool is_null_length = false;
};

} // namespace chazelle
