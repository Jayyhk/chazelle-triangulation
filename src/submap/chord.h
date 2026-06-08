#pragma once

// [C91 §2.4(ii)]: Normal-form chord — "a record describing the chord plus
// pointers to the arc-structures of the two, three, or four arcs adjacent
// to it."

#include "../polygon/point.h"
#include "../polygon/perturbation.h"
#include "../common.h"

#include <cstddef>

namespace chazelle {

struct Chord {
    // [§2.4(i)]: The two regions this chord separates.
    std::size_t region[2] = {NONE, NONE};

    // [§2.4(ii)]: Per-endpoint adj arcs.  Count is 1 at polygon-vertex
    // endpoints, 2 at non-vertex endpoints (§2.2 tex 94: non-vertex
    // endpoints merge on chord removal).
    struct AdjArcs {
        std::size_t arcs[2] = {NONE, NONE};
        std::size_t count = 0;
    };
    AdjArcs left_adj;
    AdjArcs right_adj;

    // Endpoint positions on ∂C: edge index in the input table + ∂C side.
    std::size_t left_edge  = NONE;
    std::size_t right_edge = NONE;
    Side left_side  = LEFT;
    Side right_side = RIGHT;

    // Symbolic y of the (horizontal) chord.
    double y     = 0.0;
    std::size_t y_tag = NONE;

    SymbolicY symbolic_y() const noexcept { return {y, y_tag}; }

    // True iff both endpoints sit at the same ∂C point.
    bool is_null_length = false;

    // [§3.3]: Tombstone for O(1) removal; stripped by compact().
    bool dead = false;
};

} // namespace chazelle
