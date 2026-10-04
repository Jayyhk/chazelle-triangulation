#pragma once

#include "../common.h"
#include "../polygon/perturbation.h"
#include "../polygon/point.h"

#include <cstddef>

namespace chazelle {

struct Chord {
    std::size_t region[2] = {NONE, NONE};

    struct AdjArcs {
        std::size_t arcs[2] = {NONE, NONE};
        std::size_t count = 0;
    };
    AdjArcs left_adj;
    AdjArcs right_adj;

    std::size_t left_edge = NONE;
    std::size_t right_edge = NONE;
    Side left_side = LEFT;
    Side right_side = RIGHT;

    Exact y = 0.0;
    std::size_t y_tag = NONE;

    SymbolicY symbolic_y() const noexcept {
        return {y, y_tag};
    }

    bool is_null_length = false;

    bool dead = false;
};

}
