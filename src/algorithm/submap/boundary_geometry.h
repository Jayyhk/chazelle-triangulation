#pragma once

#include "../polygon/polygon.h"
#include "submap.h"

#include <array>
#include <cassert>
#include <cstddef>

namespace chazelle {

struct RegionArcs {
    static constexpr std::size_t MAX = 4;
    std::array<std::size_t, MAX> arcs = {};
    std::size_t count = 0;
    void push(std::size_t arc_idx) {
        assert(count < MAX && "[C91 §2.3 tex 114]: conformal region has ≤ 4 arcs");
        arcs[count++] = arc_idx;
    }
    const std::size_t* begin() const {
        return arcs.data();
    }
    const std::size_t* end() const {
        return arcs.data() + count;
    }
};

RegionArcs collect_region_arcs(const Submap& submap, std::size_t region);

Side shooting_direction(std::size_t edge, Side side, const Polygon& curve);

bool chord_runs_through_infinity(const Polygon& curve, const Chord& c);

bool arc_starts_at_chord_slot(const Submap& submap, const Polygon& curve, const Chord& c,
                              bool left_slot, std::size_t arc_idx);

}
