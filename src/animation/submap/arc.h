#pragma once

#include "../common.h"
#include "../polygon/perturbation.h"
#include "../polygon/point.h"
#include "../polygon/polygon.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>

namespace chazelle::animation {

struct ArcSideRange {
    Side side = LEFT;
    std::size_t first_edge = NONE;
    std::size_t last_edge = NONE;
};

struct Arc {
    std::size_t first_edge = NONE;
    Side first_side = LEFT;
    std::size_t last_edge = NONE;
    Side last_side = LEFT;

    std::size_t region_node = NONE;

    std::size_t edge_count = 0;

    bool wraps() const noexcept {
        if (first_side != last_side)
            return true;

        if (first_side == LEFT)
            return first_edge > last_edge;
        return first_edge < last_edge;
    }

    bool wraps_end() const noexcept {
        if (first_side == LEFT && last_side == RIGHT)
            return true;
        return first_side == last_side && wraps();
    }

    bool wraps_start() const noexcept {
        if (first_side == RIGHT && last_side == LEFT)
            return true;
        return first_side == last_side && wraps();
    }

    std::size_t side_ranges(std::size_t c_start, std::size_t c_end,
                            ArcSideRange out[3]) const noexcept {
        assert(c_end >= 1 && "[C91 §2.1]: C has ≥ 2 vertices");
        const std::size_t last_c_edge = c_end - 1;
        if (first_side == last_side) {
            if (!wraps()) {
                out[0] = ArcSideRange{first_side, std::min(first_edge, last_edge),
                                      std::max(first_edge, last_edge)};
                return 1;
            }

            if (first_side == LEFT) {
                assert(first_edge > last_edge);
                out[0] = ArcSideRange{LEFT, first_edge, last_c_edge};
                out[1] = ArcSideRange{RIGHT, c_start, last_c_edge};
                out[2] = ArcSideRange{LEFT, c_start, last_edge};
            } else {
                assert(first_edge < last_edge);
                out[0] = ArcSideRange{RIGHT, c_start, first_edge};
                out[1] = ArcSideRange{LEFT, c_start, last_c_edge};
                out[2] = ArcSideRange{RIGHT, last_edge, last_c_edge};
            }
            return 3;
        }
        if (first_side == LEFT) {
            assert(first_edge <= last_c_edge && last_edge <= last_c_edge);
            out[0] = ArcSideRange{LEFT, first_edge, last_c_edge};
            out[1] = ArcSideRange{RIGHT, last_edge, last_c_edge};
            return 2;
        }

        assert(first_edge >= c_start && last_edge >= c_start);
        out[0] = ArcSideRange{RIGHT, c_start, first_edge};
        out[1] = ArcSideRange{LEFT, c_start, last_edge};
        return 2;
    }

    bool covers(std::size_t edge, Side side, std::size_t c_start,
                std::size_t c_end) const noexcept {
        ArcSideRange ranges[3];
        std::size_t n = side_ranges(c_start, c_end, ranges);
        for (std::size_t i = 0; i < n; ++i)
            if (ranges[i].side == side && ranges[i].first_edge <= edge &&
                edge <= ranges[i].last_edge)
                return true;
        return false;
    }

    std::pair<std::size_t, std::size_t> underlying_edge_range(std::size_t c_start,
                                                              std::size_t c_end) const noexcept {
        assert(c_end >= 1 && "wrap ranges require c_end >= 1 (C has ≥ 2 vertices)");
        if (first_side == last_side) {
            if (!wraps())
                return {std::min(first_edge, last_edge), std::max(first_edge, last_edge)};

            return {c_start, c_end - 1};
        }
        if (first_side == LEFT) {
            assert(first_edge <= c_end - 1 && last_edge <= c_end - 1 &&
                   "[C91 §2.4 tex 142]: end wrap pieces ≤ C's last edge");
            return {std::min(first_edge, last_edge), c_end - 1};
        }

        assert(first_edge >= c_start && last_edge >= c_start &&
               "[C91 §2.4 tex 142]: start wrap pieces ≥ C's first edge");
        return {c_start, std::max(first_edge, last_edge)};
    }

    bool dead = false;
};

inline std::size_t arc_boundary_edge_count(const Arc& a, const Polygon& curve, std::size_t c_start,
                                           std::size_t c_end, const SymbolicY& start_y,
                                           const SymbolicY& end_y) {
    ArcSideRange ranges[3];
    const std::size_t n = a.side_ranges(c_start, c_end, ranges);
    std::size_t total = 0;
    for (std::size_t i = 0; i < n; ++i)
        total += curve.count_nonnull_edges(ranges[i].first_edge, ranges[i].last_edge);
    if (n == 1)
        return total;

    {
        const auto& fe = curve.edge(a.first_edge);
        const SymbolicY f_exit =
            symbolic_y_of(curve.vertex(a.first_side == LEFT ? fe.end_idx : fe.start_idx));
        if (symbolic_y_equal(start_y, f_exit))
            total -= curve.count_nonnull_edges(a.first_edge, a.first_edge);
        const auto& le = curve.edge(a.last_edge);
        const SymbolicY l_entry =
            symbolic_y_of(curve.vertex(a.last_side == LEFT ? le.start_idx : le.end_idx));
        if (symbolic_y_equal(end_y, l_entry))
            total -= curve.count_nonnull_edges(a.last_edge, a.last_edge);
    }
    return total;
}

}
