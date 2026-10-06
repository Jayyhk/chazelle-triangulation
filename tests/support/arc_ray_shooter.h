#pragma once

#include "algorithm/merge/oracle.h"

inline chazelle::Subarc test_full_subarc(const chazelle::Polygon& curve, std::size_t first,
                                         chazelle::Side side, std::size_t last) {
    using namespace chazelle;
    return Subarc{first,
                  side,
                  last,
                  side,
                  symbolic_y_of(curve.vertex(side == LEFT ? first : first + 1)),
                  symbolic_y_of(curve.vertex(side == LEFT ? last + 1 : last))};
}

inline void test_set_cut_endpoints(const chazelle::Polygon& curve, const chazelle::Subarc& target,
                                   std::vector<chazelle::ArcPiece>& pieces) {
    for (auto& p : pieces) {
        auto full =
            test_full_subarc(curve, p.subarc.first_edge, p.subarc.first_side, p.subarc.last_edge);
        p.subarc.first_y = full.first_y;
        p.subarc.last_y = full.last_y;
    }
    assert(!pieces.empty());
    pieces.front().subarc.first_y = target.first_y;
    pieces.back().subarc.last_y = target.last_y;
}

struct TestArcRayShooter : chazelle::RayShootingOracle {
    const chazelle::Polygon* input_curve;
    explicit TestArcRayShooter(const chazelle::Polygon* c) : input_curve(c) {}
    TestArcRayShooter(const chazelle::Submap&, const chazelle::Polygon& c, std::size_t)
        : input_curve(&c) {}

    static bool crossing_x(const chazelle::Polygon& curve, std::size_t e,
                           const chazelle::SymbolicY& sy, chazelle::Exact* x) {
        return chazelle::edge_crossing_x(curve, e, sy, x);
    }

    chazelle::RayHit
    shoot(chazelle::Point p, chazelle::Side dir, std::size_t, const chazelle::Subarc& target,
          chazelle::SourceOffset source_offset = chazelle::SOURCE_OFFSET_NONE) const override {
        using namespace chazelle;
        const SymbolicY sy{p.y, p.index};
        ArcSideRange ranges[3];
        const std::size_t nl =
            subarc_side_ranges(target, 0, input_curve->num_vertices() - 1, ranges);
        RayHit best;
        chazelle::Exact best_d = 0.0;
        for (std::size_t g = 0; g < nl; ++g) {
            for (std::size_t e = ranges[g].first_edge; e <= ranges[g].last_edge; ++e) {
                chazelle::Exact x;
                if (!edge_crossing_x(*input_curve, e, sy, &x))
                    continue;

                if (!subarc_contains_point(target, *input_curve, e, ranges[g].side, sy, 0,
                                           input_curve->num_vertices() - 1))
                    continue;
                const auto ed = input_curve->edge(e);
                const bool asc = point_y_below(input_curve->vertex(ed.start_idx),
                                               input_curve->vertex(ed.end_idx));
                const Side minus_x = asc ? LEFT : RIGHT;
                const Side struck = dir == RIGHT ? minus_x : (minus_x == LEFT ? RIGHT : LEFT);
                const chazelle::Exact d = dir == RIGHT ? x - p.x : p.x - x;
                const bool wrapped =
                    d < 0.0 ||
                    (d == 0.0 && !perturbed_hit_forward(*input_curve, sy, dir, source_offset, e));
                const bool better =
                    !best.hit ||
                    (wrapped != best.wrapped
                         ? !wrapped
                         : (d != best_d ? d < best_d
                                        : ray_contact_precedes(*input_curve, sy, dir, e, struck,
                                                               best.edge, best.side)));
                if (better) {
                    best = RayHit{true, x, p.y, e, struck, wrapped, NONE};
                    best_d = d;
                }
            }
        }
        return best;
    }
};
