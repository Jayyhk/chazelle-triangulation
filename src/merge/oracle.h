#pragma once

#include "../polygon/point.h"
#include "../polygon/polygon.h"
#include "../submap/submap.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

namespace chazelle {

struct Subarc {
    std::size_t first_edge;
    Side first_side;
    std::size_t last_edge;
    Side last_side;

    SymbolicY first_y{};
    SymbolicY last_y{};
};

inline void assert_subarc_clockwise(const Subarc&) {}

inline std::size_t subarc_side_ranges(const Subarc& s, std::size_t c_start, std::size_t c_end,
                                      ArcSideRange out[3]) {
    Arc a;
    a.first_edge = s.first_edge;
    a.first_side = s.first_side;
    a.last_edge = s.last_edge;
    a.last_side = s.last_side;
    a.region_node = 0;
    return a.side_ranges(c_start, c_end, out);
}

inline bool subarc_contains_point(const Subarc& s, const Polygon& curve, std::size_t edge,
                                  Side side, const SymbolicY& y, std::size_t c_start,
                                  std::size_t c_end) {
    ArcSideRange ranges[3];
    std::size_t n = subarc_side_ranges(s, c_start, c_end, ranges);
    auto face_trav_ascends = [&](std::size_t e_, Side s_) {
        const auto& ed = curve.edge(e_);
        bool asc = symbolic_y_less(symbolic_y_of(curve.vertex(ed.start_idx)),
                                   symbolic_y_of(curve.vertex(ed.end_idx)));
        return (s_ == LEFT) ? asc : !asc;
    };
    for (std::size_t g = 0; g < n; ++g) {
        if (ranges[g].side != side || edge < ranges[g].first_edge || edge > ranges[g].last_edge)
            continue;
        if (g == 0 && edge == s.first_edge && side == s.first_side) {
            assert(s.first_y.tag != SOS_NONE &&
                   "[C91 §3.0(i) tex 169]: shot subarcs carry exact endpoints");
            bool ok = face_trav_ascends(edge, side) ? symbolic_y_geq(y, s.first_y)
                                                    : symbolic_y_leq(y, s.first_y);
            if (!ok)
                continue;
        }
        if (g + 1 == n && edge == s.last_edge && side == s.last_side) {
            assert(s.last_y.tag != SOS_NONE &&
                   "[C91 §3.0(i) tex 169]: shot subarcs carry exact endpoints");
            bool ok = face_trav_ascends(edge, side) ? symbolic_y_leq(y, s.last_y)
                                                    : symbolic_y_geq(y, s.last_y);
            if (!ok)
                continue;
        }
        return true;
    }
    return false;
}

struct RayHit {
    bool hit = false;
    Exact x = 0.0;
    Exact y = 0.0;
    std::size_t edge = 0;
    Side side = LEFT;

    bool wrapped = false;

    std::size_t hit_arc_idx = NONE;
};

struct ArcPiece {
    Subarc subarc;
    const Submap* submap = nullptr;
    const Polygon* curve = nullptr;
    bool is_boundary_piece = false;

    std::size_t granularity = 0;
};

struct RayShootingOracle {
    virtual ~RayShootingOracle() = default;

    virtual RayHit shoot(Point p, Side direction, std::size_t arc_idx, const Subarc& target,
                         SourceOffset source_x_offset = SOURCE_OFFSET_NONE) const = 0;
};

struct ArcCuttingOracle {
    virtual ~ArcCuttingOracle() = default;

    virtual std::vector<ArcPiece> cut(std::size_t arc_idx, const Subarc& target) const = 0;
};

inline void assert_cut_postconditions([[maybe_unused]] const Polygon& input_curve,
                                      [[maybe_unused]] const Subarc& target, const ArcPiece* pieces,
                                      std::size_t count, [[maybe_unused]] std::size_t max_pieces,
                                      [[maybe_unused]] std::size_t h_gamma) {
    assert(count >= 1 && "[C91 §3.0(ii) tex 170]: cut() must produce ≥1 piece");
    assert(count <= max_pieces && "[C91 §3.0(ii) tex 170]: cut() must produce ≤ g(γᵢ) pieces");

    assert(pieces[0].subarc.first_edge == target.first_edge &&
           pieces[0].subarc.first_side == target.first_side &&
           "[C91 §3.0(ii) tex 170]: first piece must start at α'.first");
    assert(pieces[count - 1].subarc.last_edge == target.last_edge &&
           pieces[count - 1].subarc.last_side == target.last_side &&
           "[C91 §3.0(ii) tex 170]: last piece must end at α'.last");

    for (std::size_t j = 0; j < count; ++j) {
        const ArcPiece& p = pieces[j];

        assert(p.subarc.first_side == p.subarc.last_side &&
               "[C91 §3.0(ii)(2) tex 170]: each piece must lie on one side of C");

        if (p.subarc.first_side == LEFT) {
            assert(p.subarc.first_edge <= p.subarc.last_edge &&
                   "[C91 §3.0(ii)(1) tex 170]: LEFT piece must ascend in edge");
        } else {
            assert(p.subarc.first_edge >= p.subarc.last_edge &&
                   "[C91 §3.0(ii)(1) tex 170]: RIGHT piece must descend in edge");
        }

        const bool at_endpoint = (j == 0 || j + 1 == count);
        if (!at_endpoint) {
            assert(!p.is_boundary_piece &&
                   "[C91 §3.0(ii)(3) tex 170]: only first/last pieces may be boundary");
        }

        if (!p.is_boundary_piece) {
            assert(p.submap != nullptr &&
                   "[C91 §3.0(ii)(3) tex 170]: non-boundary piece requires a submap");
            assert(p.curve != nullptr &&
                   "[C91 §3.0(ii)(3) tex 170]: non-boundary piece requires its curve");

            {
                [[maybe_unused]] std::size_t lo = std::min(p.subarc.first_edge, p.subarc.last_edge);
                [[maybe_unused]] std::size_t hi = std::max(p.subarc.first_edge, p.subarc.last_edge);
                assert(p.curve->num_edges() == hi - lo + 1 &&
                       "[C91 §3.0(ii)(3) tex 170]: non-boundary piece is "
                       "vertex-to-vertex; p.curve must cover exactly the "
                       "piece's polygon edge range");

                assert(p.curve->vertex(0).index == input_curve.vertex(lo).index &&
                       p.curve->vertex(p.curve->num_vertices() - 1).index ==
                           input_curve.vertex(hi + 1).index &&
                       "[C91 §3.0(ii)(3) tex 170]: ᾱⱼ must be the "
                       "vertex-to-vertex subchain of Cᵢ at the piece's "
                       "edge range");
            }

            assert(!p.submap->tree_decomposition().empty() &&
                   "[C91 §2.4(iv) tex 139]: normal-form conformal submap "
                   "needs its tree decomposition");

            assert(p.granularity >= 1 && p.granularity <= h_gamma &&
                   "[C91 §4.1 tex 343]: piece granularity is at most h(γᵢ)");
#ifdef CHAZELLE_EXPENSIVE_ASSERTS
            p.submap->check_invariants(*p.curve);
            assert(p.submap->is_conformal() &&
                   "[C91 §3.0(ii)(3) tex 170]: non-boundary piece must be conformal");
            assert(p.submap->is_granular(p.granularity, *p.curve) &&
                   "[C91 §3.0(ii)(3) tex 170]: non-boundary piece must be "
                   "γⱼ-granular for its declared γⱼ");
            assert(p.submap->is_semigranular(h_gamma) &&
                   "[C91 §3.2 tex 248]: O(h(γᵢ)) vertices per piece region "
                   "— weights bounded by h(γᵢ)");
#endif
        } else {
            assert(p.subarc.first_edge == p.subarc.last_edge &&
                   "[C91 §3.0(ii)(3) tex 170]: boundary piece must be single-edge");
        }
    }

#ifndef NDEBUG
    assert(symbolic_y_equal(pieces[0].subarc.first_y, target.first_y) &&
           symbolic_y_equal(pieces[count - 1].subarc.last_y, target.last_y));
    for (std::size_t j = 0; j < count; ++j) {
        const Subarc& s = pieces[j].subarc;
        assert(s.first_y.tag != SOS_NONE && s.last_y.tag != SOS_NONE);
        const auto& first = input_curve.edge(s.first_edge);
        const auto& last = input_curve.edge(s.last_edge);
        const SymbolicY f0 = symbolic_y_of(input_curve.vertex(first.start_idx));
        const SymbolicY f1 = symbolic_y_of(input_curve.vertex(first.end_idx));
        const SymbolicY l0 = symbolic_y_of(input_curve.vertex(last.start_idx));
        const SymbolicY l1 = symbolic_y_of(input_curve.vertex(last.end_idx));
        auto between = [](const SymbolicY& y, const SymbolicY& a, const SymbolicY& b) {
            return (symbolic_y_leq(a, y) && symbolic_y_leq(y, b)) ||
                   (symbolic_y_leq(b, y) && symbolic_y_leq(y, a));
        };
        assert(between(s.first_y, f0, f1) && between(s.last_y, l0, l1));
        if (!pieces[j].is_boundary_piece) {
            assert(symbolic_y_equal(s.first_y, s.first_side == LEFT ? f0 : f1));
            assert(symbolic_y_equal(s.last_y, s.last_side == LEFT ? l1 : l0));
        }
        if (s.first_edge == s.last_edge) {
            const bool up = (symbolic_y_less(f0, f1) == (s.first_side == LEFT));
            assert(up ? symbolic_y_leq(s.first_y, s.last_y) : symbolic_y_leq(s.last_y, s.first_y));
        }
        if (j + 1 == count)
            continue;
        const Subarc& next = pieces[j + 1].subarc;
        assert(symbolic_y_equal(s.last_y, next.first_y));
        if (s.last_edge == next.first_edge && s.last_side == next.first_side)
            continue;
        const std::size_t lo = std::min(s.last_edge, next.first_edge);
        const std::size_t hi = std::max(s.last_edge, next.first_edge);
        if (s.last_side == next.first_side) {
            assert(hi == lo + 1 &&
                   symbolic_y_equal(s.last_y, symbolic_y_of(input_curve.vertex(hi))) &&
                   ((s.last_side == LEFT) == (s.last_edge < next.first_edge)));
        } else {
            assert(lo == hi && (lo == 0 || hi + 1 == input_curve.num_edges()) &&
                   symbolic_y_equal(
                       s.last_y, symbolic_y_of(input_curve.vertex(
                                     s.last_side == LEFT ? input_curve.num_vertices() - 1 : 0))));
        }
    }
#endif
}

}
