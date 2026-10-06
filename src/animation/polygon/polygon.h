#pragma once

#include "../common.h"
#include "edge.h"
#include "perturbation.h"
#include "point.h"

#include <cassert>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

namespace chazelle::animation {

class Polygon {
public:
    explicit Polygon(std::vector<Point> vertices);

    static Polygon symbolic_curve(std::vector<Point> vertices, Exact horizontal_shift = 0);

    Polygon(const Polygon& c1, const Polygon& c2);

    Polygon subchain(std::size_t first, std::size_t count) const;

    Polygon subchain(std::size_t first, std::size_t count, std::size_t minimum,
                     std::size_t maximum) const;

    Polygon reversed() const;

    Exact edge_horizontal_shift(std::size_t edge) const;

    Exact edge_horizontal_delta(std::size_t edge) const {
        return vertex(edge + 1).x - vertex(edge).x + edge_horizontal_shift(edge);
    }

    bool previous_branch_left(std::size_t vertex) const;

    std::size_t local_index_of_tag(std::size_t tag) const noexcept;

    std::size_t num_vertices() const noexcept {
        return len_;
    }
    std::size_t num_edges() const noexcept {
        return len_ - 1;
    }

    const Point& vertex(std::size_t i) const noexcept {
        assert(i < len_);
        if (pieces_.empty())
            return table_->vertices[offset_ + i];
        for (const Piece& piece : pieces_) {
            if (i < piece.count) {
                const std::size_t position = piece.reversed ? piece.count - 1 - i : i;
                return piece.table->vertices[piece.offset + position];
            }
            i -= piece.count - 1;
        }
        assert(false && "[C91 §4.2 tex 367]: region boundary pieces cover the curve");
        return pieces_.back().table->vertices[pieces_.back().offset];
    }

    Edge edge(std::size_t i) const noexcept {
        assert(i + 1 < len_);
        return Edge{i, i + 1};
    }

    bool is_endpoint(std::size_t vertex_index) const noexcept {
        assert(vertex_index < len_ && "[C91 §2.1]: invalid vertex");
        return vertex_index == 0 || vertex_index == len_ - 1;
    }

    bool is_y_extremum(std::size_t vertex_index) const noexcept {
        assert(vertex_index < len_ && "[C91 §2.1]: invalid vertex");
        if (is_endpoint(vertex_index))
            return false;
        return is_local_y_extremum(vertex(vertex_index - 1), vertex(vertex_index),
                                   vertex(vertex_index + 1));
    }

    std::size_t count_nonnull_edges(std::size_t lo, std::size_t hi) const noexcept;

    bool edge_is_null(std::size_t i) const noexcept {
        return count_nonnull_edges(i, i) == 0;
    }

    std::size_t next_nonnull_edge(std::size_t i) const noexcept;

    std::size_t table_offset() const noexcept {
        assert(pieces_.empty() && "[C91 §4.1]: original portions share the input table");
        return offset_;
    }

    std::size_t max_y_vertex() const noexcept {
        return max_y_vertex_;
    }
    std::size_t min_y_vertex() const noexcept {
        return min_y_vertex_;
    }

private:
    friend class AnimationTrace;
    struct Table {
        mutable std::size_t animation_session = 0;
        mutable std::size_t animation_table = NONE;
        std::vector<Point> vertices;

        std::vector<std::size_t> nonnull_prefix;

        std::vector<std::size_t> next_nonnull;

        std::vector<std::size_t> previous_nonnull;

        std::size_t tag_stride = 1;
        Exact horizontal_shift = 0;
    };

    struct Piece {
        std::shared_ptr<const Table> table;
        std::size_t offset;
        std::size_t count;
        bool reversed;
    };

    Polygon() = default;

    void find_y_extremes();

    void initialize(std::vector<Point> vertices, bool original, Exact horizontal_shift = 0);

    void append_pieces(const Polygon& curve);

    std::shared_ptr<const Table> table_;
    std::vector<Piece> pieces_;
    std::size_t offset_ = 0;
    std::size_t len_ = 0;
    std::size_t max_y_vertex_ = 0;
    std::size_t min_y_vertex_ = 0;
    mutable std::size_t animation_session_ = 0;
    mutable std::size_t animation_curve_ = NONE;
};

inline Exact edge_t_at_y(const Polygon& curve, std::size_t edge_idx, const SymbolicY& target_y) {
    const auto e = curve.edge(edge_idx);
    const Point& vs = curve.vertex(e.start_idx);
    const Point& ve = curve.vertex(e.end_idx);
    if (symbolic_y_equal(target_y, symbolic_y_of(vs)))
        return 0.0;
    if (symbolic_y_equal(target_y, symbolic_y_of(ve)))
        return 1.0;

    assert(((symbolic_y_less(symbolic_y_of(vs), target_y) &&
             symbolic_y_less(target_y, symbolic_y_of(ve))) ||
            (symbolic_y_less(symbolic_y_of(ve), target_y) &&
             symbolic_y_less(target_y, symbolic_y_of(vs)))) &&
           "[C91 §2 tex 47]: query must lie strictly inside the edge's "
           "perturbed y-range (no crossing otherwise)");

    assert(vs.y != ve.y && "[C91 §2 tex 47]: horizontal edge requires SoS tag-match at "
                           "one endpoint; strictly-between target is unreachable");
    return (target_y.y - vs.y) / (ve.y - vs.y);
}

inline Exact edge_x_at_y(const Polygon& curve, std::size_t edge_idx, const SymbolicY& target_y) {
    const auto e = curve.edge(edge_idx);
    const Point& vs = curve.vertex(e.start_idx);
    const Point& ve = curve.vertex(e.end_idx);
    Exact t = edge_t_at_y(curve, edge_idx, target_y);

    if (t == 0.0)
        return vs.x;
    if (t == 1.0)
        return ve.x;
    Exact x = vs.x + t * curve.edge_horizontal_delta(edge_idx);
    const Exact shift = curve.edge_horizontal_shift(edge_idx);
    if (shift != 0) {
        const Exact period = shift < 0 ? -shift : shift;
        if (x > period / 2)
            x -= period;
        else if (x < -period / 2)
            x += period;
    }
    return x;
}

inline bool edge_crossing_x(const Polygon& curve, std::size_t e, const SymbolicY& sy, Exact* x) {
    const auto& ed = curve.edge(e);
    const Point& vs = curve.vertex(ed.start_idx);
    const Point& ve = curve.vertex(ed.end_idx);
    SymbolicY y0 = symbolic_y_of(vs);
    SymbolicY y1 = symbolic_y_of(ve);
    if (symbolic_y_equal(sy, y0)) {
        *x = vs.x;
        return true;
    }
    if (symbolic_y_equal(sy, y1)) {
        *x = ve.x;
        return true;
    }
    bool between = (symbolic_y_less(y0, sy) && symbolic_y_less(sy, y1)) ||
                   (symbolic_y_less(y1, sy) && symbolic_y_less(sy, y0));
    if (!between)
        return false;

    assert(vs.y != ve.y && "[C91 §2 tex 47]: horizontal edge admits no strictly-interior "
                           "crossing");
    *x = edge_x_at_y(curve, e, sy);
    return true;
}

inline Exact perturbed_x_offset(const Polygon& curve, const SymbolicY& sy, std::size_t edge) {
    const auto& ed = curve.edge(edge);
    const Point& a = curve.vertex(ed.start_idx);
    const Point& b = curve.vertex(ed.end_idx);
    if (symbolic_y_equal(sy, symbolic_y_of(a)) || symbolic_y_equal(sy, symbolic_y_of(b)))
        return 0.0;
    const bool a_tied = (a.y == sy.y);
    const bool b_tied = (b.y == sy.y);

    if (!a_tied && !b_tied)
        return 0.0;
    assert(a_tied != b_tied && "[C91 §2 tex 47]: raw-horizontal edges admit no strict "
                               "crossing (consecutive endpoint tags)");
    const Point& at = a_tied ? a : b;
    const Point& aw = a_tied ? b : a;
    const Exact delta =
        a_tied ? curve.edge_horizontal_delta(edge) : -curve.edge_horizontal_delta(edge);
    const Exact slope = delta / (aw.y - at.y);
    return symbolic_y_less(sy, symbolic_y_of(at)) ? -slope : slope;
}

using SourceOffset = std::optional<Exact>;
inline constexpr auto SOURCE_OFFSET_NONE = std::nullopt;

inline bool perturbed_hit_forward(const Polygon& curve, const SymbolicY& sy, Side dir,
                                  const SourceOffset& source_x_offset, std::size_t hit_edge) {
    const Exact hit_off = perturbed_x_offset(curve, sy, hit_edge);
    const Exact source_off = source_x_offset.value_or(Exact{});
    if (source_off == hit_off)
        return false;
    return (dir == RIGHT) ? (hit_off > source_off) : (hit_off < source_off);
}

inline bool extremum_prev_branch_left(const Point& u, const Point& v, const Point& w) {
    assert(is_local_y_extremum(u, v, w) &&
           "[C91 §2.1 tex 72]: inside companion sides are defined at y-extrema only");
    Exact nu = u.x - v.x, du = u.y > v.y ? u.y - v.y : v.y - u.y;
    Exact nw = w.x - v.x, dw = w.y > v.y ? w.y - v.y : v.y - w.y;
    if (nu == 0.0 && du == 0.0)
        du = 1.0;
    if (nw == 0.0 && dw == 0.0)
        dw = 1.0;
    if (nu != 0.0 && du == 0.0) {
        nu = (nu > 0.0) ? 1.0 : -1.0;
    }
    if (nw != 0.0 && dw == 0.0) {
        nw = (nw > 0.0) ? 1.0 : -1.0;
    }
    const Exact lhs = nu * dw;
    const Exact rhs = nw * du;
    assert(lhs != rhs && "[C91 §2 tex 47]: extremum branches have distinct x-offsets "
                         "(equal offsets ⟹ overlapping edges, non-simple P)");
    return lhs < rhs;
}

inline bool is_inside_companion(const Polygon& curve, std::size_t edge, Side side,
                                std::size_t vidx) {
    if (vidx == 0 || vidx + 1 >= curve.num_vertices())
        return false;
    const Point& u = curve.vertex(vidx - 1);
    const Point& v = curve.vertex(vidx);
    const Point& w = curve.vertex(vidx + 1);
    if (!is_local_y_extremum(u, v, w))
        return false;

    const bool prev_left = curve.previous_branch_left(vidx);

    auto minus_x_face = [&](std::size_t e) -> Side {
        const auto& ed = curve.edge(e);
        bool asc = symbolic_y_less(symbolic_y_of(curve.vertex(ed.start_idx)),
                                   symbolic_y_of(curve.vertex(ed.end_idx)));
        return asc ? LEFT : RIGHT;
    };
    auto plus_x_face = [&](std::size_t e) -> Side {
        return minus_x_face(e) == LEFT ? RIGHT : LEFT;
    };

    Side inside_next = prev_left ? minus_x_face(vidx) : plus_x_face(vidx);
    Side inside_prev = prev_left ? plus_x_face(vidx - 1) : minus_x_face(vidx - 1);
    if (edge == vidx && side == inside_next)
        return true;
    if (edge == vidx - 1 && side == inside_prev)
        return true;
    return false;
}

inline bool ray_contact_precedes(const Polygon& curve, const SymbolicY& sy, Side dir,
                                 std::size_t e_new, Side s_new, std::size_t e_old, Side s_old) {
    if (e_new == e_old && s_new == s_old)
        return false;
    auto tag_matched_vertex = [&](std::size_t e) -> std::size_t {
        const auto& ed = curve.edge(e);
        if (symbolic_y_equal(sy, symbolic_y_of(curve.vertex(ed.start_idx))))
            return ed.start_idx;
        if (symbolic_y_equal(sy, symbolic_y_of(curve.vertex(ed.end_idx))))
            return ed.end_idx;
        return NONE;
    };
    const std::size_t vn = tag_matched_vertex(e_new);
    const std::size_t vo = tag_matched_vertex(e_old);
    if (vn != NONE && vo != NONE) {
        assert(vn == vo && "[C91 §2 tex 47]: a symbolic y names a unique vertex");
        return is_inside_companion(curve, e_old, s_old, vo) &&
               !is_inside_companion(curve, e_new, s_new, vn);
    }
    if (vn == NONE && vo == NONE) {
        const std::size_t vidx = std::max(e_new, e_old);
        assert(vidx == std::min(e_new, e_old) + 1 &&
               "[C91 §2.1 tex 68]: coincident strict crossings flank "
               "one shared vertex");
        const Point& vb = curve.vertex(vidx);
        assert(vb.y == sy.y && vb.index != sy.tag &&
               "[C91 §2 tex 47]: the flanked vertex shares the ray's "
               "raw level under a different tag");
        auto dxdy = [&](std::size_t e) -> Exact {
            const auto& ed = curve.edge(e);
            const Point& a = curve.vertex(ed.start_idx);
            const Point& b = curve.vertex(ed.end_idx);
            assert(a.y != b.y && "[C91 §2 tex 47]: raw-horizontal edges admit no "
                                 "strict crossing");
            return curve.edge_horizontal_delta(e) / (b.y - a.y);
        };
        const Exact dsign = symbolic_y_less(sy, symbolic_y_of(vb)) ? -1.0 : 1.0;
        const Exact cn = dxdy(e_new) * dsign;
        const Exact co = dxdy(e_old) * dsign;
        if (cn == co)
            return false;
        return (dir == RIGHT) ? (cn < co) : (cn > co);
    }

    const bool new_is_strict = (vn == NONE);
    const std::size_t es = new_is_strict ? e_new : e_old;
    const auto& ed = curve.edge(es);
    const Point& vs = curve.vertex(ed.start_idx);
    const Point& ve = curve.vertex(ed.end_idx);

    assert((vs.y == sy.y) != (ve.y == sy.y) &&
           "[C91 §2 tex 47]: exactly one endpoint of the strict edge "
           "lies at the tied raw level");
    const Point& at = (vs.y == sy.y) ? vs : ve;
    const Point& away = (vs.y == sy.y) ? ve : vs;
    assert(away.y != at.y);

    const bool below = symbolic_y_less(sy, symbolic_y_of(at));
    const Exact dx =
        (vs.y == sy.y) ? curve.edge_horizontal_delta(es) : -curve.edge_horizontal_delta(es);
    int slope_sign = (dx == 0) ? 0 : (((dx > 0) == (away.y > at.y)) ? 1 : -1);
    int corr = below ? -slope_sign : slope_sign;

    const bool strict_first = (dir == RIGHT) ? (corr <= 0) : (corr >= 0);
    return new_is_strict == strict_first;
}

}
