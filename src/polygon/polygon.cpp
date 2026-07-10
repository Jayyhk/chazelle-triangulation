// src/polygon/polygon.cpp

#include "polygon.h"

#include <cassert>
#include <cmath>
#include <unordered_set>

namespace chazelle {

Polygon::Polygon(std::vector<Point> vertices) {
    assert(vertices.size() >= 2 && "[C91 §2.1]: curve needs ≥ 2 vertices");

    // [C91 §2 SoS / Edelsbrunner [10]]: vertex indices must be distinct.
    // This also implies "nonclosed" ([C91 §2.1]): front.index != back.index
    // ⟹ symbolically distinct, regardless of geometric coords.
    // [C91 §2 tex 66]: spherical-plane reference point at (0, +∞) — the
    // curve must avoid that point (Jordan-curve orientation argument).
    {
        std::unordered_set<std::size_t> seen;
        for (const auto& v : vertices) {
            assert(seen.insert(v.index).second &&
                   "[C91 §2 (SoS)]: vertex indices must be distinct");
            assert(!(v.x == 0.0 && std::isinf(v.y) && v.y > 0.0) &&
                   "[C91 §2 tex 66]: curve must avoid the reference point (0, +∞)");
        }
        // [C91 §2.4 tex 133]: every curve is a contiguous subchain of the
        // input P (whose vertices carry consecutive SoS indices v₁..vₙ),
        // so a curve's indices are consecutive from vertices[0].index.
        // Callers rely on this to translate SoS tags → table positions.
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            assert(vertices[i].index == vertices[0].index + i &&
                   "[C91 §2.4 tex 133]: curve must be a contiguous "
                   "subchain of P (consecutive SoS indices)");
        }
    }

    auto t = std::make_shared<Table>();
    t->vertices = std::move(vertices);

    // [C91 §2.2]: an edge has nonzero length iff its endpoints differ in
    // at least one coordinate.  One prefix table for the whole input
    // table; every subchain view reads it with offset arithmetic.
    const std::size_t m = t->vertices.size() - 1;
    t->nonnull_prefix.resize(m + 1, 0);
    for (std::size_t i = 0; i < m; ++i) {
        const auto& s = t->vertices[i];
        const auto& e = t->vertices[i + 1];
        t->nonnull_prefix[i + 1] = t->nonnull_prefix[i]
            + ((s.x != e.x || s.y != e.y) ? 1 : 0);
    }

    // Nonnull-successor array — one O(n) backward pass at input setup
    // ([C91 §2.2 tex 106]: weights count nonnull edges only, so scans
    // must be able to step over null runs in O(1)).
    t->next_nonnull.resize(m, m);
    for (std::size_t i = m; i-- > 0; ) {
        t->next_nonnull[i] =
            (t->nonnull_prefix[i + 1] > t->nonnull_prefix[i])
                ? i
                : (i + 1 < m ? t->next_nonnull[i + 1] : m);
    }

    table_ = std::move(t);
    offset_ = 0;
    len_ = table_->vertices.size();
    find_y_extremes();
}

Polygon::Polygon(const Polygon& c1, const Polygon& c2) {
    // [C91 §2.4 tex 133]: "The input table is read-only: it is never to
    // be modified or even copied" — the union is a wider view into the
    // SAME table, formed in O(1).
    assert(c1.table_ == c2.table_ &&
           "[C91 §2.4 tex 133]: C₁ and C₂ must be subchains of the same "
           "input table");
    // [C91 §3 tex 160]: "we assume that C₁ ∩ C₂ is a vertex of P".
    assert(c2.offset_ == c1.offset_ + c1.len_ - 1 &&
           "[C91 §3 tex 160]: C₁ ∩ C₂ must be exactly the shared vertex "
           "(C₁'s last = C₂'s first)");

    table_ = c1.table_;
    offset_ = c1.offset_;
    len_ = c1.len_ + c2.len_ - 1;

    // [C91 §2 tex 47]: the union's y-extremes are the SoS-extreme of
    // the two halves' extremes — O(1), no rescan.
    const std::size_t shift = c1.len_ - 1;
    max_y_vertex_ = point_y_above(c2.vertex(c2.max_y_vertex_),
                                  c1.vertex(c1.max_y_vertex_))
        ? c2.max_y_vertex_ + shift : c1.max_y_vertex_;
    min_y_vertex_ = point_y_below(c2.vertex(c2.min_y_vertex_),
                                  c1.vertex(c1.min_y_vertex_))
        ? c2.min_y_vertex_ + shift : c1.min_y_vertex_;
}

Polygon Polygon::subchain(std::size_t first, std::size_t count) const {
    assert(count >= 2 && "[C91 §2.1]: curve needs ≥ 2 vertices");
    assert(first + count <= len_ &&
           "[C91 §2.4 tex 133]: subchain must lie within the curve");
    Polygon p;
    p.table_ = table_;
    p.offset_ = offset_ + first;
    p.len_ = count;
    p.find_y_extremes();
    return p;
}

void Polygon::find_y_extremes() {
    // [C91 §3.4 tex 306] + [C91 §2 tex 47]: global y-extreme vertices
    // under the perturbed (SoS) order; one O(n) pass at input setup.
    max_y_vertex_ = 0;
    min_y_vertex_ = 0;
    for (std::size_t i = 1; i < len_; ++i) {
        if (point_y_above(vertex(i), vertex(max_y_vertex_)))
            max_y_vertex_ = i;
        if (point_y_below(vertex(i), vertex(min_y_vertex_)))
            min_y_vertex_ = i;
    }
}

} // namespace chazelle
