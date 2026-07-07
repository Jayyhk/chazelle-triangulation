// src/visibility/chain.cpp

#include "chain.h"

namespace chazelle {

namespace {

// Minimal p with 2^p + 1 ≥ n ([C91 §4 tex 316]: pad "if necessary" —
// no more than necessary).  n ≥ 2 ⟹ well defined, and p < 8·sizeof
// (std::size_t) − 1 for any representable n, so the shift is safe.
std::size_t minimal_p(std::size_t n) {
    assert(n >= 2 && "[C91 §2.1]: curve needs ≥ 2 vertices");
    std::size_t p = 0;
    while ((std::size_t{1} << p) + 1 < n) ++p;
    return p;
}

} // namespace

std::vector<Point> pad_curve(std::vector<Point> vertices) {
    // [C91 §2.4 tex 133]: the padded curve becomes THE input table; its
    // SoS indices must be consecutive (polygon.cpp asserts this), so
    // padding continues the caller's index sequence.
    const std::size_t p = minimal_p(vertices.size());
    const std::size_t target = (std::size_t{1} << p) + 1;

    // Minimality of p: with one fewer bit the curve would not fit.
    assert((p == 0 || (std::size_t{1} << (p - 1)) + 1 < vertices.size()) &&
           "[C91 §4 tex 316]: pad only 'if necessary' — p is minimal");
    // Fewer than n − 1 vertices are added: target − 1 = 2^p < 2(n − 1),
    // so padding is O(n).
    assert(target < 2 * vertices.size() &&
           "[C91 §4 tex 316]: minimal padding less than doubles the curve");

    const Point last = vertices.back();
    while (vertices.size() < target) {
        // Null-length edge: same coordinates, next SoS index
        // ([C91 §2.1 tex 72], [C91 §2.2 tex 106]: weight-free).
        vertices.push_back(Point{last.x, last.y,
                                 vertices.back().index + 1});
    }
    return vertices;
}

GradedCurve::GradedCurve(std::vector<Point> vertices) {
    Polygon P(pad_curve(std::move(vertices)));

    const std::size_t n = P.num_vertices();
    p_ = minimal_p(n);
    assert((std::size_t{1} << p_) + 1 == n &&
           "[C91 §4 tex 316]: padded curve has n = 2^p + 1 vertices");

    // [C91 §4 tex 320 (iii)]: p + 1 grades.
    chains_.resize(p_ + 1);

    // Grade 0: 2^p single-edge chains ([C91 §4 tex 318 (i)]: 2^0 + 1 =
    // 2 vertices each), O(1) subchain views.
    chains_[0].reserve(std::size_t{1} << p_);
    for (std::size_t i = 0; i < (std::size_t{1} << p_); ++i)
        chains_[0].push_back(P.subchain(i, 2));

    // Grade λ from grade λ − 1: chain i of grade λ is v_a..v_b with
    // a − 1 = i·2^λ, b − a = 2^λ ([C91 §4 tex 316]), i.e. exactly the
    // union of chains 2i and 2i + 1 of grade λ − 1, which share the
    // vertex v_{a + 2^{λ−1}} of P.  Each union is O(1) ([C91 §2.4 tex
    // 133]: the input table is never copied), so the whole grid costs
    // Σ_{0 ≤ λ ≤ p} 2^{p−λ} = 2^{p+1} − 1 = O(n).
    for (std::size_t grade = 1; grade <= p_; ++grade) {
        const std::size_t count = std::size_t{1} << (p_ - grade);
        chains_[grade].reserve(count);
        for (std::size_t i = 0; i < count; ++i)
            chains_[grade].push_back(Polygon(chains_[grade - 1][2 * i],
                                             chains_[grade - 1][2 * i + 1]));
    }

    // [C91 §4 tex 317–321]: the "obviously" facts, verified.
    for (std::size_t grade = 0; grade <= p_; ++grade) {
        assert(chains_[grade].size() == (std::size_t{1} << (p_ - grade)) &&
               "[C91 §4 tex 319 (ii)]: 2^{p−λ} chains in grade λ");
        for (const Polygon& c : chains_[grade])
            assert(c.num_vertices() == (std::size_t{1} << grade) + 1 &&
                   "[C91 §4 tex 318 (i)]: a grade-λ chain has 2^λ + 1 "
                   "vertices");
    }

    // The single grade-p chain is all of P: same span, and its O(1)
    // union-combined y-extremes ([C91 §2 tex 47]: SoS makes the order
    // total, so the extreme vertex is unique) agree with P's scan.
    assert(chains_[p_].size() == 1 &&
           "[C91 §4 tex 319 (ii)]: grade p has a single chain");
    assert(chains_[p_][0].num_vertices() == n &&
           &chains_[p_][0].vertex(0) == &P.vertex(0) &&
           "[C91 §4 tex 316]: the grade-p chain is the whole curve P");
    assert(chains_[p_][0].max_y_vertex() == P.max_y_vertex() &&
           chains_[p_][0].min_y_vertex() == P.min_y_vertex() &&
           "[C91 §2 tex 47]: combined y-extremes match the direct scan");
}

} // namespace chazelle
