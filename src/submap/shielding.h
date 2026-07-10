#pragma once

// [C91 §2.5 Lemma 2.4]: piece structure of the curve A under a chord.
//
// Setup: closed disk D with diametrical chord ab.  Points a, c, b are
// clockwise on ∂D.  Curve A runs c→d inside D.  The intersections
// a' (first point of ab∩A from a) and b' (last) subdivide A into
// "a total of one, two, or three connected curves" ([C91 §2.5
// tex 153]), each shielded from one of B₁/B₂.
//
// This header supplies only the PIECE COUNT.  WHICH B each piece is
// shielded from is configuration-dependent — [C91 §2.5 tex 156]:
// "Which one can be determined by simple examination of the relative
// order of the points a, b, c, d, a', b' around the boundary of R" —
// and is computed at the point of use from the actual geometry
// (conformality.cpp::descend_step anchors the side parity on the
// chord itself; a fixed fig-2.8.1-style table is wrong whenever the
// configuration differs from the figure, e.g. when c coincides
// with a).

#include <cassert>
#include <cstddef>

namespace chazelle {

// [C91 §2.5 Lemma 2.4 tex 153]: the number of pieces a' and b'
// subdivide A into.
//
// @param a_prime  true if a' exists (ab ∩ A is nonempty).
// @param b_prime  true if b' exists.
// @param a_prime_eq_b_prime  true if a' = b' (single touch, no
//                            crossing); only read when both exist.
//
// Preconditions (per the lemma's setup, [C91 §2.5 tex 150–156]):
//   - (a_prime=false, b_prime=true) is impossible — a' is the FIRST
//     and b' the LAST point of ab∩A, so they exist together.
//   - Both B₁ AND B₂ are nonempty.  Callers must enforce this before
//     invoking.
inline std::size_t shielding_piece_count(bool a_prime, bool b_prime,
                                         bool a_prime_eq_b_prime
                                             = false) noexcept {
    assert(!(b_prime && !a_prime) &&
           "[C91 §2.5]: b' cannot exist without a'");
    std::size_t n;
    if (!a_prime && !b_prime)
        n = 1;                          // A avoids ab entirely
    else if (a_prime && !b_prime)
        n = 2;                          // one crossing: [c,a'], [a',d]
    else if (a_prime_eq_b_prime)
        n = 2;                          // touch: [c,a'], [a',d]
    else
        n = 3;                          // [c,a'], [a',b'], [b',d]
    // [C91 §2.5 Lemma 2.4 tex 153]: "subdivide A into a total of one,
    // two, or three connected curves."
    assert(n >= 1 && n <= 3);
    return n;
}

} // namespace chazelle
