#pragma once

// [C91 §2.5 Lemma 2.4]: Shielding classification.
//
// Setup: closed disk D with diametrical chord ab.  Points a, c, b are
// clockwise on ∂D.  H₁ = cw arc a→b (contains c), H₂ = cw arc b→a.
// Curve A runs c→d inside D.  T = cw arc d→c.  B₁ = closure(T ∩ H₁),
// B₂ = closure(T ∩ H₂).
//
// Intersections a' (first of ab∩A from a) and b' (last) subdivide A into
// 1–3 pieces, each shielded from some B_j.  Classification depends only
// on which of a', b' exist and whether they coincide.

#include <cassert>
#include <cstddef>

namespace chazelle {

// B₁ = c-side of ∂D ∩ cw d→c arc.  B₂ = the other side.
enum class BoundarySide : int { B1 = 1, B2 = 2 };

// [C91 §2.5 Lemma 2.4] result.  Pieces in order along A from c to d.
struct ShieldingResult {
    static constexpr std::size_t MAX_PIECES = 3;
    BoundarySide shielded_from[MAX_PIECES] = {};
    std::size_t num_pieces = 0;
};

// [C91 §2.5 Lemma 2.4]: Classify shielding for curve A.
//
// @param a_prime  true if a' exists (cw d→c arc passes through a).
// @param b_prime  true if b' exists (cw d→c arc passes through b).
// @param a_prime_eq_b_prime  true if a' = b' (single touch, no crossing);
//                            only read when both exist.
//
// Precondition: (a_prime=false, b_prime=true) is impossible — cw d→c
// arc cannot reach b without passing a (a, c, b clockwise, c ∈ H₁).
inline ShieldingResult classify_shielding(bool a_prime,
                                           bool b_prime,
                                           bool a_prime_eq_b_prime
                                               = false) noexcept {
    // §2.5: "the third case ... was eliminated earlier, since it
    // corresponds to a situation where one of the B_i is empty."
    assert(!(b_prime && !a_prime) &&
           "§2.5: b' cannot exist without a'");

    ShieldingResult r;

    if (!a_prime && !b_prime) {
        // §2.5: "else the lemma is trivially correct."  A doesn't cross
        // ab, so it stays entirely on the c-side (B₁-side), shielded
        // from B₂.
        r.num_pieces = 1;
        r.shielded_from[0] = BoundarySide::B2;
    }
    else if (a_prime && !b_prime) {
        // A crosses ab once.  d is on H₂ (opposite side from c).
        // [c, a']: still near c (B₁-side) → shielded from B₂.
        // [a', d]: crossed to B₂-side → shielded from B₁.
        r.num_pieces = 2;
        r.shielded_from[0] = BoundarySide::B2;
        r.shielded_from[1] = BoundarySide::B1;
    }
    else if (a_prime && b_prime && a_prime_eq_b_prime) {
        // Both intersections at the same point (a' = b' — A touches ab
        // without fully crossing).  d on H₁ (same side as c).  Both
        // pieces stay on the B₁-side, shielded from B₂.
        r.num_pieces = 2;
        r.shielded_from[0] = BoundarySide::B2;
        r.shielded_from[1] = BoundarySide::B2;
    }
    else {
        // A crosses ab twice (a' ≠ b').  d on H₁ (same side as c).
        // [c, a']:  B₁-side                → shielded from B₂.
        // [a', b']: crossed to B₂-side     → shielded from B₁.
        // [b', d]:  crossed back to B₁-side → shielded from B₂.
        assert(a_prime && b_prime && !a_prime_eq_b_prime);
        r.num_pieces = 3;
        r.shielded_from[0] = BoundarySide::B2;
        r.shielded_from[1] = BoundarySide::B1;
        r.shielded_from[2] = BoundarySide::B2;
    }

    // [C91 §2.5 Lemma 2.4] tex 153: "subdivide A into a total of one,
    // two, or three connected curves."
    assert(r.num_pieces >= 1 && r.num_pieces <= 3 &&
           "§2.5 Lemma 2.4: A must subdivide into 1, 2, or 3 pieces");

    return r;
}

} // namespace chazelle
