#pragma once

/// [C91 §2.5]: Shielding classification (Lemma 2.4).
///
/// Setup: closed disk D with boundary ∂D, diametrical chord ab.
/// Points a, c, b are clockwise on ∂D.  H₁ = cw arc a→b (contains c),
/// H₂ = cw arc b→a.  Curve A runs from c to d inside D.
/// T = cw arc d→c.  B₁ = closure(T ∩ H₁), B₂ = closure(T ∩ H₂).
///
/// Lemma 2.4: intersection points a' (first of ab∩A from a) and b'
/// (last of ab∩A from a) subdivide A into 1–3 pieces, each shielded
/// from some B_j.  The classification depends only on which of a', b'
/// exist and whether they coincide.

#include <cassert>
#include <cstddef>

namespace chazelle {

/// Identifies which boundary arc a piece of A is shielded from.
/// B1 = closure of (cw d→c arc) ∩ H₁.  Contains the side of ∂D near c.
/// B2 = closure of (cw d→c arc) ∩ H₂.  The other side.
enum class BoundarySide : int { B1 = 1, B2 = 2 };

/// Result of Lemma 2.4: A is divided into 1–3 pieces,
/// each shielded from some B_j.
///
/// Pieces are in order along A from c to d:
///   piece[0] = [c, first_cut)
///   piece[1] = [first_cut, second_cut)  (only if num_pieces == 3)
///   piece[num_pieces-1] = [..., d]
struct ShieldingResult {
    static constexpr std::size_t MAX_PIECES = 3;

    /// Which B_j each piece is shielded from.
    BoundarySide shielded_from[MAX_PIECES] = {};

    /// Number of pieces (1, 2, or 3).
    std::size_t num_pieces = 0;
};

/// [C91 §2.5 Lemma 2.4]: Classify shielding for curve A.
///
/// The curve A runs from c to d inside disk D.  Chord ab divides D;
/// a, c, b are in clockwise order on ∂D.
///
/// @param a_prime  true if a' exists (A intersects ab, and
///                 a ∈ B₁ ∪ B₂, i.e., the cw d→c arc passes
///                 through a).
/// @param b_prime  true if b' exists (A intersects ab, and
///                 b ∈ B₁ ∪ B₂, i.e., the cw d→c arc passes
///                 through b).
/// @param a_prime_eq_b_prime  true if a' = b' (single intersection
///                            point); only meaningful when both exist.
///
/// Preconditions (from the geometry — see §2.5 proof):
///   - (a_prime=false, b_prime=true) is impossible: the cw d→c arc
///     cannot pass through b without passing through a, since
///     a, c, b are clockwise and c ∈ H₁.
///   - a_prime_eq_b_prime is only read when both a_prime and b_prime
///     are true.
///
/// Returns the shielding classification for each piece of A.
inline ShieldingResult classify_shielding(bool a_prime,
                                           bool b_prime,
                                           bool a_prime_eq_b_prime
                                               = false) noexcept {
    // [C91 §2.5]: "the third case, where the counterclockwise
    // traversal from c to d avoids both a and b, was eliminated
    // earlier, since it corresponds to a situation where one of
    // the B_i is empty."  This means b cannot be on the cw d→c
    // arc without a also being on it.
    assert(!(b_prime && !a_prime) &&
           "§2.5: b' cannot exist without a' "
           "(cw d→c arc cannot reach b without passing a)");

    ShieldingResult r;

    if (!a_prime && !b_prime) {
        // [C91 §2.5]: "else the lemma is trivially correct."
        // A does not cross ab.  Entirely on B₁-side (c's side).
        // Shielded from B₂.
        r.num_pieces = 1;
        r.shielded_from[0] = BoundarySide::B2;
    }
    else if (a_prime && !b_prime) {
        // Only a' exists.  d is on H₂ (opposite side from c).
        // A crosses ab once.
        // [c, a']: near c (B₁-side) → shielded from B₂.
        // [a', d]: crossed to B₂-side → shielded from B₁.
        r.num_pieces = 2;
        r.shielded_from[0] = BoundarySide::B2;
        r.shielded_from[1] = BoundarySide::B1;
    }
    else if (a_prime && b_prime && a_prime_eq_b_prime) {
        // Both exist, single intersection point (a' = b').
        // d is on H₁ (same side as c).  A touches ab without
        // fully crossing.  Both pieces stay on B₁-side.
        // Shielded from B₂.
        r.num_pieces = 2;
        r.shielded_from[0] = BoundarySide::B2;
        r.shielded_from[1] = BoundarySide::B2;
    }
    else {
        // Both a' and b' exist, a' ≠ b'.  d is on H₁ (same side
        // as c).  A crosses ab twice.
        // [c, a']: B₁-side → shielded from B₂.
        // [a', b']: crossed to B₂-side → shielded from B₁.
        // [b', d]: crossed back to B₁-side → shielded from B₂.
        assert(a_prime && b_prime && !a_prime_eq_b_prime);
        r.num_pieces = 3;
        r.shielded_from[0] = BoundarySide::B2;
        r.shielded_from[1] = BoundarySide::B1;
        r.shielded_from[2] = BoundarySide::B2;
    }

    return r;
}

} // namespace chazelle
