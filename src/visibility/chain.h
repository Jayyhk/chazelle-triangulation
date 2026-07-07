#pragma once

// [C91 §4 tex 314–323]: The Visibility Algorithm — chains, grades, and
// the padded input curve.
//
// [C91 §4 tex 316]: "Let P be a simple nonclosed polygonal curve with n
// vertices v₁, ..., vₙ.  By padding the curve with additional vertices,
// if necessary, we can assume that n = 2^p + 1.  Any subcurve of P of
// the form v_a, ..., v_b, where a − 1 is a multiple of 2^λ and
// b − a = 2^λ is called a chain in grade λ."
//
// [C91 §4 tex 317–321]: "Obviously,
//   (i)   a grade-λ chain has 2^λ + 1 vertices,
//   (ii)  there are 2^{p−λ} chains in grade λ, and
//   (iii) there are p + 1 grades: 0, 1, ..., p."
//
// [C91 §4 tex 323]: the up-phase works bottom-up "computing conformal
// submaps of granularity roughly m^β where m is the size of the
// underlying curve"; the down-phase then refines top-down reusing "data
// structures left behind during the top-phase (the submaps for the
// chains and their ray-shooting structures)".  The phases themselves
// are [C91 §4.1]/[C91 §4.2]; this header owns the grid they run on.

#include "../polygon/polygon.h"

#include <cassert>
#include <cstddef>
#include <vector>

namespace chazelle {

// [C91 §4 tex 323]: "β is some small enough positive constant; we set
// β = 1/5, but to make the complexity analysis more explicit we leave
// β as a parameter in most of the calculations."  Kept as an exact
// rational: every use in the algorithm is of the form ⌈β·k⌉ for an
// integer k (e.g. the canonical granularity 2^{⌈β⌈log m⌉⌉} of
// [C91 §4.1 tex 327]), which must be computed in exact integer
// arithmetic, never floating point.
inline constexpr std::size_t BETA_NUM = 1;
inline constexpr std::size_t BETA_DEN = 5;

// ⌈β·k⌉ over the exact rational β = BETA_NUM / BETA_DEN.
inline constexpr std::size_t ceil_beta(std::size_t k) noexcept {
    return (k * BETA_NUM + BETA_DEN - 1) / BETA_DEN;
}

// [C91 §4 tex 316]: "By padding the curve with additional vertices, if
// necessary, we can assume that n = 2^p + 1" — p minimal.  The padding
// vertices duplicate the last vertex's coordinates with fresh
// consecutive SoS indices ([C91 §2 tex 47] keeps them symbolically
// distinct), so every added edge is of null length: null-length edges,
// arcs, and chords are first-class citizens of the paper's framework
// ([C91 §2.1 tex 72], [C91 §2.2 tex 96]), and a null-length edge
// contributes no weight ([C91 §2.2 tex 106]: weight counts "nonnull
// length edges" only), leaving V(P)'s geometry and every granularity
// budget unchanged.  Since 2^{p−1} + 1 < n for the minimal p, fewer
// than n − 1 vertices are added: O(n) time and space.
std::vector<Point> pad_curve(std::vector<Point> vertices);

// [C91 §4 tex 316–321]: the padded input curve P (n = 2^p + 1) and the
// complete grid of its chains, indexed by (grade, position).  All
// chains in all grades are materialized as views into the one input
// table ([C91 §2.4 tex 133]: "The input table is read-only: it is
// never to be modified or even copied"): there are
// Σ_{0 ≤ λ ≤ p} 2^{p−λ} = 2^{p+1} − 1 = O(n) of them, built bottom-up —
// grade-0 chains as single-edge subchain views, and each grade-λ chain
// as the O(1) union ([C91 §3 tex 160]) of its two grade-(λ−1) halves —
// for O(n) total construction time, as required by the overall linear
// budget ([C91 §4.1 tex 357]).
class GradedCurve {
public:
    // Pads ([C91 §4 tex 316]) and builds every chain of every grade.
    // O(n) time and space.
    explicit GradedCurve(std::vector<Point> vertices);

    // The padded input curve P with n = 2^p + 1 vertices ([C91 §4 tex
    // 316]); identical to chain(p, 0), the single chain in grade p.
    const Polygon& curve() const noexcept { return chains_[p_][0]; }

    // n = 2^p + 1 ([C91 §4 tex 316]).
    std::size_t p() const noexcept { return p_; }

    // [C91 §4 tex 320 (iii)]: "there are p + 1 grades: 0, 1, ..., p".
    std::size_t num_grades() const noexcept { return p_ + 1; }

    // [C91 §4 tex 319 (ii)]: "there are 2^{p−λ} chains in grade λ".
    std::size_t num_chains(std::size_t grade) const noexcept {
        assert(grade <= p_ &&
               "[C91 §4 tex 320 (iii)]: grades are 0, 1, ..., p");
        return std::size_t{1} << (p_ - grade);
    }

    // [C91 §4 tex 316]: chain i of grade λ is the subcurve
    // v_a, ..., v_b with a − 1 = i·2^λ and b − a = 2^λ (paper's
    // 1-based vertex numbering); by fact (i) it has 2^λ + 1 vertices.
    // Consecutive chains of a grade share exactly one vertex of P.
    const Polygon& chain(std::size_t grade, std::size_t i) const noexcept {
        assert(grade <= p_ &&
               "[C91 §4 tex 320 (iii)]: grades are 0, 1, ..., p");
        assert(i < num_chains(grade) &&
               "[C91 §4 tex 319 (ii)]: only 2^{p−λ} chains in grade λ");
        return chains_[grade][i];
    }

private:
    std::size_t p_ = 0;
    // chains_[λ][i] — the full grid ([C91 §4 tex 316–321]).
    std::vector<std::vector<Polygon>> chains_;
};

} // namespace chazelle
