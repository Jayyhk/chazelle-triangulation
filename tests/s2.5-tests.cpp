// tests/s2.5-tests.cpp — Tests for [C91 §2.5]: Lemma 2.4's piece
// structure.
//
// The SIDE assignment (which B each piece is shielded from) is
// configuration-dependent ([C91 §2.5 tex 156]) and computed from real
// geometry in conformality.cpp::descend_step — its behavior is pinned
// by the §3.2/§4.1 suites.  This suite pins the piece COUNTS and the
// lemma's structural preconditions.

#include "submap/shielding.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════
//  1. No intersection — trivial case
// ════════════════════════════════════════════════════════════════

static void test_no_intersection() {
    // [C91 §2.5]: "else the lemma is trivially correct."  A does not
    // cross ab: one piece.
    assert(shielding_piece_count(false, false) == 1);
    std::printf("  [PASS] no_intersection\n");
}

// ════════════════════════════════════════════════════════════════
//  2. Only a' — single crossing (Fig 2.8.2)
// ════════════════════════════════════════════════════════════════

static void test_a_prime_only() {
    // A crosses ab once at a': pieces [c, a'] and [a', d].
    assert(shielding_piece_count(true, false) == 2);
    std::printf("  [PASS] a_prime_only\n");
}

// ════════════════════════════════════════════════════════════════
//  3. Both a' and b', distinct — double crossing (Fig 2.8.1)
// ════════════════════════════════════════════════════════════════

static void test_both_distinct() {
    // Pieces [c, a'], [a', b'], [b', d].
    assert(shielding_piece_count(true, true, false) == 3);
    std::printf("  [PASS] both_distinct\n");
}

// ════════════════════════════════════════════════════════════════
//  4. Both a' and b', coincident — single touch point
// ════════════════════════════════════════════════════════════════

static void test_both_coincident() {
    // A touches ab at one point (a' = b', no crossing): pieces
    // [c, a'] and [a', d].
    assert(shielding_piece_count(true, true, true) == 2);
    std::printf("  [PASS] both_coincident\n");
}

// ════════════════════════════════════════════════════════════════
//  5. Structural properties over all admissible inputs
// ════════════════════════════════════════════════════════════════

static void test_properties() {
    // [C91 §2.5 Lemma 2.4 tex 153]: "subdivide A into a total of one,
    // two, or three connected curves" — over every admissible input,
    // the piece count is 1 + the number of DISTINCT subdivision
    // points of ab ∩ A (a' and b' coinciding contribute one point).
    for (int a = 0; a <= 1; ++a) {
        for (int b = 0; b <= 1; ++b) {
            if (b && !a) continue;  // impossible: b' requires a'
            for (int eq = 0; eq <= 1; ++eq) {
                if (!(a && b) && eq) continue;  // eq needs both
                const std::size_t n =
                    shielding_piece_count(a != 0, b != 0, eq != 0);
                const std::size_t points =
                    !a ? 0u : (!b ? 1u : (eq ? 1u : 2u));
                assert(n == 1 + points);
            }
        }
    }
    std::printf("  [PASS] properties\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::printf("[C91 §2.5 tests]:\n");
    test_no_intersection();
    test_a_prime_only();
    test_both_distinct();
    test_both_coincident();
    test_properties();
    std::printf("All §2.5 tests passed.\n");
    return 0;
}
