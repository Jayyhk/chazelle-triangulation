/// tests/s2.5-tests.cpp — Tests for §2.5: Shielding (Lemma 2.4).

#include "submap/shielding.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

// ════════════════════════════════════════════════════════════════
//  1. No intersection — trivial case
// ════════════════════════════════════════════════════════════════

static void test_no_intersection() {
    // [C91 §2.5]: "else the lemma is trivially correct."
    // A does not cross ab.  1 piece, shielded from B₂.
    auto r = classify_shielding(false, false);
    assert(r.num_pieces == 1);
    assert(r.shielded_from[0] == BoundarySide::B2);

    std::printf("  [PASS] no_intersection\n");
}

// ════════════════════════════════════════════════════════════════
//  2. Only a' — single crossing (d on opposite side)
// ════════════════════════════════════════════════════════════════

static void test_a_prime_only() {
    // d on H₂.  A crosses ab once at a'.
    // [c, a'] shielded from B₂.
    // [a', d] shielded from B₁.
    auto r = classify_shielding(true, false);
    assert(r.num_pieces == 2);
    assert(r.shielded_from[0] == BoundarySide::B2);
    assert(r.shielded_from[1] == BoundarySide::B1);

    std::printf("  [PASS] a_prime_only\n");
}

// ════════════════════════════════════════════════════════════════
//  3. Both a' and b', distinct — double crossing
// ════════════════════════════════════════════════════════════════

static void test_both_distinct() {
    // d on H₁ (between c and b).  A crosses ab twice.
    // [c, a'] from B₂, [a', b'] from B₁, [b', d] from B₂.
    auto r = classify_shielding(true, true, false);
    assert(r.num_pieces == 3);
    assert(r.shielded_from[0] == BoundarySide::B2);
    assert(r.shielded_from[1] == BoundarySide::B1);
    assert(r.shielded_from[2] == BoundarySide::B2);

    std::printf("  [PASS] both_distinct\n");
}

// ════════════════════════════════════════════════════════════════
//  4. Both a' and b', coincident — single touch point
// ════════════════════════════════════════════════════════════════

static void test_both_coincident() {
    // d on H₁.  A touches ab at one point (a' = b').
    // Both pieces stay on B₁-side.
    // [c, a'] from B₂, [a', d] from B₂.
    auto r = classify_shielding(true, true, true);
    assert(r.num_pieces == 2);
    assert(r.shielded_from[0] == BoundarySide::B2);
    assert(r.shielded_from[1] == BoundarySide::B2);

    std::printf("  [PASS] both_coincident\n");
}

// ════════════════════════════════════════════════════════════════
//  5. Shielding properties — consistency checks
// ════════════════════════════════════════════════════════════════

static void test_properties() {
    // Property: first piece is ALWAYS shielded from B₂ (c is on B₁-side).
    for (int a = 0; a <= 1; ++a) {
        for (int b = 0; b <= 1; ++b) {
            if (b && !a) continue; // impossible
            for (int eq = 0; eq <= 1; ++eq) {
                if (!(a && b) && eq) continue; // eq only meaningful when both
                auto r = classify_shielding(a != 0, b != 0, eq != 0);
                assert(r.num_pieces >= 1 && r.num_pieces <= 3);
                assert(r.shielded_from[0] == BoundarySide::B2 &&
                       "§2.5: first piece (from c) must always be "
                       "shielded from B₂");
            }
        }
    }

    // Property: when 3 pieces, middle piece is shielded from B₁
    // (it crossed to the B₂-side of ab).
    {
        auto r = classify_shielding(true, true, false);
        assert(r.num_pieces == 3);
        assert(r.shielded_from[1] == BoundarySide::B1);
    }

    // Property: when 3 pieces, first and last are both from B₂
    // (symmetric — A crosses ab and crosses back).
    {
        auto r = classify_shielding(true, true, false);
        assert(r.shielded_from[0] == BoundarySide::B2);
        assert(r.shielded_from[2] == BoundarySide::B2);
    }

    std::printf("  [PASS] properties\n");
}

// ════════════════════════════════════════════════════════════════
//  6. Fig 2.8.1 from the paper
// ════════════════════════════════════════════════════════════════

static void test_fig281() {
    // Fig 2.8.1: A crosses ab at both a' and b'.
    // "the piece of A running from c to a' is shielded from B₂"
    // "the piece from a' to b' is shielded from B₁"
    auto r = classify_shielding(true, true, false);
    assert(r.num_pieces == 3);
    assert(r.shielded_from[0] == BoundarySide::B2); // c to a'
    assert(r.shielded_from[1] == BoundarySide::B1); // a' to b'
    assert(r.shielded_from[2] == BoundarySide::B2); // b' to d

    std::printf("  [PASS] fig281\n");
}

// ════════════════════════════════════════════════════════════════
//  7. Fig 2.8.2 — only a' exists
// ════════════════════════════════════════════════════════════════

static void test_fig282() {
    // Fig 2.8.2: d on opposite side, A crosses ab once.
    // [c, a'] shielded from B₂, [a', d] shielded from B₁.
    auto r = classify_shielding(true, false);
    assert(r.num_pieces == 2);
    assert(r.shielded_from[0] == BoundarySide::B2);
    assert(r.shielded_from[1] == BoundarySide::B1);

    std::printf("  [PASS] fig282\n");
}

// ════════════════════════════════════════════════════════════════
//  8. BoundarySide enum values
// ════════════════════════════════════════════════════════════════

static void test_enum() {
    assert(static_cast<int>(BoundarySide::B1) == 1);
    assert(static_cast<int>(BoundarySide::B2) == 2);
    assert(BoundarySide::B1 != BoundarySide::B2);

    std::printf("  [PASS] enum\n");
}

// ════════════════════════════════════════════════════════════════

int main() {
    std::printf("§2.5 tests:\n");
    test_no_intersection();
    test_a_prime_only();
    test_both_distinct();
    test_both_coincident();
    test_properties();
    test_fig281();
    test_fig282();
    test_enum();
    std::printf("All §2.5 tests passed.\n");
    return 0;
}
