#include "submap/shielding.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

static void test_no_intersection() {
    assert(shielding_piece_count(false, false) == 1);
    std::printf("  [PASS] no_intersection\n");
}

static void test_a_prime_only() {
    assert(shielding_piece_count(true, false) == 2);
    std::printf("  [PASS] a_prime_only\n");
}

static void test_both_distinct() {
    assert(shielding_piece_count(true, true, false) == 3);
    std::printf("  [PASS] both_distinct\n");
}

static void test_both_coincident() {
    assert(shielding_piece_count(true, true, true) == 2);
    std::printf("  [PASS] both_coincident\n");
}

static void test_properties() {
    for (int a = 0; a <= 1; ++a) {
        for (int b = 0; b <= 1; ++b) {
            if (b && !a)
                continue;
            for (int eq = 0; eq <= 1; ++eq) {
                if (!(a && b) && eq)
                    continue;
                const std::size_t n = shielding_piece_count(a != 0, b != 0, eq != 0);
                const std::size_t points = !a ? 0u : (!b ? 1u : (eq ? 1u : 2u));
                assert(n == 1 + points);
            }
        }
    }
    std::printf("  [PASS] properties\n");
}

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
