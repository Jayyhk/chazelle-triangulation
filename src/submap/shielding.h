#pragma once

#include <cassert>
#include <cstddef>

namespace chazelle {

inline std::size_t shielding_piece_count(bool a_prime, bool b_prime,
                                         bool a_prime_eq_b_prime = false) noexcept {
    assert(!(b_prime && !a_prime) && "[C91 §2.5]: b' cannot exist without a'");
    std::size_t n;
    if (!a_prime && !b_prime)
        n = 1;
    else if (a_prime && !b_prime)
        n = 2;
    else if (a_prime_eq_b_prime)
        n = 2;
    else
        n = 3;

    assert(n >= 1 && n <= 3);
    return n;
}

}
