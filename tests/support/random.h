#pragma once

#include <cstdint>

namespace chazelle::test {

class DeterministicRandomGenerator {
public:
    explicit DeterministicRandomGenerator(std::uint64_t seed = 0)
        : state_(seed ^ 0x9e3779b97f4a7c15ULL) {}

    std::uint64_t next() {
        state_ ^= state_ << 13;
        state_ ^= state_ >> 7;
        state_ ^= state_ << 17;
        return state_;
    }

    double uniform(double lo, double hi) {
        return lo + (hi - lo) * static_cast<double>(next() % 1000000ULL) / 1000000.0;
    }

private:
    std::uint64_t state_;
};

}
