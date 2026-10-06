#pragma once

#include "exact.h"

#include <cstddef>

namespace chazelle {

enum Side : unsigned char { LEFT, RIGHT };

struct Point {
    Exact x = 0.0;
    Exact y = 0.0;

    std::size_t index = 0;
};

}
