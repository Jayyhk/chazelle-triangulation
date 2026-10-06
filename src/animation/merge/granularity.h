#pragma once

#include "../polygon/polygon.h"
#include "../submap/submap.h"

#include <cstddef>

namespace chazelle::animation {

void enforce_granularity(Submap& submap, const Polygon& curve, std::size_t granularity);

}
