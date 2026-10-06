#pragma once

#include "../algorithm/triangulation/triangulation.h"
#include "polygon/point.h"

#include <iosfwd>
#include <vector>

namespace chazelle::animation {

chazelle::Triangulation triangulate_with_trace(const std::vector<Point>& vertices,
                                               std::ostream& output);

}
