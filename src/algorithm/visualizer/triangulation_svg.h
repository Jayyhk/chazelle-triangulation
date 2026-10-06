#pragma once

#include "../triangulation/triangulation.h"

#include <iosfwd>
#include <span>

namespace chazelle {

void write_triangulation_svg(std::ostream& output, std::span<const Point> vertices,
                             std::span<const Triangle> triangles);

}
