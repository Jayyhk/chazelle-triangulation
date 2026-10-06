#pragma once

#include "algorithm/polygon/point.h"

namespace chazelle::test {

inline Exact triangle_height(const Point& point) {
    return point.y - Exact(point.index) * Exact::infinitesimal(3);
}

inline Exact triangle_orientation(const Point& a, const Point& b, const Point& c) {
    return (b.x - a.x) * (triangle_height(c) - triangle_height(a)) -
           (triangle_height(b) - triangle_height(a)) * (c.x - a.x);
}

inline bool proper_triangle_crossing(const Point& a, const Point& b, const Point& c,
                                     const Point& d) {
    const Exact abc = triangle_orientation(a, b, c);
    const Exact abd = triangle_orientation(a, b, d);
    const Exact cda = triangle_orientation(c, d, a);
    const Exact cdb = triangle_orientation(c, d, b);
    return ((abc < 0 && abd > 0) || (abc > 0 && abd < 0)) &&
           ((cda < 0 && cdb > 0) || (cda > 0 && cdb < 0));
}

}
