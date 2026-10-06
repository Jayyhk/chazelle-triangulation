#include "triangulate.h"
#include "trace.h"
#include "triangulation/triangulation.h"

#include <utility>

namespace chazelle::animation {

chazelle::Triangulation triangulate_with_trace(const std::vector<Point>& vertices,
                                               std::ostream& output) {
    AnimationTrace trace(output, vertices);
    auto result = triangulate_polygon(vertices);
    trace.finish();
    chazelle::Triangulation output_result;
    output_result.triangles.reserve(result.triangles.size());
    for (const auto& triangle : result.triangles)
        output_result.triangles.push_back({triangle.vertices});
    output_result.vertex_triangles = std::move(result.vertex_triangles);
    output_result.work = {result.work.polygon_vertices, result.work.convexity_tests,
                          result.work.forward_steps,    result.work.backward_steps,
                          result.work.removed_vertices, result.work.triangle_incidents};
    return output_result;
}

}
