#include "triangulate/convert.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

#include "geometry/perturbation.h"

namespace chazelle {

ConvertResult convert_submap_to_fm(const Submap& vp,
                                   const Polygon& polygon) {
    ConvertResult result;

    const std::size_t nv = polygon.original_size();
    if (nv < 3) return result;

    // ── Map padded vertex indices to original vertices ──────────
    //
    // After padding, the polygon has 2^p + 1 vertices:
    //   indices 0..nv-2     → original vertices (unchanged)
    //   indices nv-1..N-2   → collinear padding on last original edge
    //   index   N-1         → the last original vertex (was index nv-1)
    //
    // Map any padded index back to the nearest original vertex.
    // Padded vertices lie on the edge from original vertex nv-2 to
    // nv-1, so they map to one of those two endpoints based on
    // which is geometrically closer (by perturbed y-comparison).
    const std::size_t N = polygon.num_vertices();
    auto map_to_original = [&](std::size_t vi) -> std::size_t {
        if (vi < nv - 1) return vi;        // original vertex, unchanged
        if (vi == N - 1) return nv - 1;    // last original vertex
        // Padded vertex on last original edge → nv-2 or nv-1.
        // Use the edge endpoint mapping: the padded vertices subdivide
        // the edge from vertex nv-2 to nv-1.  Map to the closer
        // original endpoint.  Since these vertices are collinear,
        // use parametric position: indices closer to nv-1 map to nv-2,
        // indices closer to N-1 map to nv-1.
        std::size_t mid = (nv - 1 + N - 1) / 2;
        return (vi <= mid) ? nv - 2 : nv - 1;
    };

    // ── Build vertex-node linked list ───────────────────────────
    // One node per original polygon vertex, in boundary order.
    // After padding, vertices 0..nv-2 keep their indices, but the
    // last original vertex (index nv-1) moves to the end of the
    // padded array.  We must use the original coordinates.
    std::vector<Point> pts;
    pts.reserve(nv);
    for (std::size_t i = 0; i + 1 < nv; ++i) {
        pts.push_back(polygon.vertex(i));
    }
    // The last original vertex is at the end of the padded polygon.
    Point last_pt = polygon.vertex(N - 1);
    last_pt.index = nv - 1; // restore original index
    pts.push_back(last_pt);
    result.nodes = build_vertex_list(pts);
    result.first = 0;
    result.last  = nv - 1;

    // ── Build trapezoids from V(P) submap regions ───────────────
    // Each non-degenerate region produces a trapezoid whose top/bottom
    // vertices are the region's y-extreme boundary vertices.
    for (std::size_t ri = 0; ri < vp.num_nodes(); ++ri) {
        const auto& nd = vp.node(ri);
        if (nd.deleted) continue;

        std::size_t top_v = FM_NONE;
        std::size_t bot_v = FM_NONE;
        std::size_t left_e = FM_NONE;
        std::size_t right_e = FM_NONE;

        // Scan chords incident on this region for top/bottom.
        for (std::size_t ci : nd.incident_chords) {
            const auto& c = vp.chord(ci);
            if (c.left_vertex == NONE && c.right_vertex == NONE) continue;

            for (std::size_t cv_raw : {c.left_vertex, c.right_vertex}) {
                if (cv_raw == NONE) continue;
                std::size_t cv = map_to_original(cv_raw);
                Point cp = pts[cv];
                if (top_v == FM_NONE || perturbed_y_compare(cp, pts[top_v]) > 0) {
                    top_v = cv;
                }
                if (bot_v == FM_NONE || perturbed_y_compare(cp, pts[bot_v]) < 0) {
                    bot_v = cv;
                }
            }
        }

        // Scan arcs for vertex extrema and left/right edges.
        for (std::size_t ai : nd.arcs) {
            const auto& a = vp.arc(ai);
            if (a.first_edge == NONE) continue;

            std::size_t lo = std::min(a.first_edge, a.last_edge);
            std::size_t hi = std::max(a.first_edge, a.last_edge);

            if (left_e == FM_NONE || lo < left_e) left_e = lo;
            if (right_e == FM_NONE || hi > right_e) right_e = hi;

            for (std::size_t ei : {lo, hi}) {
                if (ei >= polygon.num_edges()) continue;
                const auto& edge = polygon.edge(ei);
                for (std::size_t vi_raw : {edge.start_idx, edge.end_idx}) {
                    std::size_t vi = map_to_original(vi_raw);
                    Point vp = pts[vi];
                    if (top_v == FM_NONE || perturbed_y_compare(vp, pts[top_v]) > 0) {
                        top_v = vi;
                    }
                    if (bot_v == FM_NONE || perturbed_y_compare(vp, pts[bot_v]) < 0) {
                        bot_v = vi;
                    }
                }
            }
        }

        if (top_v == FM_NONE || bot_v == FM_NONE) continue;
        if (top_v == bot_v) continue; // degenerate

        Trapezoid t;
        t.top_vertex    = top_v;
        t.bottom_vertex = bot_v;
        t.left_edge     = left_e;
        t.right_edge    = right_e;

        std::size_t ti = result.traps.size();
        result.traps.push_back(t);

        if (top_v < result.nodes.size() &&
            result.nodes[top_v].trapezoid_idx == std::numeric_limits<std::size_t>::max()) {
            result.nodes[top_v].trapezoid_idx = ti;
        }
        if (bot_v < result.nodes.size() &&
            result.nodes[bot_v].trapezoid_idx == std::numeric_limits<std::size_t>::max()) {
            result.nodes[bot_v].trapezoid_idx = ti;
        }
    }

    return result;
}

} // namespace chazelle