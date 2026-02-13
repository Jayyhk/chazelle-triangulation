#include "chazelle/fusion.h"
#include "geometry/perturbation.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <set>
#include <string>
#include <vector>

namespace chazelle {

namespace {

/// A discovered chord from one fusion pass.
/// Per §3.1 Remark 1: "the chords of S reflect the visibility of the
/// chord endpoints of S₁ and S₂ [and] they need not be incident upon
/// any vertex of ∂C."  Endpoints are stored as (edge, y) pairs.
struct DiscoveredChord {
    double y;                  ///< Y-coordinate of the chord.
    std::size_t left_edge;     ///< Edge where the left endpoint lies.
    std::size_t right_edge;    ///< Edge where the right endpoint lies.
    bool is_null_length = false;
};

/// Given a ray hit on an arc of a submap, resolve which polygon edge
/// the hit point lies on.  Returns the edge index.
///
/// Per §3.1 Remark 1, chord endpoints are arbitrary points on ∂C,
/// identified by their edge index.  The y-coordinate of the chord
/// determines the exact position on that edge.
std::size_t resolve_hit_edge(const Submap& submap,
                              const Polygon& polygon,
                              std::size_t arc_idx,
                              double ray_y,
                              double ray_hit_x) {
    const auto& arc = submap.arc(arc_idx);
    if (arc.first_edge == NONE) return NONE;

    std::size_t lo = std::min(arc.first_edge, arc.last_edge);
    std::size_t hi = std::max(arc.first_edge, arc.last_edge);

    std::size_t best_edge = NONE;
    double best_x_dist = std::numeric_limits<double>::infinity();

    for (std::size_t ei = lo; ei <= hi && ei < polygon.num_edges(); ++ei) {
        const auto& edge = polygon.edge(ei);
        const auto& p1 = polygon.vertex(edge.start_idx);
        const auto& p2 = polygon.vertex(edge.end_idx);

        double y_lo = std::min(p1.y, p2.y);
        double y_hi = std::max(p1.y, p2.y);

        if (ray_y < y_lo - 1e-12 || ray_y > y_hi + 1e-12) continue;
        if (std::abs(y_hi - y_lo) < 1e-15) continue;

        double t = (ray_y - p1.y) / (p2.y - p1.y);
        double edge_x = p1.x + t * (p2.x - p1.x);

        double dx = std::abs(edge_x - ray_hit_x);
        if (dx < best_x_dist) {
            best_x_dist = dx;
            best_edge = ei;
        }
    }

    return best_edge;
}

/// Compute the x-coordinate of a point on ∂C at height y, on the
/// polygon edge incident on vertex v.
double origin_x_for_vertex(const Polygon& polygon, std::size_t v) {
    if (v >= polygon.num_vertices()) return 0.0;
    const auto& pt = polygon.vertex(v);
    double x = pt.x;
    double y = pt.y;
    std::size_t ei = (v > 0 && v <= polygon.num_edges()) ? v - 1 : 0;
    if (ei < polygon.num_edges()) {
        const auto& edge = polygon.edge(ei);
        const auto& p1 = polygon.vertex(edge.start_idx);
        const auto& p2 = polygon.vertex(edge.end_idx);
        if (std::abs(p2.y - p1.y) > 1e-15) {
            double t = (y - p1.y) / (p2.y - p1.y);
            x = p1.x + t * (p2.x - p1.x);
        }
    }
    return x;
}

/// Do local shooting from point (origin_x, y) within a specific region
/// of a submap.  Checks all arcs (at most 4 by conformality) and chords
/// bounding the region.  Returns the nearest hit.
///
/// This is the core primitive of §3.1: "local shooting".
RayHit local_shoot_in_region(const Submap& submap,
                              const Polygon& polygon,
                              std::size_t region_idx,
                              double origin_x, double y,
                              bool shoot_right) {
    const auto& nd = submap.node(region_idx);
    RayHit best;
    best.start_region = region_idx;
    double best_dist = std::numeric_limits<double>::infinity();

    for (std::size_t ci : nd.incident_chords) {
        const auto& c = submap.chord(ci);
        if (std::abs(c.y - y) > 1e-12) continue;
        if (c.left_edge != NONE && c.right_edge != NONE &&
            c.left_edge < polygon.num_edges() &&
            c.right_edge < polygon.num_edges()) {
            double lx = polygon.edge_x_at_y(c.left_edge, c.y);
            double rx = polygon.edge_x_at_y(c.right_edge, c.y);
            double chord_x = shoot_right ? std::min(lx, rx)
                                         : std::max(lx, rx);
            double dist = shoot_right ? (chord_x - origin_x)
                                      : (origin_x - chord_x);
            if (dist > -1e-12 && dist < best_dist) {
                best_dist = dist;
                best.type = RayHit::Type::CHORD;
                best.chord_idx = ci;
                best.hit_x = chord_x;
            }
        }
    }

    for (std::size_t ai : nd.arcs) {
        const auto& a = submap.arc(ai);
        if (a.first_edge == NONE) continue;

        // §4.2: Virtual arcs represent tilted exit-chord edges.
        // Compute intersection with the tilted segment instead of
        // polygon edges.
        if (a.is_virtual()) {
            // The tilted chord connects a point on edge first_edge to
            // a point on edge last_edge, at height virtual_y with a
            // tiny tilt.  Model as a segment from
            //   (x_left, virtual_y - ε) to (x_right, virtual_y + ε)
            // where ε is infinitesimal (symbolic perturbation).
            //
            // For numerical purposes, use a small tilt: the segment
            // spans from y_lo = virtual_y - 1e-8 to y_hi = virtual_y + 1e-8.
            double vy = a.virtual_y;
            constexpr double TILT = 1e-8;
            double y_lo = vy - TILT;
            double y_hi = vy + TILT;
            if (y < y_lo - 1e-12 || y > y_hi + 1e-12) continue;
            // Compute x at the two endpoints of the tilted chord.
            double x_left = (a.first_edge < polygon.num_edges())
                ? polygon.edge_x_at_y(a.first_edge, vy) : 0.0;
            double x_right = (a.last_edge < polygon.num_edges())
                ? polygon.edge_x_at_y(a.last_edge, vy) : 0.0;
            // Linearly interpolate x along the tilted segment.
            double t = (std::abs(y_hi - y_lo) > 1e-15)
                ? (y - y_lo) / (y_hi - y_lo) : 0.5;
            double x = x_left + t * (x_right - x_left);
            double dist = shoot_right ? (x - origin_x) : (origin_x - x);
            if (dist > -1e-12 && dist < best_dist) {
                best_dist = dist;
                best.type = RayHit::Type::ARC;
                best.arc_idx = ai;
                best.hit_x = x;
            }
            continue;
        }

        std::size_t alo = std::min(a.first_edge, a.last_edge);
        std::size_t ahi = std::max(a.first_edge, a.last_edge);
        for (std::size_t ei = alo; ei <= ahi && ei < polygon.num_edges(); ++ei) {
            const auto& edge = polygon.edge(ei);
            const auto& p1 = polygon.vertex(edge.start_idx);
            const auto& p2 = polygon.vertex(edge.end_idx);
            double y_lo = std::min(p1.y, p2.y);
            double y_hi = std::max(p1.y, p2.y);
            if (y < y_lo - 1e-12 || y > y_hi + 1e-12) continue;
            if (std::abs(y_hi - y_lo) < 1e-15) continue;
            double t = (y - p1.y) / (p2.y - p1.y);
            double x = p1.x + t * (p2.x - p1.x);
            double dist = shoot_right ? (x - origin_x) : (origin_x - x);
            if (dist > -1e-12 && dist < best_dist) {
                best_dist = dist;
                best.type = RayHit::Type::ARC;
                best.arc_idx = ai;
                best.hit_x = x;
            }
        }
    }

    return best;
}

/// Determine which region of dst contains the junction vertex.
///
/// §3.1 Start-Up: "using the information about the endpoints of C2
/// encoded in the normal-form representation of S2 (namely, pointers
/// to incident arcs), we can find, in constant time, in which region
/// of S2 the point a0 lies."
std::size_t find_junction_region(const Submap& dst,
                                  std::size_t junction_vertex) {
    for (std::size_t ai = 0; ai < dst.num_arcs(); ++ai) {
        const auto& a = dst.arc(ai);
        if (a.first_edge == NONE) continue;
        std::size_t lo = std::min(a.first_edge, a.last_edge);
        std::size_t hi = std::max(a.first_edge, a.last_edge);
        if (junction_vertex == lo || junction_vertex == hi + 1) {
            if (a.region_node != NONE) return a.region_node;
        }
    }
    for (std::size_t ai = 0; ai < dst.num_arcs(); ++ai) {
        const auto& a = dst.arc(ai);
        if (a.first_edge == NONE) continue;
        std::size_t lo = std::min(a.first_edge, a.last_edge);
        std::size_t hi = std::max(a.first_edge, a.last_edge);
        if (junction_vertex >= lo && junction_vertex <= hi + 1) {
            if (a.region_node != NONE) return a.region_node;
        }
    }
    for (std::size_t i = 0; i < dst.num_nodes(); ++i) {
        if (!dst.node(i).deleted) return i;
    }
    return NONE;
}

/// One pass of fusion (§3.1): fuse src (S₁) into dst (S₂).
///
/// Per §3.1: "To fuse S1 into S2 we let a variable p run through
/// ∂C₁ in clockwise order, stopping at a0, ..., a_{m+1}, as well as
/// at some other places to be specified."
///
/// Stops are at polygon vertices (arc boundaries) and chord endpoints.
/// Chord endpoints may lie at interior points of edges (§3.1 Remark 1).
void fuse_pass(const Submap& src,
               const RayShootingOracle& src_oracle,
               const Submap& dst,
               const RayShootingOracle& /*dst_oracle*/,
               const Polygon& polygon,
               std::size_t junction_vertex,
               std::vector<DiscoveredChord>& discovered) {
    std::fprintf(stderr, "  FUSE_PASS: src(%zu chords, %zu arcs) dst(%zu chords, %zu arcs) junction=%zu\n",
                 src.num_chords(), src.num_arcs(), dst.num_chords(), dst.num_arcs(), junction_vertex);
    for (std::size_t ai = 0; ai < src.num_arcs(); ++ai) {
        const auto& a = src.arc(ai);
        std::fprintf(stderr, "    src arc %zu: edges [%zu,%zu] side=[%s,%s] region=%zu ec=%zu\n",
                     ai, a.first_edge, a.last_edge,
                     a.first_side == Side::LEFT ? "L" : "R",
                     a.last_side == Side::LEFT ? "L" : "R",
                     a.region_node, a.edge_count);
    }
    for (std::size_t ai = 0; ai < dst.num_arcs(); ++ai) {
        const auto& a = dst.arc(ai);
        std::fprintf(stderr, "    dst arc %zu: edges [%zu,%zu] side=[%s,%s] region=%zu ec=%zu\n",
                     ai, a.first_edge, a.last_edge,
                     a.first_side == Side::LEFT ? "L" : "R",
                     a.last_side == Side::LEFT ? "L" : "R",
                     a.region_node, a.edge_count);
    }

    // -- Collect stop points --
    // Each stop has an edge index and y-coordinate defining its exact
    // position on ∂C.  Vertex-based stops (arc boundaries) and
    // edge-based stops (chord endpoints) coexist.
    struct StopPoint {
        std::size_t edge_idx;   ///< Edge on which the stop lies.
        double y;               ///< Y-coordinate.
        bool is_vertex;         ///< True if this stop is a polygon vertex.
        std::size_t vertex_idx; ///< Vertex index (if is_vertex), else NONE.
    };
    std::vector<StopPoint> stops;

    auto is_dup = [&](std::size_t edge, double y) -> bool {
        for (auto& s : stops) {
            if (s.edge_idx == edge && std::abs(s.y - y) < 1e-12) return true;
        }
        return false;
    };

    auto add_vertex_stop = [&](std::size_t v) {
        if (v == NONE || v >= polygon.num_vertices()) return;
        std::size_t e = v;
        if (e >= polygon.num_edges()) e = polygon.num_edges() - 1;
        double y = polygon.vertex(v).y;
        if (!is_dup(e, y))
            stops.push_back({e, y, true, v});
    };

    auto add_edge_stop = [&](std::size_t edge, double y) {
        if (edge == NONE || edge >= polygon.num_edges()) return;
        if (!is_dup(edge, y))
            stops.push_back({edge, y, false, NONE});
    };

    // Chord endpoints of src (may be interior to edges).
    for (std::size_t ci = 0; ci < src.num_chords(); ++ci) {
        const auto& c = src.chord(ci);
        if (c.left_edge != NONE) add_edge_stop(c.left_edge, c.y);
        if (c.right_edge != NONE) add_edge_stop(c.right_edge, c.y);
    }

    // Companion vertices: arc boundary vertices.
    for (std::size_t ai = 0; ai < src.num_arcs(); ++ai) {
        const auto& arc = src.arc(ai);
        if (arc.first_edge == NONE) continue;
        std::size_t lo = std::min(arc.first_edge, arc.last_edge);
        std::size_t hi = std::max(arc.first_edge, arc.last_edge);
        add_vertex_stop(lo);
        if (hi + 1 < polygon.num_vertices()) add_vertex_stop(hi + 1);
    }

    // Always include the junction vertex.
    add_vertex_stop(junction_vertex);

    // Sort by position along ∂C: primary key = edge index, then by
    // parametric position along that edge.
    std::sort(stops.begin(), stops.end(),
              [&](const StopPoint& a, const StopPoint& b) {
                  if (a.edge_idx != b.edge_idx) return a.edge_idx < b.edge_idx;
                  if (a.edge_idx < polygon.num_edges()) {
                      const auto& e = polygon.edge(a.edge_idx);
                      const auto& p1 = polygon.vertex(e.start_idx);
                      const auto& p2 = polygon.vertex(e.end_idx);
                      double dy = p2.y - p1.y;
                      if (std::abs(dy) > 1e-15) {
                          double ta = (a.y - p1.y) / dy;
                          double tb = (b.y - p1.y) / dy;
                          return ta < tb;
                      }
                  }
                  return a.y < b.y;
              });

    if (stops.empty()) return;

    std::fprintf(stderr, "    stops(%zu):", stops.size());
    for (auto& s : stops) {
        if (s.is_vertex)
            std::fprintf(stderr, " (e%zu,v%zu)", s.edge_idx, s.vertex_idx);
        else
            std::fprintf(stderr, " (e%zu,y%.3f)", s.edge_idx, s.y);
    }
    std::fprintf(stderr, "\n");

    // -- Start-up phase --
    std::size_t current_region = find_junction_region(dst, junction_vertex);
    std::fprintf(stderr, "    junction=%zu current_region=%zu\n", junction_vertex, current_region);

    // Find index of a0 (junction vertex) in the stops.
    std::size_t a0_idx = 0;
    for (std::size_t i = 0; i < stops.size(); ++i) {
        if (stops[i].is_vertex && stops[i].vertex_idx == junction_vertex) {
            a0_idx = i;
            break;
        }
    }

    // Helper: compute origin x for any stop point.
    auto stop_ox = [&](const StopPoint& s) -> double {
        if (s.is_vertex) return origin_x_for_vertex(polygon, s.vertex_idx);
        return polygon.edge_x_at_y(s.edge_idx, s.y);
    };

    // Helper: get edge index for shooting from a stop point.
    auto stop_shoot_edge = [&](const StopPoint& s) -> std::size_t {
        if (s.is_vertex) {
            std::size_t v = s.vertex_idx;
            return (v > 0) ? v - 1 : 0;
        }
        return s.edge_idx;
    };

    std::size_t p_edge = stops[a0_idx].edge_idx;
    double p_y = stops[a0_idx].y;
    std::size_t start_k = a0_idx;

    {
        double y = p_y;
        double ox = stop_ox(stops[a0_idx]);
        std::size_t shoot_e = stop_shoot_edge(stops[a0_idx]);

        // Shoot from a0 in src to find what it sees w.r.t. C1.
        double dist_src = std::numeric_limits<double>::infinity();
        for (bool dir : {true, false}) {
            RayHit h = src_oracle.shoot(shoot_e, y, Side::LEFT, dir);
            std::fprintf(stderr, "    startup src shoot edge=%zu y=%.3f dir=%s type=%d hit_x=%.3f arc=%zu\n",
                         shoot_e, y, dir?"R":"L", (int)h.type, h.hit_x, h.arc_idx);
            if (h.type == RayHit::Type::ARC) {
                double d = dir ? (h.hit_x - ox) : (ox - h.hit_x);
                if (d > 1e-9 && d < dist_src) dist_src = d;
            }
        }

        // Shoot from a0 within R of dst to find what it sees w.r.t. C2.
        RayHit hit_dst;
        double dist_dst = std::numeric_limits<double>::infinity();
        if (current_region != NONE) {
            for (bool dir : {true, false}) {
                RayHit h = local_shoot_in_region(
                    dst, polygon, current_region, ox, y, dir);
                std::fprintf(stderr, "    startup dst local_shoot region=%zu ox=%.3f y=%.3f dir=%s type=%d hit_x=%.3f arc=%zu\n",
                             current_region, ox, y, dir?"R":"L", (int)h.type, h.hit_x, h.arc_idx);
                if (h.type == RayHit::Type::ARC) {
                    double d = dir ? (h.hit_x - ox) : (ox - h.hit_x);
                    if (d > -1e-12 && d < dist_dst) {
                        dist_dst = d;
                        hit_dst = h;
                    }
                }
            }
        }

        if (hit_dst.type == RayHit::Type::ARC &&
            dist_dst <= dist_src + 1e-12) {
            // Case 1: c0 in ∂C₂.  Record cross-chain chord a0→c0.
            std::size_t re = resolve_hit_edge(
                dst, polygon, hit_dst.arc_idx, y, hit_dst.hit_x);
            std::size_t a0_edge = stops[a0_idx].edge_idx;
            std::fprintf(stderr, "    startup Case1: a0 edge=%zu sees dC2 at re=%zu\n",
                         a0_edge, re);
            if (re != NONE && re != a0_edge) {
                DiscoveredChord dc;
                dc.y = y;
                dc.left_edge = a0_edge;
                dc.right_edge = re;
                discovered.push_back(dc);
            }
        } else {
            // Case 2: c0 in ∂C₁.
            std::size_t c0_edge = NONE;
            double best_src_dist = std::numeric_limits<double>::infinity();
            for (bool dir : {true, false}) {
                RayHit h = src_oracle.shoot(shoot_e, y, Side::LEFT, dir);
                if (h.type == RayHit::Type::ARC) {
                    double d = dir ? (h.hit_x - ox) : (ox - h.hit_x);
                    if (d > 1e-9 && d < best_src_dist) {
                        best_src_dist = d;
                        c0_edge = resolve_hit_edge(
                            src, polygon, h.arc_idx, y, h.hit_x);
                    }
                }
            }

            if (c0_edge != NONE && c0_edge != stops[a0_idx].edge_idx) {
                std::fprintf(stderr, "    startup Case2: a0 edge=%zu sees dC1 at c0 edge=%zu\n",
                             p_edge, c0_edge);

                DiscoveredChord dc;
                dc.y = y;
                dc.left_edge = c0_edge;
                dc.right_edge = stops[a0_idx].edge_idx;
                discovered.push_back(dc);

                p_edge = c0_edge;
                p_y = y;

                for (std::size_t i = 0; i < stops.size(); ++i) {
                    if (stops[i].edge_idx >= c0_edge) {
                        start_k = i;
                        break;
                    }
                }
            }
        }
    }

    // -- Main loop --
    std::size_t total_stops = stops.size();
    for (std::size_t jj = 0; jj < total_stops; ++jj) {
        std::size_t j = (start_k + jj) % total_stops;
        const auto& stop_j = stops[j];
        if (current_region == NONE) break;

        double y = stop_j.y;
        double ox = stop_ox(stop_j);
        std::size_t aj_edge = stop_j.edge_idx;
        std::size_t shoot_e = stop_shoot_edge(stop_j);
        std::fprintf(stderr, "    main_loop jj=%zu edge=%zu y=%.3f ox=%.3f p_edge=%zu R=%zu\n",
                     jj, aj_edge, y, ox, p_edge, current_region);

        // --- Case (i): aj lies in R and sees ∂C₂ ---
        bool aj_in_R = false;
        RayHit best_dst_hit;
        double best_dst_dist = std::numeric_limits<double>::infinity();

        if (current_region < dst.num_nodes() &&
            !dst.node(current_region).deleted) {
            for (bool dir : {true, false}) {
                RayHit h = local_shoot_in_region(
                    dst, polygon, current_region, ox, y, dir);
                std::fprintf(stderr, "      dst local_shoot dir=%s type=%d hit_x=%.3f arc=%zu\n",
                             dir?"R":"L", (int)h.type, h.hit_x, h.arc_idx);
                if (h.type == RayHit::Type::ARC) {
                    double d = dir ? (h.hit_x - ox) : (ox - h.hit_x);
                    if (d > -1e-12 && d < best_dst_dist) {
                        best_dst_dist = d;
                        best_dst_hit = h;
                        aj_in_R = true;
                    }
                }
            }
        }

        if (aj_in_R) {
            std::fprintf(stderr, "      aj_in_R: best_dst_dist=%.6f hit_x=%.3f arc=%zu\n",
                         best_dst_dist, best_dst_hit.hit_x, best_dst_hit.arc_idx);
            double dist_src = std::numeric_limits<double>::infinity();
            for (bool dir : {true, false}) {
                RayHit h = src_oracle.shoot(shoot_e, y, Side::LEFT, dir);
                std::fprintf(stderr, "      src shoot edge=%zu y=%.3f dir=%s type=%d hit_x=%.3f arc=%zu\n",
                             shoot_e, y, dir?"R":"L", (int)h.type, h.hit_x, h.arc_idx);
                if (h.type == RayHit::Type::ARC) {
                    double d = dir ? (h.hit_x - ox) : (ox - h.hit_x);
                    if (d > 1e-9 && d < dist_src) dist_src = d;
                }
            }

            if (best_dst_dist <= dist_src + 1e-12) {
                std::size_t re = resolve_hit_edge(
                    dst, polygon, best_dst_hit.arc_idx, y,
                    best_dst_hit.hit_x);
                std::fprintf(stderr, "      Case(i) CHORD: aj edge=%zu re=%zu\n",
                             aj_edge, re);
                if (re != NONE && re != aj_edge) {
                    DiscoveredChord dc;
                    dc.y = y;
                    dc.left_edge = aj_edge;
                    dc.right_edge = re;
                    discovered.push_back(dc);
                }
                p_edge = aj_edge;
                p_y = y;
                continue;
            }
        }

        // --- Case (ii): exit chord endpoint of R sees segment past p ---
        if (current_region >= dst.num_nodes()) continue;
        const auto& nd_R = dst.node(current_region);
        if (nd_R.deleted) continue;

        bool case_ii_found = false;
        std::size_t best_p_prime_edge = NONE;
        double best_p_prime_y = 0.0;
        std::size_t best_chord_endpt_edge = NONE;
        double best_chord_endpt_y = 0.0;
        std::size_t best_new_region = NONE;

        for (std::size_t ci : nd_R.incident_chords) {
            const auto& c = dst.chord(ci);

            for (std::size_t endpt_e : {c.left_edge, c.right_edge}) {
                if (endpt_e == NONE || endpt_e >= polygon.num_edges())
                    continue;

                double y_e = c.y;

                for (bool dir : {true, false}) {
                    RayHit h = src_oracle.shoot(endpt_e, y_e, Side::LEFT, dir);
                    if (h.type != RayHit::Type::ARC) continue;

                    const auto& hit_arc = src.arc(h.arc_idx);
                    if (hit_arc.first_edge == NONE) continue;
                    std::size_t h_lo = std::min(hit_arc.first_edge,
                                                 hit_arc.last_edge);
                    std::size_t h_hi = std::max(hit_arc.first_edge,
                                                 hit_arc.last_edge);

                    std::size_t seg_lo = std::min(p_edge, aj_edge);
                    std::size_t seg_hi = std::max(p_edge, aj_edge);
                    if (h_lo > seg_hi || h_hi < seg_lo) continue;

                    std::size_t p_prime_e = resolve_hit_edge(
                        src, polygon, h.arc_idx, y_e, h.hit_x);
                    if (p_prime_e == NONE) continue;
                    if (p_prime_e <= p_edge) continue;

                    if (best_p_prime_edge == NONE || p_prime_e > best_p_prime_edge) {
                        best_p_prime_edge = p_prime_e;
                        best_p_prime_y = y_e;
                        best_chord_endpt_edge = endpt_e;
                        best_chord_endpt_y = y_e;
                        best_new_region = (c.region[0] == current_region)
                                          ? c.region[1] : c.region[0];
                        case_ii_found = true;
                    }
                }
            }
        }

        if (case_ii_found && best_p_prime_edge != NONE) {
            DiscoveredChord dc;
            dc.y = best_chord_endpt_y;
            dc.left_edge = best_chord_endpt_edge;
            dc.right_edge = best_p_prime_edge;
            discovered.push_back(dc);

            p_edge = best_p_prime_edge;
            p_y = best_p_prime_y;
            current_region = best_new_region;
            continue;
        }

        // --- Case (iii): no match --- aj sees ∂C₁, continue.
    }
}

/// Rebuild the submap in normal form from the union of all chords and arcs.
///
/// Per §3.1: "Sort the endpoints of these chords along ∂C (by edge name,
/// then y-coordinate)."
Submap rebuild_submap(const Submap& s1, const Submap& s2,
                      const std::vector<DiscoveredChord>& discovered,
                      const Polygon& /*polygon*/) {
    Submap result;

    struct ChordRecord {
        double y;
        std::size_t left_edge;
        std::size_t right_edge;
        bool is_null_length;
        std::size_t sort_key;
    };

    // §3.1: Chords are uniquely identified by (left_edge, right_edge, y).
    // Two chords on the same edge pair but at different y-coordinates
    // are distinct visibility relations and must both be kept.
    struct ChordKey {
        std::size_t lo, hi;
        double y;
        bool operator<(const ChordKey& o) const {
            if (lo != o.lo) return lo < o.lo;
            if (hi != o.hi) return hi < o.hi;
            return y < o.y;
        }
    };
    std::set<ChordKey> seen;
    std::vector<ChordRecord> all_chords;

    auto add_chord_rec = [&](double y, std::size_t le, std::size_t re,
                             bool is_null) {
        if (le == NONE || re == NONE) return;
        if (le == re && !is_null) return;
        ChordKey key{std::min(le, re), std::max(le, re), y};
        if (!seen.insert(key).second) return;
        std::size_t sk = std::min(le, re);
        all_chords.push_back({y, le, re, is_null, sk});
    };

    for (std::size_t ci = 0; ci < s1.num_chords(); ++ci) {
        const auto& c = s1.chord(ci);
        if (c.region[0] == NONE && c.region[1] == NONE) continue;
        add_chord_rec(c.y, c.left_edge, c.right_edge, c.is_null_length);
    }
    for (std::size_t ci = 0; ci < s2.num_chords(); ++ci) {
        const auto& c = s2.chord(ci);
        if (c.region[0] == NONE && c.region[1] == NONE) continue;
        add_chord_rec(c.y, c.left_edge, c.right_edge, c.is_null_length);
    }
    for (const auto& dc : discovered) {
        add_chord_rec(dc.y, dc.left_edge, dc.right_edge, dc.is_null_length);
    }

    std::sort(all_chords.begin(), all_chords.end(),
              [](const ChordRecord& a, const ChordRecord& b) {
                  if (a.sort_key != b.sort_key) return a.sort_key < b.sort_key;
                  return a.y < b.y;
              });

    std::size_t num_regions = all_chords.size() + 1;
    for (std::size_t i = 0; i < num_regions; ++i)
        result.add_node();

    struct ArcRecord {
        ArcStructure arc;
        std::size_t sort_key;
    };
    std::vector<ArcRecord> all_arcs;

    // §3.1: Arcs from S₁ and S₂ may span across newly discovered
    // chord endpoints.  We must split arcs at chord endpoint edges
    // before assigning them to regions.  Collect all chord endpoint
    // edges as split points.
    std::vector<std::size_t> split_edges;
    for (const auto& cr : all_chords) {
        if (cr.left_edge != NONE)
            split_edges.push_back(cr.left_edge);
        if (cr.right_edge != NONE)
            split_edges.push_back(cr.right_edge);
    }
    std::sort(split_edges.begin(), split_edges.end());
    split_edges.erase(std::unique(split_edges.begin(), split_edges.end()),
                      split_edges.end());

    auto split_and_add_arc = [&](ArcStructure a) {
        if (a.first_edge == NONE) {
            all_arcs.push_back({a, 0});
            return;
        }
        std::size_t lo = std::min(a.first_edge, a.last_edge);
        std::size_t hi = std::max(a.first_edge, a.last_edge);

        // Find all split points strictly interior to [lo, hi].
        auto it_lo = std::upper_bound(split_edges.begin(),
                                       split_edges.end(), lo);
        auto it_hi = std::lower_bound(split_edges.begin(),
                                       split_edges.end(), hi);

        // Collect interior split points.
        std::vector<std::size_t> splits;
        for (auto it = it_lo; it != it_hi; ++it)
            splits.push_back(*it);

        if (splits.empty()) {
            // No splits needed — arc stays whole.
            all_arcs.push_back({a, lo});
        } else {
            // Split the arc at each interior split point.
            std::size_t cur_lo = lo;
            for (std::size_t sp : splits) {
                if (sp <= cur_lo) continue;
                // Sub-arc [cur_lo, sp-1].
                ArcStructure sub = a;
                sub.first_edge = cur_lo;
                sub.last_edge  = sp - 1;
                sub.edge_count = sp - cur_lo;
                all_arcs.push_back({sub, cur_lo});
                cur_lo = sp;
            }
            // Final sub-arc [cur_lo, hi].
            ArcStructure sub = a;
            sub.first_edge = cur_lo;
            sub.last_edge  = hi;
            sub.edge_count = hi - cur_lo + 1;
            all_arcs.push_back({sub, cur_lo});
        }
    };

    for (std::size_t i = 0; i < s1.num_arcs(); ++i)
        split_and_add_arc(s1.arc(i));
    for (std::size_t i = 0; i < s2.num_arcs(); ++i)
        split_and_add_arc(s2.arc(i));

    std::sort(all_arcs.begin(), all_arcs.end(),
              [](const ArcRecord& a, const ArcRecord& b) {
                  return a.sort_key < b.sort_key;
              });

    for (auto& arec : all_arcs) {
        auto it = std::upper_bound(
            all_chords.begin(), all_chords.end(), arec.sort_key,
            [](std::size_t key, const ChordRecord& cr) {
                return key < cr.sort_key;
            });
        std::size_t region = static_cast<std::size_t>(
            std::distance(all_chords.begin(), it));
        if (region >= num_regions) region = num_regions - 1;

        arec.arc.region_node = region;
        std::size_t ai = result.add_arc(arec.arc);
        result.node(region).arcs.push_back(ai);
    }

    for (std::size_t ci = 0; ci < all_chords.size(); ++ci) {
        Chord chord;
        chord.y = all_chords[ci].y;
        chord.left_edge = all_chords[ci].left_edge;
        chord.right_edge = all_chords[ci].right_edge;
        chord.is_null_length = all_chords[ci].is_null_length;
        chord.region[0] = ci;
        chord.region[1] = ci + 1;
        result.add_chord(chord);
        std::fprintf(stderr, "    rebuild chord %zu: le=%zu re=%zu null=%d regions=[%zu,%zu]\n",
                     ci, chord.left_edge, chord.right_edge,
                     (int)chord.is_null_length, ci, ci + 1);
    }

    result.recompute_all_weights();

    // Temporary start_arc/end_arc: normalize() will recompute these
    // properly after the caller calls it, but we set reasonable
    // defaults so that any intermediate access doesn't crash.
    // Find the split between LEFT and RIGHT arcs.
    if (!all_arcs.empty()) {
        std::size_t n_arcs = result.num_arcs();
        std::size_t split = n_arcs;
        for (std::size_t i = 0; i < n_arcs; ++i) {
            if (result.arc(i).first_side == Side::RIGHT) {
                split = i;
                break;
            }
        }
        result.start_arc = 0;
        result.end_arc = (split > 0) ? split - 1
                         : (n_arcs > 0 ? n_arcs - 1 : 0);
    }

    return result;
}

} // anonymous namespace

Submap fuse(const Submap& s1, const RayShootingOracle& oracle1,
            const Submap& s2, const RayShootingOracle& oracle2,
            const Polygon& polygon,
            std::size_t junction_vertex) {
    std::vector<DiscoveredChord> discovered;

    fuse_pass(s1, oracle1, s2, oracle2, polygon, junction_vertex, discovered);
    fuse_pass(s2, oracle2, s1, oracle1, polygon, junction_vertex, discovered);

    if (junction_vertex != NONE && junction_vertex < polygon.num_vertices() &&
        polygon.is_y_extremum(junction_vertex)) {
        double jy = polygon.vertex(junction_vertex).y;
        DiscoveredChord dc;
        dc.y = jy;
        std::size_t je = junction_vertex;
        if (je >= polygon.num_edges()) je = polygon.num_edges() - 1;
        dc.left_edge = je;
        dc.right_edge = je;
        dc.is_null_length = true;
        discovered.push_back(dc);
    }

    return rebuild_submap(s1, s2, discovered, polygon);
}

} // namespace chazelle
