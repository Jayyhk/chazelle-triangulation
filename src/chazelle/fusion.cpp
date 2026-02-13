#include "chazelle/fusion.h"
#include "geometry/perturbation.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <set>
#include <vector>

namespace chazelle {

namespace {

/// A discovered chord from one fusion pass.
struct DiscoveredChord {
    double y;                  ///< Y-coordinate of the chord.
    std::size_t left_vertex;   ///< Left endpoint (polygon vertex idx).
    std::size_t right_vertex;  ///< Right endpoint (polygon vertex idx).
    bool is_null_length = false;
};

/// Given a ray hit on an arc of a submap, resolve the polygon vertex
/// at the hit point.
///
/// Per the paper, chords connect arbitrary visible points on dC — not
/// just polygon vertices.  However, our representation stores chord
/// endpoints as vertex indices.  We resolve by finding the edge that
/// the ray actually crosses and returning the vertex endpoint of that
/// edge that is nearest (in x) to the hit point.
std::size_t resolve_hit_vertex(const Submap& submap,
                                const Polygon& polygon,
                                std::size_t arc_idx,
                                double ray_y,
                                double ray_hit_x) {
    const auto& arc = submap.arc(arc_idx);
    if (arc.first_edge == NONE) return NONE;

    std::size_t lo = std::min(arc.first_edge, arc.last_edge);
    std::size_t hi = std::max(arc.first_edge, arc.last_edge);

    std::size_t best_vertex = NONE;
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
            // Pick the vertex nearest to the hit point in both x and y.
            double d1 = std::abs(p1.x - ray_hit_x) + std::abs(p1.y - ray_y);
            double d2 = std::abs(p2.x - ray_hit_x) + std::abs(p2.y - ray_y);
            best_vertex = (d1 <= d2) ? edge.start_idx : edge.end_idx;
        }
    }

    if (best_vertex == NONE) {
        double best_dist = std::numeric_limits<double>::infinity();
        for (std::size_t v = lo; v <= hi + 1 && v < polygon.num_vertices(); ++v) {
            double dx = std::abs(polygon.vertex(v).x - ray_hit_x) +
                        std::abs(polygon.vertex(v).y - ray_y);
            if (dx < best_dist) {
                best_dist = dx;
                best_vertex = v;
            }
        }
    }

    return best_vertex;
}

/// Compute the x-coordinate of a point on dC at height y, on the
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
/// This is the core primitive of S3.1: "local shooting".
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
        if (c.left_vertex != NONE && c.right_vertex != NONE &&
            c.left_vertex < polygon.num_vertices() &&
            c.right_vertex < polygon.num_vertices()) {
            double lx = polygon.vertex(c.left_vertex).x;
            double rx = polygon.vertex(c.right_vertex).x;
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
/// S3.1 Start-Up: "using the information about the endpoints of C2
/// encoded in the normal-form representation of S2 (namely, pointers
/// to incident arcs), we can find, in constant time, in which region
/// of S2 the point a0 lies."
std::size_t find_junction_region(const Submap& dst,
                                  std::size_t junction_vertex) {
    // The junction vertex is at one endpoint of dst's chain.  Find the
    // arc whose edge range touches junction_vertex and return its region.
    for (std::size_t ai = 0; ai < dst.num_arcs(); ++ai) {
        const auto& a = dst.arc(ai);
        if (a.first_edge == NONE) continue;
        std::size_t lo = std::min(a.first_edge, a.last_edge);
        std::size_t hi = std::max(a.first_edge, a.last_edge);
        if (junction_vertex == lo || junction_vertex == hi + 1) {
            if (a.region_node != NONE) return a.region_node;
        }
    }
    // Fallback: check if junction_vertex is within an arc's range.
    for (std::size_t ai = 0; ai < dst.num_arcs(); ++ai) {
        const auto& a = dst.arc(ai);
        if (a.first_edge == NONE) continue;
        std::size_t lo = std::min(a.first_edge, a.last_edge);
        std::size_t hi = std::max(a.first_edge, a.last_edge);
        if (junction_vertex >= lo && junction_vertex <= hi + 1) {
            if (a.region_node != NONE) return a.region_node;
        }
    }
    // Final fallback: first non-deleted region.
    for (std::size_t i = 0; i < dst.num_nodes(); ++i) {
        if (!dst.node(i).deleted) return i;
    }
    return NONE;
}

/// One pass of fusion (S3.1): fuse src (S1) into dst (S2).
///
/// Implements the paper's algorithm:
///   - Start-up phase to initialize p and current region R of dst
///   - Main loop with cases (i), (ii), (iii)
///   - Current-region tracking through dst
///
/// Per S3.1: "To fuse S1 into S2 we let a variable p run through
/// dC1 in clockwise order, stopping at a0, ..., a_{m+1}, as well as
/// at some other places to be specified.  We determine what p sees
/// along the way, while keeping track of the current region of S2
/// in which p lies."
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
    // -- Collect stop vertices --
    // a1..am = chord endpoints of src
    // a0, a_{m+1} = companion vertices at junction
    struct StopVertex {
        std::size_t vertex_idx;
    };
    std::vector<StopVertex> stops;
    std::set<std::size_t> stop_set;

    auto add_stop = [&](std::size_t v) {
        if (v == NONE || v >= polygon.num_vertices()) return;
        if (stop_set.insert(v).second)
            stops.push_back({v});
    };

    // Chord endpoints.
    for (std::size_t ci = 0; ci < src.num_chords(); ++ci) {
        const auto& c = src.chord(ci);
        add_stop(c.left_vertex);
        add_stop(c.right_vertex);
    }

    // Companion vertices: arc boundary vertices.
    for (std::size_t ai = 0; ai < src.num_arcs(); ++ai) {
        const auto& arc = src.arc(ai);
        if (arc.first_edge == NONE) continue;
        std::size_t lo = std::min(arc.first_edge, arc.last_edge);
        std::size_t hi = std::max(arc.first_edge, arc.last_edge);
        add_stop(lo);
        if (hi + 1 < polygon.num_vertices()) add_stop(hi + 1);
    }

    // Always include the junction vertex.
    add_stop(junction_vertex);

    // Sort by vertex index (clockwise order = increasing vertex index).
    std::sort(stops.begin(), stops.end(),
              [](const StopVertex& a, const StopVertex& b) {
                  return a.vertex_idx < b.vertex_idx;
              });

    if (stops.empty()) return;

    std::fprintf(stderr, "    stops(%zu):", stops.size());
    for (auto& s : stops) std::fprintf(stderr, " %zu", s.vertex_idx);
    std::fprintf(stderr, "\n");

    // -- Start-up phase --
    std::size_t current_region = find_junction_region(dst, junction_vertex);
    std::fprintf(stderr, "    junction=%zu current_region=%zu\n", junction_vertex, current_region);

    // Find index of a0 (junction vertex) in the stops.
    std::size_t a0_idx = 0;
    for (std::size_t i = 0; i < stops.size(); ++i) {
        if (stops[i].vertex_idx == junction_vertex) {
            a0_idx = i;
            break;
        }
    }

    // S3.1 Start-Up: shoot from a0 to find c0.
    std::size_t p_vertex = stops[a0_idx].vertex_idx;
    std::size_t start_k = a0_idx;

    {
        double y = polygon.vertex(p_vertex).y;
        double ox = origin_x_for_vertex(polygon, p_vertex);

        // Shoot from a0 in src to find what it sees w.r.t. C1.
        double dist_src = std::numeric_limits<double>::infinity();
        for (bool dir : {true, false}) {
            std::size_t ei = (p_vertex > 0) ? p_vertex - 1 : 0;
            RayHit h = src_oracle.shoot(ei, y, Side::LEFT, dir);
            std::fprintf(stderr, "    startup src shoot edge=%zu y=%.3f dir=%s type=%d hit_x=%.3f arc=%zu\n",
                         ei, y, dir?"R":"L", (int)h.type, h.hit_x, h.arc_idx);
            if (h.type == RayHit::Type::ARC) {
                double d = dir ? (h.hit_x - ox) : (ox - h.hit_x);
                // Require d > 1e-9 to exclude self-hits (ray hitting
                // the arc at the origin vertex itself, distance ≈ 0).
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
            // Case 1: c0 in dC2.  Record cross-chain chord a0->c0.
            std::size_t rv = resolve_hit_vertex(
                dst, polygon, hit_dst.arc_idx, y, hit_dst.hit_x);
            std::fprintf(stderr, "    startup Case1: a0=%zu sees dC2 at rv=%zu dist_dst=%.3f dist_src=%.3f\n",
                         p_vertex, rv, dist_dst, dist_src);
            if (rv != NONE && rv != p_vertex) {
                DiscoveredChord dc;
                dc.y = y;
                dc.left_vertex = p_vertex;
                dc.right_vertex = rv;
                discovered.push_back(dc);
            }
            // p stays at a0, R stays.
        }
        // Case 2: c0 in dC1 -- skip forward handled by main loop.
    }

    // -- Main loop --
    // Per §3.1: p runs from a₀ to a_{m+1} in clockwise order.  a₀
    // and a_{m+1} are companion vertices at the junction, so the
    // traversal wraps around the full boundary of ∂C₁.  In our
    // representation, stops are sorted by vertex index (= clockwise
    // order).  We iterate from start_k through the end, then wrap
    // around to process stops before start_k.
    std::size_t total_stops = stops.size();
    for (std::size_t jj = 0; jj < total_stops; ++jj) {
        std::size_t j = (start_k + jj) % total_stops;
        std::size_t aj = stops[j].vertex_idx;
        if (current_region == NONE) break;

        // --- Case (i): aj lies in R and sees dC2 ---
        // Per the paper: shoot from aj within R's arcs.  If we get
        // a hit on an arc of R, aj is "in R" and sees dC2.
        double y = polygon.vertex(aj).y;
        double ox = origin_x_for_vertex(polygon, aj);
        std::fprintf(stderr, "    main_loop jj=%zu aj=%zu y=%.3f ox=%.3f p_vertex=%zu R=%zu\n",
                     jj, aj, y, ox, p_vertex, current_region);

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
            // aj sees dC2 at best_dst_hit.
            // Now check if dC2 is nearer than dC1 (what aj sees w.r.t. C1).
            double dist_src = std::numeric_limits<double>::infinity();
            std::size_t ei = (aj > 0) ? aj - 1 : 0;
            for (bool dir : {true, false}) {
                RayHit h = src_oracle.shoot(ei, y, Side::LEFT, dir);
                std::fprintf(stderr, "      src shoot edge=%zu y=%.3f dir=%s type=%d hit_x=%.3f arc=%zu\n",
                             ei, y, dir?"R":"L", (int)h.type, h.hit_x, h.arc_idx);
                if (h.type == RayHit::Type::ARC) {
                    double d = dir ? (h.hit_x - ox) : (ox - h.hit_x);
                    // Require d > 1e-9 to exclude self-hits.
                    if (d > 1e-9 && d < dist_src) dist_src = d;
                }
            }

            if (best_dst_dist <= dist_src + 1e-12) {
                // aj sees dC2 w.r.t. C.  Record cross-chain chord.
                std::size_t rv = resolve_hit_vertex(
                    dst, polygon, best_dst_hit.arc_idx, y,
                    best_dst_hit.hit_x);
                std::fprintf(stderr, "      Case(i) CHORD: aj=%zu rv=%zu dist_dst=%.6f dist_src=%.6f\n",
                             aj, rv, best_dst_dist, dist_src);
                if (rv != NONE && rv != aj) {
                    DiscoveredChord dc;
                    dc.y = y;
                    dc.left_vertex = aj;
                    dc.right_vertex = rv;
                    discovered.push_back(dc);
                }
                p_vertex = aj;
                // R stays the same per case (i).
                continue;
            }
            // dC1 is nearer: aj sees dC1 -- fall through to case (ii).
        }

        // --- Case (ii): exit chord endpoint of R sees Aj past p ---
        // For each exit chord of R, check if either endpoint sees
        // the arc segment between p and aj on dC1.
        if (current_region >= dst.num_nodes()) continue;
        const auto& nd_R = dst.node(current_region);
        if (nd_R.deleted) continue;

        bool case_ii_found = false;
        std::size_t best_p_prime = NONE;
        std::size_t best_chord_endpt = NONE;
        std::size_t best_new_region = NONE;

        for (std::size_t ci : nd_R.incident_chords) {
            const auto& c = dst.chord(ci);

            for (std::size_t endpt_v : {c.left_vertex, c.right_vertex}) {
                if (endpt_v == NONE || endpt_v >= polygon.num_vertices())
                    continue;

                double y_e = polygon.vertex(endpt_v).y;
                std::size_t ei_e = (endpt_v > 0) ? endpt_v - 1 : 0;

                // Shoot from this dst chord endpoint toward dC1.
                for (bool dir : {true, false}) {
                    RayHit h = src_oracle.shoot(ei_e, y_e, Side::LEFT, dir);
                    if (h.type != RayHit::Type::ARC) continue;

                    const auto& hit_arc = src.arc(h.arc_idx);
                    if (hit_arc.first_edge == NONE) continue;
                    std::size_t h_lo = std::min(hit_arc.first_edge,
                                                 hit_arc.last_edge);
                    std::size_t h_hi = std::max(hit_arc.first_edge,
                                                 hit_arc.last_edge);

                    // Check: hit on arc segment [p_vertex, aj]?
                    std::size_t seg_lo = std::min(p_vertex, aj);
                    std::size_t seg_hi = std::max(p_vertex, aj);
                    if (h_lo > seg_hi || h_hi < seg_lo) continue;

                    // Resolve the hit point on dC1.
                    std::size_t p_prime = resolve_hit_vertex(
                        src, polygon, h.arc_idx, y_e, h.hit_x);
                    if (p_prime == NONE) continue;
                    if (p_prime <= p_vertex) continue; // must follow p

                    // Per the paper: pick the "last" p' (largest index).
                    if (best_p_prime == NONE || p_prime > best_p_prime) {
                        best_p_prime = p_prime;
                        best_chord_endpt = endpt_v;
                        best_new_region = (c.region[0] == current_region)
                                          ? c.region[1] : c.region[0];
                        case_ii_found = true;
                    }
                }
            }
        }

        if (case_ii_found && best_p_prime != NONE) {
            // Record cross-chain chord: chord_endpoint <-> p_prime
            double y_c = polygon.vertex(best_chord_endpt).y;
            DiscoveredChord dc;
            dc.y = y_c;
            dc.left_vertex = best_chord_endpt;
            dc.right_vertex = best_p_prime;
            discovered.push_back(dc);

            p_vertex = best_p_prime;
            current_region = best_new_region;
            continue;
        }

        // --- Case (iii): no match --- aj sees dC1, continue.
    }
}

/// Rebuild the submap in normal form from the union of all chords and arcs.
///
/// Per S3.1: "Sort the endpoints of these chords along dC (by edge name,
/// then y-coordinate)."
Submap rebuild_submap(const Submap& s1, const Submap& s2,
                      const std::vector<DiscoveredChord>& discovered,
                      const Polygon& /*polygon*/) {
    Submap result;

    // -- Collect ALL unique chords --
    struct ChordRecord {
        double y;
        std::size_t left_vertex;
        std::size_t right_vertex;
        bool is_null_length;
        std::size_t sort_key;
    };

    std::set<std::pair<std::size_t, std::size_t>> seen;
    std::vector<ChordRecord> all_chords;

    auto add_chord_rec = [&](double y, std::size_t lv, std::size_t rv,
                             bool is_null) {
        if (lv == NONE || rv == NONE) return;
        if (lv == rv && !is_null) return;
        auto key = std::make_pair(std::min(lv, rv), std::max(lv, rv));
        if (!seen.insert(key).second) return;
        // Sort key: minimum vertex index (consistent with arc sort keys).
        std::size_t sk = std::min(lv, rv);
        all_chords.push_back({y, lv, rv, is_null, sk});
    };

    for (std::size_t ci = 0; ci < s1.num_chords(); ++ci) {
        const auto& c = s1.chord(ci);
        if (c.region[0] == NONE && c.region[1] == NONE) continue;
        add_chord_rec(c.y, c.left_vertex, c.right_vertex, c.is_null_length);
    }
    for (std::size_t ci = 0; ci < s2.num_chords(); ++ci) {
        const auto& c = s2.chord(ci);
        if (c.region[0] == NONE && c.region[1] == NONE) continue;
        add_chord_rec(c.y, c.left_vertex, c.right_vertex, c.is_null_length);
    }
    for (const auto& dc : discovered) {
        add_chord_rec(dc.y, dc.left_vertex, dc.right_vertex, dc.is_null_length);
    }

    // Sort along dC: by edge name (sort_key), then y.
    std::sort(all_chords.begin(), all_chords.end(),
              [](const ChordRecord& a, const ChordRecord& b) {
                  if (a.sort_key != b.sort_key) return a.sort_key < b.sort_key;
                  return a.y < b.y;
              });

    // -- Build the submap tree --
    std::size_t num_regions = all_chords.size() + 1;
    for (std::size_t i = 0; i < num_regions; ++i)
        result.add_node();

    // -- Build arcs --
    struct ArcRecord {
        ArcStructure arc;
        std::size_t sort_key;
    };
    std::vector<ArcRecord> all_arcs;

    for (std::size_t i = 0; i < s1.num_arcs(); ++i) {
        auto a = s1.arc(i);
        std::size_t sk = (a.first_edge != NONE)
                         ? std::min(a.first_edge, a.last_edge) : 0;
        all_arcs.push_back({a, sk});
    }
    for (std::size_t i = 0; i < s2.num_arcs(); ++i) {
        auto a = s2.arc(i);
        std::size_t sk = (a.first_edge != NONE)
                         ? std::min(a.first_edge, a.last_edge) : 0;
        all_arcs.push_back({a, sk});
    }

    std::sort(all_arcs.begin(), all_arcs.end(),
              [](const ArcRecord& a, const ArcRecord& b) {
                  return a.sort_key < b.sort_key;
              });

    // Assign each arc to a region.
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

    // Add chords as tree edges.
    for (std::size_t ci = 0; ci < all_chords.size(); ++ci) {
        Chord chord;
        chord.y = all_chords[ci].y;
        chord.left_vertex = all_chords[ci].left_vertex;
        chord.right_vertex = all_chords[ci].right_vertex;
        chord.is_null_length = all_chords[ci].is_null_length;
        chord.region[0] = ci;
        chord.region[1] = ci + 1;
        result.add_chord(chord);
        std::fprintf(stderr, "    rebuild chord %zu: lv=%zu rv=%zu null=%d regions=[%zu,%zu]\n",
                     ci, chord.left_vertex, chord.right_vertex,
                     (int)chord.is_null_length, ci, ci + 1);
    }

    result.recompute_all_weights();

    if (!all_arcs.empty()) {
        result.start_arc = 0;
        result.end_arc = result.num_arcs() > 0 ? result.num_arcs() - 1 : 0;
    }

    return result;
}

} // anonymous namespace

Submap fuse(const Submap& s1, const RayShootingOracle& oracle1,
            const Submap& s2, const RayShootingOracle& oracle2,
            const Polygon& polygon,
            std::size_t junction_vertex) {
    std::vector<DiscoveredChord> discovered;

    // Pass 1: fuse S1 into S2 -- walk S1's canonical vertices,
    // tracking the current region of S2.
    fuse_pass(s1, oracle1, s2, oracle2, polygon, junction_vertex, discovered);

    // Pass 2: fuse S2 into S1 -- walk S2's canonical vertices,
    // tracking the current region of S1.
    fuse_pass(s2, oracle2, s1, oracle1, polygon, junction_vertex, discovered);

    // S3.1 / Lemma 3.1, category (4): If the junction vertex is a
    // local y-extremum, it creates null-length chords.
    if (junction_vertex != NONE && junction_vertex < polygon.num_vertices() &&
        polygon.is_y_extremum(junction_vertex)) {
        double jy = polygon.vertex(junction_vertex).y;
        DiscoveredChord dc;
        dc.y = jy;
        dc.left_vertex = junction_vertex;
        dc.right_vertex = junction_vertex;
        dc.is_null_length = true;
        discovered.push_back(dc);
    }

    return rebuild_submap(s1, s2, discovered, polygon);
}

} // namespace chazelle
