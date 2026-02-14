#include "chazelle/fusion.h"
#include "geometry/perturbation.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <unordered_set>
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
/// identified by their edge index.
///
/// Use strict equality/inequality for symbolic perturbation consistency.
std::size_t resolve_hit_edge(const Submap& submap,
                              const Polygon& polygon,
                              std::size_t arc_idx,
                              double ray_y) {
    const auto& arc = submap.arc(arc_idx);
    if (arc.first_edge == NONE) return NONE;

    std::size_t lo = std::min(arc.first_edge, arc.last_edge);
    std::size_t hi = std::max(arc.first_edge, arc.last_edge);

    // In the symbolic perturbation model, general position is assumed.
    // A horizontal ray at `ray_y` intersects exactly one edge of the chain
    // [lo, hi], unless ray_y coincides with a vertex y-coordinate.
    // Even then, consistency is handled by the lexicographic order in predicates.
    // Here we just need to find *which* edge covers ray_y.

    for (std::size_t ei = lo; ei <= hi && ei < polygon.num_edges(); ++ei) {
        const auto& edge = polygon.edge(ei);
        const auto& p1 = polygon.vertex(edge.start_idx);
        const auto& p2 = polygon.vertex(edge.end_idx);

        double y_lo = std::min(p1.y, p2.y);
        double y_hi = std::max(p1.y, p2.y);

        // Exact comparison.
        if (ray_y < y_lo || ray_y > y_hi) continue;
        if (p1.y == p2.y) continue; // Horizontal edges shouldn't exist in filtered input.

        return ei;
    }

    return NONE;
}

/// Compute the x-coordinate of a point on ∂C at height y, on the
/// polygon edge incident on vertex v.
double origin_x_for_vertex(const Polygon& polygon, std::size_t v) {
    if (v >= polygon.num_vertices()) return 0.0;
    const auto& pt = polygon.vertex(v);
    // Return vertex x directly.
    return pt.x;
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
        // Exact comparison for chords (they are horizontal segments).
        if (c.y != y) continue;
        
        if (c.left_edge != NONE && c.right_edge != NONE &&
            c.left_edge < polygon.num_edges() &&
            c.right_edge < polygon.num_edges()) {
            double lx = polygon.edge_x_at_y(c.left_edge, c.y);
            double rx = polygon.edge_x_at_y(c.right_edge, c.y);
            double chord_x = shoot_right ? std::min(lx, rx)
                                         : std::max(lx, rx);
            double dist = shoot_right ? (chord_x - origin_x)
                                      : (origin_x - chord_x);
            // Strict distance check > 0 for valid hit.
            if (dist > 0.0 && dist < best_dist) {
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
        // Symbolic perturbation logic:
        // A tilted chord is formally `y = virtual_y + ε * x`.
        // In the limit ε -> 0, it acts as a horizontal segment at `virtual_y`.
        if (a.is_virtual()) {
            double vy = a.virtual_y;
            // Exact y comparison.
            if (y != vy) continue;

            double x_left = (a.first_edge < polygon.num_edges())
                ? polygon.edge_x_at_y(a.first_edge, vy) : 0.0;
            double x_right = (a.last_edge < polygon.num_edges())
                ? polygon.edge_x_at_y(a.last_edge, vy) : 0.0;
            // Midpoint approximation for hit x is fine topology-wise.
            double x = (x_left + x_right) * 0.5;
            double dist = shoot_right ? (x - origin_x) : (origin_x - x);
            if (dist > 0.0 && dist < best_dist) {
                best_dist = dist;
                best.type = RayHit::Type::ARC;
                best.arc_idx = ai;
                best.hit_x = x;
            }
            continue;
        }

        // §3.1: O(1) per arc — test only the two stored endpoint edges.
        // The paper models each arc by its endpoint edge pointers; a
        // horizontal ray intersects the arc boundary at most once per
        // endpoint edge.
        for (std::size_t ei : {a.first_edge, a.last_edge}) {
            if (ei == NONE || ei >= polygon.num_edges()) continue;
            const auto& edge = polygon.edge(ei);
            const auto& p1 = polygon.vertex(edge.start_idx);
            const auto& p2 = polygon.vertex(edge.end_idx);
            
            double y_lo = std::min(p1.y, p2.y);
            double y_hi = std::max(p1.y, p2.y);

            // Exact bounds check.
            if (y < y_lo || y > y_hi) continue;
            if (p1.y == p2.y) continue; 

            // Exact intersection computation.
            double x = horizontal_ray_x_intercept(p1, p2, y);

            double dist = shoot_right ? (x - origin_x) : (origin_x - x);
            // Strict check.
            if (dist > 0.0 && dist < best_dist) {
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
               const RayShootingOracle& dst_oracle,
               const Polygon& polygon,
               std::size_t junction_vertex,
               std::vector<DiscoveredChord>& discovered) {


    // -- Collect stop points --
    // Each stop has an edge index and y-coordinate defining its exact
    // position on ∂C.  Vertex-based stops (arc boundaries) and
    // edge-based stops (chord endpoints) coexist.
    struct StopPoint {
        std::size_t edge_idx;   ///< Edge on which the stop lies.
        double y;               ///< Y-coordinate.
        bool is_vertex;         ///< True if this stop is a polygon vertex.
        std::size_t vertex_idx; ///< Vertex index (if is_vertex), else NONE.
        bool shoot_right;       ///< §3.1: assigned chord direction.
    };
    std::vector<StopPoint> stops;

    // §3.1 Fix 3: O(1) stop deduplication via hash set.
    auto stop_hash = [](std::size_t edge, double y) -> std::size_t {
        std::size_t h = edge;
        std::uint64_t y_bits;
        std::memcpy(&y_bits, &y, sizeof(y_bits));
        h ^= std::hash<std::uint64_t>{}(y_bits) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    };
    std::unordered_set<std::size_t> seen_stops;

    auto add_vertex_stop = [&](std::size_t v, bool shoot_right) {
        if (v == NONE || v >= polygon.num_vertices()) return;
        std::size_t e = v;
        if (e >= polygon.num_edges()) e = polygon.num_edges() - 1;
        double y = polygon.vertex(v).y;
        if (seen_stops.insert(stop_hash(e, y)).second)
            stops.push_back({e, y, true, v, shoot_right});
    };

    auto add_edge_stop = [&](std::size_t edge, double y, bool shoot_right) {
        if (edge == NONE || edge >= polygon.num_edges()) return;
        if (seen_stops.insert(stop_hash(edge, y)).second)
            stops.push_back({edge, y, false, NONE, shoot_right});
    };

    // Chord endpoints of src (may be interior to edges).
    // §3.1: direction is determined by the chord's position relative
    // to the arc — left_edge → shoot right, right_edge → shoot left.
    for (std::size_t ci = 0; ci < src.num_chords(); ++ci) {
        const auto& c = src.chord(ci);
        if (c.left_edge != NONE) add_edge_stop(c.left_edge, c.y, true);
        if (c.right_edge != NONE) add_edge_stop(c.right_edge, c.y, false);
    }

    // Companion vertices: arc boundary vertices.
    // §3.1: direction from the arc's side flags.
    for (std::size_t ai = 0; ai < src.num_arcs(); ++ai) {
        const auto& arc = src.arc(ai);
        if (arc.first_edge == NONE) continue;
        std::size_t lo = std::min(arc.first_edge, arc.last_edge);
        std::size_t hi = std::max(arc.first_edge, arc.last_edge);
        // first_side == LEFT → shoot right from that endpoint.
        add_vertex_stop(lo, arc.first_side == Side::LEFT);
        if (hi + 1 < polygon.num_vertices())
            add_vertex_stop(hi + 1, arc.last_side == Side::LEFT);
    }

    // Always include the junction vertex (default: shoot right).
    add_vertex_stop(junction_vertex, true);

    // Sort by position along ∂C using perturbed comparison logic.
    std::sort(stops.begin(), stops.end(),
              [&](const StopPoint& a, const StopPoint& b) {
                  if (a.edge_idx != b.edge_idx) return a.edge_idx < b.edge_idx;
                  if (a.edge_idx < polygon.num_edges()) {
                      const auto& e = polygon.edge(a.edge_idx);
                      const auto& p1 = polygon.vertex(e.start_idx);
                      const auto& p2 = polygon.vertex(e.end_idx);
                      
                      // Perturbed comparison.
                      bool p1_less_p2 = perturbed_y_less(p1, p2);
                      
                      if (a.y == b.y) return false;

                      if (p1_less_p2) {
                          return a.y < b.y;
                      } else {
                          return a.y > b.y;
                      }
                  }
                  return a.y < b.y;
              });

    if (stops.empty()) return;

    // -- Start-up phase --
    std::size_t current_region = find_junction_region(dst, junction_vertex);

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

        // §3.1: Shoot from a0 in the assigned direction.
        bool a0_dir = stops[a0_idx].shoot_right;

        // Shoot from a0 in src to find what it sees w.r.t. C1.
        double dist_src = std::numeric_limits<double>::infinity();
        {
            RayHit h = src_oracle.shoot(shoot_e, y, Side::LEFT, a0_dir);
            if (h.type == RayHit::Type::ARC) {
                double d = a0_dir ? (h.hit_x - ox) : (ox - h.hit_x);
                if (d > 0.0) dist_src = d;
            }
        }

        // Shoot from a0 within R of dst to find what it sees w.r.t. C2.
        RayHit hit_dst;
        double dist_dst = std::numeric_limits<double>::infinity();
        if (current_region != NONE) {
            RayHit h = dst_oracle.is_built()
                ? dst_oracle.shoot_from_region(current_region, ox, y, a0_dir)
                : local_shoot_in_region(dst, polygon, current_region, ox, y, a0_dir);
            if (h.type == RayHit::Type::ARC) {
                double d = a0_dir ? (h.hit_x - ox) : (ox - h.hit_x);
                if (d > 0.0) {
                    dist_dst = d;
                    hit_dst = h;
                }
            }
        }

        if (hit_dst.type == RayHit::Type::ARC && dist_dst <= dist_src) {
            std::size_t re = resolve_hit_edge(
                dst, polygon, hit_dst.arc_idx, y);
            std::size_t a0_edge = stops[a0_idx].edge_idx;
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
            {
                RayHit h = src_oracle.shoot(shoot_e, y, Side::LEFT, a0_dir);
                if (h.type == RayHit::Type::ARC) {
                    double d = a0_dir ? (h.hit_x - ox) : (ox - h.hit_x);
                    if (d > 0.0) {
                        c0_edge = resolve_hit_edge(
                            src, polygon, h.arc_idx, y);
                    }
                }
            }

            if (c0_edge != NONE && c0_edge != stops[a0_idx].edge_idx) {
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
        bool aj_dir = stop_j.shoot_right;

        // --- Case (i): aj lies in R and sees ∂C₂ ---
        bool aj_in_R = false;
        RayHit best_dst_hit;
        double best_dst_dist = std::numeric_limits<double>::infinity();

        if (current_region < dst.num_nodes() &&
            !dst.node(current_region).deleted) {
            RayHit h = dst_oracle.is_built()
                ? dst_oracle.shoot_from_region(current_region, ox, y, aj_dir)
                : local_shoot_in_region(dst, polygon, current_region, ox, y, aj_dir);
            if (h.type == RayHit::Type::ARC) {
                double d = aj_dir ? (h.hit_x - ox) : (ox - h.hit_x);
                if (d > 0.0) {
                    best_dst_dist = d;
                    best_dst_hit = h;
                    aj_in_R = true;
                }
            }
        }

        if (aj_in_R) {
            double dist_src = std::numeric_limits<double>::infinity();
            {
                RayHit h = src_oracle.shoot(shoot_e, y, Side::LEFT, aj_dir);
                if (h.type == RayHit::Type::ARC) {
                    double d = aj_dir ? (h.hit_x - ox) : (ox - h.hit_x);
                    if (d > 0.0) dist_src = d;
                }
            }

            if (best_dst_dist <= dist_src) {
                std::size_t re = resolve_hit_edge(
                    dst, polygon, best_dst_hit.arc_idx, y);
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

                {
                    RayHit h = src_oracle.shoot(endpt_e, y_e, Side::LEFT, aj_dir);
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
                        src, polygon, h.arc_idx, y_e);
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
                      const std::vector<DiscoveredChord>& discovered_chords,
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
        bool operator==(const ChordKey& o) const {
            return lo == o.lo && hi == o.hi && y == o.y;
        }
    };
    struct ChordKeyHash {
        std::size_t operator()(const ChordKey& k) const {
            std::size_t h = k.lo;
            h ^= std::hash<std::size_t>{}(k.hi) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            std::uint64_t y_bits;
            std::memcpy(&y_bits, &k.y, sizeof(y_bits));
            h ^= std::hash<std::uint64_t>{}(y_bits) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            return h;
        }
    };
    std::unordered_set<ChordKey, ChordKeyHash> seen;
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

    // Populate all_chords from S₁, S₂, and discovered chords.
    // Chords are added via add_chord_rec which handles deduplication.

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

    for (const auto& dc : discovered_chords) {
        add_chord_rec(dc.y, dc.left_edge, dc.right_edge, dc.is_null_length);
    }

    // §3.1: Merge the three chord lists in linear time O(N).
    // S₁ and S₂ chords are guaranteed sorted by Submap::normalize().
    
    // 1. Transform and sort discovered chords.
    std::vector<ChordRecord> disc_list;
    disc_list.reserve(discovered_chords.size());
    for (const auto& dc : discovered_chords) {
         disc_list.push_back({dc.y, dc.left_edge, dc.right_edge, dc.is_null_length, 
                              std::min(dc.left_edge, dc.right_edge)});
    }
    auto chord_less = [](const ChordRecord& a, const ChordRecord& b) {
        if (a.sort_key != b.sort_key) return a.sort_key < b.sort_key;
        return a.y < b.y;
    };
    std::sort(disc_list.begin(), disc_list.end(), chord_less);

    // 2. Pre-filter S1 and S2 chords to handle invalid entries linearly.
    std::vector<ChordRecord> s1_list, s2_list;
    s1_list.reserve(s1.num_chords());
    for (std::size_t k = 0; k < s1.num_chords(); ++k) {
        const auto& c = s1.chord(k);
        if (c.region[0] != NONE || c.region[1] != NONE)
             s1_list.push_back({c.y, c.left_edge, c.right_edge, c.is_null_length, std::min(c.left_edge, c.right_edge)});
    }
    s2_list.reserve(s2.num_chords());
    for (std::size_t k = 0; k < s2.num_chords(); ++k) {
        const auto& c = s2.chord(k);
        if (c.region[0] != NONE || c.region[1] != NONE)
             s2_list.push_back({c.y, c.left_edge, c.right_edge, c.is_null_length, std::min(c.left_edge, c.right_edge)});
    }

    // 3. Perform 3-Way Merge of S1, S2, and Discovered.
    all_chords.clear();
    all_chords.reserve(s1_list.size() + s2_list.size() + disc_list.size());

    std::size_t i1 = 0, n1 = s1_list.size();
    std::size_t i2 = 0, n2 = s2_list.size();
    std::size_t id = 0, nd = disc_list.size();
    
    while (i1 < n1 || i2 < n2 || id < nd) {
        // Find min of head of lists
        int source = -1; // 0=s1, 1=s2, 2=disc
        const ChordRecord* min_c = nullptr;

        if (i1 < n1) {
            min_c = &s1_list[i1];
            source = 0;
        }
        if (i2 < n2) {
            if (!min_c || chord_less(s2_list[i2], *min_c)) {
                min_c = &s2_list[i2];
                source = 1;
            }
        }
        if (id < nd) {
            if (!min_c || chord_less(disc_list[id], *min_c)) {
                min_c = &disc_list[id];
                source = 2;
            }
        }
        
        // Push min and advance
        all_chords.push_back(*min_c);
        if (source == 0) i1++;
        else if (source == 1) i2++;
        else id++;
    }

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
        if (lo >= hi) {
             // No interior points possible if lo == hi.
             all_arcs.push_back({a, lo});
             return;
        }

        auto it_lo = std::upper_bound(split_edges.begin(),
                                       split_edges.end(), lo);
        auto it_hi = std::lower_bound(split_edges.begin(),
                                       split_edges.end(), hi);

        // Collect interior split points.
        std::vector<std::size_t> splits;
        // Safety check: sequence must be valid.
        if (std::distance(it_lo, it_hi) > 0) {
            for (auto it = it_lo; it != it_hi; ++it)
                splits.push_back(*it);
        }

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

    // S₁ arcs come before S₂ arcs in ∂C order (disjoint edge ranges),
    // so concatenation preserves ∂C order — no sort needed.
    for (std::size_t i = 0; i < s1.num_arcs(); ++i)
        split_and_add_arc(s1.arc(i));
    for (std::size_t i = 0; i < s2.num_arcs(); ++i)
        split_and_add_arc(s2.arc(i));

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

} // namespace

Submap fuse(const Submap& s1, const RayShootingOracle& oracle1,
            const Submap& s2, const RayShootingOracle& oracle2,
            const Polygon& polygon,
            std::size_t junction_vertex) {
    if (s1.num_nodes() == 0) return Submap(s2);
    if (s2.num_nodes() == 0) return Submap(s1);

    std::vector<DiscoveredChord> discovered;

    fuse_pass(s1, oracle1, s2, oracle2, polygon, junction_vertex, discovered);
    fuse_pass(s2, oracle2, s1, oracle1, polygon, junction_vertex, discovered);

    return rebuild_submap(s1, s2, discovered, polygon);
}

} // namespace chazelle
