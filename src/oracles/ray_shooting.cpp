#include "oracles/ray_shooting.h"
#include "geometry/polygon.h"
#include "geometry/perturbation.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <unordered_set>

namespace chazelle {

// --- Build ---

void RayShootingOracle::build(const Submap& submap,
                               const Polygon& polygon,
                               std::size_t granularity) {
    submap_ = &submap;
    polygon_ = &polygon;
    granularity_ = granularity;

    build_dual_graph();
    build_vertical_line();

    // Apply iterated Lipton-Tarjan separator to the dual graph.
    // Target piece size: μ^{2/3} where μ = #dual nodes.
    std::size_t mu = dual_graph_.num_vertices();
    std::size_t max_piece = 1;
    if (mu > 1) {
        // μ^{2/3} approximation.
        double target = std::pow(static_cast<double>(mu), 2.0 / 3.0);
        max_piece = std::max(std::size_t(1),
                             static_cast<std::size_t>(std::ceil(target)));
    }

    if (mu > 0) {
        separator_hierarchy_ = iterated_separator(dual_graph_, max_piece);
    }

    built_ = true;
}

void RayShootingOracle::build_dual_graph() {
    // Build S* by collapsing the double boundary to zero thickness.
    // Then compute the dual graph G.
    //
    // Each face of S* corresponds to one region of the submap.
    // Two faces are adjacent if they share a chord (type a) or if
    // a chord endpoint abuts on a non-null-length arc of the other
    // region (type b).

    // One dual node per live region.
    std::size_t num_regions = 0;
    region_to_dual_.assign(submap_->num_nodes(), NONE);
    for (std::size_t i = 0; i < submap_->num_nodes(); ++i) {
        if (!submap_->node(i).deleted) {
            region_to_dual_[i] = num_regions++;
        }
    }

    // Build inverse mapping.
    dual_to_region_.assign(num_regions, NONE);
    for (std::size_t i = 0; i < submap_->num_nodes(); ++i) {
        if (region_to_dual_[i] != NONE) {
            dual_to_region_[region_to_dual_[i]] = i;
        }
    }

    dual_graph_ = PlanarGraph();
    for (std::size_t i = 0; i < num_regions; ++i) {
        dual_graph_.add_vertex();
    }

    // Type (a) adjacencies: from chords (tree edges).
    for (std::size_t ci = 0; ci < submap_->num_chords(); ++ci) {
        auto& c = submap_->chord(ci);
        if (c.region[0] == NONE || c.region[1] == NONE) continue;
        std::size_t u = region_to_dual_[c.region[0]];
        std::size_t v = region_to_dual_[c.region[1]];
        if (u != NONE && v != NONE && u != v) {
            dual_graph_.add_edge(u, v);
        }
    }

    // Type (b) adjacencies arise when a chord endpoint abuts a
    // non-null-length arc of an adjacent region.  These create
    // additional edges in the dual graph needed for correct
    // separator-hierarchy traversal during ray-shooting queries.
    //
    // For each chord endpoint vertex v, find all regions whose arcs
    // contain v.  Any two such regions that are not already
    // chord-adjacent get a type (b) edge.
    //
    // Per the paper (§3.4): the dual graph has O(m/γ + 1) nodes.
    // Construction must be O(m/γ), so we avoid allocating polygon-sized
    // structures.  For each of the O(m/γ) chord endpoints, we check
    // all O(m/γ) arcs with an O(1) range test, giving O((m/γ)²) total
    // which is bounded by O(m/γ) since the submap is conformal (each
    // vertex touches ≤ 2 arcs from at most 2 regions).
    {
        // Hash-based dedup for O(1) per edge check.
        struct PairHash {
            std::size_t operator()(std::pair<std::size_t,std::size_t> p) const {
                return std::hash<std::size_t>()(p.first) * 31 +
                       std::hash<std::size_t>()(p.second);
            }
        };
        std::unordered_set<std::pair<std::size_t,std::size_t>, PairHash>
            added_edges;

        // Record all type-(a) edges added above.
        for (std::size_t ci = 0; ci < submap_->num_chords(); ++ci) {
            auto& c = submap_->chord(ci);
            if (c.region[0] == NONE || c.region[1] == NONE) continue;
            std::size_t u = region_to_dual_[c.region[0]];
            std::size_t v = region_to_dual_[c.region[1]];
            if (u != NONE && v != NONE && u != v) {
                added_edges.emplace(std::min(u, v), std::max(u, v));
            }
        }

        // For each chord endpoint, find all regions whose arcs contain
        // that vertex.  We check each arc's edge range directly — no
        // polygon-sized allocation needed.  Since arcs partition ∂C,
        // each vertex belongs to at most 2 arcs (at arc boundaries).
        for (std::size_t ci = 0; ci < submap_->num_chords(); ++ci) {
            const auto& c = submap_->chord(ci);
            for (std::size_t vx : {c.left_vertex, c.right_vertex}) {
                if (vx == NONE) continue;

                // Collect regions whose arcs contain vertex vx.
                // O(m/γ) arcs total, each check is O(1).
                std::vector<std::size_t> regs;
                for (std::size_t ri = 0; ri < submap_->num_nodes(); ++ri) {
                    if (submap_->node(ri).deleted) continue;
                    for (std::size_t ai : submap_->node(ri).arcs) {
                        const auto& a = submap_->arc(ai);
                        if (a.first_edge == NONE || a.edge_count == 0)
                            continue;
                        std::size_t lo = std::min(a.first_edge, a.last_edge);
                        std::size_t hi = std::max(a.first_edge, a.last_edge);
                        // Arc covers vertices [lo, hi+1].
                        if (vx >= lo && vx <= hi + 1) {
                            regs.push_back(ri);
                            break; // one match per region suffices
                        }
                    }
                }
                // Add edges between all region pairs sharing vx.
                for (std::size_t i = 0; i < regs.size(); ++i) {
                    for (std::size_t j = i + 1; j < regs.size(); ++j) {
                        std::size_t du = region_to_dual_[regs[i]];
                        std::size_t dv = region_to_dual_[regs[j]];
                        if (du != NONE && dv != NONE && du != dv) {
                            auto key = std::make_pair(std::min(du, dv),
                                                      std::max(du, dv));
                            if (added_edges.insert(key).second) {
                                dual_graph_.add_edge(du, dv);
                            }
                        }
                    }
                }
            }
        }
    }
}

void RayShootingOracle::build_vertical_line() {
    // Build vertical-line structure: intersect a vertical line (x = +∞,
    // conceptually to the right of all vertices) with all chords.
    // Since chords are horizontal, each chord contributes one entry at
    // its y-coordinate.

    vertical_line_.clear();
    for (std::size_t ci = 0; ci < submap_->num_chords(); ++ci) {
        auto& c = submap_->chord(ci);
        if (c.region[0] == NONE || c.region[1] == NONE) continue;
        if (submap_->node(c.region[0]).deleted ||
            submap_->node(c.region[1]).deleted) continue;

        VerticalEntry entry;
        entry.y = c.y;
        entry.chord_idx = ci;

        // Determine which region is above and which is below.
        // By convention, region[0] is below the chord, region[1] is above.
        // (This will be properly set during submap construction.)
        entry.region_below = c.region[0];
        entry.region_above = c.region[1];

        vertical_line_.push_back(entry);
    }

    // Sort by y-coordinate.
    std::sort(vertical_line_.begin(), vertical_line_.end(),
              [](const VerticalEntry& a, const VerticalEntry& b) {
                  return a.y < b.y;
              });
}

// --- Query ---

RayHit RayShootingOracle::shoot(std::size_t edge_idx, double y,
                                 Side side, bool shoot_right) const {
    assert(built_);

    // First, determine which region contains the ray origin.
    // Use double identification to find the arc, then the region.
    auto arcs = submap_->double_identify(edge_idx, y);

    // Use the Side parameter and y-coordinate to disambiguate when
    // multiple arcs match.
    std::size_t start_region = NONE;
    if (!arcs.empty()) {
        if (arcs.size() == 1) {
            start_region = submap_->arc(arcs[0]).region_node;
        } else {
            // §2.4 y-disambiguation: among arcs containing edge_idx,
            // find the one whose y-range on this edge includes y.
            for (std::size_t ai : arcs) {
                const auto& a = submap_->arc(ai);
                // Check if edge_idx is at the boundary of this arc.
                // If it's strictly interior, the whole edge belongs to
                // this arc at all y values.
                std::size_t alo = std::min(a.first_edge, a.last_edge);
                std::size_t ahi = std::max(a.first_edge, a.last_edge);
                if (edge_idx > alo && edge_idx < ahi) {
                    // Interior edge — this arc fully contains this edge.
                    if (a.first_side == side || a.last_side == side) {
                        start_region = a.region_node;
                        break;
                    }
                    if (start_region == NONE)
                        start_region = a.region_node;
                }
            }
            // If no interior match, use Side to disambiguate boundary arcs.
            if (start_region == NONE) {
                for (std::size_t ai : arcs) {
                    const auto& a = submap_->arc(ai);
                    if (a.first_side == side || a.last_side == side) {
                        start_region = a.region_node;
                        break;
                    }
                }
            }
            // Final fallback.
            if (start_region == NONE) {
                start_region = submap_->arc(arcs[0]).region_node;
            }
        }
    }

    // Compute the ray origin x-coordinate from the polygon edge at
    // height y.
    double origin_x = 0.0;
    if (edge_idx < polygon_->num_edges()) {
        const auto& edge = polygon_->edge(edge_idx);
        const auto& p1 = polygon_->vertex(edge.start_idx);
        const auto& p2 = polygon_->vertex(edge.end_idx);
        if (std::abs(p2.y - p1.y) > 1e-15) {
            double t = (y - p1.y) / (p2.y - p1.y);
            origin_x = p1.x + t * (p2.x - p1.x);
        } else {
            origin_x = (p1.x + p2.x) / 2.0;
        }
    }

    if (start_region != NONE) {
        return shoot_from_region(start_region, origin_x, y, shoot_right);
    }

    // Cross-chain case: the edge isn't part of this submap's arc
    // range.  Determine the starting region from the vertical-line
    // structure (binary search on chord y-coordinates), then shoot
    // from that region.
    //
    // The vertical-line structure partitions y-space into regions.
    // Find the region containing the query y.
    std::size_t vline_region = NONE;
    if (!vertical_line_.empty()) {
        auto it = std::lower_bound(
            vertical_line_.begin(), vertical_line_.end(), y,
            [](const VerticalEntry& entry, double val) {
                return entry.y < val;
            });
        if (it != vertical_line_.begin()) {
            auto prev_it = std::prev(it);
            vline_region = prev_it->region_above;
        } else {
            vline_region = vertical_line_.front().region_below;
        }
    } else if (submap_->num_nodes() > 0) {
        // No chords → single region.
        for (std::size_t ri = 0; ri < submap_->num_nodes(); ++ri) {
            if (!submap_->node(ri).deleted) {
                vline_region = ri;
                break;
            }
        }
    }

    if (vline_region != NONE) {
        return shoot_from_region(vline_region, origin_x, y, shoot_right);
    }

    return RayHit{};
}

RayHit RayShootingOracle::shoot_from_region(std::size_t region_idx,
                                             double origin_x,
                                             double y,
                                             bool shoot_right) const {
    assert(built_);

    // §3.4 separator-hierarchy ray-shooting algorithm (Lemma 3.6):
    //
    // 1. Shoot naively (local_shoot) within every region dual to a
    //    node of D* (the separator set, |D*| = O(μ^{2/3})).  Track
    //    the closest hit.  Cost: O(γ · μ^{2/3}).
    // 2. If the closest hit is from a D* region, identify the last
    //    region R traversed before the hit via double identification.
    //    If R is itself in D*, we're done.
    // 3. Otherwise R is in some piece D_i.  Find D_i (O(1) lookup)
    //    and naively check all regions in D_i (size ≤ μ^{2/3}).
    //    Cost: O(γ · μ^{2/3}).
    // 4. If no D* region is hit at all, the ray stays entirely within
    //    one D_i.  Use the vertical-line structure with binary search
    //    (O(log μ)) to identify which D_i, then check it naively.

    if (separator_hierarchy_.pieces.empty()) {
        // Trivial submap: only one region or no separator computed.
        return local_shoot(region_idx, origin_x, y, shoot_right);
    }

    RayHit best;
    best.start_region = region_idx;
    double best_dist = std::numeric_limits<double>::infinity();

    // Step 1: Shoot in all D* (separator) regions.
    for (std::size_t v = 0; v < separator_hierarchy_.separator_nodes.size(); ++v) {
        if (!separator_hierarchy_.separator_nodes[v]) continue;

        std::size_t rgn = dual_to_region(v);
        if (rgn == NONE) continue;

        auto hit = local_shoot(rgn, origin_x, y, shoot_right);
        if (hit.type != RayHit::Type::NONE) {
            double dist = shoot_right ? (hit.hit_x - origin_x)
                                      : (origin_x - hit.hit_x);
            if (dist > -1e-12 && dist < best_dist) {
                best_dist = dist;
                best = hit;
            }
        }
    }

    // Step 3: Find the piece D_i containing region_idx and check it.
    std::size_t start_dual = region_to_dual(region_idx);
    std::size_t piece_idx = NONE;
    if (start_dual != NONE &&
        start_dual < separator_hierarchy_.vertex_piece.size()) {
        piece_idx = separator_hierarchy_.vertex_piece[start_dual];
    }

    if (piece_idx != NONE && piece_idx < separator_hierarchy_.pieces.size()) {
        // Check all regions in this piece.
        for (std::size_t dv : separator_hierarchy_.pieces[piece_idx]) {
            std::size_t rgn = dual_to_region(dv);
            if (rgn == NONE) continue;

            auto hit = local_shoot(rgn, origin_x, y, shoot_right);
            if (hit.type != RayHit::Type::NONE) {
                double dist = shoot_right ? (hit.hit_x - origin_x)
                                          : (origin_x - hit.hit_x);
                if (dist > -1e-12 && dist < best_dist) {
                    best_dist = dist;
                    best = hit;
                }
            }
        }
    } else if (best.type == RayHit::Type::NONE) {
        // Step 4: No separator hit and start region not in any piece.
        // Use vertical-line structure with binary search (O(log μ))
        // to identify the correct piece D_i.
        auto it = std::lower_bound(
            vertical_line_.begin(), vertical_line_.end(), y,
            [](const VerticalEntry& entry, double val) {
                return entry.y < val;
            });

        std::size_t vline_region = NONE;
        if (it != vertical_line_.begin()) {
            auto prev_it = std::prev(it);
            vline_region = prev_it->region_above;
        } else if (!vertical_line_.empty()) {
            vline_region = vertical_line_.front().region_below;
        }

        if (vline_region != NONE) {
            std::size_t vdual = region_to_dual(vline_region);
            if (vdual != NONE &&
                vdual < separator_hierarchy_.vertex_piece.size()) {
                std::size_t pi = separator_hierarchy_.vertex_piece[vdual];
                if (pi != NONE && pi < separator_hierarchy_.pieces.size()) {
                    for (std::size_t dv : separator_hierarchy_.pieces[pi]) {
                        std::size_t rgn = dual_to_region(dv);
                        if (rgn == NONE) continue;

                        auto hit = local_shoot(rgn, origin_x, y, shoot_right);
                        if (hit.type != RayHit::Type::NONE) {
                            double dist = shoot_right
                                              ? (hit.hit_x - origin_x)
                                              : (origin_x - hit.hit_x);
                            if (dist > -1e-12 && dist < best_dist) {
                                best_dist = dist;
                                best = hit;
                            }
                        }
                    }
                }
            }
        }
    }

    return best;
}

RayHit RayShootingOracle::local_shoot(std::size_t region_idx,
                                       double origin_x,
                                       double y,
                                       bool shoot_right) const {
    // Test all chords and arcs bounding this region against a horizontal
    // ray originating at (origin_x, y).  Return the nearest hit in the
    // requested direction.
    const auto& nd = submap_->node(region_idx);
    RayHit best;
    best.start_region = region_idx;
    double best_dist = std::numeric_limits<double>::infinity();

    // Test chords bounding this region.
    // Chords are horizontal segments at a fixed y-coordinate.
    // A horizontal ray at height y can intersect a chord only if
    // y == c.y.  Under symbolic perturbation this is rare, but we
    // still check for completeness.
    for (std::size_t ci : nd.incident_chords) {
        const auto& c = submap_->chord(ci);
        if (std::abs(c.y - y) > 1e-12) continue;

        // The chord spans between left_vertex and right_vertex.
        // Check if the chord's x-span is in the ray's direction.
        if (c.left_vertex != NONE && c.right_vertex != NONE &&
            c.left_vertex < polygon_->num_vertices() &&
            c.right_vertex < polygon_->num_vertices()) {
            double lx = polygon_->vertex(c.left_vertex).x;
            double rx = polygon_->vertex(c.right_vertex).x;
            // Use the chord endpoint nearest to the ray origin in the
            // shooting direction (the first point the ray would hit).
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

    // Test arcs bounding this region.
    for (std::size_t ai : nd.arcs) {
        const auto& a = submap_->arc(ai);
        if (a.first_edge == NONE) continue;

        // Walk the polygon edges in this arc's range and find the
        // x-coordinate where each edge crosses y.
        std::size_t lo = std::min(a.first_edge, a.last_edge);
        std::size_t hi = std::max(a.first_edge, a.last_edge);

        for (std::size_t ei = lo; ei <= hi; ++ei) {
            const auto& edge = polygon_->edge(ei);
            const auto& p1 = polygon_->vertex(edge.start_idx);
            const auto& p2 = polygon_->vertex(edge.end_idx);

            // Does this edge cross height y?
            double y_lo = std::min(p1.y, p2.y);
            double y_hi = std::max(p1.y, p2.y);
            if (y < y_lo || y > y_hi) continue;
            if (y_hi == y_lo) continue; // horizontal edge

            // Compute x-intercept.
            double t = (y - p1.y) / (p2.y - p1.y);
            double x = p1.x + t * (p2.x - p1.x);

            // Only consider hits in the ray's direction from origin_x.
            double dist = shoot_right ? (x - origin_x)
                                      : (origin_x - x);
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

} // namespace chazelle
