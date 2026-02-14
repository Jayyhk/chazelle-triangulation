#include "oracles/ray_shooting.h"
#include "geometry/polygon.h"
#include "geometry/perturbation.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>
#include <vector>

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

        // Build compact separator list: O(μ^{2/3}) entries.
        // This avoids scanning the full O(μ) boolean array per query.
        separator_list_.clear();
        for (std::size_t v = 0; v < separator_hierarchy_.separator_nodes.size(); ++v) {
            if (separator_hierarchy_.separator_nodes[v])
                separator_list_.push_back(v);
        }
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
    // Per the paper (§3.4): "The latter can be done by double
    // identification, as discussed in Section 2.4, followed by
    // sorting along ∂C, which takes O(μ log m) time."  We use
    // double_identify (O(log m) per endpoint) for each of the O(μ)
    // chord endpoints, giving O(μ log m) total.
    {
        // §3.4: deterministic dedup for dual-graph edges.
        // Use std::set for O(log μ) per insertion/lookup, giving
        // O(μ log μ) total — deterministic and within the paper's
        // O(μ log m) bound.
        std::set<std::pair<std::size_t,std::size_t>> added_edges;

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
        // that vertex via double identification (§2.4).  Per the paper
        // (§3.4): "The latter can be done by double identification, as
        // discussed in Section 2.4, followed by sorting along ∂C, which
        // takes O(μ log m) time."  Each double_identify call is O(log m)
        // and there are O(μ) chord endpoints, giving O(μ log m) total.
        for (std::size_t ci = 0; ci < submap_->num_chords(); ++ci) {
            const auto& c = submap_->chord(ci);
            for (std::size_t ex : {c.left_edge, c.right_edge}) {
                if (ex == NONE) continue;

                // Use double_identify on the chord's endpoint edge
                // to find all arcs passing through this point.
                // Per §3.1 Remark 1, endpoints are on edges, not
                // vertices.  O(log m) binary search.
                std::size_t edge_for_ex = ex;
                if (edge_for_ex >= polygon_->num_edges() && edge_for_ex > 0)
                    edge_for_ex = polygon_->num_edges() - 1;

                auto arcs = submap_->double_identify(edge_for_ex, c.y);

                // Also check the adjacent edge (the endpoint sits
                // at the boundary between edges).
                if (ex > 0) {
                    auto arcs2 = submap_->double_identify(ex - 1, c.y);
                    for (std::size_t ai2 : arcs2) {
                        bool dup = false;
                        for (std::size_t ai : arcs)
                            if (ai == ai2) { dup = true; break; }
                        if (!dup) arcs.push_back(ai2);
                    }
                }

                // Collect distinct regions from the identified arcs.
                std::vector<std::size_t> regs;
                for (std::size_t ai : arcs) {
                    std::size_t rn = submap_->arc(ai).region_node;
                    if (rn == NONE || submap_->node(rn).deleted) continue;
                    bool dup = false;
                    for (std::size_t r : regs)
                        if (r == rn) { dup = true; break; }
                    if (!dup) regs.push_back(rn);
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
    // §3.4: "Take a vertical line passing to the right of all the
    // vertices of P, and intersect it with the chords of the regions
    // in S.  This breaks up the line into segments, every one of which
    // falls entirely within some region; to split up the line and
    // identify the regions cut by each segment can be done by
    // traversing G and checking each chord for intersection with the
    // vertical line.  Since the regions cut correspond to nodes of G
    // lying on a path, sorting the intersections comes for free, and
    // all the work can be done in O(μ) time."
    //
    // We traverse the submap tree via DFS, rooted at the topmost
    // region (above the chord with the highest y).  At each node we
    // visit children in decreasing y-order of their connecting chord.
    // Since each conformal region has O(1) incident chords (degree ≤ 4),
    // this is O(1) per node, giving O(μ) total.

    vertical_line_.clear();

    // Collect valid chords with their y-coordinates.
    struct ChordInfo {
        std::size_t chord_idx;
        double y;
        std::size_t region_below; // region[0]
        std::size_t region_above; // region[1]
    };
    std::vector<ChordInfo> chords;
    for (std::size_t ci = 0; ci < submap_->num_chords(); ++ci) {
        auto& c = submap_->chord(ci);
        if (c.region[0] == NONE || c.region[1] == NONE) continue;
        if (submap_->node(c.region[0]).deleted ||
            submap_->node(c.region[1]).deleted) continue;
        chords.push_back({ci, c.y, c.region[0], c.region[1]});
    }

    if (chords.empty()) return;

    // Build per-region incidence lists (indices into 'chords').
    // Use region_to_dual_ for compact indexing.
    std::size_t num_dual = dual_to_region_.size();
    std::vector<std::vector<std::size_t>> incident(num_dual);
    for (std::size_t i = 0; i < chords.size(); ++i) {
        std::size_t du0 = region_to_dual(chords[i].region_below);
        std::size_t du1 = region_to_dual(chords[i].region_above);
        if (du0 != NONE) incident[du0].push_back(i);
        if (du1 != NONE) incident[du1].push_back(i);
    }

    // Sort each incidence list by decreasing y.  Since each conformal
    // region has degree ≤ 4, each list has at most 3 items → O(1).
    for (auto& inc : incident) {
        std::sort(inc.begin(), inc.end(),
                  [&](std::size_t a, std::size_t b) {
                      return chords[a].y > chords[b].y;
                  });
    }

    // Find the topmost region: above the chord with maximum y.
    std::size_t top_ci = 0;
    for (std::size_t i = 1; i < chords.size(); ++i) {
        if (chords[i].y > chords[top_ci].y)
            top_ci = i;
    }
    std::size_t top_region = chords[top_ci].region_above;

    // DFS: at each region, visit incident chords in decreasing y,
    // emitting each chord and recursing into the other region.
    // This produces chords in decreasing y-order.
    //
    // We use an explicit stack to handle deep trees without
    // risking stack overflow.  Each stack frame tracks the region
    // and its position in the sorted incidence list.
    std::vector<bool> visited(chords.size(), false);
    vertical_line_.reserve(chords.size());

    struct Frame {
        std::size_t region;
        std::size_t pos; // index into incident[dual_id]
    };
    std::vector<Frame> stack;
    stack.push_back({top_region, 0});

    while (!stack.empty()) {
        auto& fr = stack.back();
        std::size_t du = region_to_dual(fr.region);
        if (du == NONE) { stack.pop_back(); continue; }

        const auto& inc = incident[du];
        // Advance past visited chords.
        while (fr.pos < inc.size() && visited[inc[fr.pos]])
            ++fr.pos;

        if (fr.pos >= inc.size()) {
            // All incident chords visited; backtrack.
            stack.pop_back();
            continue;
        }

        std::size_t ci = inc[fr.pos];
        visited[ci] = true;
        ++fr.pos; // advance for when we return to this frame

        VerticalEntry entry;
        entry.y = chords[ci].y;
        entry.chord_idx = chords[ci].chord_idx;
        entry.region_below = chords[ci].region_below;
        entry.region_above = chords[ci].region_above;
        vertical_line_.push_back(entry);

        // Recurse into the other region.
        std::size_t other = (chords[ci].region_below == fr.region)
                                ? chords[ci].region_above
                                : chords[ci].region_below;
        stack.push_back({other, 0});
    }

    // The DFS produced chords in decreasing y-order;
    // reverse to get increasing y for binary-search queries.
    std::reverse(vertical_line_.begin(), vertical_line_.end());
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
            // §2.4: Side-based disambiguation should always find a
            // match in a well-formed normal-form submap.
            if (start_region == NONE && !arcs.empty()) {
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

RayHit RayShootingOracle::shoot_from_point(double origin_x, double y,
                                            bool shoot_right) const {
    assert(built_);

    // §4.1 oracle (i): shoot from a point that is EXTERNAL to
    // this submap's chain C_μ.  The point lies outside C_μ's double
    // boundary, so double_identify would return empty.  Instead, use
    // the vertical-line structure to identify the starting region at
    // height y, then call shoot_from_region with the caller-supplied
    // origin_x.
    //
    // The vertical-line structure partitions the y-axis by chord
    // y-coordinates.  For a point at height y, the region at
    // (x = ∞, y) is found by binary search.  Since the point is
    // exterior to C_μ, this correctly identifies the region on the
    // boundary-facing side of C_μ closest to the point.
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

    // Helper: naively check all regions in a piece, updating best.
    auto check_piece = [&](std::size_t pi) {
        if (pi == NONE || pi >= separator_hierarchy_.pieces.size()) return;
        for (std::size_t dv : separator_hierarchy_.pieces[pi]) {
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
    };

    // Step 1: Shoot in all D* (separator) regions.
    // Use the compact separator_list_ (O(μ^{2/3}) entries) instead
    // of scanning the full O(μ) boolean array.
    for (std::size_t v : separator_list_) {
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

    if (best.type != RayHit::Type::NONE) {
        // Step 2: A D* hit was found.  Per §3.4: "Let R be the last
        // region of S traversed before the first hit.  To identify R
        // can be done by double identification, followed by checking
        // the local orientation of the hit."
        //
        // R is the region on the ray-origin side of the hit boundary
        // element.  We find it by double-identifying at the hit point.
        // If the hit is on an arc, we identify which edge was hit and
        // look up the arc on the opposite side of ∂C (the side facing
        // the ray).  If on a chord, R is whichever of the chord's two
        // regions faces the ray origin.

        std::size_t R = NONE;
        std::vector<std::size_t> R_candidates;

        if (best.type == RayHit::Type::ARC && best.arc_idx != NONE) {
            // The hit arc belongs to a D* region (the "far" side).
            // R is on the "near" side — the region the ray was
            // travelling through when it hit.  Use double_identify
            // on the hit edge to find the companion arc on the
            // opposite side of ∂C.  Per the paper: "To identify R
            // can be done by double identification, followed by
            // checking the local orientation of the hit."
            // The hit edge is already tracked by local_shoot per §3
            // item (i): O(log m) via double_identify.
            if (best.hit_edge != NONE) {
                auto arcs = submap_->double_identify(best.hit_edge, y);
                // Among the identified arcs, R is the one that is
                // NOT the hit arc (i.e., the companion on the
                // other side of ∂C).
                for (std::size_t ai : arcs) {
                    if (ai == best.arc_idx) continue;
                    std::size_t candidate = submap_->arc(ai).region_node;
                    if (candidate != NONE && !submap_->node(candidate).deleted) {
                        R = candidate;
                        break;
                    }
                }
                // If double_identify returned only one arc (the
                // hit arc itself), R is that arc's own region —
                // meaning the ray started inside the same region.
                if (R == NONE && !arcs.empty()) {
                    R = submap_->arc(arcs[0]).region_node;
                }
            }
        } else if (best.type == RayHit::Type::CHORD && best.chord_idx != NONE) {
            // The chord separates region[0] and region[1].  Per the
            // paper: "To identify R can be done by double
            // identification, followed by checking the local
            // orientation of the hit."  For a horizontal ray hitting
            // a horizontal chord, the local orientation is degenerate.
            // We treat both adjacent regions as candidates for R and
            // check both their pieces below.
            const auto& c = submap_->chord(best.chord_idx);
            if (c.region[0] != NONE &&
                !submap_->node(c.region[0]).deleted)
                R_candidates.push_back(c.region[0]);
            if (c.region[1] != NONE &&
                !submap_->node(c.region[1]).deleted)
                R_candidates.push_back(c.region[1]);
        }

        if (R != NONE) R_candidates.push_back(R);
        if (R_candidates.empty()) R_candidates.push_back(region_idx);

        // For each candidate R, check if it is in D*.  If not,
        // find its piece D_i and naively check all regions in D_i.
        // Per §3.4: "If R is a region dual to a node v of D*, then
        // the starting point of the ray lies in R (otherwise an
        // earlier hit would have been detected) and we are trivially
        // done."  Otherwise: "We can find w [...] by first finding
        // D_i, which takes constant time since we know R, and then
        // naively checking all the regions dual to nodes in D_i."
        std::set<std::size_t> checked_pieces;
        for (std::size_t cand_R : R_candidates) {
            std::size_t cand_dual = region_to_dual(cand_R);
            bool cand_in_Dstar =
                (cand_dual != NONE &&
                 cand_dual < separator_hierarchy_.separator_nodes.size() &&
                 separator_hierarchy_.separator_nodes[cand_dual]);

            if (!cand_in_Dstar) {
                std::size_t cand_piece = NONE;
                if (cand_dual != NONE &&
                    cand_dual < separator_hierarchy_.vertex_piece.size()) {
                    cand_piece = separator_hierarchy_.vertex_piece[cand_dual];
                }
                if (cand_piece != NONE &&
                    checked_pieces.insert(cand_piece).second) {
                    check_piece(cand_piece);
                }
            }
        }
    } else {
        // Step 4: No D* hit at all.  Per §3.4: "Then the ray-
        // shooting takes place entirely within the regions dual to
        // the nodes of a single D_i.  To find out which one, we
        // shoot toward the vertical line and find which segment of
        // the line is hit.  This takes O(log μ) time by binary
        // search.  We can now identify the region R immediately.
        // The remainder of the algorithm is unchanged."
        //
        // Binary search the vertical-line structure to find the
        // region at height y, then check its piece D_i.
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
            std::size_t vdual = region_to_dual(vline_region);
            std::size_t vpiece = NONE;
            if (vdual != NONE &&
                vdual < separator_hierarchy_.vertex_piece.size()) {
                vpiece = separator_hierarchy_.vertex_piece[vdual];
            }
            check_piece(vpiece);
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

        // The chord spans between left_edge and right_edge.
        // Check if the chord's x-span is in the ray's direction.
        if (c.left_edge != NONE && c.right_edge != NONE &&
            c.left_edge < polygon_->num_edges() &&
            c.right_edge < polygon_->num_edges()) {
            double lx = polygon_->edge_x_at_y(c.left_edge, c.y);
            double rx = polygon_->edge_x_at_y(c.right_edge, c.y);
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

        // §4.2: Virtual arcs represent tilted exit-chord edges.
        // Symbolic perturbation: the tilted edge has infinitesimal
        // y-extent [vy - ε, vy + ε].  In the limit ε → 0, only rays
        // at y ≈ vy intersect it, and the intersection point is the
        // midpoint of the two endpoint x-coordinates.
        if (a.is_virtual()) {
            double vy = a.virtual_y;
            if (std::abs(y - vy) > 1e-9) continue;
            double x_left = (a.first_edge < polygon_->num_edges())
                ? polygon_->edge_x_at_y(a.first_edge, vy) : 0.0;
            double x_right = (a.last_edge < polygon_->num_edges())
                ? polygon_->edge_x_at_y(a.last_edge, vy) : 0.0;
            double x = (x_left + x_right) * 0.5;
            double dist = shoot_right ? (x - origin_x)
                                      : (origin_x - x);
            if (dist > -1e-12 && dist < best_dist) {
                best_dist = dist;
                best.type = RayHit::Type::ARC;
                best.arc_idx = ai;
                best.hit_x = x;
                best.hit_edge = a.first_edge;
            }
            continue;
        }

        // §3.4: Iterate ALL edges of the arc [first_edge..last_edge].
        // Each arc has at most O(γ) edges by granularity; conformality
        // limits us to ≤4 arcs per region, so total work is O(γ).
        std::size_t e_lo = std::min(a.first_edge, a.last_edge);
        std::size_t e_hi = std::max(a.first_edge, a.last_edge);
        for (std::size_t ei = e_lo; ei <= e_hi && ei < polygon_->num_edges(); ++ei) {
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
                best.hit_edge = ei;
            }
        }
    }

    return best;
}

} // namespace chazelle
