#include "chazelle/conformality.h"
#include "oracles/arc_cutting.h"
#include "visibility/tree_decomposition.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <queue>
#include <vector>

namespace chazelle {

namespace {

struct ChordCandidate {
    double y = 0.0;
    std::size_t left_edge = NONE;
    std::size_t right_edge = NONE;
    std::size_t region_above = NONE;
    std::size_t region_below = NONE;
    bool valid = false;
};

/// Test whether a hit from vertex v on A₁ lands on arc A₂.
/// We check that the hit edge index lies within A₂'s edge range.
bool hit_lands_on_arc(const RayHit& hit, const Submap& submap,
                      const ArcStructure& arc_j) {
    if (hit.type == RayHit::Type::ARC && hit.arc_idx != NONE) {
        const auto& hit_arc = submap.arc(hit.arc_idx);
        // Same edge range as arc_j means the hit is on A₂.
        std::size_t j_lo = std::min(arc_j.first_edge, arc_j.last_edge);
        std::size_t j_hi = std::max(arc_j.first_edge, arc_j.last_edge);
        std::size_t h_lo = std::min(hit_arc.first_edge, hit_arc.last_edge);
        std::size_t h_hi = std::max(hit_arc.first_edge, hit_arc.last_edge);
        // The hit arc overlaps A₂ if the ranges overlap.
        return h_lo <= j_hi && h_hi >= j_lo;
    }
    return false;
}

/// Perform local shooting from a point (origin_x, y) in region R of the
/// CURRENT-LEVEL submap (not the subarc's canonical submap).
/// Returns the first arc/chord hit and its x-coordinate.
RayHit shoot_in_region(const Submap& submap,
                       std::size_t region_idx,
                       const Polygon& polygon,
                       double origin_x, double y,
                       bool shoot_right) {
    // Naive local shoot: walk all edges of arcs/chords in this region.
    const auto& nd = submap.node(region_idx);
    RayHit best;
    best.start_region = region_idx;
    double best_dist = std::numeric_limits<double>::infinity();

    // Test chords.
    for (std::size_t ci : nd.incident_chords) {
        const auto& c = submap.chord(ci);
        if (std::abs(c.y - y) > 1e-12) continue;
        if (c.left_edge != NONE && c.right_edge != NONE &&
            c.left_edge < polygon.num_edges() &&
            c.right_edge < polygon.num_edges()) {
            double lx = polygon.edge_x_at_y(c.left_edge, c.y);
            double rx = polygon.edge_x_at_y(c.right_edge, c.y);
            // Use the chord endpoint nearest in the shooting direction.
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

    // Test arcs.
    for (std::size_t ai : nd.arcs) {
        const auto& a = submap.arc(ai);
        if (a.first_edge == NONE) continue;

        // §4.2: Virtual arcs represent tilted exit-chord edges.
        if (a.is_virtual()) {
            double vy = a.virtual_y;
            constexpr double TILT = 1e-8;
            double vy_lo = vy - TILT;
            double vy_hi = vy + TILT;
            if (y < vy_lo - 1e-12 || y > vy_hi + 1e-12) continue;
            double x_left = (a.first_edge < polygon.num_edges())
                ? polygon.edge_x_at_y(a.first_edge, vy) : 0.0;
            double x_right = (a.last_edge < polygon.num_edges())
                ? polygon.edge_x_at_y(a.last_edge, vy) : 0.0;
            double t = (std::abs(vy_hi - vy_lo) > 1e-15)
                ? (y - vy_lo) / (vy_hi - vy_lo) : 0.5;
            double x = x_left + t * (x_right - x_left);
            double dist = shoot_right ? (x - origin_x)
                                      : (origin_x - x);
            if (dist > -1e-12 && dist < best_dist) {
                best_dist = dist;
                best.type = RayHit::Type::ARC;
                best.arc_idx = ai;
                best.hit_x = x;
            }
            continue;
        }

        std::size_t lo = std::min(a.first_edge, a.last_edge);
        std::size_t hi = std::max(a.first_edge, a.last_edge);
        for (std::size_t ei = lo; ei <= hi; ++ei) {
            if (ei >= polygon.num_edges()) break;
            const auto& edge = polygon.edge(ei);
            const auto& p1 = polygon.vertex(edge.start_idx);
            const auto& p2 = polygon.vertex(edge.end_idx);
            double y_lo = std::min(p1.y, p2.y);
            double y_hi = std::max(p1.y, p2.y);
            if (y < y_lo || y > y_hi) continue;
            if (y_hi == y_lo) continue;
            double t = (y - p1.y) / (p2.y - p1.y);
            double x = p1.x + t * (p2.x - p1.x);
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

/// §3.2 Lemma 2.4 topological test.
///
/// Given centroid chord ab of Sα, with subarc α having endpoints c, d:
///   - If a or b lies on α, shoot from it; if the ray hits A₂, success.
///   - Otherwise determine which side of α (α₁ or α₂, the two pieces
///     of ∂ᾱ between a and b) is shielded from A₂.
///
/// Returns:
///   -1 → left child subtree should be explored (right side shielded)
///    0 → success: found a point that sees A₂ (candidate filled in)
///    1 → right child subtree should be explored (left side shielded)
///    2 → indeterminate (fallback to exhaustive)
int lemma24_test(
    const Submap& parent_submap,
    std::size_t region_idx,
    const Polygon& polygon,
    const Submap& cs_submap,
    const Chord& centroid_chord,
    const ArcStructure& arc_j,
    std::size_t alpha_start,  // start vertex of subarc α
    std::size_t alpha_end,    // end vertex of subarc α
    ChordCandidate& cand) {

    // Centroid chord endpoints: edge indices a, b.
    std::size_t ea = centroid_chord.left_edge;
    std::size_t eb = centroid_chord.right_edge;
    if (ea == NONE || eb == NONE) return 2;
    if (ea >= polygon.num_edges() || eb >= polygon.num_edges())
        return 2;

    // Check if a lies on α = [alpha_start, alpha_end).
    bool a_on_alpha = (ea >= alpha_start && ea < alpha_end);
    bool b_on_alpha = (eb >= alpha_start && eb < alpha_end);

    // If a is on α, shoot from a to see if it hits A₂.
    if (a_on_alpha) {
        // Compute origin_x on the polygon edge ea at centroid chord y.
        // The chord endpoint a lies on edge ea at height centroid_chord.y,
        // NOT at the vertex's own y-coordinate.
        double pa_y = centroid_chord.y;
        double origin_x = polygon.edge_x_at_y(ea, pa_y);
        // Shoot in both directions from a in the parent submap's region.
        for (bool dir : {true, false}) {
            auto hit = shoot_in_region(parent_submap, region_idx, polygon,
                                       origin_x, pa_y, dir);
            if (hit_lands_on_arc(hit, parent_submap, arc_j)) {
                cand.y = pa_y;
                cand.left_edge = ea;
                // Resolve right_edge: find the edge on the hit
                // arc whose x-intercept at y is closest to hit_x.
                cand.right_edge = NONE;
                if (hit.arc_idx != NONE && hit.arc_idx < parent_submap.num_arcs()) {
                    const auto& ha = parent_submap.arc(hit.arc_idx);
                    if (ha.first_edge != NONE) {
                        std::size_t hlo = std::min(ha.first_edge, ha.last_edge);
                        std::size_t hhi = std::max(ha.first_edge, ha.last_edge);
                        double best_xd = std::numeric_limits<double>::infinity();
                        for (std::size_t ei2 = hlo; ei2 <= hhi && ei2 < polygon.num_edges(); ++ei2) {
                            const auto& e2 = polygon.edge(ei2);
                            const auto& ep1 = polygon.vertex(e2.start_idx);
                            const auto& ep2 = polygon.vertex(e2.end_idx);
                            double ey_lo = std::min(ep1.y, ep2.y);
                            double ey_hi = std::max(ep1.y, ep2.y);
                            if (pa_y < ey_lo - 1e-12 || pa_y > ey_hi + 1e-12) continue;
                            if (std::abs(ey_hi - ey_lo) < 1e-15) continue;
                            double t2 = (pa_y - ep1.y) / (ep2.y - ep1.y);
                            double ex = ep1.x + t2 * (ep2.x - ep1.x);
                            double dx = std::abs(ex - hit.hit_x);
                            if (dx < best_xd) {
                                best_xd = dx;
                                cand.right_edge = ei2;
                            }
                        }
                        // No vertex fallback needed; right_edge is set to the hit edge.
                    }
                }
                if (cand.right_edge != NONE) {
                    cand.valid = true;
                    return 0; // success
                }
            }
        }
    }

    // If b is on α, shoot from b to see if it hits A₂.
    if (b_on_alpha) {
        // Compute origin_x on the polygon edge eb at centroid chord y.
        // The chord endpoint b lies on edge eb at height centroid_chord.y,
        // NOT at the vertex's own y-coordinate.
        double pb_y = centroid_chord.y;
        double origin_x = polygon.edge_x_at_y(eb, pb_y);
        for (bool dir : {true, false}) {
            auto hit = shoot_in_region(parent_submap, region_idx, polygon,
                                       origin_x, pb_y, dir);
            if (hit_lands_on_arc(hit, parent_submap, arc_j)) {
                cand.y = pb_y;
                cand.left_edge = eb;
                cand.right_edge = NONE;
                if (hit.arc_idx != NONE && hit.arc_idx < parent_submap.num_arcs()) {
                    const auto& ha = parent_submap.arc(hit.arc_idx);
                    if (ha.first_edge != NONE) {
                        std::size_t hlo = std::min(ha.first_edge, ha.last_edge);
                        std::size_t hhi = std::max(ha.first_edge, ha.last_edge);
                        double best_xd = std::numeric_limits<double>::infinity();
                        for (std::size_t ei2 = hlo; ei2 <= hhi && ei2 < polygon.num_edges(); ++ei2) {
                            const auto& e2 = polygon.edge(ei2);
                            const auto& ep1 = polygon.vertex(e2.start_idx);
                            const auto& ep2 = polygon.vertex(e2.end_idx);
                            double ey_lo = std::min(ep1.y, ep2.y);
                            double ey_hi = std::max(ep1.y, ep2.y);
                            if (pb_y < ey_lo - 1e-12 || pb_y > ey_hi + 1e-12) continue;
                            if (std::abs(ey_hi - ey_lo) < 1e-15) continue;
                            double t2 = (pb_y - ep1.y) / (ep2.y - ep1.y);
                            double ex = ep1.x + t2 * (ep2.x - ep1.x);
                            double dx = std::abs(ex - hit.hit_x);
                            if (dx < best_xd) {
                                best_xd = dx;
                                cand.right_edge = ei2;
                            }
                        }
                        // No vertex fallback needed; right_edge is set to the hit edge.
                    }
                }
                if (cand.right_edge != NONE) {
                    cand.valid = true;
                    return 0; // success
                }
            }
        }
    }

    // Neither a' nor b' hits A₂.  Use Lemma 2.4 to identify the
    // shielded side.  The chord ab splits ∂ᾱ into two halves:
    //   α₁ = portion of ∂ᾱ between a and b going one way
    //   α₂ = portion going the other way
    // We identify which half is "shielded" from A₂ — i.e., no point
    // of that half can see A₂ without crossing ab or A.
    //
    // The centroid chord ab splits the submap into two subtrees:
    //   left_child  (region[0] side of chord)
    //   right_child (region[1] side of chord)
    // The shielded side corresponds to the subtree we should NOT
    // explore.
    //
    // Per the paper: we determine which Bⱼ is shielded by examining
    // the relative positions of a, b, c, d (endpoints of α) and the
    // edge ranges of A₂.  The side of ab whose vertex indices are
    // on the SAME side as A₂'s edge range (in the circular ordering)
    // is the non-shielded side.

    // Determine which side of the chord contains edges overlapping
    // with A₂ (arc_j).  The chord's region[0] and region[1] give us
    // the two sides.  We check which side's arcs overlap arc_j's
    // edge range.
    std::size_t j_lo = std::min(arc_j.first_edge, arc_j.last_edge);
    std::size_t j_hi = std::max(arc_j.first_edge, arc_j.last_edge);

    // Exact topological test: in the linear ordering of edges
    // along ∂C, if ea < eb, then edges in [ea, eb) go to one
    // subtree and those outside go to the other.  If arc_j's edge
    // range overlaps [ea, eb), the shielded side is the other half,
    // so descend into the [ea, eb) side.
    //
    // This is exact (not a heuristic) because ∂C is an open curve
    // with a linear edge ordering, not circular.
    std::size_t c_lo = std::min(ea, eb);
    std::size_t c_hi = std::max(ea, eb);

    // Does A₂'s edge range overlap the [c_lo, c_hi) interval?
    bool j_overlaps_inner = (j_lo < c_hi && j_hi >= c_lo);

    if (j_overlaps_inner) {
        // A₂ faces the "inner" interval [c_lo, c_hi).
        // The "outer" interval is shielded → descend into inner.
        // If the centroid chord's left subtree corresponds to the
        // inner interval, descend left; else descend right.
        // We use the chord's region[0] arcs to determine:
        std::size_t r0 = centroid_chord.region[0];
        if (r0 != NONE && !cs_submap.node(r0).deleted) {
            bool r0_is_inner = false;
            for (std::size_t ai : cs_submap.node(r0).arcs) {
                const auto& a = cs_submap.arc(ai);
                if (a.first_edge != NONE) {
                    std::size_t rep = std::min(a.first_edge, a.last_edge);
                    if (rep >= c_lo && rep < c_hi) {
                        r0_is_inner = true;
                        break;
                    }
                }
            }
            return r0_is_inner ? -1 : 1; // descend into inner side
        }
        return -1; // default: try left
    } else {
        // A₂ faces the "outer" interval.
        // The "inner" interval is shielded → descend into outer.
        std::size_t r0 = centroid_chord.region[0];
        if (r0 != NONE && !cs_submap.node(r0).deleted) {
            bool r0_is_inner = false;
            for (std::size_t ai : cs_submap.node(r0).arcs) {
                const auto& a = cs_submap.arc(ai);
                if (a.first_edge != NONE) {
                    std::size_t rep = std::min(a.first_edge, a.last_edge);
                    if (rep >= c_lo && rep < c_hi) {
                        r0_is_inner = true;
                        break;
                    }
                }
            }
            return r0_is_inner ? 1 : -1; // descend into outer side
        }
        return 1; // default: try right
    }
}

/// Binary search through tree decomposition T of subarc α's canonical
/// submap Sα to find a vertex of α that sees A₂.  §3.2, Lemma 3.2.
///
/// At each internal node with centroid chord ab:
///   1. Test a and b via local shooting (Lemma 2.4).
///   2. If neither sees A₂, identify the shielded side and descend
///      into the non-shielded subtree.
///
/// At a leaf: exhaustively check all vertices in the corresponding
/// region that belong to α.
ChordCandidate search_tree_decomposition(
    const Submap& parent_submap,
    std::size_t region_idx,
    const Polygon& polygon,
    const CanonicalSubmap& cs,
    std::size_t alpha_start,
    std::size_t alpha_end,
    const ArcStructure& arc_j) {

    const auto& cs_submap = cs.submap;

    // Use cached tree decomposition from the canonical submap.
    // Build it on first use if not yet built.
    if (!cs.td_built) {
        const_cast<CanonicalSubmap&>(cs).tree_decomp.build(cs_submap);
        const_cast<CanonicalSubmap&>(cs).td_built = true;
    }
    const auto& td = cs.tree_decomp;

    if (td.size() == 0) return {};

    // Walk from root down the tree.
    std::size_t cur = td.root();

    while (cur != TD_NONE) {
        const auto& td_node = td.node(cur);

        if (td_node.is_leaf()) {
            // Exhaustive check: test each vertex in this leaf's region
            // that belongs to α = [alpha_start, alpha_end).
            std::size_t leaf_region = td_node.region_idx;
            if (leaf_region == TD_NONE || leaf_region >= cs_submap.num_nodes())
                break;

            const auto& leaf_nd = cs_submap.node(leaf_region);
            // Collect vertex range from this region's arcs.
            for (std::size_t ai : leaf_nd.arcs) {
                const auto& a = cs_submap.arc(ai);
                if (a.first_edge == NONE) continue;
                std::size_t lo = std::min(a.first_edge, a.last_edge);
                std::size_t hi = std::max(a.first_edge, a.last_edge);
                // Each edge ei corresponds to vertex ei and ei+1.
                for (std::size_t v = lo; v <= hi + 1; ++v) {
                    if (v < alpha_start || v >= alpha_end) continue;
                    if (v >= polygon.num_vertices()) continue;

                    const auto& pt = polygon.vertex(v);
                    std::size_t ev = (v > 0) ? v - 1 : 0;
                    double origin_x = pt.x;
                    if (ev < polygon.num_edges()) {
                        const auto& edge = polygon.edge(ev);
                        const auto& p1 = polygon.vertex(edge.start_idx);
                        const auto& p2 = polygon.vertex(edge.end_idx);
                        if (std::abs(p2.y - p1.y) > 1e-15) {
                            double t = (pt.y - p1.y) / (p2.y - p1.y);
                            origin_x = p1.x + t * (p2.x - p1.x);
                        }
                    }
                    // Shoot in both directions.
                    for (bool dir : {true, false}) {
                        auto hit = shoot_in_region(
                            parent_submap, region_idx, polygon,
                            origin_x, pt.y, dir);
                        if (hit_lands_on_arc(hit, parent_submap, arc_j)) {
                            ChordCandidate cand;
                            cand.y = pt.y;
                            cand.left_edge = v;
                            // Resolve right_edge: find the edge on
                            // the hit arc whose x-intercept is closest to hit_x.
                            cand.right_edge = NONE;
                            if (hit.arc_idx != NONE &&
                                hit.arc_idx < parent_submap.num_arcs()) {
                                const auto& ha = parent_submap.arc(hit.arc_idx);
                                if (ha.first_edge != NONE) {
                                    std::size_t hlo = std::min(ha.first_edge, ha.last_edge);
                                    std::size_t hhi = std::max(ha.first_edge, ha.last_edge);
                                    double best_xd = std::numeric_limits<double>::infinity();
                                    for (std::size_t ei2 = hlo; ei2 <= hhi && ei2 < polygon.num_edges(); ++ei2) {
                                        const auto& e2 = polygon.edge(ei2);
                                        const auto& ep1 = polygon.vertex(e2.start_idx);
                                        const auto& ep2 = polygon.vertex(e2.end_idx);
                                        double ey_lo = std::min(ep1.y, ep2.y);
                                        double ey_hi = std::max(ep1.y, ep2.y);
                                        if (pt.y < ey_lo - 1e-12 || pt.y > ey_hi + 1e-12) continue;
                                        if (std::abs(ey_hi - ey_lo) < 1e-15) continue;
                                        double t2 = (pt.y - ep1.y) / (ep2.y - ep1.y);
                                        double ex = ep1.x + t2 * (ep2.x - ep1.x);
                                        double dx = std::abs(ex - hit.hit_x);
                                        if (dx < best_xd) {
                                            best_xd = dx;
                                            cand.right_edge = ei2;
                                        }
                                    }
                                    // No vertex fallback needed; right_edge is set to the hit edge.
                                }
                            }
                            // Skip if we couldn't resolve or got a self-loop.
                            if (cand.right_edge != NONE && cand.right_edge != v) {
                                cand.valid = true;
                                return cand;
                            }
                        }
                    }
                }
            }
            break; // done with this leaf
        }

        // Internal node: apply the Lemma 2.4 test.
        assert(td_node.chord_idx != TD_NONE);
        if (td_node.chord_idx >= cs_submap.num_chords()) break;

        const auto& centroid = cs_submap.chord(td_node.chord_idx);
        ChordCandidate cand;
        int result = lemma24_test(
            parent_submap, region_idx, polygon,
            cs_submap, centroid, arc_j,
            alpha_start, alpha_end, cand);

        if (result == 0) {
            return cand; // success
        } else if (result < 0) {
            // Descend into left child.
            cur = td_node.left_child;
        } else if (result == 1) {
            // Descend into right child.
            cur = td_node.right_child;
        } else {
            // Indeterminate — the topological test could not determine
            // the shielded side (degenerate geometry).  Per the paper,
            // Lemma 2.4 should always determine the shielded side for
            // non-degenerate inputs.  As a fallback, try both children
            // but descend only one level deep into each, checking only
            // the nearest leaf.  This bounds the work to O(log r) per
            // indeterminate case instead of O(r).
            //
            // Try left child first — descend to its leftmost leaf.
            std::size_t try_left = td_node.left_child;
            std::size_t try_right = td_node.right_child;

            // Quick probe of the left subtree's nearest leaf.
            std::size_t probe = try_left;
            while (probe != TD_NONE && !td.node(probe).is_leaf()) {
                if (td.node(probe).left_child != TD_NONE)
                    probe = td.node(probe).left_child;
                else
                    probe = td.node(probe).right_child;
            }
            if (probe != TD_NONE && td.node(probe).is_leaf()) {
                std::size_t leaf_region = td.node(probe).region_idx;
                if (leaf_region != TD_NONE && leaf_region < cs_submap.num_nodes()) {
                    const auto& leaf_nd = cs_submap.node(leaf_region);
                    bool found = false;
                    for (std::size_t ai2 : leaf_nd.arcs) {
                        if (found) break;
                        const auto& a2 = cs_submap.arc(ai2);
                        if (a2.first_edge == NONE) continue;
                        std::size_t lo2 = std::min(a2.first_edge, a2.last_edge);
                        std::size_t hi2 = std::max(a2.first_edge, a2.last_edge);
                        for (std::size_t v = lo2; v <= hi2 + 1 && !found; ++v) {
                            if (v < alpha_start || v >= alpha_end) continue;
                            if (v >= polygon.num_vertices()) continue;
                            const auto& pt2 = polygon.vertex(v);
                            std::size_t ev2 = (v > 0) ? v - 1 : 0;
                            double ox2 = pt2.x;
                            if (ev2 < polygon.num_edges()) {
                                const auto& edge2 = polygon.edge(ev2);
                                const auto& pp1 = polygon.vertex(edge2.start_idx);
                                const auto& pp2 = polygon.vertex(edge2.end_idx);
                                if (std::abs(pp2.y - pp1.y) > 1e-15) {
                                    double tt = (pt2.y - pp1.y) / (pp2.y - pp1.y);
                                    ox2 = pp1.x + tt * (pp2.x - pp1.x);
                                }
                            }
                            for (bool dir2 : {true, false}) {
                                auto hit2 = shoot_in_region(
                                    parent_submap, region_idx, polygon,
                                    ox2, pt2.y, dir2);
                                if (hit_lands_on_arc(hit2, parent_submap, arc_j)) {
                                    ChordCandidate cand2;
                                    cand2.y = pt2.y;
                                    cand2.left_edge = v;
                                    cand2.right_edge = NONE;
                                    if (hit2.arc_idx != NONE &&
                                        hit2.arc_idx < parent_submap.num_arcs()) {
                                        const auto& ha2 = parent_submap.arc(hit2.arc_idx);
                                        if (ha2.first_edge != NONE) {
                                            std::size_t hlo2 = std::min(ha2.first_edge, ha2.last_edge);
                                            std::size_t hhi2 = std::max(ha2.first_edge, ha2.last_edge);
                                            double best_xd2 = std::numeric_limits<double>::infinity();
                                            for (std::size_t ei3 = hlo2; ei3 <= hhi2 && ei3 < polygon.num_edges(); ++ei3) {
                                                const auto& e3 = polygon.edge(ei3);
                                                const auto& ep1 = polygon.vertex(e3.start_idx);
                                                const auto& ep2 = polygon.vertex(e3.end_idx);
                                                double ey_lo = std::min(ep1.y, ep2.y);
                                                double ey_hi = std::max(ep1.y, ep2.y);
                                                if (pt2.y < ey_lo - 1e-12 || pt2.y > ey_hi + 1e-12) continue;
                                                if (std::abs(ey_hi - ey_lo) < 1e-15) continue;
                                                double t3 = (pt2.y - ep1.y) / (ep2.y - ep1.y);
                                                double ex = ep1.x + t3 * (ep2.x - ep1.x);
                                                double dx = std::abs(ex - hit2.hit_x);
                                                if (dx < best_xd2) {
                                                    best_xd2 = dx;
                                                    cand2.right_edge = ei3;
                                                }
                                            }
                                            // No vertex fallback needed; right_edge is set to the hit edge.
                                        }
                                    }
                                    if (cand2.right_edge != NONE && cand2.right_edge != v) {
                                        cand2.valid = true;
                                        return cand2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // Left probe failed — try the right subtree.
            cur = try_right;
        }
    }

    return {};
}

/// Find a chord to insert that splits a high-degree region.
///
/// For a region R with > 4 arcs, Lemma 3.3 guarantees two non-consecutive
/// arcs Aᵢ, Aⱼ with a vertex on Aᵢ visible to Aⱼ.
///
/// §3.2 algorithm:
///   1. Pick two non-consecutive arcs A₁, A₂ of R.
///   2. Decompose A₁ into O(log γ) subarcs via arc-cutting.
///   3. For each subarc α, retrieve its canonical submap Sα.
///   4. Binary search through the tree decomposition of Sα (Lemma 3.2):
///      at each level, test the centroid chord via local shooting and
///      Lemma 2.4 to narrow down.
///   5. Find a vertex v on A₁ that sees A₂; return chord candidate.
ChordCandidate find_splitting_chord(
    const Submap& submap,
    std::size_t region_idx,
    const GradeStorage& storage,
    const Polygon& polygon,
    std::size_t granularity) {

    const auto& nd = submap.node(region_idx);
    assert(nd.degree() > 4);

    if (nd.arcs.size() < 2) return {};

    // Sort the region's arcs by edge index to get the circular
    // boundary order along ∂C.  This is essential for correctly
    // determining which arcs are "consecutive" (adjacent in the
    // boundary sequence) vs "non-consecutive" per Lemma 3.3.
    std::vector<std::size_t> sorted_arcs;
    sorted_arcs.reserve(nd.arcs.size());
    for (std::size_t a : nd.arcs) {
        if (submap.arc(a).first_edge != NONE)
            sorted_arcs.push_back(a);
    }
    if (sorted_arcs.size() < 2) return {};

    std::sort(sorted_arcs.begin(), sorted_arcs.end(),
              [&](std::size_t a1, std::size_t a2) {
                  return std::min(submap.arc(a1).first_edge,
                                  submap.arc(a1).last_edge)
                       < std::min(submap.arc(a2).first_edge,
                                  submap.arc(a2).last_edge);
              });

    // §3.2: "apply Lemma 3.2 to every pair of nonconsecutive arcs
    // until we find a chord which we can add to S."
    //
    // Two arcs are consecutive if they are adjacent in the sorted
    // (boundary-order) sequence.  We iterate over ALL non-consecutive
    // pairs, trying each until we find a splitting chord.
    std::size_t k = sorted_arcs.size();

    auto index_distance = [&](std::size_t idx_a, std::size_t idx_b) -> std::size_t {
        std::size_t d = (idx_b >= idx_a) ? (idx_b - idx_a) : (idx_a - idx_b);
        return std::min(d, k - d);
    };

    // Collect all non-consecutive pairs.  Try (0, k/2) first as a
    // likely candidate, then all remaining non-consecutive pairs.
    struct ArcPair { std::size_t ii; std::size_t jj; };
    std::vector<ArcPair> pairs;

    // Preferred pair first (maximally separated).
    if (index_distance(0, k / 2) > 1) {
        pairs.push_back({0, k / 2});
    }

    // All other non-consecutive pairs.
    for (std::size_t ii = 0; ii < k; ++ii) {
        for (std::size_t jj = ii + 2; jj < k; ++jj) {
            if (index_distance(ii, jj) <= 1) continue;
            // Skip the preferred pair already added above.
            if (ii == 0 && jj == k / 2) continue;
            pairs.push_back({ii, jj});
        }
    }

    for (const auto& [pi, pj] : pairs) {
        std::size_t ai = sorted_arcs[pi];
        std::size_t aj = sorted_arcs[pj];

        const auto& arc_i_final = submap.arc(ai);
        const auto& arc_j_final = submap.arc(aj);
        if (arc_i_final.first_edge == NONE || arc_j_final.first_edge == NONE)
            continue;

        // Decompose A₁ into O(log γ) dyadic pieces.
        std::size_t start = std::min(arc_i_final.first_edge, arc_i_final.last_edge);
        std::size_t end   = std::max(arc_i_final.first_edge, arc_i_final.last_edge) + 1;
        auto pieces = cut_arc(start, end, polygon);

        // §3.2: For each piece α, retrieve its canonical submap and binary
        // search through the tree decomposition.
        for (auto& piece : pieces) {
            if (piece.grade >= storage.num_grades()) continue;
            if (piece.chain_index >= storage.num_chains(piece.grade)) continue;

            const auto& cs = storage.get(piece.grade, piece.chain_index);
            if (!cs.oracle.is_built()) continue;

            // Binary search through tree decomposition of Sα.
            auto cand = search_tree_decomposition(
                submap, region_idx, polygon, cs,
                piece.start_vertex, piece.end_vertex,
                arc_j_final);
            if (cand.valid) return cand;
        }

        // Also try with the roles swapped: search A₂ for a vertex
        // that sees A₁.  The paper says "apply Lemma 3.2 to every
        // pair," and the lemma is stated with A₁ containing the
        // vertex — so we must try both directions.
        std::size_t start2 = std::min(arc_j_final.first_edge, arc_j_final.last_edge);
        std::size_t end2   = std::max(arc_j_final.first_edge, arc_j_final.last_edge) + 1;
        auto pieces2 = cut_arc(start2, end2, polygon);

        for (auto& piece : pieces2) {
            if (piece.grade >= storage.num_grades()) continue;
            if (piece.chain_index >= storage.num_chains(piece.grade)) continue;

            const auto& cs = storage.get(piece.grade, piece.chain_index);
            if (!cs.oracle.is_built()) continue;

            auto cand = search_tree_decomposition(
                submap, region_idx, polygon, cs,
                piece.start_vertex, piece.end_vertex,
                arc_i_final);
            if (cand.valid) return cand;
        }
    }

    (void)granularity;
    return {};
}

} // anonymous namespace

void restore_conformality(Submap& submap,
                          const GradeStorage& storage,
                          const Polygon& polygon,
                          std::size_t granularity) {
    // Worklist-based approach: only (re-)examine regions with degree > 4.
    // When a chord is inserted, only the two affected regions need to be
    // re-checked, giving O(total insertions) work instead of O(n²).
    // §2.3: "conformal submaps [are] those with node-degree at most 4."
    // Node-degree = number of incident chords (tree edges).
    std::queue<std::size_t> worklist;
    for (std::size_t i = 0; i < submap.num_nodes(); ++i) {
        if (!submap.node(i).deleted && submap.node(i).degree() > 4) {
            worklist.push(i);
        }
    }

    // Safety bound: in a conformal submap with N arcs, we need at most
    // N chord insertions to reduce every region to ≤ 4 arcs.  If we
    // exceed this, something is wrong (degenerate chord endpoints
    // causing splits that don't reduce degree).
    std::size_t total_arcs = 0;
    for (std::size_t i = 0; i < submap.num_nodes(); ++i) {
        if (!submap.node(i).deleted)
            total_arcs += submap.node(i).arcs.size();
    }
    std::size_t max_iterations = total_arcs + submap.num_chords() + 100;
    std::size_t iterations = 0;

    while (!worklist.empty()) {
        std::size_t i = worklist.front();
        worklist.pop();

        if (++iterations > max_iterations) break; // safety bound

        if (submap.node(i).deleted) continue;
        if (submap.node(i).degree() <= 4) continue;

            // Find a chord to split this region.
            auto cand = find_splitting_chord(submap, i, storage,
                                              polygon, granularity);
            if (cand.valid) {
                // Skip degenerate self-loop chords (left == right).
                // These arise when vertex resolution fails to find a
                // distinct target on the opposite arc.
                if (cand.right_edge == NONE ||
                    cand.left_edge == cand.right_edge)
                    continue;

                // Split region i by inserting a new chord.
                // The new chord creates a new region.
                std::size_t new_region = submap.add_node();

                Chord new_chord;
                new_chord.y = cand.y;
                new_chord.left_edge = cand.left_edge;
                new_chord.right_edge = cand.right_edge;
                new_chord.region[0] = i;
                new_chord.region[1] = new_region;

                std::size_t new_ci = submap.add_chord(new_chord);
                std::fprintf(stderr, "    CONF add chord %zu: le=%zu re=%zu r0=%zu r1=%zu\n",
                             new_ci, new_chord.left_edge, new_chord.right_edge, i, new_region);

                // Partition arcs between region i and new_region.
                //
                // Per Lemma 2.1, removing the two visible points (the chord
                // endpoints) splits ∂C into two boundary pieces.  Each arc
                // belongs to exactly one piece, determined by whether the
                // arc's edge range falls on the same side of the chord's
                // vertex as the first boundary piece or the second.
                //
                // Topological classification: arcs whose first_edge index
                // lies between left_edge and right_edge (in the
                // circular edge ordering) go to one region; the rest go to
                // the other.
                std::size_t lv = cand.left_edge;
                std::size_t rv = (cand.right_edge != NONE)
                                    ? cand.right_edge : lv;

                // Normalise so lo ≤ hi in the edge index space.
                std::size_t lo = std::min(lv, rv);
                std::size_t hi = std::max(lv, rv);

                // §2.3: If a chord endpoint falls in the middle of an
                // arc, split that arc into two sub-arcs at the endpoint.
                // We must split before partitioning so each sub-arc
                // lands entirely on one side of the chord.
                {
                    // We iterate a snapshot of arcs since split_arc_at_vertex
                    // may add new arcs to the region.
                    auto arcs_snapshot = submap.node(i).arcs;
                    for (std::size_t arc_idx : arcs_snapshot) {
                        submap.split_arc_at_vertex(arc_idx, lo);
                        submap.split_arc_at_vertex(arc_idx, hi);
                    }
                    // Also split any newly created arcs at the other endpoint.
                    auto arcs_snapshot2 = submap.node(i).arcs;
                    for (std::size_t arc_idx : arcs_snapshot2) {
                        submap.split_arc_at_vertex(arc_idx, lo);
                        submap.split_arc_at_vertex(arc_idx, hi);
                    }
                }

                auto& arcs = submap.node(i).arcs;
                std::vector<std::size_t> keep_arcs;   // stay with region i
                std::vector<std::size_t> move_arcs;   // go to new_region

                for (std::size_t arc_idx : arcs) {
                    auto& a = submap.arc(arc_idx);
                    if (a.first_edge == NONE) {
                        // Null-length arc — keep with the original region.
                        keep_arcs.push_back(arc_idx);
                        continue;
                    }
                    // Exact topological test: the chord endpoints
                    // [lo, hi) partition ∂C into two boundary pieces.
                    // An arc belongs to the inner piece if its edge
                    // range overlaps [lo, hi).  Since chords split ∂C
                    // cleanly, each arc is entirely on one side.
                    std::size_t a_lo = std::min(a.first_edge, a.last_edge);
                    std::size_t a_hi = std::max(a.first_edge, a.last_edge);
                    if (a_lo < hi && a_hi >= lo) {
                        move_arcs.push_back(arc_idx);
                    } else {
                        keep_arcs.push_back(arc_idx);
                    }
                }

                // If the split is degenerate (everything on one side),
                // move the arc with the largest edge range to the other
                // side.  This ensures a meaningful split.
                if (move_arcs.empty() && keep_arcs.size() > 1) {
                    // Find the arc with the largest edge range to move.
                    std::size_t best_idx = 0;
                    std::size_t best_range = 0;
                    for (std::size_t k = 0; k < keep_arcs.size(); ++k) {
                        const auto& ka = submap.arc(keep_arcs[k]);
                        std::size_t range = 0;
                        if (ka.first_edge != NONE) {
                            range = std::max(ka.first_edge, ka.last_edge)
                                  - std::min(ka.first_edge, ka.last_edge) + 1;
                        }
                        if (range > best_range) {
                            best_range = range;
                            best_idx = k;
                        }
                    }
                    move_arcs.push_back(keep_arcs[best_idx]);
                    keep_arcs.erase(keep_arcs.begin() +
                                    static_cast<std::ptrdiff_t>(best_idx));
                }

                // Apply the partition.
                arcs = keep_arcs;
                for (std::size_t arc_idx : move_arcs) {
                    submap.arc(arc_idx).region_node = new_region;
                    submap.node(new_region).arcs.push_back(arc_idx);
                }

                // Partition incident chords: a chord belongs to the
                // new region if EITHER of its vertex endpoints falls
                // in the inner boundary piece [lo, hi).
                auto& ichords = submap.node(i).incident_chords;
                std::vector<std::size_t> keep_chords;
                std::vector<std::size_t> move_chords;

                for (std::size_t ci : ichords) {
                    if (ci == new_ci) {
                        // The newly-inserted chord must not be
                        // re-classified — its regions are already
                        // set correctly.
                        keep_chords.push_back(ci);
                        continue;
                    }
                    auto& c = submap.chord(ci);
                    // Classify by both endpoints' positions.
                    std::size_t lv = c.left_edge;
                    std::size_t rv = c.right_edge;
                    bool lv_inner = (lv != NONE && lv >= lo && lv < hi);
                    bool rv_inner = (rv != NONE && rv >= lo && rv < hi);
                    if (lv_inner || rv_inner) {
                        move_chords.push_back(ci);
                    } else {
                        keep_chords.push_back(ci);
                    }
                }

                ichords = keep_chords;
                for (std::size_t ci : move_chords) {
                    for (int s = 0; s < 2; ++s) {
                        if (submap.chord(ci).region[s] == i) {
                            submap.chord(ci).region[s] = new_region;
                        }
                    }
                    submap.node(new_region).incident_chords.push_back(ci);
                }

                submap.recompute_weight(i);
                submap.recompute_weight(new_region);

                // Re-examine both regions if they still have high degree.
                if (submap.node(i).degree() > 4) worklist.push(i);
                if (submap.node(new_region).degree() > 4) worklist.push(new_region);
            }
    }
}

} // namespace chazelle
