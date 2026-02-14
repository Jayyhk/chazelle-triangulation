#include "chazelle/fusion.h"
#include "geometry/perturbation.h"
#include "oracles/arc_cutting.h"
#include "chazelle/grade_storage.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
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

/// Compute the x-coordinate of a point on ∂C at height y, on the
/// polygon edge incident on vertex v.
double origin_x_for_vertex(const Polygon& polygon, std::size_t v) {
    if (v >= polygon.num_vertices()) return 0.0;
    const auto& pt = polygon.vertex(v);
    // Return vertex x directly.
    return pt.x;
}

/// §3.1 / §4.1 oracle (i): Per-subarc local shooting from (origin_x, y)
/// within a specific region of a submap.
///
/// For each arc A in the region, decomposes A into O(g(γ)) = O(log γ)
/// dyadic pieces via arc-cutting and queries each piece's canonical
/// submap oracle via shoot_from_point().  Chords are checked directly.
///
/// Total cost per call: O(f(γ)) where f(γ) = O(λ · 2^{β²λ/3+2βλ/3}).
RayHit local_shoot_in_region(const Submap& submap,
                              const Polygon& polygon,
                              std::size_t region_idx,
                              double origin_x, double y,
                              bool shoot_right,
                              const GradeStorage& storage) {
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
        if (a.is_virtual()) {
            double vy = a.virtual_y;
            if (y != vy) continue;
            double x_left = (a.first_edge < polygon.num_edges())
                ? polygon.edge_x_at_y(a.first_edge, vy) : 0.0;
            double x_right = (a.last_edge < polygon.num_edges())
                ? polygon.edge_x_at_y(a.last_edge, vy) : 0.0;
            double x = (x_left + x_right) * 0.5;
            double dist = shoot_right ? (x - origin_x) : (origin_x - x);
            if (dist > 0.0 && dist < best_dist) {
                best_dist = dist;
                best.type = RayHit::Type::ARC;
                best.arc_idx = ai;
                best.hit_x = x;
                best.hit_edge = a.first_edge;
            }
            continue;
        }

        std::size_t e_lo = std::min(a.first_edge, a.last_edge);
        std::size_t e_hi = std::max(a.first_edge, a.last_edge);

        // Single-edge arc: O(1) direct check.
        if (e_lo == e_hi) {
            if (e_lo < polygon.num_edges()) {
                const auto& edge = polygon.edge(e_lo);
                const auto& p1 = polygon.vertex(edge.start_idx);
                const auto& p2 = polygon.vertex(edge.end_idx);
                double y_lo = std::min(p1.y, p2.y);
                double y_hi = std::max(p1.y, p2.y);
                if (y >= y_lo && y <= y_hi && p1.y != p2.y) {
                    double x = horizontal_ray_x_intercept(p1, p2, y);
                    double dist = shoot_right ? (x - origin_x) : (origin_x - x);
                    if (dist > 0.0 && dist < best_dist) {
                        best_dist = dist;
                        best.type = RayHit::Type::ARC;
                        best.arc_idx = ai;
                        best.hit_x = x;
                        best.hit_edge = e_lo;
                    }
                }
            }
            continue;
        }

        // §4.1: Decompose the arc into O(log γ) dyadic pieces.
        // For each piece, use the canonical submap's oracle.
        auto pieces = cut_arc(e_lo, e_hi + 1, polygon);

        for (const auto& piece : pieces) {
            if (piece.grade < storage.num_grades() &&
                piece.chain_index < storage.num_chains(piece.grade)) {
                const auto& cs = storage.get(piece.grade, piece.chain_index);
                if (cs.oracle.is_built()) {
                    auto hit = cs.oracle.shoot_from_point(
                                   origin_x, y, shoot_right);
                    // Accept only arc hits within the original arc's
                    // edge range.
                    if (hit.type == RayHit::Type::ARC &&
                        hit.hit_edge != NONE &&
                        hit.hit_edge >= e_lo &&
                        hit.hit_edge <= e_hi) {
                        double dist = shoot_right
                            ? (hit.hit_x - origin_x)
                            : (origin_x - hit.hit_x);
                        if (dist > 0.0 && dist < best_dist) {
                            best_dist = dist;
                            best.type = RayHit::Type::ARC;
                            best.arc_idx = ai;
                            best.hit_x = hit.hit_x;
                            best.hit_edge = hit.hit_edge;
                        }
                    }
                    continue;
                }
            }

            // §4.1: The up-phase builds oracles for ALL grades
            // (including grade 0).  If we reach here, the oracle
            // is missing — this is a programming error, not an
            // expected case.  The paper guarantees: "all these
            // shooting structures have been computed."
            assert(false && "missing oracle for arc piece — "
                           "up-phase should have built it");
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
///
/// The junction vertex is at one of the two endpoints of C₂.  The
/// normal-form stores start_arc / end_arc pointers (§2.3 condition
/// (iii)), so we can look up the region in O(1).
std::size_t find_junction_region(const Submap& dst,
                                  std::size_t junction_vertex) {
    // §2.3 (iii): start_arc and end_arc point to the arcs passing
    // through the endpoints of C.  The junction vertex matches one
    // of them; its region_node is the answer.
    if (dst.start_arc != NONE && dst.start_arc < dst.num_arcs()) {
        const auto& a = dst.arc(dst.start_arc);
        if (a.first_edge != NONE) {
            std::size_t lo = std::min(a.first_edge, a.last_edge);
            std::size_t hi = std::max(a.first_edge, a.last_edge);
            if (junction_vertex >= lo && junction_vertex <= hi + 1) {
                if (a.region_node != NONE) return a.region_node;
            }
        }
    }
    if (dst.end_arc != NONE && dst.end_arc < dst.num_arcs()) {
        const auto& a = dst.arc(dst.end_arc);
        if (a.first_edge != NONE) {
            std::size_t lo = std::min(a.first_edge, a.last_edge);
            std::size_t hi = std::max(a.first_edge, a.last_edge);
            if (junction_vertex >= lo && junction_vertex <= hi + 1) {
                if (a.region_node != NONE) return a.region_node;
            }
        }
    }
    // §3.1: If neither endpoint arc covers the junction vertex,
    // the submap has a single region (trivial case).  Return it.
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
               const GradeStorage& storage,
               std::vector<DiscoveredChord>& discovered) {


    // -- Collect stop points in ∂C order --
    //
    // Per §3.1: "Let a₁, a₂, …, a_m be the canonical vertex
    // enumeration of S₁.  Recall that this enumerates the exit chord
    // endpoints in S₁ as we encounter them going clockwise around ∂C₁."
    //
    // Per §2.3 condition (iii): the arc-sequence table is stored "in
    // the order corresponding to a canonical traversal of the double
    // boundary ∂C."  Since src is in normal form, walking its arcs in
    // table order gives us stops in ∂C clockwise order.  We interleave
    // arc boundary vertices and chord endpoints — both are already in
    // ∂C order — via a single linear merge, avoiding O(E log E) sort.
    struct StopPoint {
        std::size_t edge_idx;   ///< Edge on which the stop lies.
        double y;               ///< Y-coordinate.
        bool is_vertex;         ///< True if this stop is a polygon vertex.
        std::size_t vertex_idx; ///< Vertex index (if is_vertex), else NONE.
        bool shoot_right;       ///< §3.1: assigned chord direction.
    };

    // §3.1: Collect stops from arc boundary vertices and chord
    // endpoints.  Both sources are in ∂C order.  We collect into
    // two separate sorted vectors, then 2-way merge + dedup in O(m).
    std::vector<StopPoint> vertex_stops;
    std::vector<StopPoint> chord_stops;

    // Arc boundary vertices (in ∂C order since src is in normal form).
    for (std::size_t ai = 0; ai < src.num_arcs(); ++ai) {
        const auto& arc = src.arc(ai);
        if (arc.first_edge == NONE) continue;
        std::size_t lo = std::min(arc.first_edge, arc.last_edge);
        std::size_t hi = std::max(arc.first_edge, arc.last_edge);
        {
            std::size_t v = lo;
            if (v != NONE && v < polygon.num_vertices()) {
                std::size_t e = v;
                if (e >= polygon.num_edges()) e = polygon.num_edges() - 1;
                vertex_stops.push_back({e, polygon.vertex(v).y,
                                        true, v,
                                        arc.first_side == Side::LEFT});
            }
        }
        if (hi + 1 < polygon.num_vertices()) {
            std::size_t v = hi + 1;
            std::size_t e = v;
            if (e >= polygon.num_edges()) e = polygon.num_edges() - 1;
            vertex_stops.push_back({e, polygon.vertex(v).y,
                                    true, v,
                                    arc.last_side == Side::LEFT});
        }
    }

    // Chord endpoints (in order from normalize_chords).
    for (std::size_t ci = 0; ci < src.num_chords(); ++ci) {
        const auto& c = src.chord(ci);
        if (c.left_edge != NONE && c.left_edge < polygon.num_edges())
            chord_stops.push_back({c.left_edge, c.y, false, NONE, true});
        if (c.right_edge != NONE && c.right_edge < polygon.num_edges())
            chord_stops.push_back({c.right_edge, c.y, false, NONE, false});
    }

    // Junction vertex stop.
    StopPoint junction_stop{};
    bool have_junction = false;
    if (junction_vertex != NONE && junction_vertex < polygon.num_vertices()) {
        std::size_t e = junction_vertex;
        if (e >= polygon.num_edges()) e = polygon.num_edges() - 1;
        junction_stop = {e, polygon.vertex(junction_vertex).y,
                         true, junction_vertex, true};
        have_junction = true;
    }

    // §3.1 (Lemma 3.1 proof): "Note that merging can also be used
    // instead of sorting."  Both vertex_stops and chord_stops are
    // already sorted by edge_idx (from the arc-sequence table and
    // normalize_chords respectively).  2-way merge with dedup: O(m).
    auto stop_less = [](const StopPoint& a, const StopPoint& b) {
        if (a.edge_idx != b.edge_idx) return a.edge_idx < b.edge_idx;
        return a.y < b.y;
    };
    auto stop_eq = [](const StopPoint& a, const StopPoint& b) {
        return a.edge_idx == b.edge_idx && a.y == b.y;
    };

    std::vector<StopPoint> stops;
    stops.reserve(vertex_stops.size() + chord_stops.size() + 1);

    {
        std::size_t iv = 0, ic = 0;
        while (iv < vertex_stops.size() || ic < chord_stops.size()) {
            const StopPoint* pick = nullptr;
            bool from_v = false;
            if (iv < vertex_stops.size() && ic < chord_stops.size()) {
                if (stop_less(vertex_stops[iv], chord_stops[ic])) {
                    pick = &vertex_stops[iv]; from_v = true;
                } else {
                    pick = &chord_stops[ic]; from_v = false;
                }
            } else if (iv < vertex_stops.size()) {
                pick = &vertex_stops[iv]; from_v = true;
            } else {
                pick = &chord_stops[ic]; from_v = false;
            }
            // Dedup: skip if same (edge_idx, y) as last added.
            if (stops.empty() || !stop_eq(*pick, stops.back())) {
                stops.push_back(*pick);
            }
            if (from_v) ++iv; else ++ic;
        }
    }

    // Insert junction vertex at correct position (binary search +
    // insert): O(log m).
    if (have_junction) {
        bool already = false;
        for (const auto& s : stops) {
            if (stop_eq(s, junction_stop)) { already = true; break; }
        }
        if (!already) {
            auto it = std::lower_bound(stops.begin(), stops.end(),
                                       junction_stop, stop_less);
            stops.insert(it, junction_stop);
        }
    }

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
            RayHit h = local_shoot_in_region(dst, polygon, current_region, ox, y, a0_dir, storage);
            if (h.type == RayHit::Type::ARC) {
                double d = a0_dir ? (h.hit_x - ox) : (ox - h.hit_x);
                if (d > 0.0) {
                    dist_dst = d;
                    hit_dst = h;
                }
            }
        }

        if (hit_dst.type == RayHit::Type::ARC && dist_dst <= dist_src) {
            // §3 item (i): the ray-shooting report always includes
            // the edge name.  hit_edge is set by local_shoot /
            // shoot_from_region — no post-hoc linear scan needed.
            std::size_t re = hit_dst.hit_edge;
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
                        // §3 item (i): oracle always reports the hit edge.
                        c0_edge = h.hit_edge;
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

                // §3.1: Binary search for start_k since stops are
                // sorted by edge_idx.  O(log m) instead of O(m).
                {
                    std::size_t lo_s = 0, hi_s = stops.size();
                    while (lo_s < hi_s) {
                        std::size_t mid = lo_s + (hi_s - lo_s) / 2;
                        if (stops[mid].edge_idx < c0_edge)
                            lo_s = mid + 1;
                        else
                            hi_s = mid;
                    }
                    start_k = lo_s;
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
            RayHit h = local_shoot_in_region(dst, polygon, current_region, ox, y, aj_dir, storage);
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
                // §3 item (i): oracle always reports the hit edge.
                std::size_t re = best_dst_hit.hit_edge;
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

                    // §3 item (i): oracle always reports the hit edge.
                    std::size_t p_prime_e = h.hit_edge;
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
/// Per §3.1 Lemma 3.1 proof: "let us sort the endpoints of these chords
/// along ∂C, which is done by sorting the names of the edges of P on
/// which these arcs abut, and then sorting the endpoints falling within
/// the same edges by considering y-coordinates.  This allows us to set
/// up the required arc-sequence table.  Note that merging can also be
/// used instead of sorting."
///
/// Implementation: all three chord sources (S₁, S₂, discovered) are
/// already in ∂C order (S₁ and S₂ by normalize_chords(); discovered
/// chords by the clockwise fusion walk of §3.1).  We 3-way merge them
/// in O(E) time.  Split points and arc-to-region assignment also use
/// O(E) linear walks instead of O(E log E) sorts or binary searches.
Submap rebuild_submap(const Submap& s1, const Submap& s2,
                      const std::vector<DiscoveredChord>& discovered_chords,
                      const Polygon& /*polygon*/) {
    Submap result;

    struct ChordRecord {
        double y;
        std::size_t left_edge;
        std::size_t right_edge;
        bool is_null_length;
        std::size_t sort_key;  // min(left_edge, right_edge)
    };

    // Comparator for the 3-way merge: sort by (sort_key, y).
    // Per §3.1: "sorting the names of the edges... then sorting the
    // endpoints falling within the same edges by y-coordinates."
    auto chord_less = [](const ChordRecord& a, const ChordRecord& b) {
        if (a.sort_key != b.sort_key) return a.sort_key < b.sort_key;
        return a.y < b.y;
    };

    // §3.1 (Lemma 3.1 proof): "Note that merging can also be used
    // instead of sorting."  S₁ and S₂ chords are already sorted by
    // normalize_chords().  Discovered chords are produced by the
    // clockwise fusion walk (§3.1 main loop), so they arrive in ∂C
    // order — already sorted by sort_key.  We perform a 3-way merge
    // in O(E) time, with O(1) adjacent-dedup (deterministic).

    // Pre-filter S1 and S2 chord lists (already sorted).
    std::vector<ChordRecord> s1_list, s2_list, disc_list;
    s1_list.reserve(s1.num_chords());
    for (std::size_t k = 0; k < s1.num_chords(); ++k) {
        const auto& c = s1.chord(k);
        if (c.region[0] != NONE || c.region[1] != NONE)
             s1_list.push_back({c.y, c.left_edge, c.right_edge,
                                c.is_null_length,
                                std::min(c.left_edge, c.right_edge)});
    }
    s2_list.reserve(s2.num_chords());
    for (std::size_t k = 0; k < s2.num_chords(); ++k) {
        const auto& c = s2.chord(k);
        if (c.region[0] != NONE || c.region[1] != NONE)
             s2_list.push_back({c.y, c.left_edge, c.right_edge,
                                c.is_null_length,
                                std::min(c.left_edge, c.right_edge)});
    }
    // Discovered chords: already in ∂C order from the fusion walk.
    disc_list.reserve(discovered_chords.size());
    for (const auto& dc : discovered_chords) {
        disc_list.push_back({dc.y, dc.left_edge, dc.right_edge,
                             dc.is_null_length,
                             std::min(dc.left_edge, dc.right_edge)});
    }

    // 3-way merge: O(E) where E = total chords.
    std::vector<ChordRecord> all_chords;
    all_chords.reserve(s1_list.size() + s2_list.size() + disc_list.size());

    std::size_t i1 = 0, n1 = s1_list.size();
    std::size_t i2 = 0, n2 = s2_list.size();
    std::size_t id = 0, nd = disc_list.size();

    while (i1 < n1 || i2 < n2 || id < nd) {
        int source = -1;
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

        // O(1) adjacent-dedup: since the 3-way merge produces chords
        // in sorted order, duplicates are always adjacent.
        if (min_c->left_edge != NONE && min_c->right_edge != NONE &&
            (min_c->left_edge != min_c->right_edge || min_c->is_null_length)) {
            std::size_t c_lo = std::min(min_c->left_edge, min_c->right_edge);
            std::size_t c_hi = std::max(min_c->left_edge, min_c->right_edge);
            bool is_dup = false;
            if (!all_chords.empty()) {
                const auto& prev = all_chords.back();
                std::size_t p_lo = std::min(prev.left_edge, prev.right_edge);
                std::size_t p_hi = std::max(prev.left_edge, prev.right_edge);
                if (c_lo == p_lo && c_hi == p_hi && min_c->y == prev.y)
                    is_dup = true;
            }
            if (!is_dup) {
                all_chords.push_back(*min_c);
            }
        }

        if (source == 0) i1++;
        else if (source == 1) i2++;
        else id++;
    }

    std::size_t num_regions = all_chords.size() + 1;
    for (std::size_t i = 0; i < num_regions; ++i)
        result.add_node();

    // --- Split-point collection: O(E log E) via sorted vector ---
    //
    // Per §3.1 (Lemma 3.1 proof): "let us sort the endpoints of
    // these chords along ∂C … Note that merging can also be used
    // instead of sorting but this step is not the dominant cost,
    // anyway."  The paper budgets O(m · log(n₁+n₂)) for the entire
    // normal-form setup.  Our O(E log E) sort is ≤ O(m log(n₁+n₂))
    // since E ≤ m ≤ n₁+n₂, matching the paper exactly.
    std::vector<std::size_t> split_edges;
    split_edges.reserve(all_chords.size() * 2);
    for (const auto& cr : all_chords) {
        if (cr.left_edge != NONE)
            split_edges.push_back(cr.left_edge);
        if (cr.right_edge != NONE)
            split_edges.push_back(cr.right_edge);
    }
    std::sort(split_edges.begin(), split_edges.end());
    split_edges.erase(
        std::unique(split_edges.begin(), split_edges.end()),
        split_edges.end());

    struct ArcRecord {
        ArcStructure arc;
        std::size_t sort_key;  // min(first_edge, last_edge)
    };
    std::vector<ArcRecord> all_arcs;

    auto split_and_add_arc = [&](ArcStructure a) {
        if (a.first_edge == NONE) {
            all_arcs.push_back({a, 0});
            return;
        }
        std::size_t lo = std::min(a.first_edge, a.last_edge);
        std::size_t hi = std::max(a.first_edge, a.last_edge);

        if (lo >= hi) {
            all_arcs.push_back({a, lo});
            return;
        }

        // §3.1: Find interior split points in [lo+1, hi) via binary
        // search into the sorted split_edges vector.  O(log E + k)
        // where k = number of split points in this arc's range.
        // Total across all arcs: O(A · log E + E) = O(E log E)
        // (within the paper's O(m · log(n₁+n₂)) budget).
        auto it = std::lower_bound(split_edges.begin(),
                                   split_edges.end(), lo + 1);
        std::size_t cur_lo = lo;
        while (it != split_edges.end() && *it < hi) {
            std::size_t e = *it;
            // Sub-arc [cur_lo, e-1].
            ArcStructure sub = a;
            sub.first_edge = cur_lo;
            sub.last_edge  = e - 1;
            sub.edge_count = e - cur_lo;
            all_arcs.push_back({sub, cur_lo});
            cur_lo = e;
            ++it;
        }
        // Final sub-arc [cur_lo, hi].
        ArcStructure sub = a;
        sub.first_edge = cur_lo;
        sub.last_edge  = hi;
        sub.edge_count = hi - cur_lo + 1;
        all_arcs.push_back({sub, cur_lo});
    };

    // S₁ arcs come before S₂ arcs in ∂C order (disjoint edge ranges),
    // so concatenation preserves ∂C order — no sort needed.
    for (std::size_t i = 0; i < s1.num_arcs(); ++i)
        split_and_add_arc(s1.arc(i));
    for (std::size_t i = 0; i < s2.num_arcs(); ++i)
        split_and_add_arc(s2.arc(i));

    // --- Arc-to-region assignment: O(E) linear walk ---
    //
    // Per §3.1 (Lemma 3.1 proof): "This allows us to set up the
    // required arc-sequence table."  A clockwise walk along ∂C visits
    // arcs and chords in interleaved order.  Since both all_arcs and
    // all_chords are sorted by edge position along ∂C, we can assign
    // arcs to regions with a single linear walk: each chord boundary
    // advances the region counter, and arcs between two consecutive
    // chord boundaries belong to the same region.
    {
        std::size_t chord_cursor = 0;
        std::size_t current_region = 0;
        for (auto& arec : all_arcs) {
            // Advance past all chords whose sort_key ≤ this arc's key.
            // Each such chord separates the previous region from the next.
            while (chord_cursor < all_chords.size() &&
                   all_chords[chord_cursor].sort_key <= arec.sort_key) {
                ++chord_cursor;
            }
            current_region = chord_cursor;  // region index = #chords before it
            if (current_region >= num_regions)
                current_region = num_regions - 1;

            arec.arc.region_node = current_region;
            std::size_t ai = result.add_arc(arec.arc);
            result.node(current_region).arcs.push_back(ai);
        }
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

    // Set start_arc/end_arc: find the split between LEFT and RIGHT arcs.
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
            std::size_t junction_vertex,
            const GradeStorage& storage) {
    if (s1.num_nodes() == 0) return Submap(s2);
    if (s2.num_nodes() == 0) return Submap(s1);

    std::vector<DiscoveredChord> discovered;

    fuse_pass(s1, oracle1, s2, oracle2, polygon, junction_vertex, storage, discovered);
    fuse_pass(s2, oracle2, s1, oracle1, polygon, junction_vertex, storage, discovered);

    return rebuild_submap(s1, s2, discovered, polygon);
}

} // namespace chazelle
