#pragma once

/// Ray-Shooting Oracle — Lemma 3.6 of Chazelle 1991.
///
/// Given a horizontal ray from a point on ∂C, find the first chord or arc
/// of a submap S that the ray hits.
///
/// Construction (preprocessing):
///   1. Build S* — a planar subdivision from S by collapsing the double
///      boundary to zero thickness. Edges of S* have positive length.
///   2. Compute dual graph G of S* (planar, O(m/γ + 1) nodes).
///   3. Apply iterated Lipton-Tarjan separator to G → separator hierarchy.
///   4. Build vertical-line structure for initial region identification.
///
/// Query: O(γ^{1/3} · m^{2/3}) time.
///   - Walk the separator hierarchy; at each level, test separator faces.
///   - Descend into the appropriate sub-piece; repeat.

#include "visibility/submap.h"
#include "separator/separator.h"

#include <cstddef>
#include <optional>
#include <vector>

namespace chazelle {

/// Result of a ray-shooting query.
struct RayHit {
    enum class Type { CHORD, ARC, NONE };

    Type type = Type::NONE;

    /// If type == CHORD: the index of the chord hit.
    std::size_t chord_idx = NONE;

    /// If type == ARC: the index of the arc hit.
    std::size_t arc_idx = NONE;

    /// The region the ray starts in (before hitting anything).
    std::size_t start_region = NONE;

    /// X-coordinate of the hit point.
    double hit_x = 0.0;

    /// The polygon edge index containing the hit point.
    /// Per §3 item (i): "the report should also include the name
    /// of the edge of P that contains [the hit point]."
    std::size_t hit_edge = NONE;
};

/// Preprocessed ray-shooting structure for a single submap.
class RayShootingOracle {
public:
    RayShootingOracle() = default;

    /// Build the oracle for a given submap.
    /// Time: O(m log m / γ + 1) where m = #vertices of the curve,
    ///       γ = granularity.
    void build(const Submap& submap, const class Polygon& polygon,
               std::size_t granularity);

    /// Shoot a horizontal ray from a point on ∂C.
    /// @param edge_idx  The polygon edge containing the ray origin.
    /// @param y         The y-coordinate of the ray origin.
    /// @param side      Which side of the double boundary.
    /// @param shoot_right  True for rightward ray, false for leftward.
    /// @return The first chord or arc hit.
    RayHit shoot(std::size_t edge_idx, double y, Side side,
                 bool shoot_right) const;

    /// Shoot from a specific region (when region is already known).
    /// origin_x is the x-coordinate of the ray's starting point.
    RayHit shoot_from_region(std::size_t region_idx, double origin_x,
                             double y, bool shoot_right) const;

    /// §4.1 oracle (i): Shoot a horizontal ray from an EXTERNAL point
    /// (origin_x, y) — i.e., a point NOT on this submap's chain ∂C_μ.
    /// Uses the vertical-line structure to identify the starting region,
    /// then delegates to shoot_from_region with the caller-supplied
    /// origin_x.  This is the cross-chain query needed by the per-
    /// subarc oracle decomposition.
    RayHit shoot_from_point(double origin_x, double y,
                            bool shoot_right) const;

    /// Is this oracle built and ready for queries?
    bool is_built() const { return built_; }

    /// Update the submap pointer after the submap has been moved in
    /// memory (e.g. after being stored in a vector).  The oracle's
    /// internal structures (dual graph, separator hierarchy, vertical
    /// line) remain valid since they don't store pointers into the
    /// submap's internal arrays.
    void rebind_submap(const Submap& submap) { submap_ = &submap; }

    /// Access the underlying submap (needed for vertex resolution).
    const Submap& submap_ref() const { return *submap_; }

private:
    bool built_ = false;

    /// Reference to the submap (non-owning).
    const Submap* submap_ = nullptr;
    const Polygon* polygon_ = nullptr;
    std::size_t granularity_ = 0;

    /// Dual graph G of S*.
    PlanarGraph dual_graph_;

    /// Separator hierarchy from iterated Lipton-Tarjan.
    IteratedSeparatorResult separator_hierarchy_;

    /// Compact list of separator dual-node indices (|D*| = O(μ^{2/3})).
    /// Pre-built from separator_hierarchy_.separator_nodes to avoid
    /// scanning the full O(μ) boolean array during queries.
    std::vector<std::size_t> separator_list_;

    /// Vertical-line structure: sorted list of chord y-coordinates
    /// for initial region identification.
    struct VerticalEntry {
        double y;              ///< Y-coordinate of chord intersection.
        std::size_t chord_idx; ///< Which chord.
        std::size_t region_above; ///< Region above this chord.
        std::size_t region_below; ///< Region below this chord.
    };
    std::vector<VerticalEntry> vertical_line_;

    /// Build the dual graph G of S*.
    void build_dual_graph();

    /// Build the vertical-line structure.
    void build_vertical_line();

    /// Local shooting: test all arcs/chords of a region against a ray.
    /// origin_x is the x-coordinate of the ray's starting point on ∂C.
    RayHit local_shoot(std::size_t region_idx, double origin_x,
                       double y, bool shoot_right) const;

    /// Mapping between submap region indices and dual-graph node indices.
    std::vector<std::size_t> region_to_dual_;  ///< region → dual node
    std::vector<std::size_t> dual_to_region_;  ///< dual node → region

    /// Look up the dual-graph node for a submap region.
    std::size_t region_to_dual(std::size_t region_idx) const {
        if (region_idx < region_to_dual_.size())
            return region_to_dual_[region_idx];
        return NONE;
    }

    /// Look up the submap region for a dual-graph node.
    std::size_t dual_to_region(std::size_t dual_idx) const {
        if (dual_idx < dual_to_region_.size())
            return dual_to_region_[dual_idx];
        return NONE;
    }
};

} // namespace chazelle
