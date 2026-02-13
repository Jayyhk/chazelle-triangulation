#pragma once

/// Submap — the central data structure of Chazelle's algorithm.
///
/// A submap is a coarsened visibility map obtained by retaining a subset
/// of chords of V(C).  It stores:
///   1. Submap tree (unrooted): nodes = regions, edges = chords.
///   2. Chord records on each tree edge.
///   3. Arc-sequence table: arcs in canonical ∂C traversal order.
///   4. Tree decomposition (if conformal).
///
/// A submap in "normal form" satisfies conditions (i)–(iv) of §2.3.
/// A submap is "conformal" if every tree node has degree ≤ 4.
/// A submap is "γ-granular" if it's conformal, γ-semigranular (every
/// region has weight ≤ γ), AND maximal (contracting any edge incident
/// on a degree-<3 node would exceed γ).
///
/// Key operations:
///   - Insert/remove chords (used by fusion, conformality, granularity).
///   - Compute region weights.
///   - Build tree decomposition (for conformal submaps).
///   - Double identification: locate arc(s) passing through a given point.

#include "arc_structure.h"
#include "chord.h"

#include <cstddef>
#include <vector>
#include <limits>

namespace chazelle {

/// A node (region) in the submap tree.
struct SubmapNode {
    /// Indices of incident chords (tree edges) in Submap::chords_.
    /// For a conformal submap, this has at most 4 entries.
    std::vector<std::size_t> incident_chords;

    /// Indices of arcs belonging to this region (in Submap::arc_sequence_).
    std::vector<std::size_t> arcs;

    /// Cached weight of this region.
    /// Weight = 0 if empty, else max edge_count over all arcs.
    std::size_t weight = 0;

    /// Is this region empty (weight 0)?
    bool is_empty() const { return weight == 0; }

    /// Degree in the submap tree.
    std::size_t degree() const { return incident_chords.size(); }

    /// Soft-delete flag.
    bool deleted = false;
};

class TreeDecomposition; // forward declaration

class Submap {
public:
    Submap() = default;

    // --- Construction ---

    /// Add a region node.  Returns its index.
    std::size_t add_node();

    /// Add an arc to the arc-sequence table.  Returns its index.
    std::size_t add_arc(ArcStructure arc);

    /// Add a chord (tree edge) connecting two regions.
    /// Also updates the incident_chords lists and adj_arcs pointers.
    /// Returns the chord index.
    std::size_t add_chord(Chord chord);

    /// Split an arc at a given vertex index.  If the arc's edge range
    /// [first_edge, last_edge] contains vertex_idx in its interior,
    /// splits it into two arcs: [first_edge, vertex_idx-1] and
    /// [vertex_idx, last_edge].  Both arcs belong to the same region.
    /// Returns true if a split was performed.
    bool split_arc_at_vertex(std::size_t arc_idx, std::size_t vertex_idx);

    /// Remove a chord (contract the tree edge, merging two regions).
    /// The region with the higher index is absorbed into the one with
    /// the lower index.  Returns the surviving region index.
    std::size_t remove_chord(std::size_t chord_idx);

    // --- Accessors ---

    std::size_t num_nodes() const { return nodes_.size(); }
    std::size_t num_chords() const { return chords_.size(); }
    std::size_t num_arcs() const { return arc_sequence_.size(); }

    const SubmapNode& node(std::size_t i) const { return nodes_[i]; }
    SubmapNode& node(std::size_t i) { return nodes_[i]; }

    const Chord& chord(std::size_t i) const { return chords_[i]; }
    Chord& chord(std::size_t i) { return chords_[i]; }

    const ArcStructure& arc(std::size_t i) const { return arc_sequence_[i]; }
    ArcStructure& arc(std::size_t i) { return arc_sequence_[i]; }

    const std::vector<SubmapNode>& nodes() const { return nodes_; }
    const std::vector<Chord>& chords() const { return chords_; }
    const std::vector<ArcStructure>& arc_sequence() const { return arc_sequence_; }

    // --- Weight ---

    /// Recompute the weight of a single region from its arcs.
    void recompute_weight(std::size_t node_idx);

    /// Recompute all region weights.
    void recompute_all_weights();

    // --- Queries ---

    /// Is this submap conformal (every node has degree ≤ 4)?
    bool is_conformal() const;

    /// Is this submap γ-semigranular (every region has weight ≤ γ)?
    bool is_semigranular(std::size_t gamma) const;

    /// Max degree of any node in the tree.
    std::size_t max_degree() const;

    // --- Double Identification (§2.4) ---

    /// Given an edge index in the input table and a y-coordinate,
    /// find the arc(s) that pass through that point.
    /// Returns indices into arc_sequence_.
    /// Binary search in the arc-sequence table: split circular sequence
    /// at the two C-endpoint arcs into two monotone subsequences, binary
    /// search each by edge name, disambiguate by y-coordinate.  O(log m).
    std::vector<std::size_t> double_identify(std::size_t edge_idx,
                                              double y) const;

    // --- Endpoint info ---

    /// Index of the arc passing through the start endpoint of C.
    std::size_t start_arc = NONE;
    /// Index of the arc passing through the end endpoint of C.
    std::size_t end_arc = NONE;

private:
    std::vector<SubmapNode> nodes_;
    std::vector<Chord> chords_;
    std::vector<ArcStructure> arc_sequence_;
};

/// Type alias: the complete visibility map is just a Submap with all chords.
using VisibilityMap = Submap;

} // namespace chazelle
