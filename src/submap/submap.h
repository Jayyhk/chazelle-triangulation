#pragma once

/// [C91 §2.4]: Normal-form submap representation.
///
/// (i)   Tree in standard edge/node adjacency fashion.
/// (ii)  Chords → adj arcs; arcs → region node.
/// (iii) Arc-sequence table in ∂C traversal order.
/// (iv)  Tree decomposition if conformal.

#include "chord.h"
#include "arc.h"
#include "tree_decomposition.h"
#include "../common.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

namespace chazelle {

/// [C91 §2.4(i)]: A node (region) of the submap tree.
struct SubmapNode {
    /// Indices of incident chords in the Submap's chord table.
    /// [C91 §2.3]: conformal ⟹ degree ≤ 4.
    std::vector<std::size_t> incident_chords;

    /// Degree = number of live (non-dead) incident chords.
    std::size_t degree() const noexcept;

    /// Tombstone flag for O(1) removal (§3.3).
    bool dead = false;
};

class Submap {
public:
    // ── Construction ────────────────────────────────────────────

    /// Add a region node.  Returns its index.
    std::size_t add_node();

    /// Add an arc to the arc-sequence table.  Returns its index.
    std::size_t add_arc(Arc arc);

    /// Add a chord connecting two regions.  Returns its index.
    /// Updates both regions' incident_chords lists.
    std::size_t add_chord(Chord chord);

    /// [C91 §2.2]: Remove a chord.  "The removal of a chord entails
    /// removing not only the chord itself but also those endpoints
    /// that are not vertices of C, and glueing back ∂C at those
    /// points."  Only merges arcs at non-vertex endpoints.
    /// Returns the surviving region index.
    ///
    /// @param polygon  The input polygon — needed to determine which
    ///                 chord endpoints are polygon vertices (and thus
    ///                 should NOT be cleaned up).
    std::size_t remove_chord(std::size_t chord_idx,
                              const class Polygon& polygon);

    // ── Invariant checks ────────────────────────────────────────

    /// [C91 §2.2]: "the dual graph of a submap is itself a tree."
    /// num_live_regions == num_live_chords + 1.
    void assert_tree_property() const;

    /// [C91 §2.2]: Check all structural invariants.
    ///   - Tree property
    ///   - Every live chord's regions are live nodes
    ///   - Every live arc's region_node is a live node
    void check_invariants() const;

    // ── Accessors ───────────────────────────────────────────────

    std::size_t num_nodes()  const noexcept { return nodes_.size(); }
    std::size_t num_chords() const noexcept { return chords_.size(); }
    std::size_t num_arcs()   const noexcept { return arc_sequence_.size(); }

    const SubmapNode& node(std::size_t i)  const noexcept {
        assert(i < nodes_.size());
        return nodes_[i];
    }
    SubmapNode& node(std::size_t i) noexcept {
        assert(i < nodes_.size());
        return nodes_[i];
    }

    const Chord& chord(std::size_t i) const noexcept {
        assert(i < chords_.size());
        return chords_[i];
    }
    Chord& chord(std::size_t i) noexcept {
        assert(i < chords_.size());
        return chords_[i];
    }

    const Arc& arc(std::size_t i) const noexcept {
        assert(i < arc_sequence_.size());
        return arc_sequence_[i];
    }
    Arc& arc(std::size_t i) noexcept {
        assert(i < arc_sequence_.size());
        return arc_sequence_[i];
    }

    // ── §2.4(iii): Endpoint pointers ────────────────────────────

    /// [C91 §2.4(iii)]: "the endpoints of C are identified by
    /// appropriate pointers into the input table as well as by
    /// pointers to the arc-structures whose corresponding arcs
    /// pass through the endpoints."
    std::size_t start_vertex = NONE; ///< Input table: first vertex of C.
    std::size_t end_vertex   = NONE; ///< Input table: last vertex of C.
    std::size_t start_arc    = NONE; ///< Arc-sequence: arc at C's start.
    std::size_t end_arc      = NONE; ///< Arc-sequence: arc at C's end.

    // ── §2.4: Double identification ───────────────────────────────

    /// [C91 §2.4]: "we call it the double identification of a point
    /// of ∂C."  Given an edge and a symbolic y-coordinate, find all
    /// arcs passing through that point.  O(log m).
    /// Returns at most 6 arcs (worst case: y-extremum vertex with
    /// chords on both sides).
    struct DoubleIdentifyResult {
        static constexpr std::size_t MAX = 6;
        std::array<std::size_t, MAX> arcs = {};
        std::size_t count = 0;
        void push(std::size_t arc_idx) {
            assert(count < MAX);
            arcs[count++] = arc_idx;
        }
        const std::size_t* begin() const { return arcs.data(); }
        const std::size_t* end()   const { return arcs.data() + count; }
    };

    DoubleIdentifyResult double_identify(std::size_t edge_idx,
                                          SymbolicY y) const;

    // ── §2.2 / §2.3: Submap properties ─────────────────────────

    /// [C91 §2.2]: Weight of a region = max nonnull edge count
    /// over its arcs.  Returns 0 for empty regions.
    std::size_t region_weight(std::size_t node_idx) const noexcept;

    /// [C91 §2.3]: "conformal submaps [are] those with node-degree
    /// at most 4."
    bool is_conformal() const noexcept;

    /// [C91 §2.3]: "If only condition (i) holds, then the submap is
    /// γ-semigranular."  Every region weight ≤ γ.
    bool is_semigranular(std::size_t gamma) const noexcept;

    /// [C91 §2.3]: γ-granular = (i) all weights ≤ γ AND (ii)
    /// "contracting any edge incident upon at least one node of
    /// degree less than 3 produces a new node whose weight exceeds γ."
    /// "By default, if (i) holds but the submap has no exit chord,
    /// it is still said to be γ-granular."
    bool is_granular(std::size_t gamma,
                      const class Polygon& polygon) const noexcept;

    /// [C91 §2.3]: Simulate contracting a chord — compute what the
    /// merged region's weight would be.  The weight may be less than
    /// the sum because "one or both endpoints of the chord might not
    /// be vertices of ∂C and might thus disappear."
    /// Requires Polygon to determine which endpoints are vertices
    /// (§2.2 tex 94: only non-vertex endpoints merge).
    std::size_t simulated_contraction_weight(
        std::size_t chord_idx,
        const class Polygon& polygon) const noexcept;

    // ── Tombstone compaction ────────────────────────────────────

    /// Strip all dead arcs, chords, and nodes.  Rebuilds index
    /// mappings in O(m).  Called once before putting S in normal
    /// form (§3.3: "We can now put S in normal form").
    void compact();

    /// Count of live (non-dead) nodes, chords, arcs.
    std::size_t num_live_nodes()  const noexcept;
    std::size_t num_live_chords() const noexcept;
    std::size_t num_live_arcs()   const noexcept;

    // ── §2.4(iv): Tree decomposition ────────────────────────────

    /// [C91 §2.4(iv)]: "If the submap is conformal, then its tree
    /// decomposition should be available."
    void build_tree_decomposition();
    const TreeDecomposition& tree_decomposition() const noexcept {
        return tree_decomp_;
    }

private:
    TreeDecomposition tree_decomp_;
    /// [C91 §2.4(i)]: Nodes (regions) of the submap tree.
    std::vector<SubmapNode> nodes_;

    /// [C91 §2.4(ii)]: Chords (edges of the submap tree).
    std::vector<Chord> chords_;

    /// [C91 §2.4(iii)]: "The arc-structures are stored in a table
    /// in the order corresponding to a canonical traversal of ∂C."
    std::vector<Arc> arc_sequence_;

    /// [C91 §2.4]: Cached index of first RIGHT arc in arc_sequence_.
    /// Enables O(1) LEFT/RIGHT split for double_identify (§2.4 tex 144:
    /// "logarithmic in the number of arcs").
    std::size_t left_right_boundary_ = 0;
};

} // namespace chazelle
