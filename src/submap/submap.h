#pragma once

// [C91 §2.4]: Normal-form submap.
// (i) tree (regions + chords), (ii) chord↔arc adjacency, (iii) arc-sequence
// in ∂C traversal order, (iv) tree decomposition if conformal.

#include "chord.h"
#include "arc.h"
#include "tree_decomposition.h"
#include "../common.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

namespace chazelle {

// [C91 §2.4(i)]: Tree node (region).
struct SubmapNode {
    // Incident chord indices.  [C91 §2.3]: conformal ⟹ degree ≤ 4.
    std::vector<std::size_t> incident_chords;

    // Number of live (non-dead) incident chords.
    std::size_t degree() const noexcept;

    bool dead = false;          // [C91 §3.3]: tombstone for O(1) removal.
};

class Submap {
public:
    // ── Construction ────────────────────────────────────────────

    std::size_t add_node();
    std::size_t add_arc(Arc arc);
    std::size_t add_chord(Chord chord);

    // [C91 §2.2 tex 94]: Remove a chord — "remove the chord and those
    // endpoints that are not vertices of C, gluing back ∂C at those
    // points."  Arcs are merged only at non-vertex endpoints.
    // @param polygon  needed to detect which endpoints are polygon vertices.
    // @return surviving region index.
    std::size_t remove_chord(std::size_t chord_idx,
                              const class Polygon& polygon);

    // ── Invariants ──────────────────────────────────────────────

    // [C91 §2.2]: "the dual graph of a submap is itself a tree" —
    // num_live_regions == num_live_chords + 1.
    void assert_tree_property() const;

    // Structural invariants: tree property; live chord regions / arc
    // region_nodes; arc-sequence order (LEFT first / ascending, RIGHT /
    // descending); start_arc / end_arc validity.
    void check_invariants() const;

    // Polygon-dependent invariants on top of `check_invariants()`:
    //   [C91 §2.4(iii) tex 138 + tex 144]: arcs sharing first_edge are
    //     start-y monotonic (required by double_identify Phase 2 bsearch).
    //   [C91 §2.2 tex 106]: every live arc's `edge_count` cache matches
    //     `polygon.count_nonnull_edges` over its underlying edge range.
    void check_invariants(const class Polygon& polygon) const;

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

    // ── [C91 §2.4(iii)]: C endpoints ──────────────────────────────────

    // [C91 §2.4(iii) tex 138]: "endpoints of C are identified by appropriate
    // pointers into the input table" + arc-structures passing through them.
    // We hold VERTEX-table indices instead of the paper's edge-table ones —
    // equivalent via the edge-vertex bijection (edge i = vertex i→i+1).
    std::size_t start_vertex = NONE;
    std::size_t end_vertex   = NONE;
    std::size_t start_arc    = NONE;
    std::size_t end_arc      = NONE;

    // [C91 §2.4]: First RIGHT arc in arc-sequence.  LEFT = [0, boundary),
    // RIGHT = [boundary, num_arcs).
    std::size_t left_right_boundary() const noexcept {
        return left_right_boundary_;
    }

    // ── [C91 §2.4]: Double identification ────────────────────────────

    // [C91 §2.4 tex 144]: "double identification of a point of ∂C" — all
    // arcs passing through (edge, symbolic_y).  O(log m).
    // At most 6 arcs (worst case: y-extremum with chords on both sides).
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
                                          SymbolicY y,
                                          const class Polygon& polygon) const;

    // ── [C91 §2.2 / §2.3]: Properties ──────────────────────────────────────

    // [C91 §2.2]: Region weight = max nonnull-edge count over its arcs;
    // 0 for empty regions.
    std::size_t region_weight(std::size_t node_idx) const noexcept;

    // [C91 §2.3]: conformal = node-degree ≤ 4 everywhere.
    bool is_conformal() const noexcept;

    // [C91 §2.3]: γ-semigranular = every region weight ≤ γ.
    bool is_semigranular(std::size_t gamma) const noexcept;

    // [C91 §2.3]: γ-granular = (i) all weights ≤ γ AND (ii) contracting
    // any edge incident on a < 3-degree node yields a region of weight
    // > γ.  By default, (i)-only with no exit chord is still γ-granular.
    bool is_granular(std::size_t gamma,
                      const class Polygon& polygon) const noexcept;

    // [C91 §2.3]: Weight of the merged region if `chord_idx` were contracted.
    // May be less than the sum: non-vertex chord endpoints disappear
    // ([C91 §2.2 tex 94]).
    std::size_t simulated_contraction_weight(
        std::size_t chord_idx,
        const class Polygon& polygon) const noexcept;

    // ── Compaction ──────────────────────────────────────────────

    // Strip dead arcs/chords/nodes; rebuild index mappings.  O(m).
    // Called once before putting S in normal form ([C91 §3.3]).
    void compact();

    std::size_t num_live_nodes()  const noexcept;
    std::size_t num_live_chords() const noexcept;
    std::size_t num_live_arcs()   const noexcept;

    // ── [C91 §2.4(iv)]: Tree decomposition ────────────────────────────

    // [C91 §2.4 tex 133]: arc-structures don't store endpoint y per the
    // paper ("we do not need to record the endpoints of the arc because
    // chords take care of that").  Derive the symbolic y of the arc's
    // starting position on ∂C on demand:
    //  - chord-bounded arcs: chord.y via the region's incident_chords
    //    adjacency ([C91 §2.4(ii) tex 137]);
    //  - polygon-vertex-bounded arcs (outside-pair companion at an
    //    interior y-extremum per [C91 §2.1 tex 72], or a junction after
    //    a non-merging chord removal per [C91 §2.2 tex 94]): the input
    //    table's vertex y ([C91 §2.4(iii) tex 138]).
    // O(degree) per call; O(1) for conformal submaps (degree ≤ 4).
    SymbolicY arc_start_symbolic_y(std::size_t arc_idx,
                                    const class Polygon& polygon) const;

    // [C91 §2.4(iv)]: conformal ⟹ tree decomposition available.
    void build_tree_decomposition();
    const TreeDecomposition& tree_decomposition() const noexcept {
        // [C91 §3.3 tex 276] requires mutators to run in O(1), so rather
        // than rebuilding the tree decomposition every time a chord or
        // arc changes, mutators just set a dirty flag.  When the flag is
        // set, this accessor returns an empty TreeDecomposition so that
        // any consumer requiring a fresh one (e.g. merge.h's
        // `!tree_decomposition().empty()` precondition) fails fast
        // instead of using the out-of-date one.
        static const TreeDecomposition empty_;
        return tree_decomp_dirty_ ? empty_ : tree_decomp_;
    }

private:
    TreeDecomposition tree_decomp_;
    // Set by every mutator; cleared by build_tree_decomposition().  The
    // out-of-date contents stay until the next build() rebuilds from
    // scratch (or the Submap destructor cleans up).
    bool tree_decomp_dirty_ = false;

    std::vector<SubmapNode> nodes_;     // [C91 §2.4(i)]: tree nodes.
    std::vector<Chord> chords_;         // [C91 §2.4(ii)]: tree edges.
    // [C91 §2.4(iii)]: arc-structures in canonical ∂C traversal order.
    std::vector<Arc> arc_sequence_;

    // [C91 §2.4]: First RIGHT arc index — enables O(1) LEFT/RIGHT split
    // for double_identify's O(log m) bound (tex 144).
    std::size_t left_right_boundary_ = 0;

    // [C91 §2.4 tex 144]: double_identify needs a compacted arc-sequence;
    // tracked as O(1) flag rather than O(m) per-call scan.
    bool compacted_ = true;
};

} // namespace chazelle
