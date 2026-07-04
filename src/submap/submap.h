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
    // points."
    //
    // Arc gluing happens at EVERY endpoint of the removed chord, not just
    // non-vertex ones: [C91 §2.2 tex 96] region boundaries ALTERNATE exit
    // chords with arcs (Fig. 2.4's region II — "a two-edge arc, an exit
    // chord, a one-edge arc, an exit chord, a three-edge arc (not a
    // four-edge arc!), ..." — passes through removed chords' vertex
    // endpoints inside single arcs), and [C91 §2.2 tex 108] "once removed,
    // a chord of zero length ceases to separate any arcs" (null-length
    // chords have only vertex endpoints).  The vertex/non-vertex
    // distinction of tex 94 governs ∂C's VERTEX set (mid-edge endpoints
    // vanish; vertices of C persist as interior arc vertices), not arc
    // structure.
    //
    // Exception: junctions at C's own endpoint companions ([C91 §2.1
    // tex 72] case 3) are NOT glued — the arc-sequence table keeps arcs
    // split at C's two endpoint wraps (the settled representation of
    // [C91 §2.4 tex 142] double-backing; §3.2's LogicalArc re-glues wrap
    // pieces logically where paper arcs are needed).
    //
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

    // [C91 §2.4(iii) tex 138]: "pointers to the arc-structures whose
    // corresponding arcs pass through the endpoints" — TWO arcs pass
    // through each endpoint of C.  start_arc/end_arc hold the LEFT-half
    // pair; tail_arc is the arc ENDING at C's start wrap (the last RIGHT
    // arc in canonical ∂C order).  Needed because that arc ends at no
    // chord, so it appears in no adjacency slot when its start is a
    // vertex endpoint: region_weight and the [C91 §2.2 tex 96] junction
    // glue would otherwise not reach it in O(1).  Maintained by add_arc
    // (canonical insertion order), remove_chord, insert_chord, compact,
    // and normalize.
    std::size_t tail_arc     = NONE;

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
    // May be less than the sum: "this weight might be less than the added
    // weight of the two nodes of the contracted edge" ([C91 §2.3 tex 123]).
    // Simulates exactly what remove_chord would produce: arc chains glued
    // at both endpoints (vertex and mid-edge; wrap-cap junctions excepted).
    std::size_t simulated_contraction_weight(
        std::size_t chord_idx,
        const class Polygon& polygon) const noexcept;

    // Same simulation, but with caller-supplied weights for the chord's
    // two regions instead of region_weight().  [C91 §3.3 tex 276] needs
    // the criterion-(ii) test in O(1) per chord while removals cascade;
    // the §3.3 pass maintains per-region weights itself and passes them
    // here (w0 for chord.region[0], w1 for chord.region[1]).
    std::size_t simulated_contraction_weight(
        std::size_t chord_idx,
        const class Polygon& polygon,
        std::size_t w0, std::size_t w1) const noexcept;

    // ── Compaction ──────────────────────────────────────────────

    // Strip dead arcs/chords/nodes; rebuild index mappings.  O(m).
    // Called once before putting S in normal form ([C91 §3.3]).
    void compact();

    // ── [C91 §3.3 tex 276]: Normal form ─────────────────────────

    // "We can now put S in normal form, which includes computing its
    // tree decomposition.  As we discussed earlier, this can be done
    // very simply in time O((n₁/γ₁ + n₂/γ₂ + 1)log(n₁ + n₂))."
    //
    // compact()s the tables, re-sorts the arc-sequence into canonical
    // ∂C traversal order ([C91 §2.4(iii) tex 138]: LEFT ascending
    // first_edge then RIGHT descending; within an edge, start-y follows
    // the clockwise traversal direction) — repairing the order broken by
    // [C91 §3.2]'s insert_chord appends — restores start_arc/end_arc and
    // the LEFT/RIGHT boundary, re-validates the [C91 §2.4] invariants,
    // and rebuilds the tree decomposition ([C91 §2.4(iv)]).
    //
    // O(M log M) comparison sort + O(M log M) tree decomposition
    // ([C91 §2.3 tex 116]), within tex 276's normal-form budget.
    //
    // Precondition (asserted): S is conformal — normal form (iv) requires
    // a tree decomposition, which is defined for conformal submaps.
    void normalize(const class Polygon& polygon);

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

    // Symmetric to `arc_start_symbolic_y`: the symbolic y of the arc's
    // ENDING position on ∂C.  [C91 §2.4 tex 133]: derived from the
    // bounding chord (slot 0 / count-1 slot = the arc ending at the
    // chord per the §2.4(ii) adjacency convention), or from the input
    // table's vertex y when the arc ends at a polygon vertex.
    // O(degree) per call.  Needed by [C91 §3.2]'s arc splitting and hit
    // attribution.
    SymbolicY arc_end_symbolic_y(std::size_t arc_idx,
                                  const class Polygon& polygon) const;

    // ── [C91 §3.2]: Chord insertion ───────────────────────────────────

    // Position of a new chord endpoint on ∂C, inside a region arc.
    struct ChordPointSpec {
        std::size_t arc;        // live arc of the region containing the point
        std::size_t edge;       // input-table edge of the point
        Side side;              // ∂C side of the point
        double x;               // exact x of the point on the edge at y
    };

    struct InsertChordResult {
        std::size_t chord_idx;
        std::size_t new_region;     // == chord(chord_idx).region[1]
        std::size_t p_after_arc;    // appended second half of p.arc
        std::size_t q_after_arc;    // appended second half of q.arc
    };

    // [C91 §3.2 tex 264]: "adding new chords into S" — insert the
    // visibility chord pq (horizontal at symbolic y) into `region`,
    // splitting it into two regions.  p and q must be mutually visible
    // with respect to C (caller establishes this via Lemma 3.2) and must
    // lie on two DISTINCT arcs of `region`.
    //
    // `cycle` = the region's live arcs in clockwise ∂C traversal order
    // (Lemma 2.2 [C91 §2.2 tex 96]: the ∂C-order subsequence of a
    // region's arcs IS its boundary cycle, counterclockwise).  The caller
    // ([C91 §3.2]'s driver) maintains this inventory.
    //
    // Splits p.arc / q.arc in two (the first half keeps the table slot,
    // the second half is APPENDED — canonical arc-sequence order is NOT
    // maintained; [C91 §3.3 tex 276] "We can now put S in normal form"
    // owns re-normalization, so `compacted_` is cleared to make
    // double_identify fail fast until then).
    //
    // Preconditions (asserted): neither p nor q coincides with an
    // existing chord endpoint of the region — a ∂C point's horizontal
    // visibility is unique ([C91 §2.1 tex 70]), so a chord endpoint's
    // visibility is already realized by its chord.
    //
    // O(|cycle| + degree) — O(1) for [C91 §3.2]'s bounded-arc regions
    // (tex 238), preserving Lemma 3.4's bound.
    InsertChordResult insert_chord(const ChordPointSpec& p,
                                    const ChordPointSpec& q,
                                    SymbolicY y,
                                    std::size_t region,
                                    const std::size_t* cycle,
                                    std::size_t cycle_len,
                                    const class Polygon& polygon);

    // [C91 §2.4(iv)]: conformal ⟹ tree decomposition available.
    void build_tree_decomposition();
    const TreeDecomposition& tree_decomposition() const noexcept {
        // [C91 §3.3 tex 276] requires mutators to run in O(1), so rather
        // than rebuilding the tree decomposition every time a chord or
        // arc changes, mutators just set a dirty flag.  When the flag is
        // set, this accessor returns an empty TreeDecomposition so that
        // callers requiring a fresh one (e.g. merge.h's
        // `!tree_decomposition().empty()` precondition) hit that assert
        // instead of reading an out-of-date one.
        static const TreeDecomposition empty_;
        return tree_decomp_dirty_ ? empty_ : tree_decomp_;
    }

private:
    // [C91 §2.2 tex 96]: locate the glue mate at a removed chord's
    // vertex-endpoint junction — the arc STARTING at the junction ∂C
    // point (want_after) or ENDING there (!want_after).  The 1-slot
    // adjacency convention records only the before-arc at vertex
    // endpoints, so the mate is found by scanning the adjacency slots of
    // the chord's two regions' incident chords plus the C-endpoint arcs
    // (start_arc/end_arc/tail_arc) — O(degree) = O(1) for conformal
    // submaps.  Among candidates a zero-length arc wins (it occupies the
    // junction point itself and precedes/follows immediately, [C91 §2.2
    // tex 96] "some arcs may be of zero length").
    // @param edge, side   the junction endpoint's ∂C position.
    // @param vertex_idx   the polygon vertex at the junction.
    // @param exclude, exclude2  arc indices to skip (already-known
    //     junction pieces; NONE if unused).
    std::size_t find_junction_arc(const Chord& c,
                                  std::size_t edge, Side side,
                                  std::size_t vertex_idx,
                                  bool want_after,
                                  std::size_t exclude,
                                  std::size_t exclude2,
                                  const class Polygon& polygon) const;

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
