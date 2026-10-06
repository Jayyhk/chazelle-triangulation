#pragma once

#include "../common.h"
#include "arc.h"
#include "chord.h"
#include "tree_decomposition.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

namespace chazelle {

struct SubmapNode {
    std::vector<std::size_t> incident_chords;

    std::size_t degree() const noexcept;

    bool dead = false;
};

class Submap {
public:
    std::size_t add_node();
    std::size_t add_arc(Arc arc);
    std::size_t add_chord(Chord chord);

    std::size_t remove_chord(std::size_t chord_idx, const class Polygon& polygon);

    void assert_tree_property() const;

    void check_invariants() const;

    void check_invariants(const class Polygon& polygon) const;

    std::size_t num_nodes() const noexcept {
        return nodes_.size();
    }
    std::size_t num_chords() const noexcept {
        return chords_.size();
    }
    std::size_t num_arcs() const noexcept {
        return arc_sequence_.size();
    }

    const SubmapNode& node(std::size_t i) const noexcept {
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

    std::size_t start_vertex = NONE;
    std::size_t end_vertex = NONE;
    std::size_t start_arc = NONE;
    std::size_t end_arc = NONE;

    std::size_t left_right_boundary() const noexcept {
        return left_right_boundary_;
    }

    struct DoubleIdentifyResult {
        static constexpr std::size_t MAX = 6;
        std::array<std::size_t, MAX> arcs = {};
        std::size_t count = 0;
        void push(std::size_t arc_idx) {
            for (std::size_t i = 0; i < count; ++i)
                if (arcs[i] == arc_idx)
                    return;
            assert(count < MAX);
            arcs[count++] = arc_idx;
        }
        const std::size_t* begin() const {
            return arcs.data();
        }
        const std::size_t* end() const {
            return arcs.data() + count;
        }
    };

    DoubleIdentifyResult double_identify(std::size_t edge_idx, SymbolicY y,
                                         const class Polygon& polygon) const;

    std::size_t region_weight(std::size_t node_idx) const noexcept;

    bool is_conformal() const noexcept;

    bool is_semigranular(std::size_t granularity) const noexcept;

    bool is_granular(std::size_t granularity, const class Polygon& polygon) const noexcept;

    std::size_t simulated_contraction_weight(std::size_t chord_idx,
                                             const class Polygon& polygon) const noexcept;

    std::size_t simulated_contraction_weight(std::size_t chord_idx, const class Polygon& polygon,
                                             std::size_t w0, std::size_t w1) const noexcept;

    void chord_regions_below_above(std::size_t chord_idx, const class Polygon& polygon,
                                   std::size_t* below, std::size_t* above) const;

    void compact();

    void normalize(const class Polygon& polygon);

    std::size_t num_live_nodes() const noexcept;
    std::size_t num_live_chords() const noexcept;
    std::size_t num_live_arcs() const noexcept;

    SymbolicY arc_start_symbolic_y(std::size_t arc_idx, const class Polygon& polygon) const;

    SymbolicY arc_end_symbolic_y(std::size_t arc_idx, const class Polygon& polygon) const;

    void refresh_arc_edge_counts(const class Polygon& polygon);

    struct ChordPointSpec {
        std::size_t arc;
        std::size_t edge;
        Side side;
        Exact x;
    };

    struct InsertChordResult {
        std::size_t chord_idx;
        std::size_t new_region;
        std::size_t p_after_arc;
        std::size_t q_after_arc;
    };

    InsertChordResult insert_chord(const ChordPointSpec& p, const ChordPointSpec& q, SymbolicY y,
                                   std::size_t region, const std::size_t* cycle,
                                   std::size_t cycle_len, const class Polygon& polygon);

    void build_tree_decomposition();
    const TreeDecomposition& tree_decomposition() const noexcept {
        static const TreeDecomposition empty_;
        return tree_decomp_dirty_ ? empty_ : tree_decomp_;
    }

private:
    bool arc_is_point(std::size_t arc_idx, const class Polygon& polygon) const;

    std::size_t find_junction_arc(const Chord& c, bool query_left, std::size_t edge, Side side,
                                  std::size_t vertex_idx, bool want_after, std::size_t exclude,
                                  std::size_t exclude2, const class Polygon& polygon) const;

    TreeDecomposition tree_decomp_;

    bool tree_decomp_dirty_ = false;

    std::vector<SubmapNode> nodes_;
    std::vector<Chord> chords_;

    std::vector<Arc> arc_sequence_;

    std::size_t left_right_boundary_ = 0;

    bool compacted_ = true;

    std::size_t live_chords_ = 0;
};

}
