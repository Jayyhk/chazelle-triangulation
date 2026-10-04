#pragma once

#include "../polygon/polygon.h"
#include "submap.h"

#include <cstddef>
#include <vector>

namespace chazelle {

struct PendingChord {
    SymbolicY y;
    std::size_t left_edge_c;
    Side left_side;
    std::size_t right_edge_c;
    Side right_side;
    bool is_null_length = false;
};

struct ChordEndpoint {
    std::size_t edge_c;
    Side side;
    SymbolicY y;
    std::size_t pending_idx;
    bool is_left_slot;
};

void canonicalize_chord(PendingChord& chord, const Polygon& curve);

bool chord_endpoint_precedes(const Polygon& curve, const std::vector<PendingChord>& chords,
                             const ChordEndpoint& first, const ChordEndpoint& second);

void build_submap_from_ordered_chords(Submap& out_S, const Polygon& curve,
                                      const std::vector<PendingChord>& chords,
                                      const std::vector<ChordEndpoint>& endpoints);

void build_submap_from_chords(Submap& out_S, const Polygon& curve,
                              std::vector<PendingChord> chords);

}
