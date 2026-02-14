#include "chazelle/down_phase.h"
#include "chazelle/conformality.h"
#include "chazelle/fusion.h"
#include "chazelle/granularity.h"
#include "common.h"
#include "oracles/arc_cutting.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>
#include <vector>

namespace chazelle {

namespace {

// ────────────────────────────────────────────────────────────────
// Helper: build a trivial submap for a single polygon edge.
//
// This is used for single-edge stubs from arc-cutting.
// The resulting submap has 1 region, 0 chords, and 1 arc
// spanning the single edge.
// ────────────────────────────────────────────────────────────────
struct TrivialPiece {
    Submap submap;
    RayShootingOracle oracle;
    std::size_t start_vertex;
    std::size_t end_vertex;
};

TrivialPiece make_trivial_submap(std::size_t edge_idx,
                                 std::size_t start_v,
                                 std::size_t end_v,
                                 const Polygon& polygon) {
    TrivialPiece tp;
    tp.start_vertex = start_v;
    tp.end_vertex = end_v;

    // 1 region, 1 arc, 0 chords.
    std::size_t r0 = tp.submap.add_node();
    ArcStructure arc;
    arc.first_edge = edge_idx;
    arc.last_edge  = edge_idx;
    arc.first_side = Side::LEFT;
    arc.last_side  = Side::LEFT;
    arc.region_node = r0;
    arc.edge_count = 1;
    std::size_t ai = tp.submap.add_arc(arc);
    tp.submap.node(r0).arcs.push_back(ai);
    tp.submap.recompute_weight(r0);
    tp.submap.start_arc = ai;
    tp.submap.end_arc   = ai;

    // Build a trivial oracle.
    tp.oracle.build(tp.submap, polygon, 2);
    return tp;
}

// ────────────────────────────────────────────────────────────────
// Build a trivial submap for an exit chord treated as a 3-edge
// polygonal curve (§4.2).
//
// Per Chazelle §4.2: "Let T₁ (resp. T₂) be canonical submaps
// for the 3-edge polygonal curve a₁'a₁b₁b₁' (resp. a₂'a₂b₂b₂')."
// The three edges are:
//   1. a₁'a₁ — stub from nearest ∂C vertex to left chord endpoint
//   2. a₁b₁  — the exit chord itself (tilted symbolically)
//   3. b₁b₁' — stub from right chord endpoint to nearest ∂C vertex
//
// The middle arc is virtual: its virtual_y records the chord
// height so that ray-shooting computes intersections with the
// tilted segment rather than polygon edges.
// ────────────────────────────────────────────────────────────────
TrivialPiece make_exit_chord_submap(std::size_t left_edge,
                                    std::size_t right_edge,
                                    double chord_y,
                                    const Polygon& polygon) {
    TrivialPiece tp;
    tp.start_vertex = left_edge;
    tp.end_vertex   = right_edge;

    // Single region.
    std::size_t r0 = tp.submap.add_node();

    // Arc 1: left stub  a₁'a₁  (portion of polygon edge left_edge).
    if (left_edge < polygon.num_edges()) {
        ArcStructure stub_left;
        stub_left.first_edge  = left_edge;
        stub_left.last_edge   = left_edge;
        stub_left.first_side  = Side::LEFT;
        stub_left.last_side   = Side::LEFT;
        stub_left.region_node = r0;
        stub_left.edge_count  = 1;
        std::size_t ai0 = tp.submap.add_arc(stub_left);
        tp.submap.node(r0).arcs.push_back(ai0);
    }

    // Arc 2: tilted chord  a₁b₁  (virtual edge at chord_y).
    ArcStructure chord_arc;
    chord_arc.first_edge  = left_edge;
    chord_arc.last_edge   = right_edge;
    chord_arc.first_side  = Side::LEFT;
    chord_arc.last_side   = Side::RIGHT;
    chord_arc.region_node = r0;
    chord_arc.edge_count  = 1;
    chord_arc.virtual_y   = chord_y;   // marks this as a virtual arc
    std::size_t ai1 = tp.submap.add_arc(chord_arc);
    tp.submap.node(r0).arcs.push_back(ai1);

    // Arc 3: right stub  b₁b₁'  (portion of polygon edge right_edge).
    if (right_edge < polygon.num_edges() && right_edge != left_edge) {
        ArcStructure stub_right;
        stub_right.first_edge  = right_edge;
        stub_right.last_edge   = right_edge;
        stub_right.first_side  = Side::RIGHT;
        stub_right.last_side   = Side::RIGHT;
        stub_right.region_node = r0;
        stub_right.edge_count  = 1;
        std::size_t ai2 = tp.submap.add_arc(stub_right);
        tp.submap.node(r0).arcs.push_back(ai2);
    }

    tp.submap.recompute_weight(r0);
    tp.submap.start_arc = 0;
    tp.submap.end_arc   = tp.submap.num_arcs() - 1;

    // Build a trivial oracle.
    tp.oracle.build(tp.submap, polygon, 2);
    return tp;
}

// ────────────────────────────────────────────────────────────────
// Walking R's boundary  (§4.2)
//
// A region R's boundary alternates arcs and exit chords in cyclic
// ∂C order.  We produce a sequence of "boundary elements" that
// represents R* — the boundary of R treated as a non-closed curve.
// ────────────────────────────────────────────────────────────────

/// One element of R*'s boundary: either an arc or an exit chord
/// treated as a curve edge.
struct BoundaryElement {
    enum Kind { ARC, EXIT_CHORD };
    Kind kind;
    std::size_t arc_idx;    // valid when kind == ARC
    std::size_t chord_idx;  // valid when kind == EXIT_CHORD
    std::size_t sort_key;   // first_edge for arcs, left_edge for chords
};

/// Walk the boundary of region R in ∂C order.
/// Returns the cyclic sequence of arcs and exit chords.
std::vector<BoundaryElement>
walk_region_boundary(const Submap& submap, std::size_t region_idx) {
    const auto& nd = submap.node(region_idx);

    // Collect arcs with sort keys.
    struct ArcEntry { std::size_t arc_idx; std::size_t sort_key; };
    std::vector<ArcEntry> arcs;
    for (std::size_t ai : nd.arcs) {
        const auto& a = submap.arc(ai);
        if (a.first_edge == NONE) continue;
        std::size_t lo = std::min(a.first_edge, a.last_edge);
        arcs.push_back({ai, lo});
    }
    std::sort(arcs.begin(), arcs.end(),
              [](const ArcEntry& a, const ArcEntry& b) {
                  return a.sort_key < b.sort_key;
              });

    // Collect chords with sort keys (using left_edge).
    struct ChordEntry { std::size_t chord_idx; std::size_t sort_key; };
    std::vector<ChordEntry> chords;
    for (std::size_t ci : nd.incident_chords) {
        const auto& c = submap.chord(ci);
        // Skip deleted / invalid chords.
        if (c.region[0] == NONE || c.region[1] == NONE) continue;
        if (submap.node(c.region[0]).deleted ||
            submap.node(c.region[1]).deleted) continue;
        std::size_t sk = (c.left_edge != NONE) ? c.left_edge : 0;
        chords.push_back({ci, sk});
    }
    std::sort(chords.begin(), chords.end(),
              [](const ChordEntry& a, const ChordEntry& b) {
                  return a.sort_key < b.sort_key;
              });

    // Merge arcs and chords in ∂C order.
    // A chord at left_edge = e goes BEFORE the arc starting at
    // first_edge = e (since the chord endpoint is on edge e, which
    // lies at the beginning of that edge's range).
    std::vector<BoundaryElement> boundary;
    std::size_t ci = 0;
    for (auto& ae : arcs) {
        // Emit any chords whose sort_key ≤ this arc's sort_key.
        while (ci < chords.size() && chords[ci].sort_key <= ae.sort_key) {
            boundary.push_back({BoundaryElement::EXIT_CHORD,
                                NONE, chords[ci].chord_idx,
                                chords[ci].sort_key});
            ++ci;
        }
        boundary.push_back({BoundaryElement::ARC,
                            ae.arc_idx, NONE,
                            ae.sort_key});
    }
    // Remaining chords (wrap-around).
    while (ci < chords.size()) {
        boundary.push_back({BoundaryElement::EXIT_CHORD,
                            NONE, chords[ci].chord_idx,
                            chords[ci].sort_key});
        ++ci;
    }
    return boundary;
}

// ────────────────────────────────────────────────────────────────
// refine_region — Lemma 4.2 of Chazelle 1991
//
// For each region R of the current submap S:
//
// 1. Walk R's boundary in ∂C order to obtain R* — the boundary of
//    R as a non-closed polygonal curve.  R*'s "edges" consist of
//    (a) the polygon arcs of R, each decomposed into prior-grade
//        chain pieces via arc-cutting,   AND
//    (b) the exit chords of R, each treated as a curve edge
//        (not as a chord — §4.2: "we treat the edges a₁b₁ and
//        a₂b₂ as part of the input curve").
//
// 2. Retrieve canonical submaps for each chain piece from grade
//    storage.  Build trivial 1-edge submaps for exit-chord edges
//    and for single-edge stubs from arc-cutting.
//
// 3. Merge all pieces via balanced binary tree (Lemma 4.1)
//    using fuse (§3.1) + restore_conformality (§3.2).
//    This is where cross-chain visibility chords are discovered.
//
// 4. Extract chords from V(R*) whose both endpoints lie on the
//    arcs (not the exit-chord edges) of R.  These are the new
//    chords internal to R.
//
// 5. Insert the extracted chords into the current submap,
//    splitting R into sub-regions.
// ────────────────────────────────────────────────────────────────
void refine_region(Submap& submap,
                   std::size_t region_idx,
                   GradeStorage& storage,
                   const Polygon& polygon,
                   std::size_t current_grade) {
    auto& nd = submap.node(region_idx);
    if (nd.deleted || nd.is_empty()) return;

    // ── Step 1: Walk R's boundary ───────────────────────────────
    auto boundary = walk_region_boundary(submap, region_idx);
    if (boundary.empty()) return;

    // Record which polygon edge ranges belong to R's actual arcs
    // (not exit-chord edges).  Used later for chord extraction.
    struct ArcRange { std::size_t lo, hi; };
    std::vector<ArcRange> region_arc_ranges;
    for (std::size_t ai : nd.arcs) {
        const auto& arc = submap.arc(ai);
        if (arc.first_edge == NONE) continue;
        std::size_t a_lo = std::min(arc.first_edge, arc.last_edge);
        std::size_t a_hi = std::max(arc.first_edge, arc.last_edge);
        region_arc_ranges.push_back({a_lo, a_hi});
    }
    auto edge_in_region = [&](std::size_t e) -> bool {
        for (const auto& r : region_arc_ranges)
            if (e >= r.lo && e <= r.hi) return true;
        return false;
    };

    // ── Step 2: Build the pieces to merge ───────────────────────
    //
    // For each boundary element, produce one or more "merge pieces"
    // (each with a submap + oracle + vertex range).
    //
    // A merge piece is either:
    //   - A canonical submap from grade storage (for a chain piece
    //     from arc-cutting), OR
    //   - A trivial 1-region submap (for a single-edge stub from
    //     arc-cutting, or an exit chord treated as a curve edge).

    struct MergePiece {
        // Exactly one of these is set:
        const CanonicalSubmap* stored = nullptr;  // from grade storage
        TrivialPiece* trivial = nullptr;          // owned locally
        std::size_t start_vertex;
        std::size_t end_vertex;
    };

    // Locally-owned trivial submaps (exit chords + single-edge stubs).
    std::vector<std::unique_ptr<TrivialPiece>> trivial_pieces;

    std::vector<MergePiece> merge_pieces;

    for (auto& be : boundary) {
        if (be.kind == BoundaryElement::ARC) {
            // Decompose this arc into chain pieces via arc-cutting.
            const auto& arc = submap.arc(be.arc_idx);
            if (arc.first_edge == NONE) continue;
            std::size_t start = std::min(arc.first_edge, arc.last_edge);
            std::size_t end   = std::max(arc.first_edge, arc.last_edge) + 1;
            auto pieces = cut_arc(start, end, polygon);

            for (auto& piece : pieces) {
                if (piece.grade < storage.num_grades() &&
                    piece.chain_index < storage.num_chains(piece.grade)) {
                    // Chain piece from grade storage.
                    auto& cs = storage.get(piece.grade, piece.chain_index);
                    std::size_t gp = compute_granularity(piece.grade);
                    if (!cs.oracle.is_built())
                        cs.oracle.build(cs.submap, polygon, gp);
                    else
                        cs.oracle.rebind_submap(cs.submap);

                    MergePiece mp;
                    mp.stored = &cs;
                    mp.start_vertex = piece.start_vertex;
                    mp.end_vertex = piece.end_vertex;
                    merge_pieces.push_back(mp);
                } else {
                    // Single-edge stub — build trivial submap.
                    auto tp = std::make_unique<TrivialPiece>(
                        make_trivial_submap(piece.start_vertex,
                                            piece.start_vertex,
                                            piece.end_vertex,
                                            polygon));
                    // Rebind oracle: make_trivial_submap built the
                    // oracle pointing to a stack temporary's submap;
                    // after move-constructing on the heap the pointer
                    // is stale.
                    tp->oracle.rebind_submap(tp->submap);
                    MergePiece mp;
                    mp.trivial = tp.get();
                    mp.start_vertex = piece.start_vertex;
                    mp.end_vertex = piece.end_vertex;
                    trivial_pieces.push_back(std::move(tp));
                    merge_pieces.push_back(mp);
                }
            }
        } else {
            // EXIT_CHORD → treat as a 3-edge curve (§4.2).
            //
            // §4.2: "we treat the edges a₁b₁ and a₂b₂ as part of
            // the input curve although they are not part of P."
            //
            // "Let T₁ (resp. T₂) be canonical submaps for the
            // 3-edge polygonal curve a₁'a₁b₁b₁' (resp. a₂'a₂b₂b₂').
            // ... T₁ and T₂ are computed directly (tilting the edges
            // a₁b₁ and a₂b₂ symbolically to keep the merging
            // algorithm from complaining later)."
            //
            // The exit chord connects left_edge to right_edge at
            // height y.  We build a 3-arc submap: left stub, tilted
            // chord (virtual arc), right stub.
            const auto& chord = submap.chord(be.chord_idx);
            std::size_t le = chord.left_edge;
            std::size_t re = chord.right_edge;
            if (le == NONE || re == NONE) continue;
            if (le == re) continue; // null-length chord → skip

            auto tp = std::make_unique<TrivialPiece>(
                make_exit_chord_submap(le, re, chord.y, polygon));
            tp->oracle.rebind_submap(tp->submap);
            MergePiece mp;
            mp.trivial = tp.get();
            mp.start_vertex = le;
            mp.end_vertex   = re;
            trivial_pieces.push_back(std::move(tp));
            merge_pieces.push_back(mp);
        }
    }

    if (merge_pieces.size() < 2) return; // nothing to merge

    // ── Step 3: Balanced binary tree merge (§4.2 / Lemma 4.1) ───
    //
    // Per the paper: "we consider a perfectly balanced binary tree
    // whose leaves are in bijection with the D_i and we merge
    // submaps bottom-up by following the tree pattern."
    //
    // We merge pairs bottom-up: (0,1), (2,3), … then merge those
    // results pairwise, etc., until one submap remains.

    auto get_submap = [](const MergePiece& mp) -> const Submap& {
        return mp.stored ? mp.stored->submap : mp.trivial->submap;
    };

    // §4.2: the pieces have granularity at most 2^⌈β⌈βλ⌉⌉ = fine_gamma.
    // The merged submap of V(R*) should be fine_gamma-granular conformal.
    std::size_t coarse_exp = (current_grade + 4) / 5;       // ⌈βλ⌉
    std::size_t fine_exp   = (coarse_exp + 4) / 5;          // ⌈β⌈βλ⌉⌉
    std::size_t gamma_merge = std::size_t(1) << fine_exp;   // 2^⌈β⌈βλ⌉⌉

    // Layer 0: copy each piece's submap+oracle into the working set.
    // Per Lemma 4.1: "trivially reset the granularity of each S_i to γ"
    // — enforce granularity with gamma_merge on each piece before merging.
    struct TreeNode {
        Submap submap;
        RayShootingOracle oracle;
        std::size_t start_vertex;
        std::size_t end_vertex;
    };
    std::vector<TreeNode> layer;
    layer.reserve(merge_pieces.size());
    for (auto& mp : merge_pieces) {
        TreeNode tn;
        tn.submap = get_submap(mp);
        tn.submap.set_chain_info(mp.start_vertex, mp.end_vertex, &polygon);
        enforce_granularity(tn.submap, gamma_merge);
        tn.oracle.build(tn.submap, polygon, gamma_merge); // Removed: Unsafe to store pointer to stack variable.
        tn.start_vertex = mp.start_vertex;
        tn.end_vertex = mp.end_vertex;
        layer.push_back(std::move(tn));
    }

    // Fix 1: Build oracles after nodes are in their stable vector location.
    for (auto& node : layer) {
        node.oracle.build(node.submap, polygon, gamma_merge);
    }

    // Bottom-up pairwise merge until one node remains.
    while (layer.size() > 1) {
        std::vector<TreeNode> next_layer;
        next_layer.reserve((layer.size() + 1) / 2); // Ensure safe push_back without reallocation

        for (std::size_t i = 0; i + 1 < layer.size(); i += 2) {
            // Junction vertex between layer[i] and layer[i+1].
            std::size_t junction = layer[i].end_vertex;

            // Ensure oracles are built.
            if (!layer[i].oracle.is_built())
                layer[i].oracle.build(layer[i].submap, polygon, gamma_merge);
            if (!layer[i+1].oracle.is_built())
                layer[i+1].oracle.build(layer[i+1].submap, polygon, gamma_merge);

            // Fuse the two submaps.
            Submap fused = fuse(layer[i].submap, layer[i].oracle,
                                layer[i+1].submap, layer[i+1].oracle,
                                polygon, junction, storage);

            // Restore conformality on the fused result.
            restore_conformality(fused, storage, polygon, gamma_merge);

            // §4.2 / Lemma 4.1: enforce granularity after conformality.
            enforce_granularity(fused, gamma_merge);

            TreeNode merged;
            merged.submap = std::move(fused);
            merged.start_vertex = layer[i].start_vertex;
            merged.end_vertex = layer[i+1].end_vertex;
            // §2.3: set chain info and normalize for correct
            // double_identify and C-vertex guards.
            merged.submap.set_chain_info(merged.start_vertex,
                                         merged.end_vertex, &polygon);
            merged.submap.normalize();
            // Build oracle for the next level.
            merged.oracle.build(merged.submap, polygon, gamma_merge);
            next_layer.push_back(std::move(merged));
        }
        // If odd number of nodes, carry the last one up unchanged.
        if (layer.size() % 2 == 1) {
            next_layer.push_back(std::move(layer.back()));
        }
        layer = std::move(next_layer);

        // Fix 2: Rebind oracles. When nodes were moved from 'next_layer' (or stack)
        // to 'layer', they changed address. The oracles internally hold pointers
        // to the old addresses. We must update them to point to the current addresses.
        for (auto& node : layer) {
            node.oracle.rebind_submap(node.submap);
        }
    }

    Submap& merged_submap = layer[0].submap;

    // ── Step 4: Extract chords from V(R*) internal to R ─────────
    //
    // Per §4.2: "extract the relevant information, i.e., the exit
    // chords falling entirely within each region R.  This involves
    // checking the exit chords of the computed submap of V(R*) and
    // keeping only those both of whose endpoints lie on the arcs
    // (in the double boundary sense) of the region R."
    //
    // A chord is internal to R if BOTH its left_edge and
    // right_edge edges fall within R's actual arcs (not within
    // the exit-chord-turned-curve-edges).

    struct InternalChord {
        double y;
        std::size_t left_edge;
        std::size_t right_edge;
    };
    std::vector<InternalChord> internal_chords;

    // Deduplication set.
    std::set<std::pair<std::size_t, std::size_t>> seen;

    for (std::size_t ci = 0; ci < merged_submap.num_chords(); ++ci) {
        const auto& chord = merged_submap.chord(ci);
        if (chord.left_edge == NONE || chord.right_edge == NONE)
            continue;

        // Skip self-loop chords (le == re) that are NOT null-length.
        // These are degenerate artefacts that can't partition arcs.
        // Genuine null-length chords (at y-extrema) are handled
        // separately by add_missing_null_length_chords().
        if (chord.left_edge == chord.right_edge &&
            !chord.is_null_length)
            continue;

        // Deduplicate by edge pair.
        auto key = std::make_pair(
            std::min(chord.left_edge, chord.right_edge),
            std::max(chord.left_edge, chord.right_edge));
        if (!seen.insert(key).second) continue;

        // Check both endpoints lie on R's actual arcs.
        bool left_in = edge_in_region(chord.left_edge);

        bool right_in = edge_in_region(chord.right_edge);

        if (left_in && right_in) {
            internal_chords.push_back(
                {chord.y, chord.left_edge, chord.right_edge});
        }
    }

    if (internal_chords.empty()) return;

    // ── Step 5: Insert internal chords into the submap ───────────

    // §4.2: "extract the relevant information, i.e., the exit chords
    // falling entirely within each region R.  This involves checking
    // the exit chords of the computed submap of V(R*) and keeping
    // only those both of whose endpoints lie on the arcs."
    //
    // The paper constructs S* directly — no sort is prescribed.
    // Non-crossing horizontal chord insertions are order-independent,
    // so we process chords in whatever order they were collected.

    // Build a set of edge pairs already present in the parent submap.
    // This prevents re-inserting chords that are still active.
    // Chords removed by granularity enforcement are fully invalidated
    // (left_edge = NONE) by remove_chord(), so they are naturally
    // excluded.  This allows them to be rediscovered at finer grades.
    std::set<std::pair<std::size_t, std::size_t>> parent_chords;
    for (std::size_t ci = 0; ci < submap.num_chords(); ++ci) {
        const auto& c = submap.chord(ci);
        if (c.left_edge == NONE || c.right_edge == NONE) continue;
        if (c.region[0] == NONE || c.region[1] == NONE) continue;
        if (submap.node(c.region[0]).deleted ||
            submap.node(c.region[1]).deleted) continue;
        parent_chords.emplace(std::min(c.left_edge, c.right_edge),
                              std::max(c.left_edge, c.right_edge));
    }

    // Track all sub-regions created from the original region_idx.
    // As chords are inserted, arcs move to new sub-regions; later
    // chords must search the correct sub-region.
    std::vector<std::size_t> live_subregions = {region_idx};

    for (auto& ic : internal_chords) {
        // Skip chords that already exist in the parent submap.
        auto key = std::make_pair(
            std::min(ic.left_edge, ic.right_edge),
            std::max(ic.left_edge, ic.right_edge));
        if (parent_chords.count(key)) continue;

        std::size_t le_edge = ic.left_edge;

        // Find which sub-region currently holds the arc containing
        // this chord's left endpoint edge.
        std::size_t target = NONE;
        for (std::size_t sr : live_subregions) {
            for (std::size_t ai : submap.node(sr).arcs) {
                const auto& a = submap.arc(ai);
                if (a.first_edge == NONE) continue;
                std::size_t a_lo = std::min(a.first_edge, a.last_edge);
                std::size_t a_hi = std::max(a.first_edge, a.last_edge);
                if (le_edge >= a_lo && le_edge <= a_hi) {
                    target = sr;
                    break;
                }
            }
            if (target != NONE) break;
        }
        if (target == NONE) continue;

        std::size_t new_region = submap.add_node();

        Chord new_chord;
        new_chord.y = ic.y;
        new_chord.left_edge = ic.left_edge;
        new_chord.right_edge = ic.right_edge;
        new_chord.region[0] = target;
        new_chord.region[1] = new_region;

        std::size_t new_chord_idx = submap.add_chord(new_chord);
        parent_chords.insert(key);

        // Null-length chords create an empty region; no arc partition.
        if (ic.left_edge == ic.right_edge) {
            submap.recompute_weight(target);
            submap.recompute_weight(new_region);
            continue;
        }

        // Partition arcs between target and new_region.
        std::size_t split_lo = (ic.left_edge != NONE)
                                   ? ic.left_edge : 0;
        std::size_t split_hi = (ic.right_edge != NONE)
                                   ? ic.right_edge : split_lo;
        if (split_lo > split_hi) std::swap(split_lo, split_hi);

        auto& arcs = submap.node(target).arcs;
        std::vector<std::size_t> keep_arcs;
        std::vector<std::size_t> move_arcs;

        for (std::size_t ai : arcs) {
            auto& a = submap.arc(ai);
            if (a.first_edge == NONE) {
                keep_arcs.push_back(ai);
                continue;
            }
            std::size_t a_lo = std::min(a.first_edge, a.last_edge);
            std::size_t a_hi = std::max(a.first_edge, a.last_edge);
            if (a_lo < split_hi && a_hi >= split_lo) {
                move_arcs.push_back(ai);
            } else {
                keep_arcs.push_back(ai);
            }
        }

        if (move_arcs.empty() && keep_arcs.size() > 1) {
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

        arcs = keep_arcs;
        for (std::size_t ai : move_arcs) {
            submap.arc(ai).region_node = new_region;
            submap.node(new_region).arcs.push_back(ai);
        }

        auto& ichords = submap.node(target).incident_chords;
        std::vector<std::size_t> keep_chords;
        std::vector<std::size_t> move_chords;

        for (std::size_t ci : ichords) {
            if (ci == new_chord_idx) {
                // The splitting chord itself stays with the target —
                // don't reassign its regions.
                keep_chords.push_back(ci);
                continue;
            }
            auto& c = submap.chord(ci);
            std::size_t cv = c.left_edge;
            if (cv != NONE && cv >= split_lo && cv < split_hi) {
                move_chords.push_back(ci);
            } else {
                keep_chords.push_back(ci);
            }
        }

        ichords = keep_chords;
        for (std::size_t ci : move_chords) {
            for (int s = 0; s < 2; ++s) {
                if (submap.chord(ci).region[s] == target) {
                    submap.chord(ci).region[s] = new_region;
                }
            }
            submap.node(new_region).incident_chords.push_back(ci);
        }

        submap.recompute_weight(target);
        submap.recompute_weight(new_region);

        // Track the new sub-region for subsequent chord insertions.
        live_subregions.push_back(new_region);
    }
}

/// Process one grade in the down-phase.
///
/// When \p is_final is true we are at the last grade and the goal is
/// the *complete* visibility map V(P), not a granular submap.  Per
/// Lemma 4.2 base case: "the regions of S have bounded size, and
/// therefore the missing chords can be provided in constant time per
/// region."  No granularity enforcement is performed at the final
/// grade — doing so would remove null-length chords at y-extrema
/// whose empty regions are needed by FM Algorithm 2 to locate all
/// Class B trapezoids and guarantee monotone sub-polygons.
void process_grade(Submap& submap,
                   GradeStorage& storage,
                   const Polygon& polygon,
                   std::size_t grade,
                   bool is_final) {
    // §4.2: the result of refinement is 2^⌈β⌈βλ⌉⌉-semigranular.
    //
    // β = 1/5.  compute_granularity(λ) = 2^⌈λ/5⌉  (the coarse γ).
    // The finer granularity after re-merge is 2^⌈β⌈βλ⌉⌉:
    //   coarse_exp = ⌈λ/5⌉
    //   fine_exp   = ⌈coarse_exp / 5⌉ = ⌈⌈λ/5⌉ / 5⌉
    //   fine_gamma = 2^fine_exp
    std::size_t coarse_exp = (grade + 4) / 5;           // ⌈λ/5⌉ = ⌈βλ⌉
    std::size_t fine_exp   = (coarse_exp + 4) / 5;      // ⌈⌈βλ⌉/5⌉ = ⌈β⌈βλ⌉⌉
    std::size_t fine_gamma = std::size_t(1) << fine_exp;

    // Snapshot the original chord count BEFORE refinement.
    // Per §4.2: "some of the chords connecting the R*'s might be
    // removable now.  We can check each of the exit chords directly."
    // Only ORIGINAL exit chords should be checked for granularity
    // removal — not the newly inserted chords from the re-merge.
    std::size_t original_chord_count = submap.num_chords();

    // Collect all live region indices (snapshot, since we'll modify).
    std::vector<std::size_t> regions;
    for (std::size_t i = 0; i < submap.num_nodes(); ++i) {
        if (!submap.node(i).deleted) {
            regions.push_back(i);
        }
    }

    // Refine each region.
    for (std::size_t r : regions) {
        if (submap.node(r).deleted) continue;
        std::size_t pre_chords = submap.num_chords();
        refine_region(submap, r, storage, polygon, grade);
        std::size_t post_chords = submap.num_chords();
        if (post_chords > pre_chords) {
            // New chords were added during refinement.
        }
    }

    // Restore conformality after chord reinterpretation (§4.2).
    std::size_t pre_conf = submap.num_chords();
    restore_conformality(submap, storage, polygon, fine_gamma);
    std::size_t post_conf = submap.num_chords();
    if (post_conf > pre_conf) {
        // Conformality restoration added new chords.
    }

    // Enforce granularity — but NOT at the final grade, where the
    // goal is the complete V(P) (Chazelle §4.2 / Theorem 4.3).
    //
    // Per §4.2: only check the ORIGINAL exit chords for removability,
    // not the newly inserted ones.  We use enforce_granularity with
    // the fine granularity and only check chords that existed before
    // the refinement step.
    if (!is_final) {
        enforce_granularity(submap, fine_gamma, original_chord_count);
    }

    // §2.3: Restore normal form.  Refinement, conformality restoration,
    // and granularity enforcement may have appended arcs out of order.
    // Normalizing restores the sorted arc-sequence table required for
    // correct double_identify at the next grade.
    submap.normalize();
}

// §4.2: Null-length chords at y-extrema are discovered naturally by
// the algorithm during the down-phase refinement.  No separate
// post-processing step is needed.

} // anonymous namespace

Submap down_phase(const Polygon& polygon, GradeStorage& storage) {
    std::size_t p = polygon.num_grades();

    // Start from the top-grade canonical submap (the whole polygon).
    assert(storage.num_grades() > 0);
    Submap current = storage.get(p, 0).submap;



    // Process grades following the paper's recursive structure (Lemma 4.2):
    //   λ → ⌈βλ⌉ → ⌈β⌈βλ⌉⌉ → … → constant.
    // At each step, the submap is 2^⌈βλ⌉-granular; after refinement
    // it becomes 2^⌈β⌈βλ⌉⌉-granular, and the next λ' = ⌈βλ⌉.
    //
    // With β = 1/5 the sequence is e.g. p=20 → 4 → 1.
    // The final step (λ reaching a small constant) produces V(P).
    std::size_t lambda = p;
    while (lambda >= 1) {
        std::size_t next_lambda = (lambda + 4) / 5;  // ⌈βλ⌉ = ⌈λ/5⌉
        bool is_final = (next_lambda == 0 || next_lambda >= lambda);
        process_grade(current, storage, polygon, lambda, is_final);
        if (is_final) break;
        lambda = next_lambda;
    }

    // §4.2: Null-length chords at y-extrema are discovered naturally
    // by the algorithm during the down-phase refinement.

    // The visibility map is complete.
    return current;
}

} // namespace chazelle
