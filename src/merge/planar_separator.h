#pragma once

// [LT79 §3 tex 520–668]: Lipton–Tarjan planar-separator machinery, needed
// by [C91 §3.4 tex 297]: "The second reward is that we can apply the
// linear-time algorithm of Lipton and Tarjan [21] to find a good
// separator."  Truth sources: papers/lipton-tarjan1979-transcribed.tex
// (cited as [LT79 ... tex N]) and papers/chazelle1991-transcribed.tex.
//
// [C91 §3.4 tex 297–302]: the separator "partitions the nodes of G into
// three subsets A, B, D, such that
//   (i)   no edge joins a node of A with a node of B,
//   (ii)  neither A nor B contains more than 2μ/3 nodes, and
//   (iii) D contains at most sqrt(8μ) nodes."
// This is [LT79 Corollary sqrt-sep tex 382–393] (2·sqrt(2n) = sqrt(8n)),
// obtained from [LT79 Theorem thm:main tex 316–323] with unit costs 1/n.
//
// [C91 §3.4 tex 304]: the recursive decomposition — "We repeat the
// procedure over with respect to each of G_A and G_B and iterate in this
// fashion until none of the graphs have more than μ^δ nodes, for some
// fixed δ (0 < δ < 1).  Let D* be the set of all separators ... We
// easily verify that |D*| = O(μ^δ) provided that δ is chosen large
// enough; for example, δ = 2/3 [22].  In O(μ log μ) time we can compute
// D* and partition the remaining nodes into subsets D₁, D₂, etc., each
// of size at most μ^{2/3}, such that no path of G can join two nodes in
// distinct subsets without passing through a node of D*."

#include "../common.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace chazelle {

// [LT79 §3 tex 524–528]: "we need a good representation of a planar
// embedding of a graph.  For this purpose we use a list structure whose
// elements correspond to the edges of the graph.  Stored with each edge
// are its endpoints and four pointers, designating the edges immediately
// clockwise and counter-clockwise around each of the endpoints."
//
// Half-edge encoding: undirected edge e has directed halves 2e (u→v)
// and 2e+1 (v→u); the four [LT79] pointers are the cyclic `next`/`prev`
// rotation links of the two halves around their origins.
//
// [LT79 §3 tex 534–538] Step 1 asks for a planar embedding via the
// planarity algorithm of Hopcroft–Tarjan; our caller ([C91 §3.4]'s dual
// graph G of S*) inherits its embedding from the geometry of S*, so the
// rotation lists are supplied by the caller instead.
class EmbeddedPlanarGraph {
public:
    // `rotations[v]` = the edge indices incident to v in clockwise
    // embedding order around v.  Every edge {u, v} of `edges` must
    // appear exactly once in rotations[u] and once in rotations[v]
    // (twice in rotations[u] == rotations[v] is impossible: loops are
    // excluded; parallel edges are allowed by the representation but
    // [C91 §3.4 tex 295] G is simple, and the separator asserts
    // simplicity where its correctness needs it).
    EmbeddedPlanarGraph(std::size_t num_vertices,
                        std::vector<std::pair<std::size_t, std::size_t>> edges,
                        const std::vector<std::vector<std::size_t>>& rotations);

    std::size_t num_vertices() const noexcept { return num_vertices_; }
    std::size_t num_edges()    const noexcept { return edges_.size(); }

    std::size_t edge_u(std::size_t e) const noexcept { return edges_[e].first; }
    std::size_t edge_v(std::size_t e) const noexcept { return edges_[e].second; }

    // Half-edge navigation.  Half h = 2e + b travels FROM origin(h) TO
    // to(h); twin(h) = h ^ 1.
    std::size_t half_to(std::size_t h) const noexcept {
        return (h & 1) ? edges_[h >> 1].first : edges_[h >> 1].second;
    }
    std::size_t half_origin(std::size_t h) const noexcept {
        return half_to(h ^ 1);
    }
    // Clockwise rotation successor / predecessor around origin(h).
    std::size_t rot_next(std::size_t h) const noexcept { return nxt_[h]; }
    std::size_t rot_prev(std::size_t h) const noexcept { return prv_[h]; }

    // Plain adjacency (rotation order), for BFS etc.
    const std::vector<std::size_t>& incident_halves(std::size_t v) const noexcept {
        return rot_[v];
    }

    // Induced subgraph on `nodes` (rotation order preserved — an induced
    // subgraph of an embedded planar graph inherits its embedding).
    // `nodes` need not be sorted.  [C91 §3.4 tex 304]: "G_A ... obtained
    // by keeping only the nodes of A and the edges of G that join only
    // nodes of A."  Also returns old-index mapping via out param.
    EmbeddedPlanarGraph induced(const std::vector<std::size_t>& nodes,
                                std::vector<std::size_t>* old_index) const;

private:
    std::size_t num_vertices_ = 0;
    std::vector<std::pair<std::size_t, std::size_t>> edges_;
    std::vector<std::size_t> nxt_, prv_;       // per half-edge
    std::vector<std::vector<std::size_t>> rot_; // per vertex: half-edge ids, cw
};

// [C91 §3.4 tex 297–302] / [LT79 Corollary sqrt-sep tex 382–388]:
// partition labels for one separator application.
enum class SepPart : std::uint8_t { A = 0, B = 1, D = 2 };

// One application of the separator theorem with unit costs
// ([LT79 §3 tex 531–650] Steps 2–10; Step 1's embedding is the input).
// Postconditions (asserted):
//   (i)   no edge joins A and B                    [C91 §3.4 tex 299]
//   (ii)  3·|A| ≤ 2n and 3·|B| ≤ 2n                [C91 §3.4 tex 300]
//   (iii) |D|² ≤ 8n                                [C91 §3.4 tex 301]
// Time: O(n + m) ([LT79 §3 tex 652–653]: "All steps except Step 9
// obviously run in O(n) time" + the Step 9 proof at tex 656–667).
std::vector<SepPart> planar_separator(const EmbeddedPlanarGraph& g);

// [C91 §3.4 tex 304]: D* and the D_i partition, δ = 2/3.
struct SeparatorDecomposition {
    // node → subset index (D_i), or NONE for nodes of D*.
    std::vector<std::size_t> subset;
    std::size_t num_subsets = 0;
    std::size_t dstar_size = 0;
};

// Recurse with planar_separator until no piece has more than μ^{2/3}
// nodes ([C91 §3.4 tex 304]; integer test: piece³ ≤ μ²).  Postconditions
// (asserted): every |D_i| satisfies |D_i|³ ≤ μ², and
// |D*| ≤ 19·μ^{2/3} — the paper states |D*| = O(μ^{2/3}) [C91 §3.4 tex
// 304]; the constant 19 is proven in the .cpp (geometric size-class sum
// over the recursion tree, √8·√(3/2)/(1 − √(2/3)) < 19).
// The "no path of G joins two nodes in distinct subsets without passing
// through D*" property holds structurally: every edge leaving a
// recursion piece ends in an ancestor separator ([LT79 tex 198–199]
// property (i) inductively).
// Time: O(μ log μ) ([C91 §3.4 tex 304]).
SeparatorDecomposition build_separator_decomposition(const EmbeddedPlanarGraph& g);

} // namespace chazelle
