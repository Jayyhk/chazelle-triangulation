#pragma once

/// DCEL-style planar graph for the Lipton-Tarjan separator.
///
/// Representation follows Lipton & Tarjan 1979, §3: each undirected edge
/// stores its two endpoints and four link pointers (CW and CCW around
/// each endpoint).  Each vertex stores one incident edge.
///
/// This is NOT a half-edge data structure — each undirected edge is a
/// single record with four pointers.

#include "common.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace chazelle {

/// An undirected edge in the planar graph.
/// `endpoint[s]` is the vertex on side `s` (s ∈ {0,1}).
/// `cw[s]` / `ccw[s]` are the next edges clockwise / counter-clockwise
/// around `endpoint[s]` in the planar embedding.
struct PlanarEdge {
    std::size_t endpoint[2] = {NONE, NONE};
    std::size_t cw[2]       = {NONE, NONE}; ///< Edge indices.
    std::size_t ccw[2]      = {NONE, NONE}; ///< Edge indices.

    bool is_tree  = false; ///< Marked during BFS spanning tree construction.
    bool deleted  = false; ///< Soft-delete flag for shrink operations.

    /// Return the side (0 or 1) that vertex `v` is on.
    /// Precondition: v == endpoint[0] or v == endpoint[1].
    int side_of(std::size_t v) const {
        if (endpoint[0] == v) return 0;
        assert(endpoint[1] == v);
        return 1;
    }

    /// Return the other endpoint.
    std::size_t other(std::size_t v) const {
        return endpoint[1 - side_of(v)];
    }
};

/// A vertex in the planar graph.
struct PlanarVertex {
    std::size_t some_edge = NONE; ///< Index of any incident edge.
    double      cost      = 1.0; ///< Vertex cost (for the separator).
    bool        deleted   = false;

    // BFS data (populated by separator algorithm):
    std::size_t level      = NONE; ///< BFS level.
    std::size_t parent     = NONE; ///< Parent vertex in BFS tree.
    std::size_t parent_edge= NONE; ///< Tree edge to parent.
    double      desc_cost  = 0.0;  ///< Total cost of descendants (incl self).
};

/// A mutable planar graph in edge-list representation.
class PlanarGraph {
public:
    PlanarGraph() = default;

    /// Pre-allocate `nv` vertices and `ne` edges.
    PlanarGraph(std::size_t nv, std::size_t ne);

    // ── Sizes ───────────────────────────────────────────────────

    std::size_t num_vertices() const { return vertices_.size(); }
    std::size_t num_edges()    const { return edges_.size(); }

    // ── Element access ──────────────────────────────────────────

    PlanarVertex&       vertex(std::size_t i)       { return vertices_[i]; }
    const PlanarVertex& vertex(std::size_t i) const { return vertices_[i]; }

    PlanarEdge&       edge(std::size_t i)       { return edges_[i]; }
    const PlanarEdge& edge(std::size_t i) const { return edges_[i]; }

    // ── Construction helpers ────────────────────────────────────

    /// Add a vertex with the given cost.  Returns vertex index.
    std::size_t add_vertex(double cost = 1.0);

    /// Add an undirected edge between u and v, splicing it into the
    /// CW/CCW orbits at both endpoints.  If `after_u` / `after_v`
    /// are not NONE, the new edge is inserted immediately CCW of
    /// `after_u` around u (and similarly for v) in the embedding.
    /// Returns edge index.
    std::size_t add_edge(std::size_t u, std::size_t v,
                         std::size_t after_u = NONE,
                         std::size_t after_v = NONE);

    /// Insert edge `e` into the CW/CCW orbit of vertex `v`,
    /// immediately after (CCW of) edge `after` in the orbit.
    /// If `after` is NONE and vertex has no edges, initialise orbit.
    void splice_into_orbit(std::size_t e, std::size_t v,
                           std::size_t after);

    // ── Traversal helpers ───────────────────────────────────────

    /// Next edge CW around vertex v, starting from edge e.
    std::size_t next_cw(std::size_t e, std::size_t v) const {
        return edges_[e].cw[edges_[e].side_of(v)];
    }

    /// Next edge CCW around vertex v, starting from edge e.
    std::size_t next_ccw(std::size_t e, std::size_t v) const {
        return edges_[e].ccw[edges_[e].side_of(v)];
    }

    /// Iterate over all edges incident to vertex v in CW order,
    /// calling f(edge_index) for each.  Skips deleted edges.
    /// Includes a traversal guard to avoid infinite loops if the
    /// CW chain passes through deleted edges.
    template <typename F>
    void for_each_edge_cw(std::size_t v, F&& f) const {
        std::size_t start = vertices_[v].some_edge;
        if (start == NONE) return;
        // If the starting edge is deleted, walk CW to find a live one.
        // If all edges are deleted, we give up after a full orbit.
        if (edges_[start].deleted) {
            std::size_t e = next_cw(start, v);
            std::size_t g = edges_.size();
            while (e != start && edges_[e].deleted && g-- > 0) {
                e = next_cw(e, v);
            }
            if (edges_[e].deleted) return; // all edges deleted
            start = e;
        }
        std::size_t e = start;
        std::size_t guard = edges_.size() + 1; // max iterations
        do {
            if (!edges_[e].deleted) {
                f(e);
            }
            e = next_cw(e, v);
            if (--guard == 0) break; // safety: avoid infinite loop
        } while (e != start);
    }

    /// Face traversal: given a directed edge (from → to via edge e),
    /// return the next directed edge around the face (CCW-next at `to`).
    /// Returns {next_edge, next_destination}.
    struct DirEdge {
        std::size_t edge_idx;
        std::size_t dest; ///< The vertex we're heading towards.
    };
    DirEdge face_next(std::size_t edge_idx, std::size_t dest) const;

    // ── Raw access for algorithms ───────────────────────────────

    std::vector<PlanarVertex>& vertices() { return vertices_; }
    std::vector<PlanarEdge>&   edges()    { return edges_; }
    const std::vector<PlanarVertex>& vertices() const { return vertices_; }
    const std::vector<PlanarEdge>&   edges()    const { return edges_; }

private:
    std::vector<PlanarVertex> vertices_;
    std::vector<PlanarEdge>   edges_;
};

} // namespace chazelle
