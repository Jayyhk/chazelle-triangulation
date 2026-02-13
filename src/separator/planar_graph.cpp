#include "separator/planar_graph.h"

#include <cassert>

namespace chazelle {

PlanarGraph::PlanarGraph(std::size_t nv, std::size_t ne) {
    vertices_.reserve(nv);
    edges_.reserve(ne);
}

std::size_t PlanarGraph::add_vertex(double cost) {
    std::size_t idx = vertices_.size();
    vertices_.push_back(PlanarVertex{.cost = cost});
    return idx;
}

std::size_t PlanarGraph::add_edge(std::size_t u, std::size_t v,
                                  std::size_t after_u,
                                  std::size_t after_v) {
    std::size_t idx = edges_.size();
    edges_.push_back(PlanarEdge{});
    auto& e = edges_[idx];
    e.endpoint[0] = u;
    e.endpoint[1] = v;

    splice_into_orbit(idx, u, after_u);
    splice_into_orbit(idx, v, after_v);

    return idx;
}

void PlanarGraph::splice_into_orbit(std::size_t e, std::size_t v,
                                    std::size_t after) {
    auto& vert = vertices_[v];
    int s = edges_[e].side_of(v);

    if (vert.some_edge == NONE) {
        // First edge at this vertex — point to self.
        vert.some_edge = e;
        edges_[e].cw[s]  = e;
        edges_[e].ccw[s] = e;
        return;
    }

    if (after == NONE) {
        // Insert immediately before the vertex's "some_edge" in CW order.
        after = vert.some_edge;
        // Find the edge just CW-before `some_edge` — i.e., the edge whose
        // CW-next is some_edge.  That's the CCW-prev of some_edge.
        int as = edges_[after].side_of(v);
        after = edges_[after].ccw[as]; // the edge just before some_edge in CW
    }

    // Insert `e` immediately CCW of `after` around vertex v.
    int as = edges_[after].side_of(v);
    std::size_t next = edges_[after].cw[as]; // was CW-next of after

    // Link: after → e → next
    edges_[after].cw[as] = e;
    edges_[e].ccw[s]     = after;
    edges_[e].cw[s]      = next;

    int ns = edges_[next].side_of(v);
    edges_[next].ccw[ns] = e;
}

PlanarGraph::DirEdge
PlanarGraph::face_next(std::size_t edge_idx, std::size_t dest) const {
    // We arrived at `dest` via `edge_idx`.
    // The next edge around the face is the CCW-next edge at `dest`,
    // and we travel along it away from `dest`.
    // Skip deleted edges to maintain correct face tracing.
    int s = edges_[edge_idx].side_of(dest);
    std::size_t nxt = edges_[edge_idx].ccw[s];
    // Walk past deleted edges.
    std::size_t guard = edges_.size();
    while (nxt != edge_idx && edges_[nxt].deleted && guard-- > 0) {
        int ns = edges_[nxt].side_of(dest);
        nxt = edges_[nxt].ccw[ns];
    }
    std::size_t new_dest = edges_[nxt].other(dest);
    return {nxt, new_dest};
}

} // namespace chazelle
