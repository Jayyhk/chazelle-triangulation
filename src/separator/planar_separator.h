#pragma once

#include "../common.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace chazelle {

class EmbeddedPlanarGraph {
public:
    EmbeddedPlanarGraph(std::size_t num_vertices,
                        std::vector<std::pair<std::size_t, std::size_t>> edges,
                        const std::vector<std::vector<std::size_t>>& rotations);

    std::size_t num_vertices() const noexcept {
        return num_vertices_;
    }
    std::size_t num_edges() const noexcept {
        return edges_.size();
    }

    std::size_t edge_u(std::size_t e) const noexcept {
        return edges_[e].first;
    }
    std::size_t edge_v(std::size_t e) const noexcept {
        return edges_[e].second;
    }

    std::size_t half_to(std::size_t h) const noexcept {
        return (h & 1) ? edges_[h >> 1].first : edges_[h >> 1].second;
    }
    std::size_t half_origin(std::size_t h) const noexcept {
        return half_to(h ^ 1);
    }

    std::size_t rot_next(std::size_t h) const noexcept {
        return nxt_[h];
    }
    std::size_t rot_prev(std::size_t h) const noexcept {
        return prv_[h];
    }

    const std::vector<std::size_t>& incident_halves(std::size_t v) const noexcept {
        return rot_[v];
    }

    EmbeddedPlanarGraph induced(const std::vector<std::size_t>& nodes,
                                std::vector<std::size_t>* old_index) const;

private:
    std::size_t num_vertices_ = 0;
    std::vector<std::pair<std::size_t, std::size_t>> edges_;
    std::vector<std::size_t> nxt_, prv_;
    std::vector<std::vector<std::size_t>> rot_;
};

enum class SepPart : std::uint8_t { A = 0, B = 1, D = 2 };

std::vector<SepPart> planar_separator(const EmbeddedPlanarGraph& g);

struct SeparatorDecomposition {
    std::vector<std::size_t> subset;
    std::size_t num_subsets = 0;
    std::size_t dstar_size = 0;
};

SeparatorDecomposition build_separator_decomposition(const EmbeddedPlanarGraph& g);

}
