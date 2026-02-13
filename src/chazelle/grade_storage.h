#pragma once

/// Grade storage — per-grade canonical submaps + ray-shooting structures.
///
/// Organization (§4.1):
///   grade_storage[λ] has (n−1)/2^λ entries.
///   Each entry stores the canonical submap for one chain at grade λ,
///   plus its preprocessed ray-shooting oracle.
///
/// Total storage is O(n) by geometric series.

#include "visibility/submap.h"
#include "visibility/tree_decomposition.h"
#include "oracles/ray_shooting.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace chazelle {

/// Stored data for one chain's canonical submap at a given grade.
struct CanonicalSubmap {
    Submap submap;
    RayShootingOracle oracle;
    TreeDecomposition tree_decomp;  ///< Cached tree decomposition (§2.3).
    bool td_built = false;          ///< Whether tree_decomp has been built.

    /// Chain metadata.
    std::size_t grade = 0;
    std::size_t chain_index = 0;
    std::size_t start_vertex = 0;
    std::size_t end_vertex = 0;
};

/// All grade storage: grade_storage[λ][chain_idx].
class GradeStorage {
public:
    GradeStorage() = default;

    /// Initialize for a polygon with p grades.
    void init(std::size_t num_grades, std::size_t num_vertices);

    /// Store a canonical submap for a chain.
    void store(std::size_t grade, std::size_t chain_index,
               CanonicalSubmap submap);

    /// Retrieve a canonical submap for a chain.
    const CanonicalSubmap& get(std::size_t grade,
                                std::size_t chain_index) const;
    CanonicalSubmap& get(std::size_t grade,
                          std::size_t chain_index);

    /// Number of chains at a given grade.
    std::size_t num_chains(std::size_t grade) const;

    /// Total number of grades.
    std::size_t num_grades() const { return storage_.size(); }

private:
    /// storage_[grade][chain_index]
    std::vector<std::vector<CanonicalSubmap>> storage_;
};

} // namespace chazelle
