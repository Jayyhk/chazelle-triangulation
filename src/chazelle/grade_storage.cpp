#include "chazelle/grade_storage.h"

#include <cassert>
#include <stdexcept>

namespace chazelle {

void GradeStorage::init(std::size_t num_grades, std::size_t num_vertices) {
    storage_.resize(num_grades);
    // At grade λ, number of chains = (num_vertices - 1) / 2^λ.
    // num_vertices = 2^p + 1, so num_edges = 2^p.
    // At grade λ: 2^p / 2^λ = 2^(p-λ) chains.
    std::size_t num_edges = num_vertices - 1;
    for (std::size_t lambda = 0; lambda < num_grades; ++lambda) {
        std::size_t nc = num_edges >> lambda; // 2^(p - λ)
        storage_[lambda].resize(nc);
    }
}

void GradeStorage::store(std::size_t grade, std::size_t chain_index,
                          CanonicalSubmap submap) {
    assert(grade < storage_.size());
    assert(chain_index < storage_[grade].size());
    storage_[grade][chain_index] = std::move(submap);
}

const CanonicalSubmap& GradeStorage::get(std::size_t grade,
                                          std::size_t chain_index) const {
    assert(grade < storage_.size());
    assert(chain_index < storage_[grade].size());
    return storage_[grade][chain_index];
}

CanonicalSubmap& GradeStorage::get(std::size_t grade,
                                    std::size_t chain_index) {
    assert(grade < storage_.size());
    assert(chain_index < storage_[grade].size());
    return storage_[grade][chain_index];
}

std::size_t GradeStorage::num_chains(std::size_t grade) const {
    assert(grade < storage_.size());
    return storage_[grade].size();
}

} // namespace chazelle
