#include "polygon.h"

#include <algorithm>
#include <cassert>

namespace chazelle {

Polygon::Polygon(std::vector<Point> vertices) {
    initialize(std::move(vertices), true);
}

Polygon Polygon::symbolic_curve(std::vector<Point> vertices, Exact horizontal_shift) {
    Polygon curve;
    curve.initialize(std::move(vertices), false, std::move(horizontal_shift));
    return curve;
}

void Polygon::initialize(std::vector<Point> vertices, bool original, Exact horizontal_shift) {
    assert(vertices.size() >= 2 && "[C91 §2.1]: curve needs at least two vertices");
    assert(vertices[0].index != SOS_NONE && "[C91 §2]: a vertex requires a symbolic tag");
    assert(horizontal_shift == 0 || vertices.size() == 2);
    if (horizontal_shift != 0) {
        const Exact period = horizontal_shift < 0 ? -horizontal_shift : horizontal_shift;
        for ([[maybe_unused]] const Point& point : vertices)
            assert(
                -period / 2 < point.x && point.x < period / 2 &&
                "[C91 §2 tex 59, §4.2 tex 371]: wrapping edges use the horizontal identification");
    }
    auto table = std::make_shared<Table>();
    table->horizontal_shift = std::move(horizontal_shift);
    table->tag_stride =
        vertices[1].index > vertices[0].index ? vertices[1].index - vertices[0].index : 0;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        assert(vertices[i].index != SOS_NONE);
        if (original) {
            assert(i < SOS_NONE - vertices[0].index && vertices[i].index == vertices[0].index + i &&
                   "[C91 §2.4 tex 133]: original tags follow input table order");
        }
        if (table->tag_stride != 0 &&
            (i > (SOS_NONE - 1 - vertices[0].index) / table->tag_stride ||
             vertices[i].index != vertices[0].index + i * table->tag_stride))
            table->tag_stride = 0;
    }
    assert(table->tag_stride != 0 || vertices.size() <= 32);
    table->vertices = std::move(vertices);
    const std::size_t count = table->vertices.size() - 1;
    table->nonnull_prefix.resize(count + 1, 0);
    for (std::size_t i = 0; i < count; ++i) {
        const Point& first = table->vertices[i];
        const Point& last = table->vertices[i + 1];
        table->nonnull_prefix[i + 1] =
            table->nonnull_prefix[i] +
            ((first.x != last.x || first.y != last.y || table->horizontal_shift != 0) ? 1 : 0);
    }
    table->next_nonnull.resize(count, count);
    table->previous_nonnull.resize(count, NONE);
    for (std::size_t i = count; i-- > 0;)
        table->next_nonnull[i] = table->nonnull_prefix[i + 1] > table->nonnull_prefix[i]
                                     ? i
                                     : (i + 1 < count ? table->next_nonnull[i + 1] : count);
    for (std::size_t i = 0; i < count; ++i)
        table->previous_nonnull[i] = table->nonnull_prefix[i + 1] > table->nonnull_prefix[i]
                                         ? i
                                         : (i > 0 ? table->previous_nonnull[i - 1] : NONE);
    table_ = std::move(table);
    len_ = table_->vertices.size();
    find_y_extremes();
}

void Polygon::append_pieces(const Polygon& curve) {
    if (curve.pieces_.empty())
        pieces_.push_back({curve.table_, curve.offset_, curve.len_, false});
    else
        pieces_.insert(pieces_.end(), curve.pieces_.begin(), curve.pieces_.end());
    assert(pieces_.size() <= 64 &&
           "[C91 §4.2 tex 367–372]: a region has a constant number of boundary pieces");
}

Polygon::Polygon(const Polygon& first, const Polygon& second) {
    [[maybe_unused]] const Point& junction = first.vertex(first.num_vertices() - 1);
    [[maybe_unused]] const Point& next = second.vertex(0);
    assert(junction.index == next.index && junction.x == next.x && junction.y == next.y &&
           "[C91 §3 tex 160]: curves share their common endpoint");
    if (first.pieces_.empty() && second.pieces_.empty() && first.table_ == second.table_ &&
        second.offset_ == first.offset_ + first.len_ - 1) {
        table_ = first.table_;
        offset_ = first.offset_;
    } else {
        append_pieces(first);
        append_pieces(second);
    }
    len_ = first.len_ + second.len_ - 1;
    const std::size_t shift = first.len_ - 1;
    max_y_vertex_ =
        point_y_above(second.vertex(second.max_y_vertex_), first.vertex(first.max_y_vertex_))
            ? second.max_y_vertex_ + shift
            : first.max_y_vertex_;
    min_y_vertex_ =
        point_y_below(second.vertex(second.min_y_vertex_), first.vertex(first.min_y_vertex_))
            ? second.min_y_vertex_ + shift
            : first.min_y_vertex_;
}

Polygon Polygon::subchain(std::size_t first, std::size_t count, std::size_t minimum,
                          std::size_t maximum) const {
    assert(count >= 2 && first + count <= len_ && minimum < count && maximum < count);
    Polygon result;
    result.len_ = count;
    result.min_y_vertex_ = minimum;
    result.max_y_vertex_ = maximum;
    if (pieces_.empty()) {
        result.table_ = table_;
        result.offset_ = offset_ + first;
        return result;
    }
    std::size_t remaining = count - 1;
    for (const Piece& piece : pieces_) {
        const std::size_t edges = piece.count - 1;
        if (first >= edges) {
            first -= edges;
            continue;
        }
        const std::size_t taken = std::min(remaining, edges - first);
        const std::size_t offset =
            piece.reversed ? piece.offset + edges - first - taken : piece.offset + first;
        result.pieces_.push_back({piece.table, offset, taken + 1, piece.reversed});
        remaining -= taken;
        first = 0;
        if (remaining == 0)
            break;
    }
    assert(remaining == 0);
    return result;
}

Polygon Polygon::subchain(std::size_t first, std::size_t count) const {
    Polygon result = subchain(first, count, 0, 0);
    result.find_y_extremes();
    return result;
}

Polygon Polygon::reversed() const {
    Polygon result;
    result.append_pieces(*this);
    std::reverse(result.pieces_.begin(), result.pieces_.end());
    for (Piece& piece : result.pieces_)
        piece.reversed = !piece.reversed;
    result.len_ = len_;
    result.min_y_vertex_ = len_ - 1 - min_y_vertex_;
    result.max_y_vertex_ = len_ - 1 - max_y_vertex_;
    return result;
}

Exact Polygon::edge_horizontal_shift(std::size_t edge) const {
    assert(edge < num_edges());
    if (pieces_.empty())
        return table_->horizontal_shift;
    for (const Piece& piece : pieces_) {
        if (edge < piece.count - 1)
            return piece.reversed ? -piece.table->horizontal_shift : piece.table->horizontal_shift;
        edge -= piece.count - 1;
    }
    assert(false);
    return 0;
}

bool Polygon::previous_branch_left(std::size_t vertex_index) const {
    assert(!is_endpoint(vertex_index) && is_y_extremum(vertex_index));
    Point previous = vertex(vertex_index - 1);
    Point next = vertex(vertex_index + 1);
    previous.x -= edge_horizontal_shift(vertex_index - 1);
    next.x += edge_horizontal_shift(vertex_index);
    return extremum_prev_branch_left(previous, vertex(vertex_index), next);
}

std::size_t Polygon::local_index_of_tag(std::size_t tag) const noexcept {
    auto find = [&](const Piece& piece) {
        if (piece.table->tag_stride != 0) {
            const std::size_t first = piece.table->vertices[piece.offset].index;
            if (tag < first || (tag - first) % piece.table->tag_stride != 0)
                return NONE;
            const std::size_t position = (tag - first) / piece.table->tag_stride;
            if (position >= piece.count)
                return NONE;
            return piece.reversed ? piece.count - 1 - position : position;
        }
        for (std::size_t i = 0; i < piece.count; ++i)
            if (piece.table->vertices[piece.offset + i].index == tag)
                return piece.reversed ? piece.count - 1 - i : i;
        return NONE;
    };
    if (pieces_.empty())
        return find({table_, offset_, len_, false});
    std::size_t offset = 0;
    for (const Piece& piece : pieces_) {
        const std::size_t position = find(piece);
        if (position != NONE)
            return offset + position;
        offset += piece.count - 1;
    }
    return NONE;
}

std::size_t Polygon::count_nonnull_edges(std::size_t lo, std::size_t hi) const noexcept {
    assert(lo <= hi && hi < num_edges());
    if (pieces_.empty())
        return table_->nonnull_prefix[offset_ + hi + 1] - table_->nonnull_prefix[offset_ + lo];
    std::size_t count = 0;
    std::size_t offset = 0;
    for (const Piece& piece : pieces_) {
        const std::size_t edges = piece.count - 1;
        if (lo < offset + edges && hi >= offset) {
            const std::size_t first = std::max(lo, offset) - offset;
            const std::size_t last = std::min(hi, offset + edges - 1) - offset;
            const std::size_t begin = piece.offset + (piece.reversed ? edges - 1 - last : first);
            const std::size_t end = piece.offset + (piece.reversed ? edges - first : last + 1);
            count += piece.table->nonnull_prefix[end] - piece.table->nonnull_prefix[begin];
        }
        offset += edges;
    }
    return count;
}

std::size_t Polygon::next_nonnull_edge(std::size_t first) const noexcept {
    assert(first < num_edges());
    if (pieces_.empty()) {
        const std::size_t position = table_->next_nonnull[offset_ + first];
        return position >= offset_ + num_edges() ? num_edges() : position - offset_;
    }
    std::size_t offset = 0;
    for (const Piece& piece : pieces_) {
        const std::size_t edges = piece.count - 1;
        if (first < offset + edges) {
            const std::size_t local = first > offset ? first - offset : 0;
            if (piece.reversed) {
                const std::size_t position =
                    piece.table->previous_nonnull[piece.offset + edges - 1 - local];
                if (position != NONE && position >= piece.offset)
                    return offset + edges - 1 - (position - piece.offset);
            } else {
                const std::size_t position = piece.table->next_nonnull[piece.offset + local];
                if (position < piece.offset + edges)
                    return offset + position - piece.offset;
            }
        }
        offset += edges;
    }
    return num_edges();
}

void Polygon::find_y_extremes() {
    max_y_vertex_ = 0;
    min_y_vertex_ = 0;
    for (std::size_t i = 1; i < len_; ++i) {
        if (point_y_above(vertex(i), vertex(max_y_vertex_)))
            max_y_vertex_ = i;
        if (point_y_below(vertex(i), vertex(min_y_vertex_)))
            min_y_vertex_ = i;
    }
}

}
