#include "trace.h"

#include "submap/boundary_geometry.h"

#include <cassert>
#include <ostream>
#include <stdexcept>

namespace chazelle::animation {
namespace {

thread_local AnimationTrace* active_trace = nullptr;
thread_local std::size_t next_session = 0;

void write_exact(std::ostream& output, const Exact& value) {
    output << '{';
    for (const bool denominator : {false, true}) {
        output << (denominator ? ",\"d\":[" : "\"n\":[");
        bool first = true;
        value.visit_terms(
            [&](bool is_denominator, const auto& powers, const Rational& coefficient) {
                if (is_denominator != denominator)
                    return;
                if (!first)
                    output << ',';
                first = false;
                output << "[[" << powers[0] << ',' << powers[1] << ',' << powers[2] << ','
                       << powers[3] << "],\"" << coefficient.to_string() << "\"]";
            });
        output << ']';
    }
    output << '}';
}

void write_point(std::ostream& output, const Point& point) {
    output << '[';
    write_exact(output, point.x);
    output << ',';
    write_exact(output, point.y);
    output << ']';
}

void write_indices(std::ostream& output, std::span<const std::size_t> indices) {
    output << '[';
    bool first = true;
    for (const std::size_t index : indices) {
        if (!first)
            output << ',';
        first = false;
        output << index;
    }
    output << ']';
}

}

AnimationTrace::AnimationTrace(std::ostream& output, std::span<const Point> vertices)
    : output_(output), vertices_(vertices), previous_(active_trace), session_(++next_session) {
    assert(session_ != 0);
    output_ << "{\"schema\":5,\"vertices\":[";
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        if (i != 0)
            output_ << ',';
        write_point(output_, vertices[i]);
    }
    output_ << "],\"events\":[";
    active_trace = this;
}

AnimationTrace::~AnimationTrace() {
    active_trace = previous_;
}
AnimationTrace* AnimationTrace::current() noexcept {
    return active_trace;
}

void AnimationTrace::begin(std::string_view kind) {
    assert(!finished_);
    if (events_ != 0)
        output_ << ',';
    output_ << "{\"seq\":" << events_++ << ",\"kind\":\"" << kind << '"';
}

void AnimationTrace::end() {
    output_ << '}';
}

void AnimationTrace::checkpoint(std::string_view name, std::size_t parameter) {
    begin("checkpoint");
    output_ << ",\"name\":\"" << name << '"';
    if (parameter != NONE)
        output_ << ",\"parameter\":" << parameter;
    end();
}

void animation_checkpoint(std::string_view name, std::size_t parameter) {
    if (auto* trace = AnimationTrace::current())
        trace->checkpoint(name, parameter);
}

void AnimationTrace::boundary(const Polygon& curve) {
    begin("boundary");
    output_ << ",\"first_tag\":" << curve.vertex(0).index << ",\"points\":[";
    for (std::size_t i = 0; i < curve.num_vertices(); ++i) {
        if (i != 0)
            output_ << ',';
        write_point(output_, curve.vertex(i));
    }
    output_ << ']';
    end();
}

void AnimationTrace::chain(std::size_t grade, std::size_t index, const Polygon& curve) {
    const std::size_t curve_index = this->curve(curve);
    begin("chain");
    output_ << ",\"grade\":" << grade << ",\"index\":" << index << ",\"curve\":" << curve_index
            << ",\"first\":" << curve.vertex(0).index
            << ",\"last\":" << curve.vertex(curve.num_vertices() - 1).index;
    end();
}

std::size_t AnimationTrace::curve(const Polygon& polygon) {
    if (polygon.animation_session_ == session_)
        return polygon.animation_curve_;
    auto table = [&](const Polygon::Table& data) {
        if (data.animation_session == session_)
            return;
        data.animation_session = session_;
        data.animation_table = events_;
        begin("coordinate_table");
        output_ << ",\"horizontal_shift\":";
        write_exact(output_, data.horizontal_shift);
        output_ << ",\"points\":[";
        for (std::size_t i = 0; i < data.vertices.size(); ++i) {
            if (i != 0)
                output_ << ',';
            write_point(output_, data.vertices[i]);
        }
        output_ << ']';
        end();
    };
    if (polygon.pieces_.empty())
        table(*polygon.table_);
    else
        for (const auto& piece : polygon.pieces_)
            table(*piece.table);
    polygon.animation_session_ = session_;
    polygon.animation_curve_ = events_;
    begin("curve");
    output_ << ",\"pieces\":[";
    if (polygon.pieces_.empty()) {
        output_ << '[' << polygon.table_->animation_table << ',' << polygon.offset_ << ','
                << polygon.len_ << ",false]";
    } else {
        for (std::size_t i = 0; i < polygon.pieces_.size(); ++i) {
            if (i != 0)
                output_ << ',';
            const auto& piece = polygon.pieces_[i];
            output_ << '[' << piece.table->animation_table << ',' << piece.offset << ','
                    << piece.count << ',' << (piece.reversed ? "true" : "false") << ']';
        }
    }
    output_ << ']';
    end();
    return polygon.animation_curve_;
}

std::size_t AnimationTrace::map_id(const Submap& submap) const {
    assert(submap.animation_session_ == session_ && submap.animation_map_ != NONE &&
           "A submap must be recorded at construction or copy before its later operations");
    return submap.animation_map_;
}

void AnimationTrace::copy_submap(const Submap& source, Submap& copy, const Polygon& polygon) {
    const auto curve_index = curve(polygon);
    const auto source_index = map_id(source);
    copy.animation_session_ = session_;
    copy.animation_map_ = events_;
    begin("copy_submap");
    output_ << ",\"source\":" << source_index << ",\"map\":" << map_id(copy)
            << ",\"curve\":" << curve_index;
    end();
}

void AnimationTrace::merge_inputs(const Submap& first, const Submap& second) {
    begin("merge_inputs");
    output_ << ",\"first\":" << map_id(first) << ",\"second\":" << map_id(second);
    end();
}

void AnimationTrace::ray(const Polygon& polygon, const Point& origin, Side direction,
                         const RayHit& hit, std::size_t region) {
    const auto curve_index = curve(polygon);
    begin("ray");
    output_ << ",\"curve\":" << curve_index << ",\"origin\":";
    write_point(output_, origin);
    output_ << ",\"direction\":" << static_cast<unsigned>(direction)
            << ",\"hit\":" << (hit.hit ? "true" : "false");
    if (region != NONE)
        output_ << ",\"region\":" << region;
    if (hit.hit) {
        output_ << ",\"contact\":";
        write_point(output_, {hit.x, hit.y, origin.index});
        output_ << ",\"edge\":" << hit.edge << ",\"wrapped\":" << (hit.wrapped ? "true" : "false");
    }
    end();
}

void AnimationTrace::cursor(const Polygon& polygon, const Point& point) {
    const auto curve_index = curve(polygon);
    begin("fusion_cursor");
    output_ << ",\"curve\":" << curve_index << ",\"point\":";
    write_point(output_, point);
    end();
}

void AnimationTrace::invalidate(const Submap& submap, const Polygon& polygon, std::size_t index) {
    const auto curve_index = curve(polygon);
    begin("fusion_remove");
    output_ << ",\"map\":" << map_id(submap) << ",\"curve\":" << curve_index << ",\"id\":" << index;
    chord_geometry(polygon, submap.chord(index));
    end();
}

void AnimationTrace::arc_geometry(const Polygon& polygon, const Arc& arc, SymbolicY start,
                                  SymbolicY finish) {
    output_ << ",\"region\":" << arc.region_node << ",\"edge_count\":" << arc.edge_count
            << ",\"traversal\":"
            << (arc.first_side == LEFT ? arc.first_edge
                                       : 2 * polygon.num_edges() - 1 - arc.first_edge)
            << ",\"parameter\":";
    const Exact parameter = edge_t_at_y(polygon, arc.first_edge, start);
    write_exact(output_, arc.first_side == LEFT ? parameter : 1 - parameter);
    output_ << ",\"ascending\":"
            << (shooting_direction(arc.first_edge, arc.first_side, polygon) == LEFT ? "true"
                                                                                    : "false")
            << ",\"start_tag\":" << start.tag << ",\"end_tag\":" << finish.tag << ",\"start\":";
    write_point(output_, {edge_x_at_y(polygon, arc.first_edge, start), start.y, start.tag});
    output_ << ",\"end\":";
    write_point(output_, {edge_x_at_y(polygon, arc.last_edge, finish), finish.y, finish.tag});
    output_ << ",\"ranges\":[";
    ArcSideRange ranges[3];
    const auto count = arc.side_ranges(0, polygon.num_vertices() - 1, ranges);
    for (std::size_t i = 0; i < count; ++i) {
        if (i != 0)
            output_ << ',';
        const bool forward = ranges[i].side == LEFT;
        output_ << '[' << (forward ? ranges[i].first_edge : ranges[i].last_edge + 1) << ','
                << (forward ? ranges[i].last_edge + 1 : ranges[i].first_edge) << ']';
    }
    output_ << "],\"wraps\":[";
    bool first_wrap = true;
    auto wrapped_edge = [&](std::size_t edge, const Exact& shift) {
        if (shift == 0)
            return;
        const Exact horizon = shift / 2;
        for (std::size_t i = 0; i < count; ++i) {
            const auto& range = ranges[i];
            if (edge < range.first_edge || edge > range.last_edge)
                continue;
            const bool forward = range.side == LEFT;
            const SymbolicY a = i == 0 && edge == arc.first_edge
                                    ? start
                                    : symbolic_y_of(polygon.vertex(forward ? edge : edge + 1));
            const SymbolicY b = i + 1 == count && edge == arc.last_edge
                                    ? finish
                                    : symbolic_y_of(polygon.vertex(forward ? edge + 1 : edge));
            const Exact x0 = polygon.vertex(edge).x;
            const Exact delta = polygon.edge_horizontal_delta(edge);
            const Exact first = x0 + edge_t_at_y(polygon, edge, a) * delta;
            const Exact last = x0 + edge_t_at_y(polygon, edge, b) * delta;
            if (!((first < horizon && horizon < last) || (last < horizon && horizon < first)))
                continue;
            if (!first_wrap)
                output_ << ',';
            output_ << '[' << edge << ',' << static_cast<unsigned>(range.side) << ']';
            first_wrap = false;
        }
    };
    if (polygon.pieces_.empty()) {
        wrapped_edge(0, polygon.table_->horizontal_shift);
    } else {
        std::size_t offset = 0;
        for (const auto& piece : polygon.pieces_) {
            wrapped_edge(offset, piece.reversed ? -piece.table->horizontal_shift
                                                : piece.table->horizontal_shift);
            offset += piece.count - 1;
        }
    }
    output_ << ']';
}

std::size_t AnimationTrace::build_begin(Submap& submap, const Polygon& polygon) {
    const auto curve_index = curve(polygon);
    submap.animation_session_ = session_;
    submap.animation_map_ = events_;
    begin("build_begin");
    output_ << ",\"map\":" << map_id(submap) << ",\"curve\":" << curve_index;
    end();
    return map_id(submap);
}

void AnimationTrace::build_walk(const Submap& submap, const Polygon& polygon, const Arc& arc,
                                const SymbolicY& start, const SymbolicY& finish) {
    begin("build_walk");
    output_ << ",\"map\":" << map_id(submap);
    arc_geometry(polygon, arc, start, finish);
    end();
}

void AnimationTrace::build_enter(const Submap& submap, std::size_t parent, std::size_t region,
                                 std::size_t pending, const Polygon& polygon, const Chord& chord) {
    begin("build_enter");
    output_ << ",\"map\":" << map_id(submap) << ",\"parent\":" << parent << ",\"region\":" << region
            << ",\"pending\":" << pending;
    chord_geometry(polygon, chord);
    end();
}

void AnimationTrace::build_leave(const Submap& submap, std::size_t region, std::size_t parent) {
    begin("build_leave");
    output_ << ",\"map\":" << map_id(submap) << ",\"region\":" << region
            << ",\"parent\":" << parent;
    end();
}

void AnimationTrace::build_arc(const Submap& submap, const Polygon& polygon, std::size_t index,
                               const Arc& arc, const SymbolicY& start, const SymbolicY& finish) {
    begin("build_arc");
    output_ << ",\"map\":" << map_id(submap) << ",\"id\":" << index;
    arc_geometry(polygon, arc, start, finish);
    end();
}

void AnimationTrace::build_chord(const Submap& submap, const Polygon& polygon, std::size_t index,
                                 const Chord& chord) {
    begin("build_chord");
    output_ << ",\"map\":" << map_id(submap) << ",\"id\":" << index << ",\"regions\": ["
            << chord.region[0] << ',' << chord.region[1] << ']';
    chord_geometry(polygon, chord);
    chord_incidence(submap, chord);
    end();
}

void AnimationTrace::build_end(const Submap& submap) {
    begin("build_end");
    output_ << ",\"map\":" << map_id(submap)
            << ",\"root\":" << submap.arc(submap.start_arc).region_node;
    end();
}

void AnimationTrace::settled(std::string_view name, const Submap& submap) {
    begin("checkpoint");
    output_ << ",\"name\":\"" << name << "\",\"map\":" << map_id(submap)
            << ",\"root\":" << submap.arc(submap.start_arc).region_node;
    end();
}

void AnimationTrace::remove(const Submap& submap, const Polygon& polygon, std::size_t index,
                            const Chord& chord) {
    const auto curve_index = curve(polygon);
    begin("remove_chord");
    output_ << ",\"map\":" << map_id(submap) << ",\"curve\":" << curve_index << ",\"id\":" << index
            << ",\"regions\": [" << chord.region[0] << ',' << chord.region[1] << ']';
    chord_geometry(polygon, chord);
    end();
}

void AnimationTrace::split_arc(const Submap& submap, const Polygon& polygon, std::size_t before,
                               std::size_t after, const SymbolicY& start, const SymbolicY& middle,
                               const SymbolicY& finish) {
    build_arc(submap, polygon, before, submap.arc(before), start, middle);
    build_arc(submap, polygon, after, submap.arc(after), middle, finish);
}

void AnimationTrace::arc_owner(const Submap& submap, std::size_t arc, std::size_t region) {
    begin("arc_owner");
    output_ << ",\"map\":" << map_id(submap) << ",\"arc\":" << arc << ",\"region\":" << region;
    end();
}

void AnimationTrace::delete_arc(const Submap& submap, std::size_t arc) {
    begin("delete_arc");
    output_ << ",\"map\":" << map_id(submap) << ",\"arc\":" << arc;
    end();
}

void AnimationTrace::insert(const Submap& submap, const Polygon& polygon, std::size_t index,
                            const Chord& chord) {
    begin("insert_chord");
    output_ << ",\"map\":" << map_id(submap) << ",\"id\":" << index << ",\"regions\": ["
            << chord.region[0] << ',' << chord.region[1] << ']';
    chord_geometry(polygon, chord);
    chord_incidence(submap, chord);
    end();
}

void AnimationTrace::reindex(const Submap& submap, std::span<const std::size_t> regions,
                             std::span<const std::size_t> chords,
                             std::span<const std::size_t> arcs) {
    begin("reindex");
    output_ << ",\"map\":" << map_id(submap) << ",\"regions\":";
    write_indices(output_, regions);
    output_ << ",\"chords\":";
    write_indices(output_, chords);
    output_ << ",\"arcs\":";
    write_indices(output_, arcs);
    end();
}

void AnimationTrace::reindex_arcs(const Submap& submap, std::span<const std::size_t> arcs) {
    begin("reindex_arcs");
    output_ << ",\"map\":" << map_id(submap) << ",\"arcs\":";
    write_indices(output_, arcs);
    end();
}

void AnimationTrace::chord_geometry(const Polygon& curve, const Chord& chord) {
    output_ << ",\"left\":";
    write_exact(output_, edge_x_at_y(curve, chord.left_edge, chord.symbolic_y()));
    output_ << ",\"right\":";
    write_exact(output_, edge_x_at_y(curve, chord.right_edge, chord.symbolic_y()));
    output_ << ",\"y\":";
    write_exact(output_, chord.y);
    output_ << ",\"tag\":" << chord.y_tag << ",\"left_edge\":" << chord.left_edge
            << ",\"right_edge\":" << chord.right_edge
            << ",\"left_side\":" << static_cast<unsigned>(chord.left_side)
            << ",\"right_side\":" << static_cast<unsigned>(chord.right_side)
            << ",\"left_direction\":"
            << static_cast<unsigned>(shooting_direction(chord.left_edge, chord.left_side, curve))
            << ",\"null\":" << (chord.is_null_length ? "true" : "false") << ",\"infinite\":"
            << (!chord.is_null_length && chord_runs_through_infinity(curve, chord) ? "true"
                                                                                   : "false");
}

void AnimationTrace::chord_incidence(const Submap& submap, const Chord& chord) {
    assert(chord.left_adj.count > 0 &&
           "[C91 section 2.2]: every exit chord has an incoming boundary arc");
    const auto region = submap.arc(chord.left_adj.arcs[0]).region_node;
    assert((region == chord.region[0] || region == chord.region[1]) &&
           "[C91 section 2.2]: the incoming arc bounds an incident region");
    output_ << ",\"left_region\":" << region;
}

void AnimationTrace::chord(std::string_view kind, const Polygon& curve, const Chord& chord) {
    const std::size_t curve_index = this->curve(curve);
    begin(kind);
    output_ << ",\"curve\":" << curve_index;
    chord_geometry(curve, chord);
    end();
}

void AnimationTrace::submap(std::string_view name, const Polygon& curve, const Submap& submap,
                            std::size_t granularity) {
    assert(submap.is_conformal() &&
           "[C91 §§4.1–4.2]: checkpoint submaps are conformal, so arc endpoint queries take O(1)");
    const std::size_t curve_index = this->curve(curve);
    begin("map_begin");
    output_ << ",\"name\":\"" << name << "\",\"granularity\":" << granularity
            << ",\"map\":" << map_id(submap) << ",\"curve\":" << curve_index
            << ",\"first\":0,\"last\":" << curve.num_vertices() - 1
            << ",\"root\":" << submap.arc(submap.start_arc).region_node;
    end();
    for (std::size_t i = 0; i < submap.num_nodes(); ++i) {
        if (submap.node(i).dead)
            continue;
        bool bounded = !submap.node(i).incident_chords.empty();
        for (const std::size_t chord_index : submap.node(i).incident_chords) {
            const Chord& chord = submap.chord(chord_index);
            if (!chord.is_null_length && chord_runs_through_infinity(curve, chord))
                bounded = false;
        }
        begin("map_region");
        output_ << ",\"id\":" << i << ",\"weight\":" << submap.region_weight(i)
                << ",\"bounded\":" << (bounded ? "true" : "false");
        end();
    }
    for (std::size_t i = 0; i < submap.num_chords(); ++i) {
        const Chord& c = submap.chord(i);
        if (c.dead)
            continue;
        begin("map_chord");
        output_ << ",\"id\":" << i << ",\"regions\": [" << c.region[0] << ',' << c.region[1] << ']';
        chord_geometry(curve, c);
        chord_incidence(submap, c);
        end();
    }
    for (std::size_t i = 0; i < submap.num_arcs(); ++i) {
        const Arc& arc = submap.arc(i);
        if (arc.dead)
            continue;
        const SymbolicY start = submap.arc_start_symbolic_y(i, curve);
        const SymbolicY finish = submap.arc_end_symbolic_y(i, curve);
        begin("map_arc");
        output_ << ",\"id\":" << i;
        arc_geometry(curve, arc, start, finish);
        end();
    }
    begin("map_end");
    end();
}

void AnimationTrace::merge(const Polygon& first, const Polygon& second, std::size_t granularity) {
    const auto first_index = curve(first);
    const auto second_index = curve(second);
    const auto combined_index = curve(Polygon(first, second));
    begin("merge");
    output_ << ",\"granularity\":" << granularity << ",\"first_curve\":" << first_index
            << ",\"second_curve\":" << second_index << ",\"curve\":" << combined_index
            << ",\"edges\":" << first.num_edges() + second.num_edges() << ",\"points\": [";
    write_point(output_, first.vertex(0));
    output_ << ',';
    write_point(output_, second.vertex(0));
    output_ << ',';
    write_point(output_, second.vertex(second.num_vertices() - 1));
    output_ << ']';
    end();
}

void AnimationTrace::trapezoid(std::array<std::size_t, 4> vertices_and_edges) {
    begin("trapezoid");
    output_ << ",\"vertices_and_edges\":";
    write_indices(output_, vertices_and_edges);
    output_ << ",\"corners\":[";
    const auto [top, bottom, left, right] = vertices_and_edges;
    const std::array<std::array<std::size_t, 2>, 4> contacts{
        {{{left, top}}, {{right, top}}, {{right, bottom}}, {{left, bottom}}}};
    for (std::size_t i = 0; i < contacts.size(); ++i) {
        if (i != 0)
            output_ << ',';
        const auto [edge, vertex] = contacts[i];
        const Point& a = vertices_[edge];
        const Point& b = vertices_[(edge + 1) % vertices_.size()];
        const Point& level = vertices_[vertex];
        Exact x;
        if (vertex == edge)
            x = a.x;
        else if (vertex == (edge + 1) % vertices_.size())
            x = b.x;
        else {
            assert(a.y != b.y && "[FM84 Algorithm 1]: a trapezoid side meets its endpoint levels");
            x = a.x + (level.y - a.y) * (b.x - a.x) / (b.y - a.y);
        }
        write_point(output_, {x, level.y, vertex});
    }
    output_ << ']';
    end();
}

void AnimationTrace::diagonal(std::size_t first, std::size_t second) {
    begin("diagonal");
    output_ << ",\"vertices\": [" << first << ',' << second << ']';
    end();
}

void AnimationTrace::piece(std::span<const std::size_t> vertices) {
    begin("piece");
    output_ << ",\"vertices\":";
    write_indices(output_, vertices);
    end();
}

void AnimationTrace::triangle(const std::array<std::size_t, 3>& vertices) {
    begin("triangle");
    output_ << ",\"vertices\":";
    write_indices(output_, vertices);
    end();
}

void AnimationTrace::finish() {
    assert(!finished_);
    output_ << "]}\n";
    output_.flush();
    if (!output_)
        throw std::runtime_error("Failed to write the animation trace.");
    finished_ = true;
}

}
