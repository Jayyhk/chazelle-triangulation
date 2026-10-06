#pragma once

#include "merge/oracle.h"
#include "polygon/polygon.h"
#include "submap/submap.h"

#include <array>
#include <iosfwd>
#include <span>
#include <string_view>

namespace chazelle::animation {

class AnimationTrace {
public:
    AnimationTrace(std::ostream& output, std::span<const Point> vertices);
    ~AnimationTrace();
    AnimationTrace(const AnimationTrace&) = delete;
    AnimationTrace& operator=(const AnimationTrace&) = delete;

    static AnimationTrace* current() noexcept;

    void checkpoint(std::string_view name, std::size_t parameter = NONE);
    void boundary(const Polygon& curve);
    void chain(std::size_t grade, std::size_t index, const Polygon& curve);
    void submap(std::string_view name, const Polygon& curve, const Submap& submap,
                std::size_t granularity);
    void chord(std::string_view kind, const Polygon& curve, const Chord& chord);
    void merge(const Polygon& first, const Polygon& second, std::size_t granularity);
    std::size_t curve(const Polygon& curve);
    std::size_t map_id(const Submap& submap) const;
    void copy_submap(const Submap& source, Submap& copy, const Polygon& curve);
    void merge_inputs(const Submap& first, const Submap& second);
    void ray(const Polygon& curve, const Point& origin, Side direction, const RayHit& hit,
             std::size_t region = NONE);
    void cursor(const Polygon& curve, const Point& point);
    void invalidate(const Submap& submap, const Polygon& curve, std::size_t index);
    std::size_t build_begin(Submap& submap, const Polygon& curve);
    void build_walk(const Submap& submap, const Polygon& curve, const Arc& arc,
                    const SymbolicY& start, const SymbolicY& finish);
    void build_enter(const Submap& submap, std::size_t parent, std::size_t region,
                     std::size_t pending, const Polygon& curve, const Chord& chord);
    void build_leave(const Submap& submap, std::size_t region, std::size_t parent);
    void build_arc(const Submap& submap, const Polygon& curve, std::size_t index, const Arc& arc,
                   const SymbolicY& start, const SymbolicY& finish);
    void build_chord(const Submap& submap, const Polygon& curve, std::size_t index,
                     const Chord& chord);
    void build_end(const Submap& submap);
    void settled(std::string_view name, const Submap& submap);
    void remove(const Submap& submap, const Polygon& curve, std::size_t index, const Chord& chord);
    void split_arc(const Submap& submap, const Polygon& curve, std::size_t before,
                   std::size_t after, const SymbolicY& start, const SymbolicY& middle,
                   const SymbolicY& finish);
    void arc_owner(const Submap& submap, std::size_t arc, std::size_t region);
    void delete_arc(const Submap& submap, std::size_t arc);
    void insert(const Submap& submap, const Polygon& curve, std::size_t index, const Chord& chord);
    void reindex(const Submap& submap, std::span<const std::size_t> regions,
                 std::span<const std::size_t> chords, std::span<const std::size_t> arcs);
    void reindex_arcs(const Submap& submap, std::span<const std::size_t> arcs);
    void trapezoid(std::array<std::size_t, 4> vertices_and_edges);
    void diagonal(std::size_t first, std::size_t second);
    void piece(std::span<const std::size_t> vertices);
    void triangle(const std::array<std::size_t, 3>& vertices);
    void finish();

    std::size_t event_count() const noexcept {
        return events_;
    }

private:
    void begin(std::string_view kind);
    void end();
    void chord_geometry(const Polygon& curve, const Chord& chord);
    void chord_incidence(const Submap& submap, const Chord& chord);
    void arc_geometry(const Polygon& curve, const Arc& arc, SymbolicY start, SymbolicY finish);

    std::ostream& output_;
    std::span<const Point> vertices_;
    AnimationTrace* previous_;
    std::size_t events_ = 0;
    bool finished_ = false;
    std::size_t session_;
};

void animation_checkpoint(std::string_view name, std::size_t parameter = NONE);

}
