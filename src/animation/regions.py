from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import pairwise

from trace_data import Viewport, exact_limit, signed_area


@dataclass
class RegionProjection:
    cells: list
    marker: tuple


def clip_polygon(points, viewport):
    for axis, boundary, sign in (
        (0, viewport.left, 1),
        (0, viewport.right, -1),
        (1, viewport.bottom, 1),
        (1, viewport.top, -1),
    ):
        clipped = []
        for a, b in zip(points, points[1:] + points[:1], strict=True):
            inside_a = sign * (a[axis] - boundary) >= 0
            inside_b = sign * (b[axis] - boundary) >= 0
            if inside_a != inside_b:
                t = (boundary - a[axis]) / (b[axis] - a[axis])
                clipped.append(tuple(a[i] + t * (b[i] - a[i]) for i in (0, 1)))
            if inside_b:
                clipped.append(b)
        points = clipped
    return points if len(points) >= 3 and signed_area(points) > 0 else []


def frame_position(viewport, point):
    x, y = point
    dx, dy = viewport.right - viewport.left, viewport.top - viewport.bottom
    if y == viewport.bottom:
        return x - viewport.left
    if x == viewport.right:
        return dx + y - viewport.bottom
    if y == viewport.top:
        return dx + dy + viewport.right - x
    if x == viewport.left:
        return 2 * dx + dy + viewport.top - y
    return None


def frame_point(viewport, position):
    dx, dy = viewport.right - viewport.left, viewport.top - viewport.bottom
    position %= 2 * (dx + dy)
    if position <= dx:
        return viewport.left + position, viewport.bottom
    if position <= dx + dy:
        return viewport.right, viewport.bottom + position - dx
    if position <= 2 * dx + dy:
        return viewport.right - position + dx + dy, viewport.top
    return viewport.left, viewport.top - position + 2 * dx + dy


def close_at_frame(submap, viewport, segments):
    ports = defaultdict(Counter)
    for region, boundary in segments.items():
        for a, b in boundary:
            for point, sign in ((a, -1), (b, 1)):
                position = frame_position(viewport, point)
                if position is not None:
                    ports[position][region] += sign
    ports = {
        position: {region: count for region, count in counts.items() if count}
        for position, counts in ports.items()
        if any(counts.values())
    }
    dx, dy = viewport.right - viewport.left, viewport.top - viewport.bottom
    perimeter = 2 * (dx + dy)
    corners = [Fraction(0), dx, dx + dy, 2 * dx + dy, perimeter]
    if not ports:
        caps = {
            arc["region"]
            for arc in submap.arcs
            for first, second in pairwise(arc["ranges"])
            if first[1] == second[0] and first[1] in (0, len(submap.boundary) - 1)
        }
        assert len(caps) == 1, (
            "[C91 section 2.1] Without a chord crossing the frame, both polar caps belong to the exterior region."
        )
        (region,) = caps
        segments[region].extend(
            (frame_point(viewport, a), frame_point(viewport, b)) for a, b in pairwise(corners)
        )
        return
    positions = sorted(ports)
    for index, start in enumerate(positions):
        end = positions[(index + 1) % len(positions)]
        assert sorted(ports[start].values()) == [-1, 1], (
            "Projected boundaries separate exactly two regions at a frame crossing."
        )
        (region,) = (region for region, count in ports[start].items() if count > 0)
        (before,) = (region for region, count in ports[end].items() if count < 0)
        assert region == before, "A frame interval stays in the same recorded region."
        if end <= start:
            end += perimeter
        stops = (
            [start]
            + [
                corner + offset
                for offset in (0, perimeter)
                for corner in corners[:-1]
                if start < corner + offset < end
            ]
            + [end]
        )
        segments[region].extend(
            (frame_point(viewport, a), frame_point(viewport, b)) for a, b in pairwise(stops)
        )


def horizontal_cells(segments):
    heights = sorted({point[1] for segment in segments for point in segment})
    for bottom, top in pairwise(heights):
        middle = (bottom + top) / 2
        crossings = defaultdict(list)
        for a, b in segments:
            if (a[1] > middle) == (b[1] > middle):
                continue

            slope = (b[0] - a[0]) / (b[1] - a[1])
            middle_x = a[0] + (middle - a[1]) * slope
            bottom_x = a[0] + (bottom - a[1]) * slope
            top_x = a[0] + (top - a[1]) * slope
            crossings[middle_x].append((1 if b[1] < a[1] else -1, bottom_x, top_x))
        positions = sorted(crossings)
        winding = 0
        for index, left in enumerate(positions):
            winding += sum(change for change, _, _ in crossings[left])
            assert winding in (0, 1), (
                "[C91 section 2.2] An oriented region boundary encloses its interior once."
            )
            if not winding or index + 1 == len(positions):
                continue
            right = positions[index + 1]
            _, left_bottom, left_top = crossings[left][0]
            _, right_bottom, right_top = crossings[right][0]
            assert left_bottom <= right_bottom and left_top <= right_top
            yield [(left_bottom, bottom), (right_bottom, bottom), (right_top, top), (left_top, top)]
        assert winding == 0, "A clipped region boundary closes."


def region_projections(submap, viewport):
    points = [*submap.boundary, (viewport.left, viewport.bottom), (viewport.right, viewport.top)]
    for arc in submap.arcs:
        for endpoint in (arc["start"], arc["end"]):
            limits = [exact_limit(coordinate).value for coordinate in endpoint]
            if None not in limits:
                points.append(tuple(limits))
    enclosing = Viewport(points, padding=0.08)
    segments = {region: [] for region in submap.regions}
    for arc in submap.arcs:
        segments[arc["region"]].extend(
            (a, b) for a, b in submap.arc_segments(arc, enclosing) if a != b
        )
    for chord in submap.chords:
        if chord["null"]:
            continue
        first = chord["left_region"]
        assert first in chord["regions"]
        (second,) = (region for region in chord["regions"] if region != first)
        for a, b in enclosing.world_chord_segments(chord):
            if a == b:
                continue
            if chord["left_direction"] == 0:
                a, b = b, a
            segments[first].append((a, b))
            segments[second].append((b, a))
    close_at_frame(submap, enclosing, segments)
    result = {}
    for region, boundary in segments.items():
        cells = [
            clipped
            for cell in horizontal_cells(boundary)
            if (clipped := clip_polygon(cell, viewport))
        ]
        if cells:
            candidates = []
            for cell in cells:
                bottom, top = min(y for _, y in cell), max(y for _, y in cell)
                middle = (bottom + top) / 2
                intersections = [
                    a[0] + (middle - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
                    for a, b in zip(cell, cell[1:] + cell[:1], strict=True)
                    if (a[1] > middle) != (b[1] > middle)
                ]
                left, right = sorted(intersections)
                candidates.append((min(right - left, top - bottom), ((left + right) / 2, middle)))
            _, marker = max(candidates)
        else:
            visible = [
                segment
                for arc in submap.region_arcs[region]
                for segment in submap.arc_segments(arc, viewport)
            ]
            if visible:
                a, b = max(
                    visible,
                    key=lambda segment: sum((b - a) ** 2 for a, b in zip(*segment, strict=True)),
                )
                marker = tuple((a + b) / 2 for a, b in zip(a, b, strict=True))
            else:
                marker = viewport.world_symbolic_point(submap.region_arcs[region][0]["start"])
                marker = (
                    min(viewport.right, max(viewport.left, marker[0])),
                    min(viewport.top, max(viewport.bottom, marker[1])),
                )
        result[region] = RegionProjection(cells, marker)
    return result
