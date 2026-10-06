import json
import sys
from dataclasses import dataclass
from fractions import Fraction
from functools import total_ordering
from itertools import pairwise
from pathlib import Path

sys.set_int_max_str_digits(0)


@dataclass(frozen=True)
class Limit:
    value: Fraction | None
    infinity: int = 0


def exact_limit(value):
    numerator = [(tuple(powers), Fraction(coefficient)) for powers, coefficient in value["n"]]
    denominator = [(tuple(powers), Fraction(coefficient)) for powers, coefficient in value["d"]]
    numerator = [(powers, coefficient) for powers, coefficient in numerator if coefficient]
    denominator = [(powers, coefficient) for powers, coefficient in denominator if coefficient]
    if not denominator:
        raise ValueError("An exact coordinate has a zero denominator.")
    if not numerator:
        return Limit(Fraction(0))
    n_powers, n_coefficient = min(numerator)
    d_powers, d_coefficient = min(denominator)
    powers = tuple(n - d for n, d in zip(n_powers, d_powers, strict=True))
    coefficient = n_coefficient / d_coefficient
    if powers > (0, 0, 0, 0):
        return Limit(Fraction(0))
    if powers < (0, 0, 0, 0):
        return Limit(None, 1 if coefficient > 0 else -1)
    return Limit(coefficient)


def finite_point(point):
    result = tuple(exact_limit(coordinate).value for coordinate in point)
    if len(result) != 2 or None in result:
        raise ValueError("Polygon vertices must have finite coordinates.")
    return result


@total_ordering
class ExactCoordinate:
    def __init__(self, value):
        self.numerator = self.terms(value["n"])
        self.denominator = self.terms(value["d"])
        if not self.denominator:
            raise ValueError("An exact coordinate has a zero denominator.")

    @staticmethod
    def terms(terms):
        result = {}
        for powers, coefficient in terms:
            key = tuple(powers)
            result[key] = result.get(key, Fraction(0)) + Fraction(coefficient)
        return {key: coefficient for key, coefficient in result.items() if coefficient}

    def compare(self, other):
        difference = {}
        for sign, numerator, denominator in (
            (1, self.numerator, other.denominator),
            (-1, other.numerator, self.denominator),
        ):
            for a, x in numerator.items():
                for b, y in denominator.items():
                    key = tuple(i + j for i, j in zip(a, b, strict=True))
                    difference[key] = difference.get(key, Fraction(0)) + sign * x * y
        difference = {key: coefficient for key, coefficient in difference.items() if coefficient}
        if not difference:
            return 0
        sign = (
            difference[min(difference)]
            * self.denominator[min(self.denominator)]
            * other.denominator[min(other.denominator)]
        )
        return 1 if sign > 0 else -1

    def __eq__(self, other):
        if not isinstance(other, ExactCoordinate):
            return NotImplemented
        return self.compare(other) == 0

    def __lt__(self, other):
        if not isinstance(other, ExactCoordinate):
            return NotImplemented
        return self.compare(other) < 0


def signed_area(points):
    return sum(
        a[0] * b[1] - a[1] * b[0] for a, b in zip(points, points[1:] + points[:1], strict=True)
    )


def polygon_interior_point(points):
    heights = sorted({y for _, y in points})
    intervals = []
    for bottom, top in pairwise(heights):
        y = (bottom + top) / 2
        intersections = sorted(
            a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            for a, b in zip(points, points[1:] + points[:1], strict=True)
            if (a[1] > y) != (b[1] > y)
        )
        assert len(intersections) % 2 == 0
        for left, right in zip(intersections[::2], intersections[1::2], strict=True):
            if left < right:
                intervals.append((right - left, (left + right) / 2, y))
    assert intervals, "A region with positive area has a nonempty horizontal interior interval."
    _, x, y = max(intervals)
    return x, y


class Submap:
    def __init__(self, metadata, events, boundary, first_tag, conformal=True):
        self.metadata = metadata
        self.regions = {event["id"]: event for event in events if event["kind"] == "map_region"}
        self.chords = [event for event in events if event["kind"] == "map_chord"]
        self.arcs = [event for event in events if event["kind"] == "map_arc"]
        self.boundary = boundary
        self.first_tag = first_tag
        self.adjacency = {region: [] for region in self.regions}
        self.region_arcs = {region: [] for region in self.regions}
        for chord in self.chords:
            a, b = chord["regions"]
            if a == b or a not in self.regions or b not in self.regions:
                raise ValueError("A recorded chord must separate two recorded regions.")
            self.adjacency[a].append(b)
            self.adjacency[b].append(a)
        for arc in self.arcs:
            self.region_arcs[arc["region"]].append(arc)
            for a, b in arc["ranges"]:
                if min(a, b) < metadata["first"] or max(a, b) > metadata["last"]:
                    raise ValueError("An arc range leaves its recorded chain.")
        for arcs in self.region_arcs.values():
            arcs.sort(key=self.boundary_position)
        if len(self.chords) + 1 != len(self.regions):
            raise ValueError("[C91 section 2.2] The submap's dual graph must be a tree.")
        if conformal and any(len(neighbors) > 4 for neighbors in self.adjacency.values()):
            raise ValueError("[C91 section 2.3] A conformal region has at most four exit chords.")
        if any(not arcs or (conformal and len(arcs) > 4) for arcs in self.region_arcs.values()):
            raise ValueError("[C91 section 2.3] Each conformal region has one to four arcs.")
        self.root = metadata["root"]
        self.parent = {self.root: None}
        self.order = []
        pending = [self.root]
        while pending:
            node = pending.pop()
            self.order.append(node)
            for neighbor in self.adjacency[node]:
                if neighbor == self.parent[node]:
                    continue
                if neighbor in self.parent:
                    raise ValueError("The recorded dual graph contains a cycle.")
                self.parent[neighbor] = node
                pending.append(neighbor)
        if len(self.order) != len(self.regions):
            raise ValueError("The recorded dual graph is disconnected.")

    def boundary_position(self, arc):
        return (
            arc["traversal"],
            ExactCoordinate(arc["parameter"]),
            -arc["start_tag"] if arc["ascending"] else arc["start_tag"],
            arc["edge_count"] != 0,
        )

    def arc_points(self, arc):
        points = []
        for a, b in arc["ranges"]:
            step = 1 if a < b else -1
            points.extend(self.boundary[tag - self.first_tag] for tag in range(a, b + step, step))
        points[0] = finite_point(arc["start"])
        points[-1] = finite_point(arc["end"])
        return points

    def region_points(self, region):
        return [point for arc in self.region_arcs[region] for point in self.arc_points(arc)]

    def single_edge_region_order(self):
        assert len(self.boundary) == 2 and not self.boundary.edge_wrap(0)
        assert len(self.chords) == 2 and len(self.regions) == 3
        assert all(chord["infinite"] and not chord["null"] for chord in self.chords)
        low, high = sorted(
            self.chords, key=lambda chord: (ExactCoordinate(chord["y"]), -chord["tag"])
        )
        middle = set(low["regions"]) & set(high["regions"])
        assert len(middle) == 1, (
            "[C91 section 2.1] A single edge has three regions joined by its two endpoint chords."
        )
        middle = middle.pop()
        below = next(region for region in low["regions"] if region != middle)
        above = next(region for region in high["regions"] if region != middle)
        return below, middle, above

    def single_edge_regions(self, viewport):
        below, middle, above = self.single_edge_region_order()
        low, high = sorted(
            self.chords, key=lambda chord: (ExactCoordinate(chord["y"]), -chord["tag"])
        )
        heights = (
            viewport.bottom,
            exact_limit(low["y"]).value,
            exact_limit(high["y"]).value,
            viewport.top,
        )
        result = {}
        for region, bottom, top in zip(
            (below, middle, above), heights[:-1], heights[1:], strict=True
        ):
            bottom, top = max(bottom, viewport.bottom), min(top, viewport.top)
            if bottom < top:
                result[region] = [
                    (viewport.left, bottom),
                    (viewport.right, bottom),
                    (viewport.right, top),
                    (viewport.left, top),
                ]
        return result

    def arc_segments(self, arc, viewport):
        result = []
        for index, (a, b) in enumerate(arc["ranges"]):
            step = 1 if a < b else -1
            side = 0 if step > 0 else 1
            for vertex in range(a, b, step):
                edge = min(vertex, vertex + step)
                start = (
                    viewport.world_symbolic_point(arc["start"])
                    if index == 0 and vertex == a
                    else self.boundary[vertex - self.first_tag]
                )
                end = (
                    viewport.world_symbolic_point(arc["end"])
                    if index + 1 == len(arc["ranges"]) and vertex + step == b
                    else self.boundary[vertex + step - self.first_tag]
                )
                wrap = self.boundary.edge_wrap(edge) * step if [edge, side] in arc["wraps"] else 0
                result.extend(viewport.edge_segments(start, end, wrap))
        return result

    def arc_paths(self, arc, viewport):
        result = []
        append_segments(result, self.arc_segments(arc, viewport))
        return [[viewport.point(point) for point in points] for points in result]

    def tree_positions(self, width=3.25, height=1.55, center=(4.65, -0.85)):
        depth = {self.root: 0}
        children = {node: [] for node in self.order}
        for node in self.order[1:]:
            parent = self.parent[node]
            depth[node] = depth[parent] + 1
            children[parent].append(node)
        spans = {}
        leaf = 0
        for node in reversed(self.order):
            if not children[node]:
                spans[node] = (leaf, leaf)
                leaf += 1
            else:
                spans[node] = (
                    min(spans[child][0] for child in children[node]),
                    max(spans[child][1] for child in children[node]),
                )
        maximum_depth = max(depth.values())
        return {
            node: (
                center[0] + (sum(spans[node]) / max(1, leaf - 1) - 1) * width / 2
                if leaf > 1
                else center[0],
                center[1] + (0.5 - depth[node] / maximum_depth) * height
                if maximum_depth
                else center[1],
                0.0,
            )
            for node in self.order
        }


def chord_identity(event):
    return (
        event["left_edge"],
        event["left_side"],
        event["right_edge"],
        event["right_side"],
        event["tag"],
        event["null"],
        event["infinite"],
        event["left_direction"],
        json.dumps([event["left"], event["right"], event["y"]], separators=(",", ":")),
    )


class Trace:
    def __init__(self, path: Path):
        with path.open(encoding="utf-8") as source:
            data = json.load(source)
        if data.get("schema") != 5:
            raise ValueError(
                "Unsupported animation trace schema; record a fresh trace with --trace."
            )
        self.vertices = [finite_point(point) for point in data["vertices"]]
        self.events = data["events"]
        self.curves = {}
        tables = {}
        for event in self.events:
            if event["kind"] == "coordinate_table":
                tables[event["seq"]] = (
                    [finite_point(point) for point in event["points"]],
                    exact_limit(event["horizontal_shift"]),
                )
            elif event["kind"] == "curve":
                self.curves[event["seq"]] = Curve(event["pieces"], tables)
        if len(self.vertices) < 3:
            raise ValueError("A trace needs at least three polygon vertices.")
        for index, event in enumerate(self.events):
            if event["seq"] != index:
                raise ValueError("Trace events must preserve execution order.")
        triangles = [event["vertices"] for event in self.events if event["kind"] == "triangle"]
        if len(triangles) != len(self.vertices) - 2:
            raise ValueError("The recorded triangulation must contain n-2 triangles.")
        for triangle in triangles:
            if len(triangle) != 3 or len(set(triangle)) != 3:
                raise ValueError("A triangle must reference three distinct original vertices.")
            if any(vertex < 0 or vertex >= len(self.vertices) for vertex in triangle):
                raise ValueError("A triangle references an absent vertex.")
        phases = [event["name"] for event in self.events if event["kind"] == "checkpoint"]
        required = ["up_phase", "down_phase", "trapezoids", "unimonotone", "triangulate"]
        if [phase for phase in phases if phase in required] != required:
            raise ValueError("The trace does not contain the complete triangulation pipeline.")
        boundary = next(event for event in self.events if event["kind"] == "boundary")
        self.boundary = [finite_point(point) for point in boundary["points"]]
        self.first_tag = boundary["first_tag"]
        self.submaps = {}
        metadata = None
        records = []
        for event in self.events:
            if event["kind"] == "map_begin":
                if metadata is not None:
                    raise ValueError("Recorded submaps cannot overlap.")
                metadata = event
                records = []
            elif event["kind"] == "map_end":
                if metadata is None:
                    raise ValueError("A recorded submap has no beginning.")
                self.submaps[metadata["seq"]] = Submap(
                    metadata, records, self.curves[metadata["curve"]], 0
                )
                metadata = None
            elif event["kind"] in {"map_chord", "map_arc", "map_region"}:
                if metadata is None:
                    raise ValueError("A submap record is outside its checkpoint.")
                records.append(event)
        if metadata is not None:
            raise ValueError("A recorded submap is incomplete.")

    def phase(self, name):
        start = next(
            index
            for index, event in enumerate(self.events)
            if event["kind"] == "checkpoint" and event["name"] == name
        )
        phases = {"up_phase", "down_phase", "trapezoids", "unimonotone", "triangulate"}
        end = next(
            (
                index
                for index in range(start + 1, len(self.events))
                if self.events[index]["kind"] == "checkpoint"
                and self.events[index]["name"] in phases
            ),
            len(self.events),
        )
        return self.events[start + 1 : end]


class Curve:
    def __init__(self, pieces, tables):
        self.pieces = [
            (tables[table][0], offset, count, reverse, tables[table][1])
            for table, offset, count, reverse in pieces
        ]
        if any(
            count < 2 or offset < 0 or offset + count > len(table)
            for table, offset, count, _, _ in self.pieces
        ):
            raise ValueError("A curve piece leaves its shared coordinate table.")
        if any(
            shift.value != 0
            and (len(table) != 2 or count != 2 or offset != 0 or not shift.infinity)
            for table, offset, count, _, shift in self.pieces
        ):
            raise ValueError(
                "[C91 section 4.2] A wrapped auxiliary edge uses two finite endpoints and an infinite period."
            )
        self.count = sum(count - 1 for _, _, count, _, _ in self.pieces) + 1

    def __len__(self):
        return self.count

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self[position] for position in range(*index.indices(self.count))]
        if index < 0:
            index += self.count
        if index < 0 or index >= self.count:
            raise IndexError(index)
        for table, offset, count, reverse, _ in self.pieces:
            if index < count:
                return table[offset + (count - 1 - index if reverse else index)]
            index -= count - 1
        raise IndexError(index)

    def edge_wrap(self, index):
        if index < 0 or index >= self.count - 1:
            raise IndexError(index)
        for _, _, count, reverse, shift in self.pieces:
            if index < count - 1:
                return -shift.infinity if reverse else shift.infinity
            index -= count - 1
        raise IndexError(index)

    def paths(self, viewport):
        result = []
        for edge in range(self.count - 1):
            append_segments(
                result, viewport.edge_segments(self[edge], self[edge + 1], self.edge_wrap(edge))
            )
        return [[viewport.point(point) for point in points] for points in result]


def append_segments(paths, segments):
    for start, end in segments:
        if paths and paths[-1][-1] == start:
            paths[-1].append(end)
        else:
            paths.append([start, end])


class Viewport:
    def __init__(self, vertices, width=8.2, height=4.9, center=(-2.0, 0.0), padding=0):
        self.minimum_x = min(x for x, _ in vertices)
        self.maximum_x = max(x for x, _ in vertices)
        self.minimum_y = min(y for _, y in vertices)
        self.maximum_y = max(y for _, y in vertices)
        dx = self.maximum_x - self.minimum_x
        dy = self.maximum_y - self.minimum_y
        if not dx and not dy:
            raise ValueError("A viewport needs distinct points.")
        self.width = Fraction(str(width))
        self.height = Fraction(str(height))
        self.scale = min(
            [self.width / dx]
            if not dy
            else [self.height / dy]
            if not dx
            else [self.width / dx, self.height / dy]
        )
        assert 0 <= padding < 0.5
        self.scale *= 1 - 2 * Fraction(str(padding))
        self.middle_x = (self.minimum_x + self.maximum_x) / 2
        self.middle_y = (self.minimum_y + self.maximum_y) / 2
        self.center = center
        self.left = self.middle_x - self.width / (2 * self.scale)
        self.right = self.middle_x + self.width / (2 * self.scale)
        self.bottom = self.middle_y - self.height / (2 * self.scale)
        self.top = self.middle_y + self.height / (2 * self.scale)

    def point(self, point):
        x, y = point
        return (
            float((x - self.middle_x) * self.scale) + self.center[0],
            float((y - self.middle_y) * self.scale) + self.center[1],
            0.0,
        )

    def world_symbolic_point(self, point):
        x, y = (exact_limit(coordinate) for coordinate in point)
        x_value = x.value if x.value is not None else (self.left if x.infinity < 0 else self.right)
        y_value = y.value if y.value is not None else (self.bottom if y.infinity < 0 else self.top)
        return x_value, y_value

    def symbolic_point(self, point):
        x, y = self.world_symbolic_point(point)
        return self.point((min(self.right, max(self.left, x)), min(self.top, max(self.bottom, y))))

    def segment(self, start, end):
        lower, upper = Fraction(0), Fraction(1)
        for axis, minimum, maximum in ((0, self.left, self.right), (1, self.bottom, self.top)):
            delta = end[axis] - start[axis]
            if delta == 0:
                if start[axis] < minimum or start[axis] > maximum:
                    return []
                continue
            a, b = sorted(((minimum - start[axis]) / delta, (maximum - start[axis]) / delta))
            lower, upper = max(lower, a), min(upper, b)
            if lower > upper:
                return []
        return [
            (
                tuple(start[axis] + lower * (end[axis] - start[axis]) for axis in (0, 1)),
                tuple(start[axis] + upper * (end[axis] - start[axis]) for axis in (0, 1)),
            )
        ]

    def edge_segments(self, start, end, wrap=0):
        if not wrap:
            return self.segment(start, end)
        assert start[1] == end[1], (
            "[C91 section 4.2] The symbolic tilt of an exit edge collapses to a horizontal edge."
        )
        first = (self.right if wrap > 0 else self.left, start[1])
        second = (self.left if wrap > 0 else self.right, end[1])
        outgoing = start[0] <= self.right if wrap > 0 else start[0] >= self.left
        incoming = end[0] >= self.left if wrap > 0 else end[0] <= self.right
        return (self.segment(start, first) if outgoing else []) + (
            self.segment(second, end) if incoming else []
        )

    def world_chord_segments(self, event):
        left = exact_limit(event["left"])
        right = exact_limit(event["right"])
        y = exact_limit(event["y"]).value
        if y is None or y < self.bottom or y > self.top:
            return []
        x0 = (
            left.value
            if left.value is not None
            else (self.left if left.infinity < 0 else self.right)
        )
        x1 = (
            right.value
            if right.value is not None
            else (self.left if right.infinity < 0 else self.right)
        )
        intervals = [(self.left, x0), (x1, self.right)] if event["infinite"] else [(x0, x1)]
        return [
            ((max(a, self.left), y), (min(b, self.right), y))
            for a, b in intervals
            if max(a, self.left) <= min(b, self.right)
        ]

    def chord_segments(self, event):
        return [(self.point(a), self.point(b)) for a, b in self.world_chord_segments(event)]
