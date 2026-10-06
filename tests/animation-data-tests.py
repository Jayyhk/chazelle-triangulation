import json
import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "animation"))
from operations import Operations
from regions import region_projections
from trace_data import (
    Curve,
    ExactCoordinate,
    Limit,
    Submap,
    Trace,
    Viewport,
    chord_identity,
    exact_limit,
    finite_point,
    polygon_interior_point,
)

DIRECTORY = Path(sys.argv.pop(1))


def rational(value):
    return {"n": [[[0, 0, 0, 0], str(value)]], "d": [[[0, 0, 0, 0], "1"]]}


def area(points):
    return sum(
        a[0] * b[1] - a[1] * b[0] for a, b in zip(points, points[1:] + points[:1], strict=True)
    )


def inside_cell(point, cell):
    x, y = point
    return all(
        (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]) > 0
        for a, b in zip(cell, cell[1:] + cell[:1], strict=True)
        if a != b
    )


class TraceTests(unittest.TestCase):
    def assert_region_projections(self, submap, viewport):
        projections = region_projections(submap, viewport)
        self.assertEqual(set(projections), set(submap.regions))
        cells = [
            (region, cell)
            for region, projection in projections.items()
            for cell in projection.cells
        ]
        self.assertEqual(
            sum(area(cell) for _, cell in cells),
            2 * (viewport.right - viewport.left) * (viewport.top - viewport.bottom),
        )
        for region, projection in projections.items():
            x, y = projection.marker
            self.assertTrue(
                viewport.left <= x <= viewport.right and viewport.bottom <= y <= viewport.top
            )
            if projection.cells:
                owners = [owner for owner, cell in cells if inside_cell(projection.marker, cell)]
                self.assertEqual(
                    owners, [region], "Each filled marker lies strictly inside its own region."
                )
            exits = [chord for chord in submap.chords if region in chord["regions"]]
            if exits and all(chord["null"] for chord in exits):
                self.assertFalse(
                    projection.cells,
                    "[C91 section 2.1] A null chord's empty region has no filled area.",
                )
        return projections

    def test_markers_for_every_region_at_every_grade(self):
        exterior = collapsed = False
        for path in DIRECTORY.glob("fixture-*.json"):
            trace = Trace(path)
            viewport = Viewport(trace.vertices, padding=0.08)
            for submap in trace.submaps.values():
                with self.subTest(path=path.name, checkpoint=submap.metadata["seq"]):
                    projections = self.assert_region_projections(submap, viewport)
                    exterior |= len(submap.boundary) > 2 and any(
                        projection.cells and not submap.regions[region]["bounded"]
                        for region, projection in projections.items()
                    )
                    collapsed |= any(not projection.cells for projection in projections.values())
        self.assertTrue(exterior, "Later submaps retain spatial markers in exterior regions.")
        self.assertTrue(
            collapsed, "Empty and infinitesimally thin regions retain boundary markers."
        )

    def test_region_markers_after_splits_and_contractions(self):
        for path in DIRECTORY.glob("fixture-*.json"):
            trace = Trace(path)
            viewport = Viewport(trace.vertices, padding=0.08)
            replay = Operations(trace)
            for event in trace.events:
                current = replay.apply(event)
                if event["kind"] == "build_end" or (
                    event["kind"] == "checkpoint" and event["name"] in {"contract_end", "split_end"}
                ):
                    with self.subTest(path=path.name, operation=event["seq"]):
                        self.assert_region_projections(current.submap(trace.curves), viewport)

    def test_region_markers_in_cropped_view(self):
        trace = Trace(next(DIRECTORY.glob("fixture-*.json")))
        full = Viewport(trace.vertices, padding=0.08)
        dx, dy = (full.right - full.left) / 100, (full.top - full.bottom) / 100
        for x, y in (
            (full.left, full.bottom),
            (full.middle_x, full.middle_y),
            (full.right, full.top),
        ):
            viewport = Viewport([(x, y), (x + dx, y + dy)])
            for submap in trace.submaps.values():
                self.assert_region_projections(submap, viewport)

    def test_region_projection_rejects_wrong_chord_direction(self):
        trace = Trace(next(DIRECTORY.glob("fixture-*.json")))
        submap = next(iter(trace.submaps.values()))
        submap.chords[0]["left_direction"] ^= 1
        with self.assertRaises(AssertionError):
            region_projections(submap, Viewport(trace.vertices, padding=0.08))

    def test_single_edge_regions_and_tree(self):
        for path in DIRECTORY.glob("fixture-*.json"):
            trace = Trace(path)
            viewport = Viewport(trace.vertices, padding=0.08)
            for submap in trace.submaps.values():
                if len(submap.boundary) != 2 or submap.boundary.edge_wrap(0):
                    continue
                with self.subTest(path=path.name, checkpoint=submap.metadata["seq"]):
                    polygons = submap.single_edge_regions(viewport)
                    projections = region_projections(submap, viewport)
                    expected = (
                        2 * (viewport.right - viewport.left) * (viewport.top - viewport.bottom)
                    )
                    self.assertEqual(sum(area(points) for points in polygons.values()), expected)
                    low, high = sorted(
                        submap.chords,
                        key=lambda chord: (ExactCoordinate(chord["y"]), -chord["tag"]),
                    )
                    (middle,) = set(low["regions"]) & set(high["regions"])
                    self.assertEqual(len(submap.adjacency[middle]), 2)
                    for region, points in polygons.items():
                        bottom, top = min(y for _, y in points), max(y for _, y in points)
                        self.assertEqual(
                            sum(area(cell) for cell in projections[region].cells), area(points)
                        )
                        self.assertTrue(
                            all(
                                bottom <= y <= top
                                for cell in projections[region].cells
                                for _, y in cell
                            )
                        )
                        _, y = polygon_interior_point(points)
                        if region == middle:
                            self.assertTrue(
                                exact_limit(low["y"]).value < y < exact_limit(high["y"]).value
                            )
                        elif region in low["regions"]:
                            self.assertLess(y, exact_limit(low["y"]).value)
                        else:
                            self.assertGreater(y, exact_limit(high["y"]).value)

    def test_region_marker_inside_concave_face(self):
        points = [
            (Fraction(x), Fraction(y)) for x, y in [(0, 0), (4, 0), (4, 1), (1, 1), (1, 4), (0, 4)]
        ]
        x, y = polygon_interior_point(points)
        self.assertTrue(0 < x < 4 and 0 < y < 4 and (x < 1 or y < 1))

    def assert_boundary_cycles(self, submap):
        for region, record in submap.regions.items():
            if not record["bounded"]:
                continue
            arcs = submap.region_arcs[region]
            exits = [chord for chord in submap.chords if region in chord["regions"]]
            contacts = [
                (exact_limit(chord["left"]), exact_limit(chord["right"]), exact_limit(chord["y"]))
                for chord in exits
            ]
            for first, second in zip(arcs, arcs[1:] + arcs[:1], strict=True):
                end = tuple(exact_limit(coordinate) for coordinate in first["end"])
                start = tuple(exact_limit(coordinate) for coordinate in second["start"])
                self.assertEqual(
                    end[1],
                    start[1],
                    "[C91 section 2.2] Adjacent arcs meet via a horizontal exit chord",
                )
                self.assertTrue(
                    any(
                        y == end[1]
                        and (
                            (left, right) == (end[0], start[0])
                            or (right, left) == (end[0], start[0])
                        )
                        for left, right, y in contacts
                    )
                )

    def test_exact_parameter_order(self):
        epsilon = {"n": [[[0, 0, 0, 1], "1"]], "d": [[[0, 0, 0, 0], "1"]]}
        twice = {"n": [[[0, 0, 0, 1], "2"]], "d": [[[0, 0, 0, 0], "1"]]}
        self.assertEqual(exact_limit(epsilon), exact_limit(twice))
        self.assertLess(ExactCoordinate(epsilon), ExactCoordinate(twice))
        self.assertLess(ExactCoordinate(rational(0)), ExactCoordinate(epsilon))
        self.assertLess(ExactCoordinate(twice), ExactCoordinate(rational(1)))
        negative_denominator = {"n": [[[0, 0, 0, 0], "-1"]], "d": [[[0, 0, 0, 0], "-2"]]}
        self.assertEqual(ExactCoordinate(negative_denominator), ExactCoordinate(rational("1/2")))

    def test_auxiliary_edge_through_infinity(self):
        table = [(Fraction(1), Fraction(2)), (Fraction(3), Fraction(2))]
        viewport = Viewport([(Fraction(0), Fraction(0)), (Fraction(4), Fraction(4))])
        for direction in (-1, 1):
            curve = Curve([[0, 0, 2, False]], {0: (table, Limit(None, direction))})
            components = curve.paths(viewport)
            self.assertEqual(len(components), 2)
            self.assertEqual(components[0][0], viewport.point(table[0]))
            self.assertEqual(
                components[0][-1],
                viewport.point((viewport.right if direction > 0 else viewport.left, Fraction(2))),
            )
            self.assertEqual(
                components[1][0],
                viewport.point((viewport.left if direction > 0 else viewport.right, Fraction(2))),
            )
            self.assertEqual(components[1][-1], viewport.point(table[1]))
            reverse = Curve([[0, 0, 2, True]], {0: (table, Limit(None, direction))})
            self.assertEqual(
                reverse.paths(viewport), [list(reversed(points)) for points in reversed(components)]
            )
            submap = Submap.__new__(Submap)
            submap.boundary, submap.first_tag = curve, 0
            arc = {
                "ranges": [[0, 1]],
                "start": [rational(1), rational(2)],
                "end": [rational(3), rational(2)],
                "wraps": [[0, 0]],
            }
            self.assertEqual(submap.arc_paths(arc, viewport), components)
            arc["wraps"] = []
            arc["end"] = [rational(2), rational(2)]
            self.assertEqual(len(submap.arc_paths(arc, viewport)), 1)

    def test_segment_clipping_before_float(self):
        viewport = Viewport([(Fraction(0), Fraction(0)), (Fraction(4), Fraction(4))])
        large = Fraction(10**10000)
        segments = viewport.segment((-large, -large), (large, large))
        self.assertEqual(len(segments), 1)
        for point in segments[0]:
            self.assertTrue(viewport.left <= point[0] <= viewport.right)
            self.assertTrue(viewport.bottom <= point[1] <= viewport.top)
            self.assertTrue(all(math.isfinite(value) for value in viewport.point(point)))
        self.assertEqual(viewport.segment((-large, large), (-large / 2, large)), [])

    def test_exact_limits(self):
        self.assertEqual(exact_limit(rational("-7/3")).value, Fraction(-7, 3))
        formal = {"n": [[[0, 0, 0, 0], "7/3"], [[1, 0, 0, 0], "1/5"]], "d": [[[0, 0, 0, 0], "1"]]}
        self.assertEqual(exact_limit(formal).value, Fraction(7, 3))
        formal["n"] = [[[1, 0, 0, 0], "1"]]
        self.assertEqual(exact_limit(formal).value, 0)
        formal["n"] = [[[-1, 0, 0, 0], "-1"]]
        self.assertEqual(exact_limit(formal).infinity, -1)
        formal["n"] = [[[1, -500, 0, 0], "1"]]
        self.assertEqual(exact_limit(formal).value, 0)
        formal["d"] = []
        with self.assertRaises(ValueError):
            exact_limit(formal)

    def test_normalize_before_float(self):
        original = [
            (Fraction(0), Fraction(0)),
            (Fraction(4), Fraction(0)),
            (Fraction(4), Fraction(4)),
            (Fraction(0), Fraction(4)),
        ]
        expected = [Viewport(original).point(point) for point in original]
        large = Fraction(10**10000)
        for scale in [large, 1 / large]:
            points = [(large + x * scale, -large + y * scale) for x, y in original]
            actual = [Viewport(points).point(point) for point in points]
            self.assertEqual(actual, expected)
            self.assertTrue(
                all(math.isfinite(component) for point in actual for component in point)
            )

    def test_wrap_and_clip(self):
        viewport = Viewport([(Fraction(0), Fraction(0)), (Fraction(4), Fraction(4))])
        event = {
            "left": rational(1),
            "right": rational(3),
            "y": rational(2),
            "infinite": False,
            "null": False,
        }
        self.assertEqual(
            viewport.chord_segments(event),
            [
                (
                    viewport.point((Fraction(1), Fraction(2))),
                    viewport.point((Fraction(3), Fraction(2))),
                )
            ],
        )
        event["infinite"] = True
        self.assertEqual(len(viewport.chord_segments(event)), 2)
        event["y"] = rational(100)
        self.assertEqual(viewport.chord_segments(event), [])

    def test_pipeline_traces(self):
        paths = list(DIRECTORY.glob("*.json"))
        self.assertGreaterEqual(len(paths), 20)
        refined = False
        for path in paths:
            with self.subTest(path=path.name):
                trace = Trace(path)
                polygon_area = abs(area(trace.vertices))
                triangles = [event for event in trace.events if event["kind"] == "triangle"]
                triangle_area = sum(
                    abs(area([trace.vertices[i] for i in event["vertices"]])) for event in triangles
                )
                self.assertEqual(triangle_area, polygon_area)
                trapezoids = [event for event in trace.events if event["kind"] == "trapezoid"]
                self.assertEqual(
                    sum(
                        abs(area([finite_point(point) for point in event["corners"]]))
                        for event in trapezoids
                    ),
                    polygon_area,
                )
                opened = False
                viewport = Viewport(trace.vertices)
                for curve in trace.curves.values():
                    if not any(shift.infinity for _, _, _, _, shift in curve.pieces):
                        continue
                    for component in curve.paths(viewport):
                        self.assertTrue(
                            all(math.isfinite(value) for point in component for value in point)
                        )
                for submap in trace.submaps.values():
                    if any(arc["wraps"] for arc in submap.arcs):
                        self.assert_region_projections(submap, viewport)
                for event in trace.events:
                    if event["kind"] == "map_begin":
                        self.assertFalse(opened)
                        opened = True
                        refined |= event["name"] == "refined"
                    elif event["kind"] == "map_end":
                        self.assertTrue(opened)
                        opened = False
                    elif event["kind"] == "map_chord":
                        self.assertTrue(opened)
                        for segment in viewport.chord_segments(event):
                            self.assertTrue(
                                all(math.isfinite(value) for point in segment for value in point)
                            )
                self.assertFalse(opened)
        self.assertTrue(refined, "The large fixture exercises actual section 4.2 refinement")

    def test_submap_structure_and_faces(self):
        for path in DIRECTORY.glob("*.json"):
            trace = Trace(path)
            for submap in trace.submaps.values():
                with self.subTest(path=path.name, checkpoint=submap.metadata["seq"]):
                    self.assertEqual(len(submap.regions), len(submap.chords) + 1)
                    positions = submap.tree_positions()
                    self.assertEqual(set(positions), set(submap.regions))
                    self.assertTrue(
                        all(math.isfinite(value) for point in positions.values() for value in point)
                    )
                    for region, record in submap.regions.items():
                        arcs = submap.region_arcs[region]
                        self.assertLessEqual(len(arcs), 4)
                        self.assertEqual(record["weight"], max(arc["edge_count"] for arc in arcs))
                        self.assertLessEqual(record["weight"], submap.metadata["granularity"])
                    self.assert_boundary_cycles(submap)

    def test_boundary_cycles_after_every_split_and_contraction(self):
        for path in DIRECTORY.glob("fixture-*.json"):
            trace = Trace(path)
            replay = Operations(trace)
            for event in trace.events:
                current = replay.apply(event)
                if event["kind"] == "build_end" or (
                    event["kind"] == "checkpoint" and event["name"] in {"contract_end", "split_end"}
                ):
                    with self.subTest(path=path.name, operation=event["seq"]):
                        self.assert_boundary_cycles(current.submap(trace.curves))

    def test_operations_reproduce_every_submap(self):
        wrapped = False
        for path in DIRECTORY.glob("*.json"):
            with self.subTest(path=path.name):
                trace = Trace(path)
                replay = Operations(trace)
                viewport = Viewport(trace.vertices, padding=0.08)
                for event in trace.events:
                    current = replay.apply(event)
                    if event["kind"] == "build_end" and any(
                        arc["wraps"] for arc in current.arcs.values()
                    ):
                        wrapped = True
                        self.assert_region_projections(current.submap(trace.curves), viewport)
                self.assertEqual(replay.verified, len(trace.submaps))
                chains = [event for event in trace.events if event["kind"] == "chain"]
                maximum_grade = max(event["grade"] for event in chains)
                for grade in range(maximum_grade + 1):
                    self.assertEqual(
                        [event["index"] for event in chains if event["grade"] == grade],
                        list(range(2 ** (maximum_grade - grade))),
                    )
        self.assertTrue(
            wrapped,
            "The refinement trace exercises constructed regions with auxiliary edges through infinity",
        )

    def test_operations_detect_missing_change(self):
        trace = Trace(next(DIRECTORY.glob("fixture-*.json")))
        removed = next(event for event in trace.events if event["kind"] == "build_arc")
        trace.events.remove(removed)
        with self.assertRaises((AssertionError, ValueError)):
            Operations(trace).verify()

    def test_builder_detects_missing_discovered_chord(self):
        trace = Trace(next(DIRECTORY.glob("fixture-*.json")))
        removed = next(event for event in trace.events if event["kind"] == "visibility_chord")
        trace.events = [
            event
            for event in trace.events
            if not (
                event["kind"] == "visibility_chord"
                and event["curve"] == removed["curve"]
                and chord_identity(event) == chord_identity(removed)
            )
        ]
        with self.assertRaises(AssertionError):
            Operations(trace).verify()

    def test_line_chain_viewports(self):
        for points in [
            [(Fraction(0), Fraction(0)), (Fraction(4), Fraction(0))],
            [(Fraction(0), Fraction(0)), (Fraction(0), Fraction(4))],
        ]:
            viewport = Viewport(points)
            self.assertTrue(
                all(math.isfinite(value) for point in points for value in viewport.point(point))
            )

    def test_chord_identity_preserves_sides(self):
        trace = Trace(next(DIRECTORY.glob("fixture-*.json")))
        chord = next(event for event in trace.events if event["kind"] == "map_chord")
        changed = dict(chord)
        changed["left_side"] = 1 - chord["left_side"]
        self.assertNotEqual(chord_identity(chord), chord_identity(changed))

    def test_reject_incomplete_trace(self):
        path = next(DIRECTORY.glob("fixture-*.json"))
        data = json.loads(path.read_text())
        data["events"].pop()
        incomplete = DIRECTORY / "incomplete.invalid"
        incomplete.write_text(json.dumps(data))
        try:
            with self.assertRaises(ValueError):
                Trace(incomplete)
        finally:
            incomplete.unlink()


if __name__ == "__main__":
    unittest.main()
