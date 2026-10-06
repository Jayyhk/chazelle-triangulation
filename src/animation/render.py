import argparse
import colorsys
import os
import tempfile
from fractions import Fraction
from pathlib import Path

import av
import numpy as np
from manim import (
    Circle,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    MoveAlongPath,
    Polygon,
    RoundedRectangle,
    Scene,
    Transform,
    VGroup,
    VMobject,
    config,
    tempconfig,
)

from operations import Operations
from regions import region_projections
from trace_data import Trace, Viewport, append_segments, chord_identity, finite_point

BACKGROUND = "#101823"
BOUNDARY = "#ECF2F8"
BLUE = "#6EBBFF"
GREEN = "#79E2AE"
ORANGE = "#FFB86B"
RED = "#FF7F91"
MUTED = "#7891AA"
HIGHLIGHT = "#FFE49D"


def color(index):
    rgb = colorsys.hsv_to_rgb((index * 0.618033988749895) % 1, 0.43, 0.96)
    return "#" + "".join(f"{round(channel * 255):02X}" for channel in rgb)


def path(points, tint, width=2.5):
    return VMobject(color=tint, stroke_width=width).set_points_as_corners(
        [np.array(point) for point in points]
    )


def paths(components, tint, width=2.5):
    return VGroup(*(path(points, tint, width) for points in components))


def draw_chord(viewport, event, tint=BLUE):
    result = VGroup()
    for first, second in viewport.chord_segments(event):
        if np.linalg.norm(np.array(second) - first) < 0.001:
            result.add(Circle(radius=0.045, color=tint, stroke_width=1.8).move_to(first))
        elif event["infinite"]:
            result.add(
                DashedLine(
                    first, second, color=tint, stroke_width=1.2, dash_length=0.1
                ).set_opacity(0.45)
            )
        else:
            result.add(Line(first, second, color=tint, stroke_width=2.2))
    return result


class Drawing:
    def __init__(self, submap, viewport, colors, center=(4.75, 0.65), height=4.6, width=3.5):
        self.submap = submap
        self.colors = colors
        self.positions = submap.tree_positions(width=width, height=height, center=center)
        if (
            len(submap.boundary) == 2
            and not submap.boundary.edge_wrap(0)
            and len(submap.chords) == 2
        ):
            self.positions = {
                region: (center[0], center[1] + (index - 1) * height / 2, 0)
                for index, region in enumerate(submap.single_edge_region_order())
            }
        radius = max(0.009, min(0.08, 0.58 / len(self.positions) ** 0.5))
        self.nodes = {
            region: Dot(point, radius=radius, color=colors[region])
            for region, point in self.positions.items()
        }
        self.edges = {
            chord["id"]: Line(
                *(self.positions[region] for region in chord["regions"]),
                color=MUTED,
                stroke_width=1.5,
            )
            for chord in submap.chords
        }
        self.graph = VGroup(*self.edges.values(), *self.nodes.values())
        self.chords = {chord["id"]: draw_chord(viewport, chord) for chord in submap.chords}
        self.faces = {}
        self.projections = region_projections(submap, viewport)
        self.markers = {}
        collapsed = {}
        for region, projection in self.projections.items():
            point = viewport.point(projection.marker)
            if projection.cells:
                self.faces[region] = VGroup(
                    *(
                        Polygon(
                            *(viewport.point(point) for point in cell),
                            color=colors[region],
                            stroke_width=0,
                            fill_color=colors[region],
                            fill_opacity=0.13,
                        )
                        for cell in projection.cells
                    )
                )
                marker = Dot(point, radius=0.055, color=colors[region])
            else:
                count = collapsed.get(projection.marker, 0)
                collapsed[projection.marker] = count + 1
                marker = Circle(
                    radius=0.07 + count * 0.035, color=colors[region], stroke_width=2
                ).move_to(point)
            self.markers[region] = marker.set_z_index(4)
        assert set(self.markers) == set(self.nodes), (
            "[C91 section 2.1] Every region has its own matching marker and tree node."
        )
        self.arcs = {
            arc["id"]: paths(submap.arc_paths(arc, viewport), colors[arc["region"]], 2).set_stroke(
                opacity=0.55
            )
            for arc in submap.arcs
        }
        self.geometry = VGroup(
            *self.faces.values(), *self.arcs.values(), *self.chords.values(), *self.markers.values()
        )


class TriangulationScene(Scene):
    def __init__(self, trace, **kwargs):
        self.trace = trace
        self.operations = Operations(trace)
        self.built = {}
        self.chains_per_grade = {}
        verification = Operations(trace)
        for event in trace.events:
            current = verification.apply(event)
            if event["kind"] == "build_end":
                self.built[event["map"]] = current.submap(trace.curves)
            elif event["kind"] == "chain":
                grade = event["grade"]
                self.chains_per_grade[grade] = self.chains_per_grade.get(grade, 0) + 1
        assert verification.verified == len(trace.submaps)
        self.grade_spacing = min(0.29, 1.45 / max(self.chains_per_grade))
        self.viewport = Viewport(
            trace.vertices, width=8.6, height=4.9, center=(-2.05, 0.65), padding=0.08
        )
        self.points = [self.viewport.point(point) for point in trace.vertices]
        self.geometry = VGroup()
        self.graph = VGroup()
        self.active_curve = VGroup()
        self.query_curve = None
        self.cursor = None
        self.ring = None
        self.colors = {}
        self.map_thumbnail = {}
        self.shelf = VGroup()
        self.walk_arcs = VGroup()
        self.faces = {}
        self.completed_regions = set()
        self.overlay = VGroup()
        self.merge_event = None
        self.chain = None
        self.fusion_chords = {}
        self.nodes = {}
        self.edges = {}
        self.chord_shapes = {}
        self.discovered = VGroup()
        self.discovered_curve = None
        self.discovered_identities = set()
        super().__init__(**kwargs)

    def discard(self, group):
        self.remove(*group.get_family())

    def clear(self, *groups, duration=0.2):
        visible = [group for group in groups if len(group)]
        if visible:
            self.play(*(FadeOut(group) for group in visible), run_time=duration)
            for group in visible:
                self.discard(group)

    def curve_path(self, identity, tint=BOUNDARY, viewport=None):
        viewport = viewport or self.viewport
        curve = self.trace.curves[identity]
        return paths(curve.paths(viewport), tint, 2.8)

    def show_curve(self, identity):
        if self.query_curve == identity:
            return
        self.clear(self.discovered, duration=0.12)
        self.discovered = VGroup()
        self.discovered_curve = identity
        self.discovered_identities = set()
        target = self.curve_path(identity).set_z_index(2)
        if len(self.active_curve):
            self.play(FadeOut(self.active_curve), FadeIn(target), run_time=0.15)
            self.discard(self.active_curve)
            self.active_curve = VGroup(target)
        else:
            self.active_curve = VGroup(target)
            self.play(Create(target), run_time=0.2)
        self.query_curve = identity

    def map_colors(self, current):
        return self.region_colors(current.identity, current.regions)

    def region_colors(self, identity, regions):
        colors = self.colors.setdefault(identity, {})
        for region in regions:
            colors.setdefault(region, color(identity + region))
        return colors

    def query(self, event):
        self.show_curve(event["curve"])
        origin = self.viewport.symbolic_point(event["origin"])
        source = Dot(origin, radius=0.06, color=HIGHLIGHT).set_z_index(6)
        direction = -1 if event["direction"] == 0 else 1
        contact = (
            self.viewport.symbolic_point(event["contact"])
            if event["hit"]
            else (
                self.viewport.point(
                    (
                        self.viewport.left if direction < 0 else self.viewport.right,
                        finite_point(event["origin"])[1],
                    )
                )
            )
        )
        intervals = [(origin, contact)]
        if event["hit"] and event["wrapped"]:
            boundary = self.viewport.point(
                (
                    self.viewport.left if direction < 0 else self.viewport.right,
                    finite_point(event["origin"])[1],
                )
            )
            opposite = self.viewport.point(
                (
                    self.viewport.right if direction < 0 else self.viewport.left,
                    finite_point(event["origin"])[1],
                )
            )
            intervals = [(origin, boundary), (opposite, contact)]
        self.play(FadeIn(source), run_time=0.12)
        drawn = VGroup()
        for index, (first, second) in enumerate(intervals):
            if index:
                gates = VGroup(
                    *(
                        Circle(radius=0.085, color=HIGHLIGHT, stroke_width=2).move_to(point)
                        for point in (intervals[index - 1][1], first)
                    )
                )
                self.play(Create(gates), run_time=0.18)
                self.remove(gates)
            if np.linalg.norm(np.array(second) - first) < 0.002:
                pulse = Circle(radius=0.14, color=HIGHLIGHT, stroke_width=2).move_to(first)
                self.play(Create(pulse), run_time=0.22)
                drawn.add(pulse)
            else:
                segment = Line(first, second, color=HIGHLIGHT, stroke_width=2.5)
                moving = Dot(first, color=HIGHLIGHT, radius=0.045).set_z_index(7)
                self.add(moving)
                self.play(Create(segment), MoveAlongPath(moving, segment), run_time=0.32)
                self.remove(moving)
                drawn.add(segment)
        if event["hit"]:
            curve = self.trace.curves[event["curve"]]
            edge = event["edge"]
            segments = self.viewport.edge_segments(
                curve[edge], curve[edge + 1], curve.edge_wrap(edge)
            )
            struck = VGroup(
                *(
                    Line(
                        self.viewport.point(a), self.viewport.point(b), color=GREEN, stroke_width=5
                    )
                    for a, b in segments
                    if a != b
                )
            ).set_z_index(5)
            hit = Circle(radius=0.1, color=GREEN, stroke_width=2).move_to(contact).set_z_index(6)
            drawn.add(struck, hit)
            self.play(FadeIn(struck), Create(hit), run_time=0.15)
        self.play(FadeOut(source), FadeOut(drawn), run_time=0.15)
        self.discard(drawn)
        self.remove(source)

    def discover(self, event):
        self.show_curve(event["curve"])
        identity = chord_identity(event)
        if identity in self.discovered_identities:
            return
        self.discovered_identities.add(identity)
        chord = draw_chord(self.viewport, event).set_stroke(opacity=0.85).set_z_index(1)
        self.discovered.add(chord)
        self.play(Create(chord), run_time=0.28)

    def begin_build(self, event, current):
        self.clear(self.geometry, self.graph, self.walk_arcs)
        if self.cursor is not None:
            self.remove(self.cursor)
        if self.ring is not None:
            self.remove(self.ring)
        self.show_curve(current.curve)
        self.geometry = VGroup()
        self.graph = VGroup()
        self.nodes = {}
        self.edges = {}
        self.pending_edges = {}
        self.walk_arcs = VGroup()
        self.faces = {}
        self.completed_regions = set()
        self.building = Drawing(
            self.built[current.identity],
            self.viewport,
            self.region_colors(current.identity, self.built[current.identity].regions),
        )
        root = self.building.nodes[0]
        self.nodes[0] = root
        self.graph.add(root)
        first = self.trace.curves[current.curve][0]
        self.cursor = Dot(self.viewport.point(first), radius=0.065, color=HIGHLIGHT).set_z_index(6)
        self.ring = Circle(radius=0.115, color=HIGHLIGHT, stroke_width=2).move_to(root)
        self.play(FadeIn(self.cursor), run_time=0.2)
        self.add_region_node(0, root, run_time=0.55)
        self.play(Create(self.ring), run_time=0.2)

    def add_region_node(self, region, node, edge=None, run_time=0.6):
        marker = self.building.markers[region]
        self.geometry.add(marker)
        self.play(FadeIn(marker), run_time=0.2)
        moving = marker.copy()
        self.add(moving)
        animations = [Transform(moving, node.copy()), FadeIn(node)]
        if edge is not None:
            animations.append(Create(edge))
        self.play(*animations, run_time=run_time)
        self.remove(moving)

    def complete_region(self, region):
        if region in self.completed_regions:
            return
        self.completed_regions.add(region)
        animations = []
        if region in self.building.faces:
            face = self.building.faces[region].set_z_index(-1)
            self.faces[region] = face
            self.geometry.add(face)
            animations.append(FadeIn(face))
        node = self.nodes[region]
        pulses = VGroup(
            *(
                Circle(radius=0.14, color=self.building.colors[region], stroke_width=2).move_to(
                    point
                )
                for point in (node.get_center(), self.building.markers[region].get_center())
            )
        )
        self.play(*animations, Create(pulses), run_time=0.45)
        self.play(FadeOut(pulses), run_time=0.18)
        self.discard(pulses)

    def walk(self, event, current):
        components = self.building.submap.arc_paths(event, self.viewport)
        if not components:
            self.play(
                self.ring.animate.move_to(self.building.nodes[event["region"]]), run_time=0.25
            )
        for coordinates in components:
            arc = path(
                coordinates, self.colors[current.identity][event["region"]], 3.2
            ).set_z_index(3)
            traveled = path(coordinates, HIGHLIGHT, 4).set_z_index(5)
            self.cursor.move_to(coordinates[0])
            self.play(
                Create(traveled),
                MoveAlongPath(self.cursor, traveled),
                self.ring.animate.move_to(self.building.nodes[event["region"]]),
                run_time=0.25,
            )
            self.add(arc)
            self.walk_arcs.add(arc)
            self.remove(traveled)

    def enter(self, event, current):
        node = self.building.nodes[event["region"]]
        parent = self.building.nodes[event["parent"]]
        chord = draw_chord(self.viewport, event, GREEN).set_z_index(3)
        graph_edge = Line(parent.get_center(), node.get_center(), color=GREEN, stroke_width=2)
        self.play(Create(chord.set_stroke(opacity=1)), run_time=0.4)
        self.add_region_node(event["region"], node, graph_edge)
        self.wait(0.25)
        self.play(
            self.ring.animate.move_to(node),
            graph_edge.animate.set_color(MUTED).set_stroke(width=1.5),
            run_time=0.25,
        )
        self.graph.add(graph_edge, node)
        self.nodes[event["region"]] = node
        self.pending_edges[chord_identity(event)] = graph_edge
        self.geometry.add(chord)

    def leave(self, event):
        self.complete_region(event["region"])
        self.play(self.ring.animate.move_to(self.building.nodes[event["parent"]]), run_time=0.4)

    def record_arc(self, event, current):
        if current.complete:
            return
        region = event["region"]
        required = {arc["id"] for arc in self.building.submap.region_arcs[region]}
        recorded = {index for index, arc in current.arcs.items() if arc["region"] == region}
        if recorded == required and region not in self.completed_regions:
            self.complete_region(region)

    def drawing(self, current):
        return Drawing(current.submap(self.trace.curves), self.viewport, self.map_colors(current))

    def refresh(self, current, duration=0.3):
        target = self.drawing(current)
        target.geometry.set_z_index(0)
        assert set(self.nodes) == set(current.regions)
        assert set(self.edges) == set(current.chords)
        changes = [
            node.animate.move_to(target.positions[region]).set_color(target.colors[region])
            for region, node in self.nodes.items()
        ]
        changes.extend(
            self.edges[index].animate.put_start_and_end_on(
                *(target.positions[region] for region in chord["regions"])
            )
            for index, chord in current.chords.items()
        )
        self.play(FadeOut(self.geometry), FadeIn(target.geometry), *changes, run_time=duration)
        self.discard(self.geometry)
        self.geometry = target.geometry
        self.chord_shapes = target.chords
        if self.ring is not None:
            self.remove(self.ring)
            self.ring = None
        return target

    def end_build(self, current):
        self.complete_region(0)
        self.clear(self.walk_arcs, duration=0.1)
        if self.cursor is not None:
            self.remove(self.cursor)
            self.cursor = None
        self.refresh(current, 0.2)
        self.clear(self.discovered, duration=0.15)
        self.discovered = VGroup()
        self.discovered_identities = set()

    def copy_map(self, event, current):
        self.colors[current.identity] = dict(self.colors[event["source"]])
        self.show_curve(current.curve)
        self.clear(self.geometry, self.graph)
        drawing = self.drawing(current)
        self.nodes, self.edges = drawing.nodes, drawing.edges
        self.chord_shapes = drawing.chords
        self.geometry, self.graph = drawing.geometry, drawing.graph
        source = self.map_thumbnail.get(event["source"])
        if source is not None:
            selected = source.copy().set_color(HIGHLIGHT).set_z_index(4)
            self.play(FadeIn(selected), FadeIn(self.geometry), FadeIn(self.graph), run_time=0.4)
            self.remove(selected)
            self.map_thumbnail[current.identity] = source
        else:
            self.play(FadeIn(self.geometry), FadeIn(self.graph), run_time=0.3)

    def merge_inputs(self, event):
        self.clear(self.geometry, self.graph)
        shapes = VGroup()
        trees = VGroup()
        self.fusion_chords = {}
        for index, identity in enumerate((event["first"], event["second"])):
            current = self.operations.maps[identity]
            submap = current.submap(self.trace.curves)
            drawing = Drawing(
                submap,
                self.viewport,
                self.map_colors(current),
                center=(4.75, 1.9 - index * 2.5),
                height=2.0,
            )
            for chord_id, chord in drawing.chords.items():
                chord.set_color(ORANGE if index else GREEN)
                self.fusion_chords[(identity, chord_id)] = (chord, drawing.edges[chord_id])
            shapes.add(drawing.geometry)
            trees.add(drawing.graph)
        self.geometry, self.graph = shapes, trees
        self.show_curve(self.merge_event["curve"])
        self.play(FadeIn(shapes), FadeIn(trees), run_time=0.4)

    def fusion_chord(self, event):
        chord = draw_chord(self.viewport, event, HIGHLIGHT).set_z_index(4)
        self.geometry.add(chord)
        self.play(Create(chord), run_time=0.25)

    def fusion_remove(self, event):
        removed = draw_chord(self.viewport, event, RED).set_z_index(5)
        self.play(Create(removed), run_time=0.18)
        self.play(FadeOut(removed), run_time=0.18)
        self.discard(removed)
        chord, edge = self.fusion_chords[(event["map"], event["id"])]
        self.play(
            chord.animate.set_color(RED).set_stroke(opacity=0.3),
            edge.animate.set_color(RED).set_opacity(0.65),
            run_time=0.15,
        )

    def contract(self, event, current):
        removed = draw_chord(self.viewport, event, RED).set_z_index(5)
        keep, dead = event["regions"]
        edge = self.edges.pop(event["id"])
        node = self.nodes.pop(dead)
        chord = self.chord_shapes.pop(event["id"])
        self.play(Create(removed), edge.animate.set_color(RED).set_stroke(width=3), run_time=0.18)
        self.play(
            node.animate.move_to(self.nodes[keep]),
            FadeOut(edge),
            FadeOut(chord),
            FadeOut(removed),
            run_time=0.3,
        )
        self.graph.remove(node, edge)
        self.geometry.remove(chord)
        self.remove(node, edge)
        self.discard(chord)
        self.discard(removed)

    def split(self, event, current):
        tint = self.map_colors(current)
        parent, region = event["regions"]
        point = self.nodes[parent].get_center() + np.array([0.65, -0.45, 0])
        point[0] = np.clip(point[0], 3.15, 6.4)
        point[1] = np.clip(point[1], -1.75, 2.9)
        node = Dot(point, radius=self.nodes[parent].width / 2, color=tint[region])
        edge = Line(self.nodes[parent].get_center(), point, color=MUTED, stroke_width=1.5)
        chord = draw_chord(self.viewport, event, GREEN).set_z_index(4)
        self.geometry.add(chord)
        self.play(Create(chord), run_time=0.3)
        moving = chord.copy()
        self.add(moving)
        self.play(Transform(moving, edge.copy()), FadeIn(node), run_time=0.35)
        self.remove(moving)
        self.add(edge)
        self.graph.add(edge, node)
        self.nodes[region] = node
        self.edges[event["id"]] = edge
        self.chord_shapes[event["id"]] = chord

    def snapshot(self, event):
        recorded = self.trace.submaps[event["seq"]]
        if event["name"] != "canonical":
            return
        grade, index = self.chain["grade"], self.chain["index"]
        count = self.chains_per_grade[grade]
        x = -5.95 + (index + 0.5) * 11.9 / count
        y = -3.55 + grade * self.grade_spacing
        width = min(0.65, 11 / count)
        viewport = Viewport(
            recorded.boundary[:], width=width - min(0.06, width / 10), height=0.19, center=(x, y)
        )
        curve = self.curve_path(recorded.metadata["curve"], MUTED, viewport).set_stroke(width=1)
        chords = VGroup(
            *(
                draw_chord(viewport, chord).set_stroke(width=0.7)
                for chord in recorded.chords
                if not chord["infinite"] and not chord["null"]
            )
        )
        frame = RoundedRectangle(
            width=width,
            height=0.25,
            corner_radius=min(0.025, width / 4),
            color="#31485D",
            stroke_width=0.7,
        ).move_to((x, y, 0))
        thumbnail = VGroup(frame, curve, chords)
        self.shelf.add(thumbnail)
        self.map_thumbnail[recorded.metadata["map"]] = thumbnail
        self.play(FadeIn(curve), FadeIn(frame), FadeIn(chords), run_time=0.2)

    def output_event(self, event):
        kind = event["kind"]
        if kind == "trapezoid":
            shape = Polygon(
                *(self.viewport.point(finite_point(point)) for point in event["corners"]),
                color=BLUE,
                stroke_width=1,
                fill_color=BLUE,
                fill_opacity=0.16,
            )
        elif kind == "diagonal":
            shape = Line(
                *(self.points[index] for index in event["vertices"]),
                color=HIGHLIGHT,
                stroke_width=2.5,
            )
        else:
            tint = color(event["seq"])
            shape = Polygon(
                *(self.points[index] for index in event["vertices"]),
                color=tint,
                fill_color=tint,
                fill_opacity=0.3,
                stroke_width=1.5,
            )
        self.overlay.add(shape)
        self.play(Create(shape) if kind == "diagonal" else FadeIn(shape), run_time=0.4)

    def construct(self):
        config.background_color = BACKGROUND
        components = []
        for a, b in zip(
            self.trace.vertices, self.trace.vertices[1:] + self.trace.vertices[:1], strict=True
        ):
            append_segments(components, self.viewport.segment(a, b))
        outline = (
            paths(
                [[self.viewport.point(point) for point in component] for component in components],
                MUTED,
                1.2,
            )
            .set_stroke(opacity=0.25)
            .set_z_index(-2)
        )
        vertices = VGroup(
            *(
                Dot(self.viewport.point(point), radius=0.03, color=BOUNDARY).set_opacity(0.4)
                for point in self.trace.vertices
                if self.viewport.left <= point[0] <= self.viewport.right
                and self.viewport.bottom <= point[1] <= self.viewport.top
            )
        )
        self.add(outline, vertices)
        frame = (
            RoundedRectangle(
                width=float(self.viewport.width),
                height=float(self.viewport.height),
                corner_radius=0,
                color=MUTED,
                stroke_width=0.8,
            )
            .move_to((*self.viewport.center, 0))
            .set_stroke(opacity=0.25)
            .set_z_index(-3)
        )
        self.add(frame)
        self.wait(0.4)
        for event in self.trace.events:
            kind = event["kind"]
            old_colors = dict(self.colors.get(event.get("map"), {}))
            current = self.operations.apply(event)
            if kind == "boundary":
                added = set(self.trace.boundary) - set(self.trace.vertices)
                dots = VGroup(
                    *(
                        Dot(self.viewport.point(point), radius=0.035, color=BOUNDARY)
                        for point in sorted(added)
                    )
                )
                if len(dots):
                    self.play(FadeIn(dots), run_time=0.4)
            elif kind == "chain":
                self.chain = event
                self.clear(self.geometry, self.graph, self.walk_arcs)
                self.geometry, self.graph, self.walk_arcs = VGroup(), VGroup(), VGroup()
                self.show_curve(event["curve"])
            elif kind == "ray":
                self.query(event)
            elif kind == "visibility_chord":
                self.discover(event)
            elif kind == "merge":
                self.merge_event = event
            elif kind == "merge_inputs":
                self.merge_inputs(event)
            elif kind == "fusion_cursor":
                marker = Circle(radius=0.09, color=ORANGE, stroke_width=2).move_to(
                    self.viewport.symbolic_point(event["point"])
                )
                self.play(Create(marker), run_time=0.12)
                self.remove(marker)
            elif kind == "fusion_chord":
                self.fusion_chord(event)
            elif kind == "fusion_remove":
                self.fusion_remove(event)
            elif kind == "build_begin":
                self.begin_build(event, current)
            elif kind == "build_walk":
                self.walk(event, current)
            elif kind == "build_enter":
                self.enter(event, current)
            elif kind == "build_leave":
                self.leave(event)
            elif kind == "build_arc":
                self.record_arc(event, current)
            elif kind == "build_chord" and not current.complete:
                self.edges[event["id"]] = self.pending_edges[chord_identity(event)]
            elif kind == "build_end":
                self.end_build(current)
            elif kind == "copy_submap":
                self.copy_map(event, current)
            elif kind == "remove_chord":
                self.contract(event, current)
            elif kind == "insert_chord":
                self.split(event, current)
            elif kind == "reindex":
                self.colors[event["map"]] = {
                    event["regions"][region]: tint
                    for region, tint in old_colors.items()
                    if event["regions"][region] in self.operations.maps[event["map"]].regions
                }
                self.nodes = {event["regions"][region]: node for region, node in self.nodes.items()}
                self.edges = {event["chords"][index]: edge for index, edge in self.edges.items()}
                self.chord_shapes = {
                    event["chords"][index]: chord for index, chord in self.chord_shapes.items()
                }
            elif kind == "map_begin":
                self.snapshot(event)
            elif kind == "checkpoint":
                if event["name"] in {"split_end", "contract_end"}:
                    self.refresh(current)
                elif event["name"] in {"trapezoids", "unimonotone", "triangulate"}:
                    self.clear(self.geometry, self.graph, self.active_curve, self.overlay)
                    self.geometry, self.graph, self.active_curve, self.overlay = (
                        VGroup() for _ in range(4)
                    )
                    if self.ring is not None:
                        self.remove(self.ring)
                        self.ring = None
            elif kind in {"trapezoid", "diagonal", "piece", "triangle"}:
                self.output_event(event)
        self.wait(3)


def write_video(source, destination):
    with (
        av.open(str(source)) as decoder,
        av.open(str(destination), "w", options={"movflags": "+faststart"}) as encoder,
    ):
        original = decoder.streams.video[0]
        stream = encoder.add_stream("libx264", rate=24, options={"crf": "18", "preset": "veryfast"})
        stream.width, stream.height = original.width, original.height
        stream.pix_fmt = "yuv420p"
        for index, frame in enumerate(decoder.decode(video=0)):
            frame.pts = index
            frame.time_base = Fraction(1, 24)
            for packet in stream.encode(frame):
                encoder.mux(packet)
        for packet in stream.encode():
            encoder.mux(packet)


def main():
    parser = argparse.ArgumentParser(description="Replay an exact Chazelle C++ trace with Manim.")
    parser.add_argument("trace", type=Path)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    if arguments.output.suffix != ".mp4":
        parser.error("The output must be an .mp4 file.")
    trace = Trace(arguments.trace)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="chazelle-manim-") as directory:
        with tempconfig(
            {
                "media_dir": directory,
                "output_file": "animation",
                "pixel_width": 2560,
                "pixel_height": 1440,
                "frame_width": 128 / 9,
                "frame_height": 8,
                "renderer": "cairo",
                "frame_rate": 24,
                "disable_caching": True,
                "verbosity": "ERROR",
                "progress_bar": "none",
                "preview": False,
                "write_to_movie": True,
                "background_color": BACKGROUND,
            }
        ):
            scene = TriangulationScene(trace)
            scene.render()
            rendered = Path(scene.renderer.file_writer.movie_file_path)
        descriptor, temporary = tempfile.mkstemp(suffix=".mp4", dir=arguments.output.parent)
        os.close(descriptor)
        try:
            write_video(rendered, temporary)
            os.replace(temporary, arguments.output)
        finally:
            Path(temporary).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
