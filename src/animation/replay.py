import networkx as nx
from manim import (
    Arrow,
    Circle,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Square,
    SurroundingRectangle,
    Transform,
    TransformFromCopy,
    VGroup,
)

from replay_steps import NONE, ReplaySteps
from trace_data import exact_limit

BLUE = "#6EBBFF"
GREEN = "#79E2AE"
RED = "#FF7F91"
GOLD = "#FFE49D"
MUTED = "#7891AA"


class VisualReplay:
    def __init__(self, scene, drawing, chord, palette, paths):
        self.scene = scene
        self.data = ReplaySteps(scene.trace)
        self.drawing = drawing
        self.chord = chord
        self.palette = palette
        self.paths = paths
        self.panel = VGroup()
        self.panel_id = None
        self.panels = {}
        self.glyphs = {}
        self.panel_edges = {}
        self.query_stack = []
        self.query_sources = {}
        self.query_contexts = {}
        self.refinements = {}
        self.refinement_stack = []
        self.capacity = VGroup()
        self.decision = VGroup()
        self.trapezoids = {}
        self.piece_outlines = VGroup()
        self.triangles = VGroup()
        self.active_piece = None
        self.remaining = []
        self.piece_boundary = VGroup()
        self.cursor = None
        self.test = VGroup()
        self.current_vertex = None
        self.pending_change = None
        self.context = None

    def pulse(self, shape, tint=GOLD, duration=0.16):
        if not len(shape.get_all_points()):
            return
        highlight = shape.copy().set_color(tint).set_opacity(0.8).set_z_index(7)
        self.scene.play(FadeIn(highlight), run_time=duration)
        self.scene.play(FadeOut(highlight), run_time=duration)
        self.scene.discard(highlight)

    def hide_panel(self):
        self.scene.clear(self.panel, self.capacity, self.decision, duration=0.12)
        self.panel, self.capacity, self.decision = VGroup(), VGroup(), VGroup()
        self.panel_id = None
        self.scene.graph.set_opacity(1)

    def show_panel(self, identity):
        if identity == self.panel_id:
            return
        self.scene.clear(self.panel, duration=0.12)
        self.scene.graph.set_opacity(0.08)
        self.panel_id = identity
        self.panel = self.panels[identity]
        self.scene.play(FadeIn(self.panel), run_time=0.25)

    def save_context(self):
        scene = self.scene
        return (
            scene.geometry,
            scene.graph,
            scene.active_curve,
            scene.query_curve,
            scene.nodes,
            scene.edges,
            scene.chord_shapes,
            self.panel,
            self.panel_id,
            scene.active_drawing,
        )

    def restore_context(self, saved):
        scene = self.scene
        scene.clear(scene.geometry, scene.graph, scene.active_curve, self.panel, duration=0.15)
        (
            scene.geometry,
            scene.graph,
            scene.active_curve,
            scene.query_curve,
            scene.nodes,
            scene.edges,
            scene.chord_shapes,
            self.panel,
            self.panel_id,
            scene.active_drawing,
        ) = saved
        scene.play(
            *(
                FadeIn(group)
                for group in (scene.geometry, scene.graph, scene.active_curve, self.panel)
                if len(group)
            ),
            run_time=0.3,
        )

    def display_map(self, model, owner, curve):
        scene = self.scene
        scene.clear(scene.geometry, scene.graph, scene.active_curve, self.panel, duration=0.15)
        view = self.drawing(model, scene.viewport, scene.region_colors(owner, model.regions))
        scene.geometry, scene.graph = view.geometry, view.graph
        scene.active_drawing = view
        scene.nodes, scene.edges, scene.chord_shapes = view.nodes, view.edges, view.chords
        scene.active_curve = VGroup(scene.curve_path(curve))
        scene.query_curve = curve
        self.panel, self.panel_id = VGroup(), None
        source = scene.map_thumbnail.get(owner)
        if source is not None:
            moving, moving_graph = scene.geometry.copy(), scene.graph.copy()
            source_graph = view.graph.copy().scale(0.1).move_to(source)
            scene.play(
                TransformFromCopy(source, moving),
                TransformFromCopy(source_graph, moving_graph),
                Create(scene.active_curve),
                run_time=0.45,
            )
            scene.discard(moving)
            scene.discard(moving_graph)
            scene.add(scene.geometry, scene.graph)
        else:
            scene.play(
                FadeIn(scene.geometry),
                FadeIn(scene.graph),
                Create(scene.active_curve),
                run_time=0.35,
            )

    def owner_map(self, owner):
        current = self.scene.operations.maps[owner]
        return current.submap(self.scene.trace.curves)

    def arc(self, model, index, tint=GOLD):
        arc = next(arc for arc in model.arcs if arc["id"] == index)
        return self.paths(model.arc_paths(arc, self.scene.viewport), tint, 5).set_z_index(6)

    def region(self, model, owner, region):
        view = self.drawing(
            model, self.scene.viewport, self.scene.region_colors(owner, model.regions)
        )
        return VGroup(
            view.faces.get(region, VGroup()),
            view.markers[region],
            *(view.arcs[arc["id"]] for arc in model.region_arcs[region]),
        )

    def centroid_positions(self, tree):
        nodes = tree.nodes
        leaves = [node for node in nodes.values() if node["chord"] == NONE]
        depths = {tree.event["root"]: 0}
        for node in nodes.values():
            if node["parent"] != NONE:
                depths[node["index"]] = depths[node["parent"]] + 1
        maximum = max(depths.values(), default=0)
        positions = {}
        for index, node in enumerate(leaves):
            positions[node["index"]] = (
                3.25 + 3.0 * (index + 0.5) / len(leaves),
                2.8 - 4.3 * depths[node["index"]] / max(1, maximum),
                0,
            )
        for node in reversed(list(nodes.values())):
            if node["chord"] != NONE:
                positions[node["index"]] = (
                    (positions[node["left"]][0] + positions[node["right"]][0]) / 2,
                    2.8 - 4.3 * depths[node["index"]] / max(1, maximum),
                    0,
                )
        return positions

    def centroid_begin(self, event):
        self.hide_panel()
        tree = self.data.trees[event["seq"]]
        tree.positions = self.centroid_positions(tree)
        self.panels[event["seq"]] = VGroup()
        self.glyphs[event["seq"]] = {}
        self.panel_edges[event["seq"]] = {}
        self.panel_id = event["seq"]
        self.panel = self.panels[event["seq"]]
        self.scene.graph.set_opacity(0.18)

    def centroid_node(self, event):
        tree = self.data.trees[event["tree"]]
        point = tree.positions[event["index"]]
        if event["chord"] == NONE:
            tint = self.scene.region_colors(tree.event["owner"], tree.submap.regions)[
                event["region"]
            ]
            glyph = Dot(point, radius=0.045, color=tint)
            source = self.region(tree.submap, tree.event["owner"], event["region"])[1]
        else:
            glyph = VGroup(
                Circle(radius=0.08, color=BLUE, stroke_width=1.5),
                Line((-0.07, 0, 0), (0.07, 0, 0), color=BLUE, stroke_width=2),
            ).move_to(point)
            record = next(chord for chord in tree.submap.chords if chord["id"] == event["chord"])
            source = self.chord(self.scene.viewport, record, GOLD)
            self.pulse(
                self.region(tree.submap, tree.event["owner"], event["centroid"]), duration=0.1
            )
        self.panel.add(glyph)
        self.glyphs[event["tree"]][event["index"]] = glyph
        animations = [TransformFromCopy(source, glyph)]
        if event["parent"] != NONE:
            edge = Line(tree.positions[event["parent"]], point, color=MUTED, stroke_width=1)
            self.panel.add(edge)
            self.panel_edges[event["tree"]][event["index"]] = edge
            animations.append(Create(edge))
        self.scene.play(*animations, run_time=0.22)

    def search_faces(self, event):
        identity = event["structure"]
        structure = self.data.structures[identity]
        self.hide_panel()
        graph = nx.Graph()
        graph.add_nodes_from(range(len(structure.regions)))
        graph.add_edges_from(structure.edges)
        planar, embedding = nx.check_planarity(graph)
        assert planar, "[C91 section 3.4] The compressed adjacency graph G is planar."
        positions = nx.planar_layout(embedding, scale=1.35, center=(4.75, 0.85))
        colors = self.scene.region_colors(structure.event["owner"], structure.submap.regions)
        self.panels[identity] = VGroup()
        self.glyphs[identity] = {}
        self.panel_edges[identity] = []
        self.panel_id, self.panel = identity, self.panels[identity]
        self.scene.graph.set_opacity(0.12)
        for face, region in enumerate(event["regions"]):
            glyph = VGroup(
                Square(side_length=0.14, color=colors[region], fill_opacity=0.5)
            ).move_to((*positions[face], 0))
            self.glyphs[identity][face] = glyph
            self.panel.add(glyph)
            marker = self.region(structure.submap, structure.event["owner"], region)[1]
            self.scene.play(TransformFromCopy(marker, glyph), run_time=0.15)

    def search_graph_edge(self, event):
        identity = event["structure"]
        glyphs = self.glyphs[identity]
        edge = Line(
            glyphs[event["first"]].get_center(),
            glyphs[event["second"]].get_center(),
            color=MUTED,
            stroke_width=1,
        )
        self.panels[identity].add(edge)
        self.panel_edges[identity].append((event["first"], event["second"], edge))
        model = self.data.structures[identity].submap
        if event["chord"] != NONE:
            record = next(chord for chord in model.chords if chord["id"] == event["chord"])
            source = self.chord(self.scene.viewport, record, GOLD)
        else:
            source = VGroup(
                self.arc(model, event["first_arc"], GREEN),
                self.arc(model, event["second_arc"], BLUE),
            )
        self.pulse(source, duration=0.1)
        moving = edge.copy()
        self.scene.play(TransformFromCopy(source, moving), run_time=0.25)
        self.scene.discard(moving)
        self.scene.add(edge)

    def separator_partition(self, event):
        identity = event["structure"]
        glyphs = self.glyphs[identity]
        animations = []
        colors = (GREEN, BLUE, GOLD)
        for face, part in zip(event["faces"], event["parts"], strict=True):
            glyph = glyphs[face]
            ring = Circle(radius=0.13, color=colors[part], stroke_width=2).move_to(glyph)
            glyph.add(ring)
            animations.append(Create(ring))
        self.scene.play(*animations, run_time=0.3)
        self.scene.wait(0.2)

    def separator_leaf(self, event):
        members = VGroup(*(self.glyphs[event["structure"]][face] for face in event["faces"]))
        frame = VGroup(
            *(
                SurroundingRectangle(
                    member, color=self.palette(event["subset"]), buff=0.05, stroke_width=1.2
                )
                for member in members
            )
        )
        self.panels[event["structure"]].add(frame)
        self.scene.play(Create(frame), run_time=0.25)

    def search_crossing(self, event):
        limit = exact_limit(event["y"])
        viewport = self.scene.viewport
        height = (
            limit.value
            if limit.value is not None
            else (viewport.top if limit.infinity > 0 else viewport.bottom)
        )
        height = max(viewport.bottom, min(viewport.top, height))
        point = viewport.point((viewport.right, height))
        point = (6.7, point[1], 0)
        tick = Line((6.55, point[1], 0), (6.85, point[1], 0), color=BLUE, stroke_width=1.5)
        self.panels[event["structure"]].add(tick)
        structure = self.data.structures[event["structure"]]
        if not hasattr(structure, "ticks"):
            structure.ticks = {}
        structure.ticks[event["index"]] = tick
        self.scene.play(
            TransformFromCopy(self.chord(self.scene.viewport, event), tick), run_time=0.2
        )

    def search_ready(self, event):
        structure = self.data.structures[event["structure"]]
        if not structure.crossings:
            return
        points = [structure.ticks[index].get_center() for index in sorted(structure.ticks)]
        line = Line(
            points[0] + [0, -0.25, 0], points[-1] + [0, 0.25, 0], color=MUTED, stroke_width=1
        )
        self.panels[event["structure"]].add(line)
        self.scene.play(Create(line), run_time=0.2)

    def search_begin(self, event):
        identity = event["seq"]
        self.query_stack.append(identity)
        self.query_contexts[identity] = None
        owner = event.get("owner")
        if event["structure"] != NONE:
            structure = self.data.structures[event["structure"]]
            owner = structure.event["owner"]
            if self.scene.query_curve != event["curve"]:
                self.query_contexts[identity] = self.save_context()
                self.display_map(structure.submap, owner, event["curve"])
            self.show_panel(event["structure"])
        elif self.scene.query_curve != event["curve"]:
            self.query_contexts[identity] = self.save_context()
            if owner is not None:
                self.display_map(self.owner_map(owner), owner, event["curve"])
            else:
                self.scene.clear(self.scene.geometry, self.scene.graph, self.panel, duration=0.15)
                self.scene.geometry, self.scene.graph, self.panel = VGroup(), VGroup(), VGroup()
                self.panel_id = None
                self.scene.active_drawing = None
                self.scene.show_curve(event["curve"])
        self.scene.show_curve(event["curve"])
        source = Dot(
            self.scene.viewport.symbolic_point(event["point"]), radius=0.07, color=GOLD
        ).set_z_index(8)
        self.query_sources[identity] = source
        self.scene.play(FadeIn(source), run_time=0.12)

    def search_edge(self, event):
        query = self.data.queries[event["query"]]
        curve = self.scene.trace.curves[query["curve"]]
        edge = event["edge"]
        segments = self.scene.viewport.edge_segments(
            curve[edge], curve[edge + 1], curve.edge_wrap(edge)
        )
        shape = VGroup(
            *(
                Line(
                    self.scene.viewport.point(a),
                    self.scene.viewport.point(b),
                    color=GOLD,
                    stroke_width=4,
                )
                for a, b in segments
                if a != b
            )
        )
        self.pulse(shape, duration=0.08)

    def search_scan(self, event):
        query = self.data.queries[event["query"]]
        structure = self.data.structures[query["structure"]]
        region = structure.regions[event["face"]]
        shape = self.region(structure.submap, structure.event["owner"], region)
        glyph = self.glyphs[query["structure"]][event["face"]]
        connector = Line(
            shape[1].get_center(), glyph.get_center(), color=GOLD, stroke_width=1
        ).set_opacity(0.35)
        self.pulse(VGroup(shape, glyph, connector), duration=0.1)

    def search_arc(self, event):
        query = self.data.queries[event["query"]]
        if query["structure"] != NONE:
            model = self.data.structures[query["structure"]].submap
        elif "owner" in query:
            model = self.owner_map(query["owner"])
        else:
            return
        self.pulse(self.arc(model, event["arc"]), duration=0.1)

    def search_piece(self, event):
        query = self.data.queries[event["query"]]
        curve = self.scene.trace.curves[query["curve"]]
        components = []
        for edge in range(
            min(event["first"], event["last"]), max(event["first"], event["last"]) + 1
        ):
            for a, b in self.scene.viewport.edge_segments(
                curve[edge], curve[edge + 1], curve.edge_wrap(edge)
            ):
                components.append([self.scene.viewport.point(a), self.scene.viewport.point(b)])
        self.pulse(self.paths(components, GREEN, 4), duration=0.1)

    def search_candidate(self, event):
        point = self.scene.viewport.symbolic_point(event["point"])
        ring = Circle(
            radius=0.09, color=GREEN if event["accepted"] else RED, stroke_width=2
        ).move_to(point)
        self.pulse(ring, GREEN if event["accepted"] else RED, duration=0.1)

    def search_subset(self, event):
        query = self.data.queries[event["query"]]
        structure = self.data.structures[query["structure"]]
        selected = VGroup(
            *(
                self.glyphs[query["structure"]][face]
                for face, subset in enumerate(structure.subsets)
                if subset == event["subset"]
            )
        )
        self.pulse(selected, GREEN, 0.2)

    def vertical_search(self, event):
        query = self.data.queries[event["query"]]
        structure = self.data.structures[query["structure"]]
        if event["mid"] in getattr(structure, "ticks", {}):
            self.pulse(structure.ticks[event["mid"]], duration=0.15)
        x = 6.7
        y = self.query_sources[event["query"]].get_center()[1]
        marker = Dot((x, y, 0), radius=0.05, color=GOLD)
        self.pulse(marker, duration=0.1)

    def boundary_search(self, event):
        query = self.data.queries[event["query"]]
        structure = self.data.structures[query["structure"]]
        self.pulse(
            self.region(structure.submap, structure.event["owner"], event["region"]), duration=0.12
        )

    def boundary_identify(self, event):
        self.boundary_search(event)

    def search_end(self, event):
        assert self.query_stack.pop() == event["query"]
        self.scene.query(
            self.scene.trace.events[event["result"]], self.query_sources.pop(event["query"])
        )
        saved = self.query_contexts.pop(event["query"])
        if saved is not None:
            self.restore_context(saved)

    def centroid_visit(self, event):
        tree = self.data.trees[event["tree"]]
        if self.context is None:
            self.context = self.save_context()
        if self.scene.query_curve != tree.submap.metadata["curve"]:
            self.display_map(tree.submap, tree.event["owner"], tree.submap.metadata["curve"])
        self.show_panel(event["tree"])
        self.pulse(self.glyphs[event["tree"]][event["node"]], duration=0.15)
        node = tree.nodes[event["node"]]
        if node["chord"] != NONE:
            record = next(chord for chord in tree.submap.chords if chord["id"] == node["chord"])
            self.pulse(self.chord(self.scene.viewport, record), duration=0.15)
        else:
            self.pulse(self.region(tree.submap, tree.event["owner"], node["region"]), GREEN, 0.15)

    def centroid_empty_region(self, event):
        tree = self.data.trees[event["tree"]]
        self.pulse(self.region(tree.submap, tree.event["owner"], event["region"]), RED, 0.2)

    def centroid_reject(self, event):
        tree = self.data.trees[event["tree"]]
        self.pulse(self.arc(tree.submap, event["probe_arc"]), BLUE, 0.15)
        node = tree.nodes[event["node"]]
        rejected = node["right"] if event["go_left"] else node["left"]
        self.pulse(self.panel_edges[event["tree"]][rejected], RED, 0.2)

    def centroid_branch(self, event):
        if event["next"] != NONE:
            self.pulse(self.panel_edges[event["tree"]][event["next"]], GREEN, 0.15)

    def conformality_test(self, event):
        self.hide_panel()
        self.scene.clear(self.capacity)
        model = self.owner_map(event["owner"])
        slots = VGroup(
            *(
                Circle(radius=0.07, color=GREEN if index < event["arcs"] else MUTED).move_to(
                    (3.35 + 0.27 * index, 2.65, 0)
                )
                for index in range(4)
            )
        )
        overflow = VGroup(
            *(
                Dot((3.35 + 0.27 * index, 2.65, 0), radius=0.06, color=RED)
                for index in range(4, event["arcs"])
            )
        )
        self.capacity = VGroup(slots, overflow)
        self.scene.play(FadeIn(self.capacity), run_time=0.25)
        self.pulse(
            self.region(model, event["owner"], event["region"]), GREEN if event["accepted"] else RED
        )

    def arc_pair(self, event):
        model = self.owner_map(event["owner"])
        self.scene.clear(self.decision)
        self.decision = VGroup(
            self.arc(model, event["first"], GREEN), self.arc(model, event["second"], BLUE)
        )
        self.scene.play(Create(self.decision), run_time=0.3)

    def arc_pair_result(self, event):
        if self.context is not None:
            self.restore_context(self.context)
            self.context = None
        self.pulse(self.decision, GREEN if event["found"] else RED)
        self.scene.clear(self.decision, duration=0.15)
        self.decision = VGroup()

    def granularity_test(self, event):
        self.hide_panel()
        model = self.owner_map(event["owner"])
        record = next(chord for chord in model.chords if chord["id"] == event["chord"])
        self.pulse(
            self.chord(self.scene.viewport, record), GREEN if event["accepted"] else RED, 0.1
        )

    def refinement_need(self, event):
        self.hide_panel()
        model = self.owner_map(event["owner"])
        for arc, weight in zip(event["arcs"], event["weights"], strict=True):
            self.pulse(self.arc(model, arc), RED if weight > event["limit"] else GREEN, 0.12)

    def refinement_begin(self, event):
        saved = self.save_context()
        model = self.owner_map(event["owner"])
        parent = self.region(model, event["owner"], event["region"])
        ghost = parent.copy().set_opacity(0.5).scale(0.22).move_to((-5.4, 2.65, 0))
        frame = SurroundingRectangle(ghost, buff=0.12, color=GOLD, stroke_width=1)
        self.scene.play(TransformFromCopy(parent, ghost), Create(frame), run_time=0.5)
        self.refinements[event["seq"]] = {
            "saved": saved,
            "parent": parent,
            "ghost": VGroup(ghost, frame),
            "event": event,
        }
        self.refinement_stack.append(event["seq"])

    def refinement_boundary(self, event):
        refinement = self.refinements[event["refinement"]]
        target = self.scene.curve_path(event["curve"], GOLD)
        self.scene.play(TransformFromCopy(refinement["parent"], target), run_time=0.6)
        self.scene.clear(target, duration=0.2)

    def refinement_map(self, event):
        refinement = self.refinements[event["refinement"]]
        refinement["auxiliary"] = event["auxiliary"]
        self.display_map(self.owner_map(event["auxiliary"]), event["auxiliary"], event["curve"])

    def refinement_discard(self, event):
        self.pulse(self.chord(self.scene.viewport, event), RED, 0.12)

    def refinement_extract(self, event):
        refinement = self.refinements[event["refinement"]]
        auxiliary = self.owner_map(refinement["auxiliary"])
        source = next(chord for chord in auxiliary.chords if chord["id"] == event["chord"])
        target = self.chord(self.scene.viewport, event, GREEN)
        self.scene.play(
            TransformFromCopy(self.chord(self.scene.viewport, source, GOLD), target), run_time=0.5
        )
        self.scene.clear(target, duration=0.2)

    def refinement_end(self, event):
        assert self.refinement_stack.pop() == event["refinement"]
        refinement = self.refinements[event["refinement"]]
        self.restore_context(refinement["saved"])
        self.scene.clear(refinement["ghost"], duration=0.25)

    def settle(self, current):
        assert self.pending_change is not None
        self.hide_panel()
        event = self.pending_change
        target = self.scene.drawing(current)
        tint = GREEN if event["kind"] == "insert_chord" else RED
        changed = self.chord(self.scene.viewport, event, tint).set_z_index(6)
        self.scene.play(Create(changed), run_time=0.3)
        previous = self.scene.active_drawing
        assert previous is not None and previous.submap.metadata["map"] == current.identity
        animations = [FadeOut(changed)]
        for name in ("nodes", "edges", "markers", "faces", "arcs", "chords"):
            before, after = getattr(previous, name), getattr(target, name)
            for identity in before.keys() & after.keys():
                animations.append(Transform(before[identity], after[identity].copy()))
            for identity in before.keys() - after.keys():
                if name in {"nodes", "markers"} and event["kind"] == "remove_chord":
                    animations.append(
                        Transform(
                            before[identity], after[event["regions"][0]].copy().set_opacity(0)
                        )
                    )
                else:
                    animations.append(FadeOut(before[identity]))
            for identity in after.keys() - before.keys():
                animations.append(FadeIn(after[identity]))
        self.scene.play(*animations, run_time=0.7)
        self.scene.discard(self.scene.geometry)
        self.scene.discard(self.scene.graph)
        self.scene.geometry, self.scene.graph = target.geometry, target.graph
        self.scene.nodes, self.scene.edges = target.nodes, target.edges
        self.scene.chord_shapes = target.chords
        self.scene.active_drawing = target
        self.scene.add(target.geometry, target.graph)
        if self.scene.ring is not None:
            self.scene.remove(self.scene.ring)
            self.scene.ring = None
        self.pending_change = None

    def relabel(self, event, current):
        drawing = self.scene.active_drawing
        if drawing is None or drawing.submap.metadata["map"] != current.identity:
            return
        if event["kind"] == "reindex":
            for name in ("nodes", "markers", "faces"):
                setattr(
                    drawing,
                    name,
                    {
                        event["regions"][identity]: shape
                        for identity, shape in getattr(drawing, name).items()
                    },
                )
            for name in ("edges", "chords"):
                setattr(
                    drawing,
                    name,
                    {
                        event["chords"][identity]: shape
                        for identity, shape in getattr(drawing, name).items()
                    },
                )
            drawing.positions = {
                event["regions"][identity]: position
                for identity, position in drawing.positions.items()
            }
            drawing.colors = self.scene.map_colors(current)
        drawing.arcs = {event["arcs"][identity]: shape for identity, shape in drawing.arcs.items()}
        drawing.submap = current.submap(self.scene.trace.curves)
        assert set(drawing.nodes) == current.regions and set(drawing.chords) == set(current.chords)

    def checkpoint(self, event):
        if event["name"] == "trapezoids":
            self.hide_panel()
            self.scene.clear(self.scene.graph, self.scene.shelf, duration=0.35)
            self.scene.graph = VGroup()
            self.scene.geometry.fade(0.65)
        elif event["name"] == "unimonotone":
            self.scene.geometry.fade(1 - 0.18 / 0.35)
            self.scene.overlay.fade(0.65)
        elif event["name"] == "triangulate":
            self.scene.overlay.fade(1 - 0.12 / 0.35)
            self.piece_outlines.fade(0.75)
            if self.cursor is not None:
                self.scene.clear(self.cursor, duration=0.15)
            self.cursor = None
            self.current_vertex = None

    def trapezoid(self, event):
        self.scene.output_event(event)
        self.trapezoids[len(self.trapezoids)] = self.scene.overlay[-1]

    def partition_vertex(self, event):
        self.move_cursor(event["vertex"])

    def partition_return(self, event):
        self.pulse(
            Line(self.scene.points[event["top"]], self.scene.points[event["bottom"]], color=BLUE),
            BLUE,
            0.18,
        )
        self.move_cursor(event["bottom"])

    def trapezoid_test(self, event):
        self.pulse(self.trapezoids[event["trapezoid"]], GREEN if event["split"] else MUTED, 0.22)
        self.pulse(
            Dot(self.scene.points[event["bottom"]], radius=0.07),
            GREEN if event["split"] else MUTED,
            0.15,
        )

    def piece(self, event):
        shape = Polygon(
            *(self.scene.points[index] for index in event["vertices"]),
            color=BLUE,
            stroke_width=2,
            fill_opacity=0,
        )
        self.piece_outlines.add(shape)
        self.scene.play(Create(shape), run_time=0.35)

    def move_cursor(self, vertex, backward=False):
        target = self.scene.points[vertex]
        if self.cursor is None:
            self.cursor = (
                Circle(radius=0.1, color=GOLD, stroke_width=2).move_to(target).set_z_index(8)
            )
            self.scene.play(Create(self.cursor), run_time=0.2)
        else:
            if self.current_vertex is not None and self.current_vertex != vertex:
                arrow = Arrow(
                    self.scene.points[self.current_vertex],
                    target,
                    buff=0.06,
                    color=BLUE if backward else GOLD,
                    stroke_width=2,
                    max_tip_length_to_length_ratio=0.18,
                )
                self.scene.play(Create(arrow), self.cursor.animate.move_to(target), run_time=0.35)
                self.scene.clear(arrow, duration=0.12)
            else:
                self.cursor.move_to(target)
        self.current_vertex = vertex

    def triangle_piece_begin(self, event):
        self.active_piece = event["seq"]
        self.remaining = list(event["vertices"])
        self.piece_boundary = Polygon(
            *(self.scene.points[index] for index in self.remaining),
            color=BLUE,
            stroke_width=3,
            fill_opacity=0,
        )
        self.scene.play(Create(self.piece_boundary), run_time=0.45)
        for key in ("start", "end"):
            self.pulse(
                Circle(radius=0.09, color=BLUE).move_to(self.scene.points[event[key]]), BLUE, 0.16
            )
        current = self.remaining[(self.remaining.index(event["start"]) + 1) % len(self.remaining)]
        self.move_cursor(current)

    def convexity_test(self, event):
        self.scene.clear(self.test, duration=0.12)
        points = [self.scene.points[event[key]] for key in ("previous", "current", "next")]
        tint = GREEN if event["convex"] else RED
        self.test = VGroup(
            Line(points[0], points[1], color=tint, stroke_width=5),
            Line(points[1], points[2], color=tint, stroke_width=5),
            Line(points[0], points[2], color=tint, stroke_width=1.2),
            *(Dot(point, radius=0.06, color=tint) for point in points),
        )
        self.scene.play(Create(self.test), run_time=0.35)
        self.scene.wait(0.2)

    def triangle(self, event):
        tint = self.palette(self.active_piece)
        triangle = Polygon(
            *(self.scene.points[index] for index in event["vertices"]),
            color=tint,
            fill_color=tint,
            fill_opacity=0.25,
            stroke_width=2,
        )
        self.triangles.add(triangle)
        self.scene.play(TransformFromCopy(self.test, triangle), run_time=0.5)

    def vertex_remove(self, event):
        self.remaining.remove(event["vertex"])
        if len(self.remaining) >= 3:
            target = Polygon(
                *(self.scene.points[index] for index in self.remaining),
                color=BLUE,
                stroke_width=3,
                fill_opacity=0,
            )
        else:
            target = Line(
                *(self.scene.points[index] for index in self.remaining), color=BLUE, stroke_width=3
            )
        self.scene.play(Transform(self.piece_boundary, target), run_time=0.45)
        self.scene.clear(self.test, duration=0.15)
        self.test = VGroup()

    def triangle_cursor(self, event):
        self.move_cursor(event["vertex"], bool(event["backward"]))
        self.scene.clear(self.test, duration=0.12)
        self.test = VGroup()

    def triangle_piece_end(self, event):
        self.scene.clear(self.piece_boundary, self.cursor, duration=0.2)
        self.piece_boundary = VGroup()
        self.cursor = None
        self.current_vertex = None
        self.active_piece = None

    def apply(self, event):
        if event["kind"] == "checkpoint" and event["name"] not in {
            "trapezoids",
            "unimonotone",
            "triangulate",
        }:
            return False
        handler = getattr(self, event["kind"], None)
        if handler is not None:
            handler(event)
            return True
        return event["kind"] == "ray" and event["seq"] in self.data.results
