from dataclasses import dataclass, field

from operations import Operations
from trace_data import Submap

NONE = (1 << 64) - 1


@dataclass
class SearchStructure:
    event: dict
    submap: Submap
    regions: list = field(default_factory=list)
    edges: list = field(default_factory=list)
    crossings: list = field(default_factory=list)
    subsets: list = field(default_factory=list)
    nodes: dict = field(default_factory=dict)


class ReplaySteps:
    def __init__(self, trace):
        self.trace = trace
        self.structures = {}
        self.trees = {}
        self.queries = {}
        self.results = {}
        self.pieces = {}
        operations = Operations(trace)
        queries = []
        refinements = []
        active_piece = None
        tested = None
        emitted = False
        for event in trace.events:
            operations.apply(event)
            kind = event["kind"]
            identity = event["seq"]
            if kind in {"search_structure", "centroid_begin"}:
                model = operations.maps[event["owner"]].submap(trace.curves)
                structure = SearchStructure(event, model)
                if kind == "search_structure":
                    self.structures[identity] = structure
                else:
                    self.trees[identity] = structure
            elif kind == "centroid_node":
                tree = self.trees[event["tree"]]
                assert event["index"] not in tree.nodes
                tree.nodes[event["index"]] = event
            elif kind == "centroid_end":
                tree = self.trees[event["tree"]]
                assert len(tree.nodes) == 2 * len(tree.submap.chords) + 1
                assert sorted(tree.nodes) == list(range(len(tree.nodes)))
                leaves = set()
                chords = set()
                for node in tree.nodes.values():
                    if node["chord"] == NONE:
                        assert node["size"] == 1 and node["region"] in tree.submap.regions
                        leaves.add(node["region"])
                    else:
                        chords.add(node["chord"])
                        left, right = (tree.nodes[node[key]] for key in ("left", "right"))
                        assert left["parent"] == right["parent"] == node["index"]
                        assert left["size"] + right["size"] == node["size"]
                        maximum = node["size"] - 1 - (node["size"] + 2) // 4
                        assert max(left["size"], right["size"]) - 1 <= maximum, (
                            "[C91 section 2.3] A centroid edge balances the remaining edges."
                        )
                assert leaves == set(tree.submap.regions)
                assert chords == {chord["id"] for chord in tree.submap.chords}
            elif kind == "search_faces":
                self.structures[event["structure"]].regions = event["regions"]
            elif kind == "search_graph_edge":
                self.structures[event["structure"]].edges.append((event["first"], event["second"]))
            elif kind == "search_crossing":
                structure = self.structures[event["structure"]]
                assert event["index"] == len(structure.crossings)
                structure.crossings.append(event)
            elif kind == "separator_partition":
                structure = self.structures[event["structure"]]
                partition = dict(zip(event["faces"], event["parts"], strict=True))
                assert len(partition) == len(event["faces"])
                size = len(partition)
                assert all(part in {0, 1, 2} for part in partition.values())
                assert all(3 * list(partition.values()).count(part) <= 2 * size for part in (0, 1))
                assert list(partition.values()).count(2) ** 2 <= 8 * size
                assert all(
                    {partition.get(a), partition.get(b)} != {0, 1} for a, b in structure.edges
                ), "[C91 section 3.4] No separator edge joins A and B."
            elif kind == "separator_leaf":
                structure = self.structures[event["structure"]]
                assert len(event["faces"]) ** 3 <= len(structure.regions) ** 2
            elif kind == "search_ready":
                structure = self.structures[event["structure"]]
                structure.subsets = event["subsets"]
                assert len(structure.regions) == 1 or len(structure.subsets) == len(
                    structure.regions
                )
                for a, b in zip(structure.crossings, structure.crossings[1:], strict=False):
                    assert a["above"] == b["below"]
            elif kind == "search_begin":
                assert event["structure"] == NONE or event["structure"] in self.structures
                self.queries[identity] = event
                queries.append(identity)
            elif kind in {
                "search_edge",
                "search_arc",
                "search_piece",
                "search_scan",
                "search_candidate",
                "vertical_search",
                "boundary_search",
                "boundary_identify",
                "search_subset",
            }:
                assert queries and event["query"] == queries[-1], (
                    "Search visits belong to the active query."
                )
                query = self.queries[event["query"]]
                if kind == "search_edge":
                    assert 0 <= event["edge"] < len(trace.curves[query["curve"]]) - 1
                if kind == "vertical_search":
                    structure = self.structures[query["structure"]]
                    assert (
                        0 <= event["lo"] <= event["mid"] < event["hi"] <= len(structure.crossings)
                    )
                if kind == "search_scan":
                    structure = self.structures[query["structure"]]
                    assert (structure.subsets[event["face"]] == NONE) == bool(event["separator"])
            elif kind == "search_end":
                assert queries.pop() == event["query"]
                result = trace.events[event["result"]]
                query = self.queries[event["query"]]
                assert result["kind"] == "ray" and result["curve"] == query["curve"]
                assert (
                    result["origin"] == query["point"] and result["direction"] == query["direction"]
                )
                self.results[event["result"]] = event["query"]
            elif kind == "centroid_visit":
                assert event["node"] in self.trees[event["tree"]].nodes
            elif kind == "centroid_branch":
                node = self.trees[event["tree"]].nodes[event["node"]]
                assert event["found"] or event["next"] in {node["left"], node["right"]}
            elif kind == "granularity_test":
                assert bool(event["eligible"]) == (
                    min(event["first_degree"], event["second_degree"]) < 3
                )
                assert bool(event["accepted"]) == (
                    event["eligible"] and event["weight"] <= event["limit"]
                )
            elif kind == "conformality_test":
                assert bool(event["accepted"]) == (event["arcs"] <= 4)
            elif kind == "refinement_need":
                assert len(event["arcs"]) == len(event["weights"])
                assert bool(event["needed"]) == any(
                    weight > event["limit"] for weight in event["weights"]
                )
            elif kind == "refinement_begin":
                refinements.append(identity)
            elif kind in {
                "refinement_map",
                "refinement_boundary",
                "refinement_extract",
                "refinement_discard",
            }:
                assert refinements and event["refinement"] == refinements[-1]
            elif kind == "refinement_end":
                assert refinements.pop() == event["refinement"]
            elif kind == "triangle_piece_begin":
                assert active_piece is None
                active_piece = identity
                self.pieces[identity] = event
                remaining = list(event["vertices"])
                current = remaining[(remaining.index(event["start"]) + 1) % len(remaining)]
            elif kind == "convexity_test":
                assert active_piece == event["piece"] and event["current"] == current
                index = remaining.index(current)
                assert event["previous"] == remaining[index - 1]
                assert event["next"] == remaining[(index + 1) % len(remaining)]
                assert event["convex"] in (0, 1)
                tested = event
                emitted = False
            elif kind == "triangle":
                assert tested is not None and tested["convex"] and not emitted
                assert event["vertices"] == [tested[key] for key in ("previous", "current", "next")]
                emitted = True
            elif kind == "vertex_remove":
                assert active_piece == event["piece"] and emitted and tested["convex"]
                assert event["vertex"] == current
                assert event["vertex"] not in {
                    self.pieces[active_piece][key] for key in ("start", "end")
                }
                remaining.remove(event["vertex"])
                assert len(remaining) == event["remaining"]
            elif kind == "triangle_cursor":
                assert active_piece == event["piece"] and tested is not None
                assert bool(tested["convex"]) == emitted
                back = tested["convex"] and tested["previous"] != self.pieces[active_piece]["start"]
                expected = tested["previous"] if back else tested["next"]
                assert event["vertex"] == expected and bool(event["backward"]) == bool(back)
                current = expected
                tested = None
            elif kind == "triangle_piece_end":
                assert active_piece == event["piece"] and len(remaining) == 2
                active_piece = None
        assert not queries and not refinements and active_piece is None
        assert operations.verified == len(trace.submaps)
