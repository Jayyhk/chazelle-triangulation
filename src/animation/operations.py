from copy import deepcopy

from trace_data import Submap, chord_identity


class WorkingSubmap:
    def __init__(self, identity, curve):
        self.identity = identity
        self.curve = curve
        self.regions = {0}
        self.chords = {}
        self.arcs = {}
        self.walk = [0]
        self.entered = []
        self.complete = False
        self.root = 0

    def submap(self, curves):
        exits = {region: [] for region in self.regions}
        for chord in self.chords.values():
            for region in chord["regions"]:
                exits[region].append(chord)
        records = [
            {
                "kind": "map_region",
                "id": region,
                "weight": max(
                    (arc["edge_count"] for arc in self.arcs.values() if arc["region"] == region),
                    default=0,
                ),
                "bounded": bool(exits[region])
                and all(chord["null"] or not chord["infinite"] for chord in exits[region]),
            }
            for region in sorted(self.regions)
        ]
        records.extend(dict(chord, kind="map_chord") for _, chord in sorted(self.chords.items()))
        records.extend(dict(arc, kind="map_arc") for _, arc in sorted(self.arcs.items()))
        metadata = {
            "first": 0,
            "last": len(curves[self.curve]) - 1,
            "root": self.root,
            "map": self.identity,
            "curve": self.curve,
        }
        return Submap(metadata, records, curves[self.curve], 0, conformal=False)


class Operations:
    def __init__(self, trace):
        self.trace = trace
        self.maps = {}
        self.active = None
        self.snapshot = None
        self.verified = 0
        self.fusion_maps = set()
        self.fusion_invalidated = set()
        self.discovered_chords = {}

    def apply(self, event):
        kind = event["kind"]
        if "map" in event and kind not in {"map_begin", "fusion_remove"}:
            self.active = event["map"]
        if kind == "visibility_chord":
            self.discovered_chords.setdefault(event["curve"], set()).add(chord_identity(event))
        elif kind == "build_begin":
            self.active = event["map"]
            assert self.active not in self.maps
            self.maps[self.active] = WorkingSubmap(self.active, event["curve"])
        elif kind == "copy_submap":
            self.active = event["map"]
            assert self.active not in self.maps
            self.maps[self.active] = deepcopy(self.maps[event["source"]])
            self.maps[self.active].identity = self.active
            self.maps[self.active].curve = event["curve"]
        elif kind == "merge_inputs":
            self.fusion_maps = {event["first"], event["second"]}
            self.fusion_invalidated = set()
            assert all(self.maps[identity].complete for identity in self.fusion_maps)
        elif kind == "fusion_remove":
            identity = (event["map"], event["id"])
            assert event["map"] in self.fusion_maps
            assert identity not in self.fusion_invalidated
            assert chord_identity(self.maps[event["map"]].chords[event["id"]]) == chord_identity(
                event
            )
            self.fusion_invalidated.add(identity)
        elif kind == "build_enter":
            current = self.maps[event["map"]]
            assert current.walk[-1] == event["parent"]
            assert event["region"] not in current.regions
            current.regions.add(event["region"])
            current.entered.append(event)
            current.walk.append(event["region"])
        elif kind == "build_leave":
            current = self.maps[event["map"]]
            assert current.walk.pop() == event["region"]
            assert current.walk[-1] == event["parent"]
        elif kind == "build_walk":
            assert self.maps[event["map"]].walk[-1] == event["region"]
        elif kind in {"build_chord", "insert_chord"}:
            current = self.maps[event["map"]]
            current.chords[event["id"]] = dict(event)
            if kind == "insert_chord":
                current.regions.update(event["regions"])
        elif kind == "build_arc":
            self.maps[event["map"]].arcs[event["id"]] = dict(event)
        elif kind == "arc_owner":
            self.maps[event["map"]].arcs[event["arc"]]["region"] = event["region"]
        elif kind == "delete_arc":
            del self.maps[event["map"]].arcs[event["arc"]]
        elif kind == "remove_chord":
            self.active = event["map"]
            current = self.maps[self.active]
            chord = current.chords.pop(event["id"])
            assert chord_identity(chord) == chord_identity(event)
            keep, removed = event["regions"]
            current.regions.remove(removed)
            if current.root == removed:
                current.root = keep
            for chord in current.chords.values():
                chord["regions"] = [
                    keep if region == removed else region for region in chord["regions"]
                ]
                if chord["left_region"] == removed:
                    chord["left_region"] = keep
            for arc in current.arcs.values():
                if arc["region"] == removed:
                    arc["region"] = keep
        elif kind in {"reindex", "reindex_arcs"}:
            current = self.maps[event["map"]]
            current.arcs = {
                event["arcs"][index]: dict(arc, id=event["arcs"][index])
                for index, arc in current.arcs.items()
            }
            if kind == "reindex":
                current.root = event["regions"][current.root]
                current.regions = {event["regions"][region] for region in current.regions}
                current.chords = {
                    event["chords"][index]: dict(
                        chord,
                        id=event["chords"][index],
                        regions=[event["regions"][region] for region in chord["regions"]],
                        left_region=event["regions"][chord["left_region"]],
                    )
                    for index, chord in current.chords.items()
                }
                for arc in current.arcs.values():
                    arc["region"] = event["regions"][arc["region"]]
        elif kind == "build_end":
            current = self.maps[event["map"]]
            assert current.walk == [0]
            assert len(current.entered) == len(current.chords)
            assert sorted(
                (chord_identity(chord), chord["regions"]) for chord in current.chords.values()
            ) == sorted(
                (chord_identity(entered), [entered["parent"], entered["region"]])
                for entered in current.entered
            )
            current.complete = True
            current.root = event["root"]
            if current.curve in self.discovered_chords:
                assert {
                    chord_identity(chord) for chord in current.chords.values()
                } == self.discovered_chords.pop(current.curve), (
                    "[C91 section 2.1] The builder uses exactly the discovered visibility chords."
                )
            current.submap(self.trace.curves)
        elif kind == "checkpoint" and event["name"] in {"contract_end", "split_end"}:
            current = self.maps[event["map"]]
            current.root = event["root"]
            current.submap(self.trace.curves)
        elif kind == "map_begin":
            self.snapshot = event
        elif kind == "map_end":
            recorded = self.trace.submaps[self.snapshot["seq"]]
            current = self.maps[self.snapshot["map"]]
            assert current.root == self.snapshot["root"], (event["seq"], "root")
            replayed = current.submap(self.trace.curves)
            assert {
                region: (record["weight"], record["bounded"])
                for region, record in replayed.regions.items()
            } == {
                region: (record["weight"], record["bounded"])
                for region, record in recorded.regions.items()
            }, (event["seq"], "regions")
            assert {
                chord["id"]: (chord_identity(chord), chord["regions"], chord["left_region"])
                for chord in replayed.chords
            } == {
                chord["id"]: (chord_identity(chord), chord["regions"], chord["left_region"])
                for chord in recorded.chords
            }, (event["seq"], "chords")
            assert {
                arc["id"]: tuple(
                    arc[key]
                    for key in (
                        "region",
                        "edge_count",
                        "start",
                        "end",
                        "ranges",
                        "traversal",
                        "parameter",
                        "ascending",
                        "start_tag",
                        "end_tag",
                        "wraps",
                    )
                )
                for arc in replayed.arcs
            } == {
                arc["id"]: tuple(
                    arc[key]
                    for key in (
                        "region",
                        "edge_count",
                        "start",
                        "end",
                        "ranges",
                        "traversal",
                        "parameter",
                        "ascending",
                        "start_tag",
                        "end_tag",
                        "wraps",
                    )
                )
                for arc in recorded.arcs
            }, (event["seq"], "arcs")
            self.verified += 1
            self.snapshot = None
        return self.maps.get(self.active)

    def verify(self):
        for event in self.trace.events:
            self.apply(event)
        assert self.verified == len(self.trace.submaps)
