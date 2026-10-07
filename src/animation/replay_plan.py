def same_submap(first, second):
    def records(items):
        return [
            {key: value for key, value in item.items() if key not in {"seq", "map"}}
            for item in items
        ]

    return (
        first.metadata["curve"] == second.metadata["curve"]
        and first.metadata["root"] == second.metadata["root"]
        and first.regions == second.regions
        and records(first.chords) == records(second.chords)
        and records(first.arcs) == records(second.arcs)
    )


class ReplayPlan:
    def __init__(self, trace):
        self.detailed_builds = set()
        self.build_entries = {}
        self.rays = set()
        self.decisions = set()
        self.chains = {}
        self.chain_snapshots = {}
        self.chain_children = {}
        self.trees = {event["tree"] for event in trace.events if event["kind"] == "centroid_visit"}
        queries = []
        groups = {}
        group = None
        limit = 0
        seen_decisions = set()
        chain = None
        for event in trace.events:
            kind = event["kind"]
            if kind == "chain":
                chain = (event["grade"], event["index"])
                assert chain not in self.chains
                assert event["last"] - event["first"] == 1 << event["grade"], (
                    "[C91 section 4] A grade-lambda chain contains 2^lambda edges."
                )
                self.chains[chain] = event
                group = event["seq"]
                limit = 4 if event["grade"] == 0 and event["index"] == 0 else 0
            elif kind == "map_begin" and event["name"] == "canonical":
                assert chain is not None and event["curve"] == self.chains[chain]["curve"]
                self.chain_snapshots[event["seq"]] = chain
                grade, index = chain
                children = ((grade - 1, 2 * index), (grade - 1, 2 * index + 1)) if grade else ()
                self.chain_children[chain] = children
                if children:
                    first, second = (self.chains[child] for child in children)
                    parent = self.chains[chain]
                    assert (
                        parent["first"] == first["first"]
                        and first["last"] == second["first"]
                        and second["last"] == parent["last"]
                    ), "[C91 section 4] Each chain is the union of its two contiguous children."
            elif kind in {"merge", "refinement_begin"} or (
                kind == "checkpoint" and event["name"] in {"down_phase", "bounded_regions"}
            ):
                group, limit = event["seq"], 4
            elif kind == "checkpoint" and event["name"] == "fusion":
                group, limit = event["seq"], 2
            elif kind == "search_begin":
                queries.append(event["seq"])
            elif kind == "search_end":
                assert queries.pop() == event["query"]
            elif kind == "ray" and limit and len(queries) <= 1:
                groups.setdefault(group, (limit, []))[1].append((event, len(queries)))
            elif kind == "build_begin":
                self.build_entries[event["map"]] = []
                if not self.detailed_builds:
                    self.detailed_builds.add(event["map"])
            elif kind == "build_enter":
                self.build_entries[event["map"]].append(event)
            elif kind == "granularity_test" and event["eligible"]:
                decision = bool(event["accepted"])
                if decision not in seen_decisions:
                    seen_decisions.add(decision)
                    self.decisions.add(event["seq"])
        assert not queries
        assert len(self.chain_snapshots) == len(self.chains)
        for limit, candidates in groups.values():
            standalone = [event for event, depth in candidates if depth == 0]
            candidates = standalone or [event for event, _ in candidates]
            chosen = []
            for predicate in (
                lambda event: event["hit"] and not event["wrapped"],
                lambda event: event["hit"] and event["wrapped"],
                lambda event: not event["hit"],
            ):
                if len(chosen) == limit:
                    break
                match = next((event for event in candidates if predicate(event)), None)
                if match is not None:
                    chosen.append(match["seq"])
            for event in candidates:
                if len(chosen) == limit:
                    break
                if event["seq"] not in chosen:
                    chosen.append(event["seq"])
            self.rays.update(chosen)

    def shows(self, event):
        kind = event["kind"]
        if kind == "ray":
            return event["seq"] in self.rays
        if kind.startswith(("search_", "separator_")) or kind in {
            "vertical_search",
            "boundary_search",
            "boundary_identify",
            "fusion_cursor",
        }:
            return False
        if kind == "granularity_test":
            return event["seq"] in self.decisions
        if kind == "conformality_test":
            return not event["accepted"]
        if kind == "centroid_begin":
            return event["seq"] in self.trees
        if kind in {"centroid_node", "centroid_end"}:
            return event["tree"] in self.trees
        if kind.startswith("build_"):
            return kind == "build_end" or event["map"] in self.detailed_builds
        return True
