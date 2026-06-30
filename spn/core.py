from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

import numpy as np

PlaceId = int
TransitionId = int
Marking = np.ndarray

PropensityFn = Callable[[Marking], float]
GuardFn = Callable[[Marking], bool]

@dataclass
class Place:
    name: str

@dataclass
class Transition:
    name: str
    pre: Dict[PlaceId, int]
    post: Dict[PlaceId, int]
    propensity_fn: PropensityFn
    guard_fn: Optional[GuardFn] = None
    delta: Dict[PlaceId, int] = field(init=False)

    def __post_init__(self):
        self.delta = {}
        all_places = set(self.pre) | set(self.post)

        for p in all_places:
            d = self.post.get(p, 0) - self.pre.get(p, 0)
            if d != 0:
                self.delta[p] = d

    def structurally_enabled(self, marking: Marking) -> bool:
        return all(marking[p] >= w for p, w in self.pre.items())

    def guard_enabled(self, marking: Marking) -> bool:
        if self.guard_fn is None:
            return True
        return bool(self.guard_fn(marking))

    def is_enabled(self, marking: Marking) -> bool:
        return self.structurally_enabled(marking) and self.guard_enabled(marking)

    def propensity(self, marking: Marking) -> float:
        if not self.is_enabled(marking):
            return 0.0

        value = float(self.propensity_fn(marking))

        if value < 0:
            raise ValueError(f"Negative propensity in transition {self.name}: {value}")

        return value

class PetriNet:
    def __init__(self):
        self.places: List[Place] = []
        self.transitions: List[Transition] = []
        self.place_index: Dict[str, PlaceId] = {}

    def add_place(self, name: str) -> PlaceId:
        if name in self.place_index:
            return self.place_index[name]

        pid = len(self.places)
        self.places.append(Place(name))
        self.place_index[name] = pid
        return pid

    def get_place(self, name: str) -> PlaceId:
        return self.place_index[name]

    def add_transition(
        self,
        name: str,
        pre: Dict[str, int],
        post: Dict[str, int],
        propensity_fn: PropensityFn,
        guard_fn: Optional[GuardFn] = None,
    ) -> TransitionId:
        pre_ids = {self.add_place(p): w for p, w in pre.items()}
        post_ids = {self.add_place(p): w for p, w in post.items()}

        tid = len(self.transitions)

        self.transitions.append(
            Transition(
                name=name,
                pre=pre_ids,
                post=post_ids,
                propensity_fn=propensity_fn,
                guard_fn=guard_fn,
            )
        )

        return tid

    def empty_marking(self) -> Marking:
        return np.zeros(len(self.places), dtype=np.int64)

    def marking_from_dict(self, values: Dict[str, int]) -> Marking:
        marking = self.empty_marking()

        for name, value in values.items():
            pid = self.add_place(name)

            if pid >= len(marking):
                new_marking = np.zeros(len(self.places), dtype=np.int64)
                new_marking[: len(marking)] = marking
                marking = new_marking

            marking[pid] = value

        return marking

    def place_names(self) -> List[str]:
        return [p.name for p in self.places]