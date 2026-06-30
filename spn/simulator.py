from typing import Optional

import numpy as np

from .core import PetriNet, Marking, TransitionId
from .result import SimulationResult

class Simulator:
    def __init__(self, net: PetriNet, seed: Optional[int] = None):
        self.net = net
        self.rng = np.random.default_rng(seed)

    def compute_propensities(self, marking: Marking) -> np.ndarray:
        return np.array(
            [t.propensity(marking) for t in self.net.transitions],
            dtype=float,
        )

    def choose_transition(self, propensities: np.ndarray, total: float) -> TransitionId:
        u = self.rng.random() * total
        cumulative = 0.0

        for i, a in enumerate(propensities):
            cumulative += a
            if cumulative >= u:
                return i

        return len(propensities) - 1

    def fire(self, marking: Marking, transition_id: TransitionId) -> None:
        transition = self.net.transitions[transition_id]

        for p, delta in transition.delta.items():
            marking[p] += delta

        if np.any(marking < 0):
            raise RuntimeError(
                f"Negative marking after firing transition {transition.name}"
            )

    def run(
        self,
        initial_marking: Marking,
        t_max: float,
        max_steps: int = 100_000,
        record_every_step: bool = True,
    ) -> SimulationResult:
        time = 0.0
        marking = initial_marking.copy()

        times = [time]
        markings = [marking.copy()]
        fired = []

        for _ in range(max_steps):
            propensities = self.compute_propensities(marking)
            total = float(np.sum(propensities))

            if total <= 0.0:
                break

            r1 = self.rng.random()
            tau = -np.log(r1) / total

            if time + tau > t_max:
                break

            transition_id = self.choose_transition(propensities, total)

            time += tau
            self.fire(marking, transition_id)

            fired.append(self.net.transitions[transition_id].name)

            if record_every_step:
                times.append(time)
                markings.append(marking.copy())

        return SimulationResult(
            times=np.array(times),
            markings=np.array(markings),
            fired_transitions=fired,
        )