from dataclasses import dataclass
from typing import List

import numpy as np

@dataclass
class SimulationResult:
    times: np.ndarray
    markings: np.ndarray
    fired_transitions: List[str]