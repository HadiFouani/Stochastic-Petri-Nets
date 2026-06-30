import matplotlib.pyplot as plt

from .core import PetriNet
from .result import SimulationResult

def plot_result(net: PetriNet, result: SimulationResult) -> None:
    names = net.place_names()

    def col(place: str) -> int:
        return names.index(place)

    times = result.times
    M = result.markings

    plt.figure(figsize=(10, 5))
    plt.step(times, M[:, col("x1")], where="post", label="x1")
    plt.step(times, M[:, col("x2")], where="post", label="x2")
    plt.xlabel("time (h)")
    plt.ylabel("molecule count")
    plt.title("Xist dynamics")
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.step(times, M[:, col("txa1")], where="post", label="txa1")
    plt.step(times, M[:, col("txa2")], where="post", label="txa2")
    plt.step(times, M[:, col("cxr1")], where="post", label="cxr1")
    plt.step(times, M[:, col("cxr2")], where="post", label="cxr2")
    plt.xlabel("time (h)")
    plt.ylabel("molecule count")
    plt.title("Regulator dynamics")
    plt.legend()
    plt.tight_layout()
    plt.show()