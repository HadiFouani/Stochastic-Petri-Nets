import matplotlib.pyplot as plt

from spn.simulator import Simulator
from spn.models.cxr_txa import build_cxr_txa_net, initial_marking
from spn.visualizer import PetriNetVisualizer


NUM_RUNS = 10
T_MAX = 300.0
MAX_STEPS = 500000


if __name__ == "__main__":

    net = build_cxr_txa_net()
    m0 = initial_marking(net)

    plt.figure(figsize=(14, 7))

    for run in range(NUM_RUNS):

        simulator = Simulator(net, seed=run)

        result = simulator.run(
            initial_marking=m0,
            t_max=T_MAX,
            max_steps=MAX_STEPS,
            record_every_step=True,
        )

        times = result.times

        x1_index = net.place_names().index("x1")
        x2_index = net.place_names().index("x2")

        x1 = [m[x1_index] for m in result.markings]
        x2 = [m[x2_index] for m in result.markings]

        plt.plot(
            times,
            x1,
            label=f"x1_run{run+1}",
            linewidth=1.2,
        )

        plt.plot(
            times,
            x2,
            label=f"x2_run{run+1}",
            linewidth=1.2,
            linestyle="--",
        )

    plt.title("Multiple Gillespie simulations")
    plt.xlabel("time (h)")
    plt.ylabel("molecule count")

    plt.legend(fontsize=8, ncol=2)
    plt.grid(True)

    plt.show()