from spn.simulator import Simulator
from spn.plotting import plot_result
from spn.models.cxr_txa import build_cxr_txa_net, initial_marking
from spn.visualizer import PetriNetVisualizer

if __name__ == "__main__":
    net = build_cxr_txa_net()
    m0 = initial_marking(net)

    simulator = Simulator(net, seed=7)

    result = simulator.run(
        initial_marking=m0,
        t_max=400.0,
        max_steps=500_000,
        record_every_step=True,
        show_current_transition=True,
    )

    print("Places:")
    print(net.place_names())

    print("\nInitial marking:")
    print(m0)

    print("\nFinal marking:")
    print(result.markings[-1])

    print("\nNumber of fired transitions:")
    print(len(result.fired_transitions))

    print("\nLast 20 fired transitions:")
    print(result.fired_transitions[-20:])

    plot_result(net, result)
    
    final_marking = result.markings[-1]
    viz = PetriNetVisualizer(net)
    viz.draw(final_marking)
