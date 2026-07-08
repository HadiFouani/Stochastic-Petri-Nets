from spn.visualizer import PetriNetVisualizer
from spn.models.cxr_txa import build_cxr_txa_net, initial_marking

net = build_cxr_txa_net()
m0 = initial_marking(net)
viz = PetriNetVisualizer(net)
viz.draw(m0)