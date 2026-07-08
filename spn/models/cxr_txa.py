from ..core import PetriNet, Marking

import numpy as np

p21 = 200.0   # Xist scaling
p22 = 500.0   # cXR scaling
p23 = 400.0   # tXA scaling

p3 = 3.0      # Xist --| tXA Hill coefficient
p4 = 0.3     # Xist --| tXA threshold

p5 = 3.0      # Xist --| cXR Hill coefficient
p6 = 0.3     # Xist --| cXR threshold

p11 = 3.0     # tXA -> Xist Hill coefficient
p12 = 0.40   # tXA -> Xist threshold

p13 = 3.0     # cXR --| Xist Hill coefficient
p14 = 0.55    # cXR --| Xist threshold

XIST_DEGRADATION = 0.1733


def safe_hill(x: float, n: float, theta: float) -> float:
    x = max(float(x), 0.0)
    theta = max(float(theta), 1e-12)
    return x**n / (x**n + theta**n)


def f_cxr(cxr: float) -> float:
    return 1.0 - safe_hill(cxr, p13, p22 * p14)


def f_txa(txa1: float, txa2: float) -> float:
    txa = 0.5 * (txa1 + txa2)
    return safe_hill(txa, p11, p23 * p12)


def g_txa(x: float) -> float:
    return 1.0 - safe_hill(x, p3, p21 * p4)


def g_cxr(x: float) -> float:
    return 1.0 - safe_hill(x, p5, p21 * p6)


def build_cxr_txa_net() -> PetriNet:
    net = PetriNet()

    for name in [
        "x1", "!x1",
        "x2", "!x2",
        "txa1", "!txa1",
        "txa2", "!txa2",
        "cxr1", "!cxr1",
        "cxr2", "!cxr2",
    ]:
        net.add_place(name)

    X1 = net.get_place("x1")
    NX1 = net.get_place("!x1")
    X2 = net.get_place("x2")
    NX2 = net.get_place("!x2")

    TXA1 = net.get_place("txa1")
    TXA2 = net.get_place("txa2")

    CXR1 = net.get_place("cxr1")
    NCXR1 = net.get_place("!cxr1")
    CXR2 = net.get_place("cxr2")
    NCXR2 = net.get_place("!cxr2")

    # t1: Xist1 production
    net.add_transition(
        name="t1_x1_prod",
        pre={"!x1": 1, "txa1": 1, "txa2": 1, "!cxr1": 1},
        post={"x1": 1, "txa1": 1, "txa2": 1, "!cxr1": 1},
        guard_fn=lambda m: (
            m[CXR1] < p22 * p14 and
            m[TXA1] + m[TXA2] > 2 * p23 * p12
        ),
        propensity_fn=lambda m: p21 * f_cxr(m[CXR1]) * f_txa(m[TXA1], m[TXA2]),
    )

    # t2: Xist1 degradation
    net.add_transition(
        name="t2_x1_deg",
        pre={"x1": 1},
        post={"!x1": 1},
        propensity_fn=lambda m: XIST_DEGRADATION * m[X1],
    )

    # t3: tXA1 production
    net.add_transition(
        name="t3_txa1_prod",
        pre={"!x1": 1, "!txa1": 1},
        post={"!x1": 1, "txa1": 1},
        guard_fn=lambda m: m[X1] < p21 * p4,
        propensity_fn=lambda m: p23 * g_txa(m[X1]),
    )

    # t4: tXA1 degradation
    net.add_transition(
        name="t4_txa1_deg",
        pre={"txa1": 1},
        post={"!txa1": 1},
        propensity_fn=lambda m: m[TXA1],
    )

    # t5: cXR1 production
    net.add_transition(
        name="t5_cxr1_prod",
        pre={"!x1": 1, "!cxr1": 1},
        post={"!x1": 1, "cxr1": 1},
        guard_fn=lambda m: m[X1] < p21 * p6,
        propensity_fn=lambda m: p22 * g_cxr(m[X1]),
    )

    # t6: cXR1 degradation
    net.add_transition(
        name="t6_cxr1_deg",
        pre={"cxr1": 1},
        post={"!cxr1": 1},
        propensity_fn=lambda m: m[CXR1],
    )

    # t7: Xist2 production
    net.add_transition(
        name="t7_x2_prod",
        pre={"!x2": 1, "txa1": 1, "txa2": 1, "!cxr2": 1},
        post={"x2": 1, "txa1": 1, "txa2": 1, "!cxr2": 1},
        guard_fn=lambda m: (
            m[CXR2] < p22 * p14 and
            m[TXA1] + m[TXA2] > 2 * p23 * p12
        ),
        propensity_fn=lambda m: p21 * f_cxr(m[CXR2]) * f_txa(m[TXA1], m[TXA2]),
    )

    # t8: Xist2 degradation
    net.add_transition(
        name="t8_x2_deg",
        pre={"x2": 1},
        post={"!x2": 1},
        propensity_fn=lambda m: XIST_DEGRADATION * m[X2],
    )

    # t9: tXA2 production
    net.add_transition(
        name="t9_txa2_prod",
        pre={"!x2": 1, "!txa2": 1},
        post={"!x2": 1, "txa2": 1},
        guard_fn=lambda m: m[X2] < p21 * p4,
        propensity_fn=lambda m: p23 * g_txa(m[X2]),
    )

    # t10: tXA2 degradation
    net.add_transition(
        name="t10_txa2_deg",
        pre={"txa2": 1},
        post={"!txa2": 1},
        propensity_fn=lambda m: m[TXA2],
    )

    # t11: cXR2 production
    net.add_transition(
        name="t11_cxr2_prod",
        pre={"!x2": 1, "!cxr2": 1},
        post={"!x2": 1, "cxr2": 1},
        guard_fn=lambda m: m[X2] < p21 * p6,
        propensity_fn=lambda m: p22 * g_cxr(m[X2]),
    )

    # t12: cXR2 degradation
    net.add_transition(
        name="t12_cxr2_deg",
        pre={"cxr2": 1},
        post={"!cxr2": 1},
        propensity_fn=lambda m: m[CXR2],
    )

    return net

def initial_marking(net: PetriNet) -> Marking:
    return net.marking_from_dict(
        {
            "x1": 0,
            "!x1": 500,
            "x2": 0,
            "!x2": 500,

            "txa1": 250,
            "!txa1": 0,
            "txa2": 250,
            "!txa2": 0,

            "cxr1": 250,
            "!cxr1": 0,
            "cxr2": 250,
            "!cxr2": 0,
        }
    )