from ..core import PetriNet, Marking

import numpy as np


# System-size parameter
N = 10


# ODE-scale parameters: these remain fixed
p21 = 347.0   # Xist maximal production rate
p22 = 76.0    # cXR maximal production rate
p23 = 79.0    # tXA maximal production rate

p3 = 3.4      # Xist --| tXA Hill coefficient
p4 = 0.017    # Xist --| tXA normalized threshold

p5 = 2.2      # Xist --| cXR Hill coefficient
p6 = 0.019    # Xist --| cXR normalized threshold

p11 = 2.7     # tXA -> Xist Hill coefficient
p12 = 1.03    # tXA -> Xist normalized threshold

p13 = 2.6     # cXR --| Xist Hill coefficient
p14 = 0.20    # cXR --| Xist normalized threshold

XIST_DEGRADATION = 0.1733


def safe_hill(x: float, n: float, theta: float) -> float:
    """
    Increasing Hill function:

        H(x) = x^n / (x^n + theta^n)

    Here x and theta are expressed at the normalized ODE scale,
    not at the integer token-count scale.
    """
    x = max(float(x), 0.0)
    theta = max(float(theta), 1e-12)

    if x == 0.0:
        return 0.0

    # Numerically stable form, especially when N or the Hill exponent is large.
    ratio = theta / x
    return 1.0 / (1.0 + ratio**n)


def f_cxr(cxr_density: float) -> float:
    """
    Repression of Xist by cXR.

    cxr_density = cXR_count / N
    """
    return 1.0 - safe_hill(
        cxr_density,
        p13,
        p22 * p14,
    )


def f_txa(txa1_density: float, txa2_density: float) -> float:
    """
    Activation of Xist by the shared tXA pool.

    txa1_density = tXA1_count / N
    txa2_density = tXA2_count / N
    """
    txa_density = 0.5 * (txa1_density + txa2_density)

    return safe_hill(
        txa_density,
        p11,
        p23 * p12,
    )


def xist_production_share(
    cxr_count: float,
    txa1_count: float,
    txa2_count: float,
    source_txa_count: float,
) -> float:
    """
    Count-level propensity of one of the two split Xist-production
    transitions.

    The density-level production term is

        beta_x = p21 * f_cxr(cxr) * f_txa(txa1, txa2).

    The corresponding count-level intensity is

        a_x^N = N * beta_x.

    It is split proportionally between tXA1 and tXA2.
    """
    total_txa_count = txa1_count + txa2_count

    if total_txa_count <= 0.0:
        return 0.0

    # Convert integer counts to normalized ODE-scale quantities.
    cxr_density = cxr_count / N
    txa1_density = txa1_count / N
    txa2_density = txa2_count / N

    density_level_production = (
        p21
        * f_cxr(cxr_density)
        * f_txa(txa1_density, txa2_density)
    )

    source_fraction = source_txa_count / total_txa_count

    # Count-level propensity: N * beta(x)
    return N * source_fraction * density_level_production


def g_txa(x_density: float) -> float:
    """
    Repression of tXA production by Xist.

    x_density = Xist_count / N
    """
    return 1.0 - safe_hill(
        x_density,
        p3,
        p21 * p4,
    )


def g_cxr(x_density: float) -> float:
    """
    Repression of cXR production by Xist.

    x_density = Xist_count / N
    """
    return 1.0 - safe_hill(
        x_density,
        p5,
        p21 * p6,
    )


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

    # ---------------------------------------------------------
    # Allele 1
    # ---------------------------------------------------------

    # Xist1 production attributed to tXA1
    net.add_transition(
        name="t_prod_1_x1",
        pre={"!x1": 1, "txa1": 1, "!cxr1": 1},
        post={"x1": 1, "txa1": 1, "!cxr1": 1},
        propensity_fn=lambda m: xist_production_share(
            cxr_count=m[CXR1],
            txa1_count=m[TXA1],
            txa2_count=m[TXA2],
            source_txa_count=m[TXA1],
        ),
    )

    # Xist1 production attributed to tXA2
    net.add_transition(
        name="t_prod_2_x1",
        pre={"!x1": 1, "txa2": 1, "!cxr1": 1},
        post={"x1": 1, "txa2": 1, "!cxr1": 1},
        propensity_fn=lambda m: xist_production_share(
            cxr_count=m[CXR1],
            txa1_count=m[TXA1],
            txa2_count=m[TXA2],
            source_txa_count=m[TXA2],
        ),
    )

    # Xist1 degradation
    #
    # delta * X1_count
    # = N * delta * (X1_count / N)
    net.add_transition(
        name="t2_x1_deg",
        pre={"x1": 1},
        post={"!x1": 1},
        propensity_fn=lambda m: XIST_DEGRADATION * m[X1],
    )

    # tXA1 production
    #
    # Count-level intensity:
    # N * p23 * g_txa(X1_count / N)
    net.add_transition(
        name="t3_txa1_prod",
        pre={"!x1": 1, "!txa1": 1},
        post={"!x1": 1, "txa1": 1},
        propensity_fn=lambda m: (
            N
            * p23
            * g_txa(m[X1] / N)
        ),
    )

    # tXA1 degradationa
    net.add_transition(
        name="t4_txa1_deg",
        pre={"txa1": 1},
        post={"!txa1": 1},
        propensity_fn=lambda m: m[TXA1],
    )

    # cXR1 production
    #
    # Count-level intensity:
    # N * p22 * g_cxr(X1_count / N)
    net.add_transition(
        name="t5_cxr1_prod",
        pre={"!x1": 1, "!cxr1": 1},
        post={"!x1": 1, "cxr1": 1},
        propensity_fn=lambda m: (
            N
            * p22
            * g_cxr(m[X1] / N)
        ),
    )

    # cXR1 degradation
    net.add_transition(
        name="t6_cxr1_deg",
        pre={"cxr1": 1},
        post={"!cxr1": 1},
        propensity_fn=lambda m: m[CXR1],
    )

    # ---------------------------------------------------------
    # Allele 2
    # ---------------------------------------------------------

    # Xist2 production attributed to tXA1
    net.add_transition(
        name="t_prod_1_x2",
        pre={"!x2": 1, "txa1": 1, "!cxr2": 1},
        post={"x2": 1, "txa1": 1, "!cxr2": 1},
        propensity_fn=lambda m: xist_production_share(
            cxr_count=m[CXR2],
            txa1_count=m[TXA1],
            txa2_count=m[TXA2],
            source_txa_count=m[TXA1],
        ),
    )

    # Xist2 production attributed to tXA2
    net.add_transition(
        name="t_prod_2_x2",
        pre={"!x2": 1, "txa2": 1, "!cxr2": 1},
        post={"x2": 1, "txa2": 1, "!cxr2": 1},
        propensity_fn=lambda m: xist_production_share(
            cxr_count=m[CXR2],
            txa1_count=m[TXA1],
            txa2_count=m[TXA2],
            source_txa_count=m[TXA2],
        ),
    )

    # Xist2 degradation
    net.add_transition(
        name="t8_x2_deg",
        pre={"x2": 1},
        post={"!x2": 1},
        propensity_fn=lambda m: XIST_DEGRADATION * m[X2],
    )

    # tXA2 production
    net.add_transition(
        name="t9_txa2_prod",
        pre={"!x2": 1, "!txa2": 1},
        post={"!x2": 1, "txa2": 1},
        propensity_fn=lambda m: (
            N
            * p23
            * g_txa(m[X2] / N)
        ),
    )

    # tXA2 degradation
    net.add_transition(
        name="t10_txa2_deg",
        pre={"txa2": 1},
        post={"!txa2": 1},
        propensity_fn=lambda m: m[TXA2],
    )

    # cXR2 production
    net.add_transition(
        name="t11_cxr2_prod",
        pre={"!x2": 1, "!cxr2": 1},
        post={"!x2": 1, "cxr2": 1},
        propensity_fn=lambda m: (
            N
            * p22
            * g_cxr(m[X2] / N)
        ),
    )

    # cXR2 degradation
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
            "!x1": int(round(N * p21)),

            "x2": 0,
            "!x2": int(round(N * p21)),

            "txa1": int(round(N * p23)),
            "!txa1": 0,

            "txa2": int(round(N * p23)),
            "!txa2": 0,

            "cxr1": int(round(N * p22)),
            "!cxr1": 0,

            "cxr2": int(round(N * p22)),
            "!cxr2": 0,
        }
    )