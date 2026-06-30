# Stochastic Petri Net Framework for Gene Regulatory Networks

A generic framework for modeling and simulating biological regulatory
networks using Stochastic Petri Nets and the Gillespie Stochastic
Simulation Algorithm.

The framework was initially developed for the stochastic simulation of
the cXR–tXA model of X chromosome inactivation proposed by Mutzel et al.,
but is designed to be extensible to arbitrary gene regulatory networks.

---

## Table of Contents

1. Introduction
2. Definitions
    - Ordinary Differential Equations
    - Petri Nets
    - Stochastic Petri Nets
3. Stochastic Simulation: Gillespie Algorithm
4. Theoretical Background
5. The cXR–tXA Model
6. Project Structure and Contribution Guidelines
7. Future Directions
8. References

# 1. Introduction

Gene Regulatory Networks (GRNs) describe interactions between genes,
RNAs, proteins and regulatory molecules responsible for controlling
cellular behaviour.

Traditionally, these systems are modeled using systems of Ordinary
Differential Equations (ODEs), where each variable represents the
concentration of a molecular species and evolves continuously in time.

For a system containing $n$ species

$$
x_1,x_2,\dots,x_n,
$$

the dynamics are generally written as

$$
\frac{dx_i}{dt}
=
f_i(x_1,\dots,x_n),
\qquad
i=1,\dots,n.
$$

ODE models are extremely successful for describing average behaviour
over large populations of cells. However, they neglect molecular noise
and stochastic fluctuations which become important whenever the number
of molecules involved is small.

Many biological phenomena emerge precisely because of this stochasticity.

Examples include:

- cell fate decisions,
- differentiation,
- genetic switches,
- phenotypic heterogeneity,
- stochastic symmetry breaking.

A particularly important example is **random X chromosome inactivation
(XCI)** in female mammals.

Initially, the two X chromosomes are equivalent and both exhibit low
levels of Xist expression:

$$
XaXa.
$$

During differentiation, stochastic fluctuations trigger the
up-regulation of Xist on one chromosome only, leading to the stable
mono-allelic state

$$
XaXi.
$$

Because the system is initially symmetric, deterministic models cannot
spontaneously choose one chromosome over the other:

$$
x_1(0)=x_2(0)
\quad\Longrightarrow\quad
x_1(t)=x_2(t)
\qquad
\forall t.
$$

Consequently, deterministic ODE models cannot reproduce the transition

$$
XaXa
\rightarrow
XaXi
$$

without introducing an explicit asymmetry.

To capture this phenomenon, stochastic models are required.

The objective of this project is therefore to transform deterministic
gene regulatory models into **Stochastic Petri Nets (SPNs)** whose
dynamics can be simulated using the **Gillespie Stochastic Simulation
Algorithm (SSA)**.

The framework was initially developed for the cXR–tXA model proposed by
Mutzel et al., but is intentionally designed to remain independent of
any particular biological system and may be reused for arbitrary
regulatory networks.

# 2. Definitions

## 2.1 Ordinary Differential Equations

An Ordinary Differential Equation (ODE) describes the evolution of a
continuous state variable as a function of time.

For a state variable

$$
x(t)\in\mathbb{R},
$$

an ODE has the general form

$$
\frac{dx}{dt}
=
f(x,t).
$$

For biological systems involving multiple interacting molecular species,
the model becomes a system of coupled ODEs

$$
\frac{d\mathbf{x}}{dt}
=
\mathbf{f}(\mathbf{x}),
$$

where

$$
\mathbf{x}
=
(x_1,x_2,\dots,x_n)^T
$$

is the vector of molecular concentrations.

The functions

$$
f_i
$$

typically contain:

- production terms,
- degradation terms,
- activation terms,
- repression terms.

For example,

$$
\frac{dx}{dt}
=
\alpha
-
\beta x
$$

represents constant production with first-order degradation.

ODE models assume:

- continuous concentrations,
- deterministic dynamics,
- infinite divisibility of molecular species,
- absence of intrinsic noise.

## 2.2 Petri Nets

A Petri Net is a mathematical formalism for describing distributed,
concurrent and event-driven systems.

A Petri net is defined by the tuple

$$
N=(P,T,Pre,Post),
$$

where

- $P$ is the set of places,
- $T$ is the set of transitions,
- $Pre$ is the pre-incidence function,
- $Post$ is the post-incidence function.

Places represent resources or states of the system.

Transitions represent events capable of modifying the state of the
system.

The current state of the system is represented by a marking

$$
M:P\rightarrow\mathbb{N},
$$

where

$$
M(p)
$$

denotes the number of tokens contained in place $p$.

---

### Transition Enabling

A transition

$$
t\in T
$$

is enabled if all required input tokens are available:

$$
M(p)
\ge
Pre(p,t)
\qquad
\forall p\in P.
$$

---

### Transition Firing

When transition

$$
t
$$

fires, the marking evolves according to

$$
M'
=
M
+
Post(:,t)
-
Pre(:,t).
$$

Defining

$$
\Delta_t
=
Post(:,t)-Pre(:,t),
$$

the firing rule becomes

$$
M'
=
M+\Delta_t.
$$

Petri nets naturally describe:

- concurrency,
- synchronization,
- causality,
- resource consumption,
- resource production.

## 2.3 Stochastic Petri Nets

A Stochastic Petri Net (SPN) extends a classical Petri net by assigning
a stochastic firing rate to each transition.

To every transition

$$
t_i
$$

is associated a propensity function

$$
a_i(M),
$$

which defines the instantaneous probability per unit time that the
transition fires when the system is in marking

$$
M.
$$

The collection of all transition propensities defines the total firing
rate

$$
a_0(M)
=
\sum_i a_i(M).
$$

Unlike deterministic ODE models, SPNs operate on:

- discrete molecule counts,
- random event times,
- stochastic transition selection.

Consequently, two simulations starting from identical initial conditions
may evolve toward different trajectories.

This property makes SPNs particularly suitable for studying:

- molecular noise,
- stochastic switching,
- cell fate decisions,
- symmetry breaking phenomena.

# 3. Stochastic Simulation: Gillespie Algorithm

Once the biological system has been represented as a Stochastic Petri
Net, the remaining problem is to determine:

- which transition fires next,
- and when this firing occurs.

The framework uses the **Stochastic Simulation Algorithm (SSA)**
introduced by Gillespie in 1977.

Unlike deterministic ODE solvers, which evolve all variables
simultaneously through infinitesimal time increments, Gillespie's
algorithm considers the system as a sequence of discrete random events
occurring in continuous time.

At any instant, only one transition fires.

---

## 3.1 Transition Propensities

Assume that the Petri net contains

$$
m
$$

transitions

$$
T=\{t_1,t_2,\dots,t_m\}.
$$

Each transition is associated with a propensity function

$$
a_i(M),
$$

which depends on the current marking

$$
M.
$$

The quantity

$$
a_i(M)dt
$$

represents the probability that transition

$$
t_i
$$

fires during the infinitesimal interval

$$
[t,t+dt).
$$

Therefore,

$$
a_i(M)
$$

can be interpreted as the instantaneous firing rate of transition

$$
t_i.
$$

Transitions that are not enabled have zero propensity:

$$
a_i(M)=0.
$$

---

## 3.2 Total Propensity

The total rate at which any event occurs is

$$
a_0(M)
=
\sum_{i=1}^{m}a_i(M).
$$

The waiting time until the next event therefore follows an exponential
distribution:

$$
\tau
\sim
Exp(a_0).
$$

whose probability density function is

$$
f(\tau)
=
a_0e^{-a_0\tau}.
$$

Its expected value is

$$
\mathbb{E}[\tau]
=
\frac{1}{a_0}.
$$

Consequently:

- large propensities imply short waiting times,
- small propensities imply long waiting times.

---

## 3.3 Sampling the Next Event Time

A random number

$$
r_1
\sim
U(0,1)
$$

is generated.

Using the inverse transform method, the waiting time is sampled as

$$
\tau
=
-\frac{\ln(r_1)}{a_0}.
$$

The simulation time is then updated:

$$
t
\leftarrow
t+\tau.
$$

---

## 3.4 Selecting the Transition

A second random number is generated:

$$
r_2
\sim
U(0,1).
$$

The interval

$$
[0,a_0]
$$

is partitioned according to the transition propensities:

$$
[0,a_1),
$$

$$
[a_1,a_1+a_2),
$$

$$
[a_1+a_2,a_1+a_2+a_3),
$$

$$
\dots
$$

The selected transition is the first transition satisfying

$$
\sum_{j=1}^{k}a_j
>
r_2a_0.
$$

Equivalently,

$$
k
=
\min
\left\{
i:
\sum_{j=1}^{i}a_j
>
r_2a_0
\right\}.
$$

This ensures that transition

$$
t_i
$$

is selected with probability

$$
P(t_i)
=
\frac{a_i}{a_0}.
$$

---

## 3.5 Marking Update

Once the transition

$$
t_k
$$

has been selected, the marking evolves according to

$$
M'
=
M+\Delta_k,
$$

where

$$
\Delta_k
=
Post(:,k)-Pre(:,k).
$$

The simulation then continues from the new state

$$
M'.
$$

---

## 3.6 Gillespie Algorithm

The complete algorithm can therefore be summarized as follows.

---

### Algorithm

Initialize:

$$
t=0,
\qquad
M=M_0.
$$

Repeat until termination:

1. Determine the set of enabled transitions.

2. Compute all propensities

$$
a_1,\dots,a_m.
$$

3. Compute

$$
a_0=\sum_i a_i.
$$

4. If

$$
a_0=0,
$$

terminate the simulation.

5. Generate

$$
r_1,r_2
\sim
U(0,1).
$$

6. Compute

$$
\tau
=
-\frac{\ln(r_1)}{a_0}.
$$

7. Select transition

$$
t_k
$$

such that

$$
\sum_{j=1}^{k}a_j
>
r_2a_0.
$$

8. Fire transition

$$
t_k.
$$

9. Update the marking

$$
M
\leftarrow
M+\Delta_k.
$$

10. Update time

$$
t
\leftarrow
t+\tau.
$$

11. Return to step 1.

---

## 3.7 Complexity

Assuming the Petri net contains

$$
m
$$

transitions, a single Gillespie iteration requires:

- evaluation of all propensities:

$$
O(m)
$$

- transition selection:

$$
O(m)
$$

- firing:

$$
O(|\Delta_k|)
$$

where

$$
|\Delta_k|
$$

denotes the number of places modified by transition

$$
t_k.
$$

The total complexity is therefore

$$
O(m)
$$

per simulation step.

For sparse biological networks, this complexity is generally dominated
by propensity evaluations.

---

## 3.8 Why Gillespie?

The Gillespie algorithm possesses several properties that make it
particularly suitable for gene regulatory networks:

- exact simulation of the underlying stochastic process,
- preservation of molecular discreteness,
- continuous time representation,
- natural handling of rare events,
- ability to reproduce stochastic switching and symmetry breaking.

These properties are essential for modeling random X chromosome
inactivation, where the transition

$$
XaXa
\rightarrow
XaXi
$$

is driven by stochastic fluctuations rather than deterministic dynamics.

# 4. Theoretical Background

The stochastic dynamics generated by a Stochastic Petri Net can be
interpreted mathematically as a **Continuous-Time Markov Chain (CTMC)**
defined over the set of reachable markings.

This section introduces the minimal stochastic process theory required
to understand this interpretation.

---

## 4.1 Random Variables

A random variable is a measurable mapping

$$
X:\Omega\rightarrow E
$$

from a probability space

$$
(\Omega,\mathcal F,\mathbb P)
$$

to a state space

$$
E.
$$

In the context of Stochastic Petri Nets, the state space is the set of
reachable markings:

$$
E=\mathcal M.
$$

A realization of the random variable therefore corresponds to a
particular marking of the Petri net.

---

## 4.2 Stochastic Processes

A stochastic process is a family of random variables indexed by time:

$$
(X_t)_{t\ge0}.
$$

For every fixed time

$$
t,
$$

the variable

$$
X_t
$$

describes the state of the system at time

$$
t.
$$

In our framework:

$$
X_t=M(t)
$$

where

$$
M(t)
$$

denotes the marking of the Petri net at time

$$
t.
$$

The stochastic simulation therefore generates a trajectory

$$
(M(t))_{t\ge0}.
$$

Unlike deterministic ODE models, two trajectories generated from the
same initial condition are not necessarily identical:

$$
M^{(1)}(t)
\neq
M^{(2)}(t).
$$

---

## 4.3 Markov Property

The stochastic process generated by the Petri net satisfies the Markov
property.

For any times

$$
0\le t_1<t_2<\dots<t_n<t
$$

we have

$$
\mathbb P
\left(
M(t+\Delta t)=m'
\mid
M(t)=m,
M(t_n),\dots,M(t_1)
\right)
=
\mathbb P
\left(
M(t+\Delta t)=m'
\mid
M(t)=m
\right).
$$

In other words:

> The future evolution of the system depends only on the current marking
> and not on the history of the process.

The marking therefore contains all information necessary to predict
future dynamics.

---

## 4.4 Continuous-Time Markov Chains

A Continuous-Time Markov Chain is a stochastic process

$$
(X_t)_{t\ge0}
$$

satisfying:

1. the Markov property,
2. a discrete state space,
3. continuous time evolution.

In our framework:

- time is continuous,

$$
t\in\mathbb R_+
$$

- states are discrete,

$$
M\in\mathbb N^{|P|}
$$

- transitions occur instantaneously.

The stochastic Petri net therefore naturally defines a CTMC.

---

## 4.5 State Space

The state space of the CTMC is the set of reachable markings:

$$
\mathcal S
=
\{M_0,M_1,M_2,\dots\}.
$$

Each state corresponds to one possible distribution of tokens among the
places of the Petri net.

For example:

$$
M_1=(100,0,250,250)
$$

and

$$
M_2=(101,0,249,250)
$$

correspond to two different states of the Markov chain.

---

## 4.6 Transition Rates

Suppose that firing transition

$$
t_k
$$

moves the system from marking

$$
M_i
$$

to marking

$$
M_j.
$$

The transition rate between these states is given by the propensity of
the transition:

$$
q_{ij}
=
a_k(M_i).
$$

If no transition connects the two markings:

$$
q_{ij}=0.
$$

---

## 4.7 Generator Matrix

The transition rates define the infinitesimal generator matrix

$$
Q=(q_{ij}).
$$

The off-diagonal entries are

$$
q_{ij}
=
a_k(M_i)
\qquad
(i\neq j)
$$

whenever transition

$$
t_k
$$

connects

$$
M_i
\rightarrow M_j.
$$

The diagonal coefficients are chosen such that each row sums to zero:

$$
q_{ii}
=
-
\sum_{j\neq i}
q_{ij}.
$$

Therefore,

$$
\sum_j q_{ij}=0.
$$

---

## 4.8 Exponential Waiting Times

An important property of CTMCs is that the waiting time before leaving a
state follows an exponential distribution.

If the current marking is

$$
M
$$

and the total propensity is

$$
a_0(M)
=
\sum_i a_i(M),
$$

then the waiting time satisfies

$$
\tau
\sim
Exp(a_0(M)).
$$

The corresponding density is

$$
f(\tau)
=
a_0e^{-a_0\tau}.
$$

Its expectation is

$$
\mathbb E[\tau]
=
\frac1{a_0}.
$$

This property is the mathematical justification behind Gillespie's
sampling formula

$$
\tau
=
-
\frac{\ln(r_1)}{a_0}.
$$

---

## 4.9 Jump Processes

Unlike ODE trajectories, which evolve continuously,

$$
x(t)
$$

is continuous in time,

the trajectories of a CTMC are piecewise constant:

$$
M(t)
=
M_k
\qquad
t_k\le t<t_{k+1}.
$$

The state changes only when a transition fires.

Consequently, the trajectories generated by the simulator have the form

$$
M_0
\rightarrow
M_1
\rightarrow
M_2
\rightarrow
\dots
$$

with random jump times

$$
0<t_1<t_2<t_3<\dots
$$

determined by the Gillespie algorithm.

---

## 4.10 Relation Between ODEs and CTMCs

The deterministic ODE model describes the average behaviour of an
infinite population of identical cells.

The stochastic Petri net instead describes the dynamics of a single
cell.

When molecule counts become very large, stochastic fluctuations become
small relative to the system size and the CTMC converges toward the
deterministic dynamics:

$$
\frac{M(t)}{N}
\rightarrow
x(t)
\qquad
N\rightarrow\infty.
$$

The deterministic ODE model can therefore be interpreted as the
large-population limit of the stochastic model.

This relationship provides the theoretical bridge between the original
ODE formulation and the stochastic Petri net implementation developed in
this project.

# 5. The cXR–tXA Model

## 5.1 Biological Motivation

The cXR–tXA model was introduced by Mutzel et al. to explain the
mechanisms underlying the initiation and maintenance of random X
chromosome inactivation (XCI).

The model aims to reproduce three essential biological properties:

- female-specific Xist activation,
- stable mono-allelic Xist expression,
- destabilization of both the fully active and fully inactive states.

The model identifies two regulatory components as sufficient to explain
these observations:

- a cis-acting Xist repressor:

$$
cXR
$$

- a trans-acting Xist activator:

$$
tXA
$$

The interaction between these regulators generates an extended
symmetric toggle switch capable of producing stochastic symmetry
breaking and stable mono-allelic expression. :contentReference[oaicite:0]{index=0}

---

## 5.2 Original Deterministic Model

The original model consists of six dynamical variables:

$$
x_1,\;x_2,\;
cXR_1,\;cXR_2,\;
tXA_1,\;tXA_2
$$

where:

| Variable | Meaning |
|----------|---------|
| $x_1$ | Xist level on chromosome 1 |
| $x_2$ | Xist level on chromosome 2 |
| $cXR_1$ | cis-repressor on chromosome 1 |
| $cXR_2$ | cis-repressor on chromosome 2 |
| $tXA_1$ | trans-activator on chromosome 1 |
| $tXA_2$ | trans-activator on chromosome 2 |

The deterministic model uses continuous concentrations normalized to
the interval

$$
[0,1].
$$

The stochastic formulation instead uses explicit molecule counts between
approximately 50 and 500 molecules per allele. :contentReference[oaicite:1]{index=1}

---

## 5.3 ODE System

### Xist Dynamics

For chromosome 1:

$$
\frac{dx_1}{dt}
=
p_{21}
f_{cXR}(cXR_1)
f_{tXA}(tXA_1,tXA_2)
-
\delta_x x_1
$$

For chromosome 2:

$$
\frac{dx_2}{dt}
=
p_{21}
f_{cXR}(cXR_2)
f_{tXA}(tXA_1,tXA_2)
-
\delta_x x_2
$$

where

$$
\delta_x = 0.1733
\;\text{h}^{-1}
$$

corresponds to the experimentally measured Xist degradation rate. :contentReference[oaicite:2]{index=2}

---

### cXR Repression Function

The repression exerted by cXR on Xist is modeled by

$$
f_{cXR}(cXR_i)
=
1-
\frac{
cXR_i^{p_{13}}
}{
cXR_i^{p_{13}}
+
(p_{22}p_{14})^{p_{13}}
}
$$

---

### tXA Activation Function

The trans-acting activator acts globally on both chromosomes:

$$
f_{tXA}(tXA_1,tXA_2)
=
\frac{
\left(
0.5(tXA_1+tXA_2)
\right)^{p_{11}}
}{
\left(
0.5(tXA_1+tXA_2)
\right)^{p_{11}}
+
(p_{23}p_{12})^{p_{11}}
}
$$

---

### tXA Dynamics

For chromosome 1:

$$
\frac{dtXA_1}{dt}
=
p_{23}
g_{tXA}(x_1)
-
tXA_1
$$

For chromosome 2:

$$
\frac{dtXA_2}{dt}
=
p_{23}
g_{tXA}(x_2)
-
tXA_2
$$

with

$$
g_{tXA}(x_i)
=
1-
\frac{x_i^{p_3}}
{x_i^{p_3}+(p_{21}p_4)^{p_3}}
$$

---

### cXR Dynamics

For chromosome 1:

$$
\frac{dcXR_1}{dt}
=
p_{22}
g_{cXR}(x_1)
-
cXR_1
$$

For chromosome 2:

$$
\frac{dcXR_2}{dt}
=
p_{22}
g_{cXR}(x_2)
-
cXR_2
$$

with

$$
g_{cXR}(x_i)
=
1-
\frac{x_i^{p_5}}
{x_i^{p_5}+(p_{21}p_6)^{p_5}}
$$

---

## 5.4 ODE Decomposition

Each differential equation was decomposed into:

- production terms,
- degradation terms.

For example,

$$
\frac{dx_1}{dt}
=
P_{x_1}
-
D_{x_1}
$$

with

$$
P_{x_1}
=
p_{21}
f_{cXR}(cXR_1)
f_{tXA}(tXA_1,tXA_2)
$$

and

$$
D_{x_1}
=
\delta_x x_1.
$$

The same decomposition was applied to all species.

This decomposition transforms the continuous model into a collection of
elementary stochastic reactions.

---

## 5.5 Species Definition

The stochastic Petri net contains twelve places:

| Place | Meaning |
|------|---------|
| x1 | Xist molecules on chromosome 1 |
| !x1 | absence of Xist on chromosome 1 |
| x2 | Xist molecules on chromosome 2 |
| !x2 | absence of Xist on chromosome 2 |
| txa1 | trans activator on chromosome 1 |
| !txa1 | inactive trans activator 1 |
| txa2 | trans activator on chromosome 2 |
| !txa2 | inactive trans activator 2 |
| cxr1 | cis repressor on chromosome 1 |
| !cxr1 | inactive cis repressor 1 |
| cxr2 | cis repressor on chromosome 2 |
| !cxr2 | inactive cis repressor 2 |

The complementary places

$$
!x_i,\;!txa_i,\;!cxr_i
$$

do not correspond to biological species.

They are auxiliary places introduced to represent inhibitory relations
within a purely positive Petri net formalism.

---

## 5.6 Transition Definition

The model contains twelve transitions:

| Transition | Description |
|-----------|------------|
| $t_1$ | Xist1 production |
| $t_2$ | Xist1 degradation |
| $t_3$ | tXA1 production |
| $t_4$ | tXA1 degradation |
| $t_5$ | cXR1 production |
| $t_6$ | cXR1 degradation |
| $t_7$ | Xist2 production |
| $t_8$ | Xist2 degradation |
| $t_9$ | tXA2 production |
| $t_{10}$ | tXA2 degradation |
| $t_{11}$ | cXR2 production |
| $t_{12}$ | cXR2 degradation |

---

## 5.7 Structural Enabling

A transition is structurally enabled if all required tokens are present:

$$
M(p)
\ge
Pre(p,t)
\qquad
\forall p\in P
$$

which corresponds to

```python
all(marking[p] >= w for p,w in transition.pre.items())
```
## 5.8 Biological Enabling Rules

Structural enabling alone is insufficient to represent biological
regulation.

In classical Petri nets, a transition is enabled if the required input
tokens are present:

$$
M(p) \geq Pre(p,t)
\qquad
\forall p \in P
$$

However, gene regulatory networks involve additional logical conditions
such as activation thresholds and repression mechanisms.

To represent these regulatory constraints, each transition may be
associated with an additional **guard function**.

A transition is therefore enabled only if both conditions hold:

$$
Enabled(t)
=
StructuralEnabled(t)
\land
GuardEnabled(t)
$$

or equivalently,

$$
Enabled(t)
=
\left(
M(p) \geq Pre(p,t)
\quad
\forall p \in P
\right)
\land
Guard(t,M)
$$

where

$$
Guard(t,M)
\in
\{0,1\}
$$

depends on the current marking.

---

### Xist Production

The production of Xist on chromosome 1 corresponds to transition

$$
t_1.
$$

This transition requires:

1. sufficient activation by the trans-acting activator,

$$
tXA_1 + tXA_2
>
2p_{23}p_{12}
$$

2. sufficiently weak repression by the cis-acting repressor,

$$
cXR_1
<
p_{22}p_{14}
$$

The guard condition is therefore

$$
Guard_{t_1}(M)
=
\left(
cXR_1 < p_{22}p_{14}
\right)
\land
\left(
tXA_1+tXA_2>2p_{23}p_{12}
\right)
$$

---

### tXA Production

The production of the trans-acting activator on chromosome 1
corresponds to transition

$$
t_3.
$$

Since Xist represses tXA expression, production is allowed only if the
Xist level remains below the repression threshold:

$$
x_1
<
p_{21}p_4
$$

which gives

$$
Guard_{t_3}(M)
=
\left(
x_1 < p_{21}p_4
\right)
$$

---

### cXR Production

The production of the cis-acting repressor on chromosome 1
corresponds to transition

$$
t_5.
$$

Similarly, Xist represses cXR expression and therefore

$$
x_1
<
p_{21}p_6
$$

must hold.

The guard condition becomes

$$
Guard_{t_5}(M)
=
\left(
x_1 < p_{21}p_6
\right)
$$

---

Equivalent conditions are defined for chromosome 2:

$$
Guard_{t_7},
\qquad
Guard_{t_9},
\qquad
Guard_{t_{11}}.
$$

---

## 5.9 Transition Propensities

The propensity functions are directly derived from the production and
degradation terms appearing in the original ODE system.

Consequently, the stochastic model preserves the same biological
regulatory mechanisms while introducing intrinsic noise through random
transition firing.

---

### Xist Production

The propensity associated with Xist production on chromosome 1 is

$$
a_1(M)
=
p_{21}
f_{cXR}(cXR_1)
f_{tXA}(tXA_1,tXA_2)
$$

---

### Xist Degradation

The propensity associated with Xist degradation on chromosome 1 is

$$
a_2(M)
=
\delta_x x_1
=
0.1733x_1
$$

---

### tXA Production

The propensity associated with tXA production on chromosome 1 is

$$
a_3(M)
=
p_{23}
g_{tXA}(x_1)
$$

where

$$
g_{tXA}(x_1)
=
1-
\frac{x_1^{p_3}}
{x_1^{p_3}+(p_{21}p_4)^{p_3}}
$$

---

### tXA Degradation

The propensity associated with tXA degradation on chromosome 1 is

$$
a_4(M)
=
tXA_1
$$

---

### cXR Production

The propensity associated with cXR production on chromosome 1 is

$$
a_5(M)
=
p_{22}
g_{cXR}(x_1)
$$

where

$$
g_{cXR}(x_1)
=
1-
\frac{x_1^{p_5}}
{x_1^{p_5}+(p_{21}p_6)^{p_5}}
$$

---

### cXR Degradation

The propensity associated with cXR degradation on chromosome 1 is

$$
a_6(M)
=
cXR_1
$$

---

The propensities for chromosome 2 are defined symmetrically:

$$
a_7(M),a_8(M),a_9(M),a_{10}(M),a_{11}(M),a_{12}(M)
$$

using the corresponding chromosome-specific variables.

---

## 5.10 Resulting Dynamics

The resulting Stochastic Petri Net preserves the regulatory structure of
the original ODE model while introducing stochastic fluctuations at the
molecular level.

The model is therefore able to reproduce:

- stochastic symmetry breaking,
- random mono-allelic Xist activation,
- stable XaXi states,
- instability of XaXa states,
- instability of XiXi states.

These properties emerge naturally from the stochastic dynamics and
cannot be reproduced by a deterministic symmetric ODE system without
introducing an explicit asymmetry between the two chromosomes.

# 6. Project Structure and Contribution Guidelines

This project is designed with two primary objectives:

1. provide a reusable and extensible framework for stochastic Petri nets,
2. keep biological models independent from the simulation engine.

The simulator itself should remain completely agnostic to the underlying
biological system.

The cXR–tXA model is therefore implemented as a particular instance of
the framework rather than being embedded inside the core simulation
logic.

---

## 6.1 Project Structure

The project is organized as follows:

```text
stochastic-petri-nets/
│
├── spn/
│   │
│   ├── __init__.py
│   │
│   ├── core.py
│   ├── simulator.py
│   ├── result.py
│   ├── plotting.py
│   │
│   └── models/
│       │
│       ├── __init__.py
│       └── cxr_txa.py
│
├── run_cxr_txa.py
│
├── requirements.txt
├── README.md
├── .gitignore
└── .venv/
```

---

## 6.2 Core Framework

The `spn` package contains the generic stochastic Petri net
implementation.

This code should remain independent of any biological assumptions.

### `core.py`

Contains the core Petri net data structures:

- `Place`
- `Transition`
- `PetriNet`

It also defines:

- markings,
- place identifiers,
- transition identifiers,
- propensity functions,
- guard functions.

This file defines the mathematical formalism of the Petri net.

---

### `simulator.py`

Contains the implementation of the Gillespie Stochastic Simulation
Algorithm.

Responsibilities:

- compute propensities,
- sample waiting times,
- select transitions,
- update markings,
- generate trajectories.

No biological logic should ever appear in this file.

---

### `result.py`

Defines the output object returned by simulations.

Currently this includes:

- simulation times,
- marking trajectories,
- transition history.

Additional statistics may be added in the future.

Examples include:

- occupancy probabilities,
- residence times,
- first passage times,
- stationary distributions.

---

### `plotting.py`

Contains visualization utilities.

Examples:

- molecular trajectories,
- state occupancy plots,
- histograms,
- phase diagrams,
- ensemble averages.

Plotting should remain separated from the simulation engine.

---

## 6.3 Biological Models

All biological models are placed inside:

```text
spn/models/
```

Each model should provide:

1. model parameters,
2. regulatory functions,
3. network construction,
4. initial conditions.

The simulator should never contain model-specific information.

---

### Example

The cXR–tXA model is implemented in:

```text
spn/models/cxr_txa.py
```

This file contains:

- Hill functions,
- model parameters,
- transition definitions,
- initial markings.

---

## 6.4 Separation of Concerns

The project follows the principle:

> One responsibility per file.

The simulator should not know anything about:

- X chromosome inactivation,
- genes,
- proteins,
- Hill functions,
- biological thresholds.

Similarly, biological models should not implement:

- transition selection,
- event scheduling,
- random number generation,
- simulation bookkeeping.

This separation simplifies:

- maintenance,
- debugging,
- extension,
- testing.

---

## 6.5 Adding a New Biological Model

To implement a new regulatory network, contributors should only create
a new file inside:

```text
spn/models/
```

For example:

```text
spn/models/toggle_switch.py
spn/models/repressilator.py
spn/models/lac_operon.py
```

The only required components are:

### Species definition

```python
net.add_place(...)
```

### Transition definition

```python
net.add_transition(...)
```

### Initial conditions

```python
def initial_marking(...):
```

The simulator itself should not require any modifications.

---

## 6.6 Coding Guidelines

Contributors are encouraged to follow the following conventions.

### Naming conventions

Places:

```text
x1
x2
txa1
cxr1
```

Complementary places:

```text
!x1
!txa1
!cxr1
```

Transitions:

```text
t1_x1_prod
t2_x1_deg
t3_txa1_prod
...
```

This naming convention makes trajectories and debugging significantly
easier.

---

### Type hints

All public functions should include explicit type annotations.

Example:

```python
def propensity(self, marking: Marking) -> float:
```

This improves:

- readability,
- IDE support,
- static analysis.

---

### Deterministic Random Seeds

When publishing figures or reproducing experiments, contributors are
encouraged to use fixed random seeds:

```python
Simulator(net, seed=42)
```

This improves reproducibility.

---

### Avoid Hard-Coded Constants

Biological parameters should never appear directly inside transition
definitions.

Instead:

```python
p21 = 300
p22 = 300
```

should be declared once at the beginning of the model file.

This greatly simplifies parameter sweeps and sensitivity analyses.

---

## 6.7 Testing Recommendations

Contributors are encouraged to test the following properties whenever
modifying the simulator.

### Structural consistency

Verify:

$$
M(p) \geq 0
\qquad
\forall p
$$

after every firing event.

---

### Transition enabling

Verify that disabled transitions always have:

$$
a_i(M)=0
$$

---

### Conservation relations

Whenever applicable, verify invariants such as

$$
x_i+\bar{x}_i=C
$$

or

$$
cXR_i+\overline{cXR_i}=C
$$

for constant molecule pools.

---

### Reproducibility

Two simulations initialized with the same random seed should produce
identical trajectories.

---

## 6.8 Version Control

Recommended workflow:

```bash
git checkout -b feature/new-model
```

Develop and test locally:

```bash
git add .
git commit -m "Add repressilator model"
```

Push changes:

```bash
git push origin feature/new-model
```

Finally create a pull request for review.

---

## 6.9 Dependency Management

The project uses a dedicated Python virtual environment:

```bash
python -m venv .venv
```

Activation:

### Windows PowerShell

```powershell
.venv\Scripts\Activate.ps1
```

### Linux / macOS

```bash
source .venv/bin/activate
```

Dependencies are installed using:

```bash
pip install -r requirements.txt
```

---

## 6.10 Long-Term Design Philosophy

The long-term objective is to evolve from:

```text
cXR–tXA model
```

towards:

```text
generic biological stochastic Petri net framework
```

The core simulator should therefore remain:

- generic,
- modular,
- extensible,
- biologically agnostic.

The biological models should be interchangeable components that can be
plugged into the simulator without requiring modifications to the
underlying engine.

# 7. Future Directions

This project is currently centered on the stochastic Petri net
implementation of the cXR–tXA model. However, several extensions are
natural and would make the framework more general, more visual, and
more useful for biological modeling.

---

## 7.1 Petri Net Graph Visualization

A useful next step is to implement automatic graph drawing for Petri
nets.

The goal is to generate a visual representation of any model defined
with the framework.

Places could be represented as circles:

```text
x1, x2, txa1, cxr1, ...
```

Transitions could be represented as boxes:

```text
t1_x1_prod, t2_x1_deg, ...
```

Arcs would be derived automatically from the `pre` and `post`
dictionaries of each transition.

Possible tools:

```text
networkx
matplotlib
graphviz
pygraphviz
pydot
```

A future function could have the form:

```python
draw_petri_net(net)
```

or

```python
export_petri_net_graph(net, filename="cxr_txa_net.png")
```

This would help contributors verify that the implemented Petri net
matches the intended biological model.

---

## 7.2 Exporting Incidence Matrices

Another useful extension is to export the structural matrices of the
Petri net:

$$
Pre
$$

$$
Post
$$

and

$$
C = Post - Pre.
$$

This would make it easier to compare the implemented model with the
mathematical formulation.

Possible outputs:

```text
CSV
NumPy arrays
LaTeX tables
Markdown tables
```

This would also help verify structural properties such as conservation
relations and reachability.

---

## 7.3 Connection with Most Permissive Boolean Networks

The project could be extended toward **Most Permissive Boolean Networks
(MPBNs)**.

Boolean Networks provide a qualitative abstraction of gene regulatory
networks, where each variable takes values in

$$
\{0,1\}.
$$

Most Permissive semantics enriches classical Boolean dynamics by
allowing intermediate permissive states, making it possible to represent
behaviors that may be missed by more restrictive update modes.

A promising direction is to study how SPNs, CPNs, and MPBNs relate to
each other.

One possible research path is:

```text
ODE model
    ↓
Piecewise-affine / PADE abstraction
    ↓
Boolean or multi-valued abstraction
    ↓
MPBN semantics
    ↓
Petri net encoding
    ↓
stochastic or continuous Petri net simulation
```

Recent work by Haar and Kolčák studies how Continuous Petri Nets can
faithfully represent Most Permissive Boolean Networks and create a
formal bridge between MPBNs and continuous refinements such as ODE
models. :contentReference[oaicite:0]{index=0}

---

## 7.4 Deeper CTMC Analysis

Currently, the framework uses CTMC theory mainly through Gillespie
simulation.

A future direction is to analyze the CTMC induced by the stochastic
Petri net more explicitly.

For a marking

$$
M
$$

and a transition

$$
t_i
$$

with propensity

$$
a_i(M),
$$

the CTMC jump rate is

$$
q(M,M+\Delta_i)=a_i(M).
$$

The generator matrix is then

$$
Q=(q_{ij})
$$

with

$$
q_{ii}=-\sum_{j\neq i}q_{ij}.
$$

Possible future analyses include:

- absorbing states,
- recurrent classes,
- transient states,
- expected hitting times,
- first-passage probabilities,
- stationary distributions,
- quasi-stationary distributions.

This would allow the project to move from pure simulation toward
formal stochastic analysis.

---

## 7.5 Relation Between SPNs and CTMCs

The stochastic Petri net implemented here naturally defines a CTMC over
the set of reachable markings.

Future documentation could make this relation explicit by constructing:

1. the reachable marking graph,
2. the transition rate graph,
3. the infinitesimal generator matrix,
4. the corresponding Kolmogorov forward equation.

The Kolmogorov forward equation has the form

$$
\frac{d}{dt}p(t)
=
p(t)Q
$$

where

$$
p(t)
$$

is the probability distribution over markings at time

$$
t.
$$

This would provide a direct connection between:

```text
SPN simulation
CTMC theory
Chemical Master Equation
Gillespie SSA
```

---

## 7.6 Parameter Sweeps and Ensemble Simulations

The current framework simulates one stochastic trajectory at a time.

A natural extension is to run many independent simulations with
different random seeds.

This would allow estimation of:

- probability of mono-allelic Xist activation,
- probability of bi-allelic activation,
- time to first Xist up-regulation,
- distribution of final states,
- sensitivity to parameters.

For example:

```python
for seed in range(1000):
    simulator = Simulator(net, seed=seed)
    result = simulator.run(...)
```

The output could then be summarized statistically.

---

## 7.7 Model Calibration

Another future direction is parameter calibration.

Given experimental measurements, one could estimate parameters such as:

$$
p_{21},p_{22},p_{23}
$$

or Hill thresholds such as:

$$
p_4,p_6,p_{12},p_{14}.
$$

Possible methods include:

- grid search,
- random search,
- Bayesian optimization,
- approximate Bayesian computation,
- likelihood-free inference.

This would make the framework useful not only for qualitative
simulation, but also for quantitative model fitting.

---

## 7.8 Alternative Simulation Methods

The Gillespie SSA is exact but can become computationally expensive
when the number of firing events is large.

Future versions may include approximate or hybrid methods such as:

- tau-leaping,
- chemical Langevin equation,
- hybrid deterministic/stochastic simulation,
- moment closure approximations,
- finite-state projection.

These methods would be useful when molecule numbers are large or when
many trajectories must be simulated.

---

## 7.9 Support for Other Biological Networks

Although the first application is the cXR–tXA model, the framework is
intended to support other biological systems.

Possible examples include:

- toggle switches,
- repressilators,
- signaling pathways,
- differentiation networks,
- cell-cycle models,
- epidemic models,
- metabolic networks.

The core simulator should remain unchanged when adding such models.

Only a new file inside

```text
spn/models/
```

should be required.

---

## 7.10 Documentation and Reproducible Experiments

Future work should also improve reproducibility.

Possible additions:

```text
docs/theory.md
docs/cxr_txa_model.md
docs/gillespie.md
docs/mpbn.md
notebooks/
examples/
tests/
```

Each numerical experiment should specify:

- parameter values,
- initial marking,
- simulation time,
- number of runs,
- random seed,
- classification rule for final states.

This would make the repository easier to use for future contributors.

---

# 8. References

This project is based on several mathematical, biological, and
computational references.

---

## 8.1 Biological Model

- Verena Mutzel et al.  
  **A symmetric toggle switch explains the onset of random X inactivation in different mammals.**  
  *Nature Structural & Molecular Biology*, 2019.  
  This is the main biological reference for the cXR–tXA model and the
  stochastic modeling of random X chromosome inactivation. :contentReference[oaicite:1]{index=1}

- Verena Mutzel et al.  
  **Supplementary Notes: Model Description.**  
  This document contains the detailed ODE model, stochastic simulation
  setup, parameter ranges, scaling factors, and Gillespie-based
  simulations of the cXR–tXA model. :contentReference[oaicite:2]{index=2}

---

## 8.2 Gillespie Algorithm and Stochastic Simulation

- Daniel T. Gillespie.  
  **Exact Stochastic Simulation of Coupled Chemical Reactions.**  
  *The Journal of Physical Chemistry*, 1977.  
  Classical reference for the stochastic simulation algorithm.

- Gillespie algorithm — Wikipedia.  
  Useful introductory reference for the main idea, history, and
  variants of the algorithm. :contentReference[oaicite:3]{index=3}

- Alberto Policriti.  
  **Stochastic Simulation and Gillespie’s Algorithm.**  
  Lecture notes introducing exponential waiting times, CTMCs, stochastic
  chemical systems, and rate functions. :contentReference[oaicite:4]{index=4}

---

## 8.3 Petri Nets and Stochastic Petri Nets

- Tadao Murata.  
  **Petri Nets: Properties, Analysis and Applications.**  
  *Proceedings of the IEEE*, 1989.

- Wolfgang Reisig.  
  **Petri Nets: An Introduction.**  
  Springer, 1985.

- M. K. Molloy.  
  **Performance Analysis Using Stochastic Petri Nets.**  
  *IEEE Transactions on Computers*, 1982.

- M. Ajmone Marsan, G. Balbo, and G. Conte.  
  **A Class of Generalized Stochastic Petri Nets for the Performance
  Evaluation of Multiprocessor Systems.**  
  *ACM Transactions on Computer Systems*, 1984.

- Falko Bause and Pieter S. Kritzinger.  
  **Stochastic Petri Nets: An Introduction to the Theory.**

Several of these classical Petri net and SPN references are summarized
in the Gillespie/SPN notes used for this project. :contentReference[oaicite:5]{index=5}

---

## 8.4 ODE, CPN, and Petri Net Conversion

- Hadi Fouani.  
  **From Ordinary Differential Equations to Continuous Petri Nets: A
  Structural Method for ODE–CPN Conversion.**  
  Working document, 2026.  
  This document formalizes the relation

  $$
  \dot m(t)=C v(m(t))
  $$

  and explains how an ODE vector field can be decomposed into
  non-negative transition rates and integer stoichiometric directions. :contentReference[oaicite:6]{index=6}

- Laura Recalde, Serge Haddad, and Manuel Silva.  
  **Continuous Petri Nets: Expressive Power and Decidability Issues.**  
  Reference on continuous Petri nets, reachability, liveness,
  deadlock-freeness, and timed continuous Petri net semantics. :contentReference[oaicite:7]{index=7}

---

## 8.5 Boolean Networks, MPBNs, and Petri Nets

- Stefan Haar and Juri Kolčák.  
  **Continuous Petri Nets Faithfully Fluidify Most Permissive Boolean
  Networks.**  
  CMSB, 2025.  
  This work studies the formal relation between Continuous Petri Nets
  and Most Permissive Boolean Networks. :contentReference[oaicite:8]{index=8}

- Thomas Chatain, Stefan Haar, Juraj Kolčák, Loïc Paulevé, and Aalok
  Thakkar.  
  **Concurrency in Boolean Networks.**  
  *Natural Computing*, 2020.  
  Related work on Boolean network semantics and concurrency.

---

## 8.6 PADE, Piecewise-Affine Models, and Qualitative Analysis

- Hidde de Jong et al.  
  **Qualitative Simulation of Genetic Regulatory Networks.**  
  This work develops qualitative simulation methods for genetic
  regulatory networks described by piecewise-linear differential
  equations. :contentReference[oaicite:9]{index=9}

- Jean-Luc Gouzé and Madalena Chaves.  
  **Piecewise Affine Models of Regulatory Genetic Networks.**  
  This reference reviews qualitative analysis of piecewise-affine
  differential equation models for genetic regulatory networks. :contentReference[oaicite:10]{index=10}

- Laurent Tournier and collaborators.  
  **Hierarchical Analysis of Piecewise Affine Models of Gene Regulatory
  Networks.**  
  This work studies hierarchical and qualitative analysis of
  piecewise-affine systems introduced for genetic regulatory networks. :contentReference[oaicite:11]{index=11}

---

## 8.7 Additional Useful References

- Monika Heiner, David Gilbert, and Robin Donaldson.  
  **Petri Nets for Systems and Synthetic Biology.**  
  Useful reference for the application of Petri nets in biological
  modeling.

- René David and Hassane Alla.  
  **Discrete, Continuous, and Hybrid Petri Nets.**  
  Springer, 2010.

- Javier Esparza and Mogens Nielsen.  
  **Decidability Issues for Petri Nets: A Survey.**

- C. Rackoff.  
  **The Covering and Boundedness Problems for Vector Addition Systems.**

- Ernst W. Mayr.  
  **An Algorithm for the General Petri Net Reachability Problem.**