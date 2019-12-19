---
title: 'EoN (Epidemics on Networks): a fast, flexible Python package for simulation, analytic approximation, and analysis of epidemics on networks'

tags:
  - SIR
  - SIS
  - Epidemic
  - Network
  - Stochastic simulation
  - Python
  
authors:
 - name: Joel C. Miller
   orcid: 0000-0003-4426-0405
   affiliation: "1, 2"
 - name: Tony Ting
   affiliation: "2"

affiliations:
 - name: La Trobe University, Melbourne, Australia
   index: 1
 - name: Institute for Disease Modeling, Seattle, Washington, USA
   index: 2
   
date: 19 July 2019

bibliography: paper.bib


---

# Summary

EoN (EpidemicsOnNetworks) is a pure-python package designed to assist studies of
infectious processes spreading through networks.  It originally rose out of the 
book *Mathematics of Epidemics on Networks* [@kiss:EoN], and now consists of over 
100 user-functions.


EoN provides a set of tools for

- Susceptible-Infected-Susceptible (SIS) and Susceptible-Infected-Recovered 
(SIR) disease
  - Stochastic simulation of disease spread in networkx graphs
    - continuous time Markovian
    - continuous time nonMarkovian
    - discrete time
  - Numerical solution of over 20 differential equation models, including
    - individual and pair-based models
    - pairwise models
    - edge-based compartmental models
- Stochastic simulation of a wide range of Simple and Complex contagions
- Visualization and analysis of stochastic simulations

These algorithms are built on the networkx package [@hagberg2008exploring].
EoN's documentation is maintained at 
https://epidemicsonnetworks.readthedocs.io/en/latest/ 
including numerous examples at 
https://epidemicsonnetworks.readthedocs.io/en/latest/Examples.html.  In this 
paper we provide brief descriptions with examples of a few of EoN's tools.
The examples shown are intended to demonstrate the ability of the tools.  The 
online documentation gives more detail about how to use them.

We model spreading processes on a contact network.  
## SIR and SIS disease

### Stochastic simulation
The stochastic SIR and SIS simulation tools allow the user to investigate
standard SIS and SIR dynamics (SEIR/SIRS and other processes are addressed 
within the simple contagion model):

- Markovian SIS and SIR simulations  (``fast_SIS``, ``Gillespie_SIS``, ``fast_SIR``, and ``Gillespie_SIR``).
- non-Markovian SIS and SIR simulations (``fast_nonMarkovian_SIS`` and ``fast_nonMarkovian_SIR``).
- discrete time SIS and SIR simulations where infections last a single time step 
  (``basic_discrete_SIS``, ``basic_discrete_SIR``, and ``discrete_SIR``).

For both Markovian and non-Markovian methods it is possible for the transition 
rates to depend on intrinsic properties of 
individuals and of partnerships.

The continuous-time stochastic simulations have two different implementations: a 
Gillespie implementation [@gillespie1977exact; @doob1945markoff] and an Event-driven
implementation.  Both approaches are efficient.  They have similar speed if the 
dynamics are Markovian (depending on the network and disease parameters either
may be faster than the other), but the event-driven implementation can also handle 
non-Markovian dynamics.  In earlier versions, the event-driven simulations were 
consistently faster than the Gillespie simulations, and thus they are named 
`fast_SIR` and `fast_SIS`.  The Gillespie simulations have since been optimized
using ideas from [@holme2014model] and [@cota2017optimized].

The algorithms can typically handle an SIR epidemic spreading on 
hundreds of thousands of individuals in well under a minute on a laptop.  The 
SIS versions are slower because the number of events that happen is often much
larger in an SIS simulation.

### Differential Equations Models

EoN also provides a set of tools to numerically solve approximately 20 differential equations
models for SIS or SIR disease spread in networks.  The various models use different
information about the network to make deterministic predictions about
the number infected at different times.  These use the Scipy integration tools.
The derivations of the models and explanations of their simplifying assumptions
are described in [@kiss:EoN].

Depending on the model, we need different information about the network structure.
The algorithms allow us to provide the information as inputs.  However, there
is also a version of each model which takes a network as an input 
instead and then measures the network properties.


## Simple and Complex Contagions

There are other contagious processes in networks which have received attention.
Many of these fall into one of two types, "simple contagions" and "complex 
contagions".

In a "simple contagion" an individual ``u`` may be induced to change status by
an interaction with its partner ``v``.  This status change occurs with the same
rate regardless of the statuses of other partners of ``u`` (although the other
partners may cause ``u`` to change to another status first).  SIS and SIR 
diseases are special cases of simple contagions.

In a "complex contagion" however, we permit the rate at which ``u`` changes 
from one status to another to depend on the statuses of others in some more
complicated way.  Two infected individuals may cause a susceptible individual
to become infected at some higher rate than would result from them acting 
independently.  This is frequently thought to model social contagions where an
individual may only believe something if multiple partners believe it 
[@centola:cascade].  

The simple and complex contagions are currently implemented only in a 
Gillespie setting, and so they require Markovian assumptions.  Although they 
are reasonably fast, it would typically be feasible to make a bespoke algorithm
that runs significantly faster.

### Simple contagions

EoN provides a function ``Gillespie_simple_contagion`` which allows a user to 
specify the rules governing an arbitrary simple contagion.

Examples are provided in the online documentation, including

- SEIR disease (there is an exposed state before becoming infectious)
- SIRS disease (recovered individuals eventually become susceptible again)
- SIRV disease (individuals may get vaccinated) 
- Competing SIR diseases (there is cross immunity)
- Cooperative SIR diseases (infection with one disease helps spread the other)

The implementation requires the user to separate out two distinct ways that 
transitions occur: those that are intrinsic to an individual's current state
and those that are induced by a partner.  To help demonstrate, consider an 
"SEIR" epidemic, where individuals begin susceptible, but when they interact 
with infectious partners they may enter an exposed state.  They remain in that 
exposed state for some period of time before transitioning into the infectious 
state independently of the status of any partner.
They remain infectious and eventually transition into the recovered state, again
independently of the status of any partner.  Here the "E" to "I" and "I" to "R"
transitions are intrinsic to the individual's state, while the "S" to "E" 
transition is induced by a partner.  

### Complex contagions

Complex contagions are implemented through ``Gillespie_complex_contagion`` which 
allows a user to specify the rules governing a relatively arbitrary complex 
contagion.  The one criteria we note is that there is no memory - an individual 
will change from one status to another based on the current statuses of its 
neighbors, and not based on previous interactions with some neighbors who may 
have since changed status.

In the Gillespie implementation, we need a user-defined function which 
calculates the rate at which ``u`` will change status (given knowledge about 
the current state of the system) and another user-defined function which
chooses the new status of ``u`` given that it is changing status.  We finally 
need a user-defined function that will determine which other nodes have their 
rate change due to ``u``'s transition.  By knowing the rates of all nodes
the Gillespie algorithm can choose the time of the next transition and which
node transitions.  Then it finds the new state, and finally it calculates the
new rates for all nodes affected by the change.

Once these functions are defined, the Gillespie algorithm is able to perform
the complex contagion simulation.


### Visualization & Analysis

By default the simulations return numpy arrays providing the number of individuals
with each state at each time.  However if we set a flag ``return_full_data=True``,
then the simulations return a ``Simulation_Investigation`` object.  With the
``Simulation_Investigation`` object, there are methods which allow us to 
reconstruct all details of the simulation.  We can know the exact status of 
each individual at each time, as well as who infected whom.  

There are also 
methods  provided to produce output from the ``Simulation_Investigation`` 
object.  These allow us to produce a snapshot of the network at a given time.
By default the visualization also includes the time series (e.g., S, I, and R) 
plotted beside the network snapshot.  These time series plots 
can be removed, or replaced by other time series, for example we could plot
multiple time series in the same axis, or time series generated by one of the
differential equations models.  With appropriate additional packages needed for
matplotlib's animation tools, the software can produce animations as well.

For SIR outbreaks, the ``Simulation_Investigation`` object includes a 
transmission tree.  For SIS and simple contagions, it includes a directed 
multigraph showing the transmissions that occurred (this may not be a tree).
However for complex contagions, we cannot determine who
is responsible for inducing a transition, so the implementation does not provide
a transmission tree.  The transmission tree is useful for constructing synthetic
phylogenies as in [@moshiri2018favites].


## Discussion

EoN provides a number of tools for studying infectious processes spreading in
contact networks.  The examples given here are intended to demonstrate the
range of EoN, but they represent only a fraction of the possibilities.

Full documentation is available at https://epidemicsonnetworks.readthedocs.io/en/latest/

# Dependencies:

scipy
numpy
networkx
matplotlib

# Related Packages

There are several alternative software packages that allow for simulation of 
epidemics on networks.  Here we briefly review some of these.


### epydemic

Epydemic is a python package that can simulate SIS and SIR epidemics in 
networks.  It is also built on networkx.  It can handle both discrete-time 
simulations or continuous-time Markovian simulations for which it uses a 
Gillespie-style algorithm.  It can handle more processes than just SIS or SIR 
disease.  In fact it can handle any model which can be simulated using the 
``EoN.simple_contagion``.

The documentation is available at https://pyepydemic.readthedocs.io/en/latest/


### Graph-tool

Graph-tool [@peixoto_graph-tool_2014] is a python package that serves as an 
alternative to networkx.  Many of its underlying processes are written in C++, 
so it is often much faster than networkx.

Graph-tool has a number of built-in dynamic models, including the SIS, SIR, 
and SIRS models.  The disease models are currently available only in 
discrete-time versions.

The documentation for these disease models is available at 
https://graph-tool.skewed.de/static/doc/dynamics.html.

### EpiModel

EpiModel [@jenness:EpiModel] is an R package that can handle SI, SIS, and SIR 
disease spread.  It is possible to extend EpiModel to other models.  EpiModel
is built around the StatNet package.  More details about EpiModel are available
at https://www.epimodel.org/



# Funding and Support

The development of EoN has been supported by Global Good and by La Trobe 
University.

# References
